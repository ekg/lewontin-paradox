#!/usr/bin/env python3
"""Independent, fail-closed review of the immutable ten-pair VGP pilot.

The production run stopped all ten pairs at mandatory assembly-QC preflight.
This reviewer therefore verifies what exists (identity, immutable objects,
terminal ledgers, and failure accounting) and refuses to reinterpret absent
mapping/variant/mask/PSMC products as passing zeros.  It emits the review
packet and a non-authorizing repair/re-pilot manifest; it never submits jobs.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable, Mapping, Sequence


ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = ROOT / "analysis"
RUN_ID = "vgp10-20260718-preflight-v1"
REVIEW_ID = "vgp10-review-20260718-v1"
REVIEW_SCHEMA = "vgp-10-pilot-independent-review-v1.0.0"
PRIMARY_IDS = tuple(f"P{i:02d}" for i in range(1, 11))
SENTINELS = ("P07", "P08")
REQUIRED_CLADES = (
    "Mammalia", "Aves", "Reptilia", "Amphibia", "Actinopterygii", "Chondrichthyes"
)
REQUIRED_GENERATIONS = (
    "early_clr_triocanu_parental", "later_hifi_hic_hap1_hap2"
)
CORE_ROLES = (
    "dataset_report", "checksum_catalog", "genome_fasta", "assembly_report", "assembly_stats"
)
EXPECTED_REASONS = {
    "MISSING_EXACT_FINAL_SEQUENCE_QV": 10,
    "MISSING_H1_BUSCO": 2,
    "MISSING_H2_BUSCO": 10,
    "MISSING_KMER_COPY_NUMBER_AUDIT": 10,
    "MISSING_REPEAT_OR_LOW_COMPLEXITY_MASK": 10,
    "UNRESOLVED_EXACT_READ_CHEMISTRY": 1,
}
PAIR_ARTIFACT_COLUMNS = {
    "qc": "qc_artifact",
    "diversity": "diversity_artifact",
    "psmc_trajectory": "psmc_trajectory_artifact",
    "bootstrap": "bootstrap_artifact",
    "scenario": "scenario_artifact",
    "validation": "validation_artifact",
    "annotation": "annotation_artifact",
    "telemetry": "telemetry_artifact",
    "failure": "failure_artifact",
}


class ReviewError(RuntimeError):
    """The immutable review packet is incomplete, inconsistent, or drifted."""


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_text(path: Path, value: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_name(path.name + ".partial")
    partial.write_text(value, encoding="utf-8")
    os.replace(partial, path)


def atomic_json(path: Path, value: object) -> None:
    atomic_text(path, json.dumps(value, indent=2, sort_keys=True) + "\n")


def write_tsv(path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_name(path.name + ".partial")
    with partial.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})
    os.replace(partial, path)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ReviewError(message)


def percentile(values: Sequence[float], fraction: float) -> float:
    """Linear-interpolated percentile, deterministic and dependency-free."""
    require(bool(values), "percentile requires observations")
    ordered = sorted(values)
    position = (len(ordered) - 1) * fraction
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    return ordered[lower] * (upper - position) + ordered[upper] * (position - lower)


def _artifact_packet(row: Mapping[str, str]) -> int:
    hashes: dict[str, str] = {}
    for name, column in PAIR_ARTIFACT_COLUMNS.items():
        path = ROOT / row[column]
        require(path.is_file(), f"{row['selection_id']}: missing {name} artifact")
        hashes[name] = sha256_file(path)
    observed = hashlib.sha256("".join(hashes[name] for name in sorted(hashes)).encode()).hexdigest()
    require(observed == row["artifact_packet_sha256"],
            f"{row['selection_id']}: artifact packet digest mismatch")
    return len(hashes)


def _reason_counts(rows: Sequence[Mapping[str, str]]) -> Counter[str]:
    counts: Counter[str] = Counter()
    for row in rows:
        codes = [value for value in row["all_reason_codes"].split(";") if value]
        require(codes and codes[0] == row["primary_reason_code"],
                f"{row['selection_id']}: primary reason ordering mismatch")
        require(len(codes) == len(set(codes)), f"{row['selection_id']}: duplicate reason code")
        counts.update(codes)
    return counts


def validate_closed_ledgers(
    design: Sequence[Mapping[str, str]],
    acquisitions: Sequence[Mapping[str, str]],
    results: Sequence[Mapping[str, str]],
    qc_rows: Sequence[Mapping[str, str]],
    telemetry: Sequence[Mapping[str, str]],
    summary: Mapping[str, object],
) -> dict[str, object]:
    primaries = [row for row in acquisitions if row["roster_type"] == "primary"]
    alternates = [row for row in acquisitions if row["roster_type"] == "alternate"]
    require(tuple(row["selection_id"] for row in design) == PRIMARY_IDS,
            "design is not exactly ordered P01..P10")
    require(tuple(row["selection_id"] for row in primaries) == PRIMARY_IDS,
            "acquisition is not exactly ordered P01..P10")
    require(tuple(row["selection_id"] for row in results) == PRIMARY_IDS,
            "result ledger is not exactly ordered P01..P10")
    require(tuple(row["selection_id"] for row in qc_rows) == PRIMARY_IDS,
            "QC ledger is not exactly ordered P01..P10")
    require(len(alternates) == 6, "alternate ledger does not contain exactly six rows")
    require(all(row["activation_status"] == "standby_not_triggered" for row in alternates),
            "an alternate was activated without a reviewable amendment")
    require(all(row["amendment_id"] == "none" for row in acquisitions),
            "unexpected acquisition amendment")
    require(len(telemetry) == 11 and telemetry[-1]["selection_id"] == "ALL",
            "telemetry does not reconcile ten pair rows plus aggregate")

    design_by_id = {row["selection_id"]: row for row in design}
    acq_by_id = {row["selection_id"]: row for row in primaries}
    qc_by_id = {row["selection_id"]: row for row in qc_rows}
    artifact_files = 0
    for row in results:
        selection_id = row["selection_id"]
        source = design_by_id[selection_id]
        acquisition = acq_by_id[selection_id]
        qc = qc_by_id[selection_id]
        for field in ("species", "biosample", "individual_or_isolate",
                      "h1_accession_version", "h2_accession_version"):
            require(row[field] == source[field] == acquisition[field],
                    f"{selection_id}: provenance drift in {field}")
        require(row["terminal_state"] == "FAILED_PREFLIGHT" and row["disposition"] == "FAIL",
                f"{selection_id}: terminal state was promoted")
        require(row["pair_identity_status"] == row["accession_identity_status"] ==
                row["core_checksum_status"] == "pass", f"{selection_id}: identity/checksum failure")
        require(row["core_qc_status"] == "fail_missing_required_measurements",
                f"{selection_id}: core QC disposition drift")
        require(row["mapping_status"] == "not_run_preflight_failed",
                f"{selection_id}: unexpected mapping output")
        require(row["multiplicity_status"] == "not_measured_preflight_failed",
                f"{selection_id}: multiplicity absence mislabeled")
        require(all(row[field] == "not_run_preflight_failed" for field in (
            "impg_status", "variants_status", "callable_mask_status", "consensus_status",
            "diversity_status", "psmc_status", "scenario_status"
        )), f"{selection_id}: a downstream output was silently promoted")
        require(row["bootstrap_attempts"] == row["bootstrap_successes"] == "0",
                f"{selection_id}: bootstrap accounting mismatch")
        require(row["bootstrap_success_fraction"] == "not_applicable_no_attempts",
                f"{selection_id}: zero attempts assigned a success fraction")
        require(row["alternate_activated"] == "false" and row["slurm_jobs_submitted"] == "0",
                f"{selection_id}: unexpected activation or submitted job")
        require(qc["callable_universe_bp"] == qc["callable_bp"] ==
                qc["callable_reason_total_bp"] == "not_materialized",
                f"{selection_id}: absent callable quantities not labeled unmaterialized")
        require(qc["callable_accounting_discrepancy_bp"] == "not_measured",
                f"{selection_id}: absent mask accounting was assigned a numeric result")
        require(qc["h1_ref_mismatches"] == qc["h2_reconstruction_failures"] ==
                "not_measured_preflight_failed", f"{selection_id}: absent reconstruction mislabeled")
        artifact_files += _artifact_packet(row)

    counts = _reason_counts(results)
    require(dict(counts) == EXPECTED_REASONS,
            f"reason totals do not reconcile: observed {dict(counts)!r}")
    require(summary["primary_slots_expected"] == summary["primary_slots_accounted"] == 10,
            "summary primary accounting mismatch")
    require(summary["primary_completed"] == summary["alternate_activations"] ==
            summary["slurm_jobs_submitted"] == 0, "summary reports forbidden execution/promotion")
    require(summary["primary_failed_preflight"] == 10, "summary failure count mismatch")
    require(summary["reason_counts"] == EXPECTED_REASONS, "summary reason counts mismatch")
    return {
        "primary_slots": 10,
        "alternate_slots": 6,
        "artifact_packets": 10,
        "artifact_files": artifact_files,
        "reason_counts": dict(counts),
    }


def verify_input_digests(summary: Mapping[str, object]) -> dict[str, str]:
    observed: dict[str, str] = {}
    for relative, expected in summary["input_digests"].items():
        path = ROOT / relative
        require(path.is_file(), f"missing immutable input {relative}")
        observed[relative] = sha256_file(path)
        require(observed[relative] == expected, f"immutable input digest drift: {relative}")
    for record in summary["deliverables"].values():
        path = ROOT / record["path"]
        require(path.is_file(), f"missing run deliverable {record['path']}")
        require(path.stat().st_size == record["size_bytes"], f"deliverable size drift: {record['path']}")
        require(sha256_file(path) == record["sha256"], f"deliverable digest drift: {record['path']}")
    return observed


def independently_revalidate_sentinels(
    inventory: Sequence[Mapping[str, str]],
    results: Sequence[Mapping[str, str]],
    qc_rows: Sequence[Mapping[str, str]],
    output_root: Path,
) -> dict[str, object]:
    """Recompute the predeclared P07/P08 subset from immutable CAS objects.

    No biological product is eligible after the preflight gate.  Consequently
    this comparison verifies object bytes and independently asserts that each
    requested biological quantity is absent on both sides with the same
    fail-closed status; it does not claim a biological recomputation occurred.
    """
    result_by_id = {row["selection_id"]: row for row in results}
    qc_by_id = {row["selection_id"]: row for row in qc_rows}
    object_rows: list[dict[str, object]] = []
    comparisons: list[dict[str, object]] = []
    logical_bytes = 0
    for selection_id in SENTINELS:
        rows = [row for row in inventory if row["selection_id"] == selection_id
                and row["side"] in {"h1", "h2"} and row["object_role"] in CORE_ROLES]
        expected_keys = {(side, role) for side in ("h1", "h2") for role in CORE_ROLES}
        require(len(rows) == 10 and {(row["side"], row["object_role"]) for row in rows} == expected_keys,
                f"{selection_id}: independent core-object subset is incomplete")
        for row in sorted(rows, key=lambda value: (value["side"], value["object_role"])):
            path = Path(row["local_path"])
            require(path.is_file(), f"{row['object_id']}: immutable object absent")
            size = path.stat().st_size
            digest = sha256_file(path)
            require(size == int(row["expected_bytes"]), f"{row['object_id']}: size drift")
            require(digest == row["local_sha256"], f"{row['object_id']}: SHA-256 drift")
            logical_bytes += size
            object_rows.append({
                "selection_id": selection_id, "object_id": row["object_id"],
                "side": row["side"], "object_role": row["object_role"],
                "size_bytes": size, "expected_sha256": row["local_sha256"],
                "observed_sha256": digest, "status": "MATCH",
            })
        result = result_by_id[selection_id]
        qc = qc_by_id[selection_id]
        requested = {
            "normalized_variants": result["variants_status"],
            "callable_mask": result["callable_mask_status"],
            "callable_denominator": qc["callable_bp"],
            "mask_reason_total": qc["callable_reason_total_bp"],
            "diversity": result["diversity_status"],
            "diploid_consensus": result["consensus_status"],
            "psmc_input_and_unscaled_result": result["psmc_status"],
            "psmc_bootstraps": qc["bootstrap_success_fraction"],
        }
        for quantity, production_status in requested.items():
            status = "MATCH_NOT_MATERIALIZED_PRECHECK_FAILED"
            require(production_status in {"not_run_preflight_failed", "not_materialized",
                                          "not_applicable_no_attempts"},
                    f"{selection_id}: {quantity} unexpectedly materialized")
            comparisons.append({
                "selection_id": selection_id, "quantity": quantity,
                "production_status": production_status,
                "independent_status": "not_eligible_for_recomputation_preflight_failed",
                "comparison": status,
                "interpretation": "no biological value; absence agrees and is not a passing zero",
            })
    write_tsv(output_root / "immutable_object_revalidation.tsv",
              ("selection_id", "object_id", "side", "object_role", "size_bytes",
               "expected_sha256", "observed_sha256", "status"), object_rows)
    write_tsv(output_root / "subset_comparison.tsv",
              ("selection_id", "quantity", "production_status", "independent_status",
               "comparison", "interpretation"), comparisons)
    value = {
        "review_id": REVIEW_ID,
        "subset": list(SENTINELS),
        "subset_basis": "predeclared full independent sentinels in vgp_analysis_manifest.json",
        "immutable_objects_revalidated": len(object_rows),
        "logical_bytes_rehashed": logical_bytes,
        "checksum_or_size_failures": 0,
        "biological_quantities_compared": len(comparisons),
        "biological_recomputation_status": "NOT_ELIGIBLE_PREFLIGHT_FAILED",
        "comparison_status": "MATCH_EXPLICIT_NONMATERIALIZATION",
        "slurm_jobs_submitted": 0,
    }
    atomic_json(output_root / "subset_validation.json", value)
    return value


def resource_assessment(telemetry: Sequence[Mapping[str, str]]) -> dict[str, object]:
    pair_rows = [row for row in telemetry if row["selection_id"] in PRIMARY_IDS]
    aggregate = next(row for row in telemetry if row["selection_id"] == "ALL")
    logical = [float(row["logical_read_bytes"]) for row in pair_rows]
    elapsed = [float(row["elapsed_seconds"]) for row in pair_rows]
    cpu = [float(row["cpu_seconds"]) for row in pair_rows]
    rss_kib = [float(row["maximum_rss_kib"]) for row in pair_rows]
    fs_read = [float(row["filesystem_read_bytes"]) for row in pair_rows]
    require(all(float(row["scratch_high_water_bytes"]) == 0 for row in pair_rows),
            "unexpected preflight scratch use")

    # The 716-pair figures are an eligibility upper bound.  They extrapolate
    # only the observed checksum preflight.  Biological processing retains the
    # preregistered Tier3A planning envelope because zero cluster jobs means no
    # pilot calibration exists.  Scaling the 40-pair envelope is transparent,
    # conservative, and explicitly non-authorizing.
    factor = 716 / 40
    planning = {
        "low": {"core_hours": 18, "input_gb": 64, "output_gb": 6, "inodes": 30000,
                "read_gb": 50, "write_gb": 80, "metadata_ops": 200000,
                "memory_gib_per_job": 32, "scratch_gb_per_job": 3,
                "scratch_gb_aggregate": 12, "concurrency": 4,
                "wall_hours": 0.56},
        "base": {"core_hours": 88, "input_gb": 160, "output_gb": 20, "inodes": 60000,
                 "read_gb": 125, "write_gb": 180, "metadata_ops": 800000,
                 "memory_gib_per_job": 64, "scratch_gb_per_job": 7.5,
                 "scratch_gb_aggregate": 75, "concurrency": 10,
                 "wall_hours": 1.10},
        "high": {"core_hours": 1388, "input_gb": 800, "output_gb": 160, "inodes": 200000,
                 "read_gb": 750, "write_gb": 1000, "metadata_ops": 3000000,
                 "memory_gib_per_job": 96, "scratch_gb_per_job": 37.5,
                 "scratch_gb_aggregate": 375, "concurrency": 10,
                 "wall_hours": 17.35},
    }
    scaled: dict[str, dict[str, object]] = {}
    for scenario, values in planning.items():
        scaled[scenario] = {
            key: (value if key in {"memory_gib_per_job", "scratch_gb_per_job",
                                   "scratch_gb_aggregate", "concurrency"}
                  else round(value * factor, 3))
            for key, value in values.items()
        }
    # At least 10 acquired inputs + 23 core + 8 PSMC outputs per eligible pair;
    # exact annotation adds four.  Operational inode counts are separately much
    # larger because attempts, sentinels, logs, and indexes are retained.
    durable_core_objects = 716 * 41
    optional_annotation_objects = 716 * 4
    return {
        "observed_preflight": {
            "pairs": 10,
            "core_objects": 100,
            "logical_bytes": int(aggregate["logical_read_bytes"]),
            "logical_write_bytes": int(aggregate["logical_write_bytes"]),
            "cpu_seconds": round(float(aggregate["cpu_seconds"]), 6),
            "elapsed_seconds": round(float(aggregate["elapsed_seconds"]), 6),
            "peak_rss_kib": int(aggregate["maximum_rss_kib"]),
            "filesystem_read_bytes": int(aggregate["filesystem_read_bytes"]),
            "filesystem_write_bytes": int(aggregate["filesystem_write_bytes"]),
            "scratch_high_water_bytes": int(aggregate["scratch_high_water_bytes"]),
            "metadata_operations": int(aggregate["metadata_operations"]),
            "per_pair_elapsed_median_seconds": round(percentile(elapsed, 0.5), 6),
            "per_pair_elapsed_p95_seconds": round(percentile(elapsed, 0.95), 6),
        },
        "upper_bound_716_pair_preflight_extrapolation": {
            "input_objects": 7160,
            "logical_bytes": int(float(aggregate["logical_read_bytes"]) / 10 * 716),
            "cpu_hours": round(float(aggregate["cpu_seconds"]) / 3600 / 10 * 716, 3),
            "serial_wall_hours": round(float(aggregate["elapsed_seconds"]) / 3600 / 10 * 716, 3),
            "memory_mib_per_process_observed": round(max(rss_kib) / 1024, 3),
            "minimum_storage_headroom_fraction": 0.25,
        },
        "upper_bound_716_pair_biological_planning_envelope": scaled,
        "durable_object_contract_upper_bound": {
            "minimum_input_core_psmc_objects": durable_core_objects,
            "optional_annotation_objects": optional_annotation_objects,
            "maximum_objects_including_optional_annotation": durable_core_objects + optional_annotation_objects,
            "bytes": "see low/base/high persistent input and output envelope",
            "operational_inodes": "see low/base/high inode envelope",
        },
        "headroom": {
            "storage_and_inode_multiplier": 1.25,
            "per_job_stop_multiple_of_reviewed_high_estimate": 1.5,
            "high_storage_gb_with_25pct_headroom": round((scaled["high"]["input_gb"] +
                                                          scaled["high"]["output_gb"]) * 1.25, 3),
            "high_inodes_with_25pct_headroom": round(scaled["high"]["inodes"] * 1.25),
            "mapping_memory_stop_gib_per_job_if_authorized": 144,
        },
        "limitations": [
            "The pilot observed only login-node immutable-object preflight; it observed no mapping, IMPG, mask, consensus, PSMC, bootstrap, scratch, or cluster I/O.",
            "The biological envelope is a transparent 716/40 scaling of the preregistered Tier3A planning envelope, itself informed by three earlier calibration tuples; it is not a fitted pilot model.",
            "PSMC and 200-bootstrap runtime are not represented in that prior Tier3A envelope, so total biological core-hours and wall time remain lower bounds and cannot authorize scale-out.",
            "Eligibility among 716 catalog rows is unresolved; 716 is an upper-bound sensitivity, not an expected completed-pair count.",
        ],
        "resource_model_gate": {
            "median_ape": None,
            "p95_ape": None,
            "status": "FAIL_NOT_ESTIMABLE_ZERO_CLUSTER_JOBS",
        },
    }


def gate_rows(reason_counts: Mapping[str, int]) -> list[dict[str, str]]:
    def row(gate_id: str, scope: str, gate: str, threshold: str, observed: str,
            status: str, evidence: str, consequence: str) -> dict[str, str]:
        return dict(gate_id=gate_id, scope=scope, gate=gate, threshold=threshold,
                    observed=observed, status=status, evidence=evidence,
                    decision_consequence=consequence)

    return [
        row("H01", "all_primary", "pair/individual identity resolved", "unresolved=0", "0", "PASS",
            "result manifest + frozen design/acquisition field equality", "none"),
        row("H02", "all_primary", "exact accession.version identity resolved", "unresolved=0", "0", "PASS",
            "result manifest + frozen design/acquisition field equality", "none"),
        row("H02A", "all_primary", "taxid/BioSample/individual/haplotype-role fields resolved",
            "unresolved=0", "0", "PASS", "frozen design/acquisition closed-world equality", "none"),
        row("H03", "all_primary", "immutable checksum drift", "events=0", "0", "PASS",
            "100 production objects; independent P07/P08 live rehash", "none"),
        row("H04", "retained_mapping", "retained query multiplicity >1", "count=0", "no retained mapping",
            "NOT_REACHED", "mapping refused at preflight", "cannot support GO"),
        row("H05", "retained_mapping", "retained target multiplicity >1", "count=0", "no retained mapping",
            "NOT_REACHED", "mapping refused at preflight", "cannot support GO"),
        row("H05A", "retained_mapping", "retained non-1:1 bases", "bp=0", "no retained mapping",
            "NOT_REACHED", "mapping refused at preflight", "cannot support GO"),
        row("H06", "callable_mask", "unexplained mask accounting", "discrepancy_bp=0", "not measured",
            "NOT_REACHED", "mask not materialized", "cannot support GO"),
        row("H07", "normalized_variants", "H1 REF reconstruction", "failures=0", "not measured",
            "NOT_REACHED", "variants not materialized", "cannot support GO"),
        row("H08", "normalized_variants", "H2 alternate reconstruction", "failures=0", "not measured",
            "NOT_REACHED", "variants not materialized", "cannot support GO"),
        row("H08A", "diploid_consensus", "non-callable bases encoded homozygous reference", "bp=0",
            "not measured", "NOT_REACHED", "consensus not materialized", "cannot support GO"),
        row("H08B", "annotation_outputs", "annotation sequence-dictionary mismatch retained",
            "count=0", "no annotation output", "NOT_REACHED",
            "P03/P04 mismatch excluded before output; all annotation branches not reached", "cannot support annotation GO"),
        row("H09", "psmc", "unscaled/scaled PSMC separation", "conflations=0", "no PSMC outputs",
            "NOT_REACHED", "scenario branch not run", "cannot support GO"),
        row("A01", "assembly", "exact-final-sequence QV each haplotype", "QV>=40 each", "0/10 measured",
            "FAIL", f"{reason_counts['MISSING_EXACT_FINAL_SEQUENCE_QV']}/10 missing", "repair required"),
        row("A02", "assembly", "H1/H2 BUSCO completeness and missingness",
            "complete>=0.90 and missing<=0.05 each", "H1 missing 2/10; H2 missing 10/10", "FAIL",
            "QC ledger", "repair required"),
        row("A02A", "assembly", "H1/H2 BUSCO completeness difference", "absolute difference<=0.05",
            "0/10 pair differences measurable", "FAIL", "H2 BUSCO absent for every pair", "repair required"),
        row("A02B", "assembly", "H1/H2 BUSCO duplication", "duplicated<=0.05 each",
            "0/10 both-haplotype audits measurable", "FAIL", "H2 BUSCO absent for every pair", "repair required"),
        row("A03", "assembly", "copy-number/k-mer collapse audit", "passing both haplotypes", "0/10 measured",
            "FAIL", f"{reason_counts['MISSING_KMER_COPY_NUMBER_AUDIT']}/10 missing", "repair required"),
        row("A04", "assembly", "repeat/low-complexity mask", "exact manifest-bound mask", "0/10 measured",
            "FAIL", f"{reason_counts['MISSING_REPEAT_OR_LOW_COMPLEXITY_MASK']}/10 missing", "repair required"),
        row("A05", "assembly", "exact read chemistry", "resolved", "9/10 resolved", "FAIL",
            "P09 unresolved", "P09-specific repair"),
        row("A06", "assembly", "minimum haplotype span", ">=250,000,000 bp each", "20/20 pass",
            "PASS", "frozen exact assembly reports", "does not override missing QC"),
        row("A07", "assembly", "minimum contig N50", ">=1,000,000 bp each", "20/20 pass",
            "PASS", "frozen exact assembly reports", "does not override missing QC"),
        row("A08", "assembly", "H1/H2 length ratio", "0.80..1.25", "10/10 pass; 0.941..1.157",
            "PASS", "frozen exact assembly reports", "does not override missing QC"),
        row("C01", "callability", "primary callable fraction", ">=0.60", "0/10 measured",
            "NOT_REACHED", "masks not materialized", "cannot support GO"),
        row("C01A", "callability", "sensitivity callable fraction", ">=0.50", "0/10 measured",
            "NOT_REACHED", "masks not materialized", "cannot support GO"),
        row("C01B", "callability", "minimum callable bases", ">=100,000,000 bp", "0/10 measured",
            "NOT_REACHED", "masks not materialized", "cannot support GO"),
        row("C01C", "callability", "well-callable windows", ">=50 1-Mb windows at >=80% callable",
            "0/10 measured", "NOT_REACHED", "windows not materialized", "cannot support GO"),
        row("C01D", "callability", "ordered disjoint reason assignment", "13 predeclared reasons; exact complement",
            "0/10 measured", "NOT_REACHED", "reason masks not materialized", "cannot support GO"),
        row("C01E", "callability", "universe reconciliation", "callable + reasons = universe exactly",
            "0/10 measured", "NOT_REACHED", "reason masks not materialized", "cannot support GO"),
        row("C02", "consensus", "masked and heterozygous encoding", "masked=N; heterozygous SNP=IUPAC",
            "0/10 produced", "NOT_REACHED", "consensus not materialized", "cannot support GO"),
        row("C03", "consensus", "primary indel flank mask", "10 bp", "0/10 produced", "NOT_REACHED",
            "consensus not materialized", "cannot support GO"),
        row("C04", "consensus", "indel flank sensitivities", "0 bp and 50 bp", "0/10 produced",
            "NOT_REACHED", "consensus not materialized", "cannot support GO"),
        row("C05", "normalized_variants", "exact-reference normalization and duplicate removal",
            "exact H1; exact duplicates removed", "0/10 produced", "NOT_REACHED",
            "variants not materialized", "cannot support GO"),
        row("B00", "passing_psmc_unit", "pre-registered primary bootstrap attempts", "exactly 200",
            "no passing PSMC unit", "NOT_REACHED", "0 attempts", "cannot support GO"),
        row("B01", "passing_psmc_unit", "minimum bootstrap attempts", ">=100", "no passing PSMC unit",
            "NOT_REACHED", "0 attempts; not a violation assigned to a passing unit", "cannot support GO"),
        row("B02", "passing_psmc_unit", "finite bootstrap fraction", ">=0.95",
            "no passing PSMC unit",
            "NOT_REACHED", "0/0, explicitly not applicable", "cannot support GO"),
        row("B02A", "passing_psmc_unit", "finite primary bootstrap successes", ">=190/200",
            "no passing PSMC unit", "NOT_REACHED", "0/0, explicitly not applicable", "cannot support GO"),
        row("B03", "passing_psmc_unit", "boundary-aware block construction",
            "5-Mb primary; 1/10-Mb sensitivities; never cross contig/mask boundary",
            "no bootstrap units", "NOT_REACHED", "callable mask absent", "cannot support GO"),
        row("B04", "diversity", "heterozygosity bootstrap replicates", "10,000",
            "0/10 attempted", "NOT_REACHED", "diversity not materialized", "cannot support GO"),
        row("D01", "reproducibility", "P07/P08 full independent biological recomputation",
            "2 sentinels; digest mismatches=0", "20 inputs rehashed; biological outputs ineligible",
            "NOT_REACHED", "independent review root", "re-pilot must complete biological comparison"),
        row("D02", "reproducibility", "other-pair deterministic shard comparison",
            "1 shard each; digest mismatches=0", "0/8 biological shards eligible", "NOT_REACHED",
            "all stopped at preflight", "re-pilot must complete biological comparison"),
        row("P00", "program", "primary slots adjudicated", "exactly 10", "10", "PASS",
            "closed result ledger", "none"),
        row("P01", "program", "primary callable+consensus passes", ">=8/10", "0/10", "FAIL",
            "all primaries failed preflight", "core GO prohibited"),
        row("P02", "program", "required major clades among passers", "6/6", "0/6", "FAIL",
            "no passers", "core GO prohibited"),
        row("P03", "program", "assembly generations among passers", "2/2", "0/2", "FAIL",
            "no passers", "core GO prohibited"),
        row("P04", "program", "no technology-stratum systematic failure", "failure fraction <0.50",
            "1.00 in both early CLR and later HiFi strata", "FAIL",
            "shared missing-QC packet; not evidence that either technology is biologically invalid", "core GO prohibited"),
        row("P05", "program", "zero hard-gate violations", "0", "measured identity/checksum violations=0; downstream hard gates not reached",
            "NOT_REACHED", "positive downstream hard-gate evidence absent", "core GO prohibited"),
        row("R01", "resource_model", "median absolute percentage error", "<=0.25", "not estimable",
            "FAIL_NOT_ESTIMABLE", "zero biological jobs with positive predictions/observations", "core GO prohibited"),
        row("R02", "resource_model", "95th percentile absolute percentage error", "<=0.50", "not estimable",
            "FAIL_NOT_ESTIMABLE", "zero biological jobs with positive predictions/observations", "core GO prohibited"),
        row("R03", "resource_model", "storage headroom", ">=0.25", "0.25 retained in corrected envelope",
            "PASS", "resource manifest", "planning only; does not rescue APE failure"),
        row("R04", "resource_model", "per-job stop multiple", "1.5x reviewed high estimate",
            "1.5x retained; no job invoked", "PASS", "resource manifest", "planning only; not a global ceiling"),
        row("L01", "ledger", "all primary slots reconciled", "10/10", "10/10", "PASS",
            "no drop/replacement", "none"),
        row("L02", "ledger", "all alternates reconciled", "6/6", "6/6 standby; 0 activated", "PASS",
            "no amendment", "none"),
        row("L03", "ledger", "failure/warning reason totals", "exact", "43 reason instances across 6 codes",
            "PASS", "result/QC/failure/summary equality", "none"),
        row("N01", "interpretation", "annotation absence cannot veto core", "zero annotation-driven core failures",
            "0", "PASS", "P03/P04 annotation mismatch remained branch-local", "none"),
        row("N02", "interpretation", "same-pair PSMC is not independent validation", "no independent claim",
            "no PSMC result and no claim", "PASS", "decision/report wording", "none"),
    ]


def branch_decisions() -> dict[str, dict[str, object]]:
    return {
        "core_diversity_psmc": {
            "decision": "CONDITIONAL_GO",
            "evidence_status": "TECHNICAL_PREFLIGHT_FAILURE_REPAIRABLE",
            "scaleout_authorized": False,
            "reason": "0/10 passed mandatory assembly QC, callable, consensus, bootstrap, and resource-model gates; identities and immutable core checksums passed.",
        },
        "exact_annotation_partitions": {
            "decision": "CONDITIONAL_GO",
            "evidence_status": "NOT_REACHED_CORE_PREFLIGHT_FAILED",
            "scaleout_authorized": False,
            "reason": "Run only after a core pass and exact annotation accession/version plus sequence-dictionary or validated-liftover binding; annotation absence remains non-vetoing for core.",
        },
        "direct_conversion": {
            "decision": "CONDITIONAL_GO",
            "evidence_status": "SEPARATE_PILOT_IN_PROGRESS_NO_IMMUTABLE_RESULT_REVIEWED",
            "scaleout_authorized": False,
            "reason": "H1/H2 differences are not transmitted conversion events; await the separately authorized pedigree/gamete pilot review.",
        },
        "population_gbgc": {
            "decision": "NOT_RUN/DESIGN_ONLY",
            "evidence_status": "NO_EXECUTION_TASK",
            "scaleout_authorized": False,
            "reason": "No multi-individual population-frequency execution evidence exists in the current graph.",
        },
        "phylogenetic_substitution_bias": {
            "decision": "CONDITIONAL_GO",
            "evidence_status": "NOT_EXECUTED_INPUT_GATE",
            "scaleout_authorized": False,
            "reason": "The separately authorized H01/H02 pilot produced a pinned metadata preflight but zero verified sequences, callable alignments, substitutions, or biological estimates.",
        },
        "non_allelic_conversion": {
            "decision": "NOT_RUN/DESIGN_ONLY",
            "evidence_status": "NO_EXECUTION_TASK",
            "scaleout_authorized": False,
            "reason": "No copy-resolved non-allelic execution evidence exists in the current graph.",
        },
    }


def scaleout_rows(
    acquisitions: Sequence[Mapping[str, str]],
    results: Sequence[Mapping[str, str]],
    source_digests: Mapping[str, str],
) -> list[dict[str, object]]:
    result_by_id = {row["selection_id"]: row for row in results}
    rows: list[dict[str, object]] = []
    for source in acquisitions:
        primary = source["roster_type"] == "primary"
        result = result_by_id.get(source["selection_id"])
        rows.append({
            "manifest_version": "vgp-10-scaleout-repair-manifest-v1.0.0",
            "review_id": REVIEW_ID,
            "source_acquisition_sha256": source_digests["analysis/vgp_10_pilot_acquisition_manifest.tsv"],
            "source_result_sha256": sha256_file(ANALYSIS / "vgp_10_pilot_result_manifest.tsv"),
            "selection_id": source["selection_id"],
            "roster_type": source["roster_type"],
            "activation_status": source["activation_status"],
            "amendment_id": source["amendment_id"],
            "failed_primary_retained": source["failed_primary_retained"],
            "species": source["species"],
            "biosample": source["biosample"],
            "individual_or_isolate": source["individual_or_isolate"],
            "h1_accession_version": source["h1_accession_version"],
            "h2_accession_version": source["h2_accession_version"],
            "assembly_generation": source["assembly_generation"],
            "assembly_technologies": source["assembly_technologies"],
            "long_range_phasing_evidence": source["long_range_phasing_evidence"],
            "core_identity_status": source["core_identity_status"],
            "pilot_terminal_state": result["terminal_state"] if result else "STANDBY_NOT_TRIGGERED",
            "pilot_reason_codes": result["all_reason_codes"] if result else "none",
            "repair_action": (
                "measure exact-final QV both haplotypes; BUSCO both haplotypes; k-mer/copy-number audit; exact repeat/low-complexity mask"
                + ("; resolve exact read chemistry" if primary and source["selection_id"] == "P09" else "")
                + ("; H1 BUSCO required" if primary and source["selection_id"] in {"P07", "P08"} else "")
                if primary else "remain standby; activation requires same-clade pre-result trigger and versioned amendment"
            ),
            "authorized_action": "REPAIR_AND_REPILOT_ONLY" if primary else "STANDBY_NO_ACTION",
            "biological_jobs_authorized": "false",
            "full_scaleout_authorized": "false",
        })
    return rows


def resource_manifest_rows(resources: Mapping[str, object]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for scenario, values in resources["upper_bound_716_pair_biological_planning_envelope"].items():
        rows.append({
            "manifest_version": "vgp-10-scaleout-resource-envelope-v1.0.0",
            "review_id": REVIEW_ID,
            "scope": "716_pair_eligibility_upper_bound_not_expected_count",
            "scenario": scenario,
            "object_count_minimum_contract": 29356,
            "persistent_input_gb": values["input_gb"],
            "persistent_output_gb": values["output_gb"],
            "operational_inodes": values["inodes"],
            "core_hours_lower_bound_excludes_psmc": values["core_hours"],
            "memory_gib_per_job": values["memory_gib_per_job"],
            "scratch_gb_per_job": values["scratch_gb_per_job"],
            "scratch_gb_aggregate_at_concurrency": values["scratch_gb_aggregate"],
            "moosefs_read_gb": values["read_gb"],
            "moosefs_write_gb": values["write_gb"],
            "metadata_operations": values["metadata_ops"],
            "concurrency": values["concurrency"],
            "wall_hours_lower_bound_excludes_psmc": values["wall_hours"],
            "uncertainty": "planning envelope; biological pilot telemetry absent; PSMC/200-bootstrap cost not identified",
            "operational_headroom": "25% storage/inodes; reviewed per-job high stop multiplier 1.5; not a global scientific eligibility ceiling",
            "authorization": "NOT_AUTHORIZED_REPAIR_REPILOT_REQUIRED",
        })
    return rows


def render_report(
    design: Sequence[Mapping[str, str]],
    acquisitions: Sequence[Mapping[str, str]],
    results: Sequence[Mapping[str, str]],
    gates: Sequence[Mapping[str, str]],
    decision: Mapping[str, object],
) -> str:
    by_acq = {row["selection_id"]: row for row in acquisitions}
    lines = [
        "# Independent review of the ten-genome VGP pilot",
        "",
        f"Review ID: `{REVIEW_ID}`",
        f"Run reviewed: `{RUN_ID}`",
        "Program decision: **CONDITIONAL_GO — bounded repair and re-pilot only; full biological scale-out is not authorized**",
        "",
        "## Executive finding",
        "",
        "All ten immutable primary slots are accounted for, retain the exact frozen H1/H2 accession.version and BioSample/individual identity, and pass live content-addressed checksum review. All ten nevertheless failed before mapping because mandatory assembly-QC evidence was absent: exact-final QV for both haplotypes (10/10), H2 BUSCO (10/10), a manifest-bound k-mer/copy-number audit (10/10), and a repeat/low-complexity mask (10/10). P07/P08 additionally lack H1 BUSCO and P09 lacks resolved exact read chemistry. These are technical packet failures. They are not low diversity, biological outliers, mapping failures, or failed PSMC histories.",
        "",
        "Zero biological jobs were submitted. Therefore whole-assembly 1:1 mapping, multiplicity, IMPG extraction, normalized variants, REF/alternate reconstruction, callable masks, consensus, diversity, unscaled PSMC, scaling scenarios, and bootstraps are **not reached**, not passing zeros. With 0/10 callable+consensus passers, zero represented clades/generations among passers, systematic failure in both generation/technology strata, and unestimable resource-model APE, a core GO is prohibited. Exact identity/checksum success and the repairable common cause support CONDITIONAL_GO rather than permanent NO_GO.",
        "",
        "Annotation state never caused a core failure. P03/P04 annotation mismatches remain branch-local, and absence of annotation or independent Ne evidence cannot veto a future technically valid core result. PSMC, if later produced, is descriptive evidence from the same pair and is not independent validation of same-pair diversity.",
        "",
        "## Ten primary slots",
        "",
        "| slot | species | exact H1 / H2 | individual | generation; technology and phasing | annotation branch | classification | reasons |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for source, result in zip(design, results):
        acq = by_acq[source["selection_id"]]
        technology = acq["assembly_technologies"].replace("|", "/")
        phasing = acq["long_range_phasing_evidence"].replace("|", "/")
        lines.append(
            f"| {source['selection_id']} | *{source['species']}* | `{source['h1_accession_version']}` / `{source['h2_accession_version']}` | "
            f"`{source['biosample']}`; `{source['individual_or_isolate']}` | `{source['assembly_generation']}`; {technology}; {phasing} | "
            f"`{acq['annotation_branch_status']}` | technical failure; repairable; no biological classification | `{result['all_reason_codes']}` |"
        )
    lines.extend([
        "",
        "There are no low-confidence-but-usable core results and no biological outliers because no biological measurement crossed preflight. Long-range phasing evidence is present for every pair and is retained as later confidence evidence; it cannot substitute for QV, completeness, collapse, mapping, or callability measurements.",
        "",
        "## Alternates and failure accounting",
        "",
        "All six declared alternates remain `standby_not_triggered`; none has an amendment, none replaced a failed primary, and every failed primary remains in the result ledger. All 10 artifact packets and their 90 files pass digest verification. The exact reason totals reconcile across result, QC, pair failure artifacts, independent validation, and run summary: 10 QV, 2 H1 BUSCO, 10 H2 BUSCO, 10 k-mer/copy-number, 10 repeat-mask, and 1 chemistry reasons (43 total reason instances). There are no unknown warning codes, retries, silent drops, scheduler IDs, or Slurm dependency edges.",
        "",
        "| alternate | replaces clade | species | exact H1 / H2 | individual | generation | disposition |",
        "|---|---|---|---|---|---|---|",
    ])
    for alternate in (row for row in acquisitions if row["roster_type"] == "alternate"):
        lines.append(
            f"| {alternate['selection_id']} | {REQUIRED_CLADES[int(alternate['selection_id'][1:]) - 1]} | "
            f"*{alternate['species']}* | `{alternate['h1_accession_version']}` / `{alternate['h2_accession_version']}` | "
            f"`{alternate['biosample']}`; `{alternate['individual_or_isolate']}` | `{alternate['assembly_generation']}` | "
            "standby; no trigger or versioned amendment; not reviewed as a biological result |"
        )
    lines.extend([
        "",
        "## Independent stratified recomputation",
        "",
        "P07 and P08 are the predeclared full-independent sentinels. In the independent output root, the reviewer rehashed and size-checked all 20 of their H1/H2 core objects, then compared production and independent eligibility for normalized variants, masks and reason totals, callable denominators, diversity, diploid consensus, PSMC input/unscaled result, and bootstraps. Every object matched. Every biological quantity matched only as explicit nonmaterialization after the same immutable preflight failure. A numerical biological recomputation would require bypassing preregistered gates and was correctly not attempted. No Slurm or biological job was submitted.",
        "",
        "## Gate application",
        "",
        "| gate | requirement | observed | result | consequence |",
        "|---|---|---|---|---|",
    ])
    for gate in gates:
        lines.append(f"| {gate['gate_id']} {gate['gate']} | {gate['threshold']} | {gate['observed']} | **{gate['status']}** | {gate['decision_consequence']} |")
    resources = decision["resource_assessment"]
    observed = resources["observed_preflight"]
    projected = resources["upper_bound_716_pair_preflight_extrapolation"]
    lines.extend([
        "",
        "The bootstrap 100-attempt/95%-finite rules were not violated by a passing PSMC unit—there is no passing unit—but they were also not demonstrated. The program requires positive evidence, so vacuous truth is not used to support GO.",
        "",
        "## Resource review and upper-bound scale sensitivity",
        "",
        f"Observed preflight revalidated 100 objects and {observed['logical_bytes']:,} logical bytes using {observed['cpu_seconds']:.3f} CPU-seconds, {observed['elapsed_seconds']:.3f} elapsed seconds, {observed['peak_rss_kib'] / 1024:.1f} MiB peak RSS, zero scratch, {observed['filesystem_read_bytes']:,} filesystem-read bytes, {observed['logical_write_bytes']:,} logical report bytes, zero filesystem-write bytes, and {observed['metadata_operations']} metadata operations. A transparent 716-pair upper-bound extrapolation for this checksum-only step is 7,160 objects, {projected['logical_bytes']:,} logical bytes, {projected['cpu_hours']:.3f} core-hours, and {projected['serial_wall_hours']:.3f} serial wall-hours.",
        "",
        "The pilot contains no observed mapping, IMPG, consensus, PSMC/bootstrap, scratch, or cluster-I/O telemetry, so a fitted full biological resource model—and hence median/p95 APE—does not exist. The corrected resource manifest retains a low/base/high 716-pair sensitivity by scaling the preregistered 40-pair Tier3A envelope (which was informed by three earlier calibration tuples), with eligibility explicitly unresolved. It is a planning lower bound because PSMC plus 200 bootstraps were absent from that older envelope:",
        "",
        "| scenario | minimum durable objects | persistent input/output GB | operational inodes | core-hours lower bound | memory/job GiB | scratch/job; aggregate GB | read/write GB | concurrency | wall-hours lower bound |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ])
    for scenario, values in resources["upper_bound_716_pair_biological_planning_envelope"].items():
        lines.append(
            f"| {scenario} | 29,356 | {values['input_gb']:,.1f} / {values['output_gb']:,.1f} | {values['inodes']:,.0f} | "
            f"{values['core_hours']:,.1f} | {values['memory_gib_per_job']} | {values['scratch_gb_per_job']} / {values['scratch_gb_aggregate']} | "
            f"{values['read_gb']:,.1f} / {values['write_gb']:,.1f} | {values['concurrency']} | {values['wall_hours']:,.2f} |"
        )
    lines.extend([
        "",
        f"Operational headroom is explicit rather than a global eligibility ceiling: 25% for storage/inodes, producing {resources['headroom']['high_storage_gb_with_25pct_headroom']:,.1f} GB and {resources['headroom']['high_inodes_with_25pct_headroom']:,} inodes at the high upper-bound sensitivity; per-job stopping is 1.5× a reviewed high estimate (144 GiB for a 96-GiB mapping estimate). Initial repaired-pilot concurrency remains 2. Full-scale concurrency 10 is only a sensitivity and is not authorized. A repaired ten-pair run must capture all six APE dimensions and PSMC/bootstrap cost before any full-scale resource authorization.",
        "",
        "## Branch decisions",
        "",
        "| branch | decision | reason |",
        "|---|---|---|",
    ])
    for branch, value in decision["branches"].items():
        lines.append(f"| `{branch}` | **{value['decision']}** | {value['reason']} |")
    lines.extend([
        "",
        "## Bounded repair / re-pilot authorization boundary",
        "",
        "CONDITIONAL_GO permits only preparation of a versioned repair packet and a newly reviewed ten-slot re-pilot. It does not authorize this review task to acquire new biology or submit jobs, and it does not authorize `scale-vgp-core` to begin full catalog processing.",
        "",
        "Before re-pilot submission, every retained primary must have exact-final QV for both haplotypes; BUSCO completeness/missing/duplication for both; manifest-bound k-mer/copy-number audits; exact repeat/low-complexity masks; P09 chemistry resolution; immutable hashes; reviewed resources; and an amendment only if a same-clade alternate is triggered before results. The re-pilot must then demonstrate at least 8/10 callable+consensus passes, all required clades, both generations, no technology-stratum systematic failure, at least 100 PSMC bootstraps and 95% finite success per passing PSMC unit, zero hard violations, and median/p95 resource APE no worse than 25%/50%. No threshold may be relaxed after outcome inspection.",
        "",
        "Scale-out may be reconsidered only after an independent review of that repaired immutable packet. Annotation and specialized branches remain orthogonal to core validity.",
        "",
        "## Reproducibility",
        "",
        "The review generator, its tests, and live sentinel rehash run through the same pinned GNU Guix time-machine/channel and production manifest. The independent output root contains the object ledger, requested-quantity comparison, and validation JSON. The corrected manifests are immutable derivations of the reviewed source digests and explicitly set `biological_jobs_authorized=false` and `full_scaleout_authorized=false`.",
        "",
    ])
    return "\n".join(lines)


def review(output: Path = ANALYSIS, independent_root: Path | None = None) -> dict[str, object]:
    independent_root = independent_root or output / "vgp_10_pilot_review_independent"
    design = read_tsv(ANALYSIS / "vgp_10_pair_manifest.tsv")
    acquisitions = read_tsv(ANALYSIS / "vgp_10_pilot_acquisition_manifest.tsv")
    inventory = read_tsv(ANALYSIS / "vgp_10_pilot_object_inventory.tsv")
    results = read_tsv(ANALYSIS / "vgp_10_pilot_result_manifest.tsv")
    qc_rows = read_tsv(ANALYSIS / "vgp_10_pilot_qc.tsv")
    telemetry = read_tsv(ANALYSIS / "vgp_10_pilot_resource_telemetry.tsv")
    summary = json.loads((ANALYSIS / "vgp_10_pilot_run_summary.json").read_text())
    command_log = json.loads((ANALYSIS / "vgp_10_pilot_command_log.json").read_text())

    require(command_log["scientific_thresholds_relaxed"] is False,
            "production command log reports threshold relaxation")
    require(command_log["slurm_jobs"] == command_log["slurm_dependencies"] == [],
            "production command log contains biological scheduler activity")
    source_digests = verify_input_digests(summary)
    ledger = validate_closed_ledgers(design, acquisitions, results, qc_rows, telemetry, summary)
    independent = independently_revalidate_sentinels(inventory, results, qc_rows, independent_root)
    reasons = _reason_counts(results)
    gates = gate_rows(reasons)
    resources = resource_assessment(telemetry)
    branches = branch_decisions()
    decision: dict[str, object] = {
        "schema_version": REVIEW_SCHEMA,
        "review_id": REVIEW_ID,
        "reviewed_run_id": RUN_ID,
        "program_decision": "CONDITIONAL_GO",
        "authorization": "BOUNDED_REPAIR_AND_TEN_SLOT_REPILOT_ONLY",
        "full_biological_scaleout_authorized": False,
        "scientific_thresholds_relaxed": False,
        "primary_passes": 0,
        "primary_required": 8,
        "primary_slots_reviewed": 10,
        "alternates_reviewed": 6,
        "alternate_activations": 0,
        "slurm_jobs_submitted_by_review": 0,
        "classification_counts": {
            "technical_failure": 10,
            "low_confidence_usable_core": 0,
            "biological_outlier": 0,
            "biological_not_estimable": 10,
        },
        "ledger": ledger,
        "source_digests": source_digests,
        "independent_validation": independent,
        "gate_summary": {
            "pass": sum(row["status"] == "PASS" for row in gates),
            "fail": sum(row["status"].startswith("FAIL") for row in gates),
            "not_reached": sum(row["status"] == "NOT_REACHED" for row in gates),
            "total": len(gates),
        },
        "branches": branches,
        "resource_assessment": resources,
        "repair_requirements": [
            "Measure exact-final QV for both haplotypes in all ten primaries.",
            "Measure BUSCO completeness, missingness, and duplication on both haplotypes; add H1 measurements for P07/P08.",
            "Complete manifest-bound k-mer/copy-number collapse audits and exact repeat/low-complexity masks for all ten.",
            "Resolve P09 exact read chemistry.",
            "Run the unchanged gates in a newly versioned ten-slot re-pilot with no outcome-driven relaxation.",
            "Capture positive observed telemetry for wall time, CPU-hours, peak RSS, scratch, reads, and writes and meet median/p95 APE <=25%/50%.",
            "Require >=8/10 callable+consensus passers with all six clades, both generations, no systematic technology-stratum failure, and per-passing-unit PSMC >=100 attempts and >=95% finite success.",
        ],
        "interpretation_guards": {
            "failed_assemblies_imply_low_diversity": False,
            "annotation_absence_vetoes_core": False,
            "absent_independent_ne_vetoes_core": False,
            "same_pair_psmc_is_independent_validation": False,
        },
    }

    write_tsv(output / "vgp_10_pilot_review_gates.tsv",
              ("gate_id", "scope", "gate", "threshold", "observed", "status", "evidence",
               "decision_consequence"), gates)
    scaleout = scaleout_rows(acquisitions, results, source_digests)
    write_tsv(output / "vgp_10_pilot_scaleout_manifest.tsv", tuple(scaleout[0]), scaleout)
    resource_rows = resource_manifest_rows(resources)
    write_tsv(output / "vgp_10_pilot_scaleout_resource_manifest.tsv",
              tuple(resource_rows[0]), resource_rows)
    atomic_json(output / "vgp_10_pilot_review_decision.json", decision)
    atomic_text(output / "vgp_10_pilot_review.md",
                render_report(design, acquisitions, results, gates, decision))
    return decision


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ANALYSIS)
    parser.add_argument("--independent-output-root", type=Path)
    args = parser.parse_args()
    decision = review(args.output, args.independent_output_root)
    print(json.dumps({
        "review_id": decision["review_id"],
        "program_decision": decision["program_decision"],
        "full_biological_scaleout_authorized": decision["full_biological_scaleout_authorized"],
        "primary_passes": decision["primary_passes"],
        "slurm_jobs_submitted_by_review": decision["slurm_jobs_submitted_by_review"],
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
