#!/usr/bin/env python3
"""Build the closed-world VGP Freeze 1 core scale-out ledger.

This program deliberately separates *preparing* a complete scale-out ledger
from authorizing biological execution.  The reviewed ten-pair pilot currently
permits only a bounded repair and re-pilot, and the Freeze 1 mirror currently
contains no verified payload.  Consequently the committed invocation accounts
for every frozen row and every catalog-linked haplotype, but submits no Slurm
job and materializes no biological estimate.

Future execution must use a new, review-bound manifest.  This program never
turns a CONDITIONAL_GO repair packet into a full-scale GO and never consults a
moving assembly release.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import os
import re
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Mapping, Sequence


ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = ROOT / "analysis"
CATALOG = Path(
    "/moosefs/erikg/lewontin-paradox-data/vgp/phase1-freeze-1.0/manifests/"
    "VGPPhase1-freeze-1.0.commit-dc1b2af5a7741b97d66fb10cb2bce97f41765cdf.tsv"
)
CATALOG_COMMIT = "dc1b2af5a7741b97d66fb10cb2bce97f41765cdf"
CATALOG_SHA256 = "9c58420484a8b76a2d6175b7c26bf709e68bdc726a67fc7541b8c2b5a2fc13a4"
CATALOG_BYTES = 327_466
CATALOG_ROWS = 716
CATALOG_LINES = 717
ACCESSION_RE = re.compile(r"^GC[AF]_\d{9}\.\d+$")

DECISION = ANALYSIS / "vgp_10_pilot_review_decision.json"
REVIEW_SCALEOUT = ANALYSIS / "vgp_10_pilot_scaleout_manifest.tsv"
REVIEW_RESOURCES = ANALYSIS / "vgp_10_pilot_scaleout_resource_manifest.tsv"
MIRROR_SUMMARY = ANALYSIS / "vgp_freeze1_mirror_summary.json"
MIRROR_MANIFEST = ANALYSIS / "vgp_freeze1_mirror_manifest.tsv"
DESIGN = ANALYSIS / "vgp_analysis_manifest.json"
CHANNELS = ANALYSIS / "guix/vgp_10_pilot/channels.scm"
GUIX_MANIFEST = ANALYSIS / "guix/vgp_10_pilot/manifest.scm"

MANIFEST_FIELDS = (
    "manifest_version", "catalog_commit", "catalog_sha256", "record_id",
    "record_type", "catalog_row", "scientific_name", "lineage", "taxid",
    "catalog_status", "sex", "main_assembly_id", "h1_accession_version",
    "ucsc_mirror_accession_version", "linked_ordinal", "link_class",
    "linked_assembly_id", "h2_accession_version", "pair_id",
    "catalog_assembly_technology", "catalog_qv", "long_range_phase_signal",
    "exact_individual_provenance_status", "mirror_h1_state", "mirror_h2_state",
    "disposition", "primary_reason_code", "blocking_evidence_codes",
    "confidence_tier", "wave_id", "core_result_status", "psmc_status",
    "annotation_partition_status", "annotation_absence_core_veto",
    "same_pair_psmc_independent_evidence", "biological_jobs_authorized",
    "source_review_decision_sha256", "source_mirror_summary_sha256",
)

QC_FIELDS = (
    "manifest_version", "pair_id", "catalog_row", "link_class",
    "pair_identity", "same_individual_provenance", "accession_version_frozen",
    "h1_mirror_verified", "h2_mirror_verified", "assembly_generation",
    "hifi_or_base_accuracy", "qv_each_haplotype", "completeness_each_haplotype",
    "duplication_collapse", "long_range_phase_confidence", "mutually_comparable",
    "whole_assembly_1to1_multiplicity", "mask_accounting", "callable_denominator",
    "consensus", "psmc_unscaled", "psmc_bootstrap_attempts",
    "psmc_boundary_aware", "annotation_binding", "core_eligible",
    "disposition", "primary_reason_code", "annotation_absence_core_veto",
    "same_pair_psmc_independent_evidence",
)

TELEMETRY_FIELDS = (
    "telemetry_id", "scope", "source", "measurement_status", "objects",
    "logical_bytes", "cpu_seconds", "elapsed_seconds", "peak_rss_kib",
    "scratch_bytes", "filesystem_read_bytes", "filesystem_write_bytes",
    "slurm_jobs_submitted", "pairs_completed", "reason_code", "source_sha256",
)

WAVE_FIELDS = (
    "wave_manifest_version", "wave_id", "record_type", "scenario",
    "authorization", "pair_count", "max_pair_jobs", "max_array_concurrency",
    "memory_gib_per_job", "scratch_gb_per_job", "wall_hours_lower_bound_excludes_psmc",
    "per_job_wall_limit_status",
    "max_transient_retries", "max_resource_reestimate_retries",
    "max_scientific_relaxation_retries", "job_stop_multiple",
    "checkpoint_contract", "atomic_promotion_contract", "input_digest_contract",
    "output_digest_contract", "environment_digest", "stop_condition",
)


class ScaleoutError(RuntimeError):
    """A frozen-input, authorization, or closed-world invariant failed."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def atomic_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.partial-{os.getpid()}")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        handle.write(text)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def atomic_json(path: Path, value: Mapping[str, object]) -> None:
    atomic_text(path, json.dumps(value, indent=2, sort_keys=True) + "\n")


def atomic_tsv(path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, object]]) -> None:
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    atomic_text(path, buffer.getvalue())


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def split_catalog_list(value: str) -> list[str]:
    return [item.strip() for item in value.split(",") if item.strip()]


def verify_catalog(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise ScaleoutError(f"pinned catalog missing: {path}")
    observed = (path.stat().st_size, sha256_file(path), sum(1 for _ in path.open("rb")))
    expected = (CATALOG_BYTES, CATALOG_SHA256, CATALOG_LINES)
    if observed != expected:
        raise ScaleoutError(f"frozen catalog identity drift: observed={observed}, expected={expected}")
    rows = read_tsv(path)
    if len(rows) != CATALOG_ROWS:
        raise ScaleoutError(f"frozen catalog row drift: {len(rows)} != {CATALOG_ROWS}")
    return rows


def verify_authorization(decision_path: Path, review_manifest_path: Path) -> dict[str, object]:
    decision = json.loads(decision_path.read_text(encoding="utf-8"))
    core = decision.get("branches", {}).get("core_diversity_psmc", {})
    expected = (
        decision.get("program_decision"),
        decision.get("authorization"),
        decision.get("full_biological_scaleout_authorized"),
        core.get("decision"),
        core.get("scaleout_authorized"),
    )
    if expected != (
        "CONDITIONAL_GO", "BOUNDED_REPAIR_AND_TEN_SLOT_REPILOT_ONLY", False,
        "CONDITIONAL_GO", False,
    ):
        raise ScaleoutError(f"review authorization boundary drift: {expected}")
    review_rows = read_tsv(review_manifest_path)
    if len(review_rows) != 16:
        raise ScaleoutError("review scale-out roster must retain 10 primaries and 6 alternates")
    if any(
        row["authorized_action"] not in {"REPAIR_AND_REPILOT_ONLY", "STANDBY_NO_ACTION"}
        or row["biological_jobs_authorized"] != "false"
        or row["full_scaleout_authorized"] != "false"
        for row in review_rows
    ):
        raise ScaleoutError("review scale-out manifest does not uniformly refuse full execution")
    return decision


def verify_mirror(summary_path: Path, manifest_path: Path) -> tuple[dict[str, object], dict[str, str]]:
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    reconciliation = summary.get("catalog_reconciliation", {})
    if reconciliation != {
        "catalog_rows": 716, "closed_world": True, "expected_accession_roots": 581,
        "extra_roots": 0, "missing_roots": 0, "observed_accession_roots": 581,
        "released_rows": 581, "unreleased_rows": 135, "accession_or_version_drift": 0,
    }:
        raise ScaleoutError("mirror catalog reconciliation drift")
    if summary.get("bulk_launch") != {
        "launched": False, "reason": "quota_visibility_unavailable_fail_closed",
        "slurm_jobs_launched": 0,
    }:
        raise ScaleoutError("mirror launch state drift; a new reviewed handoff is required")
    rows = read_tsv(manifest_path)
    if len(rows) != 47_870 or {row["state"] for row in rows} != {"planned"}:
        raise ScaleoutError("mirror is not the reviewed all-planned 47,870-object snapshot")
    states: dict[str, str] = {}
    for row in rows:
        accession = row["accession_version"]
        state = row["state"]
        old = states.setdefault(accession, state)
        if old != state:
            raise ScaleoutError(f"mixed mirror state for {accession}")
    return summary, states


def phase_signal(technology: str) -> str:
    lower = technology.lower()
    if "trio" in lower:
        return "TRIO_CATALOG_SIGNAL_CONFIDENCE_SUPPORT_ONLY"
    if "hi-c" in lower or "hic" in lower or "omni-c" in lower:
        return "HIC_CATALOG_SIGNAL_CONFIDENCE_SUPPORT_ONLY"
    return "NO_LONG_RANGE_SIGNAL_IN_CATALOG_NOT_A_CORE_VETO"


def linked_records(row: Mapping[str, str]) -> list[tuple[str, int, str, str, bool]]:
    result: list[tuple[str, int, str, str, bool]] = []
    for link_class, id_column, accession_column in (
        ("other_high_quality", "Assembly IDs other high-quality haplotypes", "Accession #s other high-quality haplotypes"),
        ("alternate", "Assembly IDs alternate haplotypes", "Accession #s alternate haplotypes"),
    ):
        identifiers = split_catalog_list(row[id_column])
        accessions = split_catalog_list(row[accession_column])
        cardinality_match = len(identifiers) == len(accessions)
        for ordinal, accession in enumerate(accessions, 1):
            identifier = identifiers[ordinal - 1] if cardinality_match else "UNRESOLVED_CATALOG_CARDINALITY"
            result.append((link_class, ordinal, identifier, accession, cardinality_match))
    return result


def build_ledgers(
    catalog_rows: Sequence[Mapping[str, str]], mirror_states: Mapping[str, str],
    review_digest: str, mirror_digest: str,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    manifest: list[dict[str, object]] = []
    qc: list[dict[str, object]] = []
    seen_pair_ids: set[str] = set()
    linked_count = 0

    for catalog_row, row in enumerate(catalog_rows, 2):
        h1 = row["Accession # for main haplotype"].strip()
        ucsc = row["UCSC Browser main haplotype"].strip()
        main_id = row["Assembly ID main haplotype"].strip()
        links = linked_records(row)
        for accession in ([h1] if h1 else []) + [entry[3] for entry in links] + ([ucsc] if ucsc else []):
            if not ACCESSION_RE.fullmatch(accession):
                raise ScaleoutError(f"catalog row {catalog_row} has invalid accession/version {accession!r}")

        if not h1:
            row_disposition, row_reason = "excluded", "CATALOG_NO_RELEASED_MAIN_ACCESSION"
        elif not links:
            row_disposition, row_reason = "excluded", "CATALOG_NO_LINKED_HAPLOTYPE"
        elif all(h1 == link[3] for link in links):
            row_disposition, row_reason = "excluded", "CATALOG_ONLY_SELF_LINK_NOT_PAIR"
        else:
            row_disposition, row_reason = "failed", "UPSTREAM_SCALEOUT_NOT_AUTHORIZED"

        common = {
            "manifest_version": "vgp-core-freeze1-scaleout-v1.0.0",
            "catalog_commit": CATALOG_COMMIT,
            "catalog_sha256": CATALOG_SHA256,
            "catalog_row": catalog_row,
            "scientific_name": row["Scientific Name"].strip(),
            "lineage": row["Lineage"].strip(),
            "taxid": row["NCBI taxon ID"].strip(),
            "catalog_status": row["Status"].strip() or "EMPTY",
            "sex": row["sex"].strip() or "UNRESOLVED",
            "main_assembly_id": main_id or "UNRELEASED",
            "h1_accession_version": h1 or "UNRELEASED",
            "ucsc_mirror_accession_version": ucsc or "UNRELEASED",
            "catalog_assembly_technology": row["Assembly tech"].strip() or "UNRESOLVED",
            "catalog_qv": row["QV"].strip() or "UNMEASURED",
            "long_range_phase_signal": phase_signal(row["Assembly tech"]),
            "source_review_decision_sha256": review_digest,
            "source_mirror_summary_sha256": mirror_digest,
            "annotation_absence_core_veto": "false",
            "same_pair_psmc_independent_evidence": "false",
            "biological_jobs_authorized": "false",
        }
        manifest.append({
            **common, "record_id": f"ROW-{catalog_row:04d}", "record_type": "catalog_row",
            "linked_ordinal": 0, "link_class": "none", "linked_assembly_id": "",
            "h2_accession_version": "", "pair_id": "", "exact_individual_provenance_status": "NOT_AUDITED",
            "mirror_h1_state": mirror_states.get(ucsc, "not_in_frozen_mirror") if ucsc else "not_released",
            "mirror_h2_state": "not_applicable", "disposition": row_disposition,
            "primary_reason_code": row_reason,
            "blocking_evidence_codes": "MIRROR_ZERO_VERIFIED_PAYLOAD;PILOT_ZERO_CORE_PASSES",
            "confidence_tier": "X" if row_disposition == "excluded" else "UNASSIGNED",
            "wave_id": "NONE", "core_result_status": "NOT_RUN",
            "psmc_status": "NOT_RUN", "annotation_partition_status": "NOT_RUN_CORE_NOT_PASSED",
        })

        class_ordinals = Counter()
        for link_class, _within_class, linked_id, h2, cardinality_match in links:
            linked_count += 1
            class_ordinals[link_class] += 1
            tag = "HQ" if link_class == "other_high_quality" else "ALT"
            pair_id = f"R{catalog_row:04d}-{tag}{class_ordinals[link_class]:02d}"
            if pair_id in seen_pair_ids:
                raise ScaleoutError(f"duplicate pair id: {pair_id}")
            seen_pair_ids.add(pair_id)
            is_self = h1 == h2
            disposition = "excluded" if is_self else "failed"
            reason = "CATALOG_SELF_LINK_NOT_PAIR" if is_self else "UPSTREAM_SCALEOUT_NOT_AUTHORIZED"
            evidence = ["MIRROR_ZERO_VERIFIED_PAYLOAD", "PILOT_ZERO_CORE_PASSES"]
            if not cardinality_match:
                evidence.append("CATALOG_ASSEMBLY_ID_ACCESSION_CARDINALITY_MISMATCH")
            if not is_self:
                evidence.extend((
                    "EXACT_INDIVIDUAL_PROVENANCE_NOT_AUDITED", "EXACT_FINAL_QV_NOT_MEASURED",
                    "COMPLETENESS_NOT_MEASURED", "DUPLICATION_COLLAPSE_NOT_MEASURED",
                    "WHOLE_ASSEMBLY_1TO1_NOT_RUN", "CALLABILITY_NOT_MEASURED",
                ))
            manifest.append({
                **common, "record_id": f"PAIR-{pair_id}", "record_type": "linked_pair",
                "linked_ordinal": linked_count, "link_class": link_class,
                "linked_assembly_id": linked_id, "h2_accession_version": h2,
                "pair_id": pair_id, "exact_individual_provenance_status": "NOT_AUDITED",
                "mirror_h1_state": mirror_states.get(ucsc, "not_in_frozen_mirror") if ucsc else "not_released",
                "mirror_h2_state": mirror_states.get(h2, "not_in_frozen_mirror"),
                "disposition": disposition, "primary_reason_code": reason,
                "blocking_evidence_codes": ";".join(evidence),
                "confidence_tier": "X" if is_self else "UNASSIGNED",
                "wave_id": "NONE", "core_result_status": "NOT_RUN",
                "psmc_status": "NOT_RUN", "annotation_partition_status": "NOT_RUN_CORE_NOT_PASSED",
            })
            catalog_phase = phase_signal(row["Assembly tech"])
            qc.append({
                "manifest_version": "vgp-core-freeze1-scaleout-qc-v1.0.0", "pair_id": pair_id,
                "catalog_row": catalog_row, "link_class": link_class,
                "pair_identity": "FAIL_SELF_LINK" if is_self else "NOT_AUDITED_EXACT",
                "same_individual_provenance": "NOT_AUDITED_EXACT", "accession_version_frozen": "PASS_CATALOG",
                "h1_mirror_verified": "FAIL_ZERO_VERIFIED", "h2_mirror_verified": "FAIL_NOT_IN_VERIFIED_MIRROR",
                "assembly_generation": "NOT_AUDITED", "hifi_or_base_accuracy": "NOT_AUDITED_EXACT",
                "qv_each_haplotype": "NOT_MEASURED", "completeness_each_haplotype": "NOT_MEASURED",
                "duplication_collapse": "NOT_MEASURED", "long_range_phase_confidence": catalog_phase,
                "mutually_comparable": "NOT_AUDITED", "whole_assembly_1to1_multiplicity": "NOT_RUN",
                "mask_accounting": "NOT_RUN", "callable_denominator": "NOT_RUN", "consensus": "NOT_RUN",
                "psmc_unscaled": "NOT_RUN", "psmc_bootstrap_attempts": 0,
                "psmc_boundary_aware": "NOT_RUN", "annotation_binding": "NOT_AUDITED_OPTIONAL",
                "core_eligible": "false", "disposition": disposition, "primary_reason_code": reason,
                "annotation_absence_core_veto": "false", "same_pair_psmc_independent_evidence": "false",
            })

    if len(manifest) != CATALOG_ROWS + linked_count or linked_count != 569:
        raise ScaleoutError(
            f"closed-world linked ledger drift: rows={len(manifest)}, linked={linked_count}; expected 1285/569"
        )
    if len(qc) != linked_count:
        raise ScaleoutError("pair/QC multiplicity is not 1:1")
    return manifest, qc


def build_telemetry(
    mirror_summary: Mapping[str, object], resources: Sequence[Mapping[str, str]],
    input_digests: Mapping[str, str],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = [
        {
            "telemetry_id": "SCALEOUT-GENERATOR", "scope": "closed_world_manifest_only",
            "source": "this_run", "measurement_status": "OBSERVED", "objects": 1_285,
            "logical_bytes": "", "cpu_seconds": "", "elapsed_seconds": "", "peak_rss_kib": "",
            "scratch_bytes": 0, "filesystem_read_bytes": "", "filesystem_write_bytes": "",
            "slurm_jobs_submitted": 0, "pairs_completed": 0,
            "reason_code": "UPSTREAM_SCALEOUT_NOT_AUTHORIZED",
            "source_sha256": input_digests["review_decision"],
        },
        {
            "telemetry_id": "FREEZE1-MIRROR", "scope": "47,870_frozen_objects",
            "source": "analysis/vgp_freeze1_mirror_summary.json", "measurement_status": "OBSERVED",
            "objects": 47_870,
            "logical_bytes": mirror_summary["inventory_totals"]["full_release"]["bytes"],
            "cpu_seconds": "", "elapsed_seconds": "", "peak_rss_kib": "", "scratch_bytes": 0,
            "filesystem_read_bytes": 0, "filesystem_write_bytes": 0,
            "slurm_jobs_submitted": 0, "pairs_completed": 0,
            "reason_code": "MIRROR_QUOTA_VISIBILITY_UNAVAILABLE", "source_sha256": input_digests["mirror_summary"],
        },
    ]
    for resource in resources:
        rows.append({
            "telemetry_id": f"PLANNING-{resource['scenario'].upper()}",
            "scope": resource["scope"], "source": "analysis/vgp_10_pilot_scaleout_resource_manifest.tsv",
            "measurement_status": "PLANNING_NOT_OBSERVED_NOT_AUTHORIZED",
            "objects": resource["object_count_minimum_contract"],
            "logical_bytes": "", "cpu_seconds": "", "elapsed_seconds": "", "peak_rss_kib": "",
            "scratch_bytes": "", "filesystem_read_bytes": "", "filesystem_write_bytes": "",
            "slurm_jobs_submitted": 0, "pairs_completed": 0,
            "reason_code": "BIOLOGICAL_RESOURCE_MODEL_NOT_ESTIMABLE",
            "source_sha256": input_digests["review_resources"],
        })
    return rows


def build_waves(resources: Sequence[Mapping[str, str]], environment_digest: str) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = [{
        "wave_manifest_version": "vgp-core-wave-plan-v1.0.0", "wave_id": "GATE-0",
        "record_type": "administrative_stop", "scenario": "current_reviewed_state",
        "authorization": "NOT_AUTHORIZED", "pair_count": 0, "max_pair_jobs": 0,
        "max_array_concurrency": 0, "memory_gib_per_job": 0, "scratch_gb_per_job": 0,
        "wall_hours_lower_bound_excludes_psmc": 0, "per_job_wall_limit_status": "ZERO_JOB_GATE",
        "max_transient_retries": 0,
        "max_resource_reestimate_retries": 0, "max_scientific_relaxation_retries": 0,
        "job_stop_multiple": 0, "checkpoint_contract": "NO_JOB_NO_CHECKPOINT",
        "atomic_promotion_contract": "NO_PROMOTION_WITHOUT_COMPLETE_DIGESTED_PAIR_PACKET",
        "input_digest_contract": "CATALOG_REVIEW_MIRROR_ENVIRONMENT_ALL_EXACT",
        "output_digest_contract": "PAIR_PACKET_SHA256_BEFORE_AND_AFTER_ATOMIC_RENAME",
        "environment_digest": environment_digest,
        "stop_condition": "CONDITIONAL_GO_REPAIR_ONLY_OR_MIRROR_NOT_VERIFIED",
    }]
    for resource in resources:
        rows.append({
            "wave_manifest_version": "vgp-core-wave-plan-v1.0.0",
            "wave_id": f"TEMPLATE-{resource['scenario'].upper()}", "record_type": "planning_template",
            "scenario": resource["scenario"], "authorization": "NOT_AUTHORIZED_PLANNING_ONLY",
            "pair_count": 0, "max_pair_jobs": 25,
            "max_array_concurrency": resource["concurrency"],
            "memory_gib_per_job": resource["memory_gib_per_job"],
            "scratch_gb_per_job": resource["scratch_gb_per_job"],
            "wall_hours_lower_bound_excludes_psmc": resource["wall_hours_lower_bound_excludes_psmc"],
            "per_job_wall_limit_status": "REQUIRES_REVIEWED_BIOLOGICAL_TELEMETRY_BEFORE_AUTHORIZATION",
            "max_transient_retries": 2, "max_resource_reestimate_retries": 1,
            "max_scientific_relaxation_retries": 0, "job_stop_multiple": 1.5,
            "checkpoint_contract": "STAGE_DIGEST_CHECKPOINT_THEN_RESUMABLE_RETRY",
            "atomic_promotion_contract": "COMPLETE_DIGESTED_PAIR_PACKET_RENAME_ONLY",
            "input_digest_contract": "PER_OBJECT_SHA256_AND_FROZEN_ACCESSION_VERSION",
            "output_digest_contract": "ALL_OUTPUT_SHA256_PLUS_COMPLETENESS_GATE",
            "environment_digest": environment_digest,
            "stop_condition": "FINITE_25_PAIR_WAVE_OR_ANY_HARD_GATE_OR_1.5X_REVIEWED_HIGH_ESTIMATE",
        })
    return rows


def paper_tables(manifest: Sequence[Mapping[str, object]]) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    pairs = [row for row in manifest if row["record_type"] == "linked_pair"]
    pair_rows = [{
        "pair_id": row["pair_id"], "catalog_row": row["catalog_row"],
        "scientific_name": row["scientific_name"], "link_class": row["link_class"],
        "h1_accession_version": row["h1_accession_version"],
        "h2_accession_version": row["h2_accession_version"], "disposition": row["disposition"],
        "primary_reason_code": row["primary_reason_code"], "core_callable_diversity": "NOT_ESTIMABLE",
        "psmc_unscaled": "NOT_ESTIMABLE", "annotation_partitions": "NOT_RUN",
        "same_pair_psmc_independent_evidence": "false",
    } for row in pairs]
    counts = Counter((row["record_type"], row["disposition"], row["primary_reason_code"]) for row in manifest)
    summary = [{
        "record_type": key[0], "disposition": key[1], "primary_reason_code": key[2],
        "count": value, "biological_estimates": 0,
        "interpretation": "operational/metadata disposition only; no biological estimate",
        "same_pair_psmc_independent_evidence": "false",
    } for key, value in sorted(counts.items())]
    return pair_rows, summary


def sensitivity_rows() -> list[dict[str, str]]:
    return [{
        "dimension": dimension, "status": "NOT_ESTIMABLE_ZERO_APPROVED_COMPLETED_PAIRS",
        "estimate": "NA", "reason_code": "UPSTREAM_SCALEOUT_NOT_AUTHORIZED",
        "interpretation": "No threshold effect may be inferred from technical non-execution.",
    } for dimension in (
        "QV", "completeness", "duplication_or_collapse", "long_range_phasing_evidence",
        "callability", "genome_size", "assembly_generation", "bootstrap_block_size",
        "bootstrap_finite_fraction", "mutation_rate_scenario", "generation_time_scenario",
    )]


def scaling_rows() -> list[dict[str, str]]:
    return [
        {
            "scenario_id": "UNSCALED_PRIMARY", "confidence_tier": "NONE_APPROVED",
            "mutation_rate_per_generation": "NOT_APPLICABLE", "mutation_rate_source": "NOT_APPLICABLE",
            "generation_time_years": "NOT_APPLICABLE", "generation_time_source": "NOT_APPLICABLE",
            "trajectory_status": "NOT_MATERIALIZED_CORE_NOT_RUN",
            "reporting_rule": "Unscaled PSMC must remain the primary reproducible object.",
        },
        {
            "scenario_id": "SPECIES_SCENARIOS_REQUIRED_AFTER_REVIEW", "confidence_tier": "NONE_APPROVED",
            "mutation_rate_per_generation": "NOT_SELECTED", "mutation_rate_source": "REQUIRES_MANIFEST_BOUND_SOURCE",
            "generation_time_years": "NOT_SELECTED", "generation_time_source": "REQUIRES_MANIFEST_BOUND_SOURCE",
            "trajectory_status": "NOT_MATERIALIZED_NO_APPROVED_TIER",
            "reporting_rule": "Never invent or pool scaling values; report every approved scenario separately.",
        },
    ]


def independent_validation_rows() -> list[dict[str, object]]:
    return [
        {
            "validation_component": component,
            "status": status,
            "selected_eligible_pairs": 0,
            "completed_checks": 0,
            "slurm_jobs_submitted": 0,
            "reason_code": reason,
            "interpretation": interpretation,
        }
        for component, status, reason, interpretation in (
            (
                "stratified_sample_across_clade_generation_genome_size_diversity_confidence",
                "EMPTY_BY_PREDECLARED_GATE", "ZERO_ELIGIBLE_PAIRS",
                "No post-hoc sample may be drawn from unaudited or failed pair candidates.",
            ),
            (
                "independent_biological_recomputation", "NOT_PERFORMED",
                "CORE_NOT_RUN_AND_SCALEOUT_NOT_AUTHORIZED",
                "Nonmaterialization agrees with the reviewed gate; it is not a biological replication.",
            ),
            (
                "raw_read_check", "NOT_PERFORMED", "RAW_READS_NOT_IN_VERIFIED_FREEZE1_MIRROR",
                "Raw reads cannot be inferred from assembly products or silently acquired outside the manifest.",
            ),
            (
                "kmer_copy_number_check", "NOT_PERFORMED", "EXACT_READ_KMERS_NOT_MATERIALIZED",
                "Copy-number/collapse evidence remains a mandatory unmeasured core gate.",
            ),
            (
                "literature_triangulation", "NOT_PERFORMED", "ZERO_ELIGIBLE_STRATIFIED_SAMPLE",
                "Literature values must be source-bound after eligible sample selection, never selected post hoc.",
            ),
        )
    ]


def render_results(
    manifest: Sequence[Mapping[str, object]], qc: Sequence[Mapping[str, object]],
    input_digests: Mapping[str, str], output_digests: Mapping[str, str], generated: str,
) -> str:
    pairs = [row for row in manifest if row["record_type"] == "linked_pair"]
    rows = [row for row in manifest if row["record_type"] == "catalog_row"]
    pair_counts = Counter(row["disposition"] for row in pairs)
    row_counts = Counter(row["disposition"] for row in rows)
    return f"""# VGP Freeze 1 core scale-out: closed-world gated result

Generated: `{generated}`. Release: VGP Phase 1 Freeze 1 commit
`{CATALOG_COMMIT}`, catalog SHA-256 `{CATALOG_SHA256}`.

## Outcome

**No biological scale-out was authorized or run.** The independent pilot review
is `CONDITIONAL_GO` for bounded repair and ten-slot re-pilot only, explicitly
sets full scale-out to false, and reports zero core pilot passes. The Freeze 1
mirror is a complete metadata inventory but all 47,870 objects remain `planned`;
zero payload objects are verified or reusable. This run therefore submitted
zero Slurm jobs, completed zero pairs, computed zero callable-diversity or PSMC
estimates, and promoted zero pair packets. Technical non-execution is not a
low-diversity result.

## Closed-world accounting

The manifest contains all **{len(rows)} catalog rows** plus all **{len(pairs)}
catalog-linked haplotype entries** (264 labeled other-high-quality and 305
labeled alternate). Catalog-row dispositions are: {dict(sorted(row_counts.items()))}.
Linked-entry dispositions are: {dict(sorted(pair_counts.items()))}. Three linked
entries repeat the H1 accession and are excluded as self-links, not pairs. The
remaining 566 distinct pair candidates are failed operationally at the upstream
authorization gate; they are not declared scientifically ineligible. Each pair
has exactly one manifest row and one QC row ({len(qc)} = {len(pairs)}).

The catalog label “other high-quality” is discovery metadata, not proof that a
pair passes exact-individual provenance, final-sequence QV, completeness,
duplication/collapse, mutual comparability, whole-assembly 1:1, or measured
callability. Those gates remain unmeasured. Hi-C/trio catalog signals are retained
as confidence context; their absence is explicitly not a core veto. Annotation
absence is never a core veto.

## Biological outputs and interpretation

No callable denominator, diversity estimate, mask, consensus, unscaled PSMC,
or bootstrap exists. PSMC requires 200 predeclared boundary-aware replicates
(minimum 190 finite; blocks never cross contigs or mask discontinuities) after a
core pass. The scenario table retains unscaled PSMC as primary and refuses to
invent mutation or generation values. There are no approved confidence tiers,
so there are no scaled trajectories to report.

CDS, fourfold, nonsynonymous, synonymous, WS, SW, and GC3 partitions were not
run. They remain optional post-core outputs requiring an exact native annotation
accession/version and equal sequence dictionary, or a separately validated,
manifest-bound liftover. No valid non-annotated core result was deleted or
downgraded—none yet exists.

Assembly-derived PSMC from an H1/H2 pair is descriptive demographic context.
It reuses the same individual, haplotypes, variants, callable mask, and consensus
as same-pair diversity and **is not statistically independent evidence** for
that diversity.

## Independent validation and sensitivity

Independent biological recomputation, raw-read/k-mer checking, and literature
triangulation were not performed: no pair passed the core gates, the mirror has
no verified sequence payload, and the authorization forbids full processing.
The stratified biological sample is therefore empty rather than post hoc. The
sensitivity table reports QV, completeness, duplication/collapse, long-range
phasing, callability, genome size, generation, bootstrap, and scaling effects as
not estimable. It is invalid to infer any sensitivity from universal technical
non-execution.

## Resumable wave contract

`GATE-0` is the only current state and has a finite zero-job limit. The low/base/
high rows are non-authorized planning templates copied from the reviewed pilot
resource envelope, whose mapping/PSMC prediction error was not estimable. A
future GO requires a new immutable pair audit and wave manifest; the templates
cannot authorize themselves. Each future wave is capped at 25 pairs, uses finite
array concurrency and per-job resources, allows two transient retries and one
resource re-estimation retry, permits zero scientific-threshold relaxation,
checkpoints stage digests, and promotes only a complete output packet by atomic
rename after digest/completeness verification. A job stops at a hard scientific
gate or 1.5× its then-reviewed high resource estimate. These are operational
limits, never a global byte/memory eligibility ceiling. A new assembly release
requires a new manifest and cannot drift into Freeze 1.

## Digest verification

Inputs:

{chr(10).join(f'- `{name}`: `{digest}`' for name, digest in sorted(input_digests.items()))}

Outputs (computed after atomic write):

{chr(10).join(f'- `{name}`: `{digest}`' for name, digest in sorted(output_digests.items()))}
"""


def generate(output: Path = ANALYSIS, catalog_path: Path = CATALOG) -> dict[str, object]:
    catalog_rows = verify_catalog(catalog_path)
    decision = verify_authorization(DECISION, REVIEW_SCALEOUT)
    mirror_summary, mirror_states = verify_mirror(MIRROR_SUMMARY, MIRROR_MANIFEST)
    resources = read_tsv(REVIEW_RESOURCES)
    if {row["scenario"] for row in resources} != {"low", "base", "high"}:
        raise ScaleoutError("reviewed resource scenario set drift")
    input_paths = {
        "catalog": catalog_path, "review_decision": DECISION,
        "review_scaleout": REVIEW_SCALEOUT, "review_resources": REVIEW_RESOURCES,
        "mirror_summary": MIRROR_SUMMARY, "mirror_manifest": MIRROR_MANIFEST,
        "design": DESIGN, "guix_channels": CHANNELS, "guix_manifest": GUIX_MANIFEST,
    }
    input_digests = {name: sha256_file(path) for name, path in input_paths.items()}
    environment_digest = hashlib.sha256(
        f"{input_digests['guix_channels']}:{input_digests['guix_manifest']}".encode()
    ).hexdigest()
    manifest, qc = build_ledgers(
        catalog_rows, mirror_states, input_digests["review_decision"], input_digests["mirror_summary"]
    )
    telemetry = build_telemetry(mirror_summary, resources, input_digests)
    waves = build_waves(resources, environment_digest)
    pair_table, summary_table = paper_tables(manifest)
    sensitivity = sensitivity_rows()
    scaling = scaling_rows()
    independent = independent_validation_rows()

    paths = {
        "manifest": output / "vgp_core_scaleout_manifest.tsv",
        "qc": output / "vgp_core_scaleout_qc.tsv",
        "telemetry": output / "vgp_core_scaleout_telemetry.tsv",
        "waves": output / "vgp_core_scaleout_wave_manifest.tsv",
        "paper_pairs": output / "vgp_core_scaleout_paper_pairs.tsv",
        "paper_summary": output / "vgp_core_scaleout_paper_summary.tsv",
        "sensitivity": output / "vgp_core_scaleout_sensitivity.tsv",
        "scaling": output / "vgp_core_scaleout_scaling_scenarios.tsv",
        "independent_validation": output / "vgp_core_scaleout_independent_validation.tsv",
    }
    atomic_tsv(paths["manifest"], MANIFEST_FIELDS, manifest)
    atomic_tsv(paths["qc"], QC_FIELDS, qc)
    atomic_tsv(paths["telemetry"], TELEMETRY_FIELDS, telemetry)
    atomic_tsv(paths["waves"], WAVE_FIELDS, waves)
    atomic_tsv(paths["paper_pairs"], tuple(pair_table[0]), pair_table)
    atomic_tsv(paths["paper_summary"], tuple(summary_table[0]), summary_table)
    atomic_tsv(paths["sensitivity"], tuple(sensitivity[0]), sensitivity)
    atomic_tsv(paths["scaling"], tuple(scaling[0]), scaling)
    atomic_tsv(paths["independent_validation"], tuple(independent[0]), independent)
    output_digests = {name: sha256_file(path) for name, path in paths.items()}
    generated = utc_now()
    results_path = output / "vgp_core_scaleout_results.md"
    atomic_text(results_path, render_results(manifest, qc, input_digests, output_digests, generated))
    output_digests["results"] = sha256_file(results_path)

    pair_rows = [row for row in manifest if row["record_type"] == "linked_pair"]
    row_rows = [row for row in manifest if row["record_type"] == "catalog_row"]
    summary = {
        "schema_version": "vgp-core-scaleout-summary-v1.0.0", "generated_at_utc": generated,
        "catalog": {"commit": CATALOG_COMMIT, "sha256": CATALOG_SHA256, "rows": len(row_rows)},
        "linked_haplotype_entries": len(pair_rows),
        "distinct_nonself_pair_candidates": sum(row["disposition"] != "excluded" for row in pair_rows),
        "pair_dispositions": dict(sorted(Counter(row["disposition"] for row in pair_rows).items())),
        "catalog_row_dispositions": dict(sorted(Counter(row["disposition"] for row in row_rows).items())),
        "eligible_pairs": 0, "completed_pairs": 0, "biological_estimates": 0,
        "slurm_jobs_submitted": 0, "atomic_promotions": 0,
        "authorization": decision["authorization"], "full_biological_scaleout_authorized": False,
        "mirror_verified_or_reused_objects": 0,
        "annotation_absence_core_veto": False,
        "same_pair_psmc_independent_evidence": False,
        "input_digests": input_digests, "output_digests": output_digests,
    }
    atomic_json(output / "vgp_core_scaleout_summary.json", summary)
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ANALYSIS)
    parser.add_argument("--catalog", type=Path, default=CATALOG)
    args = parser.parse_args()
    try:
        if not os.environ.get("GUIX_ENVIRONMENT"):
            raise ScaleoutError("scale-out generation must run inside pinned GNU Guix")
        summary = generate(args.output, args.catalog)
    except ScaleoutError as error:
        print(f"REFUSED: {error}", file=os.sys.stderr)
        return 65
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
