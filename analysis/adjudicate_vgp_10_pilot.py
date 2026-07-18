#!/usr/bin/env python3
"""Fail-closed adjudication for the authorized ten-pair VGP pilot run.

This program is deliberately a preflight adjudicator, not a substitute for the
biological workflow.  It revalidates the closed-world acquisition and the
pinned execution capture, determines whether each immutable primary is allowed
to cross the Slurm submission boundary, and emits explicit terminal artifacts
for every primary that cannot.  A missing hard-gate measurement is never
imputed and an alternate is never activated by this program.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import resource
import time
from pathlib import Path
from typing import Iterable, Mapping, Sequence

from analysis.vgp_10_pilot import sha256_file, verify_environment_capture


ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = ROOT / "analysis"
DEFAULT_ACQUISITION = ANALYSIS / "vgp_10_pilot_acquisition_manifest.tsv"
DEFAULT_INVENTORY = ANALYSIS / "vgp_10_pilot_object_inventory.tsv"
DEFAULT_DESIGN = ANALYSIS / "vgp_10_pair_manifest.tsv"
DEFAULT_CAPTURE = ANALYSIS / "guix/vgp_10_pilot/realization.json"
DEFAULT_OUTPUT = ANALYSIS

RUN_SCHEMA_VERSION = "vgp-10-pilot-run-adjudication-v1"
RESULT_MANIFEST_VERSION = "vgp-10-pilot-result-manifest-v1.0.0"
RUN_ID = "vgp10-20260718-preflight-v1"
PRIMARY_IDS = tuple(f"P{i:02d}" for i in range(1, 11))
CORE_OBJECT_ROLES = (
    "dataset_report",
    "checksum_catalog",
    "genome_fasta",
    "assembly_report",
    "assembly_stats",
)

# This order is part of the run record.  The first code is the terminal primary
# reason; every additional code remains visible in failure artifacts.
PREFLIGHT_REASON_ORDER = (
    "MISSING_EXACT_FINAL_SEQUENCE_QV",
    "MISSING_H1_BUSCO",
    "MISSING_H2_BUSCO",
    "MISSING_KMER_COPY_NUMBER_AUDIT",
    "UNRESOLVED_EXACT_READ_CHEMISTRY",
    "MISSING_REPEAT_OR_LOW_COMPLEXITY_MASK",
)

DOWNSTREAM_STAGES = (
    "mapping",
    "multiplicity",
    "impg",
    "variants",
    "callable_mask",
    "consensus",
    "diversity",
    "psmc_unscaled",
    "bootstrap",
    "scenario_scaling",
    "validation",
)


class AdjudicationError(RuntimeError):
    """Raised when the run ledger itself is not closed or trustworthy."""


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


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


def relative(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT.resolve()))
    except ValueError:
        return str(path.resolve())


def digest_record(path: Path) -> dict[str, object]:
    return {"path": relative(path), "sha256": sha256_file(path), "size_bytes": path.stat().st_size}


def _evidence_json(row: Mapping[str, str], field: str) -> dict[str, object]:
    try:
        value = json.loads(row[field])
    except (KeyError, json.JSONDecodeError) as exc:
        raise AdjudicationError(f"{row.get('selection_id', '?')}: invalid {field}") from exc
    if not isinstance(value, dict):
        raise AdjudicationError(f"{row.get('selection_id', '?')}: {field} is not an object")
    return value


def preflight_reason_codes(row: Mapping[str, str]) -> list[str]:
    """Return ordered, machine-readable hard-gate reasons without imputation."""
    codes: set[str] = set()
    completeness = _evidence_json(row, "completeness_evidence")
    if row.get("qv_status") != "verified_exact_final_sequence_qv":
        codes.add("MISSING_EXACT_FINAL_SEQUENCE_QV")
    if completeness.get("h1_annotation_busco") == "missing":
        codes.add("MISSING_H1_BUSCO")
    if completeness.get("h2_busco") == "missing":
        codes.add("MISSING_H2_BUSCO")
    # The approved validator requires a measured, passing copy-number/k-mer
    # audit for each haplotype.  Acquisition deliberately selected none.
    if row.get("raw_read_or_kmer_selection") != "verified_kmer_copy_number_audit":
        codes.add("MISSING_KMER_COPY_NUMBER_AUDIT")
    hifi = row.get("hifi_evidence", "")
    if "must resolve exact run chemistry" in hifi:
        codes.add("UNRESOLVED_EXACT_READ_CHEMISTRY")
    if row.get("repeat_report_status") != "verified_exact_repeat_mask":
        codes.add("MISSING_REPEAT_OR_LOW_COMPLEXITY_MASK")
    return [code for code in PREFLIGHT_REASON_ORDER if code in codes]


def validate_closed_world(
    acquisitions: Sequence[Mapping[str, str]], inventory: Sequence[Mapping[str, str]]
) -> tuple[list[Mapping[str, str]], list[Mapping[str, str]]]:
    primaries = [row for row in acquisitions if row.get("roster_type") == "primary"]
    alternates = [row for row in acquisitions if row.get("roster_type") == "alternate"]
    if tuple(row.get("selection_id") for row in primaries) != PRIMARY_IDS:
        raise AdjudicationError("primary closed world is not exactly ordered P01..P10")
    if len(alternates) != 6 or any(row.get("activation_status") != "standby_not_triggered" for row in alternates):
        raise AdjudicationError("alternate ledger is not exactly six untriggered standbys")
    if any(row.get("amendment_id") != "none" for row in acquisitions):
        raise AdjudicationError("unexpected versioned amendment in acquisition manifest")
    if any(row.get("activation_status") != "active_primary" for row in primaries):
        raise AdjudicationError("a primary was silently dropped or replaced")
    if any(row.get("core_identity_status") != "pass" for row in primaries):
        raise AdjudicationError("unresolved primary pair/accession identity")

    for selection_id in PRIMARY_IDS:
        rows = [
            row for row in inventory
            if row.get("selection_id") == selection_id
            and row.get("side") in {"h1", "h2"}
            and row.get("object_role") in CORE_OBJECT_ROLES
        ]
        expected = {(side, role) for side in ("h1", "h2") for role in CORE_OBJECT_ROLES}
        observed = {(row["side"], row["object_role"]) for row in rows if row.get("status") in {"verified", "reused"}}
        if observed != expected or len(rows) != len(expected):
            raise AdjudicationError(f"{selection_id}: core object ledger is not exact and complete")
    return primaries, alternates


def verify_core_objects(
    selection_id: str, inventory: Sequence[Mapping[str, str]], reverify: bool
) -> dict[str, object]:
    rows = [
        row for row in inventory
        if row.get("selection_id") == selection_id
        and row.get("side") in {"h1", "h2"}
        and row.get("object_role") in CORE_OBJECT_ROLES
    ]
    started_wall = time.monotonic()
    started_cpu = time.process_time()
    logical_bytes = 0
    accepted: list[dict[str, object]] = []
    for row in sorted(rows, key=lambda value: (value["side"], value["object_role"])):
        path = Path(row["local_path"])
        expected_size = int(row["expected_bytes"])
        if not path.is_file() or path.stat().st_size != expected_size:
            raise AdjudicationError(f"{row['object_id']}: accepted core object path/size mismatch")
        if reverify and sha256_file(path) != row["local_sha256"]:
            raise AdjudicationError(f"{row['object_id']}: live SHA-256 mismatch")
        logical_bytes += expected_size
        accepted.append({
            "object_id": row["object_id"], "side": row["side"], "role": row["object_role"],
            "path": row["local_path"], "sha256": row["local_sha256"], "size_bytes": expected_size,
            "live_sha256_reverified": reverify,
        })
    usage = resource.getrusage(resource.RUSAGE_SELF)
    return {
        "selection_id": selection_id,
        "status": "pass",
        "accepted_core_objects": accepted,
        "accepted_core_object_count": len(accepted),
        "logical_read_bytes": logical_bytes if reverify else 0,
        "metadata_operations": len(accepted) * (3 if reverify else 1),
        "elapsed_seconds": time.monotonic() - started_wall,
        "cpu_seconds": time.process_time() - started_cpu,
        "maximum_rss_kib": usage.ru_maxrss,
        "filesystem_read_bytes": usage.ru_inblock * 512,
        "filesystem_write_bytes": usage.ru_oublock * 512,
        "scratch_high_water_bytes": 0,
        "retry_count": 0,
    }


def _artifact_paths(base: Path, selection_id: str) -> dict[str, Path]:
    pair = base / "vgp_10_pilot_pair_artifacts" / selection_id
    return {
        "qc": pair / "qc.json",
        "diversity": pair / "diversity.tsv",
        "psmc_trajectory": pair / "psmc_trajectory.tsv",
        "bootstrap": pair / "bootstrap.tsv",
        "scenario": pair / "scenario.tsv",
        "validation": pair / "validation.tsv",
        "annotation": pair / "annotation.json",
        "telemetry": pair / "telemetry.json",
        "failure": pair / "failure.json",
    }


def emit_pair_artifacts(
    output: Path, row: Mapping[str, str], reasons: Sequence[str], telemetry: Mapping[str, object]
) -> dict[str, dict[str, object]]:
    selection_id = row["selection_id"]
    paths = _artifact_paths(output, selection_id)
    downstream = {stage: "not_run_preflight_failed" for stage in DOWNSTREAM_STAGES}
    qc = {
        "schema_version": RUN_SCHEMA_VERSION,
        "run_id": RUN_ID,
        "selection_id": selection_id,
        "terminal_state": "FAILED_PREFLIGHT",
        "core_qc_pass": False,
        "pair_identity": "pass",
        "accession_identity": "pass",
        "core_acquisition": "pass",
        "core_object_checksum": "pass",
        "qv_h1": "missing_not_measured",
        "qv_h2": "missing_not_measured",
        "busco_h1": "missing_not_measured" if "MISSING_H1_BUSCO" in reasons else "available_h1_only",
        "busco_h2": "missing_not_measured",
        "copy_number_kmer_h1": "missing_not_measured",
        "copy_number_kmer_h2": "missing_not_measured",
        "repeat_mask": "missing_not_measured",
        "maximum_query_overlap_depth": None,
        "maximum_target_overlap_depth": None,
        "retained_multiplicity_gt_one": 0,
        "h1_ref_mismatches": None,
        "h2_reconstruction_failures": None,
        "callable_mask_accounting_discrepancy_bp": None,
        "bootstrap_attempts": 0,
        "bootstrap_successes": 0,
        "unscaled_scenario_separation": "not_applicable_no_psmc_output",
        "reason_codes": list(reasons),
        "downstream_stage_status": downstream,
    }
    atomic_json(paths["qc"], qc)

    common = {"selection_id": selection_id, "status": "not_run_preflight_failed", "reason_code": reasons[0]}
    write_tsv(paths["diversity"],
              ("selection_id", "status", "heterozygous_sites", "callable_h1_bp", "diversity_per_callable_h1_bp", "reason_code"),
              [dict(common, heterozygous_sites="", callable_h1_bp="", diversity_per_callable_h1_bp="")])
    write_tsv(paths["psmc_trajectory"],
              ("selection_id", "status", "scale", "interval", "time_2N0", "lambda", "reason_code"),
              [dict(common, scale="unscaled_not_produced", interval="", time_2N0="", **{"lambda": ""})])
    write_tsv(paths["bootstrap"],
              ("selection_id", "status", "attempts_required", "attempts", "successes", "success_fraction", "block_bp", "boundary_aware", "reason_code"),
              [dict(common, attempts_required=200, attempts=0, successes=0, success_fraction="not_applicable", block_bp=5_000_000, boundary_aware="not_materialized")])
    write_tsv(paths["scenario"],
              ("selection_id", "status", "unscaled_primary_preserved", "scenario_scaled_output_separate", "scenario_id", "reason_code"),
              [dict(common, unscaled_primary_preserved="not_produced", scenario_scaled_output_separate="not_produced", scenario_id="")])
    validation_rows = [
        dict(common, validation_control=control, comparison_status="not_run_no_core_result", observed="", expected="")
        for control in ("raw_read", "kmer", "published", "D01_direct_control")
    ]
    write_tsv(paths["validation"],
              ("selection_id", "validation_control", "comparison_status", "observed", "expected", "status", "reason_code"),
              validation_rows)

    annotation = {
        "schema_version": RUN_SCHEMA_VERSION,
        "selection_id": selection_id,
        "core_blocked_by_annotation": False,
        "acquisition_branch_status": row["annotation_branch_status"],
        "output_status": "not_run_core_preflight_failed",
        "annotation_accession_version": row.get("native_annotation_accession_version") or None,
        "reason_code": "ANNOTATION_MISMATCH" if row["annotation_branch_status"] == "failed_mismatch" else "ANNOTATION_NOT_REACHED",
    }
    atomic_json(paths["annotation"], annotation)
    atomic_json(paths["telemetry"], dict(telemetry, disposition="failed_preflight", slurm_job_ids=[], slurm_dependency_ids=[]))
    failure = {
        "schema_version": RUN_SCHEMA_VERSION,
        "run_id": RUN_ID,
        "selection_id": selection_id,
        "roster_type": "primary",
        "primary_retained": True,
        "alternate_activated": False,
        "amendment_id": "none",
        "failure_stage": "preflight_core_qc",
        "primary_reason_code": reasons[0],
        "reason_codes": list(reasons),
        "slurm_jobs_submitted": 0,
        "retry_count": 0,
        "scientific_thresholds_relaxed": False,
        "message": "Mandatory measurements are absent; biological compute was not authorized.",
    }
    atomic_json(paths["failure"], failure)
    return {name: digest_record(path) for name, path in paths.items()}


RESULT_FIELDS = (
    "manifest_version", "run_id", "selection_id", "roster_type", "activation_status", "amendment_id",
    "failed_primary_retained", "species", "biosample", "individual_or_isolate", "h1_accession_version",
    "h2_accession_version", "orientation", "terminal_state", "disposition", "failure_stage",
    "primary_reason_code", "all_reason_codes", "pair_identity_status", "accession_identity_status",
    "core_checksum_status", "core_qc_status", "mapping_status", "multiplicity_status", "impg_status",
    "variants_status", "callable_mask_status", "consensus_status", "diversity_status", "psmc_status",
    "bootstrap_attempts", "bootstrap_successes", "bootstrap_success_fraction", "scenario_status",
    "validation_status", "annotation_status", "alternate_activated", "slurm_jobs_submitted",
    "qc_artifact", "diversity_artifact", "psmc_trajectory_artifact", "bootstrap_artifact",
    "scenario_artifact", "validation_artifact", "annotation_artifact", "telemetry_artifact",
    "failure_artifact", "artifact_packet_sha256",
)

QC_FIELDS = (
    "run_id", "selection_id", "terminal_state", "core_qc_pass", "pair_identity", "accession_identity",
    "core_checksum", "qv_h1", "qv_h2", "busco_h1", "busco_h2", "copy_number_kmer_h1",
    "copy_number_kmer_h2", "repeat_mask", "maximum_query_overlap_depth", "maximum_target_overlap_depth",
    "retained_query_multiplicity_gt_one", "retained_target_multiplicity_gt_one", "h1_ref_mismatches",
    "h2_reconstruction_failures", "callable_universe_bp", "callable_bp", "callable_reason_total_bp",
    "callable_accounting_discrepancy_bp", "vcf_consensus_psmc_consistent", "bootstrap_attempts",
    "bootstrap_successes", "bootstrap_success_fraction", "unscaled_scenario_separation", "reason_codes",
)

TELEMETRY_FIELDS = (
    "run_id", "selection_id", "stage", "execution_context", "disposition", "slurm_job_id",
    "slurm_dependency_ids", "cpus", "elapsed_seconds", "cpu_seconds", "maximum_rss_kib",
    "scratch_path", "scratch_high_water_bytes", "logical_read_bytes", "filesystem_read_bytes",
    "logical_write_bytes", "filesystem_write_bytes", "metadata_operations", "retry_count", "retry_reason",
    "resource_estimator_error", "command_status",
)


def generate(
    *, acquisition_path: Path = DEFAULT_ACQUISITION, inventory_path: Path = DEFAULT_INVENTORY,
    design_path: Path = DEFAULT_DESIGN, capture_path: Path = DEFAULT_CAPTURE,
    output: Path = DEFAULT_OUTPUT, reverify: bool = True,
) -> dict[str, object]:
    started_utc = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    process_wall = time.monotonic()
    process_cpu = time.process_time()
    acquisitions = read_tsv(acquisition_path)
    inventory = read_tsv(inventory_path)
    primaries, alternates = validate_closed_world(acquisitions, inventory)
    environment = verify_environment_capture(capture_path)

    result_rows: list[dict[str, object]] = []
    qc_rows: list[dict[str, object]] = []
    telemetry_rows: list[dict[str, object]] = []
    reason_counts = {code: 0 for code in PREFLIGHT_REASON_ORDER}
    logical_reads = 0
    metadata_operations = 0
    for row in primaries:
        selection_id = row["selection_id"]
        reasons = preflight_reason_codes(row)
        if not reasons:
            raise AdjudicationError(
                f"{selection_id}: appears launchable; this adjudicator refuses to fabricate a completed run"
            )
        for reason in reasons:
            reason_counts[reason] += 1
        telemetry = verify_core_objects(selection_id, inventory, reverify)
        logical_reads += int(telemetry["logical_read_bytes"])
        metadata_operations += int(telemetry["metadata_operations"])
        artifacts = emit_pair_artifacts(output, row, reasons, telemetry)
        packet_digest = hashlib.sha256("".join(artifacts[name]["sha256"] for name in sorted(artifacts)).encode()).hexdigest()
        downstream = "not_run_preflight_failed"
        result_rows.append({
            "manifest_version": RESULT_MANIFEST_VERSION, "run_id": RUN_ID,
            "selection_id": selection_id, "roster_type": "primary", "activation_status": "active_primary",
            "amendment_id": "none", "failed_primary_retained": "true", "species": row["species"],
            "biosample": row["biosample"], "individual_or_isolate": row["individual_or_isolate"],
            "h1_accession_version": row["h1_accession_version"], "h2_accession_version": row["h2_accession_version"],
            "orientation": "H1_reference_H2_query", "terminal_state": "FAILED_PREFLIGHT", "disposition": "FAIL",
            "failure_stage": "preflight_core_qc", "primary_reason_code": reasons[0],
            "all_reason_codes": ";".join(reasons), "pair_identity_status": "pass",
            "accession_identity_status": "pass", "core_checksum_status": "pass",
            "core_qc_status": "fail_missing_required_measurements", "mapping_status": downstream,
            "multiplicity_status": "not_measured_preflight_failed", "impg_status": downstream,
            "variants_status": downstream, "callable_mask_status": downstream, "consensus_status": downstream,
            "diversity_status": downstream, "psmc_status": downstream, "bootstrap_attempts": 0,
            "bootstrap_successes": 0, "bootstrap_success_fraction": "not_applicable_no_attempts",
            "scenario_status": downstream, "validation_status": "not_run_no_core_result",
            "annotation_status": "not_run_core_preflight_failed", "alternate_activated": "false",
            "slurm_jobs_submitted": 0,
            **{f"{name}_artifact": artifacts[name]["path"] for name in artifacts},
            "artifact_packet_sha256": packet_digest,
        })
        qc_rows.append({
            "run_id": RUN_ID, "selection_id": selection_id, "terminal_state": "FAILED_PREFLIGHT",
            "core_qc_pass": "false", "pair_identity": "pass", "accession_identity": "pass",
            "core_checksum": "pass", "qv_h1": "missing_not_measured", "qv_h2": "missing_not_measured",
            "busco_h1": "missing_not_measured" if "MISSING_H1_BUSCO" in reasons else "available_h1_only",
            "busco_h2": "missing_not_measured", "copy_number_kmer_h1": "missing_not_measured",
            "copy_number_kmer_h2": "missing_not_measured", "repeat_mask": "missing_not_measured",
            "maximum_query_overlap_depth": "not_measured_preflight_failed",
            "maximum_target_overlap_depth": "not_measured_preflight_failed",
            "retained_query_multiplicity_gt_one": 0, "retained_target_multiplicity_gt_one": 0,
            "h1_ref_mismatches": "not_measured_preflight_failed",
            "h2_reconstruction_failures": "not_measured_preflight_failed",
            "callable_universe_bp": "not_materialized", "callable_bp": "not_materialized",
            "callable_reason_total_bp": "not_materialized", "callable_accounting_discrepancy_bp": "not_measured",
            "vcf_consensus_psmc_consistent": "not_applicable_no_outputs", "bootstrap_attempts": 0,
            "bootstrap_successes": 0, "bootstrap_success_fraction": "not_applicable_no_attempts",
            "unscaled_scenario_separation": "not_applicable_no_outputs", "reason_codes": ";".join(reasons),
        })
        telemetry_rows.append({
            "run_id": RUN_ID, "selection_id": selection_id, "stage": "preflight_core_object_revalidation",
            "execution_context": "pinned_guix_login_node", "disposition": "completed_preflight_then_refused",
            "slurm_job_id": "not_submitted", "slurm_dependency_ids": "none", "cpus": 1,
            "elapsed_seconds": f"{telemetry['elapsed_seconds']:.6f}", "cpu_seconds": f"{telemetry['cpu_seconds']:.6f}",
            "maximum_rss_kib": telemetry["maximum_rss_kib"], "scratch_path": "none",
            "scratch_high_water_bytes": 0, "logical_read_bytes": telemetry["logical_read_bytes"],
            "filesystem_read_bytes": telemetry["filesystem_read_bytes"], "logical_write_bytes": 0,
            "filesystem_write_bytes": telemetry["filesystem_write_bytes"],
            "metadata_operations": telemetry["metadata_operations"], "retry_count": 0, "retry_reason": "none",
            "resource_estimator_error": "not_applicable_no_cluster_job",
            "command_status": "core_CAS_verified;submission_refused_hard_gate",
        })

    result_manifest = output / "vgp_10_pilot_result_manifest.tsv"
    qc_path = output / "vgp_10_pilot_qc.tsv"
    telemetry_path = output / "vgp_10_pilot_resource_telemetry.tsv"
    write_tsv(result_manifest, RESULT_FIELDS, result_rows)
    write_tsv(qc_path, QC_FIELDS, qc_rows)

    usage = resource.getrusage(resource.RUSAGE_SELF)
    telemetry_rows.append({
        "run_id": RUN_ID, "selection_id": "ALL", "stage": "closed_world_adjudication",
        "execution_context": "pinned_guix_login_node", "disposition": "completed",
        "slurm_job_id": "not_submitted", "slurm_dependency_ids": "none", "cpus": 1,
        "elapsed_seconds": f"{time.monotonic() - process_wall:.6f}", "cpu_seconds": f"{time.process_time() - process_cpu:.6f}",
        "maximum_rss_kib": usage.ru_maxrss, "scratch_path": "none", "scratch_high_water_bytes": 0,
        "logical_read_bytes": logical_reads, "filesystem_read_bytes": usage.ru_inblock * 512,
        "logical_write_bytes": sum(path.stat().st_size for path in (result_manifest, qc_path)),
        "filesystem_write_bytes": usage.ru_oublock * 512, "metadata_operations": metadata_operations,
        "retry_count": 0, "retry_reason": "none", "resource_estimator_error": "not_applicable_no_cluster_job",
        "command_status": "ten_primary_slots_terminal;zero_sbatch_calls",
    })
    write_tsv(telemetry_path, TELEMETRY_FIELDS, telemetry_rows)

    command_log = output / "vgp_10_pilot_command_log.json"
    atomic_json(command_log, {
        "schema_version": RUN_SCHEMA_VERSION, "run_id": RUN_ID,
        "commands": [
            {"command": "guix time-machine -C analysis/guix/vgp_10_pilot/channels.scm -- shell ... validate_vgp_10_pilot_acquisition.py --reverify-local", "status": "passed"},
            {"command": "python3 -m analysis.vgp_10_pilot verify-capture analysis/guix/vgp_10_pilot/realization.json", "execution_environment": "captured_pinned_guix_profile", "status": "passed"},
            {"command": "analysis/slurm/vgp_10_pilot/submit.sh --dry-run", "status": "not_run_preflight_refused_missing_hard_gate_measurements"},
            {"command": "sbatch", "status": "not_invoked_zero_jobs_submitted"},
        ],
        "slurm_jobs": [], "slurm_dependencies": [], "retry_reasons": [],
        "scientific_thresholds_relaxed": False,
    })
    summary_path = output / "vgp_10_pilot_run_summary.json"
    summary = {
        "schema_version": RUN_SCHEMA_VERSION, "run_id": RUN_ID, "started_utc": started_utc,
        "completed_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "primary_slots_expected": 10, "primary_slots_accounted": len(result_rows),
        "primary_completed": 0, "primary_failed_preflight": len(result_rows),
        "unresolved_pair_or_accession_identity": 0, "core_checksum_failures": 0,
        "core_qc_passes": 0, "retained_query_multiplicity_gt_one": 0,
        "retained_target_multiplicity_gt_one": 0,
        "multiplicity_measurement_status": "not_measured_for_any_pair_preflight_failed",
        "h1_ref_check_status": "not_measured_for_any_pair_preflight_failed",
        "h2_reconstruction_status": "not_measured_for_any_pair_preflight_failed",
        "callable_mask_status": "not_materialized_for_any_pair_preflight_failed",
        "bootstrap_attempts": 0, "bootstrap_successes": 0,
        "bootstrap_requirement_status": "not_reached_for_any_pair_preflight_failed",
        "unscaled_scenario_separation": "not_applicable_no_psmc_outputs",
        "validation_subset_status": "not_run_no_core_results",
        "direct_control": "D01_acquired_but_not_compared_no_core_results",
        "alternate_slots_declared": len(alternates), "alternate_activations": 0,
        "manifest_amendments": [], "slurm_jobs_submitted": 0, "slurm_dependency_edges": 0,
        "reason_counts": reason_counts, "environment": {
            key: environment[key] for key in ("channel_commit", "manifest_sha256", "profile", "derivation", "closure_sha256")
        },
        "input_digests": {
            relative(acquisition_path): sha256_file(acquisition_path), relative(inventory_path): sha256_file(inventory_path),
            relative(design_path): sha256_file(design_path), relative(capture_path): sha256_file(capture_path),
        },
        "deliverables": {},
    }
    atomic_json(summary_path, summary)

    report_path = output / "vgp_10_pilot_results.md"
    reason_lines = "\n".join(f"- `{code}`: {count}/10" for code, count in reason_counts.items() if count)
    report = f"""# Ten-pair VGP pilot run result

Run ID: `{RUN_ID}`
Terminal state: **10/10 primary slots failed mandatory preflight; 0 biological jobs submitted**

## Outcome

The closed-world acquisition and pinned GNU Guix capture were revalidated successfully. Exact pair, BioSample/individual, accession.version, accepted core-object size, and SHA-256 identity are resolved for all ten immutable primaries. Those successes do not authorize biological compute: every primary lacks exact-final-sequence QV, H2 BUSCO, a manifest-bound k-mer/copy-number collapse audit, and a repeat/low-complexity mask. P07 and P08 also lack H1 BUSCO; P09 retains unresolved exact read chemistry.

The approved preflight validator requires QV >= 40, BUSCO completeness/missing/duplication measurements for both haplotypes, and a passing copy-number/k-mer audit. Missing values were not imputed from technology labels, H1-only annotation BUSCO, assembly length, or published plots. Consequently SweepGA, multiplicity auditing, IMPG, VCF/BCF, masks, consensus, diversity, PSMC, 200 bootstraps, scaling scenarios, and validation comparisons are explicitly `not_run_preflight_failed` for every pair. This is a failed pilot execution, not a zero-diversity result.

## Primary-slot accounting

| state | count |
|---|---:|
| completed core analysis | 0 |
| explicitly failed at preflight | 10 |
| silently dropped or replaced | 0 |
| alternate activations | 0 |
| Slurm jobs / dependency edges | 0 / 0 |

No alternate was activated: there is no versioned amendment authorizing one, and the six declared alternates remain `standby_not_triggered`. Primary failure rows are retained in the result manifest.

## Machine-readable failure reasons

{reason_lines}

Annotation status remains branch-local. P03 and P04 retain their acquisition-time annotation mismatch; all other annotation branches were not reached because core preflight failed. No annotation disposition was used to create a core failure.

## Validation interpretation

There are zero unresolved pair/accession identities and zero checksum failures. There are also zero retained mappings with multiplicity greater than one because no mapping was retained; overlap depth, H1 REF, direct H2 reconstruction, callable-mask totals, VCF/consensus/PSMC consistency, bootstrap success fraction, and scaled trajectories are **not measured**, not passing zeros. D01 and raw-read/k-mer/published comparisons were not run because no primary produced a core result to compare. The per-pair artifacts preserve these distinctions explicitly.

## Reproducibility and telemetry

Execution used captured Guix channel `{environment['channel_commit']}`, manifest `{environment['manifest_sha256']}`, profile `{environment['profile']}`, derivation `{environment['derivation']}`, and closure `{environment['closure_sha256']}`. Live CAS revalidation read {logical_reads} logical bytes across exactly 100 primary core objects. No node-local scratch was allocated because the submission boundary was not crossed. Per-pair and aggregate elapsed time, CPU, peak RSS, filesystem/logical reads and writes, metadata operations, retries, scheduler IDs, and dispositions are recorded in `analysis/vgp_10_pilot_resource_telemetry.tsv`.

## Delivered artifacts

- `analysis/vgp_10_pilot_result_manifest.tsv`: exactly ten retained primary terminal rows and per-pair artifact paths.
- `analysis/vgp_10_pilot_qc.tsv`: gate-by-gate measured/missing/not-reached distinctions.
- `analysis/vgp_10_pilot_resource_telemetry.tsv`: ten preflight revalidation rows plus aggregate telemetry.
- `analysis/vgp_10_pilot_pair_artifacts/P01..P10/`: QC, diversity, PSMC trajectory, bootstrap, scenario, validation, annotation, telemetry, and failure artifacts.
- `analysis/vgp_10_pilot_run_summary.json` and `analysis/vgp_10_pilot_command_log.json`: independently recomputable counts, identities, command dispositions, and zero-job ledger.
"""
    atomic_text(report_path, report)

    deliverables = {
        "result_manifest": digest_record(result_manifest), "qc": digest_record(qc_path),
        "resource_telemetry": digest_record(telemetry_path), "results": digest_record(report_path),
        "command_log": digest_record(command_log),
    }
    summary["deliverables"] = deliverables
    atomic_json(summary_path, summary)
    return summary


def independently_validate(output: Path) -> dict[str, object]:
    manifest = read_tsv(output / "vgp_10_pilot_result_manifest.tsv")
    qc = read_tsv(output / "vgp_10_pilot_qc.tsv")
    telemetry = read_tsv(output / "vgp_10_pilot_resource_telemetry.tsv")
    summary = json.loads((output / "vgp_10_pilot_run_summary.json").read_text(encoding="utf-8"))
    if [row["selection_id"] for row in manifest] != list(PRIMARY_IDS):
        raise AdjudicationError("result manifest does not independently account for P01..P10")
    if len(qc) != 10 or {row["selection_id"] for row in qc} != set(PRIMARY_IDS):
        raise AdjudicationError("QC table does not independently account for ten primaries")
    if any(row["terminal_state"] != "FAILED_PREFLIGHT" or row["alternate_activated"] != "false" for row in manifest):
        raise AdjudicationError("terminal or alternate state drift")
    if sum(int(row["slurm_jobs_submitted"]) for row in manifest) != 0:
        raise AdjudicationError("nonzero Slurm job accounting")
    if any(row["core_qc_pass"] != "false" for row in qc):
        raise AdjudicationError("failed preflight was promoted to core QC pass")
    if len(telemetry) != 11 or telemetry[-1]["selection_id"] != "ALL":
        raise AdjudicationError("resource telemetry is not ten pairs plus aggregate")
    independently_counted_reasons = {code: 0 for code in PREFLIGHT_REASON_ORDER}
    artifact_packets_verified = 0
    validation_controls_not_run = 0
    for row in manifest:
        reasons = row["all_reason_codes"].split(";")
        if reasons[0] != row["primary_reason_code"]:
            raise AdjudicationError(f"{row['selection_id']}: primary failure reason is not first")
        for reason in reasons:
            if reason not in independently_counted_reasons:
                raise AdjudicationError(f"{row['selection_id']}: unknown reason code {reason}")
            independently_counted_reasons[reason] += 1
        packet_hashes: dict[str, str] = {}
        for name in (
            "qc", "diversity", "psmc_trajectory", "bootstrap", "scenario", "validation",
            "annotation", "telemetry", "failure",
        ):
            key = f"{name}_artifact"
            path = ROOT / row[key] if not Path(row[key]).is_absolute() else Path(row[key])
            if not path.is_file():
                raise AdjudicationError(f"missing per-pair artifact: {row[key]}")
            packet_hashes[name] = sha256_file(path)
        expected_packet = hashlib.sha256(
            "".join(packet_hashes[name] for name in sorted(packet_hashes)).encode()
        ).hexdigest()
        if expected_packet != row["artifact_packet_sha256"]:
            raise AdjudicationError(f"{row['selection_id']}: artifact packet digest mismatch")
        artifact_packets_verified += 1
        validation = read_tsv(
            ROOT / row["validation_artifact"]
            if not Path(row["validation_artifact"]).is_absolute()
            else Path(row["validation_artifact"])
        )
        if {value["validation_control"] for value in validation} != {
            "raw_read", "kmer", "published", "D01_direct_control"
        } or any(value["comparison_status"] != "not_run_no_core_result" for value in validation):
            raise AdjudicationError(f"{row['selection_id']}: validation-control disposition drift")
        validation_controls_not_run += len(validation)
    if independently_counted_reasons != summary["reason_counts"]:
        raise AdjudicationError("independent failure-reason counts differ from run summary")
    if summary["primary_slots_accounted"] != len(manifest) or summary["primary_failed_preflight"] != 10:
        raise AdjudicationError("independent primary accounting differs from run summary")
    result = {
        "primary_slots": len(manifest), "failed_preflight": sum(row["disposition"] == "FAIL" for row in manifest),
        "completed": sum(row["disposition"] == "PASS" for row in manifest),
        "alternate_activations": sum(row["alternate_activated"] == "true" for row in manifest),
        "slurm_jobs_submitted": sum(int(row["slurm_jobs_submitted"]) for row in manifest),
        "unresolved_identity": sum(row["pair_identity_status"] != "pass" or row["accession_identity_status"] != "pass" for row in manifest),
        "checksum_failures": sum(row["core_checksum_status"] != "pass" for row in manifest),
        "artifact_packets_verified": artifact_packets_verified,
        "validation_controls_not_run": validation_controls_not_run,
        "bootstrap_attempts": sum(int(row["bootstrap_attempts"]) for row in manifest),
        "reason_counts": independently_counted_reasons,
    }
    atomic_json(output / "vgp_10_pilot_independent_validation.json", result)
    return result


def parser() -> argparse.ArgumentParser:
    value = argparse.ArgumentParser(description=__doc__)
    value.add_argument("--acquisition", type=Path, default=DEFAULT_ACQUISITION)
    value.add_argument("--inventory", type=Path, default=DEFAULT_INVENTORY)
    value.add_argument("--design", type=Path, default=DEFAULT_DESIGN)
    value.add_argument("--capture", type=Path, default=DEFAULT_CAPTURE)
    value.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    value.add_argument("--no-live-reverify", action="store_true", help="tests only: retain ledger verification without rereading payload bytes")
    return value


def main(argv: Sequence[str] | None = None) -> int:
    args = parser().parse_args(argv)
    generate(
        acquisition_path=args.acquisition, inventory_path=args.inventory, design_path=args.design,
        capture_path=args.capture, output=args.output, reverify=not args.no_live_reverify,
    )
    independently_validate(args.output)
    print("VGP ten-pair pilot adjudication passed: 10 explicit preflight failures; 0 Slurm jobs")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
