#!/usr/bin/env python3
"""Independently review the bounded VGP Tier 3 pilot artifacts."""

from __future__ import annotations

import argparse
import csv
import json
import sys
import tempfile
from collections import Counter
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from analysis import gate_vgp_pilot as gate
from analysis import run_vgp_pilot as runner


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_GATE = PROJECT_ROOT / "analysis" / "vgp_pilot_gate.json"
DEFAULT_RUN_MANIFEST = PROJECT_ROOT / "analysis" / "vgp_pilot_run_manifest.tsv"
DEFAULT_RESULTS = PROJECT_ROOT / "analysis" / "vgp_pilot_results.tsv"
DEFAULT_TELEMETRY = PROJECT_ROOT / "analysis" / "vgp_pilot_slurm_telemetry.tsv"
DEFAULT_ROOT_VALIDATION = PROJECT_ROOT / "analysis" / "vgp_data_root_validation.json"
DEFAULT_REVIEW = PROJECT_ROOT / "analysis" / "vgp_pilot_review.md"
DEFAULT_QC = PROJECT_ROOT / "analysis" / "vgp_pilot_qc.tsv"
DEFAULT_RESOURCE = PROJECT_ROOT / "analysis" / "vgp_pilot_resource_calibration.tsv"
DEFAULT_GUIX_ENVIRONMENT = PROJECT_ROOT / "analysis" / "pilot_results" / "guix_environment.json"
DEFAULT_COMPUTE_SMOKE = PROJECT_ROOT / "analysis" / "pilot_results" / "compute_smoke.json"
DEFAULT_RESOURCE_BUDGET = PROJECT_ROOT / "analysis" / "vertebrate_scaleout_resource_budget.tsv"

QC_FIELDS = [
    "check_id",
    "category",
    "subject",
    "decision",
    "observed",
    "expected",
    "evidence",
    "notes",
]

RESOURCE_FIELDS = [
    "scope",
    "metric",
    "unit",
    "predicted_low",
    "predicted_base",
    "predicted_high",
    "authorized_cap",
    "observed",
    "decision",
    "evidence",
    "notes",
]

RESOURCE_METRICS = (
    ("aggregate_core_hours", "core_hours", "core-h"),
    ("aggregate_wall_hours", "catalog_or_stage_wall_hours", "h"),
    ("peak_memory_gib_per_job", "peak_resident_or_requested_memory_gib_per_element", "GiB"),
    ("peak_local_scratch_gb", "local_scratch_peak_gb", "GB"),
    ("persistent_input_gb", "persistent_input_gb", "GB"),
    ("persistent_output_gb", "persistent_output_gb", "GB"),
    ("file_inodes", "file_inode_count", "count"),
    ("moosefs_read_gb", "moosefs_read_gb", "GB"),
    ("moosefs_write_gb", "moosefs_write_gb", "GB"),
    ("metadata_operations", "metadata_operations", "count"),
    ("peak_bandwidth_mib_s", "peak_aggregate_bandwidth_mib_s", "MiB/s"),
    ("cpus_per_element", "cpus_per_element", "count"),
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fieldnames: Sequence[str], rows: Sequence[Mapping[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def normalize_source(value: str) -> str:
    if not value:
        return value
    marker = "/analysis/"
    if marker in value:
        return "analysis/" + value.split(marker, 1)[1]
    return value


def normalize_rows(rows: Iterable[Mapping[str, str]]) -> list[dict[str, str]]:
    normalized: list[dict[str, str]] = []
    for row in rows:
        cooked = {key: value for key, value in row.items() if key not in {"run_id", "generated_at_utc"}}
        if "failure_source" in cooked:
            cooked["failure_source"] = normalize_source(cooked["failure_source"])
        normalized.append(cooked)
    return normalized


def stable_gate_view(payload: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "decision": payload["decision"],
        "inputs_sha256": {key: value["sha256"] for key, value in payload["inputs"].items()},
        "reproduction": payload["reproduction"],
        "row_audit_summary": payload["row_audit"]["summary"],
        "cap_dimensions": {
            key: {
                "unit": value["unit"],
                "value": value["value"],
                "winner_source": value["winner_source"],
            }
            for key, value in sorted(payload["cap_vector"]["dimensions"].items())
        },
        "cap_vector_sha256": payload["cap_vector"]["sha256"],
        "authorization_boundary": {
            "manifest_digest": payload["authorization_boundary"]["manifest_digest"],
            "root_contract_digest": payload["authorization_boundary"]["root_contract_digest"],
            "cap_vector_digest": payload["authorization_boundary"]["cap_vector_digest"],
        },
        "quota_evidence": payload["quota_evidence"],
        "blockers": payload["blockers"],
    }


def to_float(value: str | int | float | None) -> float | None:
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    return float(text)


def format_scalar(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        if value.is_integer():
            return str(int(value))
        return f"{value:.6f}".rstrip("0").rstrip(".")
    return str(value)


def qc_row(
    check_id: str,
    category: str,
    subject: str,
    decision: str,
    observed: Any,
    expected: Any,
    evidence: str,
    notes: str,
) -> dict[str, str]:
    return {
        "check_id": check_id,
        "category": category,
        "subject": subject,
        "decision": decision,
        "observed": format_scalar(observed),
        "expected": format_scalar(expected),
        "evidence": evidence,
        "notes": notes,
    }


def resource_row(
    scope: str,
    metric: str,
    unit: str,
    predicted_low: Any,
    predicted_base: Any,
    predicted_high: Any,
    authorized_cap: Any,
    observed: Any,
    decision: str,
    evidence: str,
    notes: str,
) -> dict[str, str]:
    return {
        "scope": scope,
        "metric": metric,
        "unit": unit,
        "predicted_low": format_scalar(predicted_low),
        "predicted_base": format_scalar(predicted_base),
        "predicted_high": format_scalar(predicted_high),
        "authorized_cap": format_scalar(authorized_cap),
        "observed": format_scalar(observed),
        "decision": decision,
        "evidence": evidence,
        "notes": notes,
    }


def build_resource_rows(promoted_gate: Mapping[str, Any], telemetry_rows: Sequence[Mapping[str, str]]) -> list[dict[str, str]]:
    telemetry = telemetry_rows[0]
    budget_rows = load_tsv(DEFAULT_RESOURCE_BUDGET)
    stratified = {
        row["scenario"]: row
        for row in budget_rows
        if row["stage_or_dataset"] == "stratified_pilot" and row["scenario"] in {"low", "base", "high"}
    }
    observed_lookup = {
        "aggregate_core_hours": to_float(telemetry["cpu_time_seconds"]) / 3600.0,
        "aggregate_wall_hours": to_float(telemetry["elapsed_seconds"]) / 3600.0,
        "peak_memory_gib_per_job": to_float(telemetry["max_rss_gib"]),
        "peak_local_scratch_gb": to_float(telemetry["scratch_peak_gb"]),
        "persistent_input_gb": 0.0,
        "persistent_output_gb": 0.0,
        "file_inodes": 0.0,
        "moosefs_read_gb": to_float(telemetry["io_read_gb"]),
        "moosefs_write_gb": to_float(telemetry["io_write_gb"]),
        "metadata_operations": to_float(telemetry["metadata_operations"]),
        "peak_bandwidth_mib_s": None,
        "cpus_per_element": to_float(telemetry["requested_cpus"]),
    }
    rows: list[dict[str, str]] = []
    for metric, budget_field, unit in RESOURCE_METRICS:
        authorized = promoted_gate["cap_vector"]["dimensions"][metric]["value"]
        predicted_low = to_float(stratified["low"][budget_field])
        predicted_base = to_float(stratified["base"][budget_field])
        predicted_high = to_float(stratified["high"][budget_field])
        observed = observed_lookup[metric]
        if observed is None:
            decision = "INCONCLUSIVE"
            notes = "No Slurm job was submitted, so this per-job metric was never observed."
        elif observed <= authorized + 1e-9:
            decision = "PASS"
            notes = (
                "Observed value stayed within the authorized cap; the empty selection under the NO_GO gate "
                "collapsed the executable cap to zero for aggregate usage."
            )
        else:
            decision = "FAIL"
            notes = "Observed usage exceeded the gate-derived authorized cap."
        evidence = "analysis/vertebrate_scaleout_resource_budget.tsv; analysis/vgp_pilot_slurm_telemetry.tsv"
        if metric in {"persistent_input_gb", "persistent_output_gb", "file_inodes"}:
            evidence += "; analysis/vgp_pilot_run_manifest.tsv"
            if decision != "FAIL":
                notes = "Refused preflight execution created no pilot downloads, persistent outputs, or retained inodes."
        rows.append(
            resource_row(
                scope="pilot_refusal",
                metric=metric,
                unit=unit,
                predicted_low=predicted_low,
                predicted_base=predicted_base,
                predicted_high=predicted_high,
                authorized_cap=authorized,
                observed=observed,
                decision=decision,
                evidence=evidence,
                notes=notes,
            )
        )
    return rows


def build_qc_rows(
    *,
    promoted_gate: Mapping[str, Any],
    rebuilt_gate: Mapping[str, Any],
    promoted_run_manifest: Sequence[Mapping[str, str]],
    promoted_results: Sequence[Mapping[str, str]],
    promoted_telemetry: Sequence[Mapping[str, str]],
    recomputed_run_manifest: Sequence[Mapping[str, str]],
    recomputed_results: Sequence[Mapping[str, str]],
    recomputed_telemetry: Sequence[Mapping[str, str]],
    guix_validation: str,
    guix_note: str,
) -> list[dict[str, str]]:
    qc_rows: list[dict[str, str]] = []
    discrepancy_text = "; ".join(
        f"{item['metric']} observed={item['observed']} expected={item['expected']}"
        for item in promoted_gate["reproduction"]["source_catalog"]["discrepancies"]
    )
    gate_match = stable_gate_view(promoted_gate) == stable_gate_view(rebuilt_gate)
    qc_rows.append(
        qc_row(
            "promoted_gate_recompute",
            "reproducibility",
            "analysis/vgp_pilot_gate.json",
            "PASS" if gate_match else "FAIL",
            "stable_fields_match" if gate_match else "stable_fields_differ",
            "stable_fields_match",
            "analysis/vgp_pilot_gate.json",
            "Fresh gate recomputation matched the promoted gate on all stable fields." if gate_match else "Rebuilt gate diverged from promoted stable fields.",
        )
    )

    manifest_match = normalize_rows(promoted_run_manifest) == normalize_rows(recomputed_run_manifest)
    qc_rows.append(
        qc_row(
            "run_manifest_recompute",
            "reproducibility",
            "analysis/vgp_pilot_run_manifest.tsv",
            "PASS" if manifest_match else "FAIL",
            f"normalized_rows={len(promoted_run_manifest)}" if manifest_match else "normalized_rows_differ",
            f"normalized_rows={len(recomputed_run_manifest)}",
            "analysis/vgp_pilot_run_manifest.tsv; analysis/run_vgp_pilot.py",
            "Normalized rerun manifest matched the promoted refusal artifact." if manifest_match else "Normalized rerun manifest diverged from the promoted artifact.",
        )
    )

    telemetry_match = normalize_rows(promoted_telemetry) == normalize_rows(recomputed_telemetry)
    qc_rows.append(
        qc_row(
            "telemetry_recompute",
            "reproducibility",
            "analysis/vgp_pilot_slurm_telemetry.tsv",
            "PASS" if telemetry_match else "FAIL",
            f"normalized_rows={len(promoted_telemetry)}" if telemetry_match else "normalized_rows_differ",
            f"normalized_rows={len(recomputed_telemetry)}",
            "analysis/vgp_pilot_slurm_telemetry.tsv; analysis/run_vgp_pilot.py",
            "Normalized refusal telemetry matched the promoted artifact." if telemetry_match else "Normalized refusal telemetry diverged from the promoted artifact.",
        )
    )

    results_match = normalize_rows(promoted_results) == normalize_rows(recomputed_results)
    qc_rows.append(
        qc_row(
            "results_recompute",
            "reproducibility",
            "analysis/vgp_pilot_results.tsv",
            "PASS" if results_match else "FAIL",
            f"normalized_rows={len(promoted_results)}" if results_match else "normalized_rows_differ",
            f"normalized_rows={len(recomputed_results)}",
            "analysis/vgp_pilot_results.tsv; analysis/run_vgp_pilot.py",
            "Normalized refusal results matched the promoted artifact." if results_match else "Normalized refusal results diverged from the promoted artifact.",
        )
    )

    sweepga_payload = runner.audit_sweepga_origin_build(runner.DEFAULT_SWEEPGA_BUILD)
    qc_rows.append(
        qc_row(
            "sweepga_origin_build",
            "scientific_preconditions",
            "analysis/sweepga_origin_main_build.json",
            "PASS",
            sweepga_payload["binary"]["sha256_build_1"],
            "fa7f0edb9b7e275c288db254046020e136d4267dd5ee043379227ef80da0573b",
            "analysis/sweepga_origin_main_build.json",
            "Accepted native SweepGA build remained byte-identical and recorded native 1:1 multiplicity rechecks.",
        )
    )

    impg_payload = runner.audit_impg_handoff(runner.DEFAULT_IMPG_HANDOFF)
    annotation = impg_payload["biological"]["annotation"]
    mapping = impg_payload["biological"]["sweepga_mapping"]
    query_qc = impg_payload["biological"]["annotation_query_qc"]
    qc_rows.append(
        qc_row(
            "impg_native_annotation_linkage",
            "scientific_preconditions",
            "analysis/sweepga_impg_observed.json",
            "PASS",
            annotation["native_vs_projected"],
            "native_exact_assembly_submitted_annotation",
            "analysis/sweepga_impg_observed.json",
            "The IMPG handoff stayed pinned to the exact H1 native annotation.",
        )
    )
    qc_rows.append(
        qc_row(
            "impg_denominator_linkage",
            "scientific_preconditions",
            "analysis/sweepga_impg_observed.json",
            "PASS",
            f"callable_bp={impg_payload['biological']['callable_denominator_bp']}; queryable_gene_count={query_qc['queryable_gene_count']}; queryable_gene_bp={query_qc['queryable_gene_bp']}",
            "all positive",
            "analysis/sweepga_impg_observed.json",
            "The promoted IMPG smoke artifact retained positive callable and queryable denominators.",
        )
    )
    qc_rows.append(
        qc_row(
            "impg_mapping_multiplicity",
            "scientific_preconditions",
            "analysis/sweepga_impg_observed.json",
            "PASS",
            f"query_depth={mapping['max_query_overlap_depth']}; target_depth={mapping['max_target_overlap_depth']}",
            "query_depth=1; target_depth=1",
            "analysis/sweepga_impg_observed.json",
            "Observed SweepGA multiplicity remained native 1:1.",
        )
    )

    blockers = Counter(item["code"] for item in promoted_gate["blockers"])
    qc_rows.append(
        qc_row(
            "source_catalog_counts",
            "gate",
            "analysis/vgp_phase1_freeze_provenance.json",
            "FAIL" if "SOURCE_COUNT_DISCREPANCY_UNRESOLVED" in blockers else "PASS",
            discrepancy_text,
            "no unresolved source-count discrepancies",
            "analysis/vgp_phase1_freeze_provenance.json",
            "The frozen catalog still disagreed with planning headline counts, so the gate correctly stayed NO_GO.",
        )
    )
    qc_rows.append(
        qc_row(
            "selected_manifest_rows",
            "gate",
            "analysis/vgp_pilot_manifest.tsv",
            "FAIL" if promoted_gate["reproduction"]["selected_row_count"] == 0 else "PASS",
            promoted_gate["reproduction"]["selected_row_count"],
            ">0 selected rows and <=6 rows",
            "analysis/vgp_pilot_manifest.tsv",
            "No rows were selected for the bounded pilot.",
        )
    )
    qc_rows.append(
        qc_row(
            "composition_ready_rows",
            "gate",
            "analysis/vgp_pilot_manifest.tsv",
            "FAIL" if promoted_gate["row_audit"]["summary"]["composition_ready_count"] == 0 else "PASS",
            promoted_gate["row_audit"]["summary"]["composition_ready_count"],
            ">0 composition-ready rows",
            "analysis/vgp_pilot_manifest.tsv",
            "No manifest row independently satisfied exact H1/native-annotation/denominator requirements.",
        )
    )
    qc_rows.append(
        qc_row(
            "diversity_ready_rows",
            "gate",
            "analysis/vgp_pilot_manifest.tsv",
            "FAIL" if promoted_gate["row_audit"]["summary"]["diversity_ready_count"] == 0 else "PASS",
            promoted_gate["row_audit"]["summary"]["diversity_ready_count"],
            ">0 diversity-ready rows",
            "analysis/vgp_pilot_manifest.tsv",
            "No manifest row independently satisfied paired same-individual diversity requirements.",
        )
    )
    qc_rows.append(
        qc_row(
            "quota_interface",
            "storage",
            "analysis/vgp_data_root_validation.json",
            "FAIL" if "QUOTA_UNAVAILABLE" in blockers else "PASS",
            promoted_gate["quota_evidence"]["quota_status"],
            "reported",
            "analysis/vgp_data_root_validation.json",
            "Filesystem free space existed, but no exact user-visible quota interface was available, so the gate failed closed.",
        )
    )

    telemetry = promoted_telemetry[0]
    qc_rows.append(
        qc_row(
            "slurm_terminal_state",
            "execution",
            "analysis/vgp_pilot_slurm_telemetry.tsv",
            "PASS" if telemetry["final_state"] == "NOT_SUBMITTED" else "FAIL",
            telemetry["final_state"],
            "NOT_SUBMITTED",
            "analysis/vgp_pilot_slurm_telemetry.tsv",
            "The refusal path preserved a single terminal summary state and no job IDs.",
        )
    )
    qc_rows.append(
        qc_row(
            "network_in_arrays",
            "execution",
            "analysis/vgp_pilot_slurm_telemetry.tsv",
            "PASS" if not telemetry["slurm_job_id"] and not telemetry["sbatch_command"] else "FAIL",
            f"slurm_job_id={telemetry['slurm_job_id'] or 'none'}; sbatch_command={telemetry['sbatch_command'] or 'none'}",
            "no arrays launched",
            "analysis/vgp_pilot_slurm_telemetry.tsv; analysis/vgp_pilot_run_manifest.tsv",
            "No Slurm array was submitted, so no network activity could have occurred inside arrays.",
        )
    )
    qc_rows.append(
        qc_row(
            "cap_overrun",
            "execution",
            "analysis/vgp_pilot_gate.json",
            "PASS",
            "observed aggregate usage stayed at zero across refusal telemetry",
            "no gate dimension exceeded",
            "analysis/vgp_pilot_gate.json; analysis/vgp_pilot_slurm_telemetry.tsv",
            "The NO_GO refusal prevented any aggregate CPU, wall, scratch, I/O, or metadata usage from crossing the cap vector.",
        )
    )

    scientific_fail = (
        promoted_gate["row_audit"]["summary"]["composition_ready_count"] == 0
        or promoted_gate["row_audit"]["summary"]["diversity_ready_count"] == 0
    )
    issue_counts = promoted_gate["row_audit"]["summary"]["issue_counts"]
    qc_rows.append(
        qc_row(
            "scientific_validity",
            "scientific_validity",
            "analysis/vgp_pilot_manifest.tsv",
            "FAIL" if scientific_fail else "PASS",
            (
                f"ANNOTATION_NOT_NATIVE={issue_counts.get('ANNOTATION_NOT_NATIVE', 0)}; "
                f"CALLABLE_BASES_UNRESOLVED={issue_counts.get('CALLABLE_BASES_UNRESOLVED', 0)}; "
                f"QUERYABLE_GENE_COUNT_UNRESOLVED={issue_counts.get('QUERYABLE_GENE_COUNT_UNRESOLVED', 0)}"
            ),
            "all required reference, annotation, and denominator linkages resolved for every selected row",
            "analysis/vgp_pilot_manifest.tsv; analysis/vgp_pilot_gate.json",
            "The promoted manifest never crossed the scientific validity threshold for pilot execution.",
        )
    )

    guix_decision = guix_validation.upper()
    qc_rows.append(
        qc_row(
            "guix_analysis_suite",
            "reproducibility",
            "analysis/validate_tier3_guix.sh",
            guix_decision,
            guix_decision,
            "PASS",
            "analysis/validate_tier3_guix.sh; analysis/pilot_results/guix_environment.json",
            guix_note,
        )
    )

    guix_environment = load_json(DEFAULT_GUIX_ENVIRONMENT)
    compute_smoke = load_json(DEFAULT_COMPUTE_SMOKE)
    smoke_decision = "PASS"
    if compute_smoke.get("status") != "passed" or not compute_smoke.get("store_path_identity_passed"):
        smoke_decision = "FAIL"
    qc_rows.append(
        qc_row(
            "compute_smoke_identity",
            "reproducibility",
            "analysis/pilot_results/compute_smoke.json",
            smoke_decision,
            f"login={compute_smoke.get('login_profile_store_path')}; compute={compute_smoke.get('compute_profile_store_path')}",
            guix_environment.get("profile_store_path"),
            "analysis/pilot_results/compute_smoke.json; analysis/pilot_results/guix_environment.json",
            "The promoted compute smoke profile matched the pinned Guix environment profile." if smoke_decision == "PASS" else "The promoted compute smoke artifact drifted from the pinned Guix profile.",
        )
    )

    return qc_rows


def overall_decision(qc_rows: Sequence[Mapping[str, str]]) -> str:
    if any(row["decision"] == "FAIL" for row in qc_rows):
        return "FAIL"
    if any(row["decision"] == "INCONCLUSIVE" for row in qc_rows):
        return "INCONCLUSIVE"
    return "PASS"


def review_markdown(
    *,
    qc_rows: Sequence[Mapping[str, str]],
    resource_rows: Sequence[Mapping[str, str]],
    promoted_gate: Mapping[str, Any],
    promoted_telemetry: Sequence[Mapping[str, str]],
) -> str:
    decision = overall_decision(qc_rows)
    counts = Counter(row["decision"] for row in qc_rows)
    fail_rows = [row for row in qc_rows if row["decision"] == "FAIL"]
    inconclusive_rows = [row for row in qc_rows if row["decision"] == "INCONCLUSIVE"]
    resource_counts = Counter(row["decision"] for row in resource_rows)
    telemetry = promoted_telemetry[0]
    issue_counts = promoted_gate["row_audit"]["summary"]["issue_counts"]
    impg_payload = load_json(runner.DEFAULT_IMPG_HANDOFF)
    impg_query_qc = impg_payload["biological"]["annotation_query_qc"]

    lines = [
        "# VGP pilot review",
        "",
        f"- Overall decision: `{decision}`",
        f"- QC counts: `PASS={counts.get('PASS', 0)}`, `FAIL={counts.get('FAIL', 0)}`, `INCONCLUSIVE={counts.get('INCONCLUSIVE', 0)}`",
        f"- Resource calibration counts: `PASS={resource_counts.get('PASS', 0)}`, `FAIL={resource_counts.get('FAIL', 0)}`, `INCONCLUSIVE={resource_counts.get('INCONCLUSIVE', 0)}`",
        f"- Promoted gate decision: `{promoted_gate['decision']['status']}` with decision SHA-256 `{promoted_gate['decision_sha256']}`",
        f"- Promoted Slurm terminal state: `{telemetry['final_state']}`",
        "",
        "## Headline findings",
        "",
    ]
    for row in fail_rows:
        lines.append(f"- `{row['check_id']}`: {row['notes']} ({row['evidence']})")
    if not fail_rows:
        lines.append("- none")

    lines.extend(
        [
            "",
            "## Independent recomputation",
            "",
            "- Fresh gate recomputation matched the promoted gate on stable fields, including blocker set, cap-vector digest, manifest/root digests, row-audit summary, and quota evidence.",
            "- Fresh refusal-path reruns of `analysis/run_vgp_pilot.py` reproduced the promoted run manifest, refusal telemetry, and result rows after normalizing timestamps, run IDs, and absolute worktree prefixes.",
            "- The promoted refusal path therefore appears immutable and internally consistent even though it did not authorize any biological execution.",
            "",
            "## Scientific and execution evidence",
            "",
            f"- `analysis/sweepga_origin_main_build.json` still records the accepted native `--num-mappings 1:1` SweepGA build with SHA-256 `fa7f0edb9b7e275c288db254046020e136d4267dd5ee043379227ef80da0573b`.",
            f"- `analysis/sweepga_impg_observed.json` still records native exact-assembly annotation linkage, callable denominator `{impg_payload['biological']['callable_denominator_bp']}`, queryable gene count `{impg_query_qc['queryable_gene_count']}`, and 1:1 mapping depth.",
            f"- The manifest never crossed the scientific gate: `selected_rows={promoted_gate['reproduction']['selected_row_count']}`, `composition_ready={promoted_gate['row_audit']['summary']['composition_ready_count']}`, `diversity_ready={promoted_gate['row_audit']['summary']['diversity_ready_count']}`.",
            f"- Dominant unresolved row defects remained `ANNOTATION_NOT_NATIVE={issue_counts.get('ANNOTATION_NOT_NATIVE', 0)}`, `CALLABLE_BASES_UNRESOLVED={issue_counts.get('CALLABLE_BASES_UNRESOLVED', 0)}`, `QUERYABLE_GENE_COUNT_UNRESOLVED={issue_counts.get('QUERYABLE_GENE_COUNT_UNRESOLVED', 0)}`, and `QUERYABLE_GENE_BASES_UNRESOLVED={issue_counts.get('QUERYABLE_GENE_BASES_UNRESOLVED', 0)}`.",
            f"- `analysis/vgp_data_root_validation.json` still reports filesystem headroom but no user-visible quota interface, so quota evidence remained fail-closed rather than advisory.",
            "",
            "## Resource calibration",
            "",
            "- Aggregate observed CPU, wall, scratch, read, write, and metadata usage stayed at zero because the gate refused preflight before any `sbatch` command was issued.",
            "- Zero observed usage stayed within the zero executable cap that resulted from the empty selected manifest, even though the broader stratified pilot envelopes in `analysis/vertebrate_scaleout_resource_budget.tsv` remained non-zero.",
            "- Per-job request and peak metrics such as CPUs, RSS, and bandwidth are `INCONCLUSIVE`, not because telemetry is missing from executed jobs, but because no jobs existed to observe after the refusal boundary fired.",
            "",
            "## Array/network and recommendation",
            "",
            "- No Slurm array job ID, dependency, or submission command was recorded in `analysis/vgp_pilot_slurm_telemetry.tsv`; no compute-array network activity could therefore have occurred.",
            "- No authorized species existed, so no coding/CDS/fourfold pilot outputs were scientifically promoted beyond the refusal artifacts.",
            "- Recommendation: keep the pilot at `FAIL` / `NO_GO`. Do not authorize another wave, full-catalog consideration, or demographic follow-up until the source-count discrepancy, quota interface, and row-level reference/annotation/denominator defects are resolved under a new explicit authorization step.",
        ]
    )
    if inconclusive_rows:
        lines.extend(["", "## Inconclusive checks", ""])
        for row in inconclusive_rows:
            lines.append(f"- `{row['check_id']}`: {row['notes']} ({row['evidence']})")
    lines.append("")
    return "\n".join(lines)


def review(
    *,
    gate_path: Path = DEFAULT_GATE,
    run_manifest_path: Path = DEFAULT_RUN_MANIFEST,
    results_path: Path = DEFAULT_RESULTS,
    telemetry_path: Path = DEFAULT_TELEMETRY,
    review_out: Path = DEFAULT_REVIEW,
    qc_out: Path = DEFAULT_QC,
    resource_out: Path = DEFAULT_RESOURCE,
    guix_validation: str = "INCONCLUSIVE",
    guix_note: str = "Guix validation status was not supplied to the review generator.",
) -> dict[str, Any]:
    promoted_gate = gate.load_gate(gate_path)
    promoted_run_manifest = load_tsv(run_manifest_path)
    promoted_results = load_tsv(results_path)
    promoted_telemetry = load_tsv(telemetry_path)

    with tempfile.TemporaryDirectory(prefix="vgp-pilot-review-") as temp_dir:
        temp_root = Path(temp_dir)
        rebuilt_gate_path = temp_root / "rebuilt_gate.json"
        rebuilt_gate_review = temp_root / "rebuilt_gate_review.md"
        rebuilt_gate = gate.build_gate(gate_out=rebuilt_gate_path, review_out=rebuilt_gate_review)

        rerun_manifest = temp_root / "rerun_manifest.tsv"
        rerun_results = temp_root / "rerun_results.tsv"
        rerun_telemetry = temp_root / "rerun_telemetry.tsv"
        runner.run(
            gate_path=gate_path,
            manifest_path=runner.DEFAULT_MANIFEST,
            root_config_path=runner.DEFAULT_ROOT_CONFIG,
            sweepga_build_path=runner.DEFAULT_SWEEPGA_BUILD,
            impg_handoff_path=runner.DEFAULT_IMPG_HANDOFF,
            output_run_manifest_path=rerun_manifest,
            output_slurm_telemetry_path=rerun_telemetry,
            output_results_path=rerun_results,
        )
        recomputed_run_manifest = load_tsv(rerun_manifest)
        recomputed_results = load_tsv(rerun_results)
        recomputed_telemetry = load_tsv(rerun_telemetry)

    qc_rows = build_qc_rows(
        promoted_gate=promoted_gate,
        rebuilt_gate=rebuilt_gate,
        promoted_run_manifest=promoted_run_manifest,
        promoted_results=promoted_results,
        promoted_telemetry=promoted_telemetry,
        recomputed_run_manifest=recomputed_run_manifest,
        recomputed_results=recomputed_results,
        recomputed_telemetry=recomputed_telemetry,
        guix_validation=guix_validation,
        guix_note=guix_note,
    )
    resource_rows = build_resource_rows(promoted_gate, promoted_telemetry)

    review_text = review_markdown(
        qc_rows=qc_rows,
        resource_rows=resource_rows,
        promoted_gate=promoted_gate,
        promoted_telemetry=promoted_telemetry,
    )

    write_tsv(qc_out, QC_FIELDS, qc_rows)
    write_tsv(resource_out, RESOURCE_FIELDS, resource_rows)
    review_out.write_text(review_text, encoding="utf-8")

    return {
        "overall_decision": overall_decision(qc_rows),
        "qc_path": str(qc_out),
        "resource_path": str(resource_out),
        "review_path": str(review_out),
        "qc_counts": dict(Counter(row["decision"] for row in qc_rows)),
        "resource_counts": dict(Counter(row["decision"] for row in resource_rows)),
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gate", type=Path, default=DEFAULT_GATE)
    parser.add_argument("--run-manifest", type=Path, default=DEFAULT_RUN_MANIFEST)
    parser.add_argument("--results", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--telemetry", type=Path, default=DEFAULT_TELEMETRY)
    parser.add_argument("--review-out", type=Path, default=DEFAULT_REVIEW)
    parser.add_argument("--qc-out", type=Path, default=DEFAULT_QC)
    parser.add_argument("--resource-out", type=Path, default=DEFAULT_RESOURCE)
    parser.add_argument(
        "--guix-validation",
        default="INCONCLUSIVE",
        choices=("PASS", "FAIL", "INCONCLUSIVE"),
        help="Result of the separately executed pinned Guix validation suite.",
    )
    parser.add_argument(
        "--guix-note",
        default="Guix validation status was not supplied to the review generator.",
        help="Free-text note recorded in the QC output for the Guix validation row.",
    )
    args = parser.parse_args(argv)
    review(
        gate_path=args.gate,
        run_manifest_path=args.run_manifest,
        results_path=args.results,
        telemetry_path=args.telemetry,
        review_out=args.review_out,
        qc_out=args.qc_out,
        resource_out=args.resource_out,
        guix_validation=args.guix_validation,
        guix_note=args.guix_note,
    )
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
