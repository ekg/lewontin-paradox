#!/usr/bin/env python3
"""Fail-closed compute entrypoint for the bounded VGP Tier 3 pilot."""

from __future__ import annotations

import argparse
import csv
import json
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Sequence

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from analysis import gate_vgp_pilot as gate
from analysis.tier3_common import Tier3ValidationError, sha256_file


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_GATE = PROJECT_ROOT / "analysis" / "vgp_pilot_gate.json"
DEFAULT_MANIFEST = PROJECT_ROOT / "analysis" / "vgp_pilot_manifest.tsv"
DEFAULT_ROOT_CONFIG = PROJECT_ROOT / "analysis" / "vgp_data_root_config.json"
DEFAULT_SWEEPGA_BUILD = PROJECT_ROOT / "analysis" / "sweepga_origin_main_build.json"
DEFAULT_IMPG_HANDOFF = PROJECT_ROOT / "analysis" / "sweepga_impg_observed.json"
DEFAULT_OUTPUT_RUN_MANIFEST = PROJECT_ROOT / "analysis" / "vgp_pilot_run_manifest.tsv"
DEFAULT_OUTPUT_SLURM_TELEMETRY = PROJECT_ROOT / "analysis" / "vgp_pilot_slurm_telemetry.tsv"
DEFAULT_OUTPUT_RESULTS = PROJECT_ROOT / "analysis" / "vgp_pilot_results.tsv"

RUN_MANIFEST_FIELDS = [
    "run_id",
    "generated_at_utc",
    "record_type",
    "status",
    "candidate_id",
    "gate_decision",
    "gate_decision_sha256",
    "manifest_digest",
    "root_contract_digest",
    "cap_vector_digest",
    "sweepga_origin_build_sha256",
    "impg_handoff_sha256",
    "failure_code",
    "failure_source",
    "failure_message",
    "notes",
]

SLURM_TELEMETRY_FIELDS = [
    "run_id",
    "generated_at_utc",
    "record_type",
    "status",
    "candidate_id",
    "sbatch_command",
    "slurm_job_id",
    "slurm_array_job_id",
    "dependency",
    "requested_cpus",
    "requested_memory_gib",
    "requested_wall_hours",
    "requested_scratch_gb",
    "requested_read_gb",
    "requested_write_gb",
    "retry_index",
    "final_state",
    "max_rss_gib",
    "elapsed_seconds",
    "cpu_time_seconds",
    "scratch_peak_gb",
    "io_read_gb",
    "io_write_gb",
    "metadata_operations",
    "failure_code",
    "notes",
]

RESULT_FIELDS = [
    "run_id",
    "generated_at_utc",
    "record_type",
    "status",
    "candidate_id",
    "scientific_name",
    "result_scope",
    "metric",
    "numerator",
    "denominator",
    "value",
    "exclusion_reason",
    "failure_code",
    "failure_source",
    "notes",
]


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def make_run_id() -> str:
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    return f"vgp-pilot-run-{stamp}"


def atomic_write_tsv(path: Path, fieldnames: Sequence[str], rows: Sequence[Mapping[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", newline="", dir=path.parent, prefix=f".{path.name}.", delete=False
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
        handle.flush()
        temp_path = Path(handle.name)
    temp_path.replace(path)


def read_json_object(path: Path, label: str) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as error:
        raise Tier3ValidationError(f"invalid {label}: {error}") from error
    if not isinstance(value, dict):
        raise Tier3ValidationError(f"{label} must be a JSON object")
    return value


def classify_failure_code(error: Exception) -> str:
    text = str(error)
    if "gate decision is NO_GO" in text:
        return "GATE_NO_GO"
    if "manifest digest mismatch" in text:
        return "MANIFEST_DIGEST_MISMATCH"
    if "root contract digest mismatch" in text:
        return "ROOT_CONTRACT_DIGEST_MISMATCH"
    if "cap vector digest mismatch" in text:
        return "CAP_VECTOR_DIGEST_MISMATCH"
    if "gate decision hash does not match" in text or "cap vector hash does not match the gate payload" in text:
        return "GATE_TAMPERED"
    if "sweepga" in text.lower() or "impg" in text.lower():
        return "ENVIRONMENT_AUDIT_FAILED"
    return "PRECHECK_FAILED"


def audit_sweepga_origin_build(path: Path) -> dict[str, Any]:
    payload = read_json_object(path, "SweepGA origin/main build record")
    if payload.get("completion_gate_passed") is not True:
        raise Tier3ValidationError("sweepga origin/main build completion gate did not pass")
    if payload.get("acceptance_status") != "accepted_native_long_option_after_user_clarification":
        raise Tier3ValidationError("sweepga origin/main build acceptance status drifted from the recorded approval")
    resolution = payload.get("resolution", {})
    if resolution.get("accepted_syntax") != "--num-mappings 1:1":
        raise Tier3ValidationError("sweepga origin/main build does not record native --num-mappings 1:1")
    if resolution.get("all_three_whole_haplotype_mappings_completed") is not True:
        raise Tier3ValidationError("sweepga origin/main build lacks completed whole-haplotype mappings")
    if resolution.get("all_three_native_multiplicity_rechecks_passed") is not True:
        raise Tier3ValidationError("sweepga origin/main build lacks native multiplicity rechecks")
    if resolution.get("all_three_impg_handoffs_completed") is not True:
        raise Tier3ValidationError("sweepga origin/main build lacks completed IMPG handoffs")
    binary = payload.get("binary", {})
    if binary.get("sha256_build_1") != "fa7f0edb9b7e275c288db254046020e136d4267dd5ee043379227ef80da0573b":
        raise Tier3ValidationError("sweepga origin/main build SHA-256 differs from the accepted binary")
    return payload


def audit_impg_handoff(path: Path) -> dict[str, Any]:
    payload = read_json_object(path, "SweepGA to IMPG handoff record")
    environment = payload.get("environment", {})
    biological = payload.get("biological", {})
    annotation = biological.get("annotation", {})
    query_qc = biological.get("annotation_query_qc", {})
    mapping = biological.get("sweepga_mapping", {})
    if environment.get("guix_profile", "").startswith("/gnu/store/") is not True:
        raise Tier3ValidationError("impg handoff record lacks a pinned Guix profile")
    if environment.get("sweepga_commit") != "018e4ce49d2c125820e0ac50dc5feaa02d423683":
        raise Tier3ValidationError("impg handoff record drifted from the approved SweepGA origin/main commit")
    if environment.get("impg_commit") != "101df81eb28a809c8fac97d297acd9fcfbbfa048":
        raise Tier3ValidationError("impg handoff record drifted from the approved IMPG commit")
    if annotation.get("native_vs_projected") != "native_exact_assembly_submitted_annotation":
        raise Tier3ValidationError("impg handoff record is not linked to native H1 annotation")
    if biological.get("callable_denominator_bp", 0) <= 0:
        raise Tier3ValidationError("impg handoff record has no positive callable denominator")
    if query_qc.get("queryable_gene_count", 0) <= 0:
        raise Tier3ValidationError("impg handoff record has no positive queryable gene count")
    if biological.get("normalized_bcf_indexed") is not True or biological.get("normalized_vcf_indexed") is not True:
        raise Tier3ValidationError("impg handoff record lacks normalized indexed VCF/BCF outputs")
    if mapping.get("max_query_overlap_depth") != 1 or mapping.get("max_target_overlap_depth") != 1:
        raise Tier3ValidationError("impg handoff record violates SweepGA native 1:1 multiplicity")
    return payload


def refusal_run_manifest_rows(
    *,
    run_id: str,
    generated_at_utc: str,
    gate_path: Path,
    gate_payload: Mapping[str, Any] | None,
    manifest_path: Path,
    root_config_path: Path,
    sweepga_build_path: Path,
    impg_handoff_path: Path,
    error: Exception,
) -> list[dict[str, str]]:
    decision_sha = ""
    gate_decision = "UNREADABLE"
    cap_digest = ""
    blockers: list[Mapping[str, str]] = []
    if gate_payload is not None:
        decision_sha = str(gate_payload.get("decision_sha256", ""))
        gate_decision = str(gate_payload.get("decision", {}).get("status", "UNREADABLE"))
        cap_digest = str(gate_payload.get("authorization_boundary", {}).get("cap_vector_digest", ""))
        blockers = list(gate_payload.get("blockers", []))
    base_row = {
        "run_id": run_id,
        "generated_at_utc": generated_at_utc,
        "record_type": "run_summary",
        "status": "refused_preflight",
        "candidate_id": "",
        "gate_decision": gate_decision,
        "gate_decision_sha256": decision_sha,
        "manifest_digest": sha256_file(manifest_path),
        "root_contract_digest": sha256_file(root_config_path),
        "cap_vector_digest": cap_digest,
        "sweepga_origin_build_sha256": sha256_file(sweepga_build_path),
        "impg_handoff_sha256": sha256_file(impg_handoff_path),
        "failure_code": classify_failure_code(error),
        "failure_source": str(gate_path),
        "failure_message": str(error),
        "notes": "no species were submitted; no Slurm command was issued",
    }
    rows = [base_row]
    for blocker in blockers:
        rows.append(
            {
                **base_row,
                "record_type": "gate_blocker",
                "failure_code": blocker["code"],
                "failure_source": blocker["source"],
                "failure_message": blocker["message"],
                "notes": "copied from gate blockers under the current refusal boundary",
            }
        )
    return rows


def refusal_slurm_telemetry_rows(
    *,
    run_id: str,
    generated_at_utc: str,
    error: Exception,
) -> list[dict[str, str]]:
    return [
        {
            "run_id": run_id,
            "generated_at_utc": generated_at_utc,
            "record_type": "run_summary",
            "status": "refused_preflight",
            "candidate_id": "",
            "sbatch_command": "",
            "slurm_job_id": "",
            "slurm_array_job_id": "",
            "dependency": "",
            "requested_cpus": "",
            "requested_memory_gib": "",
            "requested_wall_hours": "",
            "requested_scratch_gb": "",
            "requested_read_gb": "",
            "requested_write_gb": "",
            "retry_index": "0",
            "final_state": "NOT_SUBMITTED",
            "max_rss_gib": "",
            "elapsed_seconds": "0",
            "cpu_time_seconds": "0",
            "scratch_peak_gb": "0",
            "io_read_gb": "0",
            "io_write_gb": "0",
            "metadata_operations": "0",
            "failure_code": classify_failure_code(error),
            "notes": "gate refusal prevented sbatch submission and sacct telemetry collection",
        }
    ]


def refusal_result_rows(
    *,
    run_id: str,
    generated_at_utc: str,
    gate_path: Path,
    gate_payload: Mapping[str, Any] | None,
    error: Exception,
) -> list[dict[str, str]]:
    rows = [
        {
            "run_id": run_id,
            "generated_at_utc": generated_at_utc,
            "record_type": "run_summary",
            "status": "refused_preflight",
            "candidate_id": "",
            "scientific_name": "",
            "result_scope": "pilot",
            "metric": "validated_species_count",
            "numerator": "0",
            "denominator": "0",
            "value": "0",
            "exclusion_reason": str(error),
            "failure_code": classify_failure_code(error),
            "failure_source": str(gate_path),
            "notes": "no species passed the authorization boundary for compute",
        }
    ]
    blockers = list(gate_payload.get("blockers", [])) if gate_payload is not None else []
    for blocker in blockers:
        rows.append(
            {
                "run_id": run_id,
                "generated_at_utc": generated_at_utc,
                "record_type": "gate_blocker",
                "status": "refused_preflight",
                "candidate_id": "",
                "scientific_name": "",
                "result_scope": "authorization_boundary",
                "metric": "preflight_blocker",
                "numerator": "",
                "denominator": "",
                "value": "",
                "exclusion_reason": blocker["message"],
                "failure_code": blocker["code"],
                "failure_source": blocker["source"],
                "notes": "no aggregate coding/CDS/fourfold result was authorized",
            }
        )
    return rows


def run(
    *,
    gate_path: Path = DEFAULT_GATE,
    manifest_path: Path = DEFAULT_MANIFEST,
    root_config_path: Path = DEFAULT_ROOT_CONFIG,
    sweepga_build_path: Path = DEFAULT_SWEEPGA_BUILD,
    impg_handoff_path: Path = DEFAULT_IMPG_HANDOFF,
    output_run_manifest_path: Path = DEFAULT_OUTPUT_RUN_MANIFEST,
    output_slurm_telemetry_path: Path = DEFAULT_OUTPUT_SLURM_TELEMETRY,
    output_results_path: Path = DEFAULT_OUTPUT_RESULTS,
) -> dict[str, Any]:
    run_id = make_run_id()
    generated_at_utc = utc_now()
    gate_payload: Mapping[str, Any] | None
    try:
        gate_payload = read_json_object(gate_path, "VGP pilot gate")
    except Tier3ValidationError:
        gate_payload = None
    refusal_error: Exception
    try:
        audit_sweepga_origin_build(sweepga_build_path)
        audit_impg_handoff(impg_handoff_path)
        gate.authorize_gate_action(gate_path, manifest_path, root_config_path, "compute")
    except (Tier3ValidationError, FileNotFoundError, KeyError, ValueError) as error:
        refusal_error = error
    else:  # pragma: no cover - current repository state intentionally has no authorized compute execution path.
        raise Tier3ValidationError("authorized VGP pilot compute is not implemented in this repository state")

    run_rows = refusal_run_manifest_rows(
        run_id=run_id,
        generated_at_utc=generated_at_utc,
        gate_path=gate_path,
        gate_payload=gate_payload,
        manifest_path=manifest_path,
        root_config_path=root_config_path,
        sweepga_build_path=sweepga_build_path,
        impg_handoff_path=impg_handoff_path,
        error=refusal_error,
    )
    telemetry_rows = refusal_slurm_telemetry_rows(
        run_id=run_id,
        generated_at_utc=generated_at_utc,
        error=refusal_error,
    )
    result_rows = refusal_result_rows(
        run_id=run_id,
        generated_at_utc=generated_at_utc,
        gate_path=gate_path,
        gate_payload=gate_payload,
        error=refusal_error,
    )
    atomic_write_tsv(output_run_manifest_path, RUN_MANIFEST_FIELDS, run_rows)
    atomic_write_tsv(output_slurm_telemetry_path, SLURM_TELEMETRY_FIELDS, telemetry_rows)
    atomic_write_tsv(output_results_path, RESULT_FIELDS, result_rows)
    return {
        "run_id": run_id,
        "status": "refused_preflight",
        "failure_code": classify_failure_code(refusal_error),
        "run_manifest": str(output_run_manifest_path),
        "slurm_telemetry": str(output_slurm_telemetry_path),
        "results": str(output_results_path),
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gate", type=Path, default=DEFAULT_GATE)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--root-config", type=Path, default=DEFAULT_ROOT_CONFIG)
    parser.add_argument("--sweepga-build", type=Path, default=DEFAULT_SWEEPGA_BUILD)
    parser.add_argument("--impg-handoff", type=Path, default=DEFAULT_IMPG_HANDOFF)
    parser.add_argument("--run-manifest-out", type=Path, default=DEFAULT_OUTPUT_RUN_MANIFEST)
    parser.add_argument("--slurm-telemetry-out", type=Path, default=DEFAULT_OUTPUT_SLURM_TELEMETRY)
    parser.add_argument("--results-out", type=Path, default=DEFAULT_OUTPUT_RESULTS)
    args = parser.parse_args(argv)
    run(
        gate_path=args.gate,
        manifest_path=args.manifest,
        root_config_path=args.root_config,
        sweepga_build_path=args.sweepga_build,
        impg_handoff_path=args.impg_handoff,
        output_run_manifest_path=args.run_manifest_out,
        output_slurm_telemetry_path=args.slurm_telemetry_out,
        output_results_path=args.results_out,
    )
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
