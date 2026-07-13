#!/usr/bin/env python3
"""Finalize the fail-closed Tier 3a VGP execution inventory.

The frozen pilot registry currently contains no checksum-locked, eligible VGP
individual.  This command verifies that fact together with the pinned Guix and
compute-smoke records, then emits deterministic structured-missingness
artifacts.  It deliberately refuses to emit an empty result if the pilot is
later promoted: an eligible tuple must be run through the Slurm workflow.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import os
import sys
import tempfile
from pathlib import Path
from typing import Any, Dict, Mapping, Optional, Sequence

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from analysis.tier3_common import Tier3ValidationError, sha256_file
from analysis.tier3a_vgp_collect import IMPG_COMMIT, WFMASH_COMMIT


POLICY_ID = "tier3-decisions-v1"
SCHEMA_VERSION = "1.0"
RUN_ID = "tier3a-fail-closed-20260713"
PILOT_ID = "vgp_dual_modality_individual"
VGP_PHASE1_FREEZE_SHA256 = "9c58420484a8b76a2d6175b7c26bf709e68bdc726a67fc7541b8c2b5a2fc13a4"
BUFFALO_CORE_SHA256 = "df559451dad94b53ba8675e09811708107a57eeb6ffe8f72b944bcbbf3a1f2eb"
BUFFALO_PROVENANCE_SHA256 = "ed1eee9eafc0c6e8c8e9a7bced7d56f19c5e18186131ad297aacafa52d5bc27a"

# The second entry documents the only locally observed paired VGP-like payload.
# It is inventory evidence, not an eligible row: it is absent from the frozen
# Tier 3 staging manifest and does not have the required audits or denominator.
CANDIDATES = (
    {
        "candidate_id": "vgp_dual_modality_individual",
        "scientific_name": None,
        "taxon_id": None,
        "h1_accession": None,
        "h2_accession": None,
        "source": "predeclared_pilot_registry",
        "blocking_codes": (
            "missing_frozen_dataset_identity",
            "missing_checksum_locked_h1_h2_tuple",
            "missing_deposited_exact_reference_calls_and_callable_mask",
            "missing_native_h1_annotation_audit",
            "missing_phase_identity_audit",
            "missing_collapse_duplication_qc",
        ),
    },
    {
        "candidate_id": "mmyodau2.1_local_inventory",
        "scientific_name": "Myotis daubentonii",
        "taxon_id": 98922,
        "h1_accession": "GCF_963259705.1",
        "h2_accession": "GCA_963242275.1",
        "source": "local_inventory_not_frozen_for_tier3",
        "inventory_only_compressed_h1_sha256": "9e0ffb5048bf3653a338c3e885ab23057717e3dda36f96917f443f505fc5f7f3",
        "inventory_only_compressed_h2_sha256": "11997e23203198417a6b26ad468c92bd43e1ad18cf31ffbaf8e3a8b407e949b6",
        "blocking_codes": (
            "not_in_frozen_eligible_tier3_manifest",
            "missing_exact_buffalo_species_covariate",
            "missing_deposited_exact_reference_calls_and_callable_mask",
            "missing_native_h1_gff_checksum_and_cds_reconstruction_audit",
            "missing_phase_identity_audit",
            "missing_collapse_duplication_qc",
            "pinned_wfmash_not_run_after_preflight_failure",
        ),
    },
)

DATA_COLUMNS = (
    "dataset_id",
    "scientific_name",
    "taxon_id",
    "buffalo_species",
    "buffalo_pred_log10_N",
    "h1_reference_accession_version",
    "h1_fasta_sha256",
    "h2_query_accession_version",
    "h2_fasta_sha256",
    "modality",
    "statistic_label",
    "denominator_kind",
    "callable_bed_sha256",
    "unique_callable_h1_bases",
    "total_snv_numerator",
    "total_denominator",
    "individual_snv_heterozygosity",
    "individual_snv_heterozygosity_ci_low",
    "individual_snv_heterozygosity_ci_high",
    "fourfold_W_snv_numerator",
    "fourfold_W_denominator",
    "pi_W",
    "pi_W_ci_low",
    "pi_W_ci_high",
    "fourfold_S_snv_numerator",
    "fourfold_S_denominator",
    "pi_S",
    "pi_S_ci_low",
    "pi_S_ci_high",
    "pi_S_over_pi_W",
    "pi_S_over_pi_W_ci_low",
    "pi_S_over_pi_W_ci_high",
    "pi_S_over_pi_W_reference_conditioned",
    "annotation_provider",
    "annotation_release",
    "annotation_assembly_accession_version",
    "annotation_fasta_sha256",
    "annotation_gff_sha256",
    "annotation_sequence_region_contig_mapping_sha256",
    "annotation_genetic_code",
    "annotation_native_vs_projected",
    "annotation_cds_reconstruction_audit",
    "exclusion_counts_json",
    "bootstrap_policy",
    "guix_profile_store_path",
    "slurm_job_id",
)

FAILURE_COLUMNS = (
    "candidate_id",
    "scientific_name",
    "taxon_id",
    "h1_accession_version",
    "h2_accession_version",
    "source",
    "eligibility_status",
    "job_status",
    "slurm_job_id",
    "blocking_codes",
    "deposited_modality_status",
    "direct_wfmash_status",
    "impg_status",
    "decision_version",
)


def _read_json(path: Path, label: str) -> Dict[str, Any]:
    try:
        value = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, ValueError) as error:
        raise Tier3ValidationError(f"invalid {label}: {error}") from error
    if not isinstance(value, dict):
        raise Tier3ValidationError(f"{label} must be a JSON object")
    return value


def _atomic_write(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=str(path.parent))
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(data)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_name, path)
    finally:
        if os.path.exists(temporary_name):
            os.unlink(temporary_name)


def _tsv_bytes(columns: Sequence[str], rows: Sequence[Mapping[str, Any]]) -> bytes:
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=columns, dialect="excel-tab", lineterminator="\n")
    writer.writeheader()
    for row in rows:
        writer.writerow({key: row.get(key, "") for key in columns})
    return output.getvalue().encode("utf-8")


def _audit_inputs(
    pilot_registry_path: Path,
    environment_path: Path,
    compute_smoke_path: Path,
) -> Dict[str, Any]:
    registry = _read_json(pilot_registry_path, "pilot registry")
    environment = _read_json(environment_path, "Guix environment record")
    smoke = _read_json(compute_smoke_path, "compute smoke record")
    if registry.get("decision_version") != POLICY_ID or environment.get("decision_version") != POLICY_ID:
        raise Tier3ValidationError(f"pilot registry and Guix environment must use {POLICY_ID}")

    pilots = {
        item.get("pilot_id"): item
        for item in registry.get("pilots", ())
        if isinstance(item, dict) and isinstance(item.get("pilot_id"), str)
    }
    if PILOT_ID not in pilots:
        raise Tier3ValidationError(f"pilot registry is missing {PILOT_ID}")
    pilot = pilots[PILOT_ID]
    eligibility = pilot.get("eligibility")
    if not isinstance(eligibility, str) or not eligibility.startswith("ineligible_pending"):
        raise Tier3ValidationError(f"{PILOT_ID} is no longer pending and must be run")
    if pilot.get("scientific_name") is not None or pilot.get("taxon_id") is not None:
        raise Tier3ValidationError("pending VGP pilot unexpectedly has a frozen individual identity")
    if pilot.get("required_modalities") != [
        "deposited_exact_reference_variants_plus_mask",
        "direct_wfmash_extended_cigar",
    ]:
        raise Tier3ValidationError("VGP pilot modality contract differs from frozen policy")
    optional_impg = pilot.get("optional_impg", {})
    if optional_impg.get("execution_approved") is not False:
        raise Tier3ValidationError("IMPG execution is not approved by the frozen pilot registry")
    if optional_impg.get("phase_sensitive_eligible") is not False:
        raise Tier3ValidationError("IMPG phase-sensitive eligibility must remain false")

    profile = environment.get("profile_store_path")
    if not isinstance(profile, str) or not profile.startswith("/gnu/store/"):
        raise Tier3ValidationError("Guix environment has no valid profile store path")
    versions = environment.get("tool_versions", {})
    if WFMASH_COMMIT not in str(versions.get("wfmash", "")):
        raise Tier3ValidationError("environment does not record the pinned WFMASH commit")
    if environment.get("pack_fallback", {}).get("required") is not False:
        raise Tier3ValidationError("unexpected Guix pack fallback in primary shared-store run")
    if smoke.get("status") != "passed" or smoke.get("decision_version") != POLICY_ID:
        raise Tier3ValidationError("compute smoke must pass under the frozen policy")
    if smoke.get("compute_profile_store_path") != profile or smoke.get("login_profile_store_path") != profile:
        raise Tier3ValidationError("compute/login smoke profile differs from pinned Guix profile")
    if smoke.get("store_path_identity_passed") is not True:
        raise Tier3ValidationError("compute smoke did not prove store-path identity")
    wfmash_smoke = smoke.get("wfmash_extended_cigar", {})
    if wfmash_smoke.get("extended_cigar_passed") is not True:
        raise Tier3ValidationError("compute smoke did not prove extended WFMASH CIGAR")
    if wfmash_smoke.get("unique_mapping_passed") is not True:
        raise Tier3ValidationError("compute smoke did not prove unique mapping handling")
    if smoke.get("bcf_csi_passed") is not True:
        raise Tier3ValidationError("compute smoke did not prove normalized BCF/CSI support")
    return {"registry": registry, "environment": environment, "smoke": smoke, "pilot": pilot}


def _candidate_qc(candidate: Mapping[str, Any], pilot: Mapping[str, Any]) -> Dict[str, Any]:
    local_inventory = candidate["candidate_id"] == "mmyodau2.1_local_inventory"
    return {
        "candidate_id": candidate["candidate_id"],
        "scientific_name": candidate["scientific_name"],
        "taxon_id": candidate["taxon_id"],
        "source": candidate["source"],
        "inclusion_status": "excluded_preflight",
        "blocking_codes": list(candidate["blocking_codes"]),
        "buffalo_covariates": {
            "exact_species_match": False,
            "pred_log10_N": None,
            "congener_substitution_used": False,
        },
        "reference_tuple": {
            "h1_accession_version": candidate["h1_accession"],
            "h2_accession_version": candidate["h2_accession"],
            "h1_uncompressed_fasta_sha256": None,
            "h2_uncompressed_fasta_sha256": None,
            "h1_compressed_payload_sha256_inventory_only": candidate.get("inventory_only_compressed_h1_sha256"),
            "h2_compressed_payload_sha256_inventory_only": candidate.get("inventory_only_compressed_h2_sha256"),
            "exact_reference_tuple_frozen": False,
            "contig_dictionary_audit": "not_run_missing_frozen_tuple",
        },
        "annotation": {
            "primary_gc3_4d_status": "unavailable_missing_native_exact_reference_annotation_audit",
            "provider": None,
            "release": None,
            "assembly_accession_version": None,
            "fasta_sha256": None,
            "gff_sha256": None,
            "sequence_region_contig_mapping_sha256": None,
            "genetic_code": None,
            "native_vs_projected": None,
            "contig_dictionary_passed": False,
            "sampled_cds_reconstruction": "not_run",
            "projected_annotation_used_for_primary": False,
        },
        "phasing_and_collapse": {
            "same_individual_identity_frozen": False,
            "h1_h2_phase_identity_audit_passed": False,
            "collapse_duplication_qc_passed": False,
            "inventory_pair_observed": local_inventory,
        },
        "deposited_modality": {
            "status": "unavailable_missing_calls_and_callable_mask",
            "normalized_bcf_sha256": None,
            "normalized_bcf_csi_sha256": None,
            "callable_mask_sha256": None,
            "unique_callable_h1_bases": None,
        },
        "direct_wfmash_modality": {
            "status": "not_run_preflight_gate_failure",
            "wfmash_commit": WFMASH_COMMIT,
            "raw_paf_sha256": None,
            "accepted_mapping_qc_sha256": None,
            "callable_bed_sha256": None,
            "normalized_bcf_sha256": None,
            "normalized_bcf_csi_sha256": None,
            "unique_callable_h1_bases": None,
            "exclusion_counts": dict.fromkeys(candidate["blocking_codes"], 1),
        },
        "statistics": {
            "statistic_label": "individual_snv_heterozygosity",
            "population_pi": None,
            "heterozygous_snv_numerator": None,
            "total_denominator": None,
            "fourfold_W_denominator": None,
            "fourfold_S_denominator": None,
            "uncertainty": "not_run_statistic_unavailable",
        },
        "dual_modality_concordance": {
            "required": candidate["candidate_id"] == PILOT_ID,
            "status": "not_run_no_eligible_dual_modality_tuple",
            "promotion_passed": False,
        },
        "impg": {
            "commit": IMPG_COMMIT,
            "execution_approved": False,
            "status": "not_run_source_only_in_guix_and_phase_orientation_untrusted",
            "indexed_region_queries": "not_demonstrated_for_real_candidate",
            "boundary_retention": "synthetic_truth_only",
            "exact_deduplication": "synthetic_truth_only",
            "phase_sensitive_eligible": False,
        },
        "job": {
            "status": "not_submitted_gate_failure",
            "slurm_job_id": None,
            "resource_request": {
                "deposited": {"cpus": "2-4", "memory_gb": "16-32", "walltime_hours": "2-4"},
                "direct_wfmash": {"cpus": "4-8", "memory_gb": "32-64", "walltime_hours": "8-24"},
            },
            "resource_telemetry": "not_applicable_no_job_submitted",
        },
        "pilot_registry_eligibility": pilot["eligibility"] if candidate["candidate_id"] == PILOT_ID else None,
    }


def finalize_ineligible_run(
    *,
    pilot_registry_path: Path,
    environment_path: Path,
    compute_smoke_path: Path,
    data_path: Path,
    failure_ledger_path: Path,
    qc_path: Path,
) -> Dict[str, Any]:
    """Write deterministic header-only result, failure, and provenance files."""

    audited = _audit_inputs(pilot_registry_path, environment_path, compute_smoke_path)
    environment, smoke, pilot = audited["environment"], audited["smoke"], audited["pilot"]
    candidates = [_candidate_qc(candidate, pilot) for candidate in CANDIDATES]
    failures = [
        {
            "candidate_id": candidate["candidate_id"],
            "scientific_name": candidate["scientific_name"] or "",
            "taxon_id": candidate["taxon_id"] or "",
            "h1_accession_version": candidate["h1_accession"] or "",
            "h2_accession_version": candidate["h2_accession"] or "",
            "source": candidate["source"],
            "eligibility_status": "excluded_preflight",
            "job_status": "not_submitted_gate_failure",
            "slurm_job_id": "",
            "blocking_codes": ";".join(candidate["blocking_codes"]),
            "deposited_modality_status": "unavailable_missing_calls_and_callable_mask",
            "direct_wfmash_status": "not_run_preflight_gate_failure",
            "impg_status": "not_run_execution_not_approved",
            "decision_version": POLICY_ID,
        }
        for candidate in CANDIDATES
    ]
    qc_record = {
        "schema_version": SCHEMA_VERSION,
        "decision_version": POLICY_ID,
        "run_id": RUN_ID,
        "as_of_environment_record": environment.get("recorded_at"),
        "overall_status": "no_eligible_vgp_individual_tuples",
        "included_individual_count": 0,
        "candidate_count": len(candidates),
        "input_audit": {
            "frozen_analysis_manifest": "absent",
            "pilot_registry_sha256": sha256_file(pilot_registry_path),
            "guix_environment_record_sha256": sha256_file(environment_path),
            "compute_smoke_record_sha256": sha256_file(compute_smoke_path),
            "vgp_phase1_freeze_inventory_sha256": VGP_PHASE1_FREEZE_SHA256,
            "buffalo_core_sha256": BUFFALO_CORE_SHA256,
            "buffalo_provenance_sha256": BUFFALO_PROVENANCE_SHA256,
            "inventory_payloads_reverified_during_finalization": False,
            "raw_payloads_committed": False,
        },
        "environment": {
            "manager": "GNU Guix",
            "channel_commit": environment.get("channel_commit"),
            "channels_sha256": environment.get("channels_sha256"),
            "resolved_channels_sha256": environment.get("resolved_channels_sha256"),
            "manifest_sha256": environment.get("manifest_sha256"),
            "profile_store_path": environment.get("profile_store_path"),
            "profile_gc_root": environment.get("profile_gc_root"),
            "store_paths": environment.get("store_paths"),
            "tool_versions": environment.get("tool_versions"),
            "wfmash_commit": WFMASH_COMMIT,
            "impg_commit": IMPG_COMMIT,
            "impg_executable_present": False,
            "forbidden_tools_used": [],
            "compute_smoke": {
                "status": smoke.get("status"),
                "slurm_job_id": smoke.get("slurm_job_id"),
                "compute_host": smoke.get("compute_host"),
                "compute_profile_store_path": smoke.get("compute_profile_store_path"),
                "store_path_identity_passed": smoke.get("store_path_identity_passed"),
                "wfmash_extended_cigar": smoke.get("wfmash_extended_cigar"),
                "bcf_csi_passed": smoke.get("bcf_csi_passed"),
                "callable_denominator_truth": smoke.get("callable_denominator_truth"),
            },
        },
        "pilot_gate": {
            "pilot_id": PILOT_ID,
            "registry_eligibility": pilot["eligibility"],
            "dual_modality_concordance": "not_run_no_approved_checksum_locked_tuple",
            "promotion_passed": False,
            "expansion_allowed": False,
        },
        "submitted_jobs": [],
        "resource_telemetry": "no_vgp_jobs_submitted_after_preflight_gate_failure",
        "result_table": {
            "path": "analysis/tier3a_data.tsv",
            "rows": 0,
            "header_columns": list(DATA_COLUMNS),
            "interpretation": "structured_missingness; no zero heterozygosity or denominator is implied",
        },
        "candidates": candidates,
    }
    _atomic_write(data_path, _tsv_bytes(DATA_COLUMNS, ()))
    _atomic_write(failure_ledger_path, _tsv_bytes(FAILURE_COLUMNS, failures))
    _atomic_write(qc_path, (json.dumps(qc_record, indent=2, sort_keys=True) + "\n").encode("utf-8"))
    return qc_record


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pilot-registry", type=Path, required=True)
    parser.add_argument("--environment", type=Path, required=True)
    parser.add_argument("--compute-smoke", type=Path, required=True)
    parser.add_argument("--data-output", type=Path, required=True)
    parser.add_argument("--failure-ledger", type=Path, required=True)
    parser.add_argument("--qc-output", type=Path, required=True)
    args = parser.parse_args(argv)
    finalize_ineligible_run(
        pilot_registry_path=args.pilot_registry,
        environment_path=args.environment,
        compute_smoke_path=args.compute_smoke,
        data_path=args.data_output,
        failure_ledger_path=args.failure_ledger,
        qc_path=args.qc_output,
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Tier3ValidationError as error:
        raise SystemExit(f"tier3a finalization rejected: {error}")
