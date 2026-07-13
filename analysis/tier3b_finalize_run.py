#!/usr/bin/env python3
"""Finalize a fail-closed Tier 3b run when no population tuple is eligible.

This is deliberately not a downloader or variant caller.  It converts the
frozen inventory decision into small, deterministic result/QC artifacts after
verifying the pinned Guix and compute-smoke records.  If either predeclared
pilot is promoted to eligibility, this path refuses to emit an empty table:
the approved tuple must instead be executed through the Tier 3 Slurm runner.
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

if __package__ in (None, ""):  # permit ``python analysis/<script>.py``
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from analysis.tier3_common import Tier3ValidationError, sha256_file


POLICY_ID = "tier3-decisions-v1"
SCHEMA_VERSION = "1.0"
RUN_ID = "tier3b-fail-closed-20260713"
PILOT_IDS = ("population_dgrp_freeze2", "population_ag1000g_phase3")

# This is an execution inventory, not a claim that the named survey resources
# are approved.  Each entry stays excluded until a separately frozen manifest
# supplies every checksum-addressed item listed in ``missing_inputs``.
CANDIDATES = (
    {
        "candidate_id": "dgrp_freeze2",
        "pilot_id": "population_dgrp_freeze2",
        "scientific_name": "Drosophila melanogaster",
        "taxon_id": 7227,
        "source_release": "DGRP Freeze 2.0",
        "intended_design": "inbred_lines_haploidized",
        "ploidy": "one_haploid_consensus_per_inbred_line",
        "inbreeding": "established_inbred_line_panel",
        "blocking_codes": (
            "missing_frozen_analysis_manifest",
            "missing_exact_reference_tuple",
            "missing_selected_line_callable_denominator",
            "missing_native_exact_reference_annotation_audit",
            "missing_independent_approved_pilot_baseline",
        ),
    },
    {
        "candidate_id": "ag1000g_phase3_gambiae",
        "pilot_id": "population_ag1000g_phase3",
        "scientific_name": "Anopheles gambiae",
        "taxon_id": 7165,
        "source_release": "Ag1000G phase 3",
        "intended_design": "wild_diploid",
        "ploidy": "diploid_autosomes",
        "inbreeding": "wild_outbred",
        "blocking_codes": (
            "missing_frozen_analysis_manifest",
            "missing_deposited_release_exact_reference_tuple",
            "missing_selected_population_callable_denominator",
            "missing_native_exact_reference_annotation_audit",
            "missing_independent_approved_pilot_baseline",
        ),
    },
    {
        "candidate_id": "drosophila_simulans_170",
        "pilot_id": None,
        "scientific_name": "Drosophila simulans",
        "taxon_id": 7240,
        "source_release": "170-line panel (Rogers et al. 2018)",
        "intended_design": "inbred_lines_haploidized",
        "ploidy": "one_haploid_consensus_per_inbred_line",
        "inbreeding": "inbred_line_panel",
        "blocking_codes": (
            "pilot_expansion_gate_not_passed",
            "missing_frozen_analysis_manifest",
            "missing_exact_reference_tuple",
            "missing_selected_line_callable_denominator",
            "missing_native_exact_reference_annotation_audit",
        ),
    },
    {
        "candidate_id": "drosophila_pseudoobscura_panel",
        "pilot_id": None,
        "scientific_name": "Drosophila pseudoobscura",
        "taxon_id": 7237,
        "source_release": "population genomics panel (2026 preprint)",
        "intended_design": "unfrozen_pending_stable_release",
        "ploidy": "unfrozen",
        "inbreeding": "unfrozen",
        "blocking_codes": (
            "pilot_expansion_gate_not_passed",
            "missing_frozen_analysis_manifest",
            "missing_stable_deposited_release",
            "missing_exact_reference_tuple",
            "missing_selected_population_callable_denominator",
            "missing_native_exact_reference_annotation_audit",
        ),
    },
    {
        "candidate_id": "aedes_aegypti_1206",
        "pilot_id": None,
        "scientific_name": "Aedes aegypti",
        "taxon_id": 7159,
        "source_release": "1206-genome global panel (Science 2025)",
        "intended_design": "wild_diploid",
        "ploidy": "diploid_autosomes",
        "inbreeding": "wild_outbred",
        "blocking_codes": (
            "pilot_expansion_gate_not_passed",
            "missing_frozen_analysis_manifest",
            "missing_deposited_release_exact_reference_tuple",
            "missing_single_locality_structure_qc",
            "missing_selected_population_callable_denominator",
            "missing_native_exact_reference_annotation_audit",
        ),
    },
    {
        "candidate_id": "daphnia_pulex_panels",
        "pilot_id": None,
        "scientific_name": "Daphnia pulex",
        "taxon_id": 6669,
        "source_release": "Lynch-lab population genomics panels",
        "intended_design": "cyclical_parthenogen_temporal_panel",
        "ploidy": "design_not_primary_population_pi_eligible",
        "inbreeding": "clonal_cyclical_parthenogen_design",
        "blocking_codes": (
            "primary_population_design_ineligible",
            "pilot_expansion_gate_not_passed",
            "missing_frozen_analysis_manifest",
            "missing_qualifying_exact_cohort_mask",
            "missing_native_exact_reference_annotation_audit",
        ),
    },
)

DATA_COLUMNS = (
    "dataset_id",
    "scientific_name",
    "taxon_id",
    "source_release",
    "population_id",
    "sample_design",
    "sampling_units",
    "nominal_chromosomes",
    "minimum_called_chromosomes",
    "reference_accession_version",
    "fasta_sha256",
    "annotation_provider",
    "annotation_release",
    "gff_sha256",
    "annotation_status",
    "genetic_code",
    "denominator_kind",
    "invariant_sites_explicit",
    "selected_sample_list_sha256",
    "downsampling_rule",
    "downsampling_seed_digest",
    "population_pi_numerator",
    "population_pi_denominator",
    "population_pi",
    "population_pi_ci_low",
    "population_pi_ci_high",
    "pi_S_numerator",
    "pi_S_denominator",
    "pi_S",
    "pi_W_numerator",
    "pi_W_denominator",
    "pi_W",
    "pi_S_over_pi_W",
    "pi_S_over_pi_W_ci_low",
    "pi_S_over_pi_W_ci_high",
    "pi_S_over_pi_W_reference_conditioned",
    "exclusion_counts_json",
    "guix_profile_store_path",
    "slurm_job_id",
)

FAILURE_COLUMNS = (
    "candidate_id",
    "scientific_name",
    "taxon_id",
    "source_release",
    "pilot_id",
    "pilot_reproduction_status",
    "expansion_status",
    "job_status",
    "slurm_job_id",
    "blocking_codes",
    "raw_read_calling",
    "decision_version",
)


def _read_json(path: Path, label: str) -> Dict[str, Any]:
    try:
        value = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, ValueError) as error:
        raise Tier3ValidationError("invalid {}: {}".format(label, error)) from error
    if not isinstance(value, dict):
        raise Tier3ValidationError("{} must be a JSON object".format(label))
    return value


def _atomic_write(path: Path, data: bytes) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(prefix="." + path.name + ".", dir=str(path.parent))
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
        raise Tier3ValidationError("pilot registry and Guix environment must use {}".format(POLICY_ID))

    pilots = {
        item.get("pilot_id"): item
        for item in registry.get("pilots", ())
        if isinstance(item, dict) and isinstance(item.get("pilot_id"), str)
    }
    missing = sorted(set(PILOT_IDS) - set(pilots))
    if missing:
        raise Tier3ValidationError("pilot registry is missing {!r}".format(missing))
    for pilot_id in PILOT_IDS:
        eligibility = pilots[pilot_id].get("eligibility")
        if not isinstance(eligibility, str) or not eligibility.startswith("ineligible_pending"):
            raise Tier3ValidationError(
                "{} is no longer pending and must be run, not finalized as unavailable".format(pilot_id)
            )

    profile = environment.get("profile_store_path")
    if not isinstance(profile, str) or not profile.startswith("/gnu/store/"):
        raise Tier3ValidationError("Guix environment has no valid profile store path")
    if smoke.get("status") != "passed":
        raise Tier3ValidationError("compute smoke must pass before Tier 3b finalization")
    if smoke.get("compute_profile_store_path") != profile or smoke.get("login_profile_store_path") != profile:
        raise Tier3ValidationError("compute/login smoke profile differs from the pinned Guix profile")
    if smoke.get("decision_version") != POLICY_ID:
        raise Tier3ValidationError("compute smoke decision version differs from frozen policy")
    if smoke.get("store_path_identity_passed") is not True:
        raise Tier3ValidationError("compute smoke did not prove store-path identity")
    return {"registry": registry, "environment": environment, "smoke": smoke, "pilots": pilots}


def _candidate_qc(candidate: Mapping[str, Any], pilots: Mapping[str, Mapping[str, Any]]) -> Dict[str, Any]:
    pilot_id = candidate["pilot_id"]
    pilot_eligibility = pilots[pilot_id]["eligibility"] if pilot_id else None
    design = str(candidate["intended_design"])
    minimum_called = 18 if design == "inbred_lines_haploidized" else 36 if design == "wild_diploid" else None
    nominal_chromosomes = 20 if design == "inbred_lines_haploidized" else 40 if design == "wild_diploid" else None
    return {
        "candidate_id": candidate["candidate_id"],
        "scientific_name": candidate["scientific_name"],
        "taxon_id": candidate["taxon_id"],
        "source_release": candidate["source_release"],
        "inclusion_status": "excluded_preflight",
        "blocking_codes": list(candidate["blocking_codes"]),
        "pilot": {
            "pilot_id": pilot_id,
            "registry_eligibility": pilot_eligibility,
            "reproduction_status": "not_run_no_approved_baseline" if pilot_id else "not_applicable_expansion_candidate",
        },
        "reference": {
            "assembly_accession_version": None,
            "fasta_sha256": None,
            "contig_dictionary_audit": "not_run_missing_frozen_tuple",
            "vcf_ref_allele_audit": "not_run_missing_frozen_tuple",
        },
        "annotation": {
            "primary_4d_status": "unavailable_not_audited",
            "provider": None,
            "release": None,
            "assembly_accession_version": None,
            "fasta_sha256": None,
            "gff_sha256": None,
            "sequence_region_contig_mapping_sha256": None,
            "genetic_code": None,
            "native_vs_projected": None,
            "sampled_cds_reconstruction": "not_run_missing_frozen_tuple",
            "projected_annotation_used_for_primary": False,
        },
        "population_and_samples": {
            "population_id": None,
            "population_frozen": False,
            "design": design,
            "ploidy": candidate["ploidy"],
            "inbreeding": candidate["inbreeding"],
            "minimum_usable_sampling_units": 20,
            "selected_sampling_units": 0,
            "nominal_chromosomes": nominal_chromosomes,
            "minimum_called_chromosomes_per_site": minimum_called,
            "selected_sample_list_sha256": None,
            "downsampling_rule": "sha256_dataset_population_sample_take_20",
            "downsampling_seed_or_rank_digest": None,
            "exclusion_list_sha256": None,
        },
        "filters_missingness_callability": {
            "provider_filters": "not_frozen",
            "missingness_threshold": "18_of_20_line_alleles" if minimum_called == 18 else "36_of_40_chromosomes" if minimum_called == 36 else "not_applicable_ineligible_design",
            "denominator_kind": None,
            "invariant_sites_explicit": False,
            "callable_mask_sha256": None,
            "callable_sites": None,
            "exclusion_counts": None,
        },
        "statistics": {
            "population_pi": "unavailable",
            "population_pi_numerator": None,
            "population_pi_denominator": None,
            "independent_pilot_check": "not_run_no_approved_input",
            "pi_S_over_pi_W": "unavailable",
            "pi_S_over_pi_W_reference_conditioned": True,
            "block_bootstrap": "not_run_statistic_unavailable",
        },
        "job": {
            "status": "not_submitted_gate_failure",
            "slurm_job_id": None,
            "resource_request": {"cpus": "2-8", "memory_gb": "32-64", "walltime_hours": "4-12", "scratch_gb": "50-200"},
            "resource_telemetry": "not_applicable_no_job_submitted",
        },
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
    """Write deterministic empty-result, failure-ledger, and QC artifacts."""

    audited = _audit_inputs(pilot_registry_path, environment_path, compute_smoke_path)
    environment = audited["environment"]
    smoke = audited["smoke"]
    candidates = [_candidate_qc(candidate, audited["pilots"]) for candidate in CANDIDATES]
    failures = []
    for candidate, qc in zip(CANDIDATES, candidates):
        failures.append(
            {
                "candidate_id": candidate["candidate_id"],
                "scientific_name": candidate["scientific_name"],
                "taxon_id": candidate["taxon_id"],
                "source_release": candidate["source_release"],
                "pilot_id": candidate["pilot_id"] or "",
                "pilot_reproduction_status": qc["pilot"]["reproduction_status"],
                "expansion_status": "blocked_pilots_not_reproduced_exactly",
                "job_status": "not_submitted_gate_failure",
                "slurm_job_id": "",
                "blocking_codes": ";".join(candidate["blocking_codes"]),
                "raw_read_calling": "not_launched",
                "decision_version": POLICY_ID,
            }
        )

    qc_record = {
        "schema_version": SCHEMA_VERSION,
        "decision_version": POLICY_ID,
        "run_id": RUN_ID,
        "as_of_environment_record": environment.get("recorded_at"),
        "overall_status": "no_eligible_population_tuples",
        "included_population_count": 0,
        "candidate_count": len(candidates),
        "input_audit": {
            "frozen_analysis_manifest": "absent",
            "pilot_registry_sha256": sha256_file(pilot_registry_path),
            "guix_environment_record_sha256": sha256_file(environment_path),
            "compute_smoke_record_sha256": sha256_file(compute_smoke_path),
            "raw_population_payloads_committed": False,
        },
        "pilot_gate": {
            "required_order": ["population_dgrp_freeze2", "population_ag1000g_phase3", "expansion"],
            "dgrp_reproduction": "not_run_no_approved_checksum_locked_tuple_or_baseline",
            "ag1000g_reproduction": "not_run_no_approved_checksum_locked_tuple_or_baseline",
            "reproduced_exactly": False,
            "expansion_allowed": False,
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
            "compute_smoke": {
                "status": smoke.get("status"),
                "slurm_job_id": smoke.get("slurm_job_id"),
                "compute_host": smoke.get("compute_host"),
                "compute_profile_store_path": smoke.get("compute_profile_store_path"),
                "store_path_identity_passed": smoke.get("store_path_identity_passed"),
                "callable_denominator_truth": smoke.get("callable_denominator_truth"),
            },
            "forbidden_environment_managers_used": [],
        },
        "submitted_jobs": [],
        "resource_telemetry": "no_population_jobs_submitted_after_preflight_gate_failure",
        "raw_read_workflow": "not_launched_and_not_present",
        "polarization_gate": {
            "requested": False,
            "passed": False,
            "reason": "deferred_by_tier3-decisions-v1",
            "sfs_b_output": "absent",
        },
        "result_table": {
            "path": "analysis/tier3b_data.tsv",
            "rows": 0,
            "header_columns": list(DATA_COLUMNS),
            "interpretation": "structured_missingness; no diversity estimate may be synthesized",
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
        raise SystemExit("tier3b finalization rejected: {}".format(error))
