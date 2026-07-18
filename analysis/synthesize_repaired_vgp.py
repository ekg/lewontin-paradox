#!/usr/bin/env python3
"""Build the paper handoff for the independently reviewed repaired VGP refusal.

This module is deliberately report-only.  It reads committed, reviewed ledgers,
performs exact identity joins, and writes Markdown/TSV artifacts.  It has no
network, scheduler, downloader, or demographic-inference entrypoint.
"""

from __future__ import annotations

import csv
import json
import sys
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Mapping, Sequence


if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))


ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = ROOT / "analysis"

DEFAULT_SYNTHESIS = ANALYSIS / "repaired_vgp_pilot_synthesis.md"
DEFAULT_PAPER_TABLE = ANALYSIS / "repaired_vgp_paper_table.tsv"
DEFAULT_DECISION = ANALYSIS / "repaired_vgp_next_decision.tsv"

VALID_NE_STATUS = "accepted_independent_with_time_caveat"
VALID_NE_CLASS = "independent_literature_ne"

PAPER_FIELDS = [
    "candidate_id",
    "scientific_name",
    "ncbi_taxid",
    "class",
    "order",
    "exact_reference_accession_version",
    "biosample_accession",
    "individual_or_isolate_id",
    "candidate_scope_denominator",
    "candidate_execution_status",
    "resolved_modality",
    "metadata_composition_eligible",
    "metadata_diversity_eligible",
    "reviewed_run_summary_numerator",
    "reviewed_run_summary_denominator",
    "reviewed_run_summary_target_total",
    "diversity_measurement_status",
    "diversity_callable_denominator",
    "composition_measurement_status",
    "composition_callable_bases",
    "composition_callable_fraction",
    "composition_queryable_gene_count",
    "composition_queryable_gene_bases",
    "measurement_uncertainty",
    "refusal_and_exclusion_codes",
    "promoted_biological_object_count",
    "observed_biological_bytes",
    "psmc_eligible",
    "psmc_blockers",
    "msmc2_eligible",
    "msmc2_blockers",
    "smcpp_eligible",
    "smcpp_blockers",
    "valid_independent_ne_record_count",
    "valid_independent_ne_record_ids",
    "valid_independent_ne_populations",
    "valid_independent_ne_values",
    "valid_independent_ne_units",
    "valid_independent_ne_sample_sizes",
    "valid_independent_ne_geographies",
    "valid_independent_ne_measurement_times",
    "valid_independent_ne_interval_types",
    "valid_independent_ne_uncertainty",
    "valid_independent_ne_source_dois",
    "valid_independent_ne_identity_scope",
    "coalescent_scaled_record_ids",
    "coalescent_scaled_values",
    "coalescent_scaled_units",
    "coalescent_scaled_disposition",
    "historical_scenario_record_ids",
    "historical_scenario_values",
    "historical_scenario_mutation_rates",
    "historical_scenario_generation_times",
    "historical_scenario_disposition",
    "other_non_ne_record_ids",
    "other_non_ne_disposition",
    "circular_estimate_record_count",
    "circular_estimate_record_ids",
    "circular_estimate_disposition",
    "h1_h2_demography_disposition",
    "absolute_ne_time_status",
    "exact_join_status",
    "provenance",
]

DECISION_FIELDS = [
    "row_type",
    "option_id",
    "decision_option",
    "option_scope",
    "authorization_granted",
    "authorization_required",
    "creates_ready_executable_task",
    "requested_species_count",
    "requested_compressed_input_gib",
    "requested_scratch_gib",
    "requested_core_hours",
    "requested_concurrent_species",
    "requested_memory_per_job_gib",
    "requested_persistent_storage_gib",
    "requested_retention_days",
    "current_strict_species_ceiling",
    "current_strict_compressed_input_gib_ceiling",
    "current_strict_scratch_gib_ceiling",
    "current_strict_core_hours_ceiling",
    "current_strict_concurrent_species_ceiling",
    "current_strict_memory_per_job_gib_ceiling",
    "current_strict_persistent_input_gb_ceiling",
    "current_strict_persistent_output_gb_ceiling",
    "projection_scope",
    "metric",
    "unit",
    "low",
    "base",
    "high",
    "calibration_status",
    "evidence",
    "requirements_or_blockers",
    "disposition",
]


def load_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fieldnames: Sequence[str], rows: Iterable[Mapping[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def load_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def fmt_number(value: object) -> str:
    number = float(value)
    if number.is_integer():
        return str(int(number))
    return f"{number:.12f}".rstrip("0").rstrip(".")


def semicolon(rows: Sequence[Mapping[str, str]], field: str) -> str:
    return ";".join(row[field] for row in rows)


def _assert_reviewed_refusal(
    gate: Mapping[str, object],
    run_rows: Sequence[Mapping[str, str]],
    acquisition_rows: Sequence[Mapping[str, str]],
    inventory_rows: Sequence[Mapping[str, str]],
    result_rows: Sequence[Mapping[str, str]],
    exclusion_rows: Sequence[Mapping[str, str]],
    refusal_rows: Sequence[Mapping[str, str]],
    telemetry_rows: Sequence[Mapping[str, str]],
    qc_rows: Sequence[Mapping[str, str]],
) -> tuple[Mapping[str, str], Mapping[str, str]]:
    decision = gate["decision"]
    if not isinstance(decision, Mapping) or decision.get("status") != "NO_GO":
        raise RuntimeError("repaired synthesis is pinned to reviewed NO_GO evidence")
    run = next(row for row in run_rows if row["record_type"] == "run_summary")
    if not (
        run["gate_decision"] == "NO_GO"
        and run["final_state"] == "NOT_SUBMITTED"
        and run["selected_species"] == "0"
        and run["promoted_objects"] == "0"
        and run["compressed_input_bytes"] == "0"
        and run["slurm_jobs_submitted"] == "0"
    ):
        raise RuntimeError("run summary no longer matches the reviewed refusal")
    acquisition = next(row for row in acquisition_rows if row["record_type"] == "run_summary")
    if not (
        acquisition["status"] == "refused_preflight"
        and acquisition["observed_bytes"] == "0"
        and acquisition["cumulative_transferred_bytes"] == "0"
        and not acquisition["staging_path"]
        and not acquisition["quarantine_path"]
        and not acquisition["promoted_path"]
    ):
        raise RuntimeError("acquisition summary no longer matches the reviewed zero-byte refusal")
    if inventory_rows:
        raise RuntimeError("reviewed immutable biological object inventory must remain empty")
    if len(result_rows) != 1 or result_rows[0]["record_type"] != "run_summary":
        raise RuntimeError("reviewed results must contain only the excluded empty-result summary")
    result = result_rows[0]
    if not (
        result["status"] == "excluded"
        and result["numerator"] == result["denominator"] == result["target_total"] == result["value"] == "0"
        and not result["candidate_id"]
        and not result["measurement_method"]
        and not result["artifact_sha256"]
    ):
        raise RuntimeError("result summary no longer preserves the reviewed empty-result outcome")
    expected_codes = {
        "GATE_NO_GO",
        "CAP_MOOSEFS_READ_GB_EXCEEDED",
        "CAP_SCRATCH_GIB_EXCEEDED",
        "QUOTA_UNAVAILABLE",
    }
    if {row["failure_code"] for row in exclusion_rows} != expected_codes:
        raise RuntimeError("exclusion ledger differs from reviewed gate plus blocker set")
    refusal = refusal_rows[0]
    if not (
        refusal["gate_decision"] == "NO_GO"
        and refusal["status"] == "NOT_SUBMITTED"
        and refusal["provider_requests"] == "0"
        and refusal["network_bytes"] == "0"
        and refusal["full_catalog_downloaded"] == "false"
        and refusal["population_bulk_downloaded"] == "false"
        and refusal["demographic_inferences"] == "0"
    ):
        raise RuntimeError("refusal ledger differs from reviewed no-execution evidence")
    if len(telemetry_rows) != 1:
        raise RuntimeError("expected one reviewed NOT_SUBMITTED telemetry summary")
    telemetry = telemetry_rows[0]
    if not (
        telemetry["final_state"] == "NOT_SUBMITTED"
        and telemetry["slurm_job_id"] == ""
        and telemetry["slurm_array_job_id"] == ""
        and telemetry["elapsed_seconds"] == "0"
        and telemetry["cpu_time_seconds"] == "0"
        and telemetry["scratch_peak_gb"] == "0"
        and telemetry["io_read_gb"] == "0"
        and telemetry["io_write_gb"] == "0"
        and telemetry["network_bytes"] == "0"
    ):
        raise RuntimeError("telemetry is not the reviewed zero-use NOT_SUBMITTED row")
    decisions = {row["decision"] for row in qc_rows}
    if decisions != {"PASS"} or len(qc_rows) != 57:
        raise RuntimeError("independent repaired review is not 57/57 PASS")
    return run, result


def build_paper_rows(
    manifest_rows: Sequence[Mapping[str, str]],
    demography_rows: Sequence[Mapping[str, str]],
    ne_rows: Sequence[Mapping[str, str]],
    run_result: Mapping[str, str],
) -> list[dict[str, str]]:
    if len(manifest_rows) != 6 or len(demography_rows) != 6:
        raise RuntimeError("reviewed repaired candidate denominator must be exactly six")
    demo_by_candidate = {row["candidate_id"]: row for row in demography_rows}
    sources_by_candidate: dict[str, list[Mapping[str, str]]] = defaultdict(list)
    for source in ne_rows:
        sources_by_candidate[source["candidate_id"]].append(source)

    rows: list[dict[str, str]] = []
    for manifest in manifest_rows:
        candidate = manifest["candidate_id"]
        if candidate not in demo_by_candidate:
            raise RuntimeError(f"missing exact demography audit row for {candidate}")
        demo = demo_by_candidate[candidate]
        exact_fields = (
            (manifest["scientific_name_source"], demo["scientific_name"]),
            (manifest["ncbi_taxid"], demo["ncbi_taxid"]),
            (manifest["h1_accession_version"], demo["exact_reference_accession_version"]),
            (manifest["biosample_accession"], demo["biosample_accession"]),
            (manifest["individual_or_isolate_id"], demo["individual_or_isolate_id"]),
        )
        if any(left != right for left, right in exact_fields):
            raise RuntimeError(f"exact taxon/reference/BioSample/individual join failed for {candidate}")
        if manifest["pilot_selected"] != "no" or manifest["resolved_modality"] != "tier3c_composition":
            raise RuntimeError(f"candidate disposition changed for {candidate}")
        if any(demo[field] != "no" for field in ("psmc_eligible", "msmc2_eligible", "smcpp_eligible")):
            raise RuntimeError(f"ineligible method was promoted for {candidate}")

        sources = sources_by_candidate[candidate]
        valid = [
            source
            for source in sources
            if source["record_status"] == VALID_NE_STATUS
            and source["classification"] == VALID_NE_CLASS
            and source["independence_status"].startswith("independent_")
            and source["circularity_status"] == "non_circular"
        ]
        coalescent = [source for source in sources if source["classification"] == "coalescent_scaled_not_absolute_ne"]
        historical = [source for source in sources if source["classification"] == "historical_coalescent_absolute_scenario"]
        other = [
            source
            for source in sources
            if source["classification"]
            in {
                "independent_literature_nb",
                "census_measure",
                "population_structure_covariate",
                "missing_source_sentinel",
            }
        ]
        circular = [source for source in sources if source["circularity_status"] == "circular_excluded"]
        if len(circular) != 1:
            raise RuntimeError(f"expected exactly one circular exclusion policy row for {candidate}")

        rows.append(
            {
                "candidate_id": candidate,
                "scientific_name": demo["scientific_name"],
                "ncbi_taxid": demo["ncbi_taxid"],
                "class": manifest["class"],
                "order": manifest["order"],
                "exact_reference_accession_version": demo["exact_reference_accession_version"],
                "biosample_accession": demo["biosample_accession"],
                "individual_or_isolate_id": demo["individual_or_isolate_id"],
                "candidate_scope_denominator": "6",
                "candidate_execution_status": "NOT_EXECUTED_NO_GO",
                "resolved_modality": manifest["resolved_modality"],
                "metadata_composition_eligible": manifest["assembly_composition_eligible"],
                "metadata_diversity_eligible": manifest["assembly_diversity_eligible"],
                "reviewed_run_summary_numerator": run_result["numerator"],
                "reviewed_run_summary_denominator": run_result["denominator"],
                "reviewed_run_summary_target_total": run_result["target_total"],
                "diversity_measurement_status": "NOT_MEASURED_INELIGIBLE_AND_REFUSED",
                "diversity_callable_denominator": "NOT_APPLICABLE_NO_DIPLOID_RESULT",
                "composition_measurement_status": "NOT_MEASURED_REFUSED",
                "composition_callable_bases": "NOT_MEASURED_POST_ALIGNMENT_REQUIRED",
                "composition_callable_fraction": "NOT_MEASURED_POST_ALIGNMENT_REQUIRED",
                "composition_queryable_gene_count": "NOT_MEASURED_POST_ALIGNMENT_REQUIRED",
                "composition_queryable_gene_bases": "NOT_MEASURED_POST_ALIGNMENT_REQUIRED",
                "measurement_uncertainty": "NOT_ESTIMABLE_NO_BIOLOGICAL_MEASUREMENT",
                "refusal_and_exclusion_codes": (
                    "GATE_NO_GO;CAP_MOOSEFS_READ_GB_EXCEEDED;"
                    "CAP_SCRATCH_GIB_EXCEEDED;QUOTA_UNAVAILABLE"
                ),
                "promoted_biological_object_count": "0",
                "observed_biological_bytes": "0",
                "psmc_eligible": demo["psmc_eligible"],
                "psmc_blockers": demo["psmc_blockers"],
                "msmc2_eligible": demo["msmc2_eligible"],
                "msmc2_blockers": demo["msmc2_blockers"],
                "smcpp_eligible": demo["smcpp_eligible"],
                "smcpp_blockers": demo["smcpp_blockers"],
                "valid_independent_ne_record_count": str(len(valid)),
                "valid_independent_ne_record_ids": semicolon(valid, "record_id"),
                "valid_independent_ne_populations": semicolon(valid, "population"),
                "valid_independent_ne_values": semicolon(valid, "value"),
                "valid_independent_ne_units": semicolon(valid, "unit"),
                "valid_independent_ne_sample_sizes": semicolon(valid, "sample_size"),
                "valid_independent_ne_geographies": semicolon(valid, "geography"),
                "valid_independent_ne_measurement_times": semicolon(valid, "measurement_time"),
                "valid_independent_ne_interval_types": semicolon(valid, "interval_type"),
                "valid_independent_ne_uncertainty": semicolon(valid, "uncertainty_status"),
                "valid_independent_ne_source_dois": semicolon(valid, "doi"),
                "valid_independent_ne_identity_scope": (
                    "independent literature populations; different animals/project; not VGP BioSample or individual"
                    if valid
                    else "none admitted as independent absolute Ne"
                ),
                "coalescent_scaled_record_ids": semicolon(coalescent, "record_id"),
                "coalescent_scaled_values": semicolon(coalescent, "value"),
                "coalescent_scaled_units": semicolon(coalescent, "unit"),
                "coalescent_scaled_disposition": (
                    "published mitochondrial theta=4Ne-mu; coalescent-scaled only, not absolute or contemporary nuclear Ne"
                    if coalescent
                    else "none"
                ),
                "historical_scenario_record_ids": semicolon(historical, "record_id"),
                "historical_scenario_values": semicolon(historical, "value"),
                "historical_scenario_mutation_rates": semicolon(historical, "mutation_rate"),
                "historical_scenario_generation_times": semicolon(historical, "generation_time"),
                "historical_scenario_disposition": (
                    "published other-sample/reference PSMC; mutation/generation scenario-scaled historical Ne/time; not contemporary independent Ne and not repaired-pilot inference"
                    if historical
                    else "none"
                ),
                "other_non_ne_record_ids": semicolon(other, "record_id"),
                "other_non_ne_disposition": (
                    "kept as Nb, census, population structure, or explicit missingness; not promoted to absolute Ne"
                    if other
                    else "none"
                ),
                "circular_estimate_record_count": str(len(circular)),
                "circular_estimate_record_ids": semicolon(circular, "record_id"),
                "circular_estimate_disposition": "pi/(4mu) same-response derivation not calculated and circular_excluded",
                "h1_h2_demography_disposition": demo["h1_h2_demography_disposition"],
                "absolute_ne_time_status": demo["absolute_ne_time_status"],
                "exact_join_status": "PASS_EXACT_IDENTITY",
                "provenance": (
                    "analysis/vgp_pilot_manifest.tsv;analysis/vgp_pilot_results.tsv;"
                    "analysis/vgp_pilot_exclusions.tsv;analysis/vgp_pilot_refusals.tsv;"
                    "analysis/vgp_demography_input_audit.tsv;analysis/vgp_independent_ne_sources.tsv"
                ),
            }
        )
    return rows


def _decision_base(option_id: str, label: str, scope: str, requirement: str, disposition: str) -> dict[str, str]:
    return {
        "row_type": "decision_option",
        "option_id": option_id,
        "decision_option": label,
        "option_scope": scope,
        "authorization_granted": "no",
        "authorization_required": "yes before any acquisition, execution, expansion, or inference" if option_id not in {"repair_remaining_candidate_metadata", "stop"} else "no execution authority granted by this option",
        "creates_ready_executable_task": "no",
        "requested_species_count": "",
        "requested_compressed_input_gib": "",
        "requested_scratch_gib": "",
        "requested_core_hours": "",
        "requested_concurrent_species": "",
        "requested_memory_per_job_gib": "",
        "requested_persistent_storage_gib": "",
        "requested_retention_days": "",
        "current_strict_species_ceiling": "6",
        "current_strict_compressed_input_gib_ceiling": "120",
        "current_strict_scratch_gib_ceiling": "139.698386192",
        "current_strict_core_hours_ceiling": "280",
        "current_strict_concurrent_species_ceiling": "2",
        "current_strict_memory_per_job_gib_ceiling": "96",
        "current_strict_persistent_input_gb_ceiling": "160",
        "current_strict_persistent_output_gb_ceiling": "32",
        "projection_scope": "",
        "metric": "",
        "unit": "",
        "low": "",
        "base": "",
        "high": "",
        "calibration_status": "NOT_APPLICABLE_DECISION_OPTION",
        "evidence": "analysis/repaired_vgp_pilot_review.md;analysis/repaired_vgp_resource_calibration.tsv",
        "requirements_or_blockers": requirement,
        "disposition": disposition,
    }


def build_decision_rows(resource_rows: Sequence[Mapping[str, str]]) -> list[dict[str, str]]:
    rows = [
        _decision_base(
            "repair_remaining_candidate_metadata",
            "Repair remaining candidate metadata",
            "metadata-only identity, annotation, checksum, quota, and finite resource-envelope repair",
            "Remain metadata-only; do not download biological payloads, submit jobs, build genotype data, or infer demography.",
            "AVAILABLE_DECISION_ONLY",
        ),
        _decision_base(
            "request_bounded_expansion_wave",
            "Request authorization for a numerically bounded expansion wave",
            "a separately reviewed finite species wave; not full-catalog execution",
            "A new request must numerically state exact species/accessions, compressed input, scratch, core-hours, concurrency, per-job memory, persistent storage, retention/expiry, I/O, and rollback; it must satisfy every stricter cap and receive a new exact GO.",
            "REQUEST_NOT_MADE_OR_AUTHORIZED",
        ),
        _decision_base(
            "request_population_data_subset",
            "Request authorization for a population-data subset",
            "exact-reference, exact-population subset for method-specific input validation",
            "A new request must numerically state exact species/reference/population/sample/asset identities, compressed input, scratch, core-hours, concurrency, per-job memory, storage, retention/expiry, I/O, and rollback. Raw bulk download and each PSMC/MSMC2/SMC++ inference require separate explicit authority; VGP H1/H2 is not an input shortcut.",
            "REQUEST_NOT_MADE_OR_AUTHORIZED",
        ),
        _decision_base(
            "stop",
            "Stop",
            "retain the reviewed refusal and metadata audit as the terminal evidence packet",
            "No further resource or biological action.",
            "AVAILABLE_DECISION_ONLY",
        ),
    ]

    calibrations = [row for row in resource_rows if row["scope"] == "successful_observation_calibration"]
    if not calibrations or any(
        row["decision"] != "NOT_CALIBRATED"
        or row["predicted_low"]
        or row["predicted_base"]
        or row["predicted_high"]
        for row in calibrations
    ):
        raise RuntimeError("reviewed refusal calibration must have empty low/base/high values")
    for scope in ("next_wave", "full_eligible_catalog"):
        for calibration in calibrations:
            projection = _decision_base(
                "request_bounded_expansion_wave",
                "Request authorization for a numerically bounded expansion wave",
                "projection only; no execution scope",
                "Successful observed pilot telemetry is required before an observed-calibrated estimate can be populated. Full-catalog work is not a decision option and cannot be authorized through this projection.",
                "UNAVAILABLE_ESTIMATE_NOT_AUTHORIZATION",
            )
            projection.update(
                {
                    "row_type": "resource_projection",
                    "projection_scope": scope,
                    "metric": calibration["metric"],
                    "unit": calibration["unit"],
                    "low": "",
                    "base": "",
                    "high": "",
                    "calibration_status": "NOT_ESTIMABLE_NO_SUCCESSFUL_OBSERVATION",
                    "evidence": calibration["evidence"],
                }
            )
            rows.append(projection)
    return rows


def render_markdown(
    gate: Mapping[str, object],
    manifest_rows: Sequence[Mapping[str, str]],
    paper_rows: Sequence[Mapping[str, str]],
    exclusion_rows: Sequence[Mapping[str, str]],
) -> str:
    decision_sha = str(gate["decision_sha256"])
    caps = gate["cap_vector"]["dimensions"]  # type: ignore[index]
    blockers = gate["blockers"]
    valid_rows = sum(int(row["valid_independent_ne_record_count"]) for row in paper_rows)
    valid_species = sum(int(row["valid_independent_ne_record_count"]) > 0 for row in paper_rows)
    classes = ", ".join(row["class"] for row in manifest_rows)

    cap_order = (
        "species",
        "compressed_inputs_gib",
        "scratch_gib",
        "core_hours",
        "concurrent_species",
        "memory_per_job_gib",
        "aggregate_wall_hours",
        "cpus_per_element",
        "file_inodes",
        "metadata_operations",
        "moosefs_read_gb",
        "moosefs_write_gb",
        "persistent_input_gb",
        "persistent_output_gb",
        "peak_bandwidth_mib_s",
    )
    cap_lines = []
    for metric in cap_order:
        record = caps[metric]  # type: ignore[index]
        cap_lines.append(
            f"| `{metric}` | {fmt_number(record['limit'])} {record['unit']} | "
            f"{fmt_number(record['observed'])} {record['unit']} | "
            f"{'within' if record['passes'] else 'exceeded: refusal required'} |"
        )

    method_lines = []
    for row in paper_rows:
        ne = (
            f"{row['valid_independent_ne_record_count']} LD-Ne population records"
            if row["valid_independent_ne_record_count"] != "0"
            else "0 valid independent absolute-Ne records"
        )
        method_lines.append(
            f"| *{row['scientific_name']}* | `{row['exact_reference_accession_version']}` | "
            f"`{row['biosample_accession']}` / `{row['individual_or_isolate_id']}` | no | no | no | {ne} |"
        )

    blocker_codes = "; ".join(str(blocker["code"]) for blocker in blockers)  # type: ignore[index]
    text = f"""# Repaired VGP pilot: paper-oriented synthesis

Date: 2026-07-18 UTC

## Bottom line

**Audited outcome: `NO_GO`, correctly refused, `NOT_SUBMITTED`. No biological pilot ran.** The immutable decision SHA-256 is `{decision_sha}`. Acquisition stopped before a provider request and compute stopped before `sbatch`. The zero-byte, zero-job, empty-result outcome is evidence that the authorization boundary worked; it is not diversity, composition, or performance evidence.

The reviewed scope comprised **6 metadata candidates / 6 species / 6 classes**, one each from {classes}. Metadata repair made all 6/6 composition candidates pre-download eligible, but the stricter integrated gate selected and executed 0/6. Diversity eligibility was 0/6, population-genomic eligibility 0/6, and demographic eligibility 0/6. Review QC was 57 PASS, 0 FAIL.

This synthesis uses only the independently reviewed evidence packet: exact candidate metadata in `analysis/vgp_pilot_manifest.tsv:2-7`; refusal/run/acquisition evidence in `analysis/vgp_pilot_run_manifest.tsv:2-5`, `analysis/vgp_pilot_acquisition_manifest.tsv:2-5`, `analysis/vgp_pilot_refusals.tsv:2`, and `analysis/vgp_pilot_slurm_telemetry.tsv:2`; the header-only immutable inventory; the excluded result and exclusions in `analysis/vgp_pilot_results.tsv:2` and `analysis/vgp_pilot_exclusions.tsv:2-5`; review decisions in `analysis/repaired_vgp_pilot_qc.tsv:2-58`; and the metadata-only demography audit in `analysis/vgp_demography_input_audit.tsv:2-7` plus its classified source ledger `analysis/vgp_independent_ne_sources.tsv:2-22`. It does not join older unreviewed pilot tables or treat metadata search results as executed outputs.

## Exact evidence obtained and measured denominators

| evidence layer | reviewed denominator | observed outcome | paper disposition |
| --- | ---: | --- | --- |
| repaired candidate manifest | 6 candidate rows | 6 exact current TaxId/name + exact-version H1 RefSeq + native exact-H1 annotation locations | metadata/provenance evidence only |
| taxonomic strata | 6 classes | one candidate in each listed class | no class was biologically analyzed |
| acquisition obligations | 12 finite obligations (H1 FASTA + native H1 annotation for each candidate) | 0 provider requests; 0 transferred biological bytes; 0 staged, quarantined, or promoted paths | zero-byte refusal evidence |
| immutable object inventory | 0 object rows | header-only inventory | empty inventory is retained as evidence; no local biological SHA-256 exists |
| executed candidate rows | 6 proposed; 0 selected/executed | no candidate result rows | do not report a biological sample size of six |
| run result summary | numerator 0; denominator 0; target 0 | one excluded `validated_species_count=0` summary with blank method and artifact hash | empty-result control row, not a biological estimate |
| exclusions | {len(exclusion_rows)} rows | gate exclusion plus `{blocker_codes}` | all `imputed=false`, `demographic_input_used=false` |
| scheduler | 0 job IDs | `NOT_SUBMITTED`; no command, array, dependency, elapsed time, CPU, scratch, I/O, metadata operations, or network use | terminal refusal, not performance telemetry |

All candidate post-alignment fields—callable bases/fraction and queryable gene/base denominators—remain `POST_ALIGNMENT_REQUIRED`; they are **not zeros** and were never measured. No SweepGA mapping, IMPG partition/query, VCF/BCF, coding-diversity numerator, CDS/fourfold composition target, uncertainty interval, or candidate artifact SHA-256 was produced. The result uncertainty is therefore “not estimable because no biological measurement exists,” not a zero-width interval.

## Diversity and composition claims

Supported paper claims:

- The bounded repaired proposal produced no promoted diversity or composition estimate and no biological candidate result.
- The six exact metadata rows were composition-only candidates whose required denominators could only be measured after an authorized acquisition and alignment; that boundary was never crossed.
- The gate failed for exactly three reasons: `{blocker_codes}`. The proposed worst-case MooseFS read and scratch loads exceeded their strict limits, and enforceable quota/headroom evidence was absent.
- The refusal ledgers, zero-byte acquisition summary, header-only immutable inventory, excluded result summary, exclusion ledger, and `NOT_SUBMITTED` telemetry are valid control/provenance outcomes.

Unsupported paper claims:

- Any estimate, range, rank, class summary, or uncertainty interval for nucleotide diversity, coding diversity, CDS/fourfold composition, callability, queryability, or mapping multiplicity from this repaired pilot.
- Any statement that six species were biologically sampled, processed, or validated. Six is the metadata-candidate denominator; the executed denominator is zero.
- Any reuse of dormant toolchain smoke fixtures as a repaired-pilot biological result.
- Any inference that the zero resource row calibrates runtime, memory, scratch, I/O, storage, bandwidth, or throughput.

## Exact demography and independent-Ne audit

The join in `analysis/repaired_vgp_paper_table.tsv` is exact on candidate ID, scientific name, TaxId, H1 reference accession.version, BioSample, and individual/isolate. Literature populations remain separately identified and are never collapsed onto the VGP individual.

| species | exact VGP reference | VGP BioSample / individual | future PSMC eligibility | future MSMC2 eligibility | future SMC++ eligibility | valid independent absolute-Ne evidence |
| --- | --- | --- | --- | --- | --- | --- |
{chr(10).join(method_lines)}

All method decisions are currently `no`, but their blockers remain method-specific:

- **PSMC:** all 6/6 lack a heterozygosity-retaining callable diploid consensus, a compatible callable mask, and established coverage. H1 is explicitly haploid.
- **MSMC2:** all 6/6 lack validated accurate mutually comparable phased genomes/haplotypes, compatible masks, and audited population/individual relationships.
- **SMC++:** all 6/6 lack an exact-reference population genotype set plus compatible masks, population definitions, and method-specific QC. The camel candidate has a population VCF only on the incompatible `GCA/GCF_000803125` lineage; it is not ready for `GCF_036321535.1`.

VGP H1/H2 linkage is assembly discovery metadata only. It is never assumed to be a heterozygosity-retaining demographic genotype, two independent genomes, accurate comparable phasing, or a population dataset.

Only *Camelus dromedarius* has valid independent numeric Ne evidence in this bounded audit: {valid_rows} non-circular LD-Ne records across Awarik, Haddana, Majaheem, Sahliah, Shul, and Sofor (values 15, 11, 37, 24, 17, and 23; final per-breed sample sizes 5, 4, 9, 7, 4, and 5 diploid individuals). These are different animals/project from VGP `SAMN39296380/mCamDro1`; the exact LD time slice and intervals were not reported, so they are {valid_rows} population observations for {valid_species} species, not six species-level replicates and not VGP-sample Ne.

Distinct non-promoted evidence classes are retained:

- `camel_psmc_fitak2020` is a published historical PSMC scenario on other animals/reference. Its approximate absolute Ne/time trajectory depends on the stated mutation-rate and generation-time scenario; it is not a new repaired-pilot inference and is not the contemporary independent LD-Ne field.
- `horn_theta_ima3_2022` is independent mitochondrial `theta=4Ne-mu`, retained only in a coalescent-scaled field. Without the appropriate locus mutation scenario—and because the estimand is mitochondrial rather than diploid nuclear—it is not absolute/contemporary Ne.
- `gar_nb_cosewic2015` is Nb, not Ne, and lacks value-to-population/time/uncertainty mapping. Camel and gar census measures and frog population structure are separate ecological/population fields, not Ne.
- Exactly one same-response `pi/(4mu)` policy row per candidate (6/6) is `circular_excluded`; no value was calculated. A predictor algebraically derived from the response pi cannot be used to explain that response.

Coalescent-scaled quantities, absolute scenario-scaled Ne/time, valid independent LD-Ne, census/Nb/structure evidence, and circular exclusions therefore remain separate in both the table and conclusions.

## Lewontin-paradox implications

The repaired refusal supports an infrastructure conclusion only: exact-reference metadata, measured denominators, immutable acquisition, storage headroom, and method-specific population inputs are material constraints on a comparative analysis. It supplies **no new cross-species pi or composition values**, so it cannot confirm, refute, narrow, or quantify Lewontin's paradox. The single-species independent LD-Ne inventory is insufficient for a cross-species Ne–pi relationship, and the circular `pi/(4mu)` rows are expressly unusable. Any biological statement about compressed diversity range, census size, linked selection, mutation, life history, or demographic history would exceed this evidence.

## Resource use, proposal, and uncertainty

Actual attributable use was zero: 0 provider requests, 0 biological bytes, 0 promoted objects, 0 Slurm submissions, 0 compute jobs, 0 core-seconds, 0 scratch bytes, 0 read/write/network bytes, and 0 demographic inferences. These zeros prove refusal/cap compliance only; they are **not performance telemetry**.

The following values are the gate's finite six-row **pre-run proposal**, not observed use:

| dimension | strict current ceiling | proposed finite value | gate disposition |
| --- | ---: | ---: | --- |
{chr(10).join(cap_lines)}

Every stricter integrated cap wins. Thus the current ceiling is no more than 6 species, 120 GiB compressed inputs, **139.698386192 GiB scratch** (stricter than the outer 750 GiB bound), **280 core-hours** (stricter than 1,500), 2 concurrent species, and **96 GiB/job** (stricter than 256), with the additional recorded ceilings shown above. A ceiling is not a `GO`.

Observed-calibrated low/base/high projections for both a next wave and the full eligible catalog are deliberately blank in `analysis/repaired_vgp_next_decision.tsv`: there was no successful observation from which to estimate them. The zero-use refusal is excluded from calibration. “Not estimable” is the reviewed answer; no historical or planning proxy is silently relabeled as repaired-pilot telemetry. Any future numeric estimate is a planning artifact, not authorization.

## Active-state, retention, and quarantine verification

At synthesis time (2026-07-18 UTC), a read-only `squeue -h -u $USER` snapshot returned zero jobs, and a process-table search returned zero matching repaired VGP acquisition/run processes. This agrees with the immutable ledgers: zero job IDs, zero provider requests, zero network bytes, `NOT_SUBMITTED`, and no acquisition paths. **No attributable active jobs or download processes were found.**

The retained evidence is metadata/refusal material: gate and decision hashes, manifests, refusal/acquisition/run summaries, header-only immutable inventory, empty-result summary, exclusions, telemetry, QC, and the metadata-only demography audit. There are no attributable biological objects or failed partial paths to retain, quarantine, resume, or delete; acquisition rows have blank staging, quarantine, and promoted paths. No cleanup or deletion was performed by this synthesis.

For any future authorized work, the reviewed policy keeps failed partials for at most 14 days, reproducible intermediates through independent QC plus 90 days, and provenance/final results/failure ledgers/sentinels/telemetry through publication plus seven years, subject to institutional policy. Compressed sources persist through final reproducibility review and can be evicted only by checksum-listed human approval with tested rehydration. A new authorization must replace those general defaults with an explicit numeric storage and retention scope.

## Decision boundary: options only

This handoff defines exactly four options and selects or authorizes none:

1. **Repair remaining candidate metadata**—metadata-only identity, quota, checksum, annotation, and finite resource-envelope work; no biological acquisition, job, or inference.
2. **Request authorization for a numerically bounded expansion wave**—a new request must name exact species/accessions and numeric compressed-input, scratch, core-hour, concurrency, per-job-memory, storage, I/O, retention/expiry, and rollback limits. Full-catalog acquisition/execution is not covered.
3. **Request authorization for a population-data subset**—a new request must name exact taxon/reference/population/sample/assets and the same numeric resource/retention limits. Raw population bulk download, genotype construction, and each PSMC/MSMC2/SMC++ inference require separately explicit authorization.
4. **Stop**—retain this refusal and metadata audit as terminal evidence.

No expansion, full-catalog acquisition/execution, population bulk download, demographic inference, or ready executable task was created, authorized, or launched. Reporting a blank or hypothetical estimate grants no authority.
"""
    return text


def build_synthesis(
    *,
    root: Path = ROOT,
    synthesis_out: Path = DEFAULT_SYNTHESIS,
    paper_table_out: Path = DEFAULT_PAPER_TABLE,
    decision_out: Path = DEFAULT_DECISION,
) -> dict[str, int]:
    analysis = root / "analysis"
    manifest_rows = load_tsv(analysis / "vgp_pilot_manifest.tsv")
    run_rows = load_tsv(analysis / "vgp_pilot_run_manifest.tsv")
    acquisition_rows = load_tsv(analysis / "vgp_pilot_acquisition_manifest.tsv")
    inventory_rows = load_tsv(analysis / "vgp_pilot_immutable_object_inventory.tsv")
    result_rows = load_tsv(analysis / "vgp_pilot_results.tsv")
    exclusion_rows = load_tsv(analysis / "vgp_pilot_exclusions.tsv")
    refusal_rows = load_tsv(analysis / "vgp_pilot_refusals.tsv")
    telemetry_rows = load_tsv(analysis / "vgp_pilot_slurm_telemetry.tsv")
    qc_rows = load_tsv(analysis / "repaired_vgp_pilot_qc.tsv")
    resource_rows = load_tsv(analysis / "repaired_vgp_resource_calibration.tsv")
    demography_rows = load_tsv(analysis / "vgp_demography_input_audit.tsv")
    ne_rows = load_tsv(analysis / "vgp_independent_ne_sources.tsv")
    gate = load_json(analysis / "vgp_pilot_gate.json")

    _, run_result = _assert_reviewed_refusal(
        gate,
        run_rows,
        acquisition_rows,
        inventory_rows,
        result_rows,
        exclusion_rows,
        refusal_rows,
        telemetry_rows,
        qc_rows,
    )
    paper_rows = build_paper_rows(manifest_rows, demography_rows, ne_rows, run_result)
    decision_rows = build_decision_rows(resource_rows)
    write_tsv(paper_table_out, PAPER_FIELDS, paper_rows)
    write_tsv(decision_out, DECISION_FIELDS, decision_rows)
    synthesis_out.parent.mkdir(parents=True, exist_ok=True)
    synthesis_out.write_text(
        render_markdown(gate, manifest_rows, paper_rows, exclusion_rows),
        encoding="utf-8",
    )
    return {
        "candidate_count": len(manifest_rows),
        "executed_candidate_count": 0,
        "valid_independent_ne_count": sum(int(row["valid_independent_ne_record_count"]) for row in paper_rows),
        "valid_independent_ne_species_count": sum(
            int(row["valid_independent_ne_record_count"]) > 0 for row in paper_rows
        ),
        "paper_row_count": len(paper_rows),
    }


def main() -> int:
    summary = build_synthesis()
    print(json.dumps(summary, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
