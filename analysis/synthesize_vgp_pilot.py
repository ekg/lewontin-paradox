#!/usr/bin/env python3
"""Synthesize the bounded VGP pilot into paper-oriented handoff artifacts."""

from __future__ import annotations

import csv
import json
import statistics
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = ROOT / "analysis" / "vgp_pilot_manifest.tsv"
DEFAULT_REJECTIONS = ROOT / "analysis" / "vgp_pilot_rejections.tsv"
DEFAULT_RESULTS = ROOT / "analysis" / "vgp_pilot_results.tsv"
DEFAULT_TELEMETRY = ROOT / "analysis" / "vgp_pilot_slurm_telemetry.tsv"
DEFAULT_QC = ROOT / "analysis" / "vgp_pilot_qc.tsv"
DEFAULT_RESOURCE = ROOT / "analysis" / "vgp_pilot_resource_calibration.tsv"
DEFAULT_REVIEW = ROOT / "analysis" / "vgp_pilot_review.md"
DEFAULT_NE_SOURCES = ROOT / "analysis" / "vgp_pilot_ne_sources.tsv"
DEFAULT_AVAILABILITY = ROOT / "analysis" / "vgp_pilot_population_data_availability.tsv"
DEFAULT_BUDGET = ROOT / "analysis" / "vertebrate_scaleout_resource_budget.tsv"
DEFAULT_DECISIONS = ROOT / "analysis" / "vertebrate_scaleout_decisions.tsv"

DEFAULT_SYNTHESIS = ROOT / "analysis" / "vgp_pilot_synthesis.md"
DEFAULT_PAPER_TABLE = ROOT / "analysis" / "vgp_pilot_paper_table.tsv"
DEFAULT_NEXT_DECISION = ROOT / "analysis" / "vgp_pilot_next_decision.tsv"

TIER3A_FULL_COUNT = 40
TIER3C_FULL_COUNT = 120
NEXT_WAVE_COUNT = 6

PAPER_TABLE_FIELDS = [
    "versioned_taxon_id",
    "candidate_id",
    "scientific_name",
    "ncbi_taxid",
    "class",
    "order",
    "seed_modalities",
    "h1_accession_version",
    "h2_accession_version",
    "same_individual_status",
    "annotation_file_status",
    "annotation_native_status",
    "callable_bases",
    "callable_fraction",
    "queryable_gene_count",
    "queryable_gene_bases",
    "assembly_composition_eligible",
    "assembly_diversity_eligible",
    "population_genomic_eligible",
    "demographic_eligible",
    "acceptance_status",
    "explicit_rejection_reason",
    "blocking_requirement_ids",
    "pilot_selected",
    "predicted_download_gb_exact",
    "predicted_core_hours_base",
    "predicted_core_hours_high",
    "predicted_peak_memory_gib_base",
    "predicted_peak_memory_gib_high",
    "predicted_wall_hours_base",
    "predicted_wall_hours_high",
    "predicted_scratch_gb_base",
    "predicted_scratch_gb_high",
    "predicted_inode_count_base",
    "predicted_inode_count_high",
    "predicted_moosefs_read_gb_base",
    "predicted_moosefs_read_gb_high",
    "predicted_moosefs_write_gb_base",
    "predicted_moosefs_write_gb_high",
    "predicted_metadata_operations_base",
    "predicted_metadata_operations_high",
    "ne_source_rows",
    "accepted_independent_ne_rows",
    "population_availability_rows",
    "callable_diploid_genotypes_availability_class",
    "population_vcf_availability_class",
    "demographic_inference_readiness",
    "join_status",
]

NEXT_DECISION_FIELDS = [
    "row_type",
    "scenario_id",
    "authorization_state",
    "species_scope",
    "species_count",
    "metric",
    "unit",
    "low",
    "base",
    "high",
    "calibration_status",
    "evidence",
    "notes",
    "recommendation",
]


@dataclass(frozen=True)
class Tier3AObservedStats:
    core_low: float
    core_base: float
    core_high: float
    wall_low: float
    wall_base: float
    wall_high: float
    input_low: float
    input_base: float
    input_high: float
    output_low: float
    output_base: float
    output_high: float
    inode_low: float
    inode_base: float
    inode_high: float
    requested_memory_gib: float


@dataclass(frozen=True)
class ScenarioMetric:
    unit: str
    low: float
    base: float
    high: float
    calibration_status: str
    evidence: str
    notes: str


def load_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fieldnames: Sequence[str], rows: Sequence[Mapping[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def to_float(value: str) -> float:
    return float(value)


def maybe_float(value: str) -> float | None:
    text = value.strip()
    if not text:
        return None
    return float(text)


def format_number(value: float | int | str) -> str:
    if isinstance(value, str):
        return value
    if isinstance(value, int):
        return str(value)
    if value.is_integer():
        return str(int(value))
    return f"{value:.6f}".rstrip("0").rstrip(".")


def bytes_to_gb(value: str) -> float:
    return float(value) / 1_000_000_000.0


def budget_lookup(rows: Sequence[Mapping[str, str]], stage: str, scenario: str) -> dict[str, str]:
    for row in rows:
        if row["stage_or_dataset"] == stage and row["scenario"] == scenario:
            return dict(row)
    raise KeyError(f"missing budget row for stage={stage} scenario={scenario}")


def summarize_tier3a_observed(rows: Sequence[Mapping[str, str]]) -> Tier3AObservedStats:
    observed = [row for row in rows if row["row_type"] == "calibration" and row["cohort_definition"] == "completed_three_species_Tier3A"]
    cores = [to_float(row["core_hours"]) for row in observed]
    walls = [to_float(row["catalog_or_stage_wall_hours"]) for row in observed]
    inputs = [to_float(row["persistent_input_gb"]) for row in observed]
    outputs = [to_float(row["persistent_output_gb"]) for row in observed]
    inodes = [to_float(row["file_inode_count"]) for row in observed]
    requested_memory = {
        float(row["peak_resident_or_requested_memory_gib_per_element"].split("_", 1)[0])
        for row in observed
    }
    if len(requested_memory) != 1:
        raise RuntimeError("Tier3A requested memory should be constant across observed calibrations")
    return Tier3AObservedStats(
        core_low=min(cores),
        core_base=statistics.mean(cores),
        core_high=max(cores),
        wall_low=min(walls),
        wall_base=statistics.mean(walls),
        wall_high=max(walls),
        input_low=min(inputs),
        input_base=statistics.mean(inputs),
        input_high=max(inputs),
        output_low=min(outputs),
        output_base=statistics.mean(outputs),
        output_high=max(outputs),
        inode_low=min(inodes),
        inode_base=statistics.mean(inodes),
        inode_high=max(inodes),
        requested_memory_gib=next(iter(requested_memory)),
    )


def scenario_metric(
    unit: str,
    low: float,
    base: float,
    high: float,
    calibration_status: str,
    evidence: str,
    notes: str,
) -> ScenarioMetric:
    return ScenarioMetric(
        unit=unit,
        low=low,
        base=base,
        high=high,
        calibration_status=calibration_status,
        evidence=evidence,
        notes=notes,
    )


def build_paper_table(
    manifest_rows: Sequence[Mapping[str, str]],
    ne_rows: Sequence[Mapping[str, str]],
    availability_rows: Sequence[Mapping[str, str]],
) -> list[dict[str, str]]:
    ne_by_candidate: dict[str, list[Mapping[str, str]]] = defaultdict(list)
    for row in ne_rows:
        candidate = row.get("pilot_candidate_id", "")
        if candidate:
            ne_by_candidate[candidate].append(row)

    availability_by_candidate: dict[str, list[Mapping[str, str]]] = defaultdict(list)
    for row in availability_rows:
        candidate = row.get("pilot_candidate_id", "")
        if candidate:
            availability_by_candidate[candidate].append(row)

    table_rows: list[dict[str, str]] = []
    for row in manifest_rows:
        candidate = row["candidate_id"]
        matched_ne = ne_by_candidate.get(candidate, [])
        matched_availability = availability_by_candidate.get(candidate, [])
        accepted_independent = [
            record
            for record in matched_ne
            if record.get("record_status") == "accepted_primary"
            and record.get("independence_tier") == "independent_primary"
        ]
        availability = matched_availability[0] if matched_availability else {}
        versioned_taxon_id = f"{candidate}|{row['h1_accession_version']}"
        table_rows.append(
            {
                "versioned_taxon_id": versioned_taxon_id,
                "candidate_id": candidate,
                "scientific_name": row["scientific_name_source"],
                "ncbi_taxid": row["ncbi_taxid"],
                "class": row["class"],
                "order": row["order"],
                "seed_modalities": row["seed_modalities"],
                "h1_accession_version": row["h1_accession_version"],
                "h2_accession_version": row["h2_accession_version"],
                "same_individual_status": row["same_individual_status"],
                "annotation_file_status": row["annotation_file_status"],
                "annotation_native_status": row["annotation_native_status"],
                "callable_bases": row["callable_bases"],
                "callable_fraction": row["callable_fraction"],
                "queryable_gene_count": row["queryable_gene_count"],
                "queryable_gene_bases": row["queryable_gene_bases"],
                "assembly_composition_eligible": row["assembly_composition_eligible"],
                "assembly_diversity_eligible": row["assembly_diversity_eligible"],
                "population_genomic_eligible": row["population_genomic_eligible"],
                "demographic_eligible": row["demographic_eligible"],
                "acceptance_status": row["acceptance_status"],
                "explicit_rejection_reason": row["explicit_acceptance_or_rejection_reason"],
                "blocking_requirement_ids": row["blocking_requirement_ids"],
                "pilot_selected": row["pilot_selected"],
                "predicted_download_gb_exact": format_number(bytes_to_gb(row["predicted_download_bytes_exact"])),
                "predicted_core_hours_base": row["predicted_core_hours_base"],
                "predicted_core_hours_high": row["predicted_core_hours_high"],
                "predicted_peak_memory_gib_base": row["predicted_peak_memory_gib_base"],
                "predicted_peak_memory_gib_high": row["predicted_peak_memory_gib_high"],
                "predicted_wall_hours_base": row["predicted_wall_hours_base"],
                "predicted_wall_hours_high": row["predicted_wall_hours_high"],
                "predicted_scratch_gb_base": row["predicted_scratch_gb_base"],
                "predicted_scratch_gb_high": row["predicted_scratch_gb_high"],
                "predicted_inode_count_base": row["predicted_inode_count_base"],
                "predicted_inode_count_high": row["predicted_inode_count_high"],
                "predicted_moosefs_read_gb_base": row["predicted_moosefs_read_gb_base"],
                "predicted_moosefs_read_gb_high": row["predicted_moosefs_read_gb_high"],
                "predicted_moosefs_write_gb_base": row["predicted_moosefs_write_gb_base"],
                "predicted_moosefs_write_gb_high": row["predicted_moosefs_write_gb_high"],
                "predicted_metadata_operations_base": row["predicted_metadata_operations_base"],
                "predicted_metadata_operations_high": row["predicted_metadata_operations_high"],
                "ne_source_rows": str(len(matched_ne)),
                "accepted_independent_ne_rows": str(len(accepted_independent)),
                "population_availability_rows": str(len(matched_availability)),
                "callable_diploid_genotypes_availability_class": availability.get(
                    "callable_diploid_genotypes_availability_class", ""
                ),
                "population_vcf_availability_class": availability.get("population_vcf_availability_class", ""),
                "demographic_inference_readiness": availability.get("demographic_inference_readiness", ""),
                "join_status": "no_selected_pilot_species_no_inventory_rows"
                if not matched_ne and not matched_availability
                else "matched_inventory_rows_present",
            }
        )
    return table_rows


def build_next_decision_rows(
    budget_rows: Sequence[Mapping[str, str]],
    telemetry_rows: Sequence[Mapping[str, str]],
) -> list[dict[str, str]]:
    telemetry = telemetry_rows[0]
    tier3a = summarize_tier3a_observed(budget_rows)
    stratified_low = budget_lookup(budget_rows, "stratified_pilot", "low")
    stratified_base = budget_lookup(budget_rows, "stratified_pilot", "base")
    stratified_high = budget_lookup(budget_rows, "stratified_pilot", "high")
    tier3c_low = budget_lookup(budget_rows, "tier3c_composition_per_species", "low")
    tier3c_base = budget_lookup(budget_rows, "tier3c_composition_per_species", "base")
    tier3c_high = budget_lookup(budget_rows, "tier3c_composition_per_species", "high")
    full_combined_low = budget_lookup(budget_rows, "full_catalog_combined", "low")
    full_combined_base = budget_lookup(budget_rows, "full_catalog_combined", "base")
    full_combined_high = budget_lookup(budget_rows, "full_catalog_combined", "high")

    current_metrics = {
        "download_gb": scenario_metric(
            "GB",
            0.0,
            0.0,
            0.0,
            "observed_refusal_zero",
            "analysis/vgp_pilot_slurm_telemetry.tsv; analysis/vgp_pilot_run_manifest.tsv",
            "The July 17, 2026 UTC NO_GO refusal submitted no jobs and staged no assets, so executable cost remains zero.",
        ),
        "aggregate_core_hours": scenario_metric(
            "core-h",
            0.0,
            0.0,
            0.0,
            "observed_refusal_zero",
            "analysis/vgp_pilot_slurm_telemetry.tsv",
            "Observed CPU time remained zero because final_state=NOT_SUBMITTED.",
        ),
        "peak_memory_gib_per_job": scenario_metric(
            "GiB",
            0.0,
            0.0,
            0.0,
            "observed_refusal_zero",
            "analysis/vgp_pilot_slurm_telemetry.tsv",
            "No executable job existed under the refusal boundary.",
        ),
        "catalog_wall_hours": scenario_metric(
            "h",
            0.0,
            0.0,
            0.0,
            "observed_refusal_zero",
            "analysis/vgp_pilot_slurm_telemetry.tsv",
            "Observed elapsed time remained zero because no sbatch command was issued.",
        ),
        "peak_local_scratch_gb": scenario_metric(
            "GB",
            0.0,
            0.0,
            0.0,
            "observed_refusal_zero",
            "analysis/vgp_pilot_slurm_telemetry.tsv",
            "Refused preflight created no scratch allocations.",
        ),
        "persistent_input_gb": scenario_metric(
            "GB",
            0.0,
            0.0,
            0.0,
            "observed_refusal_zero",
            "analysis/vgp_pilot_run_manifest.tsv; analysis/vgp_pilot_slurm_telemetry.tsv",
            "No pilot download or staging occurred.",
        ),
        "persistent_output_gb": scenario_metric(
            "GB",
            0.0,
            0.0,
            0.0,
            "observed_refusal_zero",
            "analysis/vgp_pilot_run_manifest.tsv; analysis/vgp_pilot_slurm_telemetry.tsv",
            "No pilot outputs were promoted across the authorization boundary.",
        ),
        "file_inodes": scenario_metric(
            "count",
            0.0,
            0.0,
            0.0,
            "observed_refusal_zero",
            "analysis/vgp_pilot_run_manifest.tsv; analysis/vgp_pilot_slurm_telemetry.tsv",
            "No retained pilot files were created.",
        ),
        "moosefs_read_gb": scenario_metric(
            "GB",
            0.0,
            0.0,
            0.0,
            "observed_refusal_zero",
            "analysis/vgp_pilot_slurm_telemetry.tsv",
            "No executable read traffic occurred under the refusal boundary.",
        ),
        "moosefs_write_gb": scenario_metric(
            "GB",
            0.0,
            0.0,
            0.0,
            "observed_refusal_zero",
            "analysis/vgp_pilot_slurm_telemetry.tsv",
            "No executable write traffic occurred under the refusal boundary.",
        ),
        "metadata_operations": scenario_metric(
            "count",
            0.0,
            0.0,
            0.0,
            "observed_refusal_zero",
            "analysis/vgp_pilot_slurm_telemetry.tsv",
            "No executable metadata operations were recorded.",
        ),
        "peak_bandwidth_mib_s": scenario_metric(
            "MiB/s",
            0.0,
            0.0,
            0.0,
            "observed_refusal_zero",
            "analysis/vgp_pilot_slurm_telemetry.tsv",
            "No executable data-transfer bandwidth was observed.",
        ),
    }

    next_wave_metrics = {
        "download_gb": scenario_metric(
            "GB",
            NEXT_WAVE_COUNT * tier3a.input_low,
            NEXT_WAVE_COUNT * tier3a.input_base,
            NEXT_WAVE_COUNT * tier3a.input_high,
            "observed_historical_proxy",
            "analysis/vertebrate_scaleout_resource_budget.tsv (three Tier3A observed calibrations)",
            "Uses retained staged input as the closest observed download proxy for a repaired <=6-species wave; the refusal pilot itself measured zero transfer bytes.",
        ),
        "aggregate_core_hours": scenario_metric(
            "core-h",
            NEXT_WAVE_COUNT * tier3a.core_low + NEXT_WAVE_COUNT * to_float(tier3c_low["core_hours"]),
            NEXT_WAVE_COUNT * tier3a.core_base + NEXT_WAVE_COUNT * to_float(tier3c_base["core_hours"]),
            NEXT_WAVE_COUNT * tier3a.core_high + NEXT_WAVE_COUNT * to_float(tier3c_high["core_hours"]),
            "observed_historical_proxy",
            "analysis/vertebrate_scaleout_resource_budget.tsv (three Tier3A observed calibrations; Tier3C observed low/mean/max)",
            "Contingent post-repair proxy for a <=6-species wave that would include both paired diversity and exact-native composition on the same six species.",
        ),
        "peak_memory_gib_per_job": scenario_metric(
            "GiB",
            tier3a.requested_memory_gib,
            tier3a.requested_memory_gib,
            tier3a.requested_memory_gib,
            "requested_only_no_actual_rss",
            "analysis/vertebrate_scaleout_resource_budget.tsv (three Tier3A observed calibrations)",
            "Historical Tier3A telemetry recorded 64 GiB requested per element but did not capture actual MaxRSS; do not treat this as measured peak usage.",
        ),
        "catalog_wall_hours": scenario_metric(
            "h",
            NEXT_WAVE_COUNT * tier3a.wall_low + NEXT_WAVE_COUNT * to_float(tier3c_low["wall_hours_per_element"]),
            NEXT_WAVE_COUNT * tier3a.wall_base + NEXT_WAVE_COUNT * to_float(tier3c_base["wall_hours_per_element"]),
            NEXT_WAVE_COUNT * tier3a.wall_high + NEXT_WAVE_COUNT * to_float(tier3c_high["wall_hours_per_element"]),
            "observed_historical_proxy",
            "analysis/vertebrate_scaleout_resource_budget.tsv",
            "Serial catalog wall-hours proxy; parallel execution would still require a new explicit authorization with concurrency caps.",
        ),
        "peak_local_scratch_gb": scenario_metric(
            "GB",
            to_float(stratified_low["local_scratch_peak_gb"]),
            to_float(stratified_base["local_scratch_peak_gb"]),
            to_float(stratified_high["local_scratch_peak_gb"]),
            "planning_envelope_unchanged_no_new_scratch_telemetry",
            "analysis/vertebrate_scaleout_resource_budget.tsv; analysis/vgp_pilot_resource_calibration.tsv",
            "The refusal pilot collected no scratch telemetry, so the prior stratified-pilot scratch envelope remains the only reviewed bound.",
        ),
        "persistent_input_gb": scenario_metric(
            "GB",
            NEXT_WAVE_COUNT * tier3a.input_low,
            NEXT_WAVE_COUNT * tier3a.input_base,
            NEXT_WAVE_COUNT * tier3a.input_high,
            "observed_historical_proxy",
            "analysis/vertebrate_scaleout_resource_budget.tsv (three Tier3A observed calibrations)",
            "Tier3C composition is assumed to reuse the same repaired exact-H1/native-annotation tuples for the same six species.",
        ),
        "persistent_output_gb": scenario_metric(
            "GB",
            NEXT_WAVE_COUNT * (tier3a.output_low + to_float(tier3c_low["persistent_output_gb"])),
            NEXT_WAVE_COUNT * (tier3a.output_base + to_float(tier3c_base["persistent_output_gb"])),
            NEXT_WAVE_COUNT * (tier3a.output_high + to_float(tier3c_high["persistent_output_gb"])),
            "observed_historical_proxy_plus_observed_tier3c_bounds",
            "analysis/vertebrate_scaleout_resource_budget.tsv",
            "Combines observed Tier3A promoted outputs with reviewed Tier3C per-species output bounds.",
        ),
        "file_inodes": scenario_metric(
            "count",
            NEXT_WAVE_COUNT * (tier3a.inode_low + to_float(tier3c_low["file_inode_count"])),
            NEXT_WAVE_COUNT * (tier3a.inode_base + to_float(tier3c_base["file_inode_count"])),
            NEXT_WAVE_COUNT * (tier3a.inode_high + to_float(tier3c_high["file_inode_count"])),
            "observed_historical_proxy_plus_observed_tier3c_bounds",
            "analysis/vertebrate_scaleout_resource_budget.tsv",
            "Inode totals combine observed Tier3A retained files with Tier3C per-species bounds.",
        ),
        "moosefs_read_gb": scenario_metric(
            "GB",
            to_float(stratified_low["moosefs_read_gb"]),
            to_float(stratified_base["moosefs_read_gb"]),
            to_float(stratified_high["moosefs_read_gb"]),
            "planning_envelope_unchanged_no_new_io_telemetry",
            "analysis/vertebrate_scaleout_resource_budget.tsv; analysis/vgp_pilot_resource_calibration.tsv",
            "The refusal pilot collected no new read-byte telemetry; prior stratified-pilot bounds remain unchanged.",
        ),
        "moosefs_write_gb": scenario_metric(
            "GB",
            to_float(stratified_low["moosefs_write_gb"]),
            to_float(stratified_base["moosefs_write_gb"]),
            to_float(stratified_high["moosefs_write_gb"]),
            "planning_envelope_unchanged_no_new_io_telemetry",
            "analysis/vertebrate_scaleout_resource_budget.tsv; analysis/vgp_pilot_resource_calibration.tsv",
            "The refusal pilot collected no new write-byte telemetry; prior stratified-pilot bounds remain unchanged.",
        ),
        "metadata_operations": scenario_metric(
            "count",
            to_float(stratified_low["metadata_operations"]),
            to_float(stratified_base["metadata_operations"]),
            to_float(stratified_high["metadata_operations"]),
            "planning_envelope_unchanged_no_new_metadata_telemetry",
            "analysis/vertebrate_scaleout_resource_budget.tsv; analysis/vgp_pilot_resource_calibration.tsv",
            "No executed pilot jobs existed from which to refit metadata-operation envelopes.",
        ),
        "peak_bandwidth_mib_s": scenario_metric(
            "MiB/s",
            to_float(stratified_low["peak_aggregate_bandwidth_mib_s"]),
            to_float(stratified_base["peak_aggregate_bandwidth_mib_s"]),
            to_float(stratified_high["peak_aggregate_bandwidth_mib_s"]),
            "planning_envelope_unchanged_no_new_bandwidth_telemetry",
            "analysis/vertebrate_scaleout_resource_budget.tsv; analysis/vgp_pilot_resource_calibration.tsv",
            "Bandwidth remained unresolved after the refusal pilot because no jobs executed.",
        ),
    }

    full_catalog_metrics = {
        "download_gb": scenario_metric(
            "GB",
            to_float(full_combined_low["persistent_input_gb"]),
            to_float(full_combined_base["persistent_input_gb"]),
            to_float(full_combined_high["persistent_input_gb"]),
            "planning_envelope_unchanged_due_overlap_and_missing_transfer_telemetry",
            "analysis/vertebrate_scaleout_resource_budget.tsv; analysis/vertebrate_scaleout_decisions.tsv",
            "Checksum-deduplicated full-catalog transfer/storage remains the reviewed planning envelope because the refusal pilot produced no overlap-aware transfer telemetry.",
        ),
        "aggregate_core_hours": scenario_metric(
            "core-h",
            TIER3A_FULL_COUNT * tier3a.core_low + TIER3C_FULL_COUNT * to_float(tier3c_low["core_hours"]),
            TIER3A_FULL_COUNT * tier3a.core_base + TIER3C_FULL_COUNT * to_float(tier3c_base["core_hours"]),
            TIER3A_FULL_COUNT * tier3a.core_high + TIER3C_FULL_COUNT * to_float(tier3c_high["core_hours"]),
            "observed_historical_proxy",
            "analysis/vertebrate_scaleout_resource_budget.tsv; analysis/vertebrate_scaleout_decisions.tsv",
            "Uses 40 paired Tier3A proxies plus 120 Tier3C composition proxies; this is contingent and non-authorizing.",
        ),
        "peak_memory_gib_per_job": scenario_metric(
            "GiB",
            tier3a.requested_memory_gib,
            tier3a.requested_memory_gib,
            tier3a.requested_memory_gib,
            "requested_only_no_actual_rss",
            "analysis/vertebrate_scaleout_resource_budget.tsv",
            "Historical Tier3A requested memory remained 64 GiB per element; actual RSS was not captured.",
        ),
        "catalog_wall_hours": scenario_metric(
            "h",
            TIER3A_FULL_COUNT * tier3a.wall_low + TIER3C_FULL_COUNT * to_float(tier3c_low["wall_hours_per_element"]),
            TIER3A_FULL_COUNT * tier3a.wall_base + TIER3C_FULL_COUNT * to_float(tier3c_base["wall_hours_per_element"]),
            TIER3A_FULL_COUNT * tier3a.wall_high + TIER3C_FULL_COUNT * to_float(tier3c_high["wall_hours_per_element"]),
            "observed_historical_proxy",
            "analysis/vertebrate_scaleout_resource_budget.tsv",
            "Serial wall-hours proxy only; any actual concurrency plan requires a fresh authorization packet.",
        ),
        "peak_local_scratch_gb": scenario_metric(
            "GB",
            to_float(full_combined_low["local_scratch_peak_gb"]),
            to_float(full_combined_base["local_scratch_peak_gb"]),
            to_float(full_combined_high["local_scratch_peak_gb"]),
            "planning_envelope_unchanged_no_new_scratch_telemetry",
            "analysis/vertebrate_scaleout_resource_budget.tsv; analysis/vgp_pilot_resource_calibration.tsv",
            "The refusal pilot provided no scratch telemetry to tighten the full-catalog bound.",
        ),
        "persistent_input_gb": scenario_metric(
            "GB",
            to_float(full_combined_low["persistent_input_gb"]),
            to_float(full_combined_base["persistent_input_gb"]),
            to_float(full_combined_high["persistent_input_gb"]),
            "planning_envelope_unchanged_due_overlap_and_missing_transfer_telemetry",
            "analysis/vertebrate_scaleout_resource_budget.tsv; analysis/vertebrate_scaleout_decisions.tsv",
            "Uses the reviewed checksum-deduplicated full-catalog storage envelope rather than inventing overlap from zero-transfer refusal telemetry.",
        ),
        "persistent_output_gb": scenario_metric(
            "GB",
            TIER3A_FULL_COUNT * tier3a.output_low + TIER3C_FULL_COUNT * to_float(tier3c_low["persistent_output_gb"]),
            TIER3A_FULL_COUNT * tier3a.output_base + TIER3C_FULL_COUNT * to_float(tier3c_base["persistent_output_gb"]),
            TIER3A_FULL_COUNT * tier3a.output_high + TIER3C_FULL_COUNT * to_float(tier3c_high["persistent_output_gb"]),
            "observed_historical_proxy_plus_observed_tier3c_bounds",
            "analysis/vertebrate_scaleout_resource_budget.tsv",
            "Combines observed Tier3A promoted outputs with Tier3C reviewed per-species output bounds.",
        ),
        "file_inodes": scenario_metric(
            "count",
            TIER3A_FULL_COUNT * tier3a.inode_low + TIER3C_FULL_COUNT * to_float(tier3c_low["file_inode_count"]),
            TIER3A_FULL_COUNT * tier3a.inode_base + TIER3C_FULL_COUNT * to_float(tier3c_base["file_inode_count"]),
            TIER3A_FULL_COUNT * tier3a.inode_high + TIER3C_FULL_COUNT * to_float(tier3c_high["file_inode_count"]),
            "observed_historical_proxy_plus_observed_tier3c_bounds",
            "analysis/vertebrate_scaleout_resource_budget.tsv",
            "Reported full catalog assumes 40 paired Tier3A species overlapping a 120-species Tier3C composition catalog.",
        ),
        "moosefs_read_gb": scenario_metric(
            "GB",
            to_float(full_combined_low["moosefs_read_gb"]),
            to_float(full_combined_base["moosefs_read_gb"]),
            to_float(full_combined_high["moosefs_read_gb"]),
            "planning_envelope_unchanged_no_new_io_telemetry",
            "analysis/vertebrate_scaleout_resource_budget.tsv; analysis/vgp_pilot_resource_calibration.tsv",
            "Full-catalog read traffic remains the prior reviewed envelope because the refusal pilot executed no jobs.",
        ),
        "moosefs_write_gb": scenario_metric(
            "GB",
            to_float(full_combined_low["moosefs_write_gb"]),
            to_float(full_combined_base["moosefs_write_gb"]),
            to_float(full_combined_high["moosefs_write_gb"]),
            "planning_envelope_unchanged_no_new_io_telemetry",
            "analysis/vertebrate_scaleout_resource_budget.tsv; analysis/vgp_pilot_resource_calibration.tsv",
            "Full-catalog write traffic remains unchanged after the refusal pilot.",
        ),
        "metadata_operations": scenario_metric(
            "count",
            to_float(full_combined_low["metadata_operations"]),
            to_float(full_combined_base["metadata_operations"]),
            to_float(full_combined_high["metadata_operations"]),
            "planning_envelope_unchanged_no_new_metadata_telemetry",
            "analysis/vertebrate_scaleout_resource_budget.tsv; analysis/vgp_pilot_resource_calibration.tsv",
            "No executed pilot jobs existed from which to refit full-catalog metadata load.",
        ),
        "peak_bandwidth_mib_s": scenario_metric(
            "MiB/s",
            to_float(full_combined_low["peak_aggregate_bandwidth_mib_s"]),
            to_float(full_combined_base["peak_aggregate_bandwidth_mib_s"]),
            to_float(full_combined_high["peak_aggregate_bandwidth_mib_s"]),
            "planning_envelope_unchanged_no_new_bandwidth_telemetry",
            "analysis/vertebrate_scaleout_resource_budget.tsv; analysis/vgp_pilot_resource_calibration.tsv",
            "Full-catalog bandwidth remained unresolved because the refusal pilot executed no work.",
        ),
    }

    rows: list[dict[str, str]] = [
        {
            "row_type": "recommendation",
            "scenario_id": "pilot_next_decision",
            "authorization_state": "current_review_only",
            "species_scope": "bounded_vgp_pilot",
            "species_count": "0",
            "metric": "",
            "unit": "",
            "low": "",
            "base": "",
            "high": "",
            "calibration_status": "review_conclusion",
            "evidence": "analysis/vgp_pilot_review.md; analysis/vgp_pilot_qc.tsv; analysis/vgp_pilot_ne_inventory.md; analysis/vertebrate_scaleout_decisions.tsv",
            "notes": (
                "Recommendation is stop_repair: zero species crossed the authorization boundary, no independent Ne/ecological inventory rows were populated, and any expansion/full-catalog/demographic action requires a new explicit user authorization packet analogous to A50/A60/A70/A71."
            ),
            "recommendation": "stop_repair",
        }
    ]

    def add_scenario(scenario_id: str, authorization_state: str, species_scope: str, species_count: int, metrics: Mapping[str, ScenarioMetric]) -> None:
        for metric_name, metric in metrics.items():
            rows.append(
                {
                    "row_type": "scenario_metric",
                    "scenario_id": scenario_id,
                    "authorization_state": authorization_state,
                    "species_scope": species_scope,
                    "species_count": str(species_count),
                    "metric": metric_name,
                    "unit": metric.unit,
                    "low": format_number(metric.low),
                    "base": format_number(metric.base),
                    "high": format_number(metric.high),
                    "calibration_status": metric.calibration_status,
                    "evidence": metric.evidence,
                    "notes": metric.notes,
                    "recommendation": "",
                }
            )

    add_scenario(
        scenario_id="current_no_go_executable_cost",
        authorization_state=telemetry["final_state"],
        species_scope="current_frozen_manifest",
        species_count=0,
        metrics=current_metrics,
    )
    add_scenario(
        scenario_id="contingent_repaired_next_wave_leq6",
        authorization_state="requires_new_explicit_authorization",
        species_scope="another_bounded_wave_after_repair",
        species_count=NEXT_WAVE_COUNT,
        metrics=next_wave_metrics,
    )
    add_scenario(
        scenario_id="contingent_repaired_reported_full_catalog",
        authorization_state="requires_new_explicit_authorization",
        species_scope="reported_full_catalog_after_repair",
        species_count=TIER3C_FULL_COUNT,
        metrics=full_catalog_metrics,
    )
    return rows


def summarize_classes(manifest_rows: Sequence[Mapping[str, str]]) -> list[str]:
    by_class: dict[str, dict[str, int]] = defaultdict(lambda: {"rows": 0, "tier3a_capable": 0, "tier3c_only": 0})
    for row in manifest_rows:
        class_name = row["class"]
        by_class[class_name]["rows"] += 1
        if "tier3a_seed" in row["seed_modalities"]:
            by_class[class_name]["tier3a_capable"] += 1
        else:
            by_class[class_name]["tier3c_only"] += 1
    rendered = []
    for class_name in sorted(by_class):
        counts = by_class[class_name]
        rendered.append(
            f"- `{class_name}`: {counts['rows']} candidate rows, {counts['tier3a_capable']} Tier3A-capable, {counts['tier3c_only']} Tier3C-only."
        )
    return rendered


def render_markdown(
    manifest_rows: Sequence[Mapping[str, str]],
    rejection_rows: Sequence[Mapping[str, str]],
    results_rows: Sequence[Mapping[str, str]],
    qc_rows: Sequence[Mapping[str, str]],
    resource_rows: Sequence[Mapping[str, str]],
    ne_rows: Sequence[Mapping[str, str]],
    availability_rows: Sequence[Mapping[str, str]],
    next_decision_rows: Sequence[Mapping[str, str]],
) -> str:
    row_count = len(manifest_rows)
    selected_count = sum(1 for row in manifest_rows if row["pilot_selected"] == "yes")
    class_counts = Counter(row["class"] for row in manifest_rows)
    reason_counts = Counter(row["explicit_acceptance_or_rejection_reason"] for row in manifest_rows)
    blocking_counts = Counter(row["blocking_requirement_ids"] for row in manifest_rows)
    same_individual_counts = Counter(row["same_individual_status"] for row in manifest_rows)
    qc_pass = sum(1 for row in qc_rows if row["decision"] == "PASS")
    qc_fail = sum(1 for row in qc_rows if row["decision"] == "FAIL")
    resource_pass = sum(1 for row in resource_rows if row["decision"] == "PASS")
    resource_inconclusive = sum(1 for row in resource_rows if row["decision"] == "INCONCLUSIVE")
    validated_species_row = next(row for row in results_rows if row["metric"] == "validated_species_count")
    current_cost_rows = [row for row in next_decision_rows if row["scenario_id"] == "current_no_go_executable_cost" and row["metric"]]
    current_cost = {row["metric"]: row["base"] for row in current_cost_rows}
    next_wave_core = next(
        row for row in next_decision_rows if row["scenario_id"] == "contingent_repaired_next_wave_leq6" and row["metric"] == "aggregate_core_hours"
    )
    full_core = next(
        row
        for row in next_decision_rows
        if row["scenario_id"] == "contingent_repaired_reported_full_catalog" and row["metric"] == "aggregate_core_hours"
    )

    lines = [
        "# VGP pilot synthesis handoff",
        "",
        "Date: 2026-07-17 UTC",
        "",
        "## Executive readout",
        "",
        "- Recommended next decision: `stop_repair`.",
        f"- Frozen pilot manifest rows reviewed: `{row_count}`; pilot-selected rows: `{selected_count}`; validated executable species: `{validated_species_row['value']}`.",
        "- Actual bounded pilot execution never crossed the authorization boundary: `final_state = NOT_SUBMITTED`, `sbatch_command = none`, and every executable cost metric remained zero.",
        f"- Review QC remained fail-closed (`PASS={qc_pass}`, `FAIL={qc_fail}`) and resource calibration stayed refusal-only (`PASS={resource_pass}`, `INCONCLUSIVE={resource_inconclusive}`).",
        "- Independent Ne/ecological inventory remained intentionally empty because the selected-species denominator was zero; no literature Ne, raw-read, VCF, or callable-genotype claims were populated.",
        "",
        "## Scope and provenance",
        "",
        "- This synthesis joins `analysis/vgp_pilot_manifest.tsv`, `analysis/vgp_pilot_rejections.tsv`, `analysis/vgp_pilot_results.tsv`, `analysis/vgp_pilot_qc.tsv`, `analysis/vgp_pilot_resource_calibration.tsv`, `analysis/vgp_pilot_ne_sources.tsv`, and `analysis/vgp_pilot_population_data_availability.tsv`.",
        "- The join key is versioned taxon identity `candidate_id|h1_accession_version`; every one of the 74 manifest rows appears in `analysis/vgp_pilot_paper_table.tsv` and every joined Ne/availability count is zero because both inventory files are header-only.",
        "- No additional download, Slurm submission, or demographic inference was launched while producing this handoff.",
        "",
        "## Which strata were successfully analyzed",
        "",
        "- No vertebrate class or stratum crossed the executable pilot gate. The `<=6`-species bounded pilot analyzed zero species biologically and emitted zero coding/CDS/fourfold estimates.",
        "- What did survive review was only precondition evidence: the accepted SweepGA build remained byte-identical, and the IMPG smoke artifact still showed native 1:1 mapping depth with `callable_bp = 14507`, `queryable_gene_count = 3`, and `queryable_gene_bp = 14507` on the exact-native sentinel assembly.",
        "- Those sentinel mapping/queryability/callability checks are not promoted to any pilot-species diversity, composition, population-genomic, or demographic result.",
        "",
        "### Frozen candidate strata reviewed",
        "",
        *summarize_classes(manifest_rows),
        "",
        "## Exact eligibility and failure reasons",
        "",
        "- Eligibility stayed zero for all four modalities: every row has `assembly_composition_eligible = no`, `assembly_diversity_eligible = no`, `population_genomic_eligible = no`, and `demographic_eligible = no`.",
        f"- Rejection ledger size equals manifest size: `{len(rejection_rows)}` rejected rows, `{class_counts['Mammalia']}` mammals, `{class_counts['Aves']}` birds, `{class_counts['Lepidosauria']}` lepidosaurs, `{class_counts['Amphibia']}` amphibians, `{class_counts['Actinopteri']}` actinopterygians, `{class_counts['Chondrichthyes']}` chondrichthyans, and `{class_counts['Cladistia']}` cladistian row, plus `{class_counts['UNRESOLVED']}` rows whose class label was unresolved in the frozen source.",
        f"- Blocking requirement counts were exact: `B03` alone for `{blocking_counts['B03']}` Tier3C-only rows and `B02;B03` for `{blocking_counts['B02;B03']}` Tier3A-capable rows.",
        f"- Annotation evidence failed for every row: `annotation_file_status = missing` for all 74 rows and `annotation_native_status = unresolved_or_missing` for every rejection row.",
        f"- Same-individual pairing was present for `{same_individual_counts['yes']}` rows but absent for `{same_individual_counts['no']}` rows; even the paired rows remained blocked because native exact-H1 annotation and denominator evidence were missing.",
        "",
        "### Rejection reason ledger",
        "",
    ]
    for reason, count in reason_counts.most_common():
        lines.append(f"- `{reason}`: `{count}` rows.")

    lines.extend(
        [
            "",
            "### Gate blockers copied into the refusal results",
            "",
            "- `SOURCE_COUNT_DISCREPANCY_UNRESOLVED`: the frozen raw VGP catalog still disagreed with earlier planning headline counts.",
            "- `NO_SELECTED_ROWS`: the frozen pilot manifest selected zero rows, so no bounded pilot could be authorized.",
            "- `ZERO_COMPOSITION_ELIGIBLE_ROWS`: no row independently satisfied exact-H1/native-annotation/denominator requirements.",
            "- `ZERO_DIVERSITY_ELIGIBLE_ROWS`: no row independently satisfied paired same-individual diversity requirements.",
            "- `QUOTA_UNAVAILABLE`: the environment exposed free space but no user-visible quota command, so the storage gate failed closed.",
            "",
            "## Diversity, composition, and uncertainty",
            "",
            "- Supported paper claim: the bounded pilot produced **no** promoted cross-species diversity or composition estimate. `validated_species_count = 0` is the exact result, not a missing-value placeholder.",
            "- Supported paper claim: the only numeric mapping/queryability/callability evidence is the sentinel IMPG smoke artifact described above; it remains separate from biological pilot outputs.",
            "- Unsupported paper claim: any coding diversity ratio, CDS/fourfold composition estimate, or class-level summary derived from this bounded pilot. None was authorized or computed.",
            "- Unsupported paper claim: any population-genomic or demographic inference from VGP H1/H2 pairs. No callable diploid genotype set, no raw-read cohort, no VCF, and no phased demographic input was reviewed into existence here.",
            "",
            "## Mapping, queryability, and callability behavior",
            "",
            "- Global precondition behavior remained stable: exact-native sentinel mapping was 1:1 and queryable/callable denominators stayed positive in the handoff smoke artifact.",
            "- Candidate-row behavior remained unresolved: every row in `analysis/vgp_pilot_paper_table.tsv` carries unresolved or missing `callable_bases`, `queryable_gene_count`, and `queryable_gene_bases` fields, so no row can support composition or diversity claims.",
            "- This separation matters for the paper: sentinel assembly-engineering evidence supports toolchain integrity only; it does not support taxon-level biological inference.",
            "",
            "## Independent Ne and ecological metadata inventory",
            "",
            f"- Usable independent Ne/life-history rows now: `{len(ne_rows)}`. Population-data availability rows now: `{len(availability_rows)}`.",
            "- Supported paper claim: no independent Ne estimate, life-history covariate, or ecological covariate entered the reviewed pilot because the selected-species denominator was zero and the inventory was intentionally header-only.",
            "- Supported paper claim: no circular predictor slipped in. The independent inventory explicitly rejected `pi/(4mu)`-style algebraic back-calculations, shared-sample genomic histories, and any claim that VGP H1/H2 pairs constitute callable diploid genotypes.",
            "- Unsupported paper claim: that later PSMC, MSMC2, SMC++, or population-VCF work is already feasible for any reviewed pilot species. Every such method remains unsupported here because there are zero curated population-data rows.",
            "",
            "## Supported and unsupported paper claims",
            "",
            "### Supported",
            "",
            "- The bounded VGP pilot remained a fail-closed refusal on July 17, 2026 UTC.",
            "- Zero species were selected, zero species were executed, and zero diversity/composition outputs were promoted.",
            "- All 74 reviewed vertebrate candidates failed exact eligibility for documented annotation/denominator and, where relevant, pairing reasons.",
            "- Current executable cost under the frozen gate is zero across download, compute, memory, scratch, storage, inode, and I/O dimensions.",
            "",
            "### Unsupported",
            "",
            "- Any biological estimate for the bounded pilot beyond `validated_species_count = 0`.",
            "- Any demographic or population-genomic claim from assembly pairs alone.",
            "- Any use of algebraically derived Ne or overlapping genomic histories as independent predictors.",
            "- Any interpretation that this synthesis authorizes another wave, a full catalog, bulk acquisition, or PSMC/MSMC2/SMC++ work.",
            "",
            "## Resource calibration and cost consequences",
            "",
            f"- Current executable cost is exactly zero: download `{current_cost['download_gb']}` GB, compute `{current_cost['aggregate_core_hours']}` core-h, wall `{current_cost['catalog_wall_hours']}` h, persistent input `{current_cost['persistent_input_gb']}` GB, persistent output `{current_cost['persistent_output_gb']}` GB, inodes `{current_cost['file_inodes']}`, read `{current_cost['moosefs_read_gb']}` GB, write `{current_cost['moosefs_write_gb']}` GB, metadata ops `{current_cost['metadata_operations']}`.",
            f"- Contingent post-repair <=6-species wave proxy: `{next_wave_core['low']}` / `{next_wave_core['base']}` / `{next_wave_core['high']}` core-h. Metrics that lacked executed-job telemetry remain explicitly labeled as unchanged planning envelopes in `analysis/vgp_pilot_next_decision.tsv`.",
            f"- Contingent post-repair reported full catalog proxy: `{full_core['low']}` / `{full_core['base']}` / `{full_core['high']}` core-h. Storage and I/O for the overlap-deduplicated full catalog remain planning envelopes because the refusal pilot produced zero transfer telemetry.",
            "- The key change from this bounded pilot is therefore negative but important: executable cost is now known to be zero under the present gate, while any non-zero future budget still depends on a repair step and a new explicit authorization packet.",
            "",
            "## Terminal-state and retention verification",
            "",
            "- Job/download terminal state from committed artifacts: `analysis/vgp_pilot_slurm_telemetry.tsv` records a single `run_summary` row with `final_state = NOT_SUBMITTED`, blank job identifiers, and `failure_code = GATE_NO_GO`.",
            "- No pilot retrieval or staging manifest was populated, and no biological output tree crossed the boundary; the refusal artifacts themselves are the only retained pilot outputs.",
            "- Retention and quarantine policy remains the reviewed one in `analysis/vertebrate_scaleout_decisions.tsv`: provenance, final refusal artifacts, failure ledgers, and telemetry are retained; unresolved/partial execution objects do not exist here and therefore there is nothing to quarantine or resume.",
            "",
            "## Recommendation",
            "",
            "- Choose `stop/repair`, not `another bounded expansion wave`, `full eligible-catalog consideration`, or `deferred`.",
            "- Repair means resolving the source-count discrepancy, exposing a user-visible quota interface, and producing row-level native-annotation and denominator evidence before any new selection occurs.",
            "- If a human later wants more work, it must be represented as a new explicit authorization task with numeric species scope and budgets analogous to `A50`/`A60`/`A70`/`A71`. This synthesis does not create or imply a ready executable node across that boundary.",
            "",
        ]
    )
    return "\n".join(lines)


def build_synthesis(
    manifest_path: Path = DEFAULT_MANIFEST,
    rejections_path: Path = DEFAULT_REJECTIONS,
    results_path: Path = DEFAULT_RESULTS,
    telemetry_path: Path = DEFAULT_TELEMETRY,
    qc_path: Path = DEFAULT_QC,
    resource_path: Path = DEFAULT_RESOURCE,
    ne_sources_path: Path = DEFAULT_NE_SOURCES,
    availability_path: Path = DEFAULT_AVAILABILITY,
    budget_path: Path = DEFAULT_BUDGET,
    synthesis_out: Path = DEFAULT_SYNTHESIS,
    paper_table_out: Path = DEFAULT_PAPER_TABLE,
    next_decision_out: Path = DEFAULT_NEXT_DECISION,
) -> dict[str, Any]:
    manifest_rows = load_tsv(manifest_path)
    rejection_rows = load_tsv(rejections_path)
    results_rows = load_tsv(results_path)
    telemetry_rows = load_tsv(telemetry_path)
    qc_rows = load_tsv(qc_path)
    resource_rows = load_tsv(resource_path)
    ne_rows = load_tsv(ne_sources_path)
    availability_rows = load_tsv(availability_path)
    budget_rows = load_tsv(budget_path)

    paper_rows = build_paper_table(manifest_rows, ne_rows, availability_rows)
    next_decision_rows = build_next_decision_rows(budget_rows, telemetry_rows)
    synthesis_text = render_markdown(
        manifest_rows=manifest_rows,
        rejection_rows=rejection_rows,
        results_rows=results_rows,
        qc_rows=qc_rows,
        resource_rows=resource_rows,
        ne_rows=ne_rows,
        availability_rows=availability_rows,
        next_decision_rows=next_decision_rows,
    )

    synthesis_out.write_text(synthesis_text, encoding="utf-8")
    write_tsv(paper_table_out, PAPER_TABLE_FIELDS, paper_rows)
    write_tsv(next_decision_out, NEXT_DECISION_FIELDS, next_decision_rows)

    return {
        "manifest_rows": len(manifest_rows),
        "paper_rows": len(paper_rows),
        "next_decision_rows": len(next_decision_rows),
        "selected_count": sum(1 for row in manifest_rows if row["pilot_selected"] == "yes"),
        "matched_ne_rows": len(ne_rows),
        "matched_availability_rows": len(availability_rows),
    }


def main() -> None:
    build_synthesis()


if __name__ == "__main__":
    main()
