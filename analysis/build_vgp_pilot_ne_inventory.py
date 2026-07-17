#!/usr/bin/env python3
"""Build a fail-closed Ne inventory for the frozen VGP pilot manifest."""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = REPO_ROOT / "analysis" / "vgp_pilot_manifest.tsv"
DEFAULT_PROVENANCE = REPO_ROOT / "analysis" / "vgp_phase1_freeze_provenance.json"
DEFAULT_SOURCES = REPO_ROOT / "analysis" / "vgp_pilot_ne_sources.tsv"
DEFAULT_AVAILABILITY = REPO_ROOT / "analysis" / "vgp_pilot_population_data_availability.tsv"
DEFAULT_MARKDOWN = REPO_ROOT / "analysis" / "vgp_pilot_ne_inventory.md"

SELECTED_TRUTHY = {"1", "true", "yes", "y"}

SOURCES_HEADER = [
    "record_id",
    "record_status",
    "classification",
    "exclusion_reason",
    "pilot_candidate_id",
    "scientific_name",
    "ncbi_taxid",
    "taxon_authority",
    "taxon_authority_id",
    "taxon_accepted_name",
    "taxon_scientific_name_source",
    "synonym_resolution",
    "population_id",
    "population_name_source",
    "geography_text",
    "measurement_start_date",
    "measurement_end_date",
    "biological_sampling_unit",
    "estimand_code",
    "quoted_estimand_definition",
    "method_code",
    "sample_size_individuals",
    "sample_size_families",
    "sample_size_loci",
    "estimand_value",
    "estimand_unit",
    "interval_lower",
    "interval_upper",
    "interval_type",
    "source_kind",
    "source_title",
    "source_authors",
    "source_year",
    "source_doi",
    "source_url",
    "source_locator",
    "source_retrieved_utc",
    "mutation_rate_value",
    "mutation_rate_source_record_id",
    "generation_time_value",
    "generation_time_source_record_id",
    "response_dataset_overlap",
    "independence_tier",
    "transfer_distance",
    "uncertainty_notes",
    "exclusion_notes",
]

AVAILABILITY_HEADER = [
    "pilot_candidate_id",
    "scientific_name",
    "ncbi_taxid",
    "biosample_accession",
    "individual_or_isolate_id",
    "same_individual_status",
    "same_individual_evidence",
    "paired_assembly_h1_accession_version",
    "paired_assembly_h2_accession_version",
    "raw_reads_availability_class",
    "raw_reads_locator",
    "callable_diploid_genotypes_availability_class",
    "callable_diploid_genotypes_locator",
    "phased_genomes_availability_class",
    "phased_genomes_locator",
    "population_vcf_availability_class",
    "population_vcf_locator",
    "callability_mask_availability_class",
    "callability_mask_locator",
    "reference_accession_version",
    "demographic_inference_readiness",
    "evidence_summary",
    "retrieval_date_utc",
    "notes",
]


@dataclass(frozen=True)
class InventoryContext:
    selected_rows: tuple[dict[str, str], ...]
    selected_count_manifest: int
    selected_count_provenance: int
    generated_at_utc: str


def _load_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _load_provenance(path: Path) -> dict[str, object]:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def _selected_rows(rows: Iterable[dict[str, str]]) -> list[dict[str, str]]:
    selected = []
    for row in rows:
        value = row.get("pilot_selected", "").strip().lower()
        if value in SELECTED_TRUTHY:
            selected.append(row)
    return selected


def collect_context(manifest_path: Path, provenance_path: Path) -> InventoryContext:
    manifest_rows = _load_manifest(manifest_path)
    selected_rows = tuple(_selected_rows(manifest_rows))
    provenance = _load_provenance(provenance_path)
    provenance_count = int(provenance["candidate_summary"]["selected_count"])
    generated_at_utc = str(provenance["generated_at_utc"])
    return InventoryContext(
        selected_rows=selected_rows,
        selected_count_manifest=len(selected_rows),
        selected_count_provenance=provenance_count,
        generated_at_utc=generated_at_utc,
    )


def _write_tsv(path: Path, header: list[str], rows: Iterable[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _render_markdown(context: InventoryContext) -> str:
    return """# VGP pilot independent Ne inventory

Date: 2026-07-17 UTC

## Scope

This inventory is keyed exactly to `analysis/vgp_pilot_manifest.tsv` and the
upstream freeze provenance in `analysis/vgp_phase1_freeze_provenance.json`.
On 2026-07-17 UTC the frozen manifest contained `pilot_selected = no` for every
row, and the provenance summary reported `selected_count = 0`.

Because no pilot species were selected, this task did not retrieve
species-specific literature or repository metadata. Instead it emits a
fail-closed empty inventory whose schema is ready for downstream joins and
whose validator will stop if a later manifest selects one or more species
without corresponding manual curation.

## Manifest status

| field | value |
| --- | --- |
| selected rows in manifest | {manifest_count} |
| selected rows in provenance summary | {provenance_count} |
| freeze provenance generated at | {generated_at} |
| species-specific Ne/life-history records retrieved | 0 |
| population-data availability rows emitted | 0 |

## Missing data disposition

- No selected pilot species were available for inventory.
- `analysis/vgp_pilot_ne_sources.tsv` is intentionally header-only because the
  denominator of selected taxa is empty.
- `analysis/vgp_pilot_population_data_availability.tsv` is intentionally
  header-only for the same reason.
- If a future manifest changes `pilot_selected` away from all-`no`, rerun this
  task with manual source curation; the builder and validator fail closed in
  that state.

## Circularity and authorization notes

- No Ne estimate was algebraically derived from `pi/(4mu)` or reused as an
  independent predictor.
- No raw reads, assemblies, VCFs, or other population payloads were downloaded.
- No demographic inference, PSMC, MSMC2, or SMC++ execution was attempted.
- No claim is made that VGP H1/H2 assembly pairs are callable diploid genotypes.

## Outputs

- `analysis/vgp_pilot_ne_sources.tsv`
- `analysis/vgp_pilot_population_data_availability.tsv`
- `analysis/vgp_pilot_ne_inventory.md`
""".format(
        manifest_count=context.selected_count_manifest,
        provenance_count=context.selected_count_provenance,
        generated_at=context.generated_at_utc,
    )


def build_inventory(
    manifest_path: Path = DEFAULT_MANIFEST,
    provenance_path: Path = DEFAULT_PROVENANCE,
    sources_path: Path = DEFAULT_SOURCES,
    availability_path: Path = DEFAULT_AVAILABILITY,
    markdown_path: Path = DEFAULT_MARKDOWN,
) -> InventoryContext:
    context = collect_context(manifest_path, provenance_path)
    if context.selected_count_manifest != context.selected_count_provenance:
        raise RuntimeError(
            "selected-count mismatch between manifest and freeze provenance: "
            f"{context.selected_count_manifest} vs {context.selected_count_provenance}"
        )
    if context.selected_count_manifest > 6:
        raise RuntimeError(
            f"manifest selected_count exceeds bounded pilot ceiling: {context.selected_count_manifest}"
        )
    if context.selected_count_manifest != 0:
        raise RuntimeError(
            "selected pilot species are present; manual literature curation is required "
            "before this inventory can be generated"
        )

    _write_tsv(sources_path, SOURCES_HEADER, [])
    _write_tsv(availability_path, AVAILABILITY_HEADER, [])
    markdown_path.write_text(_render_markdown(context), encoding="utf-8")
    return context


def main() -> None:
    build_inventory()


if __name__ == "__main__":
    main()
