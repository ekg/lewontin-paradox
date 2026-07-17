#!/usr/bin/env python3
"""Validate the pilot Ne inventory against the frozen VGP manifest."""

from __future__ import annotations

import csv
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from analysis.build_vgp_pilot_ne_inventory import (
    AVAILABILITY_HEADER,
    DEFAULT_AVAILABILITY,
    DEFAULT_MANIFEST,
    DEFAULT_MARKDOWN,
    DEFAULT_PROVENANCE,
    DEFAULT_SOURCES,
    SOURCES_HEADER,
    collect_context,
)

ALLOWED_CLASSIFICATIONS = {
    "independent_primary_covariate",
    "literature_ecological_ne",
    "coalescent_scaled_genomic_availability",
    "scenario_scaled_absolute_possibility",
    "excluded_circular",
}

ALLOWED_RECORD_STATUSES = {
    "candidate",
    "accepted_primary",
    "accepted_secondary",
    "accepted_scaled_genomic",
    "accepted_absolute_scenario",
    "excluded",
}

ALLOWED_OVERLAPS = {
    "none",
    "same_species_different_population",
    "same_population_different_samples",
    "shared_samples_different_loci",
    "shared_sites_or_samples",
    "derived_from_response",
    "unknown",
}

ALLOWED_INDEPENDENCE_TIERS = {
    "independent_primary",
    "independent_secondary",
    "partially_shared",
    "genomic_shared",
    "circular_excluded",
}


def _read_tsv(path: Path, expected_header: list[str]) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        observed_header = list(reader.fieldnames or [])
        if observed_header != expected_header:
            raise RuntimeError(
                f"{path} schema mismatch: expected {expected_header!r}, observed {observed_header!r}"
            )
        return list(reader)


def _assert_unique(rows: list[dict[str, str]], key: str, path: Path) -> None:
    seen: set[str] = set()
    for row in rows:
        value = row.get(key, "")
        if not value:
            continue
        if value in seen:
            raise RuntimeError(f"{path} contains duplicate {key} {value!r}")
        seen.add(value)


def validate_inventory(
    manifest_path: Path = DEFAULT_MANIFEST,
    provenance_path: Path = DEFAULT_PROVENANCE,
    sources_path: Path = DEFAULT_SOURCES,
    availability_path: Path = DEFAULT_AVAILABILITY,
    markdown_path: Path = DEFAULT_MARKDOWN,
) -> dict[str, int]:
    context = collect_context(manifest_path, provenance_path)
    if context.selected_count_manifest != context.selected_count_provenance:
        raise RuntimeError("freeze provenance selected_count does not match manifest")
    if context.selected_count_manifest > 6:
        raise RuntimeError("selected_count exceeds the bounded pilot ceiling")

    sources_rows = _read_tsv(sources_path, SOURCES_HEADER)
    availability_rows = _read_tsv(availability_path, AVAILABILITY_HEADER)
    _assert_unique(sources_rows, "record_id", sources_path)
    _assert_unique(availability_rows, "pilot_candidate_id", availability_path)

    selected_ids = {row["candidate_id"] for row in context.selected_rows}
    selected_by_id = {row["candidate_id"]: row for row in context.selected_rows}
    availability_ids = {row["pilot_candidate_id"] for row in availability_rows}
    if availability_ids != selected_ids:
        raise RuntimeError(
            "population-data availability rows do not match the selected pilot denominator"
        )

    for row in sources_rows:
        classification = row["classification"]
        if classification not in ALLOWED_CLASSIFICATIONS:
            raise RuntimeError(f"unsupported classification {classification!r}")

        record_status = row["record_status"]
        if record_status not in ALLOWED_RECORD_STATUSES:
            raise RuntimeError(f"unsupported record_status {record_status!r}")

        overlap = row["response_dataset_overlap"]
        if overlap not in ALLOWED_OVERLAPS:
            raise RuntimeError(f"unsupported response_dataset_overlap {overlap!r}")

        tier = row["independence_tier"]
        if tier not in ALLOWED_INDEPENDENCE_TIERS:
            raise RuntimeError(f"unsupported independence_tier {tier!r}")

        if tier == "independent_primary" and overlap != "none":
            raise RuntimeError(
                "circularity violation: independent_primary rows must declare overlap=none"
            )

        if row["pilot_candidate_id"] and row["pilot_candidate_id"] not in selected_ids:
            raise RuntimeError(
                f"source row references non-selected pilot candidate {row['pilot_candidate_id']!r}"
            )
        if row["pilot_candidate_id"]:
            manifest_row = selected_by_id[row["pilot_candidate_id"]]
            if row["scientific_name"] and row["scientific_name"] != manifest_row["scientific_name_source"]:
                raise RuntimeError(
                    "taxonomy mismatch between sources row and selected manifest "
                    f"for {row['pilot_candidate_id']!r}"
                )
            if row["ncbi_taxid"] and row["ncbi_taxid"] != manifest_row["ncbi_taxid"]:
                raise RuntimeError(
                    "TaxId mismatch between sources row and selected manifest "
                    f"for {row['pilot_candidate_id']!r}"
                )

    for row in availability_rows:
        manifest_row = selected_by_id[row["pilot_candidate_id"]]
        if row["scientific_name"] != manifest_row["scientific_name_source"]:
            raise RuntimeError(
                "taxonomy mismatch between availability row and selected manifest "
                f"for {row['pilot_candidate_id']!r}"
            )
        if row["ncbi_taxid"] != manifest_row["ncbi_taxid"]:
            raise RuntimeError(
                "TaxId mismatch between availability row and selected manifest "
                f"for {row['pilot_candidate_id']!r}"
            )

    markdown_text = markdown_path.read_text(encoding="utf-8")
    required_phrases = [
        "selected_count = 0",
        "species-specific Ne/life-history records retrieved | 0",
        "population-data availability rows emitted | 0",
        "No raw reads, assemblies, VCFs, or other population payloads were downloaded.",
    ]
    for phrase in required_phrases:
        if phrase not in markdown_text:
            raise RuntimeError(f"missing markdown phrase: {phrase!r}")

    if context.selected_count_manifest == 0:
        if sources_rows:
            raise RuntimeError("sources table must be empty when no pilot species are selected")
        if availability_rows:
            raise RuntimeError(
                "population-data availability table must be empty when no pilot species are selected"
            )

    return {
        "selected_count": context.selected_count_manifest,
        "source_rows": len(sources_rows),
        "availability_rows": len(availability_rows),
    }


def main() -> None:
    result = validate_inventory()
    print(
        "VGP_PILOT_NE_INVENTORY_OK"
        f" selected_count={result['selected_count']}"
        f" source_rows={result['source_rows']}"
        f" availability_rows={result['availability_rows']}"
    )


if __name__ == "__main__":
    main()
