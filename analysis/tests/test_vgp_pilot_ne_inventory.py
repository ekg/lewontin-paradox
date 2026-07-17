import csv
import json

import pytest

from analysis.build_vgp_pilot_ne_inventory import (
    AVAILABILITY_HEADER,
    SOURCES_HEADER,
    build_inventory,
)
from analysis.validate_vgp_pilot_ne_inventory import validate_inventory


def _write_manifest(path, selected_value="no"):
    header = [
        "candidate_id",
        "scientific_name_source",
        "ncbi_taxid",
        "ncbi_current_name",
        "biosample_accession",
        "individual_or_isolate_id",
        "same_individual_status",
        "same_individual_evidence",
        "h1_accession_version",
        "h2_accession_version",
        "pilot_selected",
    ]
    row = {
        "candidate_id": "candidate-1",
        "scientific_name_source": "Dasypus novemcinctus",
        "ncbi_taxid": "9361",
        "ncbi_current_name": "Dasypus novemcinctus",
        "biosample_accession": "SAMN00000001",
        "individual_or_isolate_id": "specimen-1",
        "same_individual_status": "yes",
        "same_individual_evidence": "shared_biosample",
        "h1_accession_version": "GCA_000001.1",
        "h2_accession_version": "GCA_000002.1",
        "pilot_selected": selected_value,
    }
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t")
        writer.writeheader()
        writer.writerow(row)


def _write_provenance(path, selected_count):
    payload = {
        "candidate_summary": {"selected_count": selected_count},
        "generated_at_utc": "2026-07-17T18:24:11Z",
    }
    path.write_text(json.dumps(payload), encoding="utf-8")


def test_build_inventory_writes_header_only_outputs_for_empty_pilot(tmp_path):
    manifest = tmp_path / "manifest.tsv"
    provenance = tmp_path / "provenance.json"
    sources = tmp_path / "sources.tsv"
    availability = tmp_path / "availability.tsv"
    markdown = tmp_path / "inventory.md"

    _write_manifest(manifest, selected_value="no")
    _write_provenance(provenance, selected_count=0)

    context = build_inventory(manifest, provenance, sources, availability, markdown)
    assert context.selected_count_manifest == 0
    assert sources.read_text(encoding="utf-8").splitlines() == ["\t".join(SOURCES_HEADER)]
    assert availability.read_text(encoding="utf-8").splitlines() == ["\t".join(AVAILABILITY_HEADER)]
    assert "selected_count = 0" in markdown.read_text(encoding="utf-8")


def test_build_inventory_fails_closed_when_species_are_selected(tmp_path):
    manifest = tmp_path / "manifest.tsv"
    provenance = tmp_path / "provenance.json"
    sources = tmp_path / "sources.tsv"
    availability = tmp_path / "availability.tsv"
    markdown = tmp_path / "inventory.md"

    _write_manifest(manifest, selected_value="yes")
    _write_provenance(provenance, selected_count=1)

    with pytest.raises(RuntimeError, match="manual literature curation is required"):
        build_inventory(manifest, provenance, sources, availability, markdown)


def test_validate_inventory_rejects_circular_primary_row(tmp_path):
    manifest = tmp_path / "manifest.tsv"
    provenance = tmp_path / "provenance.json"
    sources = tmp_path / "sources.tsv"
    availability = tmp_path / "availability.tsv"
    markdown = tmp_path / "inventory.md"

    _write_manifest(manifest, selected_value="no")
    _write_provenance(provenance, selected_count=0)
    build_inventory(manifest, provenance, sources, availability, markdown)

    row = {key: "" for key in SOURCES_HEADER}
    row.update(
        {
            "record_id": "bad-row",
            "record_status": "accepted_primary",
            "classification": "independent_primary_covariate",
            "response_dataset_overlap": "derived_from_response",
            "independence_tier": "independent_primary",
        }
    )
    with sources.open("a", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=SOURCES_HEADER, delimiter="\t")
        writer.writerow(row)

    with pytest.raises(RuntimeError, match="independent_primary rows must declare overlap=none"):
        validate_inventory(manifest, provenance, sources, availability, markdown)
