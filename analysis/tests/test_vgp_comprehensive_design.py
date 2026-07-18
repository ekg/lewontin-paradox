import copy
import csv
import importlib.util
import json
from pathlib import Path

import pytest


ROOT = Path(__file__).parents[2]
MODULE_PATH = ROOT / "analysis/assert_vgp_comprehensive_design.py"
SPEC = importlib.util.spec_from_file_location("assert_vgp_comprehensive_design", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
design = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(design)


def test_frozen_design_artifacts_pass_all_assertions():
    assert design.validate_design() == []


def test_exact_primary_and_alternate_counts():
    _, primary = design.read_tsv(design.PRIMARY)
    _, alternates = design.read_tsv(design.ALTERNATES)

    assert len(primary) == 10
    assert len(alternates) == 6
    assert [row["selection_id"] for row in primary] == [f"P{i:02d}" for i in range(1, 11)]
    assert [row["selection_id"] for row in alternates] == [f"A{i:02d}" for i in range(1, 7)]


def test_primary_panel_spans_required_strata():
    _, rows = design.read_tsv(design.PRIMARY)

    assert {row["clade"] for row in rows} == {
        "Mammalia", "Aves", "Reptilia", "Amphibia", "Actinopterygii", "Chondrichthyes"
    }
    assert {row["genome_size_stratum"] for row in rows} == {"small", "medium", "large"}
    assert {row["expected_diversity_stratum"] for row in rows} == {"low", "medium", "high"}
    assert any(row["assembly_generation"].startswith("early") for row in rows)
    assert any(row["assembly_generation"].startswith("later") for row in rows)


def test_no_unresolved_pair_identity_or_reused_accession():
    rows = []
    for path in (design.PRIMARY, design.ALTERNATES):
        _, current = design.read_tsv(path)
        rows.extend(current)

    accessions = [row[f"h{hap}_accession_version"] for row in rows for hap in (1, 2)]
    assert len(accessions) == len(set(accessions))
    for row in rows:
        assert row["catalog_taxid"] == row["resolved_taxid"]
        assert row["biosample"].startswith(("SAMN", "SAMEA"))
        assert "shared BioSample and isolate" in row["reciprocal_pair_evidence"]
        assert "each report links the other accession" in row["reciprocal_pair_evidence"]
        assert "UNRESOLVED" not in "\t".join(row.values())


def test_schema_rejects_annotation_as_core_gate():
    import jsonschema

    instance = json.loads(design.MANIFEST.read_text(encoding="utf-8"))
    schema = json.loads(design.SCHEMA.read_text(encoding="utf-8"))
    instance["annotation_policy"]["required_for_core"] = True

    with pytest.raises(jsonschema.ValidationError):
        jsonschema.Draft202012Validator(schema).validate(instance)


def test_schema_rejects_download_authorization():
    import jsonschema

    instance = json.loads(design.MANIFEST.read_text(encoding="utf-8"))
    schema = json.loads(design.SCHEMA.read_text(encoding="utf-8"))
    instance["authorization"]["biological_downloads_authorized"] = True

    with pytest.raises(jsonschema.ValidationError):
        jsonschema.Draft202012Validator(schema).validate(instance)


def test_validator_rejects_identity_mutation(tmp_path, monkeypatch):
    header, rows = design.read_tsv(design.PRIMARY)
    rows = copy.deepcopy(rows)
    rows[0]["h2_accession_version"] = "GCA_000000000.1"
    mutated = tmp_path / "primary.tsv"
    with mutated.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    errors = []
    _, loaded = design.read_tsv(mutated)
    design._validate_rows(mutated, loaded, design.EXPECTED_PRIMARY, errors)

    assert any("identity" in error for error in errors)


def test_quality_reviewed_graph_assertions_remain_present():
    graph_assertion = ROOT / "analysis/assert_vgp_comprehensive_wg.py"
    graph_tests = ROOT / "analysis/tests/test_assert_vgp_comprehensive_wg.py"

    assert graph_assertion.is_file()
    assert graph_tests.is_file()
    graph_spec = importlib.util.spec_from_file_location(
        "vgp_comprehensive_graph_assertions", graph_assertion
    )
    assert graph_spec is not None and graph_spec.loader is not None
    graph = importlib.util.module_from_spec(graph_spec)
    graph_spec.loader.exec_module(graph)
    assert sum(len(parents) for parents in graph.EXPECTED_AFTER.values()) == 17
    text = graph_assertion.read_text(encoding="utf-8")
    assert "population and non-allelic branches remain design-only" in text
