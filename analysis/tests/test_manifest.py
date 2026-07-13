import copy
import csv
import json
from pathlib import Path

import pytest

from analysis.tier3_common import Tier3ValidationError
from analysis.tier3_manifest import (
    BUFFALO_CORE_COUNT,
    BUFFALO_SHA256,
    buffalo_core_rows,
    load_and_validate_manifest,
    validate_manifest,
)


FIXTURES = Path(__file__).parent / "fixtures"


def manifest_fixture():
    return json.loads((FIXTURES / "manifest.valid.json").read_text())


def test_manifest_schema_and_semantic_validation_accept_truth_fixture():
    validated = validate_manifest(manifest_fixture(), verify_local_files=False)
    assert validated["datasets"][0]["diversity_eligibility"] == "eligible"


def test_manifest_rejects_assembly_coordinate_mismatch():
    manifest = manifest_fixture()
    manifest["datasets"][0]["annotation"]["assembly_accession"] = "GCA_999999.2"
    with pytest.raises(Tier3ValidationError, match="annotation assembly"):
        validate_manifest(manifest, verify_local_files=False)


def test_manifest_rejects_variant_source_without_explicit_callability():
    manifest = manifest_fixture()
    del manifest["datasets"][0]["denominator"]
    with pytest.raises(Tier3ValidationError, match="denominator"):
        validate_manifest(manifest, verify_local_files=False)


def test_manifest_rejects_projected_annotation_for_primary_4d():
    manifest = manifest_fixture()
    manifest["datasets"][0]["annotation"]["status"] = "projected"
    with pytest.raises(Tier3ValidationError, match="native"):
        validate_manifest(manifest, verify_local_files=False)


def test_manifest_rejects_duplicate_dataset_ids_and_raw_paths():
    manifest = manifest_fixture()
    manifest["datasets"].append(copy.deepcopy(manifest["datasets"][0]))
    with pytest.raises(Tier3ValidationError, match="duplicate dataset_id"):
        validate_manifest(manifest, verify_local_files=False)
    manifest = manifest_fixture()
    manifest["datasets"][0]["notes"] = "credential token=secret"
    with pytest.raises(Tier3ValidationError, match="credential"):
        validate_manifest(manifest, verify_local_files=False)


def test_manifest_local_checksum_verification_uses_file_uris(tmp_path):
    manifest = manifest_fixture()
    artifact = manifest["datasets"][0]["reference"]["fasta"]
    artifact["uri"] = (tmp_path / "absent.fa").as_uri()
    with pytest.raises(Tier3ValidationError, match="does not exist"):
        validate_manifest(manifest, verify_local_files=True)


def test_buffalo_pin_and_all_173_core_species_are_enforced(tmp_path):
    source = FIXTURES / "buffalo_core_sample.tsv"
    with pytest.raises(Tier3ValidationError, match="expected 173"):
        buffalo_core_rows(source)
    # The immutable production pin is part of the public generator API.
    assert len(BUFFALO_SHA256) == 64
    assert BUFFALO_CORE_COUNT == 173
    provenance = json.loads((FIXTURES / "buffalo_core.expected.json").read_text())
    assert provenance["core_count"] == BUFFALO_CORE_COUNT
    assert provenance["source"]["sha256"] == BUFFALO_SHA256
    names = [row["scientific_name"] for row in provenance["species"]]
    assert len(names) == len(set(names)) == BUFFALO_CORE_COUNT


def test_impg_10kb_truth_is_complete_and_remains_fail_closed():
    truth = json.loads((FIXTURES / "impg_10kb_truth.json").read_text())
    assert truth["length_bp"] == 10_000
    assert truth["reference_pansn"] == "truth#1#chr1"
    assert truth["query_pansn"] == "truth#2#chr1"
    assert [variant["kind"] for variant in truth["variants"]] == [
        "SNP",
        "insertion",
        "deletion",
    ]
    assert all(variant["expected_count"] == 1 for variant in truth["variants"])
    # Frozen policy does not permit a nominal truth fixture to turn an
    # unlocked Cargo dependency graph into an eligible executable.
    assert truth["impg_execution_approved"] is False
    assert truth["phase_sensitive_eligible"] is False


def test_all_committed_tier3_schemas_are_valid_draft_2020_12():
    import jsonschema

    for path in sorted((Path(__file__).parents[1] / "schemas").glob("*.schema.json")):
        jsonschema.Draft202012Validator.check_schema(json.loads(path.read_text()))
