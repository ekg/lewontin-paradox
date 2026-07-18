import csv
import hashlib
import json
from pathlib import Path

import pytest

from analysis import synthesize_vgp_comprehensive as synthesis


ROOT = Path(__file__).parents[2]


def rows(path):
    with Path(path).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


@pytest.fixture(scope="module")
def generated(tmp_path_factory):
    output = tmp_path_factory.mktemp("vgp-comprehensive")
    manifest = synthesis.generate(output)
    return output, manifest


def test_closed_world_reconciliation_is_exact_and_nonbiological(generated):
    _, manifest = generated
    assert manifest["closed_world"] == {
        "catalog_rows": 716,
        "catalog_unique_scientific_names": 714,
        "catalog_duplicate_name_rows": 2,
        "released_catalog_rows": 581,
        "unreleased_catalog_rows": 135,
        "mirror_inventory_objects": 47870,
        "mirror_verified_or_reused_objects": 0,
        "linked_haplotype_entries": 569,
        "distinct_nonself_pair_candidates": 566,
        "self_links_excluded": 3,
        "eligible_pairs": 0,
        "completed_pairs": 0,
        "exact_annotation_subset_pairs": 0,
        "biological_estimates": 0,
        "biological_jobs_submitted": 0,
    }
    assert manifest["core_pair_statuses"] == {
        "confidence_tiers": {"UNASSIGNED": 566, "X": 3},
        "dispositions": {"excluded": 3, "failed": 566},
        "primary_reason_codes": {
            "CATALOG_SELF_LINK_NOT_PAIR": 3,
            "UPSTREAM_SCALEOUT_NOT_AUTHORIZED": 566,
        },
    }
    assert manifest["catalog_row_statuses"] == {
        "dispositions": {"excluded": 187, "failed": 529},
        "primary_reason_codes": {
            "CATALOG_NO_LINKED_HAPLOTYPE": 49,
            "CATALOG_NO_RELEASED_MAIN_ACCESSION": 135,
            "CATALOG_ONLY_SELF_LINK_NOT_PAIR": 3,
            "UPSTREAM_SCALEOUT_NOT_AUTHORIZED": 529,
        },
    }


def test_same_pair_psmc_and_annotation_guards_are_machine_readable(generated):
    output, manifest = generated
    assert manifest["analysis_contract"]["same_pair_psmc_independent_evidence"] is False
    assert manifest["analysis_contract"]["annotation_absence_core_veto"] is False
    assert manifest["analysis_contract"]["cross_species_model_status"] == "NOT_IDENTIFIABLE_ZERO_ELIGIBLE_PAIRS"
    table = rows(output / synthesis.CORE_TABLE.name)
    by_id = {row["row_id"]: row for row in table}
    assert by_id["CORE-PSMC"]["independence"] == "SAME_PAIR_NONINDEPENDENT_DESCRIPTIVE"
    assert by_id["CORE-ANNOTATION"]["eligibility_role"] == "POST_CORE_OPTIONAL_PARTITION_ONLY"
    assert by_id["CORE-UNSCALED"]["status"] == "NOT_MATERIALIZED_CORE_NOT_RUN"
    assert by_id["CORE-SCENARIOS"]["status"] == "NOT_MATERIALIZED_NO_APPROVED_TIER"


def test_four_conversion_rows_remain_separate_and_absent_branches_are_design_only(generated):
    output, manifest = generated
    evidence = rows(output / synthesis.GENE_CONVERSION_TABLE.name)
    assert [row["branch"] for row in evidence] == [
        "direct_pedigree_or_gamete",
        "population_allele_frequency_spectrum",
        "historical_phylogenetic_substitution",
        "non_allelic_paralog",
    ]
    by_branch = {row["branch"]: row for row in evidence}
    assert by_branch["direct_pedigree_or_gamete"]["execution_state"] == "EXECUTED_PREFLIGHT_INPUT_GATE"
    assert by_branch["historical_phylogenetic_substitution"]["execution_state"] == "EXECUTED_PREFLIGHT_INPUT_GATE"
    for branch in ("population_allele_frequency_spectrum", "non_allelic_paralog"):
        assert by_branch[branch]["execution_state"] == "NOT_RUN_DESIGN_ONLY"
        assert by_branch[branch]["claim_classification"] == "design-only"
        assert by_branch[branch]["estimate"] == "NOT_ESTIMABLE"
    assert manifest["gene_conversion"]["direct"]["validated_events"] == 0
    assert manifest["gene_conversion"]["direct"]["published_candidates_excluded"] == 58
    assert manifest["gene_conversion"]["phylogenetic"]["verified_sequence_bytes"] == 0
    assert manifest["gene_conversion"]["population"]["state"] == "NOT_RUN_DESIGN_ONLY"
    assert manifest["gene_conversion"]["non_allelic"]["state"] == "NOT_RUN_DESIGN_ONLY"


def test_every_claim_has_classification_evidence_lineage_and_forbidden_extrapolation(generated):
    output, manifest = generated
    claims = rows(output / synthesis.CLAIM_LEDGER.name)
    allowed = {"supported", "bounded", "suggestive", "design-only", "not identifiable"}
    assert len(claims) >= 12
    assert {row["classification"] for row in claims} <= allowed
    assert allowed - {"suggestive"} <= {row["classification"] for row in claims}
    for row in claims:
        assert row["evidence_artifacts"]
        assert row["sampling_unit"]
        assert row["uncertainty_and_covariance"]
        assert row["data_lineage"]
        assert row["forbidden_extrapolation"]
        bound_names = {Path(path).name for path in manifest["input_digests"]}
        bound_names |= {Path(path).name for path in manifest["output_digests"]}
        assert set(row["evidence_artifacts"].split(";")) <= bound_names


def test_outputs_are_atomic_digest_bound_and_committed_outputs_reproduce(generated):
    output, manifest = generated
    for relative, digest in manifest["output_digests"].items():
        assert synthesis.sha256_file(output / Path(relative).name) == digest
    assert not list(output.glob(".*.partial-*"))
    assert synthesis.validate_outputs(output) == []

    committed = json.loads((ROOT / "analysis/vgp_comprehensive_final_manifest.json").read_text())
    assert committed == manifest
    for relative, digest in committed["output_digests"].items():
        assert hashlib.sha256((ROOT / relative).read_bytes()).hexdigest() == digest


def test_validator_rejects_promotion_of_population_or_nonallelic_branch(generated, tmp_path):
    output, _ = generated
    for filename in synthesis.OUTPUT_FILENAMES:
        (tmp_path / filename).write_bytes((output / filename).read_bytes())
    evidence = rows(tmp_path / synthesis.GENE_CONVERSION_TABLE.name)
    population = next(row for row in evidence if row["branch"] == "population_allele_frequency_spectrum")
    population["execution_state"] = "MEASURED"
    population["estimate"] = "0.1"
    with (tmp_path / synthesis.GENE_CONVERSION_TABLE.name).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(evidence[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(evidence)
    errors = synthesis.validate_outputs(tmp_path, verify_digests=False)
    assert any("population_allele_frequency_spectrum must remain NOT_RUN_DESIGN_ONLY" in error for error in errors)


def test_paper_ready_svgs_are_static_and_self_describing(generated):
    output, _ = generated
    for name in (synthesis.CLOSED_WORLD_FIGURE.name, synthesis.EVIDENCE_FIGURE.name):
        text = (output / name).read_text(encoding="utf-8")
        assert text.startswith("<svg")
        assert "<title>" in text
        assert "<desc>" in text
        assert "NOT_RUN / DESIGN_ONLY" in text or "0 eligible pairs" in text
