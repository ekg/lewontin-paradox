import copy
import csv
import importlib.util
from pathlib import Path


ROOT = Path(__file__).parents[2]
MODULE_PATH = ROOT / "analysis/assert_gene_conversion_evidence_design.py"
SPEC = importlib.util.spec_from_file_location("assert_gene_conversion_evidence_design", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
design = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(design)


def test_frozen_gene_conversion_design_passes():
    assert design.validate_design() == []


def test_exact_executable_and_design_only_branch_states():
    _, rows = design.read_tsv(design.DATASETS)
    executable = {
        row["dataset_id"]
        for row in rows
        if row["execution_state"] == "AUTHORIZED_DOWNSTREAM_PILOT"
    }
    assert executable == {"D01", "H01", "H02"}
    assert {
        row["dataset_id"]
        for row in rows
        if row["branch"] in {"population", "non_allelic"}
    } == {"P01", "P02", "N01", "N02"}
    assert all(
        row["execution_state"] == "NOT_RUN_DESIGN_ONLY"
        for row in rows
        if row["branch"] in {"population", "non_allelic"}
    )


def test_all_four_noninterchangeable_estimand_branches_present():
    header, rows = design.read_tsv(design.ESTIMANDS)
    assert {row["branch"] for row in rows} == design.BRANCHES
    assert {
        "callable_opportunity",
        "primary_model_and_null_simulation",
        "polarization_error_control",
        "mutation_bias_and_demography_control",
        "linkage_recombination_or_phylogeny_control",
        "multiple_testing",
    }.issubset(header)


def test_claim_matrix_forbids_cross_branch_substitution():
    _, rows = design.read_tsv(design.CLAIMS)
    text = "\n".join("\t".join(row.values()) for row in rows)
    assert "population B to event rate" in text
    assert "Paralog homogenization is not allelic gBGC" in text
    assert "without parents/gametes, population frequencies, or outgroups" in text
    assert "Population-frequency gBGC was not measured" in text
    assert "Non-allelic conversion was not measured" in text


def test_dataset_validation_rejects_unauthorized_population_execution(tmp_path, monkeypatch):
    header, rows = design.read_tsv(design.DATASETS)
    rows = copy.deepcopy(rows)
    next(row for row in rows if row["dataset_id"] == "P01")["execution_state"] = "AUTHORIZED_DOWNSTREAM_PILOT"
    mutated = tmp_path / "datasets.tsv"
    with mutated.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    monkeypatch.setattr(design, "DATASETS", mutated)
    errors = design.validate_design()
    assert any("unauthorized design-only branch state" in error for error in errors)


def test_dataset_validation_rejects_missing_exact_accession(tmp_path, monkeypatch):
    header, rows = design.read_tsv(design.DATASETS)
    rows = copy.deepcopy(rows)
    h01 = next(row for row in rows if row["dataset_id"] == "H01")
    h01["exact_accessions"] = h01["exact_accessions"].replace("GCA_048126635.1", "UNVERSIONED")
    mutated = tmp_path / "datasets.tsv"
    with mutated.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    monkeypatch.setattr(design, "DATASETS", mutated)
    errors = design.validate_design()
    assert any(
        "H01: missing accession GCA_048126635.1 from exact_accessions" in error
        for error in errors
    )
