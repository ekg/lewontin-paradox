import copy
import csv
import importlib.util
from pathlib import Path


ROOT = Path(__file__).parents[2]
MODULE_PATH = ROOT / "analysis/run_direct_gene_conversion_pilot.py"
SPEC = importlib.util.spec_from_file_location("run_direct_gene_conversion_pilot", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
pilot = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(pilot)


def _write_tsv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _copy_artifacts(tmp_path):
    for name in (pilot.DATASET_OUTPUT, pilot.TRACT_OUTPUT, pilot.SUMMARY_OUTPUT, pilot.REPORT_OUTPUT):
        (tmp_path / name).write_bytes((pilot.ANALYSIS / name).read_bytes())


def test_current_fail_closed_artifacts_validate():
    assert pilot.validate_artifacts() == []


def test_input_audit_reconciles_relationship_reads_reference_and_candidates():
    audit = pilot.audit_inputs()
    assert len(audit["runs"]) == 472
    assert audit["fastq_objects"] == 951
    assert audit["fastq_bytes"] == 356_984_558_868
    assert audit["bounded_objects"] == 697
    assert audit["bounded_bytes"] == 265_809_761_864
    assert len([
        row for row in audit["sample_rows"]
        if row["Sample type"] in {"Tetrad", "Tetrad (shallow sequenced)"}
    ]) == 52
    assert len(audit["marker_rows"]) == 137_339
    assert audit["nuclear_bases"] == 119_146_348
    assert len(audit["co_rows"]) == 71
    assert len({row["NCO ID"] for row in audit["nco_marker_rows"]}) == 14


def test_write_requires_exact_pinned_guix_profile_and_slurm(monkeypatch):
    monkeypatch.delenv("SLURM_JOB_ID", raising=False)
    monkeypatch.delenv("DIRECT_GC_PINNED_GUIX", raising=False)
    monkeypatch.delenv("GUIX_PROFILE", raising=False)
    errors = pilot.production_environment_errors()
    assert any("pinned channel" in error for error in errors)
    assert any("require Slurm" in error for error in errors)
    assert any("immutable Guix store profile" in error for error in errors)


def test_all_published_candidates_are_excluded_and_interval_censored():
    audit = pilot.audit_inputs()
    tracts = pilot.tract_rows(audit, "TEST_JOB")
    assert len(tracts) == 58
    assert sum(row["association_class"] == "CROSSOVER_ASSOCIATED" for row in tracts) == 44
    assert sum(row["association_class"] == "NON_CROSSOVER" for row in tracts) == 14
    assert all(row["analysis_status"] == pilot.PUBLISHED_ONLY for row in tracts)
    assert all(row["direct_rate_inclusion"] == "EXCLUDED" for row in tracts)
    assert all(0 < int(row["inner_tract_bp"]) <= int(row["outer_tract_bp"]) for row in tracts)


def test_validator_rejects_promoted_candidate(tmp_path):
    _copy_artifacts(tmp_path)
    _, rows = pilot.read_tsv(tmp_path / pilot.TRACT_OUTPUT)
    rows = copy.deepcopy(rows)
    rows[0]["analysis_status"] = "VALIDATED_DIRECT_EVENT"
    rows[0]["direct_rate_inclusion"] = "INCLUDED"
    _write_tsv(tmp_path / pilot.TRACT_OUTPUT, rows)
    assert any("promoted" in error for error in pilot.validate_artifacts(tmp_path))


def test_validator_rejects_numeric_blocked_rate(tmp_path):
    _copy_artifacts(tmp_path)
    _, rows = pilot.read_tsv(tmp_path / pilot.SUMMARY_OUTPUT)
    row = next(row for row in rows if row["estimand"] == "D_EVT_RATE_PER_BASE")
    row["estimate"] = "1e-9"
    row["ci_lower"] = "1e-10"
    row["ci_upper"] = "1e-8"
    _write_tsv(tmp_path / pilot.SUMMARY_OUTPUT, rows)
    assert any("numeric estimate" in error for error in pilot.validate_artifacts(tmp_path))


def test_validator_rejects_alternate_activation(tmp_path):
    _copy_artifacts(tmp_path)
    _, rows = pilot.read_tsv(tmp_path / pilot.DATASET_OUTPUT)
    d02 = next(row for row in rows if row["dataset_id"] == "D02")
    d02["activation_status"] = "ACTIVATED"
    _write_tsv(tmp_path / pilot.DATASET_OUTPUT, rows)
    assert any("without an amendment" in error for error in pilot.validate_artifacts(tmp_path))


def test_out_of_scope_estimands_are_not_measured():
    _, rows = pilot.read_tsv(pilot.ANALYSIS / pilot.SUMMARY_OUTPUT)
    for name in (
        "POPULATION_FREQUENCY_GBGC",
        "HISTORICAL_PHYLOGENETIC_GBGC",
        "NON_ALLELIC_CONVERSION",
        "CROSS_VERTEBRATE_TRANSFER",
    ):
        selected = [row for row in rows if row["estimand"] == name]
        assert len(selected) == 1
        assert selected[0]["status"] == "NOT_MEASURED"
        assert selected[0]["estimate"] == "NOT_MEASURED"
