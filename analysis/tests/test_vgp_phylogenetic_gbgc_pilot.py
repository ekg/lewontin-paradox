import copy
import csv
import importlib.util
from pathlib import Path


ROOT = Path(__file__).parents[2]
MODULE_PATH = ROOT / "analysis/run_vgp_phylogenetic_gbgc_pilot.py"
SPEC = importlib.util.spec_from_file_location("run_vgp_phylogenetic_gbgc_pilot", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
pilot = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(pilot)


def _write_tsv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def test_frozen_blocked_pilot_artifacts_validate():
    assert pilot.validate_artifacts() == []


def test_exact_h01_h02_accessions_and_roles_are_frozen():
    assert len(pilot.expected_accessions()) == 10
    for taxa in pilot.PANELS.values():
        roles = {taxon[1] for taxon in taxa}
        assert roles == {"focal_ingroup", "close_ingroup", "outgroup_1", "outgroup_2"}
        assert sum(taxon[1] == "close_ingroup" for taxon in taxa) == 2


def test_upstream_audit_proves_zero_verified_sequence():
    audit = pilot.audit_upstream()
    assert {row["accession_version"] for row in audit["source_sequences"]} == {
        "GCA_048126635.1", "GCA_048301445.1", "GCA_963210335.1", "GCA_964276395.1"
    }
    assert len(audit["source_sequences"]) == 8
    assert all(row["state"] == "planned" and row["observed_bytes"] == "0" for row in audit["mirror_sequences"])


def test_write_requires_pinned_guix_and_slurm(monkeypatch):
    monkeypatch.delenv("SLURM_JOB_ID", raising=False)
    monkeypatch.delenv("VGP_GBGC_PINNED_GUIX", raising=False)
    monkeypatch.delenv("GUIX_PROFILE", raising=False)
    errors = pilot.production_environment_errors()
    assert any("pinned channel" in error for error in errors)
    assert any("require Slurm" in error for error in errors)
    assert any("immutable Guix store path" in error for error in errors)


def test_validator_rejects_invented_biological_result(tmp_path):
    for name in (pilot.CLADE_OUTPUT, pilot.QC_OUTPUT, pilot.RESULT_OUTPUT, pilot.REPORT_OUTPUT):
        (tmp_path / name).write_bytes((pilot.ANALYSIS / name).read_bytes())
    _, rows = pilot.read_tsv(tmp_path / pilot.RESULT_OUTPUT)
    historical = next(row for row in rows if row["panel_id"] == "H01")
    historical["status"] = "MEASURED"
    historical["estimate_value"] = "0.42"
    _write_tsv(tmp_path / pilot.RESULT_OUTPUT, rows)
    errors = pilot.validate_artifacts(tmp_path)
    assert any("promoted to a biological result" in error for error in errors)


def test_validator_rejects_missing_control(tmp_path):
    for name in (pilot.CLADE_OUTPUT, pilot.QC_OUTPUT, pilot.RESULT_OUTPUT, pilot.REPORT_OUTPUT):
        (tmp_path / name).write_bytes((pilot.ANALYSIS / name).read_bytes())
    _, rows = pilot.read_tsv(tmp_path / pilot.QC_OUTPUT)
    rows = [row for row in rows if not (row["panel_id"] == "H02" and row["qc_id"] == "C02_NULL_SIMULATION")]
    _write_tsv(tmp_path / pilot.QC_OUTPUT, rows)
    errors = pilot.validate_artifacts(tmp_path)
    assert any("H02: missing QC/control rows" in error for error in errors)


def test_validator_rejects_local_digest_claim(tmp_path):
    for name in (pilot.CLADE_OUTPUT, pilot.QC_OUTPUT, pilot.RESULT_OUTPUT, pilot.REPORT_OUTPUT):
        (tmp_path / name).write_bytes((pilot.ANALYSIS / name).read_bytes())
    _, rows = pilot.read_tsv(tmp_path / pilot.CLADE_OUTPUT)
    rows = copy.deepcopy(rows)
    rows[0]["local_sequence_sha256"] = "0" * 64
    _write_tsv(tmp_path / pilot.CLADE_OUTPUT, rows)
    errors = pilot.validate_artifacts(tmp_path)
    assert any("invented local digest" in error for error in errors)
