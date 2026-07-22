import csv
import json
from pathlib import Path

from analysis import build_vgp_repair_base as repair_base


def test_corrective_repair_base_is_complete_and_job_free():
    summary = repair_base.check()
    assert summary["source_commits"] >= 90
    assert summary["selected_patches"] > summary["skipped_conclusion_patches"]
    assert summary["skipped_conclusion_patches"] >= 4
    assert summary["status_bound_artifacts"] >= 25
    assert summary["submitted_biological_jobs"] == 0


def test_claim_boundaries_and_historical_statuses_are_machine_readable():
    manifest = json.loads(repair_base.MANIFEST.read_text())
    forbidden = set(manifest["selection_policy"]["forbidden_inferences"])
    assert "P07 is assembly-invalid" in forbidden
    assert "technical pipeline failure is a biological exclusion" in forbidden
    assert "two completed result packets constitute full VGP scale-out" in forbidden

    rows = list(csv.DictReader(repair_base.STATUS_LEDGER.open(), delimiter="\t"))
    by_path = {row["path"]: row for row in rows}
    assert by_path["analysis/vgp_read_validation_report_v1.md"]["status"] == "SUPERSEDED"
    assert by_path["analysis/vgp_real_scaleout_v1/results.md"]["status"] == "SUPERSEDED"
    assert by_path["analysis/vgp_real_synthesis_v1/report.md"]["status"] == "SUPERSEDED"
    assert by_path["analysis/vgp_freeze1_mirror_manifest.tsv"]["status"] == "HISTORICAL"
    assert all(len(row["sha256"]) == 64 for row in rows)


def test_required_paths_are_repository_relative_and_present():
    manifest = json.loads(repair_base.MANIFEST.read_text())
    for path_text in manifest["required_contracts"]:
        path = Path(path_text)
        assert not path.is_absolute()
        assert (repair_base.ROOT / path).is_file()
