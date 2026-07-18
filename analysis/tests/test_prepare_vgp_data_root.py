from pathlib import Path

import pytest

from analysis import prepare_vgp_data_root as vgp_root


def _contract(tmp_path: Path) -> dict:
    root = tmp_path / "vgp-root"
    return {
        "task_id": "prepare-vgp-data",
        "root": str(root),
        "default_directory_mode": "2770",
        "layout": {
            "manifests": "manifests",
            "immutable_objects": "objects/sha256",
            "accession_views": "views/accession",
            "version_views": "views/version",
            "staging": "staging",
            "staging_acquisition": "staging/acquisition",
            "staging_partials": "staging/partials",
            "staging_outputs": "staging/outputs",
            "quarantine": "quarantine",
            "logs": "logs",
            "locks": "locks",
            "pilot": "pilot",
            "pilot_manifests": "pilot/manifests",
            "pilot_runs": "pilot/runs",
            "pilot_outputs": "pilot/outputs",
        },
        "transfer_contract": {
            "acquisition_staging_scope": "staging/acquisition",
            "acquisition_partial_scope": "staging/partials",
            "immutable_promotion_target": "objects/sha256",
            "accession_view_scope": "views/accession",
            "version_view_scope": "views/version",
            "compute_input_policy": "verified inputs only",
            "compute_workdir_policy": "copy to SLURM_TMPDIR",
            "compute_network_policy": "compute arrays never download",
            "output_validation_policy": "validate then atomically promote",
            "cleanup_scope_policy": "cleanup stays inside root",
        },
    }


def test_provision_layout_and_smoke_tests(tmp_path):
    contract = _contract(tmp_path)
    layout = vgp_root.build_layout(contract)
    records = vgp_root.provision_layout(layout)
    assert layout.root.is_dir()
    assert len(records) == len(contract["layout"])
    for record in records:
        assert Path(record["path"]).is_dir()
        assert record["mode_octal"] == "0o2770"

    smoke = vgp_root.smoke_test_storage_contract(layout)
    assert smoke["file_fsync"]["status"] == "pass"
    assert smoke["atomic_promotion"]["status"] == "pass"
    assert smoke["checksum_verification"]["status"] == "pass"
    assert smoke["lock_behavior"]["status"] == "pass"
    assert smoke["cleanup"]["status"] == "pass"
    assert not (layout.directories["immutable_objects"] / ".validation").exists()


def test_cleanup_guard_rejects_paths_outside_root(tmp_path):
    contract = _contract(tmp_path)
    layout = vgp_root.build_layout(contract)
    vgp_root.provision_layout(layout)
    outside = tmp_path / "elsewhere"
    outside.mkdir()
    with pytest.raises(ValueError):
        vgp_root.ensure_within_root(layout.root, outside)


def test_missing_quota_helper_is_nonblocking_but_atomic_promotion_is_hard():
    blockers = vgp_root.collect_blockers(
        system={
            "quota_state": {"status": "blocked"},
            "inode_state": {"status": "reported"},
        },
        smoke={"atomic_promotion": {"status": "blocked"}},
    )
    codes = {item["code"] for item in blockers}
    assert codes == {"ATOMIC_PROMOTION_UNVERIFIED"}
