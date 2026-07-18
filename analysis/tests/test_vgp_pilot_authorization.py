import copy
import gzip
import hashlib
import json
from pathlib import Path

import pytest

from analysis import vgp_pilot_authorization as auth
from analysis import migrate_vgp_data_root as migration


ROOT = Path(__file__).parents[2]


def test_production_authorization_reconciles_all_ten_and_selects_smallest_canary():
    value = auth.build_authorization()
    assert auth.validate_authorization(value) == {
        "authorized_pairs": 10, "canary": "P07", "optional_qc_vetoes": 0,
    }
    assert value["authorized_selection_ids"] == [f"P{i:02d}" for i in range(1, 11)]
    assert {row["core_authorization"] for row in value["pairs"]} == {"AUTHORIZED"}
    assert all(row["missing_confidence_covariates"] for row in value["pairs"])
    assert all(row["missing_covariates_block_core"] is False for row in value["pairs"])
    assert value["missing_optional_qc_global_job_count_effect"] == 0
    sizes = {
        row["selection_id"]: row["resources"]["measurement_basis"]["compressed_input_bytes"]
        for row in value["pairs"]
    }
    assert min(sizes, key=sizes.get) == "P07"


def test_committed_manifest_is_deterministic_and_has_real_generous_packets():
    committed = json.loads(auth.AUTHORIZATION.read_text())
    assert committed == auth.build_authorization()
    canary = next(row for row in committed["pairs"] if row["selection_id"] == "P07")
    whole = canary["resources"]["whole_canary_packet"]
    assert whole["initial_memory_gib"] >= 128
    assert whole["oom_retry_memory_gib"] == [256, 512]
    assert whole["node_local_scratch_bytes_high"] > canary["resources"]["measurement_basis"]["compressed_input_bytes"]
    for relative in committed["canary"]["job_packets"].values():
        text = (ROOT / relative).read_text()
        assert "#SBATCH --mem=" in text and "#SBATCH --partition=highmem" in text
        assert "run_canary.sh" in text
    worker = (auth.PACKET_DIR / "run_canary.sh").read_text()
    assert "SLURM_TMPDIR" in worker and "checkpoint" in worker


def test_active_authorization_and_execution_paths_use_canonical_shared_root():
    legacy = "/moosefs/erikg/lewontin-paradox-data/vgp"
    canonical = "/moosefs/erikg/vgp"
    value = json.loads(auth.AUTHORIZATION.read_text())
    for pair in value["pairs"]:
        for side in ("h1", "h2"):
            assert pair[side]["cas_path"].startswith(canonical + "/objects/sha256/")
            assert legacy not in pair[side]["cas_path"]
            assert Path(pair[side]["cas_path"]).resolve().is_relative_to(Path(canonical))
    active = [
        auth.AUTHORIZATION,
        auth.PREFLIGHT,
        ROOT / "analysis/vgp_data_root_validation.json",
        ROOT / "analysis/vgp_analysis_manifest.json",
        ROOT / "analysis/vgp_phase1_freeze_provenance.json",
        ROOT / "analysis/vgp_10_pilot_acquisition.py",
        ROOT / "analysis/mirror_vgp_freeze1.py",
        ROOT / "analysis/resolve_vgp_candidates.py",
        ROOT / "analysis/assert_vgp_comprehensive_design.py",
        ROOT / "analysis/vgp_core_scaleout.py",
        *auth.PACKET_DIR.glob("*"),
    ]
    assert all(legacy not in path.read_text() for path in active if path.is_file())
    config = json.loads((ROOT / "analysis/vgp_data_root_config.json").read_text())
    assert config["root"] == canonical
    assert config["migration_input_only"].startswith(legacy)


def test_applied_migration_reuses_verified_objects_without_redownload():
    report = json.loads(migration.OUTPUT.read_text())
    assert report["canonical_root"] == "/moosefs/erikg/vgp"
    assert report["legacy_root_role"] == "migration_input_only"
    assert report["unique_verified_objects"] == 140
    assert report["logical_verified_bytes"] == 12_692_829_704
    assert report["downloaded_bytes"] == report["network_requests"] == 0
    assert report["source_objects_removed"] == 0
    assert report["active_output_paths_under_legacy_root"] == 0
    assert all(row["canonical_path"].startswith("/moosefs/erikg/vgp/") for row in report["objects"])
    assert all(row["same_inode"] is True for row in report["objects"])


def test_missing_optional_qc_never_recreates_zero_jobs():
    value = auth.build_authorization()
    for pair in value["pairs"]:
        pair["missing_confidence_covariates"] = list(auth.CONFIDENCE_COVARIATES)
        pair["confidence_tier_at_authorization"] = "C"
    assert auth.validate_authorization(value)["authorized_pairs"] == 10
    assert value["authorized_selection_ids"] == [f"P{i:02d}" for i in range(1, 11)]


def test_preflight_fixture_checks_gzip_digest_and_storage_without_biology(tmp_path, monkeypatch):
    value = auth.build_authorization()
    fasta = tmp_path / "tiny.fa.gz"
    with gzip.open(fasta, "wt") as handle:
        handle.write(">fixture\nACGTACGTNN\n")
    digest = auth.sha256_file(fasta)
    for pair in value["pairs"]:
        for side in ("h1", "h2"):
            pair[side].update(cas_path=str(fasta), compressed_bytes=fasta.stat().st_size, sha256=digest)
        pair["resources"]["whole_canary_packet"]["node_local_scratch_bytes_high"] = 1
    monkeypatch.setattr(auth, "verify_environment_capture", lambda _: {"profile": "/gnu/store/fixture"})
    result = auth.preflight(value, tmp_path, run_scheduler_checks=False)
    assert result["global_authorization_result"] == "GO_10_OF_10"
    assert result["authorization_payload_sha256"] == hashlib.sha256(
        auth.canonical_json(value).encode("utf-8")
    ).hexdigest()
    assert result["biological_execution"] is False
    assert result["slurm_jobs_submitted"] == 0
    broken = copy.deepcopy(value)
    broken["pairs"][0]["h1"]["sha256"] = "0" * 64
    with pytest.raises(auth.AuthorizationError, match="digest gate"):
        auth.preflight(broken, tmp_path, run_scheduler_checks=False)


def test_shell_packets_are_syntax_valid():
    import subprocess
    for script in auth.PACKET_DIR.glob("*.sh"):
        subprocess.run(["bash", "-n", str(script)], check=True)
    for script in auth.PACKET_DIR.glob("*.sbatch"):
        subprocess.run(["bash", "-n", str(script)], check=True)
