import csv
import hashlib
import json
import stat
from pathlib import Path

import pytest

from analysis import acquire_vgp_pilot as acquire
from analysis import gate_vgp_pilot as gate
from analysis.tier3_common import Tier3ValidationError


ROOT = Path(__file__).parents[2]


def _read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _outputs(tmp_path: Path) -> dict[str, Path]:
    return {
        "output_manifest_path": tmp_path / "acquisition_manifest.tsv",
        "output_report_path": tmp_path / "acquisition_report.md",
        "output_inventory_path": tmp_path / "immutable_inventory.tsv",
        "refusal_evidence_path": tmp_path / "refusal.json",
    }


def _write_rehashed_gate(tmp_path: Path, mutate) -> Path:
    payload = json.loads((ROOT / "analysis/vgp_pilot_gate.json").read_text(encoding="utf-8"))
    mutate(payload)
    payload["decision_sha256"] = gate.sha256_json(
        {key: value for key, value in payload.items() if key != "decision_sha256"}
    )
    path = tmp_path / "altered_gate.json"
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path


def _run_with_spy(tmp_path: Path, gate_path: Path):
    calls = []

    def downloader(*args):
        calls.append(args)
        raise AssertionError("downloader crossed a refused preflight boundary")

    result = acquire.run(
        gate_path=gate_path,
        manifest_path=ROOT / "analysis/vgp_pilot_manifest.tsv",
        root_config_path=ROOT / "analysis/vgp_data_root_config.json",
        downloader=downloader,
        **_outputs(tmp_path),
    )
    return result, calls


def test_run_writes_refusal_manifest_and_report_for_current_nogo_gate(tmp_path):
    outputs = _outputs(tmp_path)
    result = acquire.run(
        gate_path=ROOT / "analysis/vgp_pilot_gate.json",
        manifest_path=ROOT / "analysis/vgp_pilot_manifest.tsv",
        root_config_path=ROOT / "analysis/vgp_data_root_config.json",
        **outputs,
    )
    assert result["status"] == "refused_preflight"
    # The independently regenerated gate binds the repaired manifest and then
    # refuses because the strict storage/cap contract remains NO_GO.
    assert result["failure_code"] == "GATE_NO_GO"
    rows = _read_manifest(outputs["output_manifest_path"])
    assert rows[0]["record_type"] == "run_summary"
    assert rows[0]["status"] == "refused_preflight"
    assert rows[0]["observed_bytes"] == "0"
    blocker_codes = {row["failure_code"] for row in rows[1:]}
    assert "QUOTA_UNAVAILABLE" in blocker_codes
    assert "CAP_SCRATCH_GIB_EXCEEDED" in blocker_codes
    report = outputs["output_report_path"].read_text(encoding="utf-8")
    assert "Gate decision: `NO_GO`" in report
    assert "Refused before first biological byte: `true`" in report
    assert "Exact-reference/native-annotation linkage validation under pinned GNU Guix was not re-run" in report
    inventory = _read_manifest(outputs["output_inventory_path"])
    assert inventory == []
    evidence = json.loads(outputs["refusal_evidence_path"].read_text(encoding="utf-8"))
    assert evidence["provider_requests_attempted"] == 0
    assert evidence["biological_payload_bytes_transferred"] == 0
    assert evidence["partial_files_created"] == 0
    assert evidence["objects_promoted"] == 0


def test_run_refuses_tampered_gate_and_records_gate_tampered_failure(tmp_path):
    gate_out = tmp_path / "tampered_gate.json"
    gate_payload = json.loads((ROOT / "analysis/vgp_pilot_gate.json").read_text(encoding="utf-8"))
    gate_payload["decision"]["status"] = "GO"
    gate_out.write_text(json.dumps(gate_payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    outputs = _outputs(tmp_path)
    result = acquire.run(
        gate_path=gate_out,
        manifest_path=ROOT / "analysis/vgp_pilot_manifest.tsv",
        root_config_path=ROOT / "analysis/vgp_data_root_config.json",
        **outputs,
    )
    assert result["status"] == "refused_preflight"
    assert result["failure_code"] == "GATE_TAMPERED"
    rows = _read_manifest(outputs["output_manifest_path"])
    assert rows[0]["failure_code"] == "GATE_TAMPERED"
    assert "gate decision hash does not match the gate payload" in rows[0]["failure_message"]


@pytest.mark.parametrize(
    ("name", "mutate", "expected_code"),
    [
        (
            "unknown_decision",
            lambda payload: payload["decision"].__setitem__("status", "UNKNOWN"),
            "GATE_DECISION_UNKNOWN",
        ),
        (
            "altered_bound_digest",
            lambda payload: payload["authorization_boundary"].__setitem__("manifest_digest", "0" * 64),
            "BOUND_DIGEST_MISMATCH",
        ),
        (
            "altered_approved_url",
            lambda payload: payload["retrieval_audit"]["rows"][0]["obligations"][0].__setitem__(
                "url", "https://ftp.ncbi.nlm.nih.gov/not-approved/GCF_036321535.1_genomic.fna.gz"
            ),
            "RETRIEVAL_CONTRACT_MISMATCH",
        ),
        (
            "altered_accession_version",
            lambda payload: payload["retrieval_audit"]["rows"][0]["obligations"][0].__setitem__(
                "accession_version", "GCF_036321535.2"
            ),
            "RETRIEVAL_CONTRACT_MISMATCH",
        ),
        (
            "relaxed_cap",
            lambda payload: payload["cap_vector"]["dimensions"]["species"].__setitem__("limit", 7.0),
            "CAP_VECTOR_MISMATCH",
        ),
        (
            "relaxed_storage",
            lambda payload: payload["storage_audit"].__setitem__("adequate", True),
            "ROOT_STORAGE_CONTRACT_MISMATCH",
        ),
    ],
)
def test_run_refuses_gate_contract_mutations_before_downloader(
    tmp_path, name, mutate, expected_code
):
    gate_path = _write_rehashed_gate(tmp_path, mutate)
    result, calls = _run_with_spy(tmp_path, gate_path)
    assert result["status"] == "refused_preflight", name
    assert result["failure_code"] == expected_code, name
    assert result["transferred_bytes"] == 0
    assert result["provider_requests_attempted"] == 0
    assert calls == []
    assert not list(tmp_path.rglob("*.part"))
    inventory = _read_manifest(tmp_path / "immutable_inventory.tsv")
    assert inventory == []


def test_run_refuses_missing_gate_before_downloader(tmp_path):
    result, calls = _run_with_spy(tmp_path, tmp_path / "missing.json")
    assert result["status"] == "refused_preflight"
    assert result["failure_code"] == "BOUND_INPUT_MISSING"
    assert result["transferred_bytes"] == 0
    assert calls == []


def test_live_storage_identity_is_cross_host_portable_but_same_filesystem_bound():
    payload = json.loads((ROOT / "analysis/vgp_pilot_gate.json").read_text(encoding="utf-8"))
    storage = payload["storage_audit"]
    storage["adequate"] = True
    storage["enforceable_allocation"]["status"] = "known"
    storage["enforceable_allocation"]["headroom_pass"] = True
    root = Path(storage["root"])

    # Numeric st_dev is intentionally absent because it varies by host mount
    # namespace; the stable path/inode/ownership/mode and catalog co-location
    # remain exact live gates.
    assert "device" not in storage["live_identity"]
    acquire._validate_live_storage(payload, root)

    storage["live_identity"]["catalog_on_same_filesystem"] = False
    with pytest.raises(Tier3ValidationError, match="catalog filesystem identity is not bound"):
        acquire._validate_live_storage(payload, root)


def test_promotion_rejects_staged_content_change_before_atomic_promotion(tmp_path):
    part = tmp_path / "staging" / "payload.part"
    part.parent.mkdir()
    part.write_bytes(b"approved-payload")
    object_root = tmp_path / "objects" / "sha256"

    def mutate(path: Path, first_sha256: str) -> None:
        assert first_sha256 == hashlib.sha256(b"approved-payload").hexdigest()
        path.write_bytes(b"altered-payload!")

    with pytest.raises(Tier3ValidationError, match="changed before promotion"):
        acquire.promote_verified_part(
            part,
            object_root,
            expected_size=len(b"approved-payload"),
            before_reverify=mutate,
        )
    assert part.read_bytes() == b"altered-payload!"
    assert not object_root.exists()


def test_promotion_verifies_official_md5_and_creates_read_only_sha256_object(tmp_path):
    content = b"small deterministic fixture"
    part = tmp_path / "staging" / "payload.part"
    part.parent.mkdir()
    part.write_bytes(content)
    object_root = tmp_path / "objects" / "sha256"
    result = acquire.promote_verified_part(
        part,
        object_root,
        expected_size=len(content),
        expected_source_checksum={"algorithm": "md5", "value": hashlib.md5(content).hexdigest()},
    )
    destination = Path(result["object_path"])
    assert destination.read_bytes() == content
    assert result["local_sha256"] == hashlib.sha256(content).hexdigest()
    assert result["source_checksum_verified"] is True
    assert stat.S_IMODE(destination.stat().st_mode) == (stat.S_IRUSR | stat.S_IRGRP)
    assert not part.exists()


def test_promotion_rejects_official_checksum_mismatch_without_promotion(tmp_path):
    part = tmp_path / "payload.part"
    part.write_bytes(b"payload")
    object_root = tmp_path / "objects"
    with pytest.raises(Tier3ValidationError, match="official MD5 mismatch"):
        acquire.promote_verified_part(
            part,
            object_root,
            expected_size=7,
            expected_source_checksum={"algorithm": "md5", "value": "0" * 32},
        )
    assert part.exists()
    assert not object_root.exists()
