from __future__ import annotations

import hashlib
import gzip
import json
from pathlib import Path

import pytest

from analysis import acquire_vgp_validation_reads as reads


def object_spec(payload: bytes) -> dict[str, object]:
    return {
        "object_id": "P07:SRRTEST:test.fastq.gz",
        "selection_id": "P07",
        "run_accession": "SRRTEST",
        "filename": "test.fastq.gz",
        "source_url": "https://example.invalid/test.fastq.gz",
        "expected_bytes": len(payload),
        "upstream_checksum_algorithm": "md5",
        "upstream_checksum": hashlib.md5(payload).hexdigest(),  # noqa: S324
    }


def test_acquire_resumes_verifies_promotes_and_reuses(tmp_path: Path) -> None:
    payload = gzip.compress(b"@read\nACGT\n+\nIIII\n" * 31, mtime=0)
    spec = object_spec(payload)
    partial = reads.partial_path(tmp_path, spec)
    partial.parent.mkdir(parents=True)
    partial.write_bytes(payload[:19])

    offsets: list[int] = []

    def downloader(_spec: dict[str, object], path: Path, offset: int) -> int:
        offsets.append(offset)
        with path.open("ab") as handle:
            handle.write(payload[offset:])
        return len(payload) - offset

    first = reads.acquire_object(spec, tmp_path, downloader=downloader)
    assert offsets == [19]
    assert first["status"] == "verified"
    assert first["resume_from_bytes"] == 19
    assert first["transferred_bytes"] == len(payload)
    assert first["invocation_transferred_bytes"] == len(payload) - 19
    target = Path(str(first["local_path"]))
    assert target.read_bytes() == payload
    assert target.name == hashlib.sha256(payload).hexdigest()
    assert target.stat().st_mode & 0o222 == 0
    view = Path(str(first["accession_view_path"]))
    assert view.is_symlink()
    assert view.resolve() == target.resolve()

    second = reads.acquire_object(spec, tmp_path, downloader=downloader)
    assert second["status"] == "reused"
    assert second["transferred_bytes"] == 0
    assert offsets == [19]


def test_checksum_mismatch_is_quarantined(tmp_path: Path) -> None:
    expected = b"expected"
    actual = b"wrong!!!"
    spec = object_spec(expected)

    def downloader(_spec: dict[str, object], path: Path, _offset: int) -> int:
        path.write_bytes(actual)
        return len(actual)

    result = reads.acquire_object(spec, tmp_path, downloader=downloader)
    assert result["status"] == "quarantined"
    assert not Path(str(result["local_path"] or "/missing")).exists()
    assert Path(str(result["quarantine_path"])).read_bytes() == actual


def test_recover_unmanifested_promotion_without_redownload(tmp_path: Path) -> None:
    payload = gzip.compress(b"@recovered\nACGT\n+\nIIII\n", mtime=0)
    spec = object_spec(payload)
    digest = hashlib.sha256(payload).hexdigest()
    target = reads.cas_path(tmp_path, digest)
    target.parent.mkdir(parents=True)
    target.write_bytes(payload)
    target.chmod(0o440)
    recovered = reads.recover_unmanifested_cas_sha(tmp_path, spec, "objects/sha256")
    assert recovered == digest

    def forbidden(*_args: object) -> int:
        raise AssertionError("recovery must not redownload")

    result = reads.acquire_object(
        spec,
        tmp_path,
        known_sha256=recovered,
        recovered_unmanifested_promotion=True,
        downloader=forbidden,
    )
    assert result["status"] == "verified"
    assert result["transferred_bytes"] == len(payload)
    assert result["resume_from_bytes"] == len(payload)
    assert Path(str(result["accession_view_path"])).resolve() == target.resolve()


def test_validate_plan_binds_every_object_to_exact_pair() -> None:
    plan = reads.load_plan()
    reads.validate_plan(plan)
    pairs = {row["selection_id"]: row for row in plan["pairs"]}
    assert set(pairs) == {"P04", "P07", "P09"}
    assert {row["selection_id"] for row in plan["objects"]} == set(pairs)
    for obj in plan["objects"]:
        pair = pairs[obj["selection_id"]]
        assert obj["biosample"] == pair["biosample"]
        assert obj["individual_or_isolate"] == pair["individual_or_isolate"]
        assert obj["h1_accession_version"] == pair["h1_accession_version"]
        assert obj["h2_accession_version"] == pair["h2_accession_version"]


def test_summary_accounts_every_object_and_byte() -> None:
    rows = [
        {"status": "verified", "expected_bytes": 10, "transferred_bytes": 10},
        {"status": "reused", "expected_bytes": 20, "transferred_bytes": 0},
        {"status": "missing", "expected_bytes": 30, "transferred_bytes": 4},
        {"status": "quarantined", "expected_bytes": 40, "transferred_bytes": 40},
    ]
    summary = reads.summarize(rows)
    assert summary["planned"] == {"objects": 4, "bytes": 100}
    assert summary["transferred"] == {"objects": 3, "bytes": 54}
    assert summary["verified"] == {"objects": 2, "bytes": 30}
    assert summary["reused"] == {"objects": 1, "bytes": 20}
    assert summary["missing"] == {"objects": 1, "bytes": 30}
    assert summary["quarantined"] == {"objects": 1, "bytes": 40}
    assert summary["accounting"]["objects_reconciled"] is True
    assert summary["accounting"]["bytes_reconciled"] is True


def test_root_is_loaded_from_repository_config(tmp_path: Path) -> None:
    config = tmp_path / "config.json"
    canonical = tmp_path / "canonical"
    config.write_text(
        json.dumps({"root": str(canonical), "layout": {"raw_reads": "raw/reads"}}),
        encoding="utf-8",
    )
    assert reads.configured_root(config) == canonical
    with pytest.raises(reads.AcquisitionError, match="absolute"):
        config.write_text(
            json.dumps({"root": "relative", "layout": {"raw_reads": "raw/reads"}}),
            encoding="utf-8",
        )
        reads.configured_root(config)


def test_verify_manifest_rehashes_payload_view_and_accounting(tmp_path: Path) -> None:
    payload = gzip.compress(b"@r\nAC\n+\nII\n", mtime=0)
    root = tmp_path / "vgp"
    config = tmp_path / "config.json"
    config.write_text(
        json.dumps(
            {
                "root": str(root),
                "layout": {
                    "immutable_objects": "objects/sha256",
                    "raw_reads": "raw/reads",
                    "pilot_manifests": "pilot/manifests",
                },
            }
        ),
        encoding="utf-8",
    )
    spec = object_spec(payload)

    def downloader(_spec: dict[str, object], path: Path, _offset: int) -> int:
        path.write_bytes(payload)
        return len(payload)

    row = reads.acquire_object(spec, root, downloader=downloader)
    manifest = {
        "schema_version": "vgp-validation-reads-manifest-v1.0.0",
        "canonical_root": str(root),
        "objects": [row],
        "summary": reads.summarize([row]),
    }
    manifest_path = tmp_path / reads.CANONICAL_MANIFEST_NAME
    canonical = root / "pilot/manifests" / reads.CANONICAL_MANIFEST_NAME
    reads.atomic_json(manifest_path, manifest)
    reads.atomic_json(canonical, manifest)
    result = reads.verify_manifest(manifest_path, config)
    assert result["verified"] is True
    assert result["verified_objects"] == 1
    assert result["verified_bytes"] == len(payload)
