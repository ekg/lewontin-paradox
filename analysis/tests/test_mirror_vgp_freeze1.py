import csv
import hashlib
import json
import sqlite3
from pathlib import Path

import pytest

from analysis import mirror_vgp_freeze1 as mirror


ROOT = Path(__file__).parents[2]


def _catalog_rows():
    rows = []
    for number in range(581):
        accession = f"GCA_{number:09d}.1"
        rows.append(mirror.CatalogRow(number + 2, f"species {number}", accession))
    rows.extend(
        mirror.CatalogRow(number + 583, f"unreleased {number}", None) for number in range(135)
    )
    return rows


def _raw_inventory(tmp_path: Path, *, omit_last=False, add_file=True):
    lines = []
    for row in _catalog_rows()[: 580 if omit_last else 581]:
        root = mirror.accession_path(row.accession)
        lines.append(f"cd+++++++++|1|2024/01/02-03:04:05|{root}/|")
        if add_file:
            lines.append(
                f">f+++++++++|7|2024/01/02-03:04:05|{root}/{row.accession}.fa.gz|"
            )
    path = tmp_path / "inventory.items"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def test_pinned_catalog_identity_and_716_row_dispositions():
    rows = mirror.verify_catalog(mirror.PINNED_CATALOG)
    assert len(rows) == 716
    assert sum(row.accession is not None for row in rows) == 581
    assert sum(row.accession is None for row in rows) == 135
    assert len({row.accession for row in rows if row.accession}) == 581


def test_accession_path_preserves_exact_version():
    assert mirror.accession_path("GCF_024166365.1") == "GCF/024/166/365/GCF_024166365.1"
    with pytest.raises(mirror.MirrorError, match="invalid accession"):
        mirror.accession_path("GCF_024166365")


def test_inventory_refuses_outside_pinned_guix_before_network(tmp_path, monkeypatch):
    monkeypatch.delenv("GUIX_ENVIRONMENT", raising=False)
    with pytest.raises(mirror.MirrorError, match="pinned GNU Guix"):
        mirror.freeze_inventory(mirror.PINNED_CATALOG, tmp_path)


def test_parse_inventory_proves_all_581_roots_and_classifies_fasta(tmp_path):
    objects = mirror.parse_rsync_inventory(_raw_inventory(tmp_path), _catalog_rows())
    assert len(objects) == 1162
    assert len({obj.accession for obj in objects}) == 581
    assert sum(obj.sequence_subset == "assembly_fasta" for obj in objects) == 581
    assert sum(obj.size for obj in objects if obj.object_type == "file") == 581 * 7


def test_parse_inventory_fails_closed_on_missing_accession_root(tmp_path):
    with pytest.raises(mirror.MirrorError, match="frozen accessions are missing"):
        mirror.parse_rsync_inventory(
            _raw_inventory(tmp_path, omit_last=True),
            _catalog_rows(),
        )


def test_parse_inventory_rejects_duplicate_release_object(tmp_path):
    raw = _raw_inventory(tmp_path)
    first = raw.read_text().splitlines()[1]
    raw.write_text(raw.read_text() + first + "\n")
    with pytest.raises(mirror.MirrorError, match="duplicate release paths"):
        mirror.parse_rsync_inventory(raw, _catalog_rows())


def test_source_rows_emit_exact_endpoint_and_iso_mtime(tmp_path):
    objects = mirror.parse_rsync_inventory(_raw_inventory(tmp_path), _catalog_rows())
    metadata = {
        "retrieval_started_utc": "2026-01-01T00:00:00Z",
        "retrieval_completed_utc": "2026-01-01T00:01:00Z",
    }
    row = next(iter(mirror.source_rows(objects, metadata)))
    assert row["source_endpoint"].startswith(mirror.TRANSPORT_ENDPOINT + "GCA/")
    assert row["source_mtime_utc"] == "2024-01-02T03:04:05Z"


def test_promote_reverifies_sha256_and_is_atomic(tmp_path):
    content = b"verified fixture"
    part = tmp_path / "staging/payload.part"
    destination = tmp_path / "objects/payload"
    part.parent.mkdir()
    part.write_bytes(content)
    observed = mirror.promote_verified(part, destination, expected_size=len(content))
    assert observed == hashlib.sha256(content).hexdigest()
    assert destination.read_bytes() == content
    assert not part.exists()


def test_promote_detects_change_before_promotion_without_overwriting(tmp_path):
    part = tmp_path / "payload.part"
    destination = tmp_path / "durable"
    part.write_bytes(b"original")

    def alter(path, _digest):
        path.write_bytes(b"mutated!")

    with pytest.raises(mirror.MirrorError, match="changed before promotion"):
        mirror.promote_verified(
            part,
            destination,
            expected_size=8,
            before_reverify=alter,
        )
    assert not destination.exists()


def test_checksum_failure_is_quarantined_and_verified_destination_survives(tmp_path):
    durable = tmp_path / "objects/payload"
    durable.parent.mkdir()
    durable.write_bytes(b"last verified")
    part = tmp_path / "staging/payload.part"
    part.parent.mkdir()
    part.write_bytes(b"bad incoming")
    with pytest.raises(mirror.MirrorError, match="upstream MD5 mismatch"):
        mirror.promote_verified(
            part,
            durable,
            expected_size=len(b"bad incoming"),
            checksum_algorithm="md5",
            checksum="0" * 32,
        )
    quarantined = mirror.quarantine_part(part, tmp_path / "quarantine", "payload", "mismatch")
    assert quarantined.read_bytes() == b"bad incoming"
    assert durable.read_bytes() == b"last verified"


def test_state_accounting_is_mutually_exclusive_and_byte_exact():
    rows = [
        {"state": "planned", "object_type": "file", "size": 11},
        {"state": "verified", "object_type": "file", "size": 13},
        {"state": "verified", "object_type": "directory", "size": 99},
    ]
    accounting = mirror.state_accounting(rows)
    assert accounting["planned"] == {"objects": 1, "files": 1, "bytes": 11}
    assert accounting["verified"] == {"objects": 2, "files": 1, "bytes": 13}
    assert sum(value["objects"] for value in accounting.values()) == 3


def test_worker_refuses_before_rsync_when_capacity_gate_is_closed(tmp_path, monkeypatch):
    capacity = tmp_path / "capacity.json"
    capacity.write_text(json.dumps({"adequate": False, "gate_reason": "quota unavailable"}))
    monkeypatch.setenv("GUIX_ENVIRONMENT", "/gnu/store/test")
    with pytest.raises(mirror.MirrorError, match="capacity gate refused"):
        mirror.run_worker(tmp_path / "state.sqlite", tmp_path, 2, capacity)


def test_worker_rejects_slurm_environment_even_with_guix_marker(tmp_path, monkeypatch):
    capacity = tmp_path / "capacity.json"
    capacity.write_text(json.dumps({"adequate": True, "gate_reason": "verified"}))
    monkeypatch.setenv("GUIX_ENVIRONMENT", "/gnu/store/test")
    monkeypatch.setenv("SLURM_JOB_ID", "123")
    with pytest.raises(mirror.MirrorError, match="must not run as a Slurm"):
        mirror.run_worker(tmp_path / "state.sqlite", tmp_path, 2, capacity)


def test_upstream_md5_manifest_cannot_name_extra_uninventoried_object(tmp_path):
    accession = "GCA_000000000.1"
    root = mirror.accession_path(accession)
    objects_root = tmp_path / "objects"
    manifest = objects_root / root / "md5sum.txt"
    manifest.parent.mkdir(parents=True)
    manifest.write_text("0" * 32 + "  renamed-extra.fa.gz\n")
    objects = [
        mirror.InventoryObject(
            "id", 2, accession, root + "/md5sum.txt", "file", 1,
            "2024/01/02-03:04:05", "", "non_sequence_product_or_metadata"
        )
    ]
    with pytest.raises(mirror.MirrorError, match="non-inventoried object"):
        mirror.parse_upstream_md5s(objects_root, objects)


def test_database_update_accepts_only_exclusive_valid_state(tmp_path):
    database = tmp_path / "state.sqlite"
    connection = sqlite3.connect(database)
    connection.execute(
        "CREATE TABLE objects (inventory_id TEXT PRIMARY KEY, state TEXT, updated_at_utc TEXT)"
    )
    connection.execute("INSERT INTO objects VALUES ('x', 'planned', '')")
    connection.commit()
    connection.close()
    mirror.update_database(database, "x", state="verified")
    connection = sqlite3.connect(database)
    assert connection.execute("SELECT state FROM objects").fetchone()[0] == "verified"
    connection.close()
    with pytest.raises(mirror.MirrorError, match="invalid state"):
        mirror.update_database(database, "x", state="half_verified")


def test_committed_mirror_artifacts_are_closed_world_and_byte_exact():
    summary = json.loads((ROOT / "analysis/vgp_freeze1_mirror_summary.json").read_text())
    reconciliation = summary["catalog_reconciliation"]
    assert reconciliation == {
        "accession_or_version_drift": 0,
        "catalog_rows": 716,
        "closed_world": True,
        "expected_accession_roots": 581,
        "extra_roots": 0,
        "missing_roots": 0,
        "observed_accession_roots": 581,
        "released_rows": 581,
        "unreleased_rows": 135,
    }
    source_path = ROOT / "analysis/vgp_freeze1_source_inventory.tsv"
    source_objects = source_files = source_bytes = fasta_files = fasta_bytes = 0
    accessions = set()
    with source_path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            source_objects += 1
            accessions.add(row["accession_version"])
            assert row["source_endpoint"] == mirror.TRANSPORT_ENDPOINT + row["source_relative_path"]
            if row["object_type"] == "file":
                source_files += 1
                source_bytes += int(row["size_bytes"])
            if row["sequence_subset"] == "assembly_fasta":
                fasta_files += 1
                fasta_bytes += int(row["size_bytes"])
    assert (source_objects, source_files, source_bytes) == (47870, 43371, 3916877494936)
    assert (fasta_files, fasta_bytes) == (581, 388785172570)
    assert len(accessions) == 581
    manifest_states = set()
    manifest_objects = manifest_bytes = 0
    with (ROOT / "analysis/vgp_freeze1_mirror_manifest.tsv").open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            manifest_objects += 1
            manifest_states.add(row["state"])
            if row["object_type"] == "file":
                manifest_bytes += int(row["size_bytes"])
    assert (manifest_objects, manifest_bytes, manifest_states) == (
        47870,
        3916877494936,
        {"planned"},
    )
    assert summary["fixture"]["passed"] is True
    assert summary["storage"]["arbitrary_global_byte_cap"] is None
    assert summary["bulk_launch"] == {
        "launched": False,
        "reason": "quota_visibility_unavailable_fail_closed",
        "slurm_jobs_launched": 0,
    }
    handoff = (ROOT / "analysis/vgp_freeze1_mirror_handoff.md").read_text()
    assert "unverified historical planning estimates only" in handoff
    assert "no bulk payload transfer or Slurm job was launched" in handoff
