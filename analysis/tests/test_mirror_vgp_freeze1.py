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
    assert row["canonical_vgp_root"] == "/moosefs/erikg/vgp"
    assert row["mirror_root"] == "/moosefs/erikg/vgp/freeze1"
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


def test_cas_promotion_and_second_view_reuse_are_atomic(tmp_path, monkeypatch):
    monkeypatch.setattr(mirror, "VGP_DATA_ROOT", tmp_path / "vgp")
    content = b"canonical shared bytes"
    part = tmp_path / "first.part"
    part.write_bytes(content)
    first_view = tmp_path / "vgp/freeze1/objects/source/first"
    digest, cas_path = mirror.promote_to_cas(part, first_view, expected_size=len(content))
    assert digest == hashlib.sha256(content).hexdigest()
    assert cas_path == tmp_path / "vgp/objects/sha256" / digest[:2] / digest[2:4] / digest
    assert cas_path.read_bytes() == content
    assert first_view.read_bytes() == content
    assert cas_path.stat().st_ino == first_view.stat().st_ino
    second_view = tmp_path / "vgp/freeze1/objects/source/second"
    mirror.atomic_publish_view(cas_path, second_view)
    assert second_view.stat().st_ino == cas_path.stat().st_ino


def test_missing_quota_helper_is_observability_only(tmp_path, monkeypatch):
    monkeypatch.setattr(mirror.shutil, "which", lambda _command: None)
    obj = mirror.InventoryObject(
        "id", 2, "GCA_000000000.1", "GCA/000/000/000/GCA_000000000.1/file",
        "file", 7, "2024/01/02-03:04:05", "", "non_sequence_product_or_metadata"
    )
    evidence = mirror.storage_evidence(tmp_path, [obj], concurrency=2)
    assert evidence["filesystem_capacity_adequate"] is True
    assert evidence["quota_visibility_adequate"] is False
    assert evidence["quota_visibility_is_policy_gate"] is False
    assert evidence["adequate"] is True
    assert evidence["gate_reason"].endswith("quota_helper_unavailable")


def test_process_restart_recovers_uncommitted_partial_accounting(tmp_path, monkeypatch):
    monkeypatch.setattr(mirror, "VGP_DATA_ROOT", tmp_path / "vgp")
    root = tmp_path / "vgp/freeze1"
    database = root / "state/mirror.sqlite3"
    obj = mirror.InventoryObject(
        "id", 2, "GCA_000000000.1", "GCA/000/000/000/GCA_000000000.1/payload",
        "file", 7, "2024/01/02-03:04:05", "", "non_sequence_product_or_metadata"
    )
    mirror.init_database(database, [obj], root)
    row = mirror.database_rows(database)[0]
    part = Path(row["staging_path"])
    part.parent.mkdir(parents=True)
    part.write_bytes(b"part")

    def resume(_source, resumed_part, timeout=300):
        assert timeout == 300
        resumed_part.write_bytes(b"partial")
        return 3

    monkeypatch.setattr(mirror, "rsync_transfer", resume)
    mirror.process_file(database, root, row)
    completed = mirror.database_rows(database)[0]
    assert completed["state"] == "verified"
    assert completed["observed_bytes"] == 7
    assert completed["transferred_bytes"] == 7
    assert completed["attempts"] == 2


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


def test_verified_upstream_conflict_requires_reproduction_and_authoritative_alternate():
    observed_md5 = hashlib.md5(b"official current bytes").hexdigest()
    observed_sha256 = hashlib.sha256(b"official current bytes").hexdigest()
    evidence = {
        "canonical_vgp_root": "/moosefs/erikg/vgp",
        "inventory_id": "inventory-id",
        "source_relative_path": "GCA/000/000/000/GCA_000000000.1/report.txt",
        "sequence_subset": "non_sequence_product_or_metadata",
        "frozen_catalog": {"algorithm": "md5", "digest": "0" * 32, "size_bytes": 22},
        "official_source_attempts": [
            {
                "source_url": mirror.TRANSPORT_ENDPOINT + "GCA/000/000/000/GCA_000000000.1/report.txt",
                "retrieval_started_utc": "2026-01-01T00:00:00Z",
                "retrieval_completed_utc": "2026-01-01T00:00:01Z",
                "size_bytes": 22,
                "md5": observed_md5,
                "sha256": observed_sha256,
                "quarantine_path": "/quarantine/attempt-1",
            },
            {
                "source_url": mirror.TRANSPORT_ENDPOINT + "GCA/000/000/000/GCA_000000000.1/report.txt",
                "retrieval_started_utc": "2026-01-01T00:01:00Z",
                "retrieval_completed_utc": "2026-01-01T00:01:01Z",
                "size_bytes": 22,
                "md5": observed_md5,
                "sha256": observed_sha256,
                "quarantine_path": "/quarantine/attempt-2",
            },
        ],
        "authoritative_alternate": {
            "assembly_accession": "GCA_000000000.1",
            "source_url": "https://ftp.ncbi.nlm.nih.gov/report.txt",
            "checksum_catalog_url": "https://ftp.ncbi.nlm.nih.gov/md5checksums.txt",
            "retrieval_started_utc": "2026-01-01T00:02:00Z",
            "retrieval_completed_utc": "2026-01-01T00:02:01Z",
            "size_bytes": 22,
            "md5": observed_md5,
            "sha256": observed_sha256,
            "catalog_algorithm": "md5",
            "catalog_digest": observed_md5,
            "quarantine_path": "/quarantine/ncbi",
        },
        "resolution": "VERIFIED_UPSTREAM_CONFLICT",
    }
    mirror.validate_verified_upstream_conflict(evidence)
    divergent_ncbi = json.loads(json.dumps(evidence))
    divergent_ncbi["authoritative_alternate"].update(
        {"size_bytes": 21, "md5": "1" * 32, "sha256": "2" * 64, "catalog_digest": "1" * 32}
    )
    mirror.validate_verified_upstream_conflict(divergent_ncbi)
    absent_ncbi = json.loads(json.dumps(evidence))
    absent_ncbi["authoritative_alternate"] = {
        "assembly_accession": "GCA_000000000.1",
        "source_url": "https://ftp.ncbi.nlm.nih.gov/accession/",
        "checksum_catalog_url": "https://ftp.ncbi.nlm.nih.gov/accession/md5checksums.txt",
        "retrieval_started_utc": "2026-01-01T00:02:00Z",
        "retrieval_completed_utc": "2026-01-01T00:02:01Z",
        "checksum_catalog_retrieval": {"sha256": "3" * 64},
        "object_available": False,
        "catalog_entry_found": False,
        "searched_equivalent_names": ["hub.txt"],
    }
    mirror.validate_verified_upstream_conflict(absent_ncbi)
    evidence["official_source_attempts"] = evidence["official_source_attempts"][:1]
    with pytest.raises(mirror.MirrorError, match="two independent"):
        mirror.validate_verified_upstream_conflict(evidence)


def test_sequence_conflict_cannot_receive_terminal_metadata_exception():
    with pytest.raises(mirror.MirrorError, match="non-sequence metadata"):
        mirror.validate_verified_upstream_conflict(
            {
                "canonical_vgp_root": "/moosefs/erikg/vgp",
                "sequence_subset": "assembly_fasta",
            }
        )


def test_process_file_reproduces_checksum_conflict_twice_then_continues(tmp_path, monkeypatch):
    monkeypatch.setattr(mirror, "VGP_DATA_ROOT", tmp_path / "vgp")
    root = tmp_path / "vgp/freeze1"
    database = root / "state/mirror.sqlite3"
    obj = mirror.InventoryObject(
        "id", 2, "GCA_000000000.1", "GCA/000/000/000/GCA_000000000.1/report.txt",
        "file", 3, "2024/01/02-03:04:05", "", "non_sequence_product_or_metadata",
        "md5", "0" * 32,
    )
    mirror.init_database(database, [obj], root)
    calls = []

    def conflicting_source(_source, part, timeout=300):
        calls.append(timeout)
        part.parent.mkdir(parents=True, exist_ok=True)
        part.write_bytes(b"bad")
        return 3

    monkeypatch.setattr(mirror, "rsync_transfer", conflicting_source)
    monkeypatch.setattr(mirror.time, "sleep", lambda _seconds: None)
    mirror.process_file(database, root, mirror.database_rows(database)[0])
    row = mirror.database_rows(database)[0]
    assert calls == [300, 300]
    assert row["state"] == "quarantined"
    assert row["attempts"] == 2
    assert row["transferred_bytes"] == 6
    connection = sqlite3.connect(database)
    assert connection.execute("SELECT COUNT(*) FROM quarantine_events").fetchone()[0] == 2
    connection.close()


def test_terminal_metadata_exception_is_accounted_but_not_remaining(tmp_path, monkeypatch):
    monkeypatch.setattr(mirror, "VGP_DATA_ROOT", tmp_path / "vgp")
    root = tmp_path / "vgp/freeze1"
    database = root / "state/mirror.sqlite3"
    obj = mirror.InventoryObject(
        "id", 2, "GCA_000000000.1", "GCA/000/000/000/GCA_000000000.1/report.txt",
        "file", 3, "2024/01/02-03:04:05", "", "non_sequence_product_or_metadata",
    )
    mirror.init_database(database, [obj], root)
    mirror.update_database(
        database, "id", state="verified_upstream_conflict", observed_bytes=3
    )
    progress = mirror.mirror_progress(database, root)["progress"]
    assert progress["verified_upstream_conflict_files"] == 1
    assert progress["verified_upstream_conflict_logical_bytes"] == 3
    assert progress["remaining_files"] == 0
    assert progress["remaining_bytes"] == 0


def test_terminal_worker_skips_unneeded_full_cas_reindex(tmp_path, monkeypatch):
    monkeypatch.setenv("GUIX_ENVIRONMENT", "/gnu/store/test")
    monkeypatch.setattr(mirror, "VGP_DATA_ROOT", tmp_path / "vgp")
    monkeypatch.setattr(mirror, "EXCEPTION_LEDGER_OUTPUT", tmp_path / "exception-ledger.json")
    root = tmp_path / "vgp/freeze1"
    database = root / "state/mirror.sqlite3"
    obj = mirror.InventoryObject(
        "id", 2, "GCA_000000000.1", "GCA/000/000/000/GCA_000000000.1/report.txt",
        "file", 3, "2024/01/02-03:04:05", "", "non_sequence_product_or_metadata",
    )
    directory = mirror.InventoryObject(
        "directory-id", 2, "GCA_000000000.1", "GCA/000/000/000/GCA_000000000.1/",
        "directory", 0, "2024/01/02-03:04:05", "", "non_sequence_product_or_metadata",
    )
    mirror.init_database(database, [directory, obj], root)
    mirror.update_database(
        database, "id", state="verified_upstream_conflict", observed_bytes=3
    )
    mirror.update_database(database, "directory-id", state="reused")
    fixture = root / "fixture/fixture-report.json"
    fixture.parent.mkdir(parents=True)
    fixture.write_text(json.dumps({"passed": True}), encoding="utf-8")
    capacity = root / "inventory/capacity-evidence.json"
    capacity.parent.mkdir(parents=True)
    capacity.write_text(
        json.dumps(
            {
                "filesystem_capacity_adequate": True,
                "write_probe_passed": True,
                "gate_reason": "verified",
                "requirements": {"total_bytes": 0, "total_inodes": 0},
            }
        ),
        encoding="utf-8",
    )

    def unexpected_reindex(_cas_root):
        pytest.fail("terminal worker must not rehash the full canonical CAS")

    def unexpected_checksum_rebind(_objects_root, _objects, _ignored=None):
        pytest.fail("terminal worker must not rebind checksums for completed payloads")

    monkeypatch.setattr(mirror, "index_verified_cas", unexpected_reindex)
    monkeypatch.setattr(mirror, "parse_upstream_md5s", unexpected_checksum_rebind)
    original_update = mirror.update_database

    def no_terminal_directory_rewrite(database_path, inventory_id, **updates):
        if inventory_id == "directory-id":
            pytest.fail("terminal worker must not rewrite completed directory rows")
        return original_update(database_path, inventory_id, **updates)

    monkeypatch.setattr(mirror, "update_database", no_terminal_directory_rewrite)
    mirror.run_worker(database, root, 2, capacity)

    connection = sqlite3.connect(database)
    run = connection.execute(
        "SELECT outcome, detail FROM worker_runs ORDER BY run_id DESC LIMIT 1"
    ).fetchone()
    connection.close()
    assert run == ("complete", "frozen_inventory_verified_with_1_upstream_conflict_exception(s)")


def test_ncbi_assembly_urls_require_one_exact_accession_version():
    index = """<a href="GCA_005190385.2_old/">old</a>
<a href="GCA_005190385.3_NGI_Narwhal_2/">exact</a>
<a href="GCA_005190385.4_new/">new</a>"""
    source, checksums = mirror.ncbi_assembly_urls("GCA_005190385.3", index)
    base = (
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/005/190/385/"
        "GCA_005190385.3_NGI_Narwhal_2"
    )
    assert source == base + "/GCA_005190385.3_NGI_Narwhal_2_assembly_report.txt"
    assert checksums == base + "/md5checksums.txt"
    with pytest.raises(mirror.MirrorError, match="0 exact-version"):
        mirror.ncbi_assembly_urls("GCA_005190385.5", index)


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


def test_upstream_md5_normalizes_only_official_mirrordata_prefix(tmp_path):
    accession = "GCA_000000000.1"
    root = mirror.accession_path(accession)
    objects_root = tmp_path / "objects"
    manifest = objects_root / root / "md5sum.txt"
    manifest.parent.mkdir(parents=True)
    target = root + "/payload.fa.gz"
    manifest.write_text("1" * 32 + "  /mirrordata/hubs/" + target + "\n")
    objects = [
        mirror.InventoryObject("manifest", 2, accession, root + "/md5sum.txt", "file", 1,
                               "2024/01/02-03:04:05", "", "non_sequence_product_or_metadata"),
        mirror.InventoryObject("payload", 2, accession, target, "file", 7,
                               "2024/01/02-03:04:05", "", "assembly_fasta"),
    ]
    assert mirror.parse_upstream_md5s(objects_root, objects) == {
        target: ("md5", "1" * 32)
    }
    manifest.write_text("1" * 32 + "  /untrusted/hubs/" + target + "\n")
    with pytest.raises(mirror.MirrorError, match="unrecognized absolute prefix"):
        mirror.parse_upstream_md5s(objects_root, objects)


def test_stale_provider_absolute_checksum_entry_is_audited_not_inventoried(tmp_path):
    accession = "GCA_000000000.1"
    root = mirror.accession_path(accession)
    objects_root = tmp_path / "objects"
    manifest = objects_root / root / "md5sum.txt"
    manifest.parent.mkdir(parents=True)
    orphan = root + "/provider-deleted.txt"
    manifest.write_text("2" * 32 + "  /mirrordata/hubs/" + orphan + "\n")
    objects = [
        mirror.InventoryObject("manifest", 2, accession, root + "/md5sum.txt", "file", 1,
                               "2024/01/02-03:04:05", "", "non_sequence_product_or_metadata")
    ]
    ignored = []
    assert mirror.parse_upstream_md5s(objects_root, objects, ignored) == {}
    assert ignored == [orphan]


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
    assert (manifest_objects, manifest_bytes) == (47870, 3916877494936)
    assert manifest_states <= mirror.VALID_STATES
    assert manifest_states == {"verified", "reused", "verified_upstream_conflict"}
    assert summary["live_progress"]["remaining_files"] == 0
    assert summary["live_progress"]["remaining_bytes"] == 0
    assert summary["live_progress"]["currently_quarantined_files"] == 0
    assert summary["live_progress"]["verified_files"] == 42920
    assert summary["live_progress"]["verified_upstream_conflict_files"] == 451
    assert len(summary["checksum_exceptions"]) == 451
    assert all(
        entry["sequence_subset"] == "non_sequence_product_or_metadata"
        and entry["scientific_sequence_review_required"] is False
        and entry["promotion_permitted"] is False
        for entry in summary["checksum_exceptions"]
    )
    assert summary["state_accounting"]["verified"]["bytes"] > 0
    assert summary["fixture"]["passed"] is True
    assert summary["canonical_vgp_root"] == "/moosefs/erikg/vgp"
    assert summary["mirror_root"] == "/moosefs/erikg/vgp/freeze1"
    assert summary["storage"]["arbitrary_global_byte_cap"] is None
    assert summary["storage"]["quota_visibility_is_policy_gate"] is False
    assert summary["storage"]["adequate"] is True
    assert summary["bulk_launch"]["reason"].startswith("capacity_write_and_inode_headroom_verified")
    assert summary["bulk_launch"]["slurm_jobs_launched"] == 0
    handoff = (ROOT / "analysis/vgp_freeze1_mirror_handoff.md").read_text()
    assert "unverified historical planning estimates only" in handoff
    assert "quota helpers are observability only" in handoff
