import csv
import gzip
import hashlib
import json
import shutil
import struct
from pathlib import Path

import pytest

from analysis import build_vgp_freeze1_bgzf as bgzf


def _mirror_manifest(path: Path, source: Path, *, subset: str = "assembly_fasta") -> None:
    fields = [
        "canonical_vgp_root",
        "mirror_root",
        "inventory_id",
        "accession_version",
        "source_relative_path",
        "object_type",
        "size_bytes",
        "sequence_subset",
        "observed_bytes",
        "local_sha256",
        "durable_path",
        "state",
    ]
    payload = {
        "canonical_vgp_root": str(source.parents[2]),
        "mirror_root": str(source.parents[1]),
        "inventory_id": "inventory-1",
        "accession_version": "GCA_000000001.1",
        "source_relative_path": "GCA/000/000/001/GCA_000000001.1/GCA_000000001.1.fa.gz",
        "object_type": "file",
        "size_bytes": str(source.stat().st_size),
        "sequence_subset": subset,
        "observed_bytes": str(source.stat().st_size),
        "local_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "durable_path": str(source),
        "state": "verified",
    }
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerow(payload)


def test_inventory_is_closed_world_and_excludes_non_assembly(tmp_path):
    source = tmp_path / "freeze1" / "objects" / "assembly.fa.gz"
    source.parent.mkdir(parents=True)
    source.write_bytes(gzip.compress(b">chr1\nACGT\n"))
    manifest = tmp_path / "mirror.tsv"
    _mirror_manifest(manifest, source)

    rows = bgzf.load_inventory(manifest, tmp_path / "derived")

    assert len(rows) == 1
    assert rows[0].source_format == "gzip"
    assert rows[0].resource_class == "tiny"
    assert rows[0].derived_relative_path.endswith(".fa.gz")


def test_detects_uncompressed_gzip_and_bgzf(tmp_path):
    plain = tmp_path / "a.fa"
    ordinary = tmp_path / "b.fa.gz"
    block = tmp_path / "c.fa.gz"
    plain.write_bytes(b">x\nA\n")
    ordinary.write_bytes(gzip.compress(plain.read_bytes()))
    block.write_bytes(
        bytes.fromhex("1f8b08040000000000ff0600424302001b00")
        + bytes.fromhex("03000000000000000000")
    )

    assert bgzf.detect_format(plain) == "uncompressed_fasta"
    assert bgzf.detect_format(ordinary) == "gzip"
    assert bgzf.detect_format(block) == "bgzf"


@pytest.mark.parametrize(
    ("size", "expected"),
    [
        (1, "tiny"),
        (256 * 1024**2, "tiny"),
        (256 * 1024**2 + 1, "small"),
        (1024**3 + 1, "medium"),
        (4 * 1024**3 + 1, "large"),
    ],
)
def test_resource_class_boundaries(size, expected):
    assert bgzf.resource_class_for_size(size).name == expected


def test_parse_gzi_rejects_truncated_and_non_monotonic(tmp_path):
    truncated = tmp_path / "bad.gzi"
    truncated.write_bytes(struct.pack("<Q", 2) + struct.pack("<QQ", 100, 200))
    with pytest.raises(bgzf.ValidationError, match="GZI_SIZE"):
        bgzf.validate_gzi(truncated, 1000)

    unordered = tmp_path / "unordered.gzi"
    unordered.write_bytes(
        struct.pack("<Q", 2)
        + struct.pack("<QQ", 200, 300)
        + struct.pack("<QQ", 100, 400)
    )
    with pytest.raises(bgzf.ValidationError, match="GZI_ORDER"):
        bgzf.validate_gzi(unordered, 1000)


def test_finalize_accounts_for_missing_worker_as_failed(tmp_path):
    source = tmp_path / "freeze1" / "objects" / "assembly.fa.gz"
    source.parent.mkdir(parents=True)
    source.write_bytes(gzip.compress(b">chr1\nACGT\n"))
    manifest = tmp_path / "mirror.tsv"
    _mirror_manifest(manifest, source)
    derived = tmp_path / "derived"
    inventory = bgzf.write_inventory(manifest, derived, tmp_path / "inventory.tsv")

    summary = bgzf.finalize(inventory, derived, tmp_path / "shared.tsv", tmp_path / "summary.json")

    assert summary["closed_world"] is True
    assert summary["counts"] == {"failed": 1}
    assert b"\r\n" not in (tmp_path / "shared.tsv").read_bytes()
    with (tmp_path / "shared.tsv").open(newline="") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert row["state"] == "failed"
    assert row["reason_code"] == "NO_WORKER_STATUS"


def test_atomic_write_json_replaces_complete_document(tmp_path):
    target = tmp_path / "state.json"
    bgzf.atomic_write_json(target, {"state": "working"})
    bgzf.atomic_write_json(target, {"state": "promoted", "n": 2})
    assert json.loads(target.read_text()) == {"n": 2, "state": "promoted"}
    assert not list(tmp_path.glob("*.tmp"))


@pytest.mark.skipif(
    shutil.which("bgzip") is None or shutil.which("samtools") is None,
    reason="pinned htslib/samtools profile not active",
)
def test_worker_streams_validates_indexes_and_promotes_real_triplet(tmp_path):
    source = tmp_path / "mirror" / "sample.fa.gz"
    source.parent.mkdir()
    uncompressed = b">chr1 description\nACGTACGT\n>chr2\nNNNNAAA\n"
    source.write_bytes(gzip.compress(uncompressed))
    relative = "GCA/000/000/001/GCA_000000001.1/GCA_000000001.1.fa.gz"
    row = bgzf.InventoryRow(
        task_index=0,
        inventory_id="object-1",
        accession_version="GCA_000000001.1",
        source_relative_path=relative,
        source_path=str(source),
        source_format="gzip",
        source_bytes=source.stat().st_size,
        source_compressed_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        mirror_state="verified",
        derived_relative_path=f"objects/{relative}",
        resource_class="tiny",
        cpus=2,
        memory_mb=4096,
        scratch_bytes=1,
        walltime="00:10:00",
    )
    inventory = tmp_path / "inventory.tsv"
    bgzf._write_tsv(inventory, bgzf.INVENTORY_FIELDS, [row.as_dict()])
    derived = tmp_path / "derived"
    status = bgzf.run_worker(inventory, derived, 0, tmp_path / "scratch")

    payload = derived / row.derived_relative_path
    assert status["state"] == "converted"
    assert status["source_decompressed_sha256"] == hashlib.sha256(uncompressed).hexdigest()
    assert status["sequence_count"] == 2
    assert payload.is_file()
    assert Path(str(payload) + ".gzi").is_file()
    assert Path(str(payload) + ".fai").is_file()
    assert gzip.decompress(payload.read_bytes()) == uncompressed
    assert (payload.parent / "provenance.json").is_file()
