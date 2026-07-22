#!/usr/bin/env python3
"""Build and validate the immutable VGP Freeze 1 BGZF assembly view.

The source mirror is read-only.  Each worker streams an assembly through a
node-local directory, verifies source and derivative identity, promotes a
content-addressed set on the shared filesystem, and finally renames one
accession directory.  The directory rename is the triplet commit point: a
consumer can never observe only some of FASTA.gz, .gzi, and .fai.
"""

from __future__ import annotations

import argparse
import csv
import dataclasses
import datetime as dt
import errno
import fcntl
import gzip
import hashlib
import json
import os
import resource
import shutil
import struct
import subprocess
import sys
import tempfile
import time
import uuid
from collections import Counter
from pathlib import Path
from typing import Iterable, Iterator, Mapping, Sequence


DEFAULT_MIRROR_MANIFEST = Path("analysis/vgp_freeze1_mirror_manifest.tsv")
DEFAULT_DERIVED_ROOT = Path("/moosefs/erikg/vgp/derived/freeze1-bgzf")
DEFAULT_INVENTORY = Path("analysis/vgp_freeze1_bgzf_inventory.tsv")
DEFAULT_SHARED_MANIFEST = Path("analysis/vgp_freeze1_bgzf_manifest.tsv")
DEFAULT_SUMMARY = Path("analysis/vgp_freeze1_bgzf_summary.json")
CONFIG_VARIABLE = "VGP_FREEZE1_ASSEMBLY_BGZF_ROOT"
BGZF_EOF = bytes.fromhex(
    "1f8b08040000000000ff0600424302001b0003000000000000000000"
)
PROMOTED_STATES = {"converted", "reused"}
MIRROR_TERMINAL_STATES = {"verified", "reused", "verified_upstream_conflict"}


class PipelineError(RuntimeError):
    reason_code = "PIPELINE_ERROR"


class InventoryError(PipelineError):
    reason_code = "INVENTORY_INVALID"


class ValidationError(PipelineError):
    reason_code = "VALIDATION_FAILED"


class ScratchError(PipelineError):
    reason_code = "SCRATCH_INSUFFICIENT"


class ToolError(PipelineError):
    reason_code = "TOOL_FAILED"


@dataclasses.dataclass(frozen=True)
class ResourceClass:
    name: str
    max_source_bytes: int | None
    cpus: int
    memory_mb: int
    scratch_multiplier: float
    scratch_floor_bytes: int
    walltime: str
    max_parallel: int

    def scratch_bytes(self, source_bytes: int) -> int:
        return max(self.scratch_floor_bytes, int(source_bytes * self.scratch_multiplier))


RESOURCE_CLASSES = (
    ResourceClass("tiny", 256 * 1024**2, 2, 4096, 2.0, 2 * 1024**3, "02:00:00", 6),
    ResourceClass("small", 1024**3, 4, 8192, 1.8, 4 * 1024**3, "04:00:00", 8),
    ResourceClass("medium", 4 * 1024**3, 6, 16384, 1.7, 8 * 1024**3, "08:00:00", 4),
    ResourceClass("large", None, 8, 32768, 1.6, 24 * 1024**3, "18:00:00", 2),
)


INVENTORY_FIELDS = (
    "task_index",
    "inventory_id",
    "accession_version",
    "source_relative_path",
    "source_path",
    "source_format",
    "source_bytes",
    "source_compressed_sha256",
    "mirror_state",
    "derived_relative_path",
    "resource_class",
    "cpus",
    "memory_mb",
    "scratch_bytes",
    "walltime",
)

SHARED_FIELDS = INVENTORY_FIELDS + (
    "state",
    "reason_code",
    "attempts",
    "derived_path",
    "gzi_path",
    "fai_path",
    "provenance_path",
    "source_decompressed_sha256",
    "source_decompressed_bytes",
    "sequence_count",
    "names_lengths_sha256",
    "derived_bgzf_sha256",
    "derived_bgzf_bytes",
    "gzi_sha256",
    "gzi_bytes",
    "fai_sha256",
    "fai_bytes",
    "compression_ratio",
    "elapsed_seconds",
    "cpu_seconds",
    "peak_memory_bytes",
    "peak_scratch_bytes",
    "random_access_probe_count",
    "bgzf_block_count",
    "updated_at_utc",
)


@dataclasses.dataclass(frozen=True)
class InventoryRow:
    task_index: int
    inventory_id: str
    accession_version: str
    source_relative_path: str
    source_path: str
    source_format: str
    source_bytes: int
    source_compressed_sha256: str
    mirror_state: str
    derived_relative_path: str
    resource_class: str
    cpus: int
    memory_mb: int
    scratch_bytes: int
    walltime: str

    @classmethod
    def from_mapping(cls, row: Mapping[str, str]) -> "InventoryRow":
        return cls(
            task_index=int(row["task_index"]),
            inventory_id=row["inventory_id"],
            accession_version=row["accession_version"],
            source_relative_path=row["source_relative_path"],
            source_path=row["source_path"],
            source_format=row["source_format"],
            source_bytes=int(row["source_bytes"]),
            source_compressed_sha256=row["source_compressed_sha256"],
            mirror_state=row["mirror_state"],
            derived_relative_path=row["derived_relative_path"],
            resource_class=row["resource_class"],
            cpus=int(row["cpus"]),
            memory_mb=int(row["memory_mb"]),
            scratch_bytes=int(row["scratch_bytes"]),
            walltime=row["walltime"],
        )

    def as_dict(self) -> dict[str, str | int]:
        return dataclasses.asdict(self)


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat().replace("+00:00", "Z")


def sha256_file(path: Path, chunk_size: int = 8 * 1024**2) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fsync_path(path: Path) -> None:
    flags = os.O_RDONLY | getattr(os, "O_DIRECTORY", 0)
    fd = os.open(path, flags)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def atomic_write_json(path: Path, payload: Mapping[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{uuid.uuid4().hex}.tmp")
    with temporary.open("w") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)
    fsync_path(path.parent)


def detect_format(path: Path) -> str:
    with path.open("rb") as handle:
        header = handle.read(18)
    if not header.startswith(b"\x1f\x8b"):
        return "uncompressed_fasta"
    if len(header) >= 18 and header[2] == 8 and header[3] & 4:
        xlen = struct.unpack("<H", header[10:12])[0]
        if xlen >= 6 and header[12:16] == b"BC\x02\x00":
            return "bgzf"
    return "gzip"


def resource_class_for_size(size: int) -> ResourceClass:
    for item in RESOURCE_CLASSES:
        if item.max_source_bytes is None or size <= item.max_source_bytes:
            return item
    raise AssertionError("unreachable resource class")


def _safe_relative_path(value: str) -> Path:
    candidate = Path(value)
    if candidate.is_absolute() or ".." in candidate.parts or not candidate.parts:
        raise InventoryError(f"unsafe source_relative_path: {value!r}")
    return candidate


def load_inventory(mirror_manifest: Path, derived_root: Path) -> list[InventoryRow]:
    del derived_root  # Root is intentionally not baked into the relocatable inventory.
    selected: list[dict[str, str]] = []
    seen_paths: set[str] = set()
    with mirror_manifest.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "inventory_id", "accession_version", "source_relative_path", "object_type",
            "sequence_subset", "observed_bytes", "local_sha256", "durable_path", "state",
        }
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise InventoryError("terminal mirror manifest lacks required fields")
        for raw in reader:
            if raw["object_type"] != "file" or raw["sequence_subset"] != "assembly_fasta":
                continue
            relative = _safe_relative_path(raw["source_relative_path"])
            if raw["source_relative_path"] in seen_paths:
                raise InventoryError(f"duplicate assembly FASTA: {relative}")
            seen_paths.add(raw["source_relative_path"])
            if raw["state"] not in MIRROR_TERMINAL_STATES:
                raise InventoryError(f"nonterminal assembly FASTA {relative}: {raw['state']}")
            source = Path(raw["durable_path"])
            if not source.is_file():
                raise InventoryError(f"missing mirrored assembly FASTA: {source}")
            observed = int(raw["observed_bytes"])
            if source.stat().st_size != observed:
                raise InventoryError(f"source size drift for {source}")
            selected.append(raw)
    if not selected:
        raise InventoryError("zero assembly_fasta rows in terminal mirror manifest")
    selected.sort(key=lambda row: row["source_relative_path"])
    output: list[InventoryRow] = []
    for index, raw in enumerate(selected):
        size = int(raw["observed_bytes"])
        resources = resource_class_for_size(size)
        source = Path(raw["durable_path"])
        output.append(
            InventoryRow(
                task_index=index,
                inventory_id=raw["inventory_id"],
                accession_version=raw["accession_version"],
                source_relative_path=raw["source_relative_path"],
                source_path=str(source),
                source_format=detect_format(source),
                source_bytes=size,
                source_compressed_sha256=raw["local_sha256"],
                mirror_state=raw["state"],
                derived_relative_path=str(Path("objects") / raw["source_relative_path"]),
                resource_class=resources.name,
                cpus=resources.cpus,
                memory_mb=resources.memory_mb,
                scratch_bytes=resources.scratch_bytes(size),
                walltime=resources.walltime,
            )
        )
    return output


def _write_tsv(path: Path, fieldnames: Sequence[str], rows: Iterable[Mapping[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=fieldnames,
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)
    fsync_path(path.parent)


def write_inventory(mirror_manifest: Path, derived_root: Path, output: Path) -> Path:
    rows = load_inventory(mirror_manifest, derived_root)
    _write_tsv(output, INVENTORY_FIELDS, (row.as_dict() for row in rows))
    derived_root.mkdir(parents=True, exist_ok=True)
    for directory in ("objects", "cas/decompressed-sha256", "locks", "status", "staging", "logs"):
        (derived_root / directory).mkdir(parents=True, exist_ok=True)
    config = {
        "schema_version": 1,
        "config_variable": CONFIG_VARIABLE,
        "config_value": str((derived_root / "objects").resolve()),
        "derived_root": str(derived_root.resolve()),
        "mirror_manifest": str(mirror_manifest.resolve()),
        "mirror_manifest_sha256": sha256_file(mirror_manifest),
        "inventory_path": str(output.resolve()),
        "inventory_sha256": sha256_file(output),
        "assembly_count": len(rows),
        "source_bytes": sum(row.source_bytes for row in rows),
        "created_at_utc": utc_now(),
    }
    atomic_write_json(derived_root / "CONFIG.json", config)
    env_path = derived_root / "CONFIG.env"
    env_path.write_text(f"export {CONFIG_VARIABLE}={derived_root.resolve() / 'objects'}\n")
    return output


def read_inventory(path: Path) -> list[InventoryRow]:
    with path.open(newline="") as handle:
        rows = [InventoryRow.from_mapping(row) for row in csv.DictReader(handle, delimiter="\t")]
    for expected, row in enumerate(rows):
        if row.task_index != expected:
            raise InventoryError(f"task index discontinuity at {expected}: {row.task_index}")
    return rows


def _run(command: Sequence[str], *, stdin=None, stdout=None) -> subprocess.CompletedProcess:
    if stdout is None:
        stdout = subprocess.PIPE
    try:
        return subprocess.run(command, stdin=stdin, stdout=stdout, stderr=subprocess.PIPE, check=True)
    except FileNotFoundError as error:
        raise ToolError(f"tool absent: {command[0]}") from error
    except subprocess.CalledProcessError as error:
        stderr = error.stderr.decode(errors="replace") if isinstance(error.stderr, bytes) else error.stderr
        raise ToolError(f"command failed ({' '.join(command)}): {stderr[-2000:]}") from error


def _capture_source_and_compress(row: InventoryRow, work: Path, bgzip_threads: int) -> dict[str, object]:
    source_sha = work / "source.sha256"
    source_bytes = work / "source.bytes"
    source_sequences = work / "source.sequences.tsv"
    output = work / "payload.fa.gz"
    if row.source_format == "gzip":
        source_handle = gzip.open(row.source_path, "rb")
    elif row.source_format == "uncompressed_fasta":
        source_handle = Path(row.source_path).open("rb")
    else:
        raise ValidationError(f"unsupported source format for compression: {row.source_format}")

    digest = hashlib.sha256()
    byte_value = 0
    sequence_rows: list[tuple[str, int]] = []
    sequence_name: str | None = None
    sequence_length = 0
    with source_handle as source, output.open("wb") as sink:
        try:
            process = subprocess.Popen(
                ["bgzip", "-@", str(max(1, bgzip_threads)), "-c"],
                stdin=subprocess.PIPE,
                stdout=sink,
                stderr=subprocess.PIPE,
            )
        except FileNotFoundError as error:
            raise ToolError("tool absent: bgzip") from error
        assert process.stdin is not None and process.stderr is not None
        try:
            for line in source:
                digest.update(line)
                byte_value += len(line)
                process.stdin.write(line)
                if line.startswith(b">"):
                    if sequence_name is not None:
                        sequence_rows.append((sequence_name, sequence_length))
                    header = line[1:].rstrip(b"\r\n").split(None, 1)
                    if not header:
                        raise ValidationError("empty FASTA sequence name")
                    sequence_name = header[0].decode("utf-8")
                    sequence_length = 0
                else:
                    sequence_length += len(b"".join(line.split()))
        except Exception as error:
            try:
                process.stdin.close()
            except BrokenPipeError:
                pass
            if process.poll() is None:
                process.terminate()
            stderr = process.stderr.read()
            process.wait()
            if isinstance(error, PipelineError):
                raise
            raise ToolError(
                f"BGZF compression stream failed: {stderr.decode(errors='replace')[-2000:]}"
            ) from error
        process.stdin.close()
        stderr = process.stderr.read()
        returncode = process.wait()
        if returncode:
            raise ToolError(
                f"BGZF compression failed ({returncode}): {stderr.decode(errors='replace')[-2000:]}"
            )
    if sequence_name is not None:
        sequence_rows.append((sequence_name, sequence_length))
    if not sequence_rows:
        raise ValidationError("source FASTA has no sequences")

    sha_value = digest.hexdigest()
    source_sha.write_text(f"{sha_value}  -\n")
    source_bytes.write_text(f"{byte_value}\n")
    with source_sequences.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerows(sequence_rows)
    if not output.is_file() or output.stat().st_size == 0:
        raise ValidationError("BGZF_ZERO_OUTPUT: compressor produced no payload")
    return {
        "payload": output,
        "source_decompressed_sha256": sha_value,
        "source_decompressed_bytes": byte_value,
        "source_sequences": source_sequences,
    }


def _capture_existing_bgzf(row: InventoryRow, work: Path) -> dict[str, object]:
    output = work / "payload.fa.gz"
    try:
        os.link(row.source_path, output)
    except OSError as error:
        if error.errno != errno.EXDEV:
            raise
        shutil.copyfile(row.source_path, output)
    source_sha = hashlib.sha256()
    source_size = 0
    process = subprocess.Popen(["bgzip", "-dc", row.source_path], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    assert process.stdout is not None
    with (work / "source.stream").open("wb") as sink:
        for chunk in iter(lambda: process.stdout.read(8 * 1024**2), b""):
            source_sha.update(chunk)
            source_size += len(chunk)
            sink.write(chunk)
    stderr = process.communicate()[1]
    if process.returncode:
        raise ToolError(f"BGZF source decompression failed: {stderr.decode(errors='replace')}")
    # This branch is rare; bounded node-local materialization is permitted only
    # for pre-existing BGZF reuse validation, never on shared storage.
    sequence_path = work / "source.sequences.tsv"
    awk = r'''BEGIN{OFS="\t"} /^>/{if(seen)print name,seq_len; h=substr($0,2);sub(/\r$/,"",h);split(h,a,/[ \t]/);name=a[1];seq_len=0;seen=1;next}{gsub(/[ \t\r]/,"",$0);seq_len+=length($0)} END{if(seen)print name,seq_len}'''
    with (work / "source.stream").open("rb") as stream, sequence_path.open("wb") as sink:
        _run(["awk", awk], stdin=stream, stdout=sink)
    (work / "source.stream").unlink()
    return {
        "payload": output,
        "source_decompressed_sha256": source_sha.hexdigest(),
        "source_decompressed_bytes": source_size,
        "source_sequences": sequence_path,
    }


def validate_gzi(path: Path, bgzf_bytes: int) -> int:
    raw_size = path.stat().st_size
    with path.open("rb") as handle:
        header = handle.read(8)
        if len(header) != 8:
            raise ValidationError("GZI_SIZE: missing entry count")
        count = struct.unpack("<Q", header)[0]
        if raw_size != 8 + 16 * count:
            raise ValidationError(f"GZI_SIZE: count={count}, bytes={raw_size}")
        prior_compressed = prior_uncompressed = -1
        for _ in range(count):
            compressed, uncompressed = struct.unpack("<QQ", handle.read(16))
            if compressed <= prior_compressed or uncompressed <= prior_uncompressed:
                raise ValidationError("GZI_ORDER: offsets are not strictly increasing")
            if compressed >= bgzf_bytes:
                raise ValidationError("GZI_RANGE: compressed offset outside BGZF")
            prior_compressed, prior_uncompressed = compressed, uncompressed
    return count


def validate_bgzf_structure(path: Path) -> int:
    size = path.stat().st_size
    if size < len(BGZF_EOF):
        raise ValidationError("BGZF_SIZE: payload shorter than EOF marker")
    with path.open("rb") as handle:
        handle.seek(-len(BGZF_EOF), os.SEEK_END)
        if handle.read() != BGZF_EOF:
            raise ValidationError("BGZF_EOF: canonical EOF marker absent")
        offset = 0
        blocks = 0
        while offset < size:
            handle.seek(offset)
            header = handle.read(18)
            if len(header) != 18 or header[:4] != b"\x1f\x8b\x08\x04":
                raise ValidationError(f"BGZF_HEADER: malformed block at {offset}")
            xlen = struct.unpack("<H", header[10:12])[0]
            if xlen < 6 or header[12:16] != b"BC\x02\x00":
                raise ValidationError(f"BGZF_EXTRA: BC field absent at {offset}")
            block_size = struct.unpack("<H", header[16:18])[0] + 1
            if block_size < 28 or block_size > 65536 or offset + block_size > size:
                raise ValidationError(f"BGZF_BLOCK_SIZE: invalid block at {offset}")
            offset += block_size
            blocks += 1
        if offset != size:
            raise ValidationError("BGZF_BOUNDARY: terminal block exceeds file")
    return blocks


def _compare_source_sequences_to_fai(source_sequences: Path, fai: Path) -> tuple[int, str, list[tuple[str, int]]]:
    digest = hashlib.sha256()
    probes: list[tuple[str, int]] = []
    count = 0
    with source_sequences.open("rb") as source, fai.open("rb") as index:
        while True:
            left = source.readline()
            right = index.readline()
            if not left and not right:
                break
            if not left or not right:
                raise ValidationError("FAI_SEQUENCE_COUNT: source and FAI differ")
            fields = right.rstrip(b"\n").split(b"\t")
            if len(fields) < 5:
                raise ValidationError("FAI_FORMAT: fewer than five columns")
            canonical = fields[0] + b"\t" + fields[1] + b"\n"
            if left != canonical:
                raise ValidationError(
                    f"FAI_NAMES_LENGTHS: source={left[:160]!r}, fai={canonical[:160]!r}"
                )
            digest.update(left)
            count += 1
            if len(probes) < 3:
                probes.append((fields[0].decode(), int(fields[1])))
    if count == 0:
        raise ValidationError("FASTA_EMPTY: no sequences")
    # First, median-ish among the first records, and final are sufficient to
    # exercise distinct virtual offsets; exact stream identity covers content.
    all_tail: tuple[str, int] | None = None
    with fai.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            all_tail = (fields[0], int(fields[1]))
    if all_tail and all_tail not in probes:
        probes.append(all_tail)
    return count, digest.hexdigest(), probes


def _verify_decompressed(path: Path) -> tuple[str, int]:
    digest = hashlib.sha256()
    total = 0
    process = subprocess.Popen(["bgzip", "-dc", str(path)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    assert process.stdout is not None
    for chunk in iter(lambda: process.stdout.read(8 * 1024**2), b""):
        digest.update(chunk)
        total += len(chunk)
    stderr = process.communicate()[1]
    if process.returncode:
        raise ToolError(f"bgzip decompression verification failed: {stderr.decode(errors='replace')[-2000:]}")
    return digest.hexdigest(), total


def _random_access_probes(path: Path, probes: Sequence[tuple[str, int]]) -> int:
    checked = 0
    chosen = probes[:2] + probes[-1:]
    seen: set[str] = set()
    for name, length in chosen:
        if name in seen or length <= 0:
            continue
        seen.add(name)
        midpoint = max(1, (length + 1) // 2)
        start = max(1, midpoint - 50)
        end = min(length, start + 100)
        result = _run(["samtools", "faidx", str(path), f"{name}:{start}-{end}"])
        lines = result.stdout.decode().splitlines() if result.stdout else []
        bases = "".join(lines[1:])
        if len(bases) != end - start + 1:
            raise ValidationError(f"FAIDX_PROBE: {name}:{start}-{end} returned {len(bases)} bases")
        checked += 1
    if checked == 0:
        raise ValidationError("FAIDX_PROBE: no random access probe completed")
    return checked


def validate_triplet(
    payload: Path,
    source_sequences: Path,
    source_sha256: str,
    source_bytes: int,
) -> dict[str, object]:
    _run(["bgzip", "-t", str(payload)])
    blocks = validate_bgzf_structure(payload)
    gzi = Path(str(payload) + ".gzi")
    fai = Path(str(payload) + ".fai")
    gzi_entries = validate_gzi(gzi, payload.stat().st_size)
    sequence_count, names_digest, probes = _compare_source_sequences_to_fai(source_sequences, fai)
    derived_sha, derived_bytes = _verify_decompressed(payload)
    if derived_sha != source_sha256:
        raise ValidationError("DECOMPRESSED_SHA256: source and BGZF differ")
    if derived_bytes != source_bytes:
        raise ValidationError("DECOMPRESSED_BYTES: source and BGZF differ")
    probe_count = _random_access_probes(payload, probes)
    return {
        "sequence_count": sequence_count,
        "names_lengths_sha256": names_digest,
        "derived_decompressed_sha256": derived_sha,
        "derived_decompressed_bytes": derived_bytes,
        "random_access_probe_count": probe_count,
        "bgzf_block_count": blocks,
        "gzi_entry_count": gzi_entries,
    }


def _index_payload(payload: Path) -> None:
    _run(["bgzip", "-r", str(payload)])
    _run(["samtools", "faidx", str(payload)])
    for suffix in (".gzi", ".fai"):
        target = Path(str(payload) + suffix)
        if not target.is_file() or target.stat().st_size == 0:
            raise ValidationError(f"INDEX_MISSING: {target}")


def _status_path(derived_root: Path, row: InventoryRow) -> Path:
    return derived_root / "status" / f"{row.task_index:06d}.{row.accession_version}.json"


def _final_paths(derived_root: Path, row: InventoryRow) -> tuple[Path, Path, Path, Path]:
    payload = derived_root / row.derived_relative_path
    return payload, Path(str(payload) + ".gzi"), Path(str(payload) + ".fai"), payload.parent / "provenance.json"


def _existing_complete(derived_root: Path, row: InventoryRow, *, deep: bool = False) -> dict[str, object] | None:
    payload, gzi, fai, provenance = _final_paths(derived_root, row)
    if not all(path.is_file() for path in (payload, gzi, fai, provenance)):
        return None
    try:
        data = json.loads(provenance.read_text())
    except (OSError, json.JSONDecodeError):
        return None
    if (
        data.get("source_relative_path") != row.source_relative_path
        or data.get("source_compressed_sha256") != row.source_compressed_sha256
        or data.get("state") not in PROMOTED_STATES
    ):
        return None
    for label, path in (("derived_bgzf", payload), ("gzi", gzi), ("fai", fai)):
        if data.get(f"{label}_bytes") != path.stat().st_size:
            return None
        if deep and data.get(f"{label}_sha256") != sha256_file(path):
            return None
    if deep:
        _run(["bgzip", "-t", str(payload)])
    validate_bgzf_structure(payload)
    validate_gzi(gzi, payload.stat().st_size)
    return data


def _copy_file_fsynced(source: Path, destination: Path) -> None:
    with source.open("rb") as incoming, destination.open("xb") as outgoing:
        shutil.copyfileobj(incoming, outgoing, 8 * 1024**2)
        outgoing.flush()
        os.fsync(outgoing.fileno())


def _promote(
    derived_root: Path,
    row: InventoryRow,
    work: Path,
    provenance: dict[str, object],
) -> tuple[str, dict[str, object]]:
    payload = work / "payload.fa.gz"
    gzi = Path(str(payload) + ".gzi")
    fai = Path(str(payload) + ".fai")
    identity = str(provenance["source_decompressed_sha256"])
    cas = derived_root / "cas" / "decompressed-sha256" / identity[:2] / identity[2:4] / identity
    lock_path = derived_root / "locks" / f"{identity}.lock"
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    with lock_path.open("a+") as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
        if cas.exists():
            metadata = json.loads((cas / "provenance.json").read_text())
            for label, name in (("derived_bgzf", "payload.fa.gz"), ("gzi", "payload.fa.gz.gzi"), ("fai", "payload.fa.gz.fai")):
                candidate = cas / name
                if metadata.get(f"{label}_sha256") != sha256_file(candidate):
                    raise ValidationError(f"CAS_DIGEST: invalid existing {label} for {identity}")
            _run(["bgzip", "-t", str(cas / "payload.fa.gz")])
            state = "reused"
            provenance["cas_reused"] = True
        else:
            cas.parent.mkdir(parents=True, exist_ok=True)
            (derived_root / "staging").mkdir(parents=True, exist_ok=True)
            stage = derived_root / "staging" / f"cas.{identity}.{uuid.uuid4().hex}.tmp"
            stage.mkdir()
            try:
                _copy_file_fsynced(payload, stage / "payload.fa.gz")
                _copy_file_fsynced(gzi, stage / "payload.fa.gz.gzi")
                _copy_file_fsynced(fai, stage / "payload.fa.gz.fai")
                atomic_write_json(stage / "provenance.json", provenance)
                fsync_path(stage)
                os.replace(stage, cas)
                fsync_path(cas.parent)
            finally:
                if stage.exists():
                    shutil.rmtree(stage)
            state = "reused" if row.source_format == "bgzf" else "converted"
            provenance["cas_reused"] = False
    final_payload = derived_root / row.derived_relative_path
    final_dir = final_payload.parent
    parent = final_dir.parent
    parent.mkdir(parents=True, exist_ok=True)
    if final_dir.exists():
        existing = _existing_complete(derived_root, row, deep=True)
        if existing is None:
            raise ValidationError(f"PROMOTION_CONFLICT: incomplete or stale {final_dir}")
        return "reused", existing
    stage_dir = parent / f".{final_dir.name}.{uuid.uuid4().hex}.tmp"
    stage_dir.mkdir()
    try:
        os.link(cas / "payload.fa.gz", stage_dir / final_payload.name)
        os.link(cas / "payload.fa.gz.gzi", stage_dir / f"{final_payload.name}.gzi")
        os.link(cas / "payload.fa.gz.fai", stage_dir / f"{final_payload.name}.fai")
        provenance.update(
            {
                "state": state,
                "cas_path": str(cas),
                "derived_path": str(final_payload),
                "promoted_at_utc": utc_now(),
            }
        )
        atomic_write_json(stage_dir / "provenance.json", provenance)
        fsync_path(stage_dir)
        os.replace(stage_dir, final_dir)
        fsync_path(parent)
    finally:
        if stage_dir.exists():
            shutil.rmtree(stage_dir)
    return state, provenance


def _reason_code(error: BaseException) -> str:
    if isinstance(error, PipelineError):
        detail = str(error).split(":", 1)[0]
        return detail if detail.isupper() and " " not in detail else error.reason_code
    if isinstance(error, OSError):
        if error.errno == errno.ENOSPC:
            return "NO_SPACE"
        if error.errno in (errno.EIO, errno.ESTALE):
            return "SHARED_IO"
    return "UNEXPECTED"


def run_worker(inventory: Path, derived_root: Path, task_index: int, scratch_root: Path) -> dict[str, object]:
    rows = read_inventory(inventory)
    if task_index < 0 or task_index >= len(rows):
        raise InventoryError(f"task index out of range: {task_index}")
    row = rows[task_index]
    status_path = _status_path(derived_root, row)
    attempts = 1
    if status_path.exists():
        try:
            attempts = int(json.loads(status_path.read_text()).get("attempts", 0)) + 1
        except (ValueError, json.JSONDecodeError):
            pass
    started = time.monotonic()
    self_before = resource.getrusage(resource.RUSAGE_SELF)
    child_before = resource.getrusage(resource.RUSAGE_CHILDREN)
    base_status: dict[str, object] = {
        **row.as_dict(),
        "attempts": attempts,
        "started_at_utc": utc_now(),
        "slurm_job_id": os.environ.get("SLURM_JOB_ID", ""),
        "slurm_array_task_id": os.environ.get("SLURM_ARRAY_TASK_ID", ""),
        "hostname": os.uname().nodename,
    }
    existing = _existing_complete(derived_root, row, deep=True)
    if existing is not None:
        status = {**base_status, **existing, "state": "reused", "reason_code": "VERIFIED_RESUME", "updated_at_utc": utc_now()}
        atomic_write_json(status_path, status)
        return status
    scratch_root.mkdir(parents=True, exist_ok=True)
    free = shutil.disk_usage(scratch_root).free
    if free < row.scratch_bytes:
        error = ScratchError(f"SCRATCH_INSUFFICIENT: require {row.scratch_bytes}, have {free}")
        failed = {**base_status, "state": "failed", "reason_code": _reason_code(error), "error": str(error), "updated_at_utc": utc_now()}
        atomic_write_json(status_path, failed)
        raise error
    work = Path(tempfile.mkdtemp(prefix=f"vgp-bgzf-{row.task_index:06d}-", dir=scratch_root))
    try:
        if sha256_file(Path(row.source_path)) != row.source_compressed_sha256:
            raise ValidationError("SOURCE_SHA256_DRIFT: mirrored source changed")
        captured = (
            _capture_existing_bgzf(row, work)
            if row.source_format == "bgzf"
            else _capture_source_and_compress(row, work, max(1, row.cpus - 1))
        )
        payload = Path(captured["payload"])
        _index_payload(payload)
        validation = validate_triplet(
            payload,
            Path(captured["source_sequences"]),
            str(captured["source_decompressed_sha256"]),
            int(captured["source_decompressed_bytes"]),
        )
        payload_sha = sha256_file(payload)
        gzi = Path(str(payload) + ".gzi")
        fai = Path(str(payload) + ".fai")
        provenance: dict[str, object] = {
            "schema_version": 1,
            "inventory_id": row.inventory_id,
            "accession_version": row.accession_version,
            "source_relative_path": row.source_relative_path,
            "source_path": row.source_path,
            "source_format": row.source_format,
            "source_bytes": row.source_bytes,
            "source_compressed_sha256": row.source_compressed_sha256,
            "source_decompressed_sha256": captured["source_decompressed_sha256"],
            "source_decompressed_bytes": captured["source_decompressed_bytes"],
            "derived_bgzf_sha256": payload_sha,
            "derived_bgzf_bytes": payload.stat().st_size,
            "gzi_sha256": sha256_file(gzi),
            "gzi_bytes": gzi.stat().st_size,
            "fai_sha256": sha256_file(fai),
            "fai_bytes": fai.stat().st_size,
            "bgzip_version": _run(["bgzip", "--version"]).stdout.decode().splitlines()[0],
            "samtools_version": _run(["samtools", "--version"]).stdout.decode().splitlines()[0],
            **validation,
        }
        state, provenance = _promote(derived_root, row, work, provenance)
        elapsed = time.monotonic() - started
        self_after = resource.getrusage(resource.RUSAGE_SELF)
        child_after = resource.getrusage(resource.RUSAGE_CHILDREN)
        cpu_seconds = (
            self_after.ru_utime - self_before.ru_utime + self_after.ru_stime - self_before.ru_stime
            + child_after.ru_utime - child_before.ru_utime + child_after.ru_stime - child_before.ru_stime
        )
        peak_scratch = sum(path.stat().st_size for path in work.iterdir() if path.is_file())
        status = {
            **base_status,
            **provenance,
            "state": state,
            "reason_code": "SOURCE_BGZF_VERIFIED" if row.source_format == "bgzf" else "STREAM_CONVERTED_VERIFIED",
            "elapsed_seconds": round(elapsed, 6),
            "cpu_seconds": round(cpu_seconds, 6),
            "peak_memory_bytes": max(self_after.ru_maxrss, child_after.ru_maxrss) * 1024,
            "peak_scratch_bytes": peak_scratch,
            "compression_ratio": payload.stat().st_size / int(captured["source_decompressed_bytes"]),
            "updated_at_utc": utc_now(),
        }
        atomic_write_json(status_path, status)
        return status
    except BaseException as error:
        failed = {
            **base_status,
            "state": "failed",
            "reason_code": _reason_code(error),
            "error": f"{type(error).__name__}: {error}",
            "elapsed_seconds": round(time.monotonic() - started, 6),
            "peak_scratch_bytes": sum(path.stat().st_size for path in work.iterdir() if path.is_file()),
            "updated_at_utc": utc_now(),
        }
        atomic_write_json(status_path, failed)
        raise
    finally:
        shutil.rmtree(work, ignore_errors=True)


def finalize(inventory: Path, derived_root: Path, shared_manifest: Path, summary_path: Path) -> dict[str, object]:
    rows = read_inventory(inventory)
    records: list[dict[str, object]] = []
    counts: Counter[str] = Counter()
    reason_counts: Counter[str] = Counter()
    for row in rows:
        status_path = _status_path(derived_root, row)
        if status_path.is_file():
            try:
                status = json.loads(status_path.read_text())
            except json.JSONDecodeError:
                status = {"state": "failed", "reason_code": "STATUS_INVALID_JSON", "attempts": 0}
        else:
            status = {"state": "failed", "reason_code": "NO_WORKER_STATUS", "attempts": 0}
        state = str(status.get("state", "failed"))
        if state not in PROMOTED_STATES | {"failed", "excluded"}:
            state = "failed"
            status["reason_code"] = "STATUS_NONTERMINAL"
        if state in PROMOTED_STATES and _existing_complete(derived_root, row) is None:
            state = "failed"
            status["reason_code"] = "PROMOTED_TRIPLET_INVALID"
        payload, gzi, fai, provenance = _final_paths(derived_root, row)
        record = {
            **row.as_dict(),
            **status,
            "state": state,
            "derived_path": str(payload) if payload.exists() else "",
            "gzi_path": str(gzi) if gzi.exists() else "",
            "fai_path": str(fai) if fai.exists() else "",
            "provenance_path": str(provenance) if provenance.exists() else "",
        }
        records.append(record)
        counts[state] += 1
        reason_counts[str(record.get("reason_code", ""))] += 1
    _write_tsv(shared_manifest, SHARED_FIELDS, records)
    promoted = [row for row in records if row["state"] in PROMOTED_STATES]
    total_decompressed = sum(int(row.get("source_decompressed_bytes", 0)) for row in promoted)
    total_derived = sum(int(row.get("derived_bgzf_bytes", 0)) for row in promoted)
    summary: dict[str, object] = {
        "schema_version": 1,
        "generated_at_utc": utc_now(),
        "closed_world": len(records) == len(rows) and sum(counts.values()) == len(rows),
        "inventory_rows": len(rows),
        "counts": dict(sorted(counts.items())),
        "reason_counts": dict(sorted(reason_counts.items())),
        "total_source_compressed_bytes": sum(row.source_bytes for row in rows),
        "total_source_decompressed_bytes": total_decompressed,
        "total_derived_bgzf_bytes": total_derived,
        "total_index_bytes": sum(int(row.get("gzi_bytes", 0)) + int(row.get("fai_bytes", 0)) for row in promoted),
        "compression_ratio_bgzf_to_decompressed": total_derived / total_decompressed if total_decompressed else None,
        "cpu_hours": sum(float(row.get("cpu_seconds", 0)) for row in records) / 3600,
        "elapsed_hours_sum": sum(float(row.get("elapsed_seconds", 0)) for row in records) / 3600,
        "peak_memory_bytes_max": max((int(row.get("peak_memory_bytes", 0)) for row in records), default=0),
        "peak_scratch_bytes_max": max((int(row.get("peak_scratch_bytes", 0)) for row in records), default=0),
        "failures": counts.get("failed", 0),
        "retries": sum(max(0, int(row.get("attempts", 0)) - 1) for row in records),
        "config_variable": CONFIG_VARIABLE,
        "config_value": str((derived_root / "objects").resolve()),
        "inventory_path": str(inventory.resolve()),
        "inventory_sha256": sha256_file(inventory),
        "shared_manifest_path": str(shared_manifest.resolve()),
    }
    atomic_write_json(summary_path, summary)
    shared_copy = derived_root / "manifest.tsv"
    shutil.copyfile(shared_manifest, shared_copy)
    atomic_write_json(derived_root / "summary.json", summary)
    return summary


def _parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    inventory = subparsers.add_parser("inventory")
    inventory.add_argument("--mirror-manifest", type=Path, default=DEFAULT_MIRROR_MANIFEST)
    inventory.add_argument("--derived-root", type=Path, default=DEFAULT_DERIVED_ROOT)
    inventory.add_argument("--output", type=Path, default=DEFAULT_INVENTORY)
    worker = subparsers.add_parser("worker")
    worker.add_argument("--inventory", type=Path, default=DEFAULT_INVENTORY)
    worker.add_argument("--derived-root", type=Path, default=DEFAULT_DERIVED_ROOT)
    worker.add_argument("--task-index", type=int, required=True)
    worker.add_argument("--scratch-root", type=Path, required=True)
    final = subparsers.add_parser("finalize")
    final.add_argument("--inventory", type=Path, default=DEFAULT_INVENTORY)
    final.add_argument("--derived-root", type=Path, default=DEFAULT_DERIVED_ROOT)
    final.add_argument("--shared-manifest", type=Path, default=DEFAULT_SHARED_MANIFEST)
    final.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = _parse_args(argv)
    if args.command == "inventory":
        output = write_inventory(args.mirror_manifest, args.derived_root, args.output)
        print(output)
    elif args.command == "worker":
        status = run_worker(args.inventory, args.derived_root, args.task_index, args.scratch_root)
        print(json.dumps(status, sort_keys=True))
    elif args.command == "finalize":
        summary = finalize(args.inventory, args.derived_root, args.shared_manifest, args.summary)
        print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except PipelineError as error:
        print(f"{_reason_code(error)}: {error}", file=sys.stderr)
        raise SystemExit(2)
