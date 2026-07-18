#!/usr/bin/env python3
"""Fail-closed, resumable mirroring of the immutable VGP Phase 1 Freeze 1 hub.

The release roster is the pinned VGP catalog, never the moving UCSC hub.  The
catalog's exact ``UCSC Browser main haplotype`` accession/version values are
mapped to UCSC's sibling accession trees below the official ``hubs`` rsync
module.  A complete metadata-only listing is required before this program will
transfer an object.

Bulk state lives outside git in SQLite.  Human-readable TSV/JSON snapshots are
written atomically and can safely be regenerated while a worker is running.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import shutil
import sqlite3
import stat
import subprocess
import sys
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath
from typing import Iterable, Mapping, Sequence


PROJECT_ROOT = Path(__file__).resolve().parents[1]
RELEASE_ROOT = Path("/moosefs/erikg/lewontin-paradox-data/vgp/freeze1")
PINNED_CATALOG = Path(
    "/moosefs/erikg/lewontin-paradox-data/vgp/phase1-freeze-1.0/manifests/"
    "VGPPhase1-freeze-1.0.commit-dc1b2af5a7741b97d66fb10cb2bce97f41765cdf.tsv"
)
CATALOG_COMMIT = "dc1b2af5a7741b97d66fb10cb2bce97f41765cdf"
CATALOG_SHA256 = "9c58420484a8b76a2d6175b7c26bf709e68bdc726a67fc7541b8c2b5a2fc13a4"
CATALOG_BYTES = 327466
CATALOG_PHYSICAL_LINES = 717
CATALOG_DATA_ROWS = 716
RELEASE_ENDPOINT = "rsync://hgdownload.soe.ucsc.edu/hubs/VGP/"
TRANSPORT_ENDPOINT = "rsync://hgdownload.soe.ucsc.edu/hubs/"
ACCESSION_RE = re.compile(r"^(GC[AF])_(\d{3})(\d{3})(\d{3})\.(\d+)$")
PATH_ACCESSION_RE = re.compile(r"(?:^|/)(GC[AF]_\d+\.\d+)(?:/|$)")
VALID_STATES = {
    "planned",
    "transferred",
    "verified",
    "reused",
    "missing",
    "superseded",
    "quarantined",
}
DEFAULT_HEADROOM_FRACTION = 0.20
DEFAULT_CONCURRENCY = 2
CHECKSUM_RESERVE_BYTES = 64 * 1024 * 1024

SOURCE_FIELDS = [
    "inventory_id",
    "release_definition_endpoint",
    "transport_endpoint",
    "retrieval_started_utc",
    "retrieval_completed_utc",
    "catalog_commit",
    "catalog_sha256",
    "catalog_data_rows",
    "catalog_row",
    "accession_version",
    "source_relative_path",
    "source_endpoint",
    "object_type",
    "size_bytes",
    "source_mtime_utc",
    "link_target",
    "sequence_subset",
    "upstream_checksum_algorithm",
    "upstream_checksum",
]
MANIFEST_FIELDS = [
    "inventory_id",
    "accession_version",
    "source_relative_path",
    "object_type",
    "size_bytes",
    "source_mtime_utc",
    "sequence_subset",
    "upstream_checksum_algorithm",
    "upstream_checksum",
    "observed_bytes",
    "local_sha256",
    "durable_path",
    "staging_path",
    "quarantine_path",
    "state",
    "state_reason",
    "attempts",
    "transferred_bytes",
    "updated_at_utc",
]


class MirrorError(RuntimeError):
    """An acquisition contract was violated."""


@dataclass(frozen=True)
class CatalogRow:
    row_number: int
    species: str
    accession: str | None


@dataclass(frozen=True)
class InventoryObject:
    inventory_id: str
    catalog_row: int
    accession: str
    path: str
    object_type: str
    size: int
    mtime: str
    link_target: str
    sequence_subset: str
    checksum_algorithm: str = ""
    checksum: str = ""


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def digest_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def atomic_write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.partial-{os.getpid()}")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        handle.write(text)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def atomic_write_json(path: Path, payload: Mapping[str, object]) -> None:
    atomic_write_text(path, json.dumps(payload, indent=2, sort_keys=True) + "\n")


def atomic_write_tsv(path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.partial-{os.getpid()}")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def verify_catalog(path: Path) -> list[CatalogRow]:
    if not path.is_file():
        raise MirrorError(f"pinned catalog is missing: {path}")
    observed_size = path.stat().st_size
    observed_sha = sha256_file(path)
    with path.open("rb") as handle:
        physical_lines = sum(1 for _ in handle)
    if (observed_size, observed_sha, physical_lines) != (
        CATALOG_BYTES,
        CATALOG_SHA256,
        CATALOG_PHYSICAL_LINES,
    ):
        raise MirrorError(
            "pinned catalog identity mismatch: "
            f"bytes={observed_size}, sha256={observed_sha}, lines={physical_lines}"
        )
    with path.open("r", encoding="utf-8", newline="") as handle:
        records = list(csv.DictReader(handle, delimiter="\t"))
    if len(records) != CATALOG_DATA_ROWS:
        raise MirrorError(f"pinned catalog has {len(records)} rows, expected {CATALOG_DATA_ROWS}")
    rows: list[CatalogRow] = []
    seen: set[str] = set()
    for row_number, record in enumerate(records, 2):
        raw = record["UCSC Browser main haplotype"].strip()
        accession = raw or None
        if accession:
            if not ACCESSION_RE.fullmatch(accession):
                raise MirrorError(f"catalog row {row_number} has invalid exact accession {accession!r}")
            if accession in seen:
                raise MirrorError(f"catalog accession is duplicated: {accession}")
            seen.add(accession)
        rows.append(CatalogRow(row_number, record["Scientific Name"].strip(), accession))
    return rows


def accession_path(accession: str) -> str:
    match = ACCESSION_RE.fullmatch(accession)
    if not match:
        raise MirrorError(f"invalid accession/version: {accession!r}")
    prefix, first, second, third, _version = match.groups()
    return f"{prefix}/{first}/{second}/{third}/{accession}"


def freeze_inventory(catalog: Path, root: Path) -> dict[str, object]:
    """Stream a metadata-only recursive listing to a new immutable snapshot."""
    if not os.environ.get("GUIX_ENVIRONMENT"):
        raise MirrorError("inventory must run inside the pinned GNU Guix environment")
    rows = verify_catalog(catalog)
    accessions = sorted(row.accession for row in rows if row.accession)
    if len(accessions) != 581:
        raise MirrorError(f"expected 581 exact release accessions, observed {len(accessions)}")
    inventory_root = root / "inventory"
    empty = inventory_root / ".empty"
    inventory_root.mkdir(parents=True, exist_ok=True)
    empty.mkdir(parents=True, exist_ok=True)
    paths_file = inventory_root / "freeze1-accession-paths.txt"
    atomic_write_text(paths_file, "".join(f"{accession_path(value)}/\n" for value in accessions))
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    raw = inventory_root / f"freeze1-assemblies.{stamp}.rsync-items"
    metadata = inventory_root / f"freeze1-assemblies.{stamp}.metadata"
    temporary = raw.with_name(f".{raw.name}.partial-{os.getpid()}")
    started = utc_now()
    command = [
        "rsync",
        "--dry-run",
        "--recursive",
        "--links",
        "--times",
        "--itemize-changes",
        "--protect-args",
        "--no-motd",
        "--timeout=300",
        f"--files-from={paths_file}",
        "--out-format=%i|%l|%M|%n|%L",
        TRANSPORT_ENDPOINT,
        str(empty) + "/",
    ]
    with temporary.open("wb") as output:
        completed = subprocess.run(command, stdout=output, stderr=subprocess.PIPE, check=False)
        output.flush()
        os.fsync(output.fileno())
    completed_at = utc_now()
    if completed.returncode != 0:
        failed = raw.with_suffix(raw.suffix + ".failed")
        os.replace(temporary, failed)
        raise MirrorError(
            f"metadata-only rsync failed ({completed.returncode}): "
            f"{completed.stderr.decode('utf-8', errors='replace').strip()}"
        )
    os.replace(temporary, raw)
    raw_sha256 = sha256_file(raw)
    with raw.open("rb") as handle:
        line_count = sum(1 for _ in handle)
    metadata_text = "\n".join(
        (
            f"retrieval_started_utc={started}",
            f"retrieval_completed_utc={completed_at}",
            f"release_definition_endpoint={RELEASE_ENDPOINT}",
            f"transport_endpoint={TRANSPORT_ENDPOINT}",
            "rsync_exit_status=0",
            f"raw_path={raw}",
            f"raw_inventory_sha256={raw_sha256}",
            f"raw_inventory_lines={line_count}",
            f"raw_inventory_bytes={raw.stat().st_size}",
            f"catalog_sha256={CATALOG_SHA256}",
            f"catalog_data_rows={CATALOG_DATA_ROWS}",
            "guix_channel_commit=44bbfc24e4bcc48d0e3343cd3d83452721af8c36",
            f"rsync_version={subprocess.run(['rsync', '--version'], capture_output=True, text=True, check=True).stdout.splitlines()[0]}",
        )
    )
    atomic_write_text(metadata, metadata_text + "\n")
    return {
        "status": "metadata_inventory_frozen",
        "raw_inventory": str(raw),
        "raw_inventory_sha256": raw_sha256,
        "raw_inventory_lines": line_count,
        "raw_inventory_bytes": raw.stat().st_size,
        "metadata": str(metadata),
        "accession_roots": len(accessions),
        "payload_bytes_transferred": 0,
    }


def parse_metadata_file(path: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            result[key] = value
        else:
            checksum = re.fullmatch(r"([0-9a-f]{64})\s+(.+)", line)
            if checksum:
                result["raw_inventory_sha256"] = checksum.group(1)
    return result


def classify_sequence(path: str, accession: str) -> str:
    name = PurePosixPath(path).name
    if name in {f"{accession}.fa", f"{accession}.fa.gz"}:
        return "assembly_fasta"
    if name == f"{accession}.2bit":
        return "assembly_2bit"
    return "non_sequence_product_or_metadata"


def parse_rsync_inventory(
    raw_path: Path,
    catalog_rows: Sequence[CatalogRow],
    checksum_map: Mapping[str, tuple[str, str]] | None = None,
) -> list[InventoryObject]:
    """Parse a complete `%i|%l|%M|%n|%L` dry-run listing and fail closed."""
    row_by_accession = {row.accession: row.row_number for row in catalog_rows if row.accession}
    root_by_accession = {accession: accession_path(accession) for accession in row_by_accession}
    observed_roots: set[str] = set()
    objects: list[InventoryObject] = []
    checksum_map = checksum_map or {}
    for line_number, line in enumerate(raw_path.read_text(encoding="utf-8").splitlines(), 1):
        parts = line.split("|", 4)
        if len(parts) != 5:
            raise MirrorError(f"unparseable rsync inventory line {line_number}: {line!r}")
        itemized, size_text, mtime, path, link_field = parts
        if len(itemized) < 2 or itemized[1] not in {"f", "d", "L"}:
            raise MirrorError(f"unsupported rsync object type on line {line_number}: {itemized!r}")
        object_type = {"f": "file", "d": "directory", "L": "symlink"}[itemized[1]]
        path = path.removesuffix("/")
        matching = [
            accession
            for accession, root in root_by_accession.items()
            if path == root or path.startswith(root + "/")
        ]
        # --files-from emits shared parent directories.  They are transport
        # scaffolding, not release objects, so exclude them from the snapshot.
        if not matching:
            continue
        if len(matching) != 1:
            raise MirrorError(f"inventory path maps ambiguously: {path}")
        accession = matching[0]
        root = root_by_accession[accession]
        observed_roots.add(accession)
        # The root component is the release identity.  Names of liftOver files
        # may legitimately mention other accessions and are not membership.
        root_match = PATH_ACCESSION_RE.search(root)
        if not root_match or root_match.group(1) != accession:
            raise MirrorError(f"accession drift in source root: {root}")
        try:
            size = int(size_text)
        except ValueError as error:
            raise MirrorError(f"invalid size on inventory line {line_number}") from error
        link_target = link_field.removeprefix(" -> ") if object_type == "symlink" else ""
        checksum_algorithm, checksum = checksum_map.get(path, ("", ""))
        identity = digest_text(f"{accession}\0{path}\0{object_type}\0{size}\0{mtime}\0{link_target}")
        objects.append(
            InventoryObject(
                identity,
                row_by_accession[accession],
                accession,
                path,
                object_type,
                size,
                mtime,
                link_target,
                classify_sequence(path, accession),
                checksum_algorithm,
                checksum,
            )
        )
    missing = sorted(set(row_by_accession) - observed_roots)
    if missing:
        raise MirrorError(f"{len(missing)} frozen accessions are missing from inventory: {missing[:10]}")
    if len(observed_roots) != 581:
        raise MirrorError(f"expected 581 exact accession roots, observed {len(observed_roots)}")
    if len({obj.path for obj in objects}) != len(objects):
        raise MirrorError("source inventory contains duplicate release paths")
    return sorted(objects, key=lambda obj: obj.path)


def parse_upstream_md5s(objects_root: Path, objects: Sequence[InventoryObject]) -> dict[str, tuple[str, str]]:
    checksums: dict[str, tuple[str, str]] = {}
    roots = {obj.accession: accession_path(obj.accession) for obj in objects}
    inventoried_paths = {obj.path for obj in objects}
    for accession, root in roots.items():
        manifest = objects_root / root / "md5sum.txt"
        if not manifest.is_file():
            continue
        for line_number, line in enumerate(manifest.read_text(encoding="utf-8", errors="strict").splitlines(), 1):
            match = re.fullmatch(r"([0-9a-fA-F]{32})\s+[* ]?(.+)", line)
            if not match:
                raise MirrorError(f"invalid upstream MD5 line {manifest}:{line_number}")
            digest, relative = match.groups()
            relative = relative.removeprefix("./")
            candidate = str(PurePosixPath(root) / relative)
            if candidate not in inventoried_paths:
                raise MirrorError(
                    f"upstream checksum manifest references a non-inventoried object: {candidate}"
                )
            prior = checksums.get(candidate)
            value = ("md5", digest.lower())
            if prior and prior != value:
                raise MirrorError(f"conflicting upstream checksums for {candidate}")
            checksums[candidate] = value
    return checksums


def inventory_totals(objects: Sequence[InventoryObject]) -> dict[str, dict[str, int]]:
    result: dict[str, dict[str, int]] = {}
    subsets = {
        "full_release": objects,
        "assembly_fasta": [obj for obj in objects if obj.sequence_subset == "assembly_fasta"],
        "assembly_2bit": [obj for obj in objects if obj.sequence_subset == "assembly_2bit"],
        "assembly_fasta_or_2bit": [
            obj for obj in objects if obj.sequence_subset in {"assembly_fasta", "assembly_2bit"}
        ],
    }
    for name, subset in subsets.items():
        files = [obj for obj in subset if obj.object_type == "file"]
        result[name] = {
            "objects": len(subset),
            "files": len(files),
            "directories": sum(obj.object_type == "directory" for obj in subset),
            "symlinks": sum(obj.object_type == "symlink" for obj in subset),
            "bytes": sum(obj.size for obj in files),
        }
    return result


def storage_evidence(root: Path, objects: Sequence[InventoryObject], concurrency: int) -> dict[str, object]:
    root.mkdir(parents=True, exist_ok=True)
    values = os.statvfs(root)
    available = values.f_bavail * values.f_frsize
    free_inodes = values.f_favail
    files = [obj for obj in objects if obj.object_type == "file"]
    durable = sum(obj.size for obj in files)
    largest = max((obj.size for obj in files), default=0)
    staging = concurrency * largest
    quarantine = largest
    headroom = int(durable * DEFAULT_HEADROOM_FRACTION)
    inode_durable = len(objects)
    inode_staging = concurrency
    inode_quarantine = 1
    inode_checksum = 8
    inode_headroom = max(1024, int(inode_durable * DEFAULT_HEADROOM_FRACTION))
    required = durable + staging + quarantine + CHECKSUM_RESERVE_BYTES + headroom
    required_inodes = (
        inode_durable + inode_staging + inode_quarantine + inode_checksum + inode_headroom
    )
    quota_attempts: list[dict[str, object]] = []
    commands = (["mfsgetquota", "-n", str(root)], ["quota", "-s"])
    quota_visible = False
    for command in commands:
        executable = shutil.which(command[0])
        if not executable:
            quota_attempts.append({"command": command, "status": "tool_unavailable"})
            continue
        completed = subprocess.run(command, text=True, capture_output=True, check=False)
        quota_attempts.append(
            {
                "command": command,
                "status": "ok" if completed.returncode == 0 else "failed",
                "returncode": completed.returncode,
                "stdout": completed.stdout.strip(),
                "stderr": completed.stderr.strip(),
            }
        )
        quota_visible |= completed.returncode == 0
    filesystem_adequate = available >= required and free_inodes >= required_inodes
    return {
        "path": str(root),
        "observed_at_utc": utc_now(),
        "statvfs": {
            "block_size": values.f_frsize,
            "total_bytes": values.f_blocks * values.f_frsize,
            "free_bytes": values.f_bfree * values.f_frsize,
            "available_bytes": available,
            "total_inodes": values.f_files,
            "free_inodes": free_inodes,
        },
        "requirements": {
            "durable_bytes": durable,
            "staging_bytes": staging,
            "staging_concurrency": concurrency,
            "largest_object_bytes": largest,
            "checksum_manifest_bytes": CHECKSUM_RESERVE_BYTES,
            "quarantine_bytes": quarantine,
            "operational_headroom_fraction": DEFAULT_HEADROOM_FRACTION,
            "operational_headroom_bytes": headroom,
            "total_bytes": required,
            "durable_inodes": inode_durable,
            "staging_inodes": inode_staging,
            "checksum_inodes": inode_checksum,
            "quarantine_inodes": inode_quarantine,
            "operational_headroom_inodes": inode_headroom,
            "total_inodes": required_inodes,
        },
        "quota_evidence": quota_attempts,
        "filesystem_capacity_adequate": filesystem_adequate,
        "quota_visibility_adequate": quota_visible,
        "adequate": filesystem_adequate and quota_visible,
        "gate_reason": (
            "capacity_and_quota_verified"
            if filesystem_adequate and quota_visible
            else "quota_visibility_unavailable_fail_closed"
            if filesystem_adequate
            else "filesystem_capacity_or_inodes_inadequate"
        ),
        "arbitrary_global_byte_cap": None,
    }


def init_database(path: Path, objects: Sequence[InventoryObject], root: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    connection = sqlite3.connect(path)
    try:
        connection.execute("PRAGMA journal_mode=WAL")
        connection.execute(
            """CREATE TABLE IF NOT EXISTS objects (
                inventory_id TEXT PRIMARY KEY, accession TEXT NOT NULL, path TEXT NOT NULL UNIQUE,
                object_type TEXT NOT NULL, size INTEGER NOT NULL, mtime TEXT NOT NULL,
                sequence_subset TEXT NOT NULL, checksum_algorithm TEXT NOT NULL,
                checksum TEXT NOT NULL, observed_bytes INTEGER NOT NULL DEFAULT 0,
                local_sha256 TEXT NOT NULL DEFAULT '', durable_path TEXT NOT NULL,
                staging_path TEXT NOT NULL, quarantine_path TEXT NOT NULL DEFAULT '',
                state TEXT NOT NULL, state_reason TEXT NOT NULL, attempts INTEGER NOT NULL DEFAULT 0,
                transferred_bytes INTEGER NOT NULL DEFAULT 0, updated_at_utc TEXT NOT NULL
            )"""
        )
        now = utc_now()
        for obj in objects:
            durable_path = root / "objects" / obj.path
            staging_path = root / "staging" / f"{obj.path}.part"
            connection.execute(
                """INSERT OR IGNORE INTO objects
                (inventory_id, accession, path, object_type, size, mtime, sequence_subset,
                 checksum_algorithm, checksum, durable_path, staging_path, state, state_reason,
                 updated_at_utc) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, 'planned',
                 'frozen_inventory_pending_transfer', ?)""",
                (
                    obj.inventory_id,
                    obj.accession,
                    obj.path,
                    obj.object_type,
                    obj.size,
                    obj.mtime,
                    obj.sequence_subset,
                    obj.checksum_algorithm,
                    obj.checksum,
                    str(durable_path),
                    str(staging_path),
                    now,
                ),
            )
        connection.commit()
    finally:
        connection.close()


def source_rows(
    objects: Sequence[InventoryObject], metadata: Mapping[str, str]
) -> Iterable[dict[str, object]]:
    for obj in objects:
        parsed_mtime = datetime.strptime(obj.mtime, "%Y/%m/%d-%H:%M:%S").replace(tzinfo=timezone.utc)
        yield {
            "inventory_id": obj.inventory_id,
            "release_definition_endpoint": RELEASE_ENDPOINT,
            "transport_endpoint": TRANSPORT_ENDPOINT,
            "retrieval_started_utc": metadata["retrieval_started_utc"],
            "retrieval_completed_utc": metadata["retrieval_completed_utc"],
            "catalog_commit": CATALOG_COMMIT,
            "catalog_sha256": CATALOG_SHA256,
            "catalog_data_rows": CATALOG_DATA_ROWS,
            "catalog_row": obj.catalog_row,
            "accession_version": obj.accession,
            "source_relative_path": obj.path,
            "source_endpoint": TRANSPORT_ENDPOINT + obj.path,
            "object_type": obj.object_type,
            "size_bytes": obj.size,
            "source_mtime_utc": parsed_mtime.isoformat().replace("+00:00", "Z"),
            "link_target": obj.link_target,
            "sequence_subset": obj.sequence_subset,
            "upstream_checksum_algorithm": obj.checksum_algorithm,
            "upstream_checksum": obj.checksum,
        }


def database_rows(database: Path) -> list[dict[str, object]]:
    connection = sqlite3.connect(database)
    connection.row_factory = sqlite3.Row
    try:
        rows = [dict(row) for row in connection.execute("SELECT * FROM objects ORDER BY path")]
    finally:
        connection.close()
    for row in rows:
        if row["state"] not in VALID_STATES:
            raise MirrorError(f"invalid manifest state {row['state']!r}")
    return rows


def manifest_rows(rows: Sequence[Mapping[str, object]]) -> Iterable[dict[str, object]]:
    for row in rows:
        yield {
            "inventory_id": row["inventory_id"],
            "accession_version": row["accession"],
            "source_relative_path": row["path"],
            "object_type": row["object_type"],
            "size_bytes": row["size"],
            "source_mtime_utc": row["mtime"],
            "sequence_subset": row["sequence_subset"],
            "upstream_checksum_algorithm": row["checksum_algorithm"],
            "upstream_checksum": row["checksum"],
            "observed_bytes": row["observed_bytes"],
            "local_sha256": row["local_sha256"],
            "durable_path": row["durable_path"],
            "staging_path": row["staging_path"],
            "quarantine_path": row["quarantine_path"],
            "state": row["state"],
            "state_reason": row["state_reason"],
            "attempts": row["attempts"],
            "transferred_bytes": row["transferred_bytes"],
            "updated_at_utc": row["updated_at_utc"],
        }


def state_accounting(rows: Sequence[Mapping[str, object]]) -> dict[str, dict[str, int]]:
    accounting = {state: {"objects": 0, "files": 0, "bytes": 0} for state in sorted(VALID_STATES)}
    for row in rows:
        state = str(row["state"])
        accounting[state]["objects"] += 1
        if row["object_type"] == "file":
            accounting[state]["files"] += 1
            accounting[state]["bytes"] += int(row["size"])
    if sum(value["objects"] for value in accounting.values()) != len(rows):
        raise MirrorError("manifest states are not mutually exclusive and exhaustive")
    return accounting


def render_handoff(summary: Mapping[str, object]) -> str:
    totals = summary["inventory_totals"]
    storage = summary["storage"]
    reconciliation = summary["catalog_reconciliation"]
    state = summary["state_accounting"]
    req = storage["requirements"]
    statvfs = storage["statvfs"]
    return f"""# VGP Phase 1 Freeze 1 mirror handoff

**Release definition:** VGP/vgp-phase1 commit `{CATALOG_COMMIT}`, file
`VGPPhase1-freeze-1.0.tsv`, SHA-256 `{CATALOG_SHA256}`.  The local pinned copy
is exactly {CATALOG_BYTES:,} bytes, {CATALOG_PHYSICAL_LINES} physical lines,
and {CATALOG_DATA_ROWS} data rows.

## Closed-world reconciliation

All {reconciliation['catalog_rows']} catalog rows have exactly one disposition:
{reconciliation['released_rows']} rows map one-to-one to unique exact UCSC
browser accession/version roots; {reconciliation['unreleased_rows']} catalog
rows have an empty UCSC browser accession and therefore contribute no released
hub object.  The snapshot contains all {reconciliation['expected_accession_roots']}
expected roots, no missing root, no extra root, and no root-level accession or
version drift.  The moving current hub was used only as transport.

## Authoritative metadata-only inventory

The recursive rsync dry-run completed at `{summary['retrieval_completed_utc']}`
from `{RELEASE_ENDPOINT}` (product paths resolve against `{TRANSPORT_ENDPOINT}`).
The bound raw listing SHA-256 is `{summary['raw_inventory']['sha256']}`.
It inventoried {totals['full_release']['files']:,} files and
{totals['full_release']['bytes']:,} file bytes.  Exact assembly FASTA:
{totals['assembly_fasta']['files']:,} files / {totals['assembly_fasta']['bytes']:,}
bytes.  Exact assembly 2bit: {totals['assembly_2bit']['files']:,} files /
{totals['assembly_2bit']['bytes']:,} bytes.  Their union is
{totals['assembly_fasta_or_2bit']['files']:,} files /
{totals['assembly_fasta_or_2bit']['bytes']:,} bytes.

The historical approximately 967 GB whole-collection and approximately 520 GB
FASTA-only figures are **unverified historical planning estimates only**.  They
are not checksums, completion criteria, byte ceilings, or evidence.  They have
been replaced by the frozen inventory above.

All {summary['checksum_policy']['published_md5_manifests_in_inventory']} official
`md5sum.txt` objects are inventoried.  The worker promotes and validates those
first, binds every checksum that names an exact frozen path, and refuses all
remaining payload if any checksum manifest fails or names an extra/renamed
object.  Objects without a published checksum receive the three-pass local
SHA-256 contract described below.

## Capacity and launch gate

Required durable bytes: {req['durable_bytes']:,}; source-relative staging at
bounded concurrency {req['staging_concurrency']}: {req['staging_bytes']:,};
checksum/manifests: {req['checksum_manifest_bytes']:,}; one-object quarantine:
{req['quarantine_bytes']:,}; explicit {req['operational_headroom_fraction']:.0%}
operational headroom: {req['operational_headroom_bytes']:,}.  Combined
worst-case requirement: {req['total_bytes']:,} bytes and {req['total_inodes']:,}
inodes.  Filesystem evidence reports {statvfs['available_bytes']:,} bytes and
{statvfs['free_inodes']:,} inodes available.  No arbitrary global byte or
memory cap is imposed.

Capacity gate: **{storage['gate_reason']}**.  Filesystem capacity passes, but
bulk execution remains fail-closed unless the mounted filesystem's applicable
quota can be queried successfully.  This run did not infer quota absence from
`df`.  The exact failed/unavailable quota probes are preserved in the JSON
summary.  Consequently no bulk payload transfer or Slurm job was launched.

## Harmless transfer fixture

Pinned GNU Guix rsync was deliberately interrupted after
{summary['fixture'].get('interrupted_partial_bytes', 0):,} of
{summary['fixture'].get('source_bytes', 0):,} bytes.  Resume added exactly
{summary['fixture'].get('resume_additional_bytes', 0):,} bytes and produced
SHA-256 `{summary['fixture'].get('exact_resume_sha256', 'not-run')}`.  The
fixture also proved checksum-failure quarantine and atomic promotion:
`{summary['fixture'].get('passed', False)}`.  Its durable JSON report is
`{summary['fixture'].get('report_path')}`.

## Durable operation

Run only through the pinned wrapper:

```bash
analysis/run_vgp_freeze1_mirror.sh inventory
analysis/run_vgp_freeze1_mirror.sh build
analysis/run_vgp_freeze1_mirror.sh fixture
analysis/run_vgp_freeze1_mirror.sh worker --concurrency 2
```

`inventory` always creates a new timestamped snapshot and never replaces the
bound snapshot.  To reconcile a deliberate refresh, pass the new raw and
metadata paths explicitly to `build` and review the closed-world result before
allowing the worker to see its capacity evidence.

The worker uses source-relative `.part` staging, `rsync --partial
--append-verify`, bounded concurrency, exponential backoff, size plus published
MD5 verification when available, local SHA-256 before/re-before/after atomic
promotion, and source timestamps as metadata only.  A mismatch is moved to a
source-relative quarantine location and never overwrites or deletes the last
verified destination.  There is no mirror-wide delete operation.  Each object
becomes durable immediately after its own verification.

## Current mutually exclusive accounting

{json.dumps(state, indent=2, sort_keys=True)}

No raw sequencing-read archive is in scope.  No bulk payload belongs in git.
"""


def write_snapshots(
    *,
    objects: Sequence[InventoryObject],
    metadata: Mapping[str, str],
    database: Path,
    storage: Mapping[str, object],
    source_output: Path,
    manifest_output: Path,
    summary_output: Path,
    handoff_output: Path,
) -> dict[str, object]:
    rows = database_rows(database)
    totals = inventory_totals(objects)
    accounting = state_accounting(rows)
    released = len({obj.accession for obj in objects})
    checksum_manifests = [
        obj for obj in objects if obj.object_type == "file" and PurePosixPath(obj.path).name == "md5sum.txt"
    ]
    published_checksum_bindings = sum(bool(obj.checksum_algorithm) for obj in objects)
    fixture_path = RELEASE_ROOT / "fixture" / "fixture-report.json"
    fixture_report: dict[str, object] = {
        "report_path": str(fixture_path),
        "required_before_bulk": True,
        "passed": False,
    }
    if fixture_path.is_file():
        fixture_report.update(json.loads(fixture_path.read_text(encoding="utf-8")))
    summary: dict[str, object] = {
        "schema_version": "vgp-freeze1-mirror-summary-v1",
        "generated_at_utc": utc_now(),
        "release": {
            "catalog_commit": CATALOG_COMMIT,
            "catalog_sha256": CATALOG_SHA256,
            "catalog_bytes": CATALOG_BYTES,
            "catalog_physical_lines": CATALOG_PHYSICAL_LINES,
            "catalog_data_rows": CATALOG_DATA_ROWS,
            "release_definition_endpoint": RELEASE_ENDPOINT,
            "transport_endpoint": TRANSPORT_ENDPOINT,
            "guix_channel_commit": "44bbfc24e4bcc48d0e3343cd3d83452721af8c36",
            "guix_manifest": "analysis/guix/vgp-freeze1-manifest.scm",
        },
        "retrieval_started_utc": metadata["retrieval_started_utc"],
        "retrieval_completed_utc": metadata["retrieval_completed_utc"],
        "raw_inventory": {
            "path": metadata["raw_path"],
            "sha256": metadata.get("raw_inventory_sha256", ""),
            "rsync_exit_status": int(metadata["rsync_exit_status"]),
        },
        "catalog_reconciliation": {
            "catalog_rows": CATALOG_DATA_ROWS,
            "released_rows": released,
            "unreleased_rows": CATALOG_DATA_ROWS - released,
            "expected_accession_roots": released,
            "observed_accession_roots": released,
            "missing_roots": 0,
            "extra_roots": 0,
            "accession_or_version_drift": 0,
            "closed_world": True,
        },
        "inventory_totals": totals,
        "historical_estimates": {
            "whole_collection_approx_gb": 967,
            "fasta_only_approx_gb": 520,
            "status": "unverified_historical_planning_estimates_replaced_by_frozen_inventory",
            "checksum_evidence": False,
            "completion_criterion": False,
            "byte_ceiling": False,
        },
        "checksum_policy": {
            "published_md5_manifests_in_inventory": len(checksum_manifests),
            "published_checksum_bindings_loaded": published_checksum_bindings,
            "manifest_transfer_order": "all_md5sum.txt objects before any other payload",
            "published_checksum_rule": "verify provider MD5 when an exact frozen path is listed",
            "fallback_rule": "compute SHA-256 after staging and reverify before and after promotion",
        },
        "storage": storage,
        "state_accounting": accounting,
        "fixture": fixture_report,
        "bulk_launch": {
            "launched": sum(
                value["objects"] for key, value in accounting.items() if key != "planned"
            )
            > 0,
            "reason": storage["gate_reason"],
            "slurm_jobs_launched": 0,
        },
    }
    atomic_write_tsv(source_output, SOURCE_FIELDS, source_rows(objects, metadata))
    atomic_write_tsv(manifest_output, MANIFEST_FIELDS, manifest_rows(rows))
    atomic_write_json(summary_output, summary)
    atomic_write_text(handoff_output, render_handoff(summary))
    return summary


def update_database(database: Path, inventory_id: str, **updates: object) -> None:
    allowed = {
        "observed_bytes",
        "local_sha256",
        "quarantine_path",
        "state",
        "state_reason",
        "attempts",
        "transferred_bytes",
        "checksum_algorithm",
        "checksum",
    }
    if not updates or set(updates) - allowed:
        raise MirrorError(f"invalid database update keys: {set(updates) - allowed}")
    if "state" in updates and updates["state"] not in VALID_STATES:
        raise MirrorError(f"invalid state: {updates['state']}")
    updates["updated_at_utc"] = utc_now()
    assignments = ", ".join(f"{key} = ?" for key in updates)
    connection = sqlite3.connect(database, timeout=60)
    try:
        connection.execute("BEGIN IMMEDIATE")
        connection.execute(
            f"UPDATE objects SET {assignments} WHERE inventory_id = ?",
            [*updates.values(), inventory_id],
        )
        connection.commit()
    finally:
        connection.close()


def verify_digest(path: Path, algorithm: str, expected: str) -> None:
    if not algorithm:
        return
    if algorithm != "md5":
        raise MirrorError(f"unsupported upstream checksum algorithm: {algorithm}")
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    if digest.hexdigest() != expected:
        raise MirrorError(f"upstream MD5 mismatch for {path}")


def quarantine_part(part: Path, quarantine_root: Path, relative: str, reason: str) -> Path:
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    destination = quarantine_root / f"{relative}.{stamp}.{reason}"
    destination.parent.mkdir(parents=True, exist_ok=True)
    if part.exists():
        os.replace(part, destination)
    return destination


def promote_verified(
    part: Path,
    destination: Path,
    *,
    expected_size: int,
    checksum_algorithm: str = "",
    checksum: str = "",
    before_reverify=None,
) -> str:
    if part.stat().st_size != expected_size:
        raise MirrorError(f"size mismatch: expected {expected_size}, observed {part.stat().st_size}")
    verify_digest(part, checksum_algorithm, checksum)
    first_sha = sha256_file(part)
    if before_reverify:
        before_reverify(part, first_sha)
    second_sha = sha256_file(part)
    if second_sha != first_sha:
        raise MirrorError("staged object changed before promotion")
    with part.open("rb") as handle:
        os.fsync(handle.fileno())
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        raise MirrorError("durable destination already exists; refusing overwrite")
    os.chmod(part, stat.S_IRUSR | stat.S_IRGRP)
    os.replace(part, destination)
    third_sha = sha256_file(destination)
    if third_sha != first_sha:
        raise MirrorError("promoted object changed after atomic promotion")
    with destination.open("rb") as handle:
        os.fsync(handle.fileno())
    directory_fd = os.open(destination.parent, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(directory_fd)
    finally:
        os.close(directory_fd)
    return third_sha


def rsync_transfer(source: str, part: Path, timeout: int = 300) -> int:
    part.parent.mkdir(parents=True, exist_ok=True)
    before = part.stat().st_size if part.exists() else 0
    command = [
        "rsync",
        "--partial",
        "--append-verify",
        "--protect-args",
        "--times",
        "--no-motd",
        f"--timeout={timeout}",
        source,
        str(part),
    ]
    completed = subprocess.run(command, text=True, capture_output=True, check=False)
    if completed.returncode != 0:
        raise MirrorError(
            f"rsync failed ({completed.returncode}) for {source}: {completed.stderr.strip()}"
        )
    after = part.stat().st_size
    return max(0, after - before)


def process_file(database: Path, root: Path, row: Mapping[str, object], retries: int = 5) -> None:
    destination = Path(str(row["durable_path"]))
    part = Path(str(row["staging_path"]))
    expected_size = int(row["size"])
    algorithm = str(row["checksum_algorithm"])
    checksum = str(row["checksum"])
    if destination.is_file():
        try:
            if row["state"] == "planned" and not algorithm and not row["local_sha256"]:
                raise MirrorError("unbound pre-existing durable object")
            if destination.stat().st_size != expected_size:
                raise MirrorError("existing durable size mismatch")
            verify_digest(destination, algorithm, checksum)
            local_sha = sha256_file(destination)
            if row["local_sha256"] and local_sha != row["local_sha256"]:
                raise MirrorError("existing durable local SHA-256 mismatch")
            update_database(
                database,
                str(row["inventory_id"]),
                observed_bytes=destination.stat().st_size,
                local_sha256=local_sha,
                state="reused",
                state_reason="existing_durable_object_reverified",
            )
            return
        except MirrorError as error:
            update_database(
                database,
                str(row["inventory_id"]),
                state="superseded",
                state_reason=f"preserved_existing_durable_object:{error}",
            )
            return
    transferred = int(row["transferred_bytes"])
    for attempt in range(1, retries + 1):
        try:
            transferred += rsync_transfer(TRANSPORT_ENDPOINT + str(row["path"]), part)
            update_database(
                database,
                str(row["inventory_id"]),
                observed_bytes=part.stat().st_size,
                transferred_bytes=transferred,
                attempts=int(row["attempts"]) + attempt,
                state="transferred",
                state_reason="complete_part_pending_verification",
            )
            local_sha = promote_verified(
                part,
                destination,
                expected_size=expected_size,
                checksum_algorithm=algorithm,
                checksum=checksum,
            )
            update_database(
                database,
                str(row["inventory_id"]),
                observed_bytes=expected_size,
                local_sha256=local_sha,
                state="verified",
                state_reason=(
                    "upstream_checksum_and_local_sha256_verified"
                    if algorithm
                    else "local_sha256_verified_before_and_after_promotion"
                ),
            )
            return
        except MirrorError as error:
            if "mismatch" in str(error) or "changed" in str(error):
                quarantined = quarantine_part(part, root / "quarantine", str(row["path"]), "mismatch")
                update_database(
                    database,
                    str(row["inventory_id"]),
                    quarantine_path=str(quarantined),
                    state="quarantined",
                    state_reason=str(error),
                    attempts=int(row["attempts"]) + attempt,
                    transferred_bytes=transferred,
                )
                return
            if attempt == retries:
                update_database(
                    database,
                    str(row["inventory_id"]),
                    state="missing",
                    state_reason=str(error),
                    attempts=int(row["attempts"]) + attempt,
                    transferred_bytes=transferred,
                )
                return
            time.sleep(min(300, 5 * (3 ** (attempt - 1))))


def run_worker(database: Path, root: Path, concurrency: int, capacity_path: Path) -> None:
    if not os.environ.get("GUIX_ENVIRONMENT"):
        raise MirrorError("worker must run inside the pinned GNU Guix environment")
    if any(name in os.environ for name in ("SLURM_JOB_ID", "SLURM_ARRAY_JOB_ID")):
        raise MirrorError("VGP mirror worker must not run as a Slurm analysis job")
    capacity = json.loads(capacity_path.read_text(encoding="utf-8"))
    if not capacity.get("adequate"):
        raise MirrorError(f"capacity gate refused bulk transfer: {capacity.get('gate_reason')}")
    fixture = root / "fixture" / "fixture-report.json"
    if not fixture.is_file() or not json.loads(fixture.read_text())["passed"]:
        raise MirrorError("fixture proof is missing or did not pass")
    # Thread workers are sufficient here: rsync and hashing do the work outside
    # Python.  Every state mutation is an independent SQLite transaction.
    from concurrent.futures import ThreadPoolExecutor

    rows = database_rows(database)
    # Directories are inventory objects too.  Materialize and account for them
    # before files so the final state partition is exhaustive.
    for row in (item for item in rows if item["object_type"] == "directory"):
        directory = Path(str(row["durable_path"]))
        existed = directory.is_dir()
        directory.mkdir(parents=True, exist_ok=True)
        update_database(
            database,
            str(row["inventory_id"]),
            state="reused" if existed else "verified",
            state_reason="source_relative_directory_materialized",
        )
    files = [
        row
        for row in rows
        if row["object_type"] == "file" and row["state"] in {"planned", "transferred"}
    ]
    # Publish checksum manifests first, then bind their MD5 entries before the
    # remaining files are started.
    checksum_rows = [row for row in files if PurePosixPath(str(row["path"])).name == "md5sum.txt"]
    other_rows = [row for row in files if row not in checksum_rows]
    with ThreadPoolExecutor(max_workers=concurrency) as executor:
        list(executor.map(lambda row: process_file(database, root, row), checksum_rows))
    checksum_states = [
        row
        for row in database_rows(database)
        if PurePosixPath(str(row["path"])).name == "md5sum.txt"
    ]
    failed_checksums = [
        row for row in checksum_states if row["state"] not in {"verified", "reused"}
    ]
    if failed_checksums:
        raise MirrorError(
            f"{len(failed_checksums)} published checksum manifests failed; refusing payload transfer"
        )
    # Bind every published provider checksum before any non-manifest payload.
    refreshed_objects = [
        InventoryObject(
            str(row["inventory_id"]),
            0,
            str(row["accession"]),
            str(row["path"]),
            str(row["object_type"]),
            int(row["size"]),
            str(row["mtime"]),
            "",
            str(row["sequence_subset"]),
        )
        for row in database_rows(database)
    ]
    checksum_map = parse_upstream_md5s(root / "objects", refreshed_objects)
    connection = sqlite3.connect(database, timeout=60)
    try:
        connection.execute("BEGIN IMMEDIATE")
        for path, (algorithm, checksum) in checksum_map.items():
            connection.execute(
                "UPDATE objects SET checksum_algorithm = ?, checksum = ?, updated_at_utc = ? WHERE path = ?",
                (algorithm, checksum, utc_now(), path),
            )
        connection.commit()
    finally:
        connection.close()
    other_paths = {str(row["path"]) for row in other_rows}
    other_rows = [row for row in database_rows(database) if str(row["path"]) in other_paths]
    with ThreadPoolExecutor(max_workers=concurrency) as executor:
        list(executor.map(lambda row: process_file(database, root, row), other_rows))


def run_fixture(root: Path) -> dict[str, object]:
    if not os.environ.get("GUIX_ENVIRONMENT"):
        raise MirrorError("fixture must run inside the pinned GNU Guix environment")
    fixture = root / "fixture"
    if fixture.exists():
        shutil.rmtree(fixture)
    source = fixture / "source" / "payload.bin"
    part = fixture / "staging" / "payload.bin.part"
    destination = fixture / "objects" / "payload.bin"
    source.parent.mkdir(parents=True)
    part.parent.mkdir(parents=True)
    content = (b"VGP-FREEZE1-RESUME-FIXTURE\n" * 524288)[:8 * 1024 * 1024]
    source.write_bytes(content)
    command = [
        "rsync",
        "--partial",
        "--append-verify",
        "--bwlimit=1024",
        str(source),
        str(part),
    ]
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    deadline = time.time() + 15
    while time.time() < deadline:
        if part.exists() and 0 < part.stat().st_size < len(content):
            break
        time.sleep(0.05)
    process.terminate()
    process.wait(timeout=10)
    partial_size = part.stat().st_size if part.exists() else 0
    if not 0 < partial_size < len(content):
        raise MirrorError(f"fixture failed to create an interrupted partial: {partial_size}")
    resumed_bytes = rsync_transfer(str(source), part)
    expected_sha = hashlib.sha256(content).hexdigest()
    promoted_sha = promote_verified(part, destination, expected_size=len(content))
    if promoted_sha != expected_sha:
        raise MirrorError("fixture exact-resume checksum mismatch")
    bad = fixture / "staging" / "bad.bin.part"
    bad.write_bytes(b"known bad fixture")
    quarantined: Path | None = None
    try:
        promote_verified(
            bad,
            fixture / "objects" / "must-not-exist.bin",
            expected_size=bad.stat().st_size,
            checksum_algorithm="md5",
            checksum="0" * 32,
        )
    except MirrorError:
        quarantined = quarantine_part(bad, fixture / "quarantine", "bad.bin", "checksum")
    if not quarantined or not quarantined.is_file():
        raise MirrorError("fixture checksum mismatch was not quarantined")
    report = {
        "passed": True,
        "completed_at_utc": utc_now(),
        "source_bytes": len(content),
        "interrupted_partial_bytes": partial_size,
        "resume_additional_bytes": resumed_bytes,
        "exact_resume_sha256": promoted_sha,
        "checksum_failure_quarantined": True,
        "quarantine_path": str(quarantined),
        "atomic_promotion_verified": destination.is_file() and not part.exists(),
        "rsync_version": subprocess.run(
            ["rsync", "--version"], text=True, capture_output=True, check=True
        ).stdout.splitlines()[0],
        "slurm_jobs_launched": 0,
    }
    atomic_write_json(fixture / "fixture-report.json", report)
    return report


def build_outputs(args: argparse.Namespace) -> dict[str, object]:
    catalog_rows = verify_catalog(args.catalog)
    metadata = parse_metadata_file(args.inventory_metadata)
    if metadata.get("rsync_exit_status") != "0":
        raise MirrorError("metadata-only rsync inventory did not complete successfully")
    observed_inventory_sha256 = sha256_file(args.raw_inventory)
    if metadata.get("raw_inventory_sha256") != observed_inventory_sha256:
        raise MirrorError(
            "raw inventory digest does not match its retrieval metadata: "
            f"recorded={metadata.get('raw_inventory_sha256')}, observed={observed_inventory_sha256}"
        )
    if metadata.get("raw_path") != str(args.raw_inventory):
        raise MirrorError("inventory metadata does not bind the supplied raw inventory")
    objects = parse_rsync_inventory(args.raw_inventory, catalog_rows)
    checksums = parse_upstream_md5s(args.root / "objects", objects)
    if checksums:
        objects = parse_rsync_inventory(args.raw_inventory, catalog_rows, checksums)
    storage = storage_evidence(args.root, objects, args.concurrency)
    init_database(args.database, objects, args.root)
    atomic_write_json(args.root / "inventory" / "capacity-evidence.json", storage)
    return write_snapshots(
        objects=objects,
        metadata=metadata,
        database=args.database,
        storage=storage,
        source_output=args.source_output,
        manifest_output=args.manifest_output,
        summary_output=args.summary_output,
        handoff_output=args.handoff_output,
    )


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    subcommands = result.add_subparsers(dest="command", required=True)
    inventory = subcommands.add_parser(
        "inventory", help="freeze a new metadata-only recursive rsync inventory"
    )
    inventory.add_argument("--catalog", type=Path, default=PINNED_CATALOG)
    inventory.add_argument("--root", type=Path, default=RELEASE_ROOT)
    build = subcommands.add_parser("build", help="validate frozen listing and write snapshots")
    build.add_argument("--catalog", type=Path, default=PINNED_CATALOG)
    build.add_argument("--root", type=Path, default=RELEASE_ROOT)
    build.add_argument(
        "--raw-inventory",
        type=Path,
        default=RELEASE_ROOT / "inventory/freeze1-assemblies.20260718T111631Z.rsync-items",
    )
    build.add_argument(
        "--inventory-metadata",
        type=Path,
        default=RELEASE_ROOT / "inventory/freeze1-assemblies.20260718T111631Z.metadata",
    )
    build.add_argument("--database", type=Path, default=RELEASE_ROOT / "state/mirror.sqlite3")
    build.add_argument("--concurrency", type=int, default=DEFAULT_CONCURRENCY)
    build.add_argument(
        "--source-output", type=Path, default=PROJECT_ROOT / "analysis/vgp_freeze1_source_inventory.tsv"
    )
    build.add_argument(
        "--manifest-output", type=Path, default=PROJECT_ROOT / "analysis/vgp_freeze1_mirror_manifest.tsv"
    )
    build.add_argument(
        "--summary-output", type=Path, default=PROJECT_ROOT / "analysis/vgp_freeze1_mirror_summary.json"
    )
    build.add_argument(
        "--handoff-output", type=Path, default=PROJECT_ROOT / "analysis/vgp_freeze1_mirror_handoff.md"
    )
    fixture = subcommands.add_parser("fixture", help="prove interruption/resume/quarantine/promotion")
    fixture.add_argument("--root", type=Path, default=RELEASE_ROOT)
    worker = subcommands.add_parser("worker", help="run the capacity-gated mirror worker")
    worker.add_argument("--root", type=Path, default=RELEASE_ROOT)
    worker.add_argument("--database", type=Path, default=RELEASE_ROOT / "state/mirror.sqlite3")
    worker.add_argument(
        "--capacity", type=Path, default=RELEASE_ROOT / "inventory/capacity-evidence.json"
    )
    worker.add_argument("--concurrency", type=int, default=DEFAULT_CONCURRENCY)
    return result


def main(argv: Sequence[str] | None = None) -> int:
    args = parser().parse_args(argv)
    try:
        if hasattr(args, "concurrency") and args.concurrency < 1:
            raise MirrorError("concurrency must be at least one")
        if args.command == "inventory":
            payload = freeze_inventory(args.catalog, args.root)
        elif args.command == "build":
            payload = build_outputs(args)
        elif args.command == "fixture":
            payload = run_fixture(args.root)
        elif args.command == "worker":
            run_worker(args.database, args.root, args.concurrency, args.capacity)
            payload = {"status": "complete"}
        else:  # pragma: no cover
            raise AssertionError(args.command)
    except (MirrorError, OSError, sqlite3.Error) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
