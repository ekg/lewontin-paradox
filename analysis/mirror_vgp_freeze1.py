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
import urllib.parse
import urllib.request
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath
from typing import Iterable, Mapping, Sequence


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATA_ROOT_CONFIG = PROJECT_ROOT / "analysis/vgp_data_root_config.json"
EXCEPTION_LEDGER_OUTPUT = PROJECT_ROOT / "analysis/vgp_freeze1_exception_ledger.json"


def configured_vgp_root(config: Path = DATA_ROOT_CONFIG) -> Path:
    """Resolve every active VGP path through the repository root contract."""
    payload = json.loads(config.read_text(encoding="utf-8"))
    root = Path(str(payload["root"]))
    if not root.is_absolute():
        raise RuntimeError(f"configured VGP root is not absolute: {root}")
    return root


VGP_DATA_ROOT = configured_vgp_root()
RELEASE_ROOT = VGP_DATA_ROOT / "freeze1"
PINNED_CATALOG = VGP_DATA_ROOT / "manifests/VGPPhase1-freeze-1.0.commit-dc1b2af5a7741b97d66fb10cb2bce97f41765cdf.tsv"
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
    "verified_upstream_conflict",
}
DEFAULT_HEADROOM_FRACTION = 0.20
DEFAULT_CONCURRENCY = 2
CHECKSUM_RESERVE_BYTES = 64 * 1024 * 1024

SOURCE_FIELDS = [
    "canonical_vgp_root",
    "mirror_root",
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
    "canonical_vgp_root",
    "mirror_root",
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
    "cas_path",
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


class TransferError(MirrorError):
    """A resumable transfer attempt failed after possibly adding payload bytes."""

    def __init__(self, message: str, additional_bytes: int):
        super().__init__(message)
        self.additional_bytes = additional_bytes


def validate_verified_upstream_conflict(evidence: Mapping[str, object]) -> None:
    """Fail closed unless evidence satisfies the metadata-only exception policy."""
    if evidence.get("canonical_vgp_root") != str(VGP_DATA_ROOT):
        raise MirrorError("upstream-conflict evidence does not name the canonical VGP root")
    if evidence.get("sequence_subset") != "non_sequence_product_or_metadata":
        raise MirrorError("terminal exception is restricted to non-sequence metadata")
    if evidence.get("resolution") != "VERIFIED_UPSTREAM_CONFLICT":
        raise MirrorError("upstream-conflict evidence has no terminal resolution")
    frozen = evidence.get("frozen_catalog")
    if not isinstance(frozen, Mapping) or frozen.get("algorithm") != "md5":
        raise MirrorError("upstream-conflict evidence lacks the frozen MD5 binding")
    attempts = evidence.get("official_source_attempts")
    if not isinstance(attempts, list) or len(attempts) < 2:
        raise MirrorError("upstream conflict requires two independent official-source attempts")
    required_attempt_fields = {
        "source_url", "retrieval_started_utc", "retrieval_completed_utc",
        "size_bytes", "md5", "sha256", "quarantine_path",
    }
    for attempt in attempts:
        if not isinstance(attempt, Mapping) or required_attempt_fields - set(attempt):
            raise MirrorError("official-source conflict attempt is incomplete")
    observed = {
        (int(attempt["size_bytes"]), str(attempt["md5"]), str(attempt["sha256"]))
        for attempt in attempts
    }
    if len(observed) != 1:
        raise MirrorError("official-source conflict attempts are not reproducible")
    observed_size, observed_md5, observed_sha256 = next(iter(observed))
    if observed_size != int(frozen.get("size_bytes", -1)):
        raise MirrorError("official-source conflict changed the frozen object size")
    if observed_md5 == str(frozen.get("digest", "")):
        raise MirrorError("official-source bytes satisfy the frozen digest; no conflict exists")
    if not re.fullmatch(r"[0-9a-f]{32}", observed_md5) or not re.fullmatch(
        r"[0-9a-f]{64}", observed_sha256
    ):
        raise MirrorError("upstream-conflict evidence contains malformed digests")
    alternate = evidence.get("authoritative_alternate")
    if not isinstance(alternate, Mapping):
        raise MirrorError("upstream conflict lacks an authoritative alternate resolution")
    required_alternate_fields = {
        "assembly_accession", "source_url", "checksum_catalog_url",
        "retrieval_started_utc", "retrieval_completed_utc",
    }
    if required_alternate_fields - set(alternate):
        raise MirrorError("authoritative alternate resolution is incomplete")
    if urllib.parse.urlparse(str(alternate["source_url"])).hostname != "ftp.ncbi.nlm.nih.gov":
        raise MirrorError("authoritative alternate source is not official NCBI FTP")
    if urllib.parse.urlparse(str(alternate["checksum_catalog_url"])).hostname != "ftp.ncbi.nlm.nih.gov":
        raise MirrorError("authoritative alternate checksum catalog is not official NCBI FTP")
    if alternate.get("object_available") is False:
        if alternate.get("catalog_entry_found") is not False:
            raise MirrorError("NCBI no-equivalent resolution did not prove catalog-entry absence")
        retrieval = alternate.get("checksum_catalog_retrieval")
        if not isinstance(retrieval, Mapping) or not re.fullmatch(
            r"[0-9a-f]{64}", str(retrieval.get("sha256", ""))
        ):
            raise MirrorError("NCBI no-equivalent resolution lacks the retrieved catalog digest")
        if not alternate.get("searched_equivalent_names"):
            raise MirrorError("NCBI no-equivalent resolution records no searched names")
        return
    object_fields = {
        "size_bytes", "md5", "sha256", "catalog_algorithm", "catalog_digest",
        "quarantine_path",
    }
    if object_fields - set(alternate):
        raise MirrorError("authoritative NCBI object resolution is incomplete")
    alternate_md5 = str(alternate["md5"])
    alternate_sha256 = str(alternate["sha256"])
    if int(alternate["size_bytes"]) <= 0 or not re.fullmatch(r"[0-9a-f]{32}", alternate_md5):
        raise MirrorError("authoritative NCBI alternate has malformed size or MD5 evidence")
    if not re.fullmatch(r"[0-9a-f]{64}", alternate_sha256):
        raise MirrorError("authoritative NCBI alternate has malformed SHA-256 evidence")
    if alternate["catalog_algorithm"] != "md5" or alternate["catalog_digest"] != alternate_md5:
        raise MirrorError("authoritative NCBI checksum catalog does not bind alternate bytes")


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
            f"canonical_vgp_root={VGP_DATA_ROOT}",
            f"mirror_root={root}",
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


def parse_upstream_md5s(
    objects_root: Path,
    objects: Sequence[InventoryObject],
    ignored_provider_orphans: list[str] | None = None,
) -> dict[str, tuple[str, str]]:
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
            provider_prefix = "/mirrordata/hubs/"
            provider_absolute = relative.startswith(provider_prefix)
            if provider_absolute:
                candidate = relative.removeprefix(provider_prefix)
            elif PurePosixPath(relative).is_absolute():
                raise MirrorError(
                    f"upstream checksum manifest uses an unrecognized absolute prefix: {relative}"
                )
            else:
                candidate = str(PurePosixPath(root) / relative)
            if candidate not in inventoried_paths:
                if provider_absolute and candidate.startswith(root + "/"):
                    if ignored_provider_orphans is not None:
                        ignored_provider_orphans.append(candidate)
                    continue
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
    probe = root / f".write-probe-{os.getpid()}"
    try:
        with probe.open("xb") as handle:
            handle.write(b"vgp-freeze1-write-probe\n")
            handle.flush()
            os.fsync(handle.fileno())
        writable = True
    except OSError:
        writable = False
    finally:
        probe.unlink(missing_ok=True)
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
    filesystem_adequate = available >= required and free_inodes >= required_inodes and writable
    return {
        "canonical_vgp_root": str(VGP_DATA_ROOT),
        "mirror_root": str(root),
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
        "write_probe_passed": writable,
        "filesystem_capacity_adequate": filesystem_adequate,
        "quota_visibility_adequate": quota_visible,
        "quota_visibility_is_policy_gate": False,
        "adequate": filesystem_adequate,
        "gate_reason": (
            "capacity_write_and_inode_headroom_verified"
            if filesystem_adequate and quota_visible
            else "capacity_write_and_inode_headroom_verified_quota_helper_unavailable"
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
        connection.execute(
            """CREATE TABLE IF NOT EXISTS worker_runs (
                run_id INTEGER PRIMARY KEY AUTOINCREMENT, started_at_utc TEXT NOT NULL,
                started_epoch REAL NOT NULL, starting_transferred_bytes INTEGER NOT NULL,
                completed_at_utc TEXT NOT NULL DEFAULT '', outcome TEXT NOT NULL DEFAULT 'running',
                detail TEXT NOT NULL DEFAULT '', completed_transferred_bytes INTEGER NOT NULL DEFAULT -1
            )"""
        )
        run_columns = {
            row[1] for row in connection.execute("PRAGMA table_info(worker_runs)")
        }
        if "completed_transferred_bytes" not in run_columns:
            connection.execute(
                "ALTER TABLE worker_runs ADD COLUMN completed_transferred_bytes INTEGER NOT NULL DEFAULT -1"
            )
        connection.execute(
            """CREATE TABLE IF NOT EXISTS quarantine_events (
                event_id INTEGER PRIMARY KEY AUTOINCREMENT, inventory_id TEXT NOT NULL,
                quarantine_path TEXT NOT NULL UNIQUE, observed_bytes INTEGER NOT NULL,
                reason TEXT NOT NULL, created_at_utc TEXT NOT NULL
            )"""
        )
        connection.execute(
            """CREATE TABLE IF NOT EXISTS upstream_conflicts (
                inventory_id TEXT PRIMARY KEY, resolution TEXT NOT NULL,
                evidence_json TEXT NOT NULL, ledger_path TEXT NOT NULL,
                created_at_utc TEXT NOT NULL
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
        quarantine_root = root / "quarantine"
        if quarantine_root.is_dir():
            suffix = re.compile(r"\.\d{8}T\d{6}(?:\d{6})?Z\.([^.]+)$")
            for quarantined in quarantine_root.rglob("*"):
                if not quarantined.is_file():
                    continue
                relative = str(quarantined.relative_to(quarantine_root))
                match = suffix.search(relative)
                if not match:
                    continue
                source_path = relative[: match.start()]
                record = connection.execute(
                    "SELECT inventory_id FROM objects WHERE path = ?", (source_path,)
                ).fetchone()
                if record:
                    connection.execute(
                        """INSERT OR IGNORE INTO quarantine_events
                           (inventory_id, quarantine_path, observed_bytes, reason, created_at_utc)
                           VALUES (?, ?, ?, ?, ?)""",
                        (
                            record[0], str(quarantined), quarantined.stat().st_size,
                            match.group(1), utc_now(),
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
            "canonical_vgp_root": str(VGP_DATA_ROOT),
            "mirror_root": str(RELEASE_ROOT),
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
        digest = str(row["local_sha256"])
        yield {
            "canonical_vgp_root": str(VGP_DATA_ROOT),
            "mirror_root": str(RELEASE_ROOT),
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
            "cas_path": (
                str(VGP_DATA_ROOT / "objects/sha256" / digest[:2] / digest[2:4] / digest)
                if digest
                else ""
            ),
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


def mirror_progress(database: Path, root: Path) -> dict[str, object]:
    """Return exact restart-safe accounting without traversing the object tree."""
    rows = database_rows(database)
    files = [row for row in rows if row["object_type"] == "file"]
    total_bytes = sum(int(row["size"]) for row in files)
    verified = [row for row in files if row["state"] in {"verified", "reused"}]
    quarantined = [row for row in files if row["state"] == "quarantined"]
    exceptions = [row for row in files if row["state"] == "verified_upstream_conflict"]
    transferred_bytes = sum(int(row["transferred_bytes"]) for row in files)
    attempts = sum(int(row["attempts"]) for row in files)
    attempted_files = sum(int(row["attempts"]) > 0 for row in files)
    retries = max(0, attempts - attempted_files)
    connection = sqlite3.connect(database)
    connection.row_factory = sqlite3.Row
    try:
        quarantine_totals = connection.execute(
            "SELECT COUNT(*), COALESCE(SUM(observed_bytes), 0) FROM quarantine_events"
        ).fetchone()
        runs = connection.execute(
            "SELECT * FROM worker_runs ORDER BY run_id DESC LIMIT 5"
        ).fetchall()
    finally:
        connection.close()
    run_history: list[dict[str, object]] = []
    for run in runs:
        ended_epoch = time.time()
        if run["completed_at_utc"]:
            ended_epoch = datetime.fromisoformat(
                str(run["completed_at_utc"]).replace("Z", "+00:00")
            ).timestamp()
        elapsed = max(0.0, ended_epoch - float(run["started_epoch"]))
        ending_transferred = int(run["completed_transferred_bytes"])
        if ending_transferred < 0:
            next_run = min(
                (newer for newer in runs if int(newer["run_id"]) > int(run["run_id"])),
                key=lambda newer: int(newer["run_id"]),
                default=None,
            )
            ending_transferred = (
                int(next_run["starting_transferred_bytes"])
                if next_run is not None
                else transferred_bytes
            )
        added = max(0, ending_transferred - int(run["starting_transferred_bytes"]))
        run_history.append({
            "run_id": int(run["run_id"]),
            "started_at_utc": run["started_at_utc"],
            "completed_at_utc": run["completed_at_utc"],
            "outcome": run["outcome"],
            "detail": run["detail"],
            "elapsed_seconds": elapsed,
            "transferred_bytes": added,
            "throughput_bytes_per_second": (added / elapsed if elapsed else 0.0),
        })
    return {
        "schema_version": "vgp-freeze1-live-progress-v1",
        "generated_at_utc": utc_now(),
        "canonical_vgp_root": str(VGP_DATA_ROOT),
        "mirror_root": str(root),
        "inventory": {
            "objects": len(rows),
            "files": len(files),
            "bytes": total_bytes,
        },
        "progress": {
            "verified_files": len(verified),
            "verified_bytes": sum(int(row["size"]) for row in verified),
            "quarantine_events": int(quarantine_totals[0]),
            "quarantined_bytes": int(quarantine_totals[1]),
            "currently_quarantined_files": len(quarantined),
            "currently_quarantined_logical_bytes": sum(
                int(row["observed_bytes"]) for row in quarantined
            ),
            "verified_upstream_conflict_files": len(exceptions),
            "verified_upstream_conflict_logical_bytes": sum(
                int(row["size"]) for row in exceptions
            ),
            "transferred_network_bytes": transferred_bytes,
            "attempts": attempts,
            "retries": retries,
            "remaining_files": len(files) - len(verified) - len(exceptions),
            "remaining_bytes": total_bytes
            - sum(int(row["size"]) for row in verified)
            - sum(int(row["size"]) for row in exceptions),
            "state_accounting": state_accounting(rows),
        },
        "latest_worker_run": run_history[0] if run_history else None,
        "worker_run_history": run_history,
    }


def write_progress(database: Path, root: Path, output: Path | None = None) -> dict[str, object]:
    # Reconcile a SIGKILL/host-loss run only when no live worker owns the lock.
    import fcntl
    lock_path = root / "state/worker.lock"
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    with lock_path.open("a+", encoding="utf-8") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            pass
        else:
            rows = database_rows(database)
            transferred = sum(
                int(row["transferred_bytes"])
                for row in rows if row["object_type"] == "file"
            )
            connection = sqlite3.connect(database, timeout=60)
            try:
                connection.execute(
                    """UPDATE worker_runs SET completed_at_utc = ?, outcome = 'interrupted',
                       detail = 'reconciled_after_process_exit', completed_transferred_bytes = ?
                       WHERE outcome = 'running'""",
                    (utc_now(), transferred),
                )
                connection.commit()
            finally:
                connection.close()
    payload = mirror_progress(database, root)
    atomic_write_json(output or root / "state/progress.json", payload)
    return payload


def render_handoff(summary: Mapping[str, object]) -> str:
    totals = summary["inventory_totals"]
    storage = summary["storage"]
    reconciliation = summary["catalog_reconciliation"]
    state = summary["state_accounting"]
    req = storage["requirements"]
    statvfs = storage["statvfs"]
    progress = summary["live_progress"]
    blocker = (
        "The mirror is incomplete; transfer continues around a reproducible official-source "
        f"checksum mismatch: {progress['currently_quarantined_files']} current file / "
        f"{progress['currently_quarantined_logical_bytes']:,} logical bytes, across "
        f"{progress['quarantine_events']} quarantine events / "
        f"{progress['quarantined_bytes']:,} physical quarantined bytes.  See "
        "`analysis/vgp_freeze1_exception_ledger.json` for any completed adjudications."
        if progress["currently_quarantined_files"]
        else (
            f"The exception ledger contains {progress['verified_upstream_conflict_files']} "
            "exhaustively reproduced non-sequence VERIFIED_UPSTREAM_CONFLICT object(s); "
            "they remain quarantined and do not block unrelated or completed transfer."
            if progress["verified_upstream_conflict_files"]
            else "No current checksum, capacity, permission, or network blocker is recorded."
        )
    )
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
remaining payload if any checksum manifest fails, escapes its exact accession,
or uses an unrecognized prefix.  It records and ignores
{summary['checksum_policy']['ignored_provider_absolute_orphan_entries']} stale
provider-absolute checksum entries that remain in official catalogs while the
named file is absent from both the frozen inventory and source.  Such entries
do not expand the closed world.  Objects without a published checksum receive
the repeated local SHA-256 contract described below.

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

Capacity gate: **{storage['gate_reason']}**.  Direct filesystem byte/inode
headroom and an fsync-backed write probe are the operational authorization
evidence.  User-visible quota helpers are observability only: an unavailable
helper neither implies a quota nor blocks an explicitly authorized transfer.
Real write, ENOSPC, inode, network, and checksum failures remain hard errors.

## Live transfer checkpoint

Verified: {progress['verified_files']:,} files / {progress['verified_bytes']:,}
bytes.  Network payload: {progress['transferred_network_bytes']:,} bytes across
{progress['attempts']:,} attempts and {progress['retries']:,} retries.  Remaining:
{progress['remaining_files']:,} files / {progress['remaining_bytes']:,} bytes.
{blocker}

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
analysis/run_vgp_freeze1_mirror.sh status
```

`inventory` always creates a new timestamped snapshot and never replaces the
bound snapshot.  To reconcile a deliberate refresh, pass the new raw and
metadata paths explicitly to `build` and review the closed-world result before
allowing the worker to see its capacity evidence.

The worker uses source-relative `.part` staging, `rsync --partial
--append-verify`, bounded concurrency, exponential backoff, size plus published
MD5 verification when available, and repeated local SHA-256 validation.  It
atomically inserts the object at the shared digest-derived CAS path and then
atomically hard-links the source-relative mirror view.  Already verified CAS
objects are revalidated and reused by exact size plus provider MD5 without a
redownload.  A mismatch is moved to source-relative quarantine and never
overwrites or deletes the last verified object.  There is no mirror-wide delete
operation.  Each object becomes available immediately after its own verification.

`status` atomically refreshes `state/progress.json` with exact inventory,
verified/quarantined/network bytes, attempts, retries, per-state accounting,
elapsed time, and run throughput.  A stopped worker leaves `.part` files and
SQLite transactions durable; invoking `worker` again continues with
`--append-verify` and revalidates any view published before an interruption.

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
    exception_ledger_output: Path,
    ignored_provider_checksum_orphans: Sequence[str] = (),
) -> dict[str, object]:
    rows = database_rows(database)
    totals = inventory_totals(objects)
    accounting = state_accounting(rows)
    released = len({obj.accession for obj in objects})
    checksum_manifests = [
        obj for obj in objects if obj.object_type == "file" and PurePosixPath(obj.path).name == "md5sum.txt"
    ]
    published_checksum_bindings = sum(bool(obj.checksum_algorithm) for obj in objects)
    progress = mirror_progress(database, Path(str(storage["path"])))
    exceptions = conflict_entries(database)
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
        "canonical_vgp_root": str(VGP_DATA_ROOT),
        "mirror_root": str(RELEASE_ROOT),
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
            "ignored_provider_absolute_orphan_entries": len(ignored_provider_checksum_orphans),
            "ignored_provider_absolute_orphan_paths": sorted(ignored_provider_checksum_orphans),
            "isolated_conflicts_continue_unrelated_transfer": True,
            "terminal_non_sequence_exception": "VERIFIED_UPSTREAM_CONFLICT",
            "sequence_conflicts_require_scientific_review": True,
        },
        "checksum_exceptions": exceptions,
        "storage": storage,
        "state_accounting": accounting,
        "operational_errors": [
            {
                "source_relative_path": row["path"],
                "state": row["state"],
                "reason": row["state_reason"],
                "observed_bytes": row["observed_bytes"],
                "quarantine_path": row["quarantine_path"],
            }
            for row in rows if row["state"] in {"missing", "quarantined", "superseded"}
        ],
        "live_progress": progress["progress"],
        "worker_run_history": progress["worker_run_history"],
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
    write_exception_ledger(database, Path(str(storage["path"])), exception_ledger_output)
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


def record_quarantine_event(
    database: Path, inventory_id: str, path: Path, reason: str
) -> None:
    connection = sqlite3.connect(database, timeout=60)
    try:
        connection.execute(
            """INSERT OR IGNORE INTO quarantine_events
               (inventory_id, quarantine_path, observed_bytes, reason, created_at_utc)
               VALUES (?, ?, ?, ?, ?)""",
            (inventory_id, str(path), path.stat().st_size, reason, utc_now()),
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
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    destination = quarantine_root / f"{relative}.{stamp}.{reason}"
    destination.parent.mkdir(parents=True, exist_ok=True)
    if part.exists():
        os.replace(part, destination)
    return destination


def file_digests(path: Path) -> dict[str, object]:
    md5 = hashlib.md5()
    sha256 = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            md5.update(block)
            sha256.update(block)
    return {"size_bytes": path.stat().st_size, "md5": md5.hexdigest(), "sha256": sha256.hexdigest()}


def fetch_https(url: str, destination: Path) -> dict[str, object]:
    parsed = urllib.parse.urlparse(url)
    if parsed.scheme != "https":
        raise MirrorError(f"authoritative alternate URL must use HTTPS: {url}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.unlink(missing_ok=True)
    started = utc_now()
    request = urllib.request.Request(url, headers={"User-Agent": "VGP-Freeze1-mirror/1.0"})
    try:
        with urllib.request.urlopen(request, timeout=300) as response, destination.open("xb") as out:
            shutil.copyfileobj(response, out, length=8 * 1024 * 1024)
            out.flush()
            os.fsync(out.fileno())
    except Exception as error:
        destination.unlink(missing_ok=True)
        raise MirrorError(f"authoritative alternate retrieval failed for {url}: {error}") from error
    completed = utc_now()
    return {
        "source_url": url,
        "retrieval_started_utc": started,
        "retrieval_completed_utc": completed,
        **file_digests(destination),
    }


def ensure_conflict_table(database: Path) -> None:
    connection = sqlite3.connect(database, timeout=60)
    try:
        connection.execute(
            """CREATE TABLE IF NOT EXISTS upstream_conflicts (
                inventory_id TEXT PRIMARY KEY, resolution TEXT NOT NULL,
                evidence_json TEXT NOT NULL, ledger_path TEXT NOT NULL,
                created_at_utc TEXT NOT NULL
            )"""
        )
        connection.commit()
    finally:
        connection.close()


def conflict_entries(database: Path) -> list[dict[str, object]]:
    ensure_conflict_table(database)
    connection = sqlite3.connect(database)
    try:
        return [
            json.loads(row[0])
            for row in connection.execute(
                "SELECT evidence_json FROM upstream_conflicts ORDER BY inventory_id"
            )
        ]
    finally:
        connection.close()


def write_exception_ledger(database: Path, root: Path, output: Path) -> dict[str, object]:
    payload = {
        "schema_version": "vgp-freeze1-upstream-conflict-ledger-v1",
        "canonical_vgp_root": str(VGP_DATA_ROOT),
        "mirror_root": str(root),
        "generated_at_utc": utc_now(),
        "entries": conflict_entries(database),
    }
    atomic_write_json(root / "state/upstream-conflict-ledger.json", payload)
    atomic_write_json(output, payload)
    return payload


def record_verified_upstream_conflict(
    database: Path, root: Path, evidence: Mapping[str, object]
) -> None:
    validate_verified_upstream_conflict(evidence)
    inventory_id = str(evidence["inventory_id"])
    ledger_path = root / "state/upstream-conflicts" / f"{inventory_id}.json"
    atomic_write_json(ledger_path, evidence)
    ensure_conflict_table(database)
    connection = sqlite3.connect(database, timeout=60)
    try:
        connection.execute("BEGIN IMMEDIATE")
        connection.execute(
            """INSERT OR REPLACE INTO upstream_conflicts
               (inventory_id, resolution, evidence_json, ledger_path, created_at_utc)
               VALUES (?, 'VERIFIED_UPSTREAM_CONFLICT', ?, ?, ?)""",
            (inventory_id, json.dumps(evidence, sort_keys=True), str(ledger_path), utc_now()),
        )
        connection.execute(
            """UPDATE objects SET state = 'verified_upstream_conflict',
               state_reason = 'VERIFIED_UPSTREAM_CONFLICT: reproducible frozen UCSC catalog/source conflict; authoritative NCBI bytes and catalog agree; bytes remain quarantined',
               updated_at_utc = ? WHERE inventory_id = ?""",
            (utc_now(), inventory_id),
        )
        connection.commit()
    finally:
        connection.close()


def reproduce_official_source_conflict(
    database: Path, root: Path, row: Mapping[str, object]
) -> list[dict[str, object]]:
    """Retrieve and quarantine two independent current official-source copies."""
    source_url = TRANSPORT_ENDPOINT + str(row["path"])
    part = Path(str(row["staging_path"]))
    attempts: list[dict[str, object]] = []
    transferred = int(row["transferred_bytes"])
    attempt_count = int(row["attempts"])
    for number in (1, 2):
        part.unlink(missing_ok=True)
        started = utc_now()
        added = rsync_transfer(source_url, part)
        completed = utc_now()
        digests = file_digests(part)
        if int(digests["size_bytes"]) != int(row["size"]):
            raise MirrorError("official-source conflict reproduction changed frozen object size")
        if digests["md5"] == row["checksum"]:
            raise MirrorError("fresh official-source retrieval now satisfies frozen checksum")
        quarantined = quarantine_part(part, root / "quarantine", str(row["path"]), "mismatch")
        record_quarantine_event(database, str(row["inventory_id"]), quarantined, f"independent_reproduction_{number}")
        transferred += added
        attempt_count += 1
        update_database(
            database,
            str(row["inventory_id"]),
            observed_bytes=int(digests["size_bytes"]),
            quarantine_path=str(quarantined),
            state="quarantined",
            state_reason=f"independent upstream checksum conflict reproduction {number}",
            attempts=attempt_count,
            transferred_bytes=transferred,
        )
        attempts.append(
            {
                "attempt": number,
                "source_url": source_url,
                "retrieval_started_utc": started,
                "retrieval_completed_utc": completed,
                **digests,
                "advertised_algorithm": row["checksum_algorithm"],
                "advertised_digest": row["checksum"],
                "quarantine_path": str(quarantined),
            }
        )
    return attempts


def adjudicate_upstream_conflict(
    database: Path,
    root: Path,
    inventory_id: str,
    alternate_source: str,
    alternate_checksum_source: str,
    ledger_output: Path,
) -> dict[str, object]:
    """Reproduce a quarantined metadata mismatch twice and resolve it against NCBI."""
    if not os.environ.get("GUIX_ENVIRONMENT"):
        raise MirrorError("conflict adjudication must run inside the pinned GNU Guix environment")
    matching = [row for row in database_rows(database) if row["inventory_id"] == inventory_id]
    if len(matching) != 1:
        raise MirrorError(f"inventory object not found: {inventory_id}")
    row = matching[0]
    if row["object_type"] != "file" or row["sequence_subset"] != "non_sequence_product_or_metadata":
        raise MirrorError("terminal conflict adjudication is restricted to non-sequence metadata files")
    if row["checksum_algorithm"] != "md5" or not row["checksum"]:
        raise MirrorError("conflict adjudication requires a frozen provider MD5 binding")
    if urllib.parse.urlparse(alternate_source).hostname != "ftp.ncbi.nlm.nih.gov" or urllib.parse.urlparse(
        alternate_checksum_source
    ).hostname != "ftp.ncbi.nlm.nih.gov":
        raise MirrorError("alternate source and checksum catalog must be official NCBI FTP HTTPS URLs")

    source_url = TRANSPORT_ENDPOINT + str(row["path"])
    attempts = reproduce_official_source_conflict(database, root, row)

    audit_root = root / "state/conflict-resolution" / inventory_id
    alternate_part = audit_root / "alternate.part"
    alternate = fetch_https(alternate_source, alternate_part)
    catalog_part = audit_root / "checksum-catalog.part"
    catalog = fetch_https(alternate_checksum_source, catalog_part)
    target_name = PurePosixPath(urllib.parse.urlparse(alternate_source).path).name
    catalog_digest = ""
    for line in catalog_part.read_text(encoding="utf-8").splitlines():
        match = re.fullmatch(r"([0-9a-fA-F]{32})\s+[* ]?(?:\./)?(.+)", line.strip())
        if match and PurePosixPath(match.group(2)).name == target_name:
            catalog_digest = match.group(1).lower()
            break
    if not catalog_digest:
        raise MirrorError(f"authoritative NCBI checksum catalog has no entry for {target_name}")
    alternate_quarantine = quarantine_part(
        alternate_part, root / "quarantine", str(row["path"]), "ncbi-authoritative-alternate"
    )
    record_quarantine_event(
        database, inventory_id, alternate_quarantine, "authoritative_ncbi_alternate"
    )
    catalog_final = audit_root / "ncbi-md5checksums.txt"
    os.replace(catalog_part, catalog_final)
    accession_match = PATH_ACCESSION_RE.search(str(row["path"]))
    evidence: dict[str, object] = {
        "schema_version": "vgp-freeze1-upstream-conflict-v1",
        "canonical_vgp_root": str(VGP_DATA_ROOT),
        "mirror_root": str(root),
        "inventory_id": inventory_id,
        "source_relative_path": row["path"],
        "sequence_subset": row["sequence_subset"],
        "frozen_catalog": {
            "source_url": source_url,
            "algorithm": row["checksum_algorithm"],
            "digest": row["checksum"],
            "size_bytes": int(row["size"]),
            "source_mtime_utc": row["mtime"],
        },
        "official_source_attempts": attempts,
        "authoritative_alternate": {
            "assembly_accession": accession_match.group(1) if accession_match else row["accession"],
            **alternate,
            "checksum_catalog_url": alternate_checksum_source,
            "checksum_catalog_retrieval": {
                **catalog,
                "local_audit_path": str(catalog_final),
            },
            "catalog_algorithm": "md5",
            "catalog_digest": catalog_digest,
            "comparison_to_official_source": (
                "identical"
                if (
                    int(alternate["size_bytes"]), alternate["md5"], alternate["sha256"]
                ) == (
                    int(attempts[0]["size_bytes"]),
                    attempts[0]["md5"],
                    attempts[0]["sha256"],
                )
                else "different_authoritative_current_record"
            ),
            "quarantine_path": str(alternate_quarantine),
        },
        "promotion_permitted": False,
        "scientific_sequence_review_required": False,
        "resolution": "VERIFIED_UPSTREAM_CONFLICT",
        "resolved_at_utc": utc_now(),
    }
    record_verified_upstream_conflict(database, root, evidence)
    ledger = write_exception_ledger(database, root, ledger_output)
    return {"evidence": evidence, "exception_ledger": ledger}


def adjudicate_ncbi_catalog_absence(
    database: Path,
    root: Path,
    inventory_id: str,
    alternate_directory: str,
    alternate_checksum_source: str,
    ledger_output: Path,
) -> dict[str, object]:
    """Document that NCBI publishes no exact equivalent for UCSC-only metadata."""
    if not os.environ.get("GUIX_ENVIRONMENT"):
        raise MirrorError("conflict adjudication must run inside the pinned GNU Guix environment")
    matching = [row for row in database_rows(database) if row["inventory_id"] == inventory_id]
    if len(matching) != 1:
        raise MirrorError(f"inventory object not found: {inventory_id}")
    row = matching[0]
    if row["object_type"] != "file" or row["sequence_subset"] != "non_sequence_product_or_metadata":
        raise MirrorError("terminal conflict adjudication is restricted to non-sequence metadata files")
    if row["checksum_algorithm"] != "md5" or not row["checksum"]:
        raise MirrorError("conflict adjudication requires a frozen provider MD5 binding")
    if urllib.parse.urlparse(alternate_directory).hostname != "ftp.ncbi.nlm.nih.gov" or urllib.parse.urlparse(
        alternate_checksum_source
    ).hostname != "ftp.ncbi.nlm.nih.gov":
        raise MirrorError("alternate directory and checksum catalog must be official NCBI URLs")

    attempts = reproduce_official_source_conflict(database, root, row)
    audit_root = root / "state/conflict-resolution" / inventory_id
    catalog_part = audit_root / "checksum-catalog.part"
    catalog = fetch_https(alternate_checksum_source, catalog_part)
    catalog_names: list[str] = []
    for line in catalog_part.read_text(encoding="utf-8").splitlines():
        match = re.fullmatch(r"[0-9a-fA-F]{32}\s+[* ]?(?:\./)?(.+)", line.strip())
        if match:
            catalog_names.append(PurePosixPath(match.group(1)).name)
    searched = [PurePosixPath(str(row["path"])).name]
    found = sorted(set(searched) & set(catalog_names))
    if found:
        raise MirrorError(f"NCBI checksum catalog unexpectedly contains equivalent object(s): {found}")
    catalog_final = audit_root / "ncbi-md5checksums.txt"
    os.replace(catalog_part, catalog_final)
    evidence: dict[str, object] = {
        "schema_version": "vgp-freeze1-upstream-conflict-v1",
        "canonical_vgp_root": str(VGP_DATA_ROOT),
        "mirror_root": str(root),
        "inventory_id": inventory_id,
        "source_relative_path": row["path"],
        "sequence_subset": row["sequence_subset"],
        "frozen_catalog": {
            "source_url": TRANSPORT_ENDPOINT + str(row["path"]),
            "algorithm": row["checksum_algorithm"],
            "digest": row["checksum"],
            "size_bytes": int(row["size"]),
            "source_mtime_utc": row["mtime"],
        },
        "official_source_attempts": attempts,
        "authoritative_alternate": {
            "assembly_accession": row["accession"],
            "source_url": alternate_directory,
            "checksum_catalog_url": alternate_checksum_source,
            "retrieval_started_utc": catalog["retrieval_started_utc"],
            "retrieval_completed_utc": catalog["retrieval_completed_utc"],
            "checksum_catalog_retrieval": {
                **catalog,
                "local_audit_path": str(catalog_final),
            },
            "object_available": False,
            "catalog_entry_found": False,
            "searched_equivalent_names": searched,
            "catalog_entry_count": len(catalog_names),
            "resolution_kind": "ncbi_catalog_no_equivalent_ucsc_metadata_object",
        },
        "promotion_permitted": False,
        "scientific_sequence_review_required": False,
        "resolution": "VERIFIED_UPSTREAM_CONFLICT",
        "resolved_at_utc": utc_now(),
    }
    record_verified_upstream_conflict(database, root, evidence)
    ledger = write_exception_ledger(database, root, ledger_output)
    return {"evidence": evidence, "exception_ledger": ledger}


def ncbi_assembly_urls(accession: str, directory_index: str) -> tuple[str, str]:
    """Resolve one exact accession/version directory from an official NCBI index."""
    candidates = sorted(
        {
            urllib.parse.unquote(match)
            for match in re.findall(r'href="([^"/]+/)"', directory_index)
            if urllib.parse.unquote(match).startswith(accession + "_")
        }
    )
    if len(candidates) != 1:
        raise MirrorError(
            f"NCBI directory index resolves {accession} to {len(candidates)} exact-version entries"
        )
    directory = candidates[0].rstrip("/")
    accession_root = accession_path(accession)
    base = f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{accession_root.rsplit('/', 1)[0]}/{directory}"
    return f"{base}/{directory}_assembly_report.txt", f"{base}/md5checksums.txt"


def adjudicate_current_metadata_conflicts(
    database: Path, root: Path, ledger_output: Path
) -> dict[str, object]:
    """Adjudicate eligible assembly-report conflicts without stopping on one failure."""
    results: list[dict[str, object]] = []
    errors: list[dict[str, str]] = []
    eligible = [
        row
        for row in database_rows(database)
        if row["state"] == "quarantined"
        and row["sequence_subset"] == "non_sequence_product_or_metadata"
        and int(row["attempts"]) >= 2
    ]
    accession_cache: dict[str, tuple[str, str, dict[str, object]]] = {}
    for row in eligible:
        inventory_id = str(row["inventory_id"])
        accession = str(row["accession"])
        try:
            if accession not in accession_cache:
                accession_root = accession_path(accession).rsplit("/", 1)[0]
                index_url = f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{accession_root}/"
                audit_root = root / "state/ncbi-resolution" / accession
                index_part = audit_root / "directory-index.part"
                index_retrieval = fetch_https(index_url, index_part)
                source_url, checksum_url = ncbi_assembly_urls(
                    accession, index_part.read_text(encoding="utf-8")
                )
                index_final = audit_root / "directory-index.html"
                os.replace(index_part, index_final)
                accession_cache[accession] = (
                    source_url,
                    checksum_url,
                    {**index_retrieval, "local_audit_path": str(index_final)},
                )
            source_url, checksum_url, index_evidence = accession_cache[accession]
            if PurePosixPath(str(row["path"])).name.endswith("_assembly_report.txt"):
                result = adjudicate_upstream_conflict(
                    database, root, inventory_id, source_url, checksum_url, ledger_output
                )
            else:
                result = adjudicate_ncbi_catalog_absence(
                    database,
                    root,
                    inventory_id,
                    checksum_url.rsplit("/", 1)[0] + "/",
                    checksum_url,
                    ledger_output,
                )
            results.append(
                {
                    "inventory_id": inventory_id,
                    "accession": accession,
                    "ncbi_directory_index": index_evidence,
                    "resolution": result["evidence"]["resolution"],
                }
            )
        except MirrorError as error:
            errors.append(
                {"inventory_id": inventory_id, "accession": accession, "error": str(error)}
            )
    return {
        "canonical_vgp_root": str(VGP_DATA_ROOT),
        "mirror_root": str(root),
        "eligible_conflicts": len(eligible),
        "resolved_conflicts": len(results),
        "errors": errors,
        "results": results,
    }


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


def atomic_publish_view(cas_path: Path, destination: Path) -> None:
    """Atomically hard-link a verified CAS object into a source-relative view."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        if not destination.is_file() or not os.path.samefile(cas_path, destination):
            raise MirrorError(f"durable view conflicts with verified CAS object: {destination}")
        return
    temporary = destination.with_name(f".{destination.name}.publish-{os.getpid()}-{time.time_ns()}")
    try:
        os.link(cas_path, temporary)
        os.replace(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)
    directory_fd = os.open(destination.parent, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(directory_fd)
    finally:
        os.close(directory_fd)


def promote_to_cas(
    part: Path,
    destination: Path,
    *,
    expected_size: int,
    checksum_algorithm: str = "",
    checksum: str = "",
) -> tuple[str, Path]:
    """Verify staging, atomically create the shared CAS object, then publish its view."""
    if part.stat().st_size != expected_size:
        raise MirrorError(f"size mismatch: expected {expected_size}, observed {part.stat().st_size}")
    verify_digest(part, checksum_algorithm, checksum)
    first_sha = sha256_file(part)
    if sha256_file(part) != first_sha:
        raise MirrorError("staged object changed before CAS promotion")
    cas_path = VGP_DATA_ROOT / "objects/sha256" / first_sha[:2] / first_sha[2:4] / first_sha
    cas_path.parent.mkdir(parents=True, exist_ok=True)
    os.chmod(part, stat.S_IRUSR | stat.S_IRGRP)
    try:
        # link(2), unlike replace(2), is a no-clobber atomic CAS insertion.
        os.link(part, cas_path)
        created = True
    except FileExistsError:
        created = False
    if cas_path.stat().st_size != expected_size or sha256_file(cas_path) != first_sha:
        raise MirrorError(f"content-addressed destination conflicts: {cas_path}")
    part.unlink()
    atomic_publish_view(cas_path, destination)
    if sha256_file(destination) != first_sha:
        raise MirrorError("published source-relative view changed after CAS promotion")
    if created:
        directory_fd = os.open(cas_path.parent, os.O_RDONLY | os.O_DIRECTORY)
        try:
            os.fsync(directory_fd)
        finally:
            os.close(directory_fd)
    return first_sha, cas_path


def index_verified_cas(cas_root: Path) -> dict[tuple[int, str], tuple[Path, str]]:
    """Revalidate canonical CAS objects once and index them for provider-MD5 reuse."""
    result: dict[tuple[int, str], tuple[Path, str]] = {}
    if not cas_root.is_dir():
        return result
    for candidate in sorted(cas_root.glob("??/??/" + "?" * 64)):
        if not candidate.is_file() or not re.fullmatch(r"[0-9a-f]{64}", candidate.name):
            continue
        local_sha = sha256_file(candidate)
        if local_sha != candidate.name:
            raise MirrorError(f"canonical CAS path failed SHA-256 identity: {candidate}")
        digest = hashlib.md5()
        with candidate.open("rb") as handle:
            for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
                digest.update(block)
        key = (candidate.stat().st_size, digest.hexdigest())
        prior = result.get(key)
        if prior and prior[1] != local_sha:
            raise MirrorError(f"ambiguous size/MD5 collision in canonical CAS: {candidate}")
        result[key] = (candidate, local_sha)
    return result


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
    after = part.stat().st_size if part.exists() else 0
    if completed.returncode != 0:
        raise TransferError(
            f"rsync failed ({completed.returncode}) for {source}: {completed.stderr.strip()}",
            max(0, after - before),
        )
    return max(0, after - before)


def process_file(
    database: Path,
    root: Path,
    row: Mapping[str, object],
    cas_by_md5: Mapping[tuple[int, str], tuple[Path, str]] | None = None,
    retries: int = 5,
) -> None:
    destination = Path(str(row["durable_path"]))
    part = Path(str(row["staging_path"]))
    expected_size = int(row["size"])
    algorithm = str(row["checksum_algorithm"])
    checksum = str(row["checksum"])
    if destination.is_file():
        try:
            if destination.stat().st_size != expected_size:
                raise MirrorError("existing durable size mismatch")
            verify_digest(destination, algorithm, checksum)
            local_sha = sha256_file(destination)
            if row["local_sha256"] and local_sha != row["local_sha256"]:
                raise MirrorError("existing durable local SHA-256 mismatch")
            cas_path = VGP_DATA_ROOT / "objects/sha256" / local_sha[:2] / local_sha[2:4] / local_sha
            if not cas_path.is_file() or not os.path.samefile(cas_path, destination):
                raise MirrorError("pre-existing durable object is not a canonical CAS view")
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
    attempts_before_run = int(row["attempts"])
    if part.is_file() and part.stat().st_size > int(row["observed_bytes"]):
        # A process can die after rsync has flushed its partial but before the
        # parent records the failed attempt.  Recover that durable evidence so
        # network-byte and retry accounting remains exact across SIGKILL/host loss.
        transferred = max(transferred, part.stat().st_size)
        attempts_before_run += 1
        update_database(
            database,
            str(row["inventory_id"]),
            observed_bytes=part.stat().st_size,
            transferred_bytes=transferred,
            attempts=attempts_before_run,
            state="transferred",
            state_reason="recovered_uncommitted_partial_after_process_restart",
        )
    if algorithm == "md5" and cas_by_md5:
        reusable = cas_by_md5.get((expected_size, checksum))
        if reusable:
            cas_path, local_sha = reusable
            verify_digest(cas_path, algorithm, checksum)
            atomic_publish_view(cas_path, destination)
            update_database(
                database,
                str(row["inventory_id"]),
                observed_bytes=expected_size,
                local_sha256=local_sha,
                state="reused",
                state_reason="verified_canonical_cas_md5_reuse_without_redownload",
                attempts=attempts_before_run,
                transferred_bytes=transferred,
            )
            part.unlink(missing_ok=True)
            return
    for attempt in range(1, retries + 1):
        try:
            transferred += rsync_transfer(TRANSPORT_ENDPOINT + str(row["path"]), part)
            update_database(
                database,
                str(row["inventory_id"]),
                observed_bytes=part.stat().st_size,
                transferred_bytes=transferred,
                attempts=attempts_before_run + attempt,
                state="transferred",
                state_reason="complete_part_pending_verification",
            )
            local_sha, _cas_path = promote_to_cas(
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
            if isinstance(error, TransferError):
                transferred += error.additional_bytes
                update_database(
                    database,
                    str(row["inventory_id"]),
                    observed_bytes=part.stat().st_size if part.exists() else 0,
                    attempts=attempts_before_run + attempt,
                    transferred_bytes=transferred,
                    state="transferred",
                    state_reason=f"resumable_partial_after_error:{error}",
                )
            if "mismatch" in str(error) or "changed" in str(error):
                quarantined = quarantine_part(part, root / "quarantine", str(row["path"]), "mismatch")
                record_quarantine_event(
                    database, str(row["inventory_id"]), quarantined, str(error)
                )
                update_database(
                    database,
                    str(row["inventory_id"]),
                    quarantine_path=str(quarantined),
                    state="quarantined",
                    state_reason=str(error),
                    attempts=attempts_before_run + attempt,
                    transferred_bytes=transferred,
                )
                if attempt < min(retries, 2):
                    time.sleep(1)
                    continue
                return
            if attempt == retries:
                update_database(
                    database,
                    str(row["inventory_id"]),
                    state="missing",
                    state_reason=str(error),
                    attempts=attempts_before_run + attempt,
                    transferred_bytes=transferred,
                )
                return
            time.sleep(min(300, 5 * (3 ** (attempt - 1))))


def recheck_live_capacity(root: Path, capacity: Mapping[str, object]) -> dict[str, object]:
    requirements = capacity["requirements"]
    values = os.statvfs(root)
    available = values.f_bavail * values.f_frsize
    free_inodes = values.f_favail
    required_bytes = int(requirements["total_bytes"])
    required_inodes = int(requirements["total_inodes"])
    probe = root / f".worker-write-probe-{os.getpid()}"
    try:
        with probe.open("xb") as handle:
            handle.write(b"vgp-freeze1-worker-write-probe\n")
            handle.flush()
            os.fsync(handle.fileno())
    finally:
        probe.unlink(missing_ok=True)
    if available < required_bytes:
        raise MirrorError(
            f"live filesystem capacity inadequate: {available} available < {required_bytes} required"
        )
    if free_inodes < required_inodes:
        raise MirrorError(
            f"live filesystem inode capacity inadequate: {free_inodes} available < {required_inodes} required"
        )
    return {
        "checked_at_utc": utc_now(),
        "available_bytes": available,
        "required_bytes": required_bytes,
        "free_inodes": free_inodes,
        "required_inodes": required_inodes,
        "write_probe_passed": True,
    }


def run_worker(database: Path, root: Path, concurrency: int, capacity_path: Path) -> None:
    if not os.environ.get("GUIX_ENVIRONMENT"):
        raise MirrorError("worker must run inside the pinned GNU Guix environment")
    if any(name in os.environ for name in ("SLURM_JOB_ID", "SLURM_ARRAY_JOB_ID")):
        raise MirrorError("VGP mirror worker must not run as a Slurm analysis job")
    import fcntl
    lock_path = root / "state/worker.lock"
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    worker_lock = lock_path.open("a+", encoding="utf-8")
    try:
        fcntl.flock(worker_lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError as error:
        raise MirrorError(f"another VGP Freeze 1 worker holds {lock_path}") from error
    worker_lock.seek(0)
    worker_lock.truncate()
    worker_lock.write(f"pid={os.getpid()} started_at_utc={utc_now()}\n")
    worker_lock.flush()
    os.fsync(worker_lock.fileno())
    capacity = json.loads(capacity_path.read_text(encoding="utf-8"))
    if not capacity.get("filesystem_capacity_adequate") or not capacity.get("write_probe_passed"):
        raise MirrorError(f"capacity gate refused bulk transfer: {capacity.get('gate_reason')}")
    live_capacity = recheck_live_capacity(root, capacity)
    atomic_write_json(root / "state/live-capacity.json", {
        "canonical_vgp_root": str(VGP_DATA_ROOT),
        "mirror_root": str(root),
        **live_capacity,
    })
    fixture = root / "fixture" / "fixture-report.json"
    if not fixture.is_file() or not json.loads(fixture.read_text())["passed"]:
        raise MirrorError("fixture proof is missing or did not pass")
    # Thread workers are sufficient here: rsync and hashing do the work outside
    # Python.  Every state mutation is an independent SQLite transaction.
    from concurrent.futures import ThreadPoolExecutor

    rows = database_rows(database)
    starting_transferred = sum(
        int(row["transferred_bytes"]) for row in rows if row["object_type"] == "file"
    )
    connection = sqlite3.connect(database, timeout=60)
    try:
        connection.execute(
            """UPDATE worker_runs SET completed_at_utc = ?, outcome = 'interrupted',
               detail = 'superseded_by_restart_after_process_exit',
               completed_transferred_bytes = ?
               WHERE outcome = 'running'""",
            (utc_now(), starting_transferred),
        )
        cursor = connection.execute(
            "INSERT INTO worker_runs (started_at_utc, started_epoch, starting_transferred_bytes) VALUES (?, ?, ?)",
            (utc_now(), time.time(), starting_transferred),
        )
        run_id = int(cursor.lastrowid)
        connection.commit()
    finally:
        connection.close()
    write_progress(database, root)
    # Directories are inventory objects too.  Materialize and account for them
    # before files so the final state partition is exhaustive.
    for row in (
        item
        for item in rows
        if item["object_type"] == "directory"
        and item["state"] not in {"verified", "reused"}
    ):
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
        row for row in rows
        if row["object_type"] == "file"
        and row["state"] in {"planned", "transferred", "missing", "quarantined"}
    ]
    # Publish checksum manifests first, then bind their MD5 entries before the
    # remaining files are started.
    checksum_rows = [row for row in files if PurePosixPath(str(row["path"])).name == "md5sum.txt"]
    other_rows = [row for row in files if row not in checksum_rows]
    with ThreadPoolExecutor(max_workers=concurrency) as executor:
        list(executor.map(lambda row: process_file(database, root, row), checksum_rows))
    write_progress(database, root)
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
    other_paths = {str(row["path"]) for row in other_rows}
    cas_by_md5: dict[tuple[int, str], tuple[Path, str]] = {}
    if other_paths:
        # Bind every published provider checksum before any non-manifest
        # payload.  A terminal run has no consumer for these bindings and must
        # not rewrite all 47,870 durable rows merely to record completion.
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
            now = utc_now()
            for path, (algorithm, checksum) in checksum_map.items():
                connection.execute(
                    "UPDATE objects SET checksum_algorithm = ?, checksum = ?, updated_at_utc = ? WHERE path = ?",
                    (algorithm, checksum, now, path),
                )
            connection.commit()
        finally:
            connection.close()
        other_rows = [
            row for row in database_rows(database) if str(row["path"]) in other_paths
        ]
        # Revalidating the shared CAS is intentionally expensive: it reads
        # every byte twice (SHA-256 plus MD5).  Build the reuse index only when
        # this run actually has a checksum-bound payload that could consume it.
        if any(row["checksum_algorithm"] == "md5" for row in other_rows):
            cas_by_md5 = index_verified_cas(VGP_DATA_ROOT / "objects/sha256")
    import threading
    progress_lock = threading.Lock()
    last_progress = [time.monotonic()]

    def process_and_checkpoint(row: Mapping[str, object]) -> None:
        process_file(database, root, row, cas_by_md5)
        with progress_lock:
            now = time.monotonic()
            if now - last_progress[0] >= 60:
                write_progress(database, root)
                last_progress[0] = now

    outcome = "complete"
    detail = "frozen_inventory_verified"
    try:
        with ThreadPoolExecutor(max_workers=concurrency) as executor:
            list(executor.map(process_and_checkpoint, other_rows))
        incomplete = [
            row for row in database_rows(database)
            if row["object_type"] == "file"
            and row["state"] not in {"verified", "reused", "verified_upstream_conflict"}
        ]
        if incomplete:
            raise MirrorError(
                f"{len(incomplete)} frozen files remain unverified; first error: "
                f"{incomplete[0]['path']}: {incomplete[0]['state_reason']}"
            )
        exceptions = [
            row for row in database_rows(database)
            if row["state"] == "verified_upstream_conflict"
        ]
        if any(row["sequence_subset"] != "non_sequence_product_or_metadata" for row in exceptions):
            raise MirrorError("sequence object was incorrectly assigned a terminal conflict exception")
        if exceptions:
            detail = f"frozen_inventory_verified_with_{len(exceptions)}_upstream_conflict_exception(s)"
    except BaseException as error:
        outcome = "failed"
        detail = str(error)
        raise
    finally:
        connection = sqlite3.connect(database, timeout=60)
        try:
            completed_transferred = sum(
                int(row["transferred_bytes"])
                for row in database_rows(database)
                if row["object_type"] == "file"
            )
            connection.execute(
                """UPDATE worker_runs SET completed_at_utc = ?, outcome = ?, detail = ?,
                   completed_transferred_bytes = ? WHERE run_id = ?""",
                (utc_now(), outcome, detail, completed_transferred, run_id),
            )
            connection.commit()
        finally:
            connection.close()
        write_progress(database, root)
        write_exception_ledger(database, root, EXCEPTION_LEDGER_OUTPUT)


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
        "canonical_vgp_root": str(VGP_DATA_ROOT),
        "mirror_root": str(root),
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
    ignored_provider_orphans: list[str] = []
    checksums = parse_upstream_md5s(
        args.root / "objects", objects, ignored_provider_orphans
    )
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
        exception_ledger_output=args.exception_ledger_output,
        ignored_provider_checksum_orphans=ignored_provider_orphans,
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
    build.add_argument(
        "--exception-ledger-output",
        type=Path,
        default=PROJECT_ROOT / "analysis/vgp_freeze1_exception_ledger.json",
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
    adjudicate = subcommands.add_parser(
        "adjudicate-conflict",
        help="independently reproduce a quarantined metadata conflict and resolve against NCBI",
    )
    adjudicate.add_argument("--root", type=Path, default=RELEASE_ROOT)
    adjudicate.add_argument("--database", type=Path, default=RELEASE_ROOT / "state/mirror.sqlite3")
    adjudicate.add_argument("--inventory-id", required=True)
    adjudicate.add_argument("--alternate-source", required=True)
    adjudicate.add_argument("--alternate-checksum-source", required=True)
    adjudicate.add_argument("--ledger-output", type=Path, default=EXCEPTION_LEDGER_OUTPUT)
    adjudicate_all = subcommands.add_parser(
        "adjudicate-conflicts",
        help="resolve all eligible quarantined assembly reports through official NCBI indexes",
    )
    adjudicate_all.add_argument("--root", type=Path, default=RELEASE_ROOT)
    adjudicate_all.add_argument(
        "--database", type=Path, default=RELEASE_ROOT / "state/mirror.sqlite3"
    )
    adjudicate_all.add_argument("--ledger-output", type=Path, default=EXCEPTION_LEDGER_OUTPUT)
    status = subcommands.add_parser("status", help="write exact restart-safe live progress")
    status.add_argument("--root", type=Path, default=RELEASE_ROOT)
    status.add_argument("--database", type=Path, default=RELEASE_ROOT / "state/mirror.sqlite3")
    status.add_argument("--output", type=Path, default=RELEASE_ROOT / "state/progress.json")
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
        elif args.command == "adjudicate-conflict":
            payload = adjudicate_upstream_conflict(
                args.database,
                args.root,
                args.inventory_id,
                args.alternate_source,
                args.alternate_checksum_source,
                args.ledger_output,
            )
        elif args.command == "adjudicate-conflicts":
            payload = adjudicate_current_metadata_conflicts(
                args.database, args.root, args.ledger_output
            )
        elif args.command == "status":
            payload = write_progress(args.database, args.root, args.output)
        else:  # pragma: no cover
            raise AssertionError(args.command)
    except (MirrorError, OSError, sqlite3.Error) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
