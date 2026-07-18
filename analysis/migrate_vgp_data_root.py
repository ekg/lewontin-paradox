#!/usr/bin/env python3
"""Relink verified legacy VGP objects into the canonical shared CAS.

The legacy project-named root is read-only migration input.  This command never
downloads and never removes a source object.  On the shared MooseFS device it
creates hard links at canonical digest-derived paths, after size and SHA-256
verification, so already acquired bytes are reused without duplication.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import sys
from pathlib import Path
from typing import Any, Sequence


ROOT = Path(__file__).resolve().parents[1]
INVENTORY = ROOT / "analysis/vgp_10_pilot_object_inventory.tsv"
OUTPUT = ROOT / "analysis/vgp_data_root_migration_v1.json"
CANONICAL = Path("/moosefs/erikg/vgp")
LEGACY = Path("/moosefs/erikg/lewontin-paradox-data/vgp/phase1-freeze-1.0")
CATALOG_NAME = "VGPPhase1-freeze-1.0.commit-dc1b2af5a7741b97d66fb10cb2bce97f41765cdf.tsv"
CATALOG_SHA256 = "9c58420484a8b76a2d6175b7c26bf709e68bdc726a67fc7541b8c2b5a2fc13a4"


class MigrationError(RuntimeError):
    """Verified bytes cannot be safely relinked into the canonical root."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verified_sources(inventory: Path = INVENTORY) -> list[dict[str, Any]]:
    with inventory.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    by_digest: dict[str, dict[str, Any]] = {}
    for row in rows:
        if row["status"] not in {"reused", "verified"}:
            continue
        digest = row["local_sha256"]
        source = Path(row["local_path"])
        if len(digest) != 64 or not source.is_absolute():
            raise MigrationError(f"invalid verified object record: {row['object_id']}")
        try:
            source.resolve().relative_to(LEGACY.resolve())
        except ValueError as error:
            raise MigrationError(f"migration source escapes legacy root: {source}") from error
        record = {
            "sha256": digest,
            "bytes": int(row["observed_bytes"]),
            "source": source,
            "object_ids": [row["object_id"]],
        }
        if digest in by_digest:
            if by_digest[digest]["bytes"] != record["bytes"] or by_digest[digest]["source"] != source:
                raise MigrationError(f"conflicting records for digest: {digest}")
            by_digest[digest]["object_ids"].append(row["object_id"])
        else:
            by_digest[digest] = record
    return [by_digest[digest] for digest in sorted(by_digest)]


def migrate(apply: bool, verify_hashes: bool = True) -> dict[str, Any]:
    records = verified_sources()
    linked = reused = bytes_total = 0
    output_rows = []
    for row in records:
        source: Path = row["source"]
        digest = row["sha256"]
        if not source.is_file() or source.stat().st_size != row["bytes"]:
            raise MigrationError(f"legacy verified source missing/size drift: {source}")
        if verify_hashes and sha256_file(source) != digest:
            raise MigrationError(f"legacy verified source digest drift: {source}")
        destination = CANONICAL / "objects/sha256" / digest[:2] / digest[2:4] / digest
        disposition = "planned_hardlink"
        same_inode = False
        if destination.exists():
            if not destination.is_file() or destination.stat().st_size != row["bytes"]:
                raise MigrationError(f"canonical destination conflicts: {destination}")
            if verify_hashes and sha256_file(destination) != digest:
                raise MigrationError(f"canonical destination digest drift: {destination}")
            same_inode = source.stat().st_ino == destination.stat().st_ino and source.stat().st_dev == destination.stat().st_dev
            disposition = "reused_existing_hardlink" if same_inode else "reused_existing_verified_copy"
            reused += 1
        elif apply:
            destination.parent.mkdir(parents=True, exist_ok=True)
            if source.stat().st_dev != destination.parent.stat().st_dev:
                raise MigrationError("legacy and canonical CAS are not on the same filesystem; hardlink refused")
            os.link(source, destination)
            same_inode = source.stat().st_ino == destination.stat().st_ino
            if not same_inode:
                raise MigrationError(f"hardlink identity check failed: {destination}")
            disposition = "hardlinked_verified"
            linked += 1
        bytes_total += row["bytes"]
        output_rows.append({
            "sha256": digest,
            "bytes": row["bytes"],
            "canonical_path": str(destination),
            "legacy_migration_source": str(source),
            "disposition": disposition,
            "same_inode": same_inode,
            "object_ids": sorted(row["object_ids"]),
        })

    catalog_source = LEGACY / "manifests" / CATALOG_NAME
    catalog_destination = CANONICAL / "manifests" / CATALOG_NAME
    if not catalog_source.is_file() or sha256_file(catalog_source) != CATALOG_SHA256:
        raise MigrationError("pinned catalog migration source failed digest gate")
    catalog_disposition = "planned_hardlink"
    if catalog_destination.exists():
        if sha256_file(catalog_destination) != CATALOG_SHA256:
            raise MigrationError("canonical pinned catalog digest drift")
        catalog_disposition = "reused_existing"
    elif apply:
        catalog_destination.parent.mkdir(parents=True, exist_ok=True)
        os.link(catalog_source, catalog_destination)
        catalog_disposition = "hardlinked_verified"

    return {
        "schema_version": "vgp-data-root-migration-v1.0.0",
        "canonical_root": str(CANONICAL),
        "legacy_root_role": "migration_input_only",
        "legacy_migration_input": str(LEGACY),
        "operation": "verified_same-filesystem_hardlink",
        "applied": apply,
        "network_requests": 0,
        "downloaded_bytes": 0,
        "source_objects_removed": 0,
        "unique_verified_objects": len(records),
        "logical_verified_bytes": bytes_total,
        "new_hardlinks": linked,
        "reused_canonical_objects": reused,
        "objects": output_rows,
        "catalog": {
            "sha256": CATALOG_SHA256,
            "canonical_path": str(catalog_destination),
            "legacy_migration_source": str(catalog_source),
            "disposition": catalog_disposition,
        },
        "active_output_root": str(CANONICAL),
        "active_output_paths_under_legacy_root": 0,
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--apply", action="store_true")
    parser.add_argument("--skip-rehash", action="store_true")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args(argv)
    try:
        result = migrate(args.apply, not args.skip_rehash)
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print(json.dumps({key: result[key] for key in (
            "canonical_root", "unique_verified_objects", "logical_verified_bytes",
            "new_hardlinks", "reused_canonical_objects", "downloaded_bytes",
        )}, sort_keys=True))
        return 0
    except (MigrationError, OSError, ValueError) as error:
        print(f"migration error: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
