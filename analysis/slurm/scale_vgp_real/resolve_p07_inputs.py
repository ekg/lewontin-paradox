#!/usr/bin/env python3
"""Fail-closed resolution of exact P07 CAS inputs and pinned FastGA tools."""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path


def digest(path: Path) -> str:
    value = hashlib.sha256(path.read_bytes()).hexdigest()
    return value


manifest_path, capture_path, amendment_path = map(Path, sys.argv[1:])
manifest = json.loads(manifest_path.read_text())
capture = json.loads(capture_path.read_text())
amendment = json.loads(amendment_path.read_text())
if manifest.get("canonical_vgp_root") != "/moosefs/erikg/vgp":
    raise SystemExit("noncanonical P07 input manifest")
if manifest.get("selection_id") != "P07":
    raise SystemExit("wrong selection")
for side in ("h1", "h2"):
    row = manifest["inputs"][side]
    path = Path(row["canonical_cas_path"])
    if not str(path).startswith("/moosefs/erikg/vgp/") or not path.is_file():
        raise SystemExit(f"noncanonical or absent {side} CAS object")
    if digest(path) != row["compressed_sha256"] or path.stat().st_size != row["compressed_bytes"]:
        raise SystemExit(f"{side} CAS integrity failure")
rows = {row["name"]: row for row in capture["executables"]}
sweepga = Path(rows["sweepga"]["path"])
if not sweepga.is_file() or digest(sweepga) != rows["sweepga"]["sha256"]:
    raise SystemExit("SweepGA capture mismatch")
for name, row in amendment["companions"].items():
    path = Path(row["path"])
    if not path.is_file() or digest(path) != row["sha256"]:
        raise SystemExit(f"FastGA amendment mismatch: {name}")
print(manifest["inputs"]["h1"]["canonical_cas_path"])
print(manifest["inputs"]["h2"]["canonical_cas_path"])
print(manifest["inputs"]["h1"]["decompressed_sha256"])
print(manifest["inputs"]["h2"]["decompressed_sha256"])
print(sweepga)
print(digest(capture_path))
print(digest(amendment_path))
