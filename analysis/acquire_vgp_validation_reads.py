#!/usr/bin/env python3
"""Acquire a checksum-locked, exact-individual VGP validation-read subset.

The immutable plan is repository data.  This program resolves every storage
path from ``analysis/vgp_data_root_config.json``, retains interrupted transfers
under the configured staging tree, verifies ENA's object MD5 and byte count,
computes a local SHA-256, and atomically promotes into the shared VGP CAS.  It
never writes biological data beneath the project-named legacy root.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import sys
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Mapping, Sequence


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
ROOT_CONFIG = REPOSITORY_ROOT / "analysis/vgp_data_root_config.json"
PLAN = REPOSITORY_ROOT / "analysis/vgp_validation_read_plan_v1.json"
OUTPUT = REPOSITORY_ROOT / "analysis/vgp_validation_reads_manifest_v1.json"
CANONICAL_MANIFEST_NAME = "vgp_validation_reads_manifest_v1.json"
HEX32 = re.compile(r"^[0-9a-f]{32}$")
HEX64 = re.compile(r"^[0-9a-f]{64}$")
VERSIONED_ASSEMBLY = re.compile(r"^GCA_[0-9]+\.[0-9]+$")
SAFE_ACCESSION = re.compile(r"^SRR[0-9]+$")
DEFAULT_RESERVE_BYTES = 5 * 1024**3


class AcquisitionError(RuntimeError):
    """The exact plan, filesystem, or verified content contract failed."""


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def load_json(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise AcquisitionError(f"expected JSON object: {path}")
    return value


def load_plan(path: Path = PLAN) -> dict[str, Any]:
    return load_json(path)


def configured_root(config_path: Path = ROOT_CONFIG) -> Path:
    config = load_json(config_path)
    root = Path(str(config.get("root", "")))
    if not root.is_absolute():
        raise AcquisitionError("configured VGP data root must be absolute")
    layout = config.get("layout", {})
    if not isinstance(layout, dict) or not (layout.get("raw_reads") or layout.get("accession_views")):
        raise AcquisitionError("root configuration lacks a raw-read/accession view layout")
    return root


def atomic_json(path: Path, value: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.partial-{os.getpid()}")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)
    fsync_directory(path.parent)


def fsync_directory(path: Path) -> None:
    descriptor = os.open(path, os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def hash_file(path: Path) -> tuple[str, str, int]:
    sha256 = hashlib.sha256()
    md5 = hashlib.md5()  # noqa: S324 - ENA transport/content checksum only
    observed = 0
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            sha256.update(block)
            md5.update(block)
            observed += len(block)
    return sha256.hexdigest(), md5.hexdigest(), observed


def cas_path(
    data_root: Path, digest: str, immutable_objects_relative: str = "objects/sha256"
) -> Path:
    if not HEX64.fullmatch(digest):
        raise AcquisitionError(f"invalid SHA-256 identity: {digest!r}")
    return data_root / immutable_objects_relative / digest[:2] / digest[2:4] / digest


def _safe_name(value: object) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value))


def unlink_if_exists(path: Path) -> None:
    """Path.unlink(missing_ok=...) compatibility for cluster Python 3.7."""
    try:
        path.unlink()
    except FileNotFoundError:
        pass


def partial_path(
    data_root: Path, spec: Mapping[str, object], staging_partials_relative: str = "staging/partials"
) -> Path:
    return data_root / staging_partials_relative / "vgp-validation-reads-v1" / f"{_safe_name(spec['object_id'])}.partial"


def read_view_path(
    data_root: Path, spec: Mapping[str, object], raw_reads_relative: str = "raw/reads"
) -> Path:
    return data_root / raw_reads_relative / str(spec["selection_id"]) / str(spec["run_accession"]) / str(spec["filename"])


def quarantine(
    path: Path, data_root: Path, reason: str, quarantine_relative: str = "quarantine"
) -> Path:
    directory = data_root / quarantine_relative / "vgp-validation-reads-v1" / utc_now().replace(":", "")
    directory.mkdir(parents=True, exist_ok=True)
    target = directory / f"{path.name}.{_safe_name(reason)[:96]}"
    os.replace(path, target)
    fsync_directory(directory)
    return target


def _atomic_symlink(target: Path, view: Path) -> None:
    view.parent.mkdir(parents=True, exist_ok=True)
    temporary = view.with_name(f".{view.name}.partial-{os.getpid()}")
    unlink_if_exists(temporary)
    relative_target = os.path.relpath(target, view.parent)
    temporary.symlink_to(relative_target)
    os.replace(temporary, view)
    fsync_directory(view.parent)


def validate_plan(plan: Mapping[str, Any]) -> None:
    if plan.get("schema_version") != "vgp-validation-read-plan-v1.0.0":
        raise AcquisitionError("unexpected validation-read plan version")
    if plan.get("canonical_root") != "/moosefs/erikg/vgp":
        raise AcquisitionError("plan canonical root drift")
    pairs_raw = plan.get("pairs")
    objects_raw = plan.get("objects")
    if not isinstance(pairs_raw, list) or not isinstance(objects_raw, list):
        raise AcquisitionError("plan pairs/objects must be arrays")
    pairs = {str(row.get("selection_id")): row for row in pairs_raw if isinstance(row, dict)}
    if set(pairs) != {"P04", "P07", "P09"} or len(pairs) != len(pairs_raw):
        raise AcquisitionError("plan must contain exact unique P04/P07/P09 pair strata")
    seen: set[str] = set()
    for spec in objects_raw:
        if not isinstance(spec, dict):
            raise AcquisitionError("object record is not a JSON object")
        object_id = str(spec.get("object_id", ""))
        if not object_id or object_id in seen:
            raise AcquisitionError(f"duplicate or blank object id: {object_id!r}")
        seen.add(object_id)
        selection = str(spec.get("selection_id", ""))
        if selection not in pairs:
            raise AcquisitionError(f"{object_id}: selection is outside exact pair roster")
        pair = pairs[selection]
        for field in (
            "species", "taxid", "biosample", "individual_or_isolate",
            "h1_accession_version", "h2_accession_version",
        ):
            if spec.get(field) != pair.get(field):
                raise AcquisitionError(f"{object_id}: {field} does not bind exact pilot pair")
        if not SAFE_ACCESSION.fullmatch(str(spec.get("run_accession", ""))):
            raise AcquisitionError(f"{object_id}: invalid run accession")
        for field in ("h1_accession_version", "h2_accession_version"):
            if not VERSIONED_ASSEMBLY.fullmatch(str(spec.get(field, ""))):
                raise AcquisitionError(f"{object_id}: assembly accession is not exact/versioned")
        if spec.get("platform") not in {"PACBIO_SMRT", "ILLUMINA"}:
            raise AcquisitionError(f"{object_id}: unsupported platform")
        if spec.get("library_strategy") != "WGS" or spec.get("library_source") != "GENOMIC":
            raise AcquisitionError(f"{object_id}: validation payload is not genomic WGS")
        if spec.get("upstream_checksum_algorithm") != "md5" or not HEX32.fullmatch(
            str(spec.get("upstream_checksum", ""))
        ):
            raise AcquisitionError(f"{object_id}: exact ENA MD5 is absent")
        if int(spec.get("expected_bytes") or 0) <= 0:
            raise AcquisitionError(f"{object_id}: expected bytes are not positive")
        if not str(spec.get("source_url", "")).startswith("https://ftp.sra.ebi.ac.uk/"):
            raise AcquisitionError(f"{object_id}: source is not the authoritative ENA HTTPS archive")


# Keep this runtime-evaluated alias compatible with the cluster host Python;
# function annotations themselves are postponed by ``__future__.annotations``.
Downloader = Callable[[dict, Path, int], int]


def http_downloader(spec: dict[str, object], path: Path, offset: int) -> int:
    """Download with byte-range resume and bounded retry, retaining partials."""
    expected = int(spec["expected_bytes"])
    initial_offset = offset
    attempts_without_progress = 0
    while offset < expected:
        before_attempt = offset
        request = urllib.request.Request(
            str(spec["source_url"]),
            headers={
                "User-Agent": "vgp-validation-read-acquisition/1.0",
                "Accept-Encoding": "identity",
                "Range": f"bytes={offset}-",
            },
        )
        try:
            with urllib.request.urlopen(request, timeout=300) as response:  # noqa: S310 - frozen exact ENA URL
                status = int(getattr(response, "status", response.getcode()))
                if offset and status != 206:
                    unlink_if_exists(path)
                    offset = 0
                mode = "ab" if offset else "wb"
                with path.open(mode) as handle:
                    while True:
                        block = response.read(8 * 1024 * 1024)
                        if not block:
                            break
                        handle.write(block)
                        offset += len(block)
                    handle.flush()
                    os.fsync(handle.fileno())
        except (OSError, TimeoutError, urllib.error.URLError) as error:
            current = path.stat().st_size if path.exists() else 0
            if current > before_attempt:
                attempts_without_progress = 0
            else:
                attempts_without_progress += 1
            offset = current
            if attempts_without_progress >= 8:
                raise AcquisitionError(
                    f"ENA transfer failed after bounded retries at {offset}/{expected} bytes: {error}"
                ) from error
            time.sleep(min(2 ** attempts_without_progress, 30))
            continue
        current = path.stat().st_size
        offset = current
        if offset == before_attempt:
            attempts_without_progress += 1
            if attempts_without_progress >= 8:
                raise AcquisitionError(f"ENA transfer made no progress at {offset}/{expected} bytes")
            time.sleep(min(2 ** attempts_without_progress, 30))
        else:
            attempts_without_progress = 0
    return max(0, offset - initial_offset)


def _verified_target(
    target: Path, spec: Mapping[str, object]
) -> tuple[bool, str, int]:
    if not target.is_file():
        return False, "", 0
    sha256, md5, size = hash_file(target)
    valid = (
        target.name == sha256
        and size == int(spec["expected_bytes"])
        and md5 == str(spec["upstream_checksum"])
    )
    return valid, sha256, size


def acquire_object(
    object_spec: Mapping[str, object],
    data_root: Path,
    *,
    downloader: Downloader = http_downloader,
    known_sha256: str = "",
    raw_reads_relative: str = "raw/reads",
    immutable_objects_relative: str = "objects/sha256",
    staging_partials_relative: str = "staging/partials",
    quarantine_relative: str = "quarantine",
    reused_source: str = "",
    recovered_unmanifested_promotion: bool = False,
) -> dict[str, object]:
    spec: dict[str, object] = dict(object_spec)
    expected = int(spec["expected_bytes"])
    result: dict[str, object] = dict(spec)
    result.update(
        status="planned", status_detail="", observed_bytes=0, transferred_bytes=0,
        invocation_transferred_bytes=0,
        resume_from_bytes=0, local_sha256="", local_path="", accession_view_path="",
        quarantine_path="", reused_source=reused_source, verified_utc="",
        verification={"expected_bytes_match": False, "ena_md5_match": False,
                      "local_sha256_path_match": False, "gzip_magic": False,
                      "post_promotion_reverified": False},
    )
    view = read_view_path(data_root, spec, raw_reads_relative)

    candidate_targets: list[Path] = []
    if HEX64.fullmatch(known_sha256):
        candidate_targets.append(cas_path(data_root, known_sha256, immutable_objects_relative))
    if view.is_symlink():
        candidate_targets.append(view.resolve())
    for target in candidate_targets:
        valid, sha256, size = _verified_target(target, spec)
        if valid:
            _atomic_symlink(target, view)
            recovered = recovered_unmanifested_promotion
            result.update(
                status="verified" if recovered else "reused",
                status_detail=(
                    "recovered, rehashed, and published an atomically promoted CAS object after a manifest/view interruption"
                    if recovered else "canonical CAS object rehashed and reused"
                ),
                observed_bytes=size, local_sha256=sha256, local_path=str(target),
                accession_view_path=str(view), verified_utc=utc_now(),
                transferred_bytes=size if recovered else 0,
                resume_from_bytes=size if recovered else 0,
                verification={"expected_bytes_match": True, "ena_md5_match": True,
                              "local_sha256_path_match": True, "gzip_magic": _gzip_magic(target),
                              "post_promotion_reverified": True},
            )
            return result
        canonical_cas = data_root / immutable_objects_relative
        inside_cas = target == canonical_cas or canonical_cas in target.parents
        if target.exists() and inside_cas:
            bad = quarantine(target, data_root, "canonical-revalidation-mismatch", quarantine_relative)
            result["quarantine_path"] = str(bad)

    partial = partial_path(data_root, spec, staging_partials_relative)
    partial.parent.mkdir(parents=True, exist_ok=True)
    offset = partial.stat().st_size if partial.exists() else 0
    result["resume_from_bytes"] = offset
    if offset > expected:
        bad = quarantine(partial, data_root, "partial-oversize", quarantine_relative)
        result.update(
            status="quarantined", observed_bytes=offset, quarantine_path=str(bad),
            status_detail=f"resumable partial exceeds ENA bytes: {offset} > {expected}",
        )
        return result
    try:
        invocation_transferred = downloader(spec, partial, offset)
        result["invocation_transferred_bytes"] = invocation_transferred
        # ``transferred_bytes`` is cumulative unique archive payload present
        # for this acquisition object, including an earlier resumable
        # invocation or an external multi-range fill of the same canonical
        # partial.  The invocation-only count remains separately auditable.
        result["transferred_bytes"] = 0 if reused_source else partial.stat().st_size
    except Exception as error:  # transport failures remain resumable, never promoted
        observed = partial.stat().st_size if partial.exists() else 0
        result.update(
            status="missing", observed_bytes=observed,
            transferred_bytes=0 if reused_source else observed,
            status_detail=f"transfer incomplete; partial retained for resume: {error}",
        )
        return result

    sha256, md5, size = hash_file(partial)
    checks = {
        "expected_bytes_match": size == expected,
        "ena_md5_match": md5 == str(spec["upstream_checksum"]),
        "local_sha256_path_match": True,
        "gzip_magic": _gzip_magic(partial),
        "post_promotion_reverified": False,
    }
    result.update(observed_bytes=size, local_sha256=sha256, verification=checks)
    failures = [name for name, passed in checks.items() if name != "post_promotion_reverified" and not passed]
    if failures:
        bad = quarantine(
            partial, data_root, "content-mismatch-" + "-".join(failures), quarantine_relative
        )
        result.update(
            status="quarantined", quarantine_path=str(bad),
            status_detail="failed pre-promotion verification: " + ", ".join(failures),
        )
        return result

    target = cas_path(data_root, sha256, immutable_objects_relative)
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists():
        valid, _, _ = _verified_target(target, spec)
        if not valid:
            existing_bad = quarantine(
                target, data_root, "CAS-collision-revalidation-mismatch", quarantine_relative
            )
            staged_bad = quarantine(
                partial, data_root, "CAS-collision-staged-source", quarantine_relative
            )
            result.update(
                status="quarantined", quarantine_path=f"{existing_bad};{staged_bad}",
                status_detail="existing digest path failed byte/MD5/SHA revalidation",
            )
            return result
        partial.unlink()
        disposition = "reused"
        detail = "concurrent canonical CAS object rehashed and reused"
    else:
        os.replace(partial, target)
        os.chmod(target, 0o440)
        fsync_directory(target.parent)
        disposition = "reused" if reused_source else "verified"
        detail = (
            f"verified legacy object reused from {reused_source} and atomically promoted"
            if reused_source
            else "ENA bytes staged, verified, and atomically promoted read-only"
        )
    valid, promoted_sha, promoted_size = _verified_target(target, spec)
    if not valid:
        bad = quarantine(
            target, data_root, "post-promotion-revalidation-mismatch", quarantine_relative
        )
        result.update(
            status="quarantined", local_path="", quarantine_path=str(bad),
            status_detail="post-promotion size/ENA-MD5/local-SHA revalidation failed",
        )
        return result
    _atomic_symlink(target, view)
    checks["post_promotion_reverified"] = True
    result.update(
        status=disposition, status_detail=detail, observed_bytes=promoted_size,
        local_sha256=promoted_sha, local_path=str(target), accession_view_path=str(view),
        verified_utc=utc_now(), verification=checks,
    )
    return result


def _gzip_magic(path: Path) -> bool:
    with path.open("rb") as handle:
        return handle.read(2) == b"\x1f\x8b"


def summarize(rows: Sequence[Mapping[str, object]]) -> dict[str, Any]:
    def bucket(selected: Sequence[Mapping[str, object]], byte_field: str = "expected_bytes") -> dict[str, int]:
        return {
            "objects": len(selected),
            "bytes": sum(int(row.get(byte_field) or 0) for row in selected),
        }

    all_rows = list(rows)
    transferred = [row for row in all_rows if int(row.get("transferred_bytes") or 0) > 0]
    locally_verified = [row for row in all_rows if row.get("status") in {"verified", "reused"}]
    summary: dict[str, Any] = {
        "planned": bucket(all_rows),
        "transferred": bucket(transferred, "transferred_bytes"),
        "verified": bucket(locally_verified),
        "newly_promoted": bucket([row for row in all_rows if row.get("status") == "verified"]),
    }
    for status in ("reused", "missing", "quarantined"):
        summary[status] = bucket([row for row in all_rows if row.get("status") == status])
    summary["pending"] = bucket([row for row in all_rows if row.get("status") == "planned"])
    terminal = summary["newly_promoted"]["objects"] + sum(
        summary[status]["objects"] for status in ("reused", "missing", "quarantined", "pending")
    )
    terminal_bytes = summary["newly_promoted"]["bytes"] + sum(
        summary[status]["bytes"] for status in ("reused", "missing", "quarantined", "pending")
    )
    summary["accounting"] = {
        "objects_reconciled": terminal == summary["planned"]["objects"],
        "bytes_reconciled": terminal_bytes == summary["planned"]["bytes"],
        "status_partition_objects": terminal,
        "status_partition_bytes": terminal_bytes,
    }
    return summary


def inventory_legacy(config: Mapping[str, Any], objects: Sequence[Mapping[str, object]]) -> dict[str, Any]:
    legacy_value = str(config.get("migration_input_only", ""))
    legacy = Path(legacy_value) if legacy_value else Path("/nonexistent")
    wanted = {str(row["filename"]): row for row in objects}
    candidates: list[dict[str, object]] = []
    if legacy.is_dir():
        for directory, _subdirs, filenames in os.walk(legacy):
            for filename in set(filenames) & set(wanted):
                path = Path(directory) / filename
                candidates.append({"path": str(path), "bytes": path.stat().st_size, "filename": filename})
    return {
        "canonical_root": str(config.get("root", "")),
        "legacy_root_role": "migration_input_only",
        "legacy_root": legacy_value,
        "matching_candidate_objects": len(candidates),
        "matching_candidate_bytes": sum(int(row["bytes"]) for row in candidates),
        "candidates": candidates,
    }


def _seed_verified_legacy(
    legacy_inventory: Mapping[str, Any], spec: Mapping[str, object], data_root: Path,
    staging_partials_relative: str = "staging/partials",
) -> str:
    partial = partial_path(data_root, spec, staging_partials_relative)
    if partial.exists():
        return ""
    for candidate in legacy_inventory.get("candidates", []):
        path = Path(str(candidate.get("path", "")))
        if path.name != spec["filename"] or not path.is_file():
            continue
        _sha, md5, size = hash_file(path)
        if size != int(spec["expected_bytes"]) or md5 != spec["upstream_checksum"]:
            continue
        partial.parent.mkdir(parents=True, exist_ok=True)
        try:
            os.link(path, partial)
        except OSError:
            shutil.copyfile(path, partial)
        return str(path)
    return ""


def prior_sha256_by_id(path: Path) -> dict[str, str]:
    if not path.is_file():
        return {}
    prior = load_json(path)
    result: dict[str, str] = {}
    for row in prior.get("objects", []):
        digest = str(row.get("local_sha256", ""))
        if HEX64.fullmatch(digest):
            result[str(row.get("object_id", ""))] = digest
    return result


def recover_unmanifested_cas_sha(
    data_root: Path, spec: Mapping[str, object], immutable_objects_relative: str
) -> str:
    """Recover a verified promotion if interruption preceded view/manifest."""
    root = data_root / immutable_objects_relative
    expected = int(spec["expected_bytes"])
    if not root.is_dir():
        return ""
    pattern = "[0-9a-f][0-9a-f]/[0-9a-f][0-9a-f]/[0-9a-f]*"
    for candidate in root.glob(pattern):
        if not candidate.is_file() or not HEX64.fullmatch(candidate.name):
            continue
        if candidate.stat().st_size != expected:
            continue
        valid, digest, _size = _verified_target(candidate, spec)
        if valid:
            return digest
    return ""


def capacity_record(data_root: Path) -> dict[str, int]:
    usage = shutil.disk_usage(data_root)
    return {"total_bytes": usage.total, "used_bytes": usage.used, "available_bytes": usage.free}


def verify_manifest(
    manifest_path: Path = OUTPUT, config_path: Path = ROOT_CONFIG
) -> dict[str, Any]:
    """Independently rehash the published acquisition and its full ledger."""
    config = load_json(config_path)
    data_root = configured_root(config_path)
    layout = config["layout"]
    manifest = load_json(manifest_path)
    errors: list[str] = []
    if manifest.get("schema_version") != "vgp-validation-reads-manifest-v1.0.0":
        errors.append("manifest schema version mismatch")
    if manifest.get("canonical_root") != str(data_root):
        errors.append("manifest canonical root does not match repository configuration")
    rows = manifest.get("objects", [])
    if not isinstance(rows, list):
        raise AcquisitionError("manifest objects is not an array")
    recomputed = summarize(rows)
    if manifest.get("summary") != recomputed:
        errors.append("manifest status/object/byte summary does not recompute exactly")
    verified_objects = 0
    verified_bytes = 0
    allowed = {"verified", "reused", "missing", "quarantined", "planned"}
    cas_root = (data_root / str(layout["immutable_objects"])).resolve()
    raw_relative = str(layout.get("raw_reads") or layout["accession_views"])
    raw_root = (data_root / raw_relative).resolve()
    for row in rows:
        object_id = str(row.get("object_id", "<missing>"))
        status = str(row.get("status", ""))
        if status not in allowed:
            errors.append(f"{object_id}: invalid status {status!r}")
            continue
        if status not in {"verified", "reused"}:
            continue
        target = Path(str(row.get("local_path", "")))
        view = Path(str(row.get("accession_view_path", "")))
        if not target.is_file():
            errors.append(f"{object_id}: verified CAS target is absent")
            continue
        resolved_target = target.resolve()
        if cas_root not in resolved_target.parents:
            errors.append(f"{object_id}: verified target escapes configured CAS")
            continue
        if not view.is_symlink() or view.resolve() != resolved_target:
            errors.append(f"{object_id}: configured raw-read accession view is absent or wrong")
        elif raw_root not in view.resolve(strict=False).parents and raw_root not in view.parent.resolve().parents:
            # The symlink resolves into CAS; its own parent, not its target,
            # establishes containment in the configured raw-read view.
            errors.append(f"{object_id}: accession view escapes configured raw-read tree")
        sha256, md5, size = hash_file(target)
        if sha256 != row.get("local_sha256") or target.name != sha256:
            errors.append(f"{object_id}: local SHA-256 or CAS identity mismatch")
        if md5 != row.get("upstream_checksum"):
            errors.append(f"{object_id}: ENA MD5 mismatch")
        if size != int(row.get("expected_bytes") or 0) or size != int(row.get("observed_bytes") or 0):
            errors.append(f"{object_id}: expected/observed byte mismatch")
        if not _gzip_magic(target):
            errors.append(f"{object_id}: gzip magic mismatch")
        if target.stat().st_mode & 0o222:
            errors.append(f"{object_id}: immutable CAS target is writable")
        verified_objects += 1
        verified_bytes += size
    if verified_objects == 0 or verified_bytes == 0:
        errors.append("zero real raw-read payloads verified")
    canonical_manifest = data_root / str(layout["pilot_manifests"]) / CANONICAL_MANIFEST_NAME
    if not canonical_manifest.is_file():
        errors.append("canonical downstream manifest copy is absent")
    elif hashlib.sha256(canonical_manifest.read_bytes()).hexdigest() != hashlib.sha256(
        manifest_path.read_bytes()
    ).hexdigest():
        errors.append("repository and canonical manifest copies differ")
    result = {
        "verified": not errors,
        "verified_objects": verified_objects,
        "verified_bytes": verified_bytes,
        "errors": errors,
        "canonical_root": str(data_root),
        "manifest": str(manifest_path),
        "verified_utc": utc_now(),
    }
    if errors:
        raise AcquisitionError("; ".join(errors))
    return result


def run(
    *, output: Path = OUTPUT, only: set[str] | None = None, delay_seconds: float = 2.0,
    reserve_bytes: int = DEFAULT_RESERVE_BYTES,
) -> dict[str, Any]:
    started_utc = utc_now()
    config = load_json(ROOT_CONFIG)
    data_root = configured_root(ROOT_CONFIG)
    plan = load_plan()
    validate_plan(plan)
    if str(data_root) != str(plan["canonical_root"]):
        raise AcquisitionError("repository root config does not match plan canonical root")
    data_root.mkdir(parents=True, exist_ok=True)
    layout = config["layout"]
    for key in ("immutable_objects", "staging_partials", "quarantine", "pilot_manifests"):
        relative = str(layout.get(key, ""))
        if not relative or Path(relative).is_absolute() or ".." in Path(relative).parts:
            raise AcquisitionError(f"unsafe or absent configured layout path: {key}")
        (data_root / relative).mkdir(parents=True, exist_ok=True)
    raw_reads_relative = str(layout.get("raw_reads") or layout.get("accession_views") or "")
    if (
        not raw_reads_relative
        or Path(raw_reads_relative).is_absolute()
        or ".." in Path(raw_reads_relative).parts
    ):
        raise AcquisitionError("unsafe or absent configured raw-read/accession view path")
    (data_root / raw_reads_relative).mkdir(parents=True, exist_ok=True)

    selected_specs = [
        dict(row) for row in plan["objects"]
        if only is None or str(row["object_id"]) in only or str(row["selection_id"]) in only
    ]
    if only is not None and not selected_specs:
        raise AcquisitionError("--only matched no exact object id or selection id")
    capacity_before = capacity_record(data_root)
    planned_transfer_bytes = sum(int(row["expected_bytes"]) for row in selected_specs)
    if capacity_before["available_bytes"] < planned_transfer_bytes + reserve_bytes:
        raise AcquisitionError(
            "live filesystem capacity is insufficient: "
            f"available={capacity_before['available_bytes']} planned={planned_transfer_bytes} reserve={reserve_bytes}"
        )

    prior = prior_sha256_by_id(output)
    legacy = inventory_legacy(config, selected_specs)
    rows: list[dict[str, object]] = []
    selected_ids = {str(row["object_id"]) for row in selected_specs}
    for index, plan_row in enumerate(plan["objects"]):
        spec = dict(plan_row)
        if str(spec["object_id"]) not in selected_ids:
            row = dict(spec)
            row.update(
                status="planned", status_detail="not selected in this incremental invocation",
                observed_bytes=0, transferred_bytes=0, resume_from_bytes=0,
                invocation_transferred_bytes=0,
                local_sha256="", local_path="", accession_view_path="", quarantine_path="",
                reused_source="", verified_utc="",
            )
            rows.append(row)
            continue
        known_sha256 = prior.get(str(spec["object_id"]), "")
        recovered_unmanifested = False
        if not known_sha256:
            known_sha256 = recover_unmanifested_cas_sha(
                data_root, spec, str(layout["immutable_objects"])
            )
            recovered_unmanifested = bool(known_sha256)
        reused_source = _seed_verified_legacy(
            legacy, spec, data_root, str(layout["staging_partials"])
        )
        row = acquire_object(
            spec, data_root, known_sha256=known_sha256,
            raw_reads_relative=raw_reads_relative,
            immutable_objects_relative=str(layout["immutable_objects"]),
            staging_partials_relative=str(layout["staging_partials"]),
            quarantine_relative=str(layout["quarantine"]), reused_source=reused_source,
            recovered_unmanifested_promotion=recovered_unmanifested,
        )
        rows.append(row)
        if delay_seconds and index + 1 < len(plan["objects"]):
            time.sleep(delay_seconds)

    summary = summarize(rows)
    manifest: dict[str, Any] = {
        "schema_version": "vgp-validation-reads-manifest-v1.0.0",
        "task_id": "acquire-vgp-validation-reads",
        "canonical_root": str(data_root),
        "root_config": str(ROOT_CONFIG.relative_to(REPOSITORY_ROOT)),
        "legacy_root_role": "migration_input_only",
        "started_utc": started_utc,
        "completed_utc": utc_now(),
        "transfer_mode": "sequential HTTPS range-resume; 8 MiB streaming chunks; bounded retry; atomic CAS promotion",
        "primary_assembly_run_is_not_gated": True,
        "plan_sha256": hashlib.sha256(PLAN.read_bytes()).hexdigest(),
        "capacity_preflight": {
            **capacity_before,
            "planned_selected_transfer_bytes": planned_transfer_bytes,
            "reserved_free_bytes": reserve_bytes,
            "passed": True,
            "basis": "live shutil.disk_usage on the canonical MooseFS path; no quota-unavailable or global-byte policy",
        },
        "capacity_after": capacity_record(data_root),
        "legacy_inventory": legacy,
        "pairs": plan["pairs"],
        "objects": rows,
        "summary": summary,
    }
    canonical_output = data_root / str(layout["pilot_manifests"]) / CANONICAL_MANIFEST_NAME
    atomic_json(output, manifest)
    atomic_json(canonical_output, manifest)
    return manifest


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument(
        "--only", action="append", default=[],
        help="incrementally acquire an exact object_id or all objects for a selection_id; repeatable",
    )
    parser.add_argument("--delay-seconds", type=float, default=2.0)
    parser.add_argument("--reserve-bytes", type=int, default=DEFAULT_RESERVE_BYTES)
    parser.add_argument(
        "--verify-only", action="store_true",
        help="rehash the existing payloads/views and recompute the published ledger without network access",
    )
    args = parser.parse_args(argv)
    try:
        if args.verify_only:
            verification = verify_manifest(args.output)
            print(
                f"verified={verification['verified_objects']} objects/"
                f"{verification['verified_bytes']} bytes under {verification['canonical_root']}"
            )
            return 0
        manifest = run(
            output=args.output, only=set(args.only) if args.only else None,
            delay_seconds=args.delay_seconds, reserve_bytes=args.reserve_bytes,
        )
    except AcquisitionError as error:
        print(f"refusing acquisition: {error}", file=sys.stderr)
        return 2
    summary = manifest["summary"]
    print(
        f"planned={summary['planned']['objects']} objects/{summary['planned']['bytes']} bytes; "
        f"transferred={summary['transferred']['objects']} objects/{summary['transferred']['bytes']} bytes; "
        f"verified={summary['verified']['objects']} objects/{summary['verified']['bytes']} bytes; "
        f"reused={summary['reused']['objects']} objects/{summary['reused']['bytes']} bytes; "
        f"missing={summary['missing']['objects']} objects/{summary['missing']['bytes']} bytes; "
        f"quarantined={summary['quarantined']['objects']} objects/{summary['quarantined']['bytes']} bytes"
    )
    return 0 if summary["verified"]["objects"] > 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
