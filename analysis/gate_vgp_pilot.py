#!/usr/bin/env python3
"""Independently regenerate and enforce the repaired bounded VGP pilot gate.

This module is deliberately offline.  It audits committed metadata and immutable
cached official responses; it never downloads a biological payload, submits a
job, or performs demographic inference.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Sequence

from analysis.tier3_common import Tier3ValidationError, sha256_file


PROJECT_ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = PROJECT_ROOT / "analysis"
DEFAULT_MANIFEST = ANALYSIS / "vgp_pilot_manifest.tsv"
DEFAULT_REJECTIONS = ANALYSIS / "vgp_pilot_rejections.tsv"
DEFAULT_SIZE_BUDGET = ANALYSIS / "vgp_pilot_size_budget.tsv"
DEFAULT_RESOLUTION_INDEX = ANALYSIS / "vgp_resolution_cache/index.json"
DEFAULT_PRIOR_SEEDS = ANALYSIS / "vgp_resolution_cache/prior_refusal/vgp_pilot_manifest.tsv"
DEFAULT_FREEZE_PROVENANCE = ANALYSIS / "vgp_phase1_freeze_provenance.json"
DEFAULT_ROOT_CONFIG = ANALYSIS / "vgp_data_root_config.json"
DEFAULT_ROOT_VALIDATION = ANALYSIS / "vgp_data_root_validation.json"
DEFAULT_DECISIONS = ANALYSIS / "vertebrate_scaleout_decisions.tsv"
DEFAULT_EXECUTION_PLAN = ANALYSIS / "vertebrate_scaleout_execution_plan.md"
DEFAULT_RESOURCE_BUDGET = ANALYSIS / "vertebrate_scaleout_resource_budget.tsv"
DEFAULT_GUIX_CHANNELS = ANALYSIS / "guix/channels.scm"
DEFAULT_GUIX_MANIFEST = ANALYSIS / "guix/manifest.scm"
DEFAULT_GUIX_ENVIRONMENT = ANALYSIS / "pilot_results/guix_environment.json"
DEFAULT_GATE_JSON = ANALYSIS / "vgp_pilot_gate.json"
DEFAULT_GATE_REVIEW = ANALYSIS / "vgp_pilot_gate_review.md"

BYTES_PER_GIB = 1024**3
EXPECTED_CATALOG_STATISTICS = {
    "physical_lines": 717,
    "header_lines": 1,
    "data_rows": 716,
    "unique_species": 714,
    "data_row_excess_over_unique_species": 2,
    "duplicated_species": [
        {"scientific_name": "Lophostoma evotis", "multiplicity": 2},
        {"scientific_name": "Micronycteris microtis", "multiplicity": 2},
    ],
}
BASE_CAPS = {
    "species": (6.0, "count"),
    "compressed_inputs_gib": (120.0, "GiB"),
    "scratch_gib": (750.0, "GiB"),
    "core_hours": (1500.0, "core-hours"),
    "concurrent_species": (2.0, "count"),
    "memory_per_job_gib": (256.0, "GiB"),
}
EXPECTED_MEASUREMENT_THRESHOLDS = {
    "callable_bases": 10_000_000,
    "callable_fraction": 0.50,
    "queryable_gene_count": 1_000,
    "queryable_gene_bases": 1_000_000,
}
VERSIONED_ACCESSION = re.compile(r"^GC[AF]_\d+\.\d+$")
HEX32 = re.compile(r"^[0-9a-f]{32}$")
HEX40 = re.compile(r"^[0-9a-f]{40}$")
HEX64 = re.compile(r"^[0-9a-f]{64}$")
AFFIRMATIVE_PHASE_STATUSES = {
    "yes",
    "phased",
    "affirmative",
    "proven_phased_h1_h2",
    "official_phased_pair",
}
PAIR_MODALITIES = {"tier3a", "tier3a_diversity", "tier3a_paired_diversity"}


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def canonical_json_bytes(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("utf-8")


def sha256_json(value: Any) -> str:
    return hashlib.sha256(canonical_json_bytes(value)).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def load_tsv(path: Path) -> list[dict[str, str]]:
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_json(path: Path, value: Mapping[str, Any]) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _finite(value: Any, *, positive: bool = False) -> float | None:
    try:
        parsed = float(str(value).replace(",", ""))
    except (TypeError, ValueError):
        return None
    if not math.isfinite(parsed) or (positive and parsed <= 0):
        return None
    return parsed


def _integer(value: Any, *, positive: bool = False) -> int | None:
    parsed = _finite(value, positive=positive)
    if parsed is None or not parsed.is_integer():
        return None
    return int(parsed)


def _blocker(code: str, message: str, source: str) -> dict[str, str]:
    return {"code": code, "message": message, "source": source}


def _within(path: Path, root: Path) -> bool:
    try:
        path.resolve(strict=False).relative_to(root.resolve(strict=False))
        return True
    except ValueError:
        return False


def _split_modalities(value: str) -> set[str]:
    return {item.strip().lower() for item in re.split(r"[,;]", value or "") if item.strip()}


def _is_tier3a(row: Mapping[str, str]) -> bool:
    resolved = row.get("resolved_modality", "").strip().lower()
    return resolved in PAIR_MODALITIES or resolved.startswith("tier3a_")


def audit_catalog(provenance: Mapping[str, Any], configured_root: Path) -> tuple[dict[str, Any], list[dict[str, str]], list[dict[str, str]]]:
    blockers: list[dict[str, str]] = []
    source_record = provenance.get("source_catalog", {})
    source_path = Path(str(source_record.get("path", "")))
    if not source_path.is_file():
        raise Tier3ValidationError(f"frozen source catalog is unavailable: {source_path}")
    raw = source_path.read_bytes()
    text = raw.decode("utf-8")
    physical_lines = len(text.splitlines())
    with source_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        fields = reader.fieldnames or []
    if "Scientific Name" not in fields:
        raise Tier3ValidationError("frozen source catalog lacks Scientific Name")
    counts = Counter(row["Scientific Name"].strip() for row in rows)
    duplicates = [
        {"scientific_name": name, "multiplicity": multiplicity}
        for name, multiplicity in sorted(counts.items())
        if multiplicity > 1
    ]
    statistics = {
        "physical_lines": physical_lines,
        "header_lines": 1 if fields else 0,
        "data_rows": len(rows),
        "unique_species": len(counts),
        "data_row_excess_over_unique_species": len(rows) - len(counts),
        "duplicated_species": duplicates,
    }
    observed_sha = hashlib.sha256(raw).hexdigest()
    if statistics != EXPECTED_CATALOG_STATISTICS:
        blockers.append(_blocker("CATALOG_UNITS_OR_DUPLICATES_MISMATCH", "catalog lines, rows, unique species, excess, or duplicate multiplicities differ from the bounded gate", str(source_path)))
    if observed_sha != source_record.get("sha256"):
        blockers.append(_blocker("CATALOG_PROVENANCE_SHA256_MISMATCH", "catalog bytes do not match frozen provenance", str(source_path)))
    if source_record.get("line_count") != physical_lines:
        blockers.append(_blocker("CATALOG_PROVENANCE_LINE_COUNT_MISMATCH", "catalog physical-line count does not match frozen provenance", str(source_path)))
    if not HEX40.fullmatch(str(source_record.get("source_commit", ""))):
        blockers.append(_blocker("CATALOG_SOURCE_COMMIT_INVALID", "catalog provenance lacks an exact source commit", str(source_path)))
    if not _within(source_path, configured_root / "manifests"):
        blockers.append(_blocker("CATALOG_OUTSIDE_CONFIGURED_ROOT", "frozen catalog is not under the configured root manifests directory", str(source_path)))
    audit = {
        "path": str(source_path),
        "sha256": observed_sha,
        "size_bytes": len(raw),
        "source_url": source_record.get("source_url"),
        "source_commit": source_record.get("source_commit"),
        "retrieved_at_utc": source_record.get("retrieved_at_utc"),
        "statistics": statistics,
        "provenance_record_sha256": sha256_json(source_record),
    }
    return audit, blockers, rows


def audit_environment(
    channels_path: Path,
    manifest_path: Path,
    environment_path: Path,
) -> tuple[dict[str, Any], list[dict[str, str]]]:
    blockers: list[dict[str, str]] = []
    record = load_json(environment_path)
    channels_sha = sha256_file(channels_path)
    manifest_sha = sha256_file(manifest_path)
    commit = str(record.get("channel_commit", ""))
    channels_text = channels_path.read_text(encoding="utf-8")
    if not HEX40.fullmatch(commit) or commit not in channels_text:
        blockers.append(_blocker("GUIX_CHANNEL_NOT_EXACTLY_PINNED", "Guix environment channel commit is absent or not pinned by channels.scm", str(environment_path)))
    if record.get("channels_sha256") != channels_sha:
        blockers.append(_blocker("GUIX_CHANNELS_DIGEST_MISMATCH", "recorded Guix channels digest differs from the live pinned file", str(environment_path)))
    if record.get("manifest_sha256") != manifest_sha:
        blockers.append(_blocker("GUIX_MANIFEST_DIGEST_MISMATCH", "recorded Guix manifest digest differs from the live pinned file", str(environment_path)))
    if "guix time-machine" not in str(record.get("execution", "")) or "--pure" not in str(record.get("execution", "")):
        blockers.append(_blocker("GUIX_EXECUTION_NOT_PINNED_PURE", "environment record does not require pinned guix time-machine shell --pure", str(environment_path)))
    identity = {
        "channel_commit": commit,
        "channels_sha256": channels_sha,
        "manifest_sha256": manifest_sha,
        "environment_record_sha256": sha256_file(environment_path),
        "profile_store_path": record.get("profile_store_path"),
        "profile_derivation": next((item for item in record.get("derivations", []) if str(item).endswith("-profile.drv")), None),
        "resolved_channels_sha256": record.get("resolved_channels_sha256"),
        "execution": record.get("execution"),
        "tool_versions": record.get("tool_versions", {}),
    }
    identity["sha256"] = sha256_json(identity)
    return identity, blockers


def _parse_statfs(commands: Sequence[Mapping[str, Any]]) -> tuple[int | None, int | None]:
    pattern = re.compile(r"block_size=(\d+).*blocks_available=(\d+).*files_free=(\d+)")
    for command in commands:
        match = pattern.search(str(command.get("stdout", "")))
        if match:
            block_size, available_blocks, free_inodes = map(int, match.groups())
            return block_size * available_blocks, free_inodes
    return None, None


def audit_storage(
    root_config: Mapping[str, Any],
    validation: Mapping[str, Any],
    catalog_path: Path,
    required_bytes: int,
    required_inodes: int,
) -> tuple[dict[str, Any], list[dict[str, str]]]:
    blockers: list[dict[str, str]] = []
    root = Path(str(root_config.get("root", "")))
    if not root.is_absolute() or not root.is_dir():
        blockers.append(_blocker("ROOT_IDENTITY_INVALID", "configured external VGP root is not an existing absolute directory", str(root)))
        live_identity: dict[str, Any] = {"exists": False, "resolved_path": str(root.resolve(strict=False))}
    else:
        stat = root.stat()
        live_identity = {
            "exists": True,
            "resolved_path": str(root.resolve()),
            "device": stat.st_dev,
            "inode": stat.st_ino,
            "uid": stat.st_uid,
            "gid": stat.st_gid,
            "mode_octal": oct(stat.st_mode & 0o7777),
        }
    if validation.get("root") != str(root) or validation.get("root_owner", {}).get("path") != str(root):
        blockers.append(_blocker("ROOT_VALIDATION_IDENTITY_MISMATCH", "storage validation does not identify the configured external root", "analysis/vgp_data_root_validation.json"))
    owner = validation.get("root_owner", {})
    if owner.get("world_writable") is not False:
        blockers.append(_blocker("ROOT_WORLD_WRITABLE", "configured VGP root is or may be world writable", "analysis/vgp_data_root_validation.json"))
    if not _within(catalog_path, root):
        blockers.append(_blocker("CATALOG_ROOT_IDENTITY_MISMATCH", "catalog provenance is outside the configured external root", str(catalog_path)))
    required_smoke = ("file_fsync", "staging_dir_fsync", "atomic_promotion", "checksum_verification", "lock_behavior", "cleanup")
    smoke = validation.get("smoke_tests", {})
    for key in required_smoke:
        if smoke.get(key, {}).get("status") != "pass":
            code = "ROOT_ATOMIC_PROMOTION_UNSAFE" if key == "atomic_promotion" else "ROOT_SAFETY_EVIDENCE_FAILED"
            blockers.append(_blocker(code, f"external-root safety evidence {key} is not pass", "analysis/vgp_data_root_validation.json"))
    if smoke.get("same_filesystem") is not True:
        blockers.append(_blocker("ROOT_ATOMIC_FILESYSTEM_MISMATCH", "staging and immutable promotion targets are not proven on one filesystem", "analysis/vgp_data_root_validation.json"))

    system = validation.get("system_evidence", {})
    free_bytes, free_inodes = _parse_statfs(system.get("commands", []))
    headroom = 0.25
    bytes_with_headroom = math.ceil(required_bytes / (1.0 - headroom))
    inodes_with_headroom = math.ceil(required_inodes / (1.0 - headroom))
    filesystem_pass = free_bytes is not None and free_bytes >= bytes_with_headroom
    inode_pass = free_inodes is not None and free_inodes >= inodes_with_headroom
    if not filesystem_pass:
        blockers.append(_blocker("FILESYSTEM_FREE_SPACE_HEADROOM_INSUFFICIENT", "filesystem free bytes do not prove 25% worst-case headroom", "analysis/vgp_data_root_validation.json"))
    if not inode_pass:
        blockers.append(_blocker("FILESYSTEM_INODE_HEADROOM_INSUFFICIENT", "filesystem free inodes do not prove 25% worst-case headroom", "analysis/vgp_data_root_validation.json"))

    quota = system.get("quota_state", {})
    quota_status = str(quota.get("status", "missing")).lower()
    remaining_bytes = _integer(quota.get("remaining_bytes"))
    remaining_inodes = _integer(quota.get("remaining_inodes"))
    if remaining_bytes is None:
        limit = _integer(quota.get("limit_bytes"))
        used = _integer(quota.get("used_bytes"))
        if limit is not None and used is not None:
            remaining_bytes = limit - used
    if remaining_inodes is None:
        limit = _integer(quota.get("limit_inodes"))
        used = _integer(quota.get("used_inodes"))
        if limit is not None and used is not None:
            remaining_inodes = limit - used
    enforceable_known = quota_status in {"available", "reported"} and remaining_bytes is not None and remaining_inodes is not None
    allocation_pass = bool(enforceable_known and remaining_bytes >= bytes_with_headroom and remaining_inodes >= inodes_with_headroom)
    if not enforceable_known:
        blockers.append(_blocker("QUOTA_UNAVAILABLE", str(quota.get("message", "no exact enforceable allocation/quota evidence is recorded")), "analysis/vgp_data_root_validation.json"))
    elif not allocation_pass:
        blockers.append(_blocker("QUOTA_HEADROOM_INSUFFICIENT", "enforceable allocation has less than 25% worst-case byte or inode headroom", "analysis/vgp_data_root_validation.json"))

    audit = {
        "root": str(root),
        "live_identity": live_identity,
        "headroom_fraction_required": headroom,
        "worst_case": {
            "required_bytes": required_bytes,
            "required_inodes": required_inodes,
            "required_bytes_with_headroom": bytes_with_headroom,
            "required_inodes_with_headroom": inodes_with_headroom,
        },
        "filesystem": {
            "free_bytes": free_bytes,
            "free_inodes": free_inodes,
            "byte_headroom_pass": filesystem_pass,
            "inode_headroom_pass": inode_pass,
            "is_not_an_enforceable_user_allocation": True,
        },
        "enforceable_allocation": {
            "status": "known" if enforceable_known else "unknown",
            "source_status": quota_status,
            "available_interfaces": quota.get("available_interfaces", []),
            "remaining_bytes": remaining_bytes,
            "remaining_inodes": remaining_inodes,
            "headroom_pass": allocation_pass,
        },
        "safety_smoke_pass": not any(item["code"].startswith("ROOT_") for item in blockers),
        "adequate": not blockers,
        "policy": "filesystem capacity and enforceable per-user/allocation capacity are independent; every stricter limit wins",
    }
    audit["sha256"] = sha256_json(audit)
    return audit, blockers


def _cache_response_path(record: Mapping[str, Any]) -> Path:
    path = Path(str(record.get("response_path", "")))
    return path if path.is_absolute() else PROJECT_ROOT / path


def audit_cache(index: Mapping[str, Any]) -> tuple[dict[str, Mapping[str, Any]], list[dict[str, str]]]:
    blockers: list[dict[str, str]] = []
    by_key: dict[str, Mapping[str, Any]] = {}
    for record in index.get("responses", []):
        key = str(record.get("request_key", ""))
        if not HEX64.fullmatch(key) or key in by_key:
            blockers.append(_blocker("CACHE_REQUEST_KEY_INVALID", "official response has an invalid or duplicate immutable request key", "analysis/vgp_resolution_cache/index.json"))
            continue
        by_key[key] = record
        path = _cache_response_path(record)
        # HEAD responses intentionally cache an empty body (zero bytes); the
        # payload's finite positive size is carried by Content-Length and is
        # audited separately for every acquisition obligation.
        size = _integer(record.get("response_size_bytes"))
        if not path.is_file() or size is None or size < 0:
            blockers.append(_blocker("CACHE_RESPONSE_MISSING", f"immutable official response is absent or has no finite size: {key}", str(path)))
            continue
        if path.stat().st_size != size or sha256_file(path) != record.get("response_sha256"):
            blockers.append(_blocker("CACHE_RESPONSE_DIGEST_OR_SIZE_MISMATCH", f"immutable official response bytes disagree with its index: {key}", str(path)))
        request = record.get("request", {})
        if request.get("endpoint") != record.get("endpoint") or not str(record.get("endpoint", "")).startswith("https://"):
            blockers.append(_blocker("CACHE_OFFICIAL_ENDPOINT_INVALID", f"response endpoint is not exact HTTPS official provenance: {key}", str(path)))
        if not request.get("source") or not request.get("source_version") or not record.get("retrieved_at_utc"):
            blockers.append(_blocker("CACHE_RETRIEVAL_PROVENANCE_INCOMPLETE", f"response lacks source/version/retrieval time: {key}", str(path)))
    return by_key, blockers


def _parse_md5_catalog(path: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        parts = line.strip().split(None, 1)
        if len(parts) == 2 and HEX32.fullmatch(parts[0].lower()):
            result[parts[1].lstrip("./")] = parts[0].lower()
    return result


def audit_retrieval_row(row: Mapping[str, str], cache: Mapping[str, Mapping[str, Any]]) -> dict[str, Any]:
    issues: list[dict[str, str]] = []
    candidate = row.get("candidate_id", "")

    def issue(code: str, message: str) -> None:
        issues.append({"code": code, "message": message})

    try:
        obligations = json.loads(row.get("acquisition_obligations", ""))
    except (TypeError, json.JSONDecodeError):
        obligations = []
        issue("RETRIEVAL_OBLIGATIONS_INVALID", "acquisition obligations are not valid JSON")
    if not isinstance(obligations, list) or not obligations:
        obligations = []
        issue("RETRIEVAL_OBLIGATIONS_EMPTY", "no staged acquisition obligations are bound")

    request_keys = [key for key in row.get("metadata_cache_request_keys", "").split(";") if key]
    request_records = [cache[key] for key in request_keys if key in cache]
    if len(request_records) != len(request_keys) or not request_keys:
        issue("RETRIEVAL_CACHE_LINKAGE_MISSING", "row does not link exclusively to present immutable official responses")
    endpoint_records = {str(record.get("endpoint")): record for record in request_records}
    md5_catalogs: list[dict[str, str]] = []
    for record in request_records:
        if record.get("request", {}).get("parameters", {}).get("object_type") == "md5_catalog":
            md5_catalogs.append(_parse_md5_catalog(_cache_response_path(record)))

    expected_assets: dict[str, dict[str, Any]] = {
        row.get("h1_fasta_url", ""): {
            "accession": row.get("h1_accession_version", ""),
            "size": _integer(row.get("h1_fasta_compressed_bytes"), positive=True),
            "md5": row.get("h1_provider_md5", "").lower(),
            "role": "h1_fasta",
        },
        row.get("annotation_gff_url", ""): {
            "accession": row.get("h1_accession_version", ""),
            "size": _integer(row.get("annotation_gff_compressed_bytes"), positive=True),
            "md5": row.get("annotation_provider_md5", "").lower(),
            "role": "native_h1_annotation",
        },
    }
    if _is_tier3a(row):
        expected_assets[row.get("h2_fasta_url", "")] = {
            "accession": row.get("h2_accession_version", ""),
            "size": _integer(row.get("h2_fasta_compressed_bytes"), positive=True),
            "md5": row.get("h2_provider_md5", "").lower(),
            "role": "h2_fasta",
        }
    expected_assets.pop("", None)
    observed_urls: set[str] = set()
    normalized: list[dict[str, Any]] = []
    required_common_steps = {"stage_full_payload", "compute_local_sha256", "reverify_local_sha256_before_promotion", "atomic_promote_read_only"}
    for obligation in obligations:
        if not isinstance(obligation, dict):
            issue("RETRIEVAL_OBLIGATION_INVALID", "an acquisition obligation is not an object")
            continue
        url = str(obligation.get("url", ""))
        observed_urls.add(url)
        expected = expected_assets.get(url)
        size = _integer(obligation.get("expected_size_bytes"), positive=True)
        accession = str(obligation.get("accession_version", ""))
        steps = [str(step) for step in obligation.get("steps", [])]
        if expected is None:
            issue("RETRIEVAL_ASSET_NOT_DECLARED", f"obligation URL is not an exact row asset: {url}")
        else:
            if accession != expected["accession"] or accession not in url or not VERSIONED_ACCESSION.fullmatch(accession):
                issue("RETRIEVAL_ACCESSION_VERSION_MISMATCH", f"obligation does not preserve exact accession.version: {url}")
            if size is None or size != expected["size"]:
                issue("RETRIEVAL_SIZE_INVALID", f"obligation size is missing, non-finite, or disagrees with official metadata: {url}")
        if obligation.get("pre_download_eligible") is not True or obligation.get("remote_checksum_required_for_pre_download") is not False:
            issue("RETRIEVAL_STAGING_POLICY_INVALID", f"obligation does not permit exact-version staged local hashing semantics: {url}")
        if not required_common_steps.issubset(steps) or size is None or f"verify_expected_size_bytes:{size}" not in steps:
            issue("RETRIEVAL_LOCAL_SHA256_OBLIGATION_MISSING", f"size/local-SHA-256/reverification/atomic-promotion steps are incomplete: {url}")
        head = endpoint_records.get(url)
        official_size = None
        if head is not None:
            official_size = _integer(head.get("response_headers", {}).get("content-length"), positive=True)
            request = head.get("request", {})
            if request.get("method") != "HEAD" or request.get("source") != "NCBI genomes FTP" or request.get("source_version") != accession:
                issue("RETRIEVAL_OFFICIAL_VERSION_PROVENANCE_INVALID", f"asset lacks exact-version official HEAD provenance: {url}")
        if head is None or official_size != size:
            issue("RETRIEVAL_OFFICIAL_SIZE_LINKAGE_MISSING", f"asset size is not linked to the immutable official exact-version HEAD response: {url}")
        basename = url.rsplit("/", 1)[-1]
        official_md5 = next((catalog[basename] for catalog in md5_catalogs if basename in catalog), "")
        provider_md5 = expected["md5"] if expected else ""
        checksum_verified = False
        if official_md5:
            checksum_verified = provider_md5 == official_md5
            if not checksum_verified or obligation.get("official_checksum_verification_required") is not True or f"verify_official_md5:{official_md5}" not in steps:
                issue("RETRIEVAL_SOURCE_CHECKSUM_NOT_VERIFIED", f"existing official MD5 is not verified and required at acquisition: {url}")
        elif provider_md5 or obligation.get("official_checksum_verification_required") is not False:
            issue("RETRIEVAL_CHECKSUM_POLICY_INCONSISTENT", f"absent source checksum is not represented correctly: {url}")
        normalized.append(
            {
                "role": expected["role"] if expected else "unknown",
                "url": url,
                "accession_version": accession,
                "expected_size_bytes": size,
                "source_checksum": {"algorithm": "md5", "value": official_md5} if official_md5 else None,
                "source_checksum_verified_against_official_catalog": checksum_verified if official_md5 else None,
                "local_sha256_after_staging_required": "compute_local_sha256" in steps,
                "local_sha256_reverification_required": "reverify_local_sha256_before_promotion" in steps,
                "remote_checksum_required_for_pre_download": obligation.get("remote_checksum_required_for_pre_download"),
                "official_response_request_key": head.get("request_key") if head else None,
                "steps": steps,
            }
        )
    if observed_urls != set(expected_assets):
        issue("RETRIEVAL_ASSET_SET_MISMATCH", "acquisition obligations do not exactly cover required H1/native-annotation/Tier3A-H2 assets")
    return {
        "candidate_id": candidate,
        "pre_download_ready": not issues,
        "metadata_cache_request_keys": request_keys,
        "obligations": normalized,
        "issues": issues,
    }


def audit_manifest_row(row: Mapping[str, str], retrieval: Mapping[str, Any]) -> dict[str, Any]:
    issues: list[dict[str, str]] = list(retrieval.get("issues", []))

    def issue(code: str, message: str) -> None:
        issues.append({"code": code, "message": message})

    h1 = row.get("h1_accession_version", "")
    if not VERSIONED_ACCESSION.fullmatch(h1):
        issue("H1_ACCESSION_INVALID", "H1 accession is missing or not exact accession.version")
    if row.get("h1_exact_version_status") != "official_current_exact_version":
        issue("H1_EXACT_VERSION_NOT_OFFICIAL", "H1 exact-version status is not official")
    if row.get("taxon_identity_status") != "exact_taxid_match" or _integer(row.get("ncbi_taxid"), positive=True) is None:
        issue("TAXON_IDENTITY_NOT_EXACT", "H1 taxon identity is not exact")
    if row.get("annotation_reference_accession_version") != h1:
        issue("ANNOTATION_REFERENCE_MISMATCH", "native annotation reference does not exactly equal H1 accession.version")
    if row.get("annotation_native_status") != "official_ncbi_native_exact_h1":
        issue("ANNOTATION_NOT_NATIVE_H1", "annotation is not official native annotation for exact H1")
    if row.get("annotation_sequence_region_linkage_status") != "proven_official_exact_h1":
        issue("ANNOTATION_LINKAGE_NOT_PROVEN", "annotation sequence-region linkage is not proven for exact H1")
    if not row.get("annotation_accession_version") or not row.get("annotation_release_version_or_date"):
        issue("ANNOTATION_VERSION_MISSING", "annotation release identity is not exact and versioned")
    if row.get("callability_reference_accession_version") != h1:
        issue("CALLABILITY_REFERENCE_MISMATCH", "post-alignment measurement reference is not exact H1")

    tier3a = _is_tier3a(row)
    if tier3a:
        h2 = row.get("h2_accession_version", "")
        if not VERSIONED_ACCESSION.fullmatch(h2):
            issue("H2_ACCESSION_INVALID", "Tier3A H2 accession is missing or not exact accession.version")
        if row.get("same_individual_status", "").lower() != "yes":
            issue("PAIR_NOT_SAME_INDIVIDUAL", "Tier3A H1/H2 pair lacks affirmative same-individual evidence")
        if row.get("phase_evidence_status", "").lower() not in AFFIRMATIVE_PHASE_STATUSES:
            issue("PAIR_PHASE_EVIDENCE_MISSING", "Tier3A H1/H2 pair lacks affirmative phasing evidence")
        if not str(row.get("pair_evidence_url", "")).startswith("https://"):
            issue("PAIR_EVIDENCE_URL_MISSING", "Tier3A pair lacks an exact HTTPS evidence locator")
    exact_size = _integer(row.get("predicted_download_bytes_exact"), positive=True)
    obligation_sum = sum(item.get("expected_size_bytes") or 0 for item in retrieval.get("obligations", []))
    if exact_size is None or exact_size != obligation_sum:
        issue("DECLARED_DOWNLOAD_SIZE_INVALID", "finite exact row download bytes do not equal required asset sizes")
    for field in (
        "predicted_core_hours_high",
        "predicted_peak_memory_gib_high",
        "predicted_scratch_gb_high",
        "predicted_inode_count_high",
        "predicted_moosefs_read_gb_high",
        "predicted_moosefs_write_gb_high",
        "predicted_metadata_operations_high",
    ):
        if _finite(row.get(field), positive=True) is None:
            issue("ROW_RESOURCE_NOT_FINITE", f"row resource bound {field} is absent, non-positive, or non-finite")
    metadata_ready = not issues
    return {
        "candidate_id": row.get("candidate_id"),
        "scientific_name": row.get("scientific_name_source"),
        "catalog_row_number": _integer(row.get("catalog_row_number"), positive=True),
        "resolved_modality": row.get("resolved_modality"),
        "tier3a_required": tier3a,
        "metadata_ready": metadata_ready,
        "tier3a_ready": metadata_ready and tier3a,
        "observed_pilot_selected": row.get("pilot_selected") == "yes",
        "issues": issues,
    }


def audit_measurement_contract(index: Mapping[str, Any], manifest_rows: Sequence[Mapping[str, str]]) -> tuple[dict[str, Any], list[dict[str, str]]]:
    blockers: list[dict[str, str]] = []
    contract = copy_json(index.get("post_alignment_measurement_contract", {}))
    if contract.get("contract_id") != "vgp_post_alignment_denominators_v1":
        blockers.append(_blocker("MEASUREMENT_CONTRACT_ID_INVALID", "post-alignment denominator contract ID is not the approved version", "analysis/vgp_resolution_cache/index.json"))
    if contract.get("phase") != "post_alignment_pre_result_acceptance" or contract.get("pre_download_prerequisite") is not False:
        blockers.append(_blocker("MEASUREMENT_CONTRACT_PHASE_INVALID", "denominators must be measured after alignment, before result acceptance, and not before download", "analysis/vgp_resolution_cache/index.json"))
    if contract.get("minimum_thresholds") != EXPECTED_MEASUREMENT_THRESHOLDS:
        blockers.append(_blocker("MEASUREMENT_THRESHOLDS_INVALID", "denominator acceptance thresholds differ from the approved executable contract", "analysis/vgp_resolution_cache/index.json"))
    command = str(contract.get("executable_command", ""))
    if "resolve_vgp_candidates.py measure" not in command or "--metrics-json" not in command or "--output-json" not in command:
        blockers.append(_blocker("MEASUREMENT_COMMAND_NOT_EXECUTABLE", "denominator contract lacks its executable evaluator command", "analysis/vgp_resolution_cache/index.json"))
    if contract.get("failure_disposition") != "exclude_downstream_result":
        blockers.append(_blocker("MEASUREMENT_FAILURE_DISPOSITION_INVALID", "missing or inadequate denominators must exclude the affected downstream result", "analysis/vgp_resolution_cache/index.json"))
    for row in manifest_rows:
        if row.get("post_alignment_measurement_contract_id") != contract.get("contract_id") or row.get("post_alignment_result_disposition") != "not_measured_exclude_until_adequate":
            blockers.append(_blocker("ROW_MEASUREMENT_CONTRACT_MISMATCH", f"row {row.get('candidate_id')} does not bind the post-alignment exclusion contract", "analysis/vgp_pilot_manifest.tsv"))
    contract["sha256"] = sha256_json(contract)
    return contract, blockers


def copy_json(value: Any) -> Any:
    return json.loads(json.dumps(value))


def evaluate_post_alignment_measurements(contract: Mapping[str, Any], measurements: Mapping[str, Any]) -> dict[str, Any]:
    failures: list[str] = []
    thresholds = contract.get("minimum_thresholds", {})
    for field in ("callable_bases", "callable_fraction", "queryable_gene_count", "queryable_gene_bases"):
        observed = _finite(measurements.get(field), positive=True)
        minimum = _finite(thresholds.get(field), positive=True)
        if observed is None:
            failures.append(f"{field}_missing")
        elif minimum is None or observed < minimum:
            failures.append(f"{field}_below_minimum")
    return {
        "contract_id": contract.get("contract_id"),
        "accepted": not failures,
        "failed_thresholds": failures,
        "result_disposition": "accept_downstream_result" if not failures else "exclude_downstream_result",
        "measurements": dict(measurements),
    }


def _parse_integrated_caps(decision_rows: Sequence[Mapping[str, str]], plan_text: str, resource_rows: Sequence[Mapping[str, str]]) -> dict[str, list[dict[str, Any]]]:
    by_id = {row["decision_id"]: row for row in decision_rows}
    d011 = by_id["D011"]["resolution_hard_gate"]
    d012 = by_id["D012"]["resolution_hard_gate"]
    d018 = by_id["D018"]["resolution_hard_gate"]
    small = re.search(r"\*\*small cap:\*\*\s*([^|]+)", plan_text)
    if not small:
        raise Tier3ValidationError("integrated execution plan lacks the small cap clause")
    body = small.group(1)

    def number(pattern: str, text: str, label: str) -> float:
        match = re.search(pattern, text, flags=re.IGNORECASE)
        if not match:
            raise Tier3ValidationError(f"cannot parse integrated cap {label}")
        return float(match.group(1))

    high = next((row for row in resource_rows if row.get("stage_or_dataset") == "stratified_pilot" and row.get("scenario") == "high"), None)
    if high is None:
        raise Tier3ValidationError("integrated high resource envelope is missing")
    plan = {
        "species": number(r"([0-9.]+)\s+slots", body, "species"),
        "compressed_inputs_gib": number(r"([0-9.]+)\s+GB input/download", body, "inputs") * 1_000_000_000 / BYTES_PER_GIB,
        "scratch_gib": number(r"([0-9.]+)\s+GB scratch", body, "scratch") * 1_000_000_000 / BYTES_PER_GIB,
        "core_hours": number(r"([0-9.]+)\s+core-h", body, "core hours"),
        "concurrent_species": number(r"Tier3A\s+([0-9.]+)", body, "concurrency"),
        "memory_per_job_gib": number(r"([0-9.]+)\s+GiB/job", body, "memory"),
        "file_inodes": number(r"([0-9.]+)k\s+files", body, "files") * 1000,
        "moosefs_read_gb": number(r"([0-9.]+)\s+GB read", body, "read"),
        "moosefs_write_gb": number(r"/([0-9.]+)\s+GB write", body, "write"),
        "metadata_operations": number(r"([0-9.]+)k\s+ops", body, "metadata operations") * 1000,
        "peak_bandwidth_mib_s": number(r"([0-9.]+)\s+MiB/s", body, "bandwidth"),
        "aggregate_wall_hours": number(r"([0-9.]+)\s+h", body, "wall hours"),
        "persistent_input_gb": number(r"([0-9.]+)\s+GB input/download", body, "persistent input"),
        "persistent_output_gb": number(r"\+\s*([0-9.]+)\s+GB output", body, "persistent output"),
        "cpus_per_element": number(r"([0-9.]+)\s+CPU/element", body, "CPUs per element"),
    }
    d011_rows = re.search(r"exact\s+(eight|8)\s+rows", d011, flags=re.IGNORECASE)
    if not d011_rows:
        raise Tier3ValidationError("cannot parse integrated cap D011 rows")
    decisions = {
        "species": 8.0,
        "compressed_inputs_gib": number(r"<=([0-9.]+)\s+GB download/input", d011, "D011 inputs") * 1_000_000_000 / BYTES_PER_GIB,
        "scratch_gib": number(r"([0-9.]+)\s+GB scratch", d012, "D012 scratch") * 1_000_000_000 / BYTES_PER_GIB,
        "core_hours": number(r"<=([0-9.]+)\s+core-h", d012, "D012 core hours"),
        "concurrent_species": number(r"transfers <=([0-9.]+)", d018, "D018 transfers"),
        "memory_per_job_gib": number(r"([0-9.]+)\s+GiB/job", d012, "D012 memory"),
        "aggregate_wall_hours": number(r"([0-9.]+)\s+h", d012, "D012 wall hours"),
        "persistent_input_gb": number(r"<=([0-9.]+)\s+GB download/input", d011, "D011 persistent input"),
        "persistent_output_gb": number(r"([0-9.]+)\s+GB output", d012, "D012 persistent output"),
    }
    resources = {
        "species": float(high["species_count"]),
        "compressed_inputs_gib": float(high["persistent_input_gb"]) * 1_000_000_000 / BYTES_PER_GIB,
        "scratch_gib": float(high["local_scratch_peak_gb"]) * 1_000_000_000 / BYTES_PER_GIB,
        "core_hours": float(high["core_hours"]),
        "concurrent_species": float(high["concurrency_cap"]),
        "memory_per_job_gib": float(high["peak_resident_or_requested_memory_gib_per_element"]),
        "file_inodes": float(high["file_inode_count"]),
        "moosefs_read_gb": float(high["moosefs_read_gb"]),
        "moosefs_write_gb": float(high["moosefs_write_gb"]),
        "metadata_operations": float(high["metadata_operations"]),
        "peak_bandwidth_mib_s": float(high["peak_aggregate_bandwidth_mib_s"]),
        "aggregate_wall_hours": float(high["catalog_or_stage_wall_hours"]),
        "persistent_input_gb": float(high["persistent_input_gb"]),
        "persistent_output_gb": float(high["persistent_output_gb"]),
        "cpus_per_element": float(high["cpus_per_element"]),
    }
    sources: dict[str, list[dict[str, Any]]] = {}
    for name, (value, unit) in BASE_CAPS.items():
        sources.setdefault(name, []).append({"source": "regate_task_bound", "limit": value, "unit": unit})
    units = {
        "file_inodes": "count",
        "moosefs_read_gb": "GB",
        "moosefs_write_gb": "GB",
        "metadata_operations": "count",
        "peak_bandwidth_mib_s": "MiB/s",
        "aggregate_wall_hours": "hours",
        "persistent_input_gb": "GB",
        "persistent_output_gb": "GB",
        "cpus_per_element": "count",
    }
    for source_name, values in (("integrated_execution_plan_small_cap", plan), ("integrated_decision_hard_gates", decisions), ("integrated_resource_envelope_high", resources)):
        for name, value in values.items():
            unit = BASE_CAPS[name][1] if name in BASE_CAPS else units[name]
            sources.setdefault(name, []).append({"source": source_name, "limit": value, "unit": unit})
    return sources


def build_cap_vector(
    manifest_rows: Sequence[Mapping[str, str]],
    audits: Sequence[Mapping[str, Any]],
    budget_rows: Sequence[Mapping[str, str]],
    decision_rows: Sequence[Mapping[str, str]],
    plan_text: str,
    resource_rows: Sequence[Mapping[str, str]],
) -> tuple[dict[str, Any], list[dict[str, str]], dict[str, Any]]:
    blockers: list[dict[str, str]] = []
    ready_ids = {row["candidate_id"] for row in audits if row["metadata_ready"]}
    candidates = [row for row in budget_rows if row.get("row_type") == "candidate"]
    candidate_by_id = {row["candidate_id"]: row for row in candidates}
    budget_issues: list[dict[str, str]] = []
    ready_budget: list[Mapping[str, str]] = []
    for row in manifest_rows:
        budget = candidate_by_id.get(row["candidate_id"])
        if budget is None:
            budget_issues.append({"candidate_id": row["candidate_id"], "code": "SIZE_BUDGET_ROW_MISSING"})
            continue
        exact = _integer(budget.get("download_bytes_exact"), positive=True)
        if exact is None or exact != _integer(row.get("predicted_download_bytes_exact"), positive=True):
            budget_issues.append({"candidate_id": row["candidate_id"], "code": "SIZE_BUDGET_EXACT_BYTES_MISMATCH"})
        for field in ("compressed_inputs_gib", "core_hours_high", "peak_memory_gib_high", "scratch_gib_high"):
            if _finite(budget.get(field), positive=True) is None:
                budget_issues.append({"candidate_id": row["candidate_id"], "code": "SIZE_BUDGET_NONFINITE", "field": field})
        if row["candidate_id"] in ready_ids:
            ready_budget.append(budget)
    aggregate = next((row for row in budget_rows if row.get("row_type") == "aggregate_proposed_metadata_eligible"), None)
    if aggregate is None:
        budget_issues.append({"candidate_id": "aggregate", "code": "SIZE_BUDGET_PROPOSED_AGGREGATE_MISSING"})
    proposed = {
        "species": float(len(ready_budget)),
        "compressed_inputs_gib": sum(float(row["compressed_inputs_gib"]) for row in ready_budget),
        "scratch_gib": sum(sorted((float(row["scratch_gib_high"]) for row in ready_budget), reverse=True)[:2]),
        "core_hours": sum(float(row["core_hours_high"]) for row in ready_budget),
        "concurrent_species": float(min(2, len(ready_budget))),
        "memory_per_job_gib": max([0.0] + [float(row["peak_memory_gib_high"]) for row in ready_budget]),
        "file_inodes": sum(float(row["predicted_inode_count_high"]) for row in manifest_rows if row["candidate_id"] in ready_ids),
        "moosefs_read_gb": sum(float(row["predicted_moosefs_read_gb_high"]) for row in manifest_rows if row["candidate_id"] in ready_ids),
        "moosefs_write_gb": sum(float(row["predicted_moosefs_write_gb_high"]) for row in manifest_rows if row["candidate_id"] in ready_ids),
        "metadata_operations": sum(float(row["predicted_metadata_operations_high"]) for row in manifest_rows if row["candidate_id"] in ready_ids),
        "peak_bandwidth_mib_s": 120.0 if ready_budget else 0.0,
        "aggregate_wall_hours": sum(float(row["predicted_wall_hours_high"]) for row in manifest_rows if row["candidate_id"] in ready_ids),
        "persistent_input_gb": sum(int(row["download_bytes_exact"]) for row in ready_budget) / 1_000_000_000,
        "persistent_output_gb": sum(
            max(0, int(row["predicted_persistent_storage_bytes_exact"]) - int(row["predicted_download_bytes_exact"]))
            for row in manifest_rows if row["candidate_id"] in ready_ids
        ) / 1_000_000_000,
        "cpus_per_element": 8.0 if ready_budget else 0.0,
    }
    if aggregate is not None:
        comparisons = {
            "download_bytes_exact": sum(int(row["download_bytes_exact"]) for row in ready_budget),
            "compressed_inputs_gib": proposed["compressed_inputs_gib"],
            "core_hours_high": proposed["core_hours"],
            "peak_memory_gib_high": proposed["memory_per_job_gib"],
            "scratch_gib_high": proposed["scratch_gib"],
        }
        for field, expected in comparisons.items():
            observed = _finite(aggregate.get(field))
            if observed is None or not math.isclose(observed, float(expected), rel_tol=0, abs_tol=1e-6):
                budget_issues.append({"candidate_id": "aggregate", "code": "SIZE_BUDGET_PROPOSED_AGGREGATE_MISMATCH", "field": field})
    for issue in budget_issues:
        blockers.append(_blocker(issue["code"], f"size budget audit failed for {issue['candidate_id']}{'.' + issue.get('field', '') if issue.get('field') else ''}", "analysis/vgp_pilot_size_budget.tsv"))

    sources = _parse_integrated_caps(decision_rows, plan_text, resource_rows)
    dimensions: dict[str, Any] = {}
    for name, candidates_sources in sources.items():
        for source in candidates_sources:
            if _finite(source["limit"], positive=True) is None:
                blockers.append(_blocker("CAP_NOT_FINITE", f"cap {name} from {source['source']} is not finite and positive", "integrated cap sources"))
        winner = min(candidates_sources, key=lambda item: float(item["limit"]))
        observed = proposed[name]
        passed = math.isfinite(observed) and observed <= float(winner["limit"]) + 1e-9
        dimensions[name] = {
            "unit": winner["unit"],
            "limit": float(winner["limit"]),
            # Compatibility alias for prior consumers; `limit` is the
            # authoritative schema-2 name.
            "value": float(winner["limit"]),
            "observed": observed,
            "passes": passed,
            "winner_source": winner["source"],
            "sources": candidates_sources,
        }
        if not passed:
            blockers.append(_blocker(f"CAP_{name.upper()}_EXCEEDED", f"proposed finite worst-case {name}={observed} exceeds strictest limit {winner['limit']} {winner['unit']}", "strict cap vector"))
    vector = {
        "policy": "minimum finite limit wins across the task bound and every integrated numeric/storage/I/O cap",
        "operational_thresholds": {
            "quota_headroom_fraction_minimum": 0.25,
            "pause_at_cap_fraction": 0.80,
            "stop_resumption_at_cap_fraction": 0.95,
        },
        "dimensions": dimensions,
        "proposed_metadata_ready": proposed,
    }
    vector["sha256"] = sha256_json(vector)
    budget_audit = {
        "candidate_row_count": len(candidates),
        "metadata_ready_candidate_count": len(ready_budget),
        "issues": budget_issues,
        "aggregate_proposed_row": dict(aggregate) if aggregate else None,
    }
    return vector, blockers, budget_audit


def audit_dispositions(
    prior_rows: Sequence[Mapping[str, str]],
    manifest_rows: Sequence[Mapping[str, str]],
    rejection_rows: Sequence[Mapping[str, str]],
    catalog_rows: Sequence[Mapping[str, str]],
    row_audits: Sequence[Mapping[str, Any]],
    caps_pass: bool,
    storage_pass: bool,
) -> tuple[dict[str, Any], list[dict[str, str]], list[str]]:
    blockers: list[dict[str, str]] = []
    manifest = {row["candidate_id"]: row for row in manifest_rows}
    rejections = {row["candidate_id"]: row for row in rejection_rows}
    audits = {row["candidate_id"]: row for row in row_audits}
    if len(manifest) != len(manifest_rows) or len(rejections) != len(rejection_rows):
        blockers.append(_blocker("CANDIDATE_ID_DUPLICATED", "manifest or rejection ledger contains duplicate candidate IDs", "analysis/vgp_pilot_manifest.tsv; analysis/vgp_pilot_rejections.tsv"))
    rows: list[dict[str, Any]] = []
    for prior in prior_rows:
        candidate = prior["candidate_id"]
        catalog_line = _integer(prior.get("catalog_row_number"), positive=True)
        issues: list[str] = []
        if catalog_line is None or catalog_line < 2 or catalog_line > len(catalog_rows) + 1:
            issues.append("CATALOG_PHYSICAL_LINE_INVALID")
            catalog = {}
        else:
            catalog = catalog_rows[catalog_line - 2]
            if catalog.get("Scientific Name") != prior.get("scientific_name_source"):
                issues.append("CATALOG_SPECIES_MISMATCH")
        observed_rejection = rejections.get(candidate)
        if observed_rejection is None:
            issues.append("REJECTION_ROW_MISSING")
        prioritized = candidate in manifest
        expected_status = "blocked_stricter_cap" if prioritized else "rejected_not_prioritized"
        if observed_rejection is not None:
            if observed_rejection.get("acceptance_status") != expected_status:
                issues.append("REJECTION_STATUS_MISMATCH")
            if observed_rejection.get("catalog_row_number") != prior.get("catalog_row_number"):
                issues.append("REJECTION_CATALOG_LINE_MISMATCH")
            if observed_rejection.get("scientific_name") != prior.get("scientific_name_source"):
                issues.append("REJECTION_SPECIES_MISMATCH")
            if prioritized:
                current = manifest[candidate]
                for field in ("h1_accession_version", "h2_accession_version", "annotation_reference_accession_version", "resolved_modality"):
                    if observed_rejection.get(field, "") != current.get(field, ""):
                        issues.append("REJECTION_PRIORITIZED_METADATA_MISMATCH")
                        break
            elif observed_rejection.get("blocking_requirement_ids") != "PILOT_LIMIT":
                issues.append("REJECTION_NONPRIORITIZED_REASON_MISMATCH")
        rows.append(
            {
                "candidate_id": candidate,
                "catalog_physical_line": catalog_line,
                "scientific_name": prior.get("scientific_name_source"),
                "prioritized": prioritized,
                "expected_rejection_status": expected_status,
                "observed_rejection_status": observed_rejection.get("acceptance_status") if observed_rejection else None,
                "matches": not issues,
                "issues": sorted(set(issues)),
            }
        )
    prior_ids = {row["candidate_id"] for row in prior_rows}
    for extra in sorted(set(rejections) - prior_ids):
        rows.append({"candidate_id": extra, "matches": False, "issues": ["UNEXPECTED_REJECTION_ROW"]})
    if len(prior_rows) != 74 or len(prior_ids) != 74 or len(rejection_rows) != 74:
        blockers.append(_blocker("SEED_OR_REJECTION_ROW_COUNT_MISMATCH", "independent seed/rejection closure is not exactly 74 unique seed rows and 74 explicit disposition rows", "analysis/vgp_resolution_cache/prior_refusal; analysis/vgp_pilot_rejections.tsv"))
    for item in rows:
        if not item["matches"]:
            blockers.append(_blocker("REJECTION_DISPOSITION_MISMATCH", f"candidate {item['candidate_id']} disposition failed: {','.join(item['issues'])}", "analysis/vgp_pilot_rejections.tsv"))
    independently_selected = sorted(
        candidate for candidate, audit in audits.items()
        if audit["metadata_ready"] and caps_pass and storage_pass
    )
    observed_selected = sorted(row["candidate_id"] for row in manifest_rows if row.get("pilot_selected") == "yes")
    if observed_selected != independently_selected:
        blockers.append(_blocker("SELECTED_ROWS_DISAGREE_WITH_INDEPENDENT_GATE", "manifest pilot_selected rows differ from independently metadata/cap/storage-ready rows", "analysis/vgp_pilot_manifest.tsv"))
    audit = {
        "seed_row_count": len(prior_rows),
        "unique_seed_candidate_count": len(prior_ids),
        "rejection_row_count": len(rejection_rows),
        "all_rows_match": all(row["matches"] for row in rows),
        "rows": rows,
    }
    audit["sha256"] = sha256_json(audit)
    return audit, blockers, independently_selected


def _pair_evidence_payload(rows: Sequence[Mapping[str, str]]) -> dict[str, Any]:
    fields = (
        "candidate_id", "resolved_modality", "h1_accession_version", "h2_accession_version",
        "biosample_accession", "individual_or_isolate_id", "h1_h2_relationship",
        "haplotype_contig_map_sha256", "haplotype_contig_relationship_audit",
        "pair_evidence_url", "pair_evidence_retrieved_at_utc", "same_individual_evidence",
        "same_individual_status", "phase_evidence_status", "linked_h2_accessions_ncbi",
    )
    payload = {"rows": [{field: row.get(field, "") for field in fields} for row in rows]}
    payload["sha256"] = sha256_json(payload)
    return payload


def _input_record(path: Path) -> dict[str, Any]:
    return {"path": str(path), "sha256": sha256_file(path), "size_bytes": path.stat().st_size}


def _evaluate_gate(
    *,
    manifest_path: Path,
    rejections_path: Path,
    size_budget_path: Path,
    resolution_index_path: Path,
    prior_seed_manifest_path: Path,
    freeze_provenance_path: Path,
    root_config_path: Path,
    root_validation_path: Path,
    decisions_path: Path,
    execution_plan_path: Path,
    resource_budget_path: Path,
    guix_channels_path: Path,
    guix_manifest_path: Path,
    guix_environment_path: Path,
    generated_at_utc: str | None = None,
) -> dict[str, Any]:
    paths = {
        "manifest": Path(manifest_path), "rejections": Path(rejections_path), "size_budget": Path(size_budget_path),
        "resolution_index": Path(resolution_index_path), "prior_seed_manifest": Path(prior_seed_manifest_path),
        "freeze_provenance": Path(freeze_provenance_path), "root_config": Path(root_config_path),
        "root_validation": Path(root_validation_path), "decisions": Path(decisions_path),
        "execution_plan": Path(execution_plan_path), "resource_budget": Path(resource_budget_path),
        "guix_channels": Path(guix_channels_path), "guix_manifest": Path(guix_manifest_path),
        "guix_environment": Path(guix_environment_path),
    }
    inputs = {name: _input_record(path) for name, path in paths.items()}
    manifest_rows = load_tsv(paths["manifest"])
    rejection_rows = load_tsv(paths["rejections"])
    budget_rows = load_tsv(paths["size_budget"])
    prior_rows = load_tsv(paths["prior_seed_manifest"])
    index = load_json(paths["resolution_index"])
    provenance = load_json(paths["freeze_provenance"])
    root_config = load_json(paths["root_config"])
    root_validation = load_json(paths["root_validation"])
    decision_rows = load_tsv(paths["decisions"])
    resource_rows = load_tsv(paths["resource_budget"])

    blockers: list[dict[str, str]] = []
    catalog, found, catalog_rows = audit_catalog(provenance, Path(root_config.get("root", "")))
    blockers.extend(found)
    environment, found = audit_environment(paths["guix_channels"], paths["guix_manifest"], paths["guix_environment"])
    blockers.extend(found)
    cache, found = audit_cache(index)
    blockers.extend(found)
    retrieval_rows = [audit_retrieval_row(row, cache) for row in manifest_rows]
    retrieval = {
        "ready_count": sum(row["pre_download_ready"] for row in retrieval_rows),
        "rows": retrieval_rows,
        "policy": "official exact-version staging; finite size; verify every available source checksum; compute and reverify local SHA-256 before atomic read-only promotion",
    }
    retrieval["sha256"] = sha256_json(retrieval)
    row_audits = [audit_manifest_row(row, retrieval_row) for row, retrieval_row in zip(manifest_rows, retrieval_rows)]
    for audit in row_audits:
        for issue in audit["issues"]:
            blockers.append(_blocker(f"ROW_{issue['code']}", f"{audit['candidate_id']}: {issue['message']}", "analysis/vgp_pilot_manifest.tsv"))
    if len(manifest_rows) > 6:
        blockers.append(_blocker("MANIFEST_SPECIES_CAP_EXCEEDED", "repaired pilot manifest contains more than six candidates", "analysis/vgp_pilot_manifest.tsv"))
    measurement, found = audit_measurement_contract(index, manifest_rows)
    blockers.extend(found)
    cap_vector, found, budget_audit = build_cap_vector(
        manifest_rows, row_audits, budget_rows, decision_rows,
        paths["execution_plan"].read_text(encoding="utf-8"), resource_rows,
    )
    blockers.extend(found)
    cap_pass = all(item["passes"] for item in cap_vector["dimensions"].values())
    ready_ids = {audit["candidate_id"] for audit in row_audits if audit["metadata_ready"]}
    required_bytes = sum(_integer(row.get("predicted_persistent_storage_bytes_exact"), positive=True) or 0 for row in manifest_rows if row["candidate_id"] in ready_ids)
    required_inodes = int(sum(_finite(row.get("predicted_inode_count_high"), positive=True) or 0 for row in manifest_rows if row["candidate_id"] in ready_ids))
    storage, found = audit_storage(root_config, root_validation, Path(catalog["path"]), required_bytes, required_inodes)
    blockers.extend(found)
    dispositions, found, independently_selected = audit_dispositions(
        prior_rows, manifest_rows, rejection_rows, catalog_rows, row_audits, cap_pass, storage["adequate"]
    )
    blockers.extend(found)
    pair_evidence = _pair_evidence_payload(manifest_rows)
    selection = {
        "metadata_ready_candidate_ids": sorted(ready_ids),
        "independently_expected_selected_candidate_ids": independently_selected,
        "observed_selected_candidate_ids": sorted(row["candidate_id"] for row in manifest_rows if row.get("pilot_selected") == "yes"),
        "selection_rule": "exact metadata AND strict cap vector AND external-root/storage contract; repair eligibility labels are not inputs",
    }
    row_audit = {
        "summary": {
            "manifest_row_count": len(manifest_rows),
            "metadata_ready_count": sum(row["metadata_ready"] for row in row_audits),
            "tier3a_ready_count": sum(row["tier3a_ready"] for row in row_audits),
            "composition_ready_count": sum(row["metadata_ready"] for row in row_audits),
            "diversity_ready_count": sum(row["tier3a_ready"] for row in row_audits),
            "issue_counts": dict(sorted(Counter(issue["code"] for row in row_audits for issue in row["issues"]).items())),
        },
        "rows": row_audits,
    }
    catalog_component = {
        "catalog_audit": catalog,
        "freeze_source_record": provenance.get("source_catalog", {}),
    }
    root_component = {
        "root_config_sha256": inputs["root_config"]["sha256"],
        "root_validation_sha256": inputs["root_validation"]["sha256"],
        "storage_audit": storage,
    }
    boundary = {
        "catalog_provenance_digest": sha256_json(catalog_component),
        "data_root_storage_contract_digest": sha256_json(root_component),
        "environment_digest": environment["sha256"],
        "cap_vector_digest": cap_vector["sha256"],
        "retrieval_checksum_obligations_digest": retrieval["sha256"],
        "pair_evidence_digest": pair_evidence["sha256"],
        "measurement_contract_digest": measurement["sha256"],
        "row_dispositions_digest": dispositions["sha256"],
        "manifest_digest": inputs["manifest"]["sha256"],
        "root_contract_digest": sha256_json(root_component),
        "input_bundle_digest": sha256_json(inputs),
        "go_token": "GO",
    }
    boundary["authorization_tuple_digest"] = sha256_json(boundary)
    status = "GO" if not blockers else "NO_GO"
    gate: dict[str, Any] = {
        "schema_version": "2.0",
        "task_id": "regate-vgp-pilot",
        "generated_at_utc": generated_at_utc or utc_now(),
        "authorization_boundary": boundary,
        "authorization_scope": {
            "full_catalog_download_authorized": False,
            "raw_population_bulk_download_authorized": False,
            "demographic_inference_authorized": False,
            "biological_payloads_downloaded_by_gate": 0,
            "jobs_launched_by_gate": 0,
        },
        "inputs": inputs,
        "reproduction": {
            "source_catalog": {**catalog, "discrepancies": [] if catalog["statistics"] == EXPECTED_CATALOG_STATISTICS else [{"metric": "catalog_statistics", "observed": catalog["statistics"], "expected": EXPECTED_CATALOG_STATISTICS}]},
            "manifest_candidate_count": len(manifest_rows),
            "selected_row_count": len(selection["observed_selected_candidate_ids"]),
            "selected_candidate_ids": selection["observed_selected_candidate_ids"],
        },
        "quota_evidence": {
            "quota_status": storage["enforceable_allocation"]["status"],
            "free_bytes": storage["filesystem"]["free_bytes"],
            "free_inodes": storage["filesystem"]["free_inodes"],
        },
        "catalog_audit": catalog,
        "environment": environment,
        "storage_audit": storage,
        "cap_vector": cap_vector,
        "size_budget_audit": budget_audit,
        "retrieval_audit": retrieval,
        "pair_evidence": pair_evidence,
        "measurement_contract": measurement,
        "row_audit": row_audit,
        "selection_audit": selection,
        "disposition_audit": dispositions,
        "blockers": sorted(blockers, key=lambda item: (item["code"], item["message"], item["source"])),
        "decision": {
            "status": status,
            "only_literal_go_authorizes_downstream": True,
            "authorization_tuple_digest": boundary["authorization_tuple_digest"],
            "blocker_count": len(blockers),
        },
    }
    gate["decision_sha256"] = sha256_json(gate)
    return gate


def review_markdown(gate: Mapping[str, Any]) -> str:
    stats = gate["catalog_audit"]["statistics"]
    lines = [
        "# Independently regenerated VGP pilot gate review", "",
        f"- Decision: `{gate['decision']['status']}`", f"- Decision SHA-256: `{gate['decision_sha256']}`",
        f"- Authorization tuple SHA-256: `{gate['authorization_boundary']['authorization_tuple_digest']}`",
        "- Downstream rule: only the literal decision `GO`, with every bound digest reverified, authorizes acquisition or compute.", "",
        "## Catalog units and duplicate identities", "",
        f"- {stats['physical_lines']} physical lines = `{stats['header_lines']}` header + `{stats['data_rows']}` data rows.",
        f"- `{stats['unique_species']}` unique species; the data-row excess over unique species is `{stats['data_row_excess_over_unique_species']}`.",
    ]
    for item in stats["duplicated_species"]:
        lines.append(f"- `{item['scientific_name']}` multiplicity: `{item['multiplicity']}`.")
    lines.extend([
        "", "## Independent row closure", "",
        f"- Manifest rows audited: `{gate['row_audit']['summary']['manifest_row_count']}`.",
        f"- Exact H1/native-annotation metadata-ready rows: `{gate['row_audit']['summary']['metadata_ready_count']}`.",
        f"- Exact paired Tier3A-ready rows: `{gate['row_audit']['summary']['tier3a_ready_count']}`.",
        f"- Seed/rejection rows independently closed: `{gate['disposition_audit']['seed_row_count']}` / `{gate['disposition_audit']['rejection_row_count']}`; all match: `{str(gate['disposition_audit']['all_rows_match']).lower()}`.",
        f"- Independently expected selected rows: `{', '.join(gate['selection_audit']['independently_expected_selected_candidate_ids']) or 'none'}`.",
        "", "## Strictest finite cap vector", "",
    ])
    for name, payload in gate["cap_vector"]["dimensions"].items():
        lines.append(f"- `{name}`: observed `{payload['observed']}`; limit `{payload['limit']}` {payload['unit']}; pass `{str(payload['passes']).lower()}`; winner `{payload['winner_source']}`.")
    storage = gate["storage_audit"]
    lines.extend([
        "", "## External-root and storage contract", "",
        f"- Root: `{storage['root']}`; live identity digest is bound in the decision.",
        f"- Filesystem free bytes/inodes: `{storage['filesystem']['free_bytes']}` / `{storage['filesystem']['free_inodes']}`.",
        f"- Enforceable allocation/quota status: `{storage['enforceable_allocation']['status']}`.",
        "- Filesystem free space and inodes are reported independently from enforceable per-user/allocation limits. Unknown quota never counts as adequate and cannot be overridden by free space.",
        "- Required worst-case byte and inode capacity includes at least 25% headroom; every stricter integrated limit wins.",
        "", "## Retrieval, checksums, pairs, and denominators", "",
        f"- Pre-download-ready exact-version rows: `{gate['retrieval_audit']['ready_count']}`.",
        "- Each exact official payload is staged, finite-size checked, checked against every available official checksum, locally SHA-256 hashed, reverified, and atomically promoted read-only. A missing remote SHA-256/MD5 is not itself a pre-download blocker.",
        f"- Pair evidence digest: `{gate['authorization_boundary']['pair_evidence_digest']}`; Tier3A requires exact versioned H2 plus affirmative same-individual and phasing evidence.",
        f"- Denominator contract: `{gate['measurement_contract']['contract_id']}`; pre-download prerequisite `{str(gate['measurement_contract']['pre_download_prerequisite']).lower()}`.",
        "- Callable/queryable denominators are measured after alignment. Missing or sub-threshold measurements exclude the affected downstream result.",
        "", "## Bound decision components", "",
    ])
    for name, digest in gate["authorization_boundary"].items():
        if name.endswith("digest") or name.endswith("_digest"):
            lines.append(f"- `{name}`: `{digest}`")
    lines.extend(["", "## Blockers", ""])
    if gate["blockers"]:
        for item in gate["blockers"]:
            lines.append(f"- `{item['code']}`: {item['message']} ({item['source']})")
    else:
        lines.append("- none")
    lines.extend([
        "", "## Authorization exclusions", "",
        "- Full-catalog acquisition is unauthorized.", "- Raw population bulk download is unauthorized.",
        "- Demographic inference is unauthorized.", "- This gate launched zero downloads and zero jobs.", "",
    ])
    return "\n".join(lines)


def build_gate(
    *,
    manifest_path: Path = DEFAULT_MANIFEST,
    rejections_path: Path = DEFAULT_REJECTIONS,
    size_budget_path: Path = DEFAULT_SIZE_BUDGET,
    resolution_index_path: Path = DEFAULT_RESOLUTION_INDEX,
    prior_seed_manifest_path: Path = DEFAULT_PRIOR_SEEDS,
    freeze_provenance_path: Path = DEFAULT_FREEZE_PROVENANCE,
    root_config_path: Path = DEFAULT_ROOT_CONFIG,
    root_validation_path: Path = DEFAULT_ROOT_VALIDATION,
    decisions_path: Path = DEFAULT_DECISIONS,
    execution_plan_path: Path = DEFAULT_EXECUTION_PLAN,
    resource_budget_path: Path = DEFAULT_RESOURCE_BUDGET,
    guix_channels_path: Path = DEFAULT_GUIX_CHANNELS,
    guix_manifest_path: Path = DEFAULT_GUIX_MANIFEST,
    guix_environment_path: Path = DEFAULT_GUIX_ENVIRONMENT,
    gate_out: Path | None = DEFAULT_GATE_JSON,
    review_out: Path | None = DEFAULT_GATE_REVIEW,
) -> dict[str, Any]:
    gate = _evaluate_gate(
        manifest_path=manifest_path, rejections_path=rejections_path, size_budget_path=size_budget_path,
        resolution_index_path=resolution_index_path, prior_seed_manifest_path=prior_seed_manifest_path,
        freeze_provenance_path=freeze_provenance_path, root_config_path=root_config_path,
        root_validation_path=root_validation_path, decisions_path=decisions_path,
        execution_plan_path=execution_plan_path, resource_budget_path=resource_budget_path,
        guix_channels_path=guix_channels_path, guix_manifest_path=guix_manifest_path,
        guix_environment_path=guix_environment_path,
    )
    if gate_out is not None:
        write_json(gate_out, gate)
    if review_out is not None:
        Path(review_out).parent.mkdir(parents=True, exist_ok=True)
        Path(review_out).write_text(review_markdown(gate), encoding="utf-8")
    return gate


def load_gate(path: Path) -> dict[str, Any]:
    payload = load_json(path)
    expected = sha256_json({key: value for key, value in payload.items() if key != "decision_sha256"})
    if payload.get("decision_sha256") != expected:
        raise Tier3ValidationError("gate decision hash does not match the gate payload")
    if payload.get("schema_version") != "2.0":
        raise Tier3ValidationError("unsupported VGP gate schema")
    status = payload.get("decision", {}).get("status")
    if status not in {"GO", "NO_GO"}:
        raise Tier3ValidationError(f"gate decision state {status!r} is not the literal GO token")
    return payload


def authorize_gate_action(
    gate_path: Path,
    manifest_path: Path,
    root_config_path: Path,
    action: str,
    *,
    rejections_path: Path | None = None,
    size_budget_path: Path | None = None,
    resolution_index_path: Path | None = None,
    prior_seed_manifest_path: Path | None = None,
    freeze_provenance_path: Path | None = None,
    root_validation_path: Path | None = None,
    decisions_path: Path | None = None,
    execution_plan_path: Path | None = None,
    resource_budget_path: Path | None = None,
    guix_channels_path: Path | None = None,
    guix_manifest_path: Path | None = None,
    guix_environment_path: Path | None = None,
) -> dict[str, Any]:
    if action not in {"acquire", "compute"}:
        raise Tier3ValidationError(f"unsupported gate action: {action}")
    gate = load_gate(gate_path)
    current = _evaluate_gate(
        manifest_path=Path(manifest_path), rejections_path=rejections_path or DEFAULT_REJECTIONS,
        size_budget_path=size_budget_path or DEFAULT_SIZE_BUDGET,
        resolution_index_path=resolution_index_path or DEFAULT_RESOLUTION_INDEX,
        prior_seed_manifest_path=prior_seed_manifest_path or DEFAULT_PRIOR_SEEDS,
        freeze_provenance_path=freeze_provenance_path or DEFAULT_FREEZE_PROVENANCE,
        root_config_path=Path(root_config_path), root_validation_path=root_validation_path or DEFAULT_ROOT_VALIDATION,
        decisions_path=decisions_path or DEFAULT_DECISIONS,
        execution_plan_path=execution_plan_path or DEFAULT_EXECUTION_PLAN,
        resource_budget_path=resource_budget_path or DEFAULT_RESOURCE_BUDGET,
        guix_channels_path=guix_channels_path or DEFAULT_GUIX_CHANNELS,
        guix_manifest_path=guix_manifest_path or DEFAULT_GUIX_MANIFEST,
        guix_environment_path=guix_environment_path or DEFAULT_GUIX_ENVIRONMENT,
    )
    expected = gate["authorization_boundary"]
    observed = current["authorization_boundary"]
    comparisons = (
        ("catalog_provenance_digest", "catalog provenance digest mismatch"),
        ("data_root_storage_contract_digest", "root/storage contract digest mismatch"),
        ("environment_digest", "environment digest mismatch"),
        ("cap_vector_digest", "cap vector digest mismatch"),
        ("retrieval_checksum_obligations_digest", "retrieval/checksum obligations digest mismatch"),
        ("pair_evidence_digest", "pair evidence digest mismatch"),
        ("measurement_contract_digest", "measurement contract digest mismatch"),
        ("row_dispositions_digest", "row dispositions digest mismatch"),
        ("manifest_digest", "manifest digest mismatch"),
        ("input_bundle_digest", "bound input bundle digest mismatch"),
        ("authorization_tuple_digest", "authorization tuple digest mismatch"),
    )
    for key, message in comparisons:
        if observed.get(key) != expected.get(key):
            raise Tier3ValidationError(f"{message}: expected {expected.get(key)}, observed {observed.get(key)}")
    if gate["decision"]["status"] != "GO":
        raise Tier3ValidationError(f"gate decision is {gate['decision']['status']}; {action} is not authorized")
    if current["decision"]["status"] != "GO":
        raise Tier3ValidationError(f"live recomputed gate decision is {current['decision']['status']}; {action} is not authorized")
    return {
        "authorized": True, "action": action, "decision_sha256": gate["decision_sha256"],
        "authorization_tuple_digest": expected["authorization_tuple_digest"],
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    build = subparsers.add_parser("build", help="regenerate the offline signed/hashed gate")
    build.add_argument("--gate-out", type=Path, default=DEFAULT_GATE_JSON)
    build.add_argument("--review-out", type=Path, default=DEFAULT_GATE_REVIEW)
    authorize = subparsers.add_parser("authorize", help="recompute and enforce all bound gate digests")
    authorize.add_argument("--gate", type=Path, default=DEFAULT_GATE_JSON)
    authorize.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    authorize.add_argument("--root-config", type=Path, default=DEFAULT_ROOT_CONFIG)
    authorize.add_argument("--action", choices=("acquire", "compute"), required=True)
    args = parser.parse_args(argv)
    if args.command == "build":
        gate = build_gate(gate_out=args.gate_out, review_out=args.review_out)
        print(json.dumps({"decision": gate["decision"]["status"], "decision_sha256": gate["decision_sha256"], "blocker_count": len(gate["blockers"])}, sort_keys=True))
        return 0
    authorize_gate_action(args.gate, args.manifest, args.root_config, args.action)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
