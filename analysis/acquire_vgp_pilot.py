#!/usr/bin/env python3
"""Gate-bound, fail-closed acquisition for the repaired bounded VGP pilot.

The authorization check is deliberately completed before a downloader is
called or a biological staging file is opened.  The currently committed gate
is ``NO_GO`` and therefore this entrypoint only writes small refusal evidence.
The GO path is nevertheless implemented so that a later, independently
regenerated decision cannot weaken staging, checksum, cap, or promotion rules.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import stat
import sys
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Mapping, MutableMapping, Sequence

if __package__ in {None, ""}:  # Support ``python3 analysis/acquire_vgp_pilot.py``.
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from analysis import gate_vgp_pilot as gate
from analysis.tier3_common import Tier3ValidationError, sha256_file


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_GATE = PROJECT_ROOT / "analysis" / "vgp_pilot_gate.json"
DEFAULT_MANIFEST = PROJECT_ROOT / "analysis" / "vgp_pilot_manifest.tsv"
DEFAULT_ROOT_CONFIG = PROJECT_ROOT / "analysis" / "vgp_data_root_config.json"
DEFAULT_OUTPUT_MANIFEST = PROJECT_ROOT / "analysis" / "vgp_pilot_acquisition_manifest.tsv"
DEFAULT_OUTPUT_REPORT = PROJECT_ROOT / "analysis" / "vgp_pilot_acquisition_report.md"
DEFAULT_OUTPUT_INVENTORY = PROJECT_ROOT / "analysis" / "vgp_pilot_immutable_object_inventory.tsv"
DEFAULT_REFUSAL_EVIDENCE = PROJECT_ROOT / "analysis" / "vgp_pilot_acquisition_refusal.json"

GIB = 1024**3
VERSIONED_ACCESSION = re.compile(r"^GC[AF]_\d+\.\d+$")
HEX32 = re.compile(r"^[0-9a-f]{32}$")
HEX64 = re.compile(r"^[0-9a-f]{64}$")
ALLOWED_ASSET_ROLES = {
    "h1_fasta",
    "h2_fasta",
    "native_h1_annotation",
    "assembly_report",
    "sequence_report",
    "required_index",
    "required_metadata",
}
BIOLOGICAL_ROLES = {"h1_fasta", "h2_fasta", "native_h1_annotation"}
HARD_CAP_CEILINGS = {
    "species": 6.0,
    "compressed_inputs_gib": 120.0,
    "scratch_gib": 750.0,
    "core_hours": 1500.0,
    "concurrent_species": 2.0,
    "memory_per_job_gib": 256.0,
}
MANIFEST_FIELDS = [
    "run_id",
    "generated_at_utc",
    "started_at_utc",
    "completed_at_utc",
    "record_type",
    "candidate_id",
    "asset_role",
    "accession_version",
    "source_url",
    "expected_bytes",
    "observed_bytes",
    "cumulative_transferred_bytes",
    "expected_sha256",
    "observed_sha256",
    "source_checksum_algorithm",
    "source_checksum_value",
    "source_checksum_verified",
    "provider_md5",
    "http_status",
    "response_headers_json",
    "retries",
    "validation_outcomes_json",
    "failure_code",
    "failure_source",
    "failure_message",
    "staging_path",
    "quarantine_path",
    "promoted_path",
    "status",
]
INVENTORY_FIELDS = [
    "run_id",
    "promoted_at_utc",
    "candidate_id",
    "asset_role",
    "accession_version",
    "source_url",
    "bytes",
    "local_sha256",
    "source_checksum_algorithm",
    "source_checksum_value",
    "source_checksum_verified",
    "object_path",
    "mode_octal",
    "authorization_tuple_digest",
    "decision_sha256",
]
SLURM_ENV_VARS = ("SLURM_JOB_ID", "SLURM_ARRAY_JOB_ID", "SLURM_ARRAY_TASK_ID", "SLURM_CLUSTER_NAME")

DownloadFunction = Callable[[Mapping[str, Any], Path, MutableMapping[str, Any]], Mapping[str, Any]]
PromoteHook = Callable[[Path, str], None]
StagedValidator = Callable[[str, Sequence[Mapping[str, Any]]], Sequence[str]]


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def make_run_id() -> str:
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    return f"vgp-pilot-acquire-{stamp}-{os.getpid()}"


def _write_tsv(path: Path, fields: Sequence[str], rows: Sequence[Mapping[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.partial-{os.getpid()}")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", extrasaction="ignore", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def write_manifest(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    _write_tsv(path, MANIFEST_FIELDS, rows)


def write_inventory(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    _write_tsv(path, INVENTORY_FIELDS, rows)


def _write_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.partial-{os.getpid()}")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def _write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.partial-{os.getpid()}")
    temporary.write_text(text, encoding="utf-8")
    os.replace(temporary, path)


def detected_slurm_environment() -> dict[str, str]:
    return {name: os.environ[name] for name in SLURM_ENV_VARS if name in os.environ}


def _safe_gate_status(payload: Mapping[str, Any] | None) -> str:
    if not isinstance(payload, Mapping):
        return "MISSING_OR_UNREADABLE"
    decision = payload.get("decision")
    return str(decision.get("status", "UNKNOWN")) if isinstance(decision, Mapping) else "UNKNOWN"


def classify_failure_code(error: Exception) -> str:
    text = str(error).lower()
    if "gate decision hash does not match" in text:
        return "GATE_TAMPERED"
    if "gate decision is no_go" in text:
        return "GATE_NO_GO"
    if "not the literal go token" in text or "unknown gate decision" in text:
        return "GATE_DECISION_UNKNOWN"
    if "manifest digest mismatch" in text:
        return "MANIFEST_DIGEST_MISMATCH"
    if "root/storage contract digest mismatch" in text or "root contract digest mismatch" in text:
        return "ROOT_STORAGE_CONTRACT_MISMATCH"
    if "cap vector" in text:
        return "CAP_VECTOR_MISMATCH"
    if "retrieval" in text or "approved url" in text or "accession.version" in text:
        return "RETRIEVAL_CONTRACT_MISMATCH"
    if "checksum" in text or "digest" in text or "gate decision hash" in text:
        return "BOUND_DIGEST_MISMATCH"
    if "does not exist" in text or "no such file" in text or "unavailable" in text:
        return "BOUND_INPUT_MISSING"
    return "PRECHECK_FAILED"


def _canonical_digest(payload: Mapping[str, Any], digest_key: str = "sha256") -> str:
    return gate.sha256_json({key: value for key, value in payload.items() if key != digest_key})


def verify_gate_internal_bindings(payload: Mapping[str, Any]) -> None:
    """Reject a self-rehashed gate whose bound nested objects were altered."""

    boundary = payload.get("authorization_boundary")
    if not isinstance(boundary, Mapping):
        raise Tier3ValidationError("gate authorization boundary is missing")
    decision = payload.get("decision")
    if not isinstance(decision, Mapping):
        raise Tier3ValidationError("unknown gate decision: decision object is missing")
    status = decision.get("status")
    if status not in {"GO", "NO_GO"}:
        raise Tier3ValidationError(f"unknown gate decision {status!r}; decision is not the literal GO token")
    if boundary.get("go_token") != "GO":
        raise Tier3ValidationError("bound literal GO token was altered")
    tuple_observed = gate.sha256_json(
        {key: value for key, value in boundary.items() if key != "authorization_tuple_digest"}
    )
    if boundary.get("authorization_tuple_digest") != tuple_observed:
        raise Tier3ValidationError("authorization tuple digest mismatch inside gate")
    if decision.get("authorization_tuple_digest") != tuple_observed:
        raise Tier3ValidationError("decision authorization tuple digest mismatch")

    nested = (
        ("cap_vector", "cap_vector_digest", "cap vector"),
        ("retrieval_audit", "retrieval_checksum_obligations_digest", "retrieval/checksum obligations"),
        ("pair_evidence", "pair_evidence_digest", "pair evidence"),
        ("measurement_contract", "measurement_contract_digest", "measurement contract"),
        ("disposition_audit", "row_dispositions_digest", "row dispositions"),
    )
    for object_key, boundary_key, label in nested:
        obj = payload.get(object_key)
        if not isinstance(obj, Mapping):
            raise Tier3ValidationError(f"{label} object is missing")
        observed = _canonical_digest(obj)
        if obj.get("sha256") != observed or boundary.get(boundary_key) != observed:
            raise Tier3ValidationError(f"{label} digest mismatch inside gate")

    environment = payload.get("environment")
    if not isinstance(environment, Mapping) or environment.get("sha256") != _canonical_digest(environment):
        raise Tier3ValidationError("environment digest mismatch inside gate")
    if boundary.get("environment_digest") != environment.get("sha256"):
        raise Tier3ValidationError("bound environment digest mismatch")
    inputs = payload.get("inputs")
    if not isinstance(inputs, Mapping) or boundary.get("input_bundle_digest") != gate.sha256_json(inputs):
        raise Tier3ValidationError("bound input bundle digest mismatch")
    manifest_input = inputs.get("manifest")
    if not isinstance(manifest_input, Mapping) or boundary.get("manifest_digest") != manifest_input.get("sha256"):
        raise Tier3ValidationError("manifest digest mismatch inside gate")
    storage = payload.get("storage_audit")
    root_input = inputs.get("root_config")
    validation_input = inputs.get("root_validation")
    if not all(isinstance(item, Mapping) for item in (storage, root_input, validation_input)):
        raise Tier3ValidationError("root/storage contract components are missing")
    root_component = {
        "root_config_sha256": root_input["sha256"],
        "root_validation_sha256": validation_input["sha256"],
        "storage_audit": storage,
    }
    root_digest = gate.sha256_json(root_component)
    if boundary.get("root_contract_digest") != root_digest or boundary.get("data_root_storage_contract_digest") != root_digest:
        raise Tier3ValidationError("root/storage contract digest mismatch inside gate")


def verify_cap_contract(payload: Mapping[str, Any]) -> None:
    cap_vector = payload["cap_vector"]
    dimensions = cap_vector.get("dimensions", {})
    for name, hard_ceiling in HARD_CAP_CEILINGS.items():
        dimension = dimensions.get(name)
        if not isinstance(dimension, Mapping):
            raise Tier3ValidationError(f"cap vector is missing {name}")
        limit = float(dimension.get("limit"))
        source_limits = [float(source["limit"]) for source in dimension.get("sources", [])]
        if not source_limits or limit != min(source_limits):
            raise Tier3ValidationError(f"cap vector relaxed {name}; strictest integrated limit was not retained")
        if limit > hard_ceiling:
            raise Tier3ValidationError(f"cap vector relaxed {name} above hard ceiling {hard_ceiling}")
    thresholds = cap_vector.get("operational_thresholds", {})
    if float(thresholds.get("quota_headroom_fraction_minimum", -1)) < 0.25:
        raise Tier3ValidationError("cap vector relaxed required storage headroom below 25 percent")


def _load_untrusted_gate(path: Path) -> Mapping[str, Any] | None:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    return value if isinstance(value, Mapping) else None


def _selected_assets(payload: Mapping[str, Any]) -> list[dict[str, Any]]:
    selection = payload["selection_audit"]
    expected = sorted(selection.get("independently_expected_selected_candidate_ids", []))
    observed = sorted(selection.get("observed_selected_candidate_ids", []))
    reproduced = sorted(payload.get("reproduction", {}).get("selected_candidate_ids", []))
    if not expected or observed != expected or reproduced != expected:
        raise Tier3ValidationError("GO gate selected-row sets are empty or disagree")
    if len(expected) > 6:
        raise Tier3ValidationError("selected species exceed cap vector")
    row_audits = {row["candidate_id"]: row for row in payload["row_audit"]["rows"]}
    retrieval = {row["candidate_id"]: row for row in payload["retrieval_audit"]["rows"]}
    assets: list[dict[str, Any]] = []
    for candidate_id in expected:
        row = row_audits.get(candidate_id)
        retrieval_row = retrieval.get(candidate_id)
        if not row or not retrieval_row or not row.get("metadata_ready") or not retrieval_row.get("pre_download_ready"):
            raise Tier3ValidationError(f"selected row is not exact-reference retrieval-ready: {candidate_id}")
        obligations = retrieval_row.get("obligations", [])
        roles = {item.get("role") for item in obligations}
        required = {"h1_fasta", "native_h1_annotation"}
        if row.get("tier3a_required"):
            required.add("h2_fasta")
        if not required.issubset(roles):
            raise Tier3ValidationError(f"selected row lacks required exact assets: {candidate_id}")
        for obligation in obligations:
            asset = dict(obligation)
            role = str(asset.get("role", ""))
            accession = str(asset.get("accession_version", ""))
            url = str(asset.get("url", ""))
            size = asset.get("expected_size_bytes")
            if role not in ALLOWED_ASSET_ROLES:
                raise Tier3ValidationError(f"unapproved asset role {role!r}: {candidate_id}")
            if not VERSIONED_ACCESSION.fullmatch(accession) or accession not in url:
                raise Tier3ValidationError(f"approved URL/accession.version mismatch: {candidate_id} {url}")
            if not isinstance(size, int) or size <= 0:
                raise Tier3ValidationError(f"approved asset has no finite positive size: {url}")
            if asset.get("local_sha256_after_staging_required") is not True or asset.get("local_sha256_reverification_required") is not True:
                raise Tier3ValidationError(f"approved asset weakens local SHA-256 obligations: {url}")
            asset["candidate_id"] = candidate_id
            assets.append(asset)
    if sum(int(asset["expected_size_bytes"]) for asset in assets) > int(
        float(payload["cap_vector"]["dimensions"]["compressed_inputs_gib"]["limit"]) * GIB
    ):
        raise Tier3ValidationError("approved asset bytes exceed strict compressed-input cap")
    return assets


def preflight_authorize(
    gate_path: Path, manifest_path: Path, root_config_path: Path
) -> tuple[dict[str, Any], dict[str, Any], list[dict[str, Any]]]:
    payload = gate.load_gate(gate_path)
    verify_gate_internal_bindings(payload)
    verify_cap_contract(payload)
    # NO_GO is independently sufficient to refuse.  In particular, do not
    # touch the external root or evaluate anything that could reach a provider
    # after this literal decision boundary.  A GO must additionally survive
    # the complete live recomputation below.
    if payload["decision"]["status"] != "GO":
        raise Tier3ValidationError(f"gate decision is {payload['decision']['status']}; acquire is not authorized")
    authorization = gate.authorize_gate_action(gate_path, manifest_path, root_config_path, "acquire")
    if payload["authorization_scope"].get("full_catalog_download_authorized") is not False:
        raise Tier3ValidationError("full catalog acquisition scope was altered")
    if payload["authorization_scope"].get("raw_population_bulk_download_authorized") is not False:
        raise Tier3ValidationError("raw population bulk acquisition scope was altered")
    if payload["authorization_scope"].get("demographic_inference_authorized") is not False:
        raise Tier3ValidationError("demographic inference scope was altered")
    return payload, authorization, _selected_assets(payload)


def _default_downloader(
    asset: Mapping[str, Any], part_path: Path, transfer_state: MutableMapping[str, Any]
) -> Mapping[str, Any]:
    """Download one exact URL with bounded retries and resumable Range staging."""

    expected = int(asset["expected_size_bytes"])
    part_path.parent.mkdir(parents=True, exist_ok=True)
    retries = 0
    response_headers: dict[str, str] = {}
    http_status = 0
    started = utc_now()
    while True:
        offset = part_path.stat().st_size if part_path.exists() else 0
        if offset > expected:
            raise Tier3ValidationError(f"partial file exceeds expected size: {part_path}")
        request = urllib.request.Request(str(asset["url"]), headers={"User-Agent": "lewontin-vgp-pilot/2.0"})
        if offset:
            request.add_header("Range", f"bytes={offset}-")
        try:
            with urllib.request.urlopen(request, timeout=120) as response:  # noqa: S310 - exact gate-bound HTTPS URL
                http_status = int(getattr(response, "status", response.getcode()))
                response_headers = {key.lower(): value for key, value in response.headers.items()}
                if offset and http_status != 206:
                    # A source that ignored Range is safe only if we restart the
                    # local part.  Bytes from earlier attempts remain accounted.
                    offset = 0
                    part_path.unlink(missing_ok=True)
                mode = "ab" if offset else "wb"
                with part_path.open(mode) as handle:
                    while True:
                        chunk = response.read(1024 * 1024)
                        if not chunk:
                            break
                        new_total = int(transfer_state["transferred_bytes"]) + len(chunk)
                        if new_total > int(transfer_state["transfer_cap_bytes"]):
                            raise Tier3ValidationError("cumulative actual transfer bytes exceed strict compressed-input cap")
                        if handle.tell() + len(chunk) > expected:
                            raise Tier3ValidationError(f"response exceeds approved expected size: {asset['url']}")
                        handle.write(chunk)
                        transfer_state["transferred_bytes"] = new_total
                    handle.flush()
                    os.fsync(handle.fileno())
            break
        except (OSError, urllib.error.URLError) as error:
            if retries >= 3:
                raise Tier3ValidationError(f"retrieval failed after {retries + 1} attempts: {asset['url']}: {error}") from error
            retries += 1
            time.sleep(min(2**retries, 8))
    return {
        "started_at_utc": started,
        "completed_at_utc": utc_now(),
        "http_status": http_status,
        "response_headers": response_headers,
        "retries": retries,
    }


def _hashes(path: Path) -> tuple[str, str, int]:
    sha = hashlib.sha256()
    md5 = hashlib.md5()  # noqa: S324 - official provider checksum, never used as the object identity
    size = 0
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            sha.update(chunk)
            md5.update(chunk)
            size += len(chunk)
    return sha.hexdigest(), md5.hexdigest(), size


def _quarantine(part_path: Path, quarantine_root: Path, reason: str) -> Path:
    quarantine_root.mkdir(parents=True, exist_ok=True)
    target = quarantine_root / part_path.name
    suffix = 1
    while target.exists():
        target = quarantine_root / f"{part_path.name}.{suffix}"
        suffix += 1
    if part_path.exists():
        os.replace(part_path, target)
        target.chmod(stat.S_IRUSR | stat.S_IRGRP)
    reason_path = target.with_name(target.name + ".reason.json")
    _write_json(reason_path, {"quarantined_at_utc": utc_now(), "reason": reason, "payload": str(target)})
    return target


def promote_verified_part(
    part_path: Path,
    object_root: Path,
    *,
    expected_size: int,
    expected_source_checksum: Mapping[str, str] | None = None,
    before_reverify: PromoteHook | None = None,
) -> dict[str, Any]:
    """Verify a complete part twice, then atomically promote it read-only."""

    first_sha, first_md5, first_size = _hashes(part_path)
    if first_size != expected_size:
        raise Tier3ValidationError(f"staged size mismatch: expected {expected_size}, observed {first_size}")
    source_verified: bool | None = None
    if expected_source_checksum:
        algorithm = expected_source_checksum.get("algorithm")
        value = str(expected_source_checksum.get("value", "")).lower()
        if algorithm != "md5" or not HEX32.fullmatch(value):
            raise Tier3ValidationError("unsupported or invalid official source checksum")
        source_verified = first_md5 == value
        if not source_verified:
            raise Tier3ValidationError(f"official MD5 mismatch: expected {value}, observed {first_md5}")
    if before_reverify:
        before_reverify(part_path, first_sha)
    second_sha, second_md5, second_size = _hashes(part_path)
    if (second_sha, second_md5, second_size) != (first_sha, first_md5, first_size):
        raise Tier3ValidationError("staged content or local digest changed before promotion")
    object_root.mkdir(parents=True, exist_ok=True)
    destination = object_root / first_sha[:2] / first_sha
    destination.parent.mkdir(parents=True, exist_ok=True)
    if part_path.stat().st_dev != destination.parent.stat().st_dev:
        raise Tier3ValidationError("staging and immutable promotion target are not on one filesystem")
    if destination.exists():
        if sha256_file(destination) != first_sha or destination.stat().st_size != first_size:
            raise Tier3ValidationError("existing immutable object does not match its content address")
        part_path.unlink()
    else:
        os.replace(part_path, destination)
    destination.chmod(stat.S_IRUSR | stat.S_IRGRP)
    return {
        "local_sha256": first_sha,
        "local_md5": first_md5,
        "bytes": first_size,
        "source_checksum_verified": source_verified,
        "object_path": str(destination),
        "mode_octal": oct(stat.S_IMODE(destination.stat().st_mode)),
    }


def default_staged_validator(candidate_id: str, staged: Sequence[Mapping[str, Any]]) -> Sequence[str]:
    """Reassert the exact-reference role set before any candidate promotion.

    Detailed sequence-region and pair facts are gate-bound to immutable official
    responses.  If a future gate adds staged report assets, they remain part of
    this exact obligation set; this function never infers a replacement pair or
    annotation.
    """

    roles = {str(item["asset"]["role"]) for item in staged}
    if not {"h1_fasta", "native_h1_annotation"}.issubset(roles):
        raise Tier3ValidationError(f"staged H1/native-annotation linkage set incomplete: {candidate_id}")
    accessions = {
        str(item["asset"]["accession_version"])
        for item in staged
        if item["asset"]["role"] in {"h1_fasta", "native_h1_annotation"}
    }
    if len(accessions) != 1:
        raise Tier3ValidationError(f"staged H1/native annotation exact-reference linkage mismatch: {candidate_id}")
    return ["exact_accession_version_revalidated", "h1_native_annotation_reference_linkage_revalidated"]


def _empty_manifest_row() -> dict[str, Any]:
    return {field: "" for field in MANIFEST_FIELDS}


def refusal_rows(
    *,
    run_id: str,
    generated_at_utc: str,
    gate_path: Path,
    gate_payload: Mapping[str, Any] | None,
    error: Exception,
) -> list[dict[str, Any]]:
    code = classify_failure_code(error)
    decision_sha = str(gate_payload.get("decision_sha256", "")) if gate_payload else ""
    rows: list[dict[str, Any]] = []
    summary = _empty_manifest_row()
    summary.update(
        {
            "run_id": run_id,
            "generated_at_utc": generated_at_utc,
            "started_at_utc": generated_at_utc,
            "completed_at_utc": utc_now(),
            "record_type": "run_summary",
            "status": "refused_preflight",
            "asset_role": "authorization_boundary",
            "source_url": str(gate_path),
            "expected_bytes": "0",
            "observed_bytes": "0",
            "cumulative_transferred_bytes": "0",
            "expected_sha256": decision_sha,
            "observed_sha256": decision_sha,
            "retries": "0",
            "validation_outcomes_json": json.dumps(["refused_before_downloader", "zero_biological_bytes"], separators=(",", ":")),
            "failure_code": code,
            "failure_source": str(gate_path),
            "failure_message": str(error),
        }
    )
    rows.append(summary)
    blockers = gate_payload.get("blockers", []) if gate_payload else []
    for blocker in blockers if isinstance(blockers, list) else []:
        row = _empty_manifest_row()
        row.update(
            {
                "run_id": run_id,
                "generated_at_utc": generated_at_utc,
                "started_at_utc": generated_at_utc,
                "completed_at_utc": summary["completed_at_utc"],
                "record_type": "gate_blocker",
                "status": "refused_preflight",
                "asset_role": "authorization_boundary",
                "source_url": blocker.get("source", ""),
                "expected_bytes": "0",
                "observed_bytes": "0",
                "cumulative_transferred_bytes": "0",
                "retries": "0",
                "failure_code": blocker.get("code", "UNKNOWN_GATE_BLOCKER"),
                "failure_source": blocker.get("source", ""),
                "failure_message": blocker.get("message", ""),
            }
        )
        rows.append(row)
    return rows


def _gate_digest(payload: Mapping[str, Any] | None, name: str) -> str:
    boundary = payload.get("authorization_boundary", {}) if payload else {}
    return str(boundary.get(name, "UNAVAILABLE")) if isinstance(boundary, Mapping) else "UNAVAILABLE"


def render_report(
    *,
    run_id: str,
    generated_at_utc: str,
    gate_path: Path,
    gate_payload: Mapping[str, Any] | None,
    manifest_path: Path,
    root_config_path: Path,
    output_manifest_path: Path,
    output_inventory_path: Path,
    refusal_evidence_path: Path,
    status: str,
    transferred_bytes: int,
    promoted_count: int,
    quarantined_count: int,
    error: Exception | None,
) -> str:
    slurm_env = detected_slurm_environment()
    lines = [
        "# Repaired VGP pilot acquisition report",
        "",
        f"- Run ID: `{run_id}`",
        f"- Generated at: `{generated_at_utc}`",
        f"- Gate path: `{gate_path}`",
        f"- Gate decision: `{_safe_gate_status(gate_payload)}`",
        f"- Acquisition status: `{status}`",
        f"- Refused before first biological byte: `{str(status == 'refused_preflight').lower()}`",
        f"- Provider requests attempted: `{0 if status == 'refused_preflight' else 'recorded per asset'}`",
        f"- Biological payload bytes transferred: `{transferred_bytes}`",
        f"- Verified immutable objects promoted: `{promoted_count}`",
        f"- Quarantine objects written: `{quarantined_count}`",
        f"- Slurm environment detected: `{str(bool(slurm_env)).lower()}`",
        f"- Output manifest: `{output_manifest_path}`",
        f"- Immutable-object inventory: `{output_inventory_path}`",
        f"- Refusal evidence: `{refusal_evidence_path if status == 'refused_preflight' else 'not applicable'}`",
        "",
        "## Authorization boundary",
        "",
        f"- Manifest path: `{manifest_path}`",
        f"- Manifest digest: `{_gate_digest(gate_payload, 'manifest_digest')}`",
        f"- Root contract path: `{root_config_path}`",
        f"- Root/storage digest: `{_gate_digest(gate_payload, 'data_root_storage_contract_digest')}`",
        f"- Environment digest: `{_gate_digest(gate_payload, 'environment_digest')}`",
        f"- Cap-vector digest: `{_gate_digest(gate_payload, 'cap_vector_digest')}`",
        f"- Retrieval/checksum digest: `{_gate_digest(gate_payload, 'retrieval_checksum_obligations_digest')}`",
        f"- Pair-evidence digest: `{_gate_digest(gate_payload, 'pair_evidence_digest')}`",
        f"- Measurement-contract digest: `{_gate_digest(gate_payload, 'measurement_contract_digest')}`",
        f"- Authorization-tuple digest: `{_gate_digest(gate_payload, 'authorization_tuple_digest')}`",
        "",
    ]
    if error is not None:
        lines.extend(["## Refusal reason", "", f"- `{classify_failure_code(error)}`: {error}", ""])
    blockers = gate_payload.get("blockers", []) if gate_payload else []
    lines.extend(["## Gate blockers", ""])
    if blockers:
        for blocker in blockers:
            lines.append(f"- `{blocker.get('code')}`: {blocker.get('message')} ({blocker.get('source')})")
    else:
        lines.append("- None recorded, or the gate payload was unavailable.")
    lines.extend(
        [
            "",
            "## Validation notes",
            "",
            "- The gate's nested authorization tuple was verified before the downloader boundary; a literal GO would additionally trigger live recomputation of every bound input and digest.",
            "- Literal GO, exact selected rows/assets, strictest cap values, 25% storage headroom, and explicit scope exclusions are mechanical preconditions.",
            "- Every authorized payload uses resumable `.part` staging, finite-size and cumulative-byte checks, any official checksum, local SHA-256, immediate SHA-256 reverification, and same-filesystem atomic read-only promotion.",
            "- Mismatches are quarantined and excluded from the immutable inventory.",
            "- Exact-reference/native-annotation/Tier3A-pair checks are gate-bound to the pinned GNU Guix environment; on refusal no staged biological report exists to re-run.",
            "- Exact-reference/native-annotation linkage validation under pinned GNU Guix was not re-run because zero assets were authorized or acquired.",
            "- No full catalog or raw population bulk acquisition, Slurm submission, or demographic inference is performed by this entrypoint.",
        ]
    )
    if slurm_env:
        lines.extend(["", "## Detected Slurm environment", "", "```json", json.dumps(slurm_env, indent=2, sort_keys=True), "```"])
    return "\n".join(lines) + "\n"


def _refusal_evidence(
    *, run_id: str, generated_at_utc: str, gate_path: Path, payload: Mapping[str, Any] | None, error: Exception
) -> dict[str, Any]:
    return {
        "schema_version": "1.0",
        "run_id": run_id,
        "generated_at_utc": generated_at_utc,
        "status": "refused_preflight",
        "failure_code": classify_failure_code(error),
        "failure_message": str(error),
        "gate_path": str(gate_path),
        "gate_decision": _safe_gate_status(payload),
        "decision_sha256": str(payload.get("decision_sha256", "")) if payload else "",
        "authorization_tuple_digest": _gate_digest(payload, "authorization_tuple_digest"),
        "provider_requests_attempted": 0,
        "biological_payload_bytes_transferred": 0,
        "partial_files_created": 0,
        "objects_promoted": 0,
        "quarantine_objects_written": 0,
        "slurm_jobs_submitted": 0,
        "demographic_inferences_performed": 0,
        "full_catalog_downloaded": False,
        "raw_population_bulk_downloaded": False,
    }


def _layout(root_config_path: Path) -> tuple[Path, dict[str, Path]]:
    config = json.loads(root_config_path.read_text(encoding="utf-8"))
    root = Path(config["root"]).resolve(strict=True)
    paths = {name: (root / value).resolve(strict=False) for name, value in config["layout"].items()}
    for path in paths.values():
        try:
            path.relative_to(root)
        except ValueError as error:
            raise Tier3ValidationError(f"configured acquisition path escapes authorized root: {path}") from error
    return root, paths


def _validate_live_storage(payload: Mapping[str, Any], root: Path) -> None:
    storage = payload["storage_audit"]
    if storage.get("adequate") is not True:
        raise Tier3ValidationError("root/storage contract is not adequate")
    allocation = storage.get("enforceable_allocation", {})
    if allocation.get("status") in {None, "", "unknown"} or allocation.get("headroom_pass") is not True:
        raise Tier3ValidationError("root/storage quota or allocation headroom is unavailable")
    recorded = storage.get("live_identity", {})
    if recorded.get("catalog_on_same_filesystem") is not True:
        raise Tier3ValidationError("root/storage catalog filesystem identity is not bound")
    observed = root.stat()
    identity_checks = {
        "inode": observed.st_ino,
        "uid": observed.st_uid,
        "gid": observed.st_gid,
        "mode_octal": oct(stat.S_IMODE(observed.st_mode)),
        "resolved_path": str(root.resolve(strict=True)),
    }
    for key, value in identity_checks.items():
        if recorded.get(key) != value:
            raise Tier3ValidationError(f"root/storage live identity changed: {key}")
    filesystem = os.statvfs(root)
    free_bytes = filesystem.f_bavail * filesystem.f_frsize
    free_inodes = filesystem.f_favail
    worst_case = storage.get("worst_case", {})
    if free_bytes < int(worst_case.get("required_bytes_with_headroom", 0)):
        raise Tier3ValidationError("live filesystem byte headroom is below the bound storage contract")
    if free_inodes < int(worst_case.get("required_inodes_with_headroom", 0)):
        raise Tier3ValidationError("live filesystem inode headroom is below the bound storage contract")


def _transfer_byte_cap(payload: Mapping[str, Any]) -> int:
    dimensions = payload["cap_vector"]["dimensions"]
    candidates = [int(float(dimensions["compressed_inputs_gib"]["limit"]) * GIB)]
    for name in ("persistent_input_gb", "moosefs_write_gb"):
        if name in dimensions:
            candidates.append(int(float(dimensions[name]["limit"]) * 1_000_000_000))
    allocation = payload["storage_audit"].get("enforceable_allocation", {})
    remaining_bytes = allocation.get("remaining_bytes")
    if isinstance(remaining_bytes, int) and remaining_bytes > 0:
        candidates.append(remaining_bytes)
    return min(candidates)


def _asset_part_name(asset: Mapping[str, Any]) -> str:
    basename = str(asset["url"]).rsplit("/", 1)[-1]
    safe = re.sub(r"[^A-Za-z0-9._-]", "_", basename)
    return f"{asset['candidate_id']}.{asset['role']}.{safe}.part"


def _acquire_authorized(
    *,
    payload: Mapping[str, Any],
    authorization: Mapping[str, Any],
    assets: Sequence[Mapping[str, Any]],
    root_config_path: Path,
    run_id: str,
    generated_at_utc: str,
    downloader: DownloadFunction,
    before_reverify: PromoteHook | None,
    staged_validator: StagedValidator,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], int, int]:
    root, paths = _layout(root_config_path)
    _validate_live_storage(payload, root)
    part_root = paths["staging_partials"] / run_id
    quarantine_root = paths["quarantine"] / run_id
    object_root = paths["immutable_objects"]
    transfer_cap = _transfer_byte_cap(payload)
    expected_total = sum(int(asset["expected_size_bytes"]) for asset in assets)
    if expected_total > transfer_cap:
        raise Tier3ValidationError(
            f"approved expected bytes {expected_total} exceed strictest transfer/storage cap {transfer_cap}"
        )
    transfer_state: MutableMapping[str, Any] = {"transferred_bytes": 0, "transfer_cap_bytes": transfer_cap}
    staged_by_candidate: dict[str, list[dict[str, Any]]] = {}
    rows: list[dict[str, Any]] = []
    inventory: list[dict[str, Any]] = []
    quarantined = 0

    for asset in assets:
        part_path = part_root / _asset_part_name(asset)
        metadata = downloader(asset, part_path, transfer_state)
        item = {"asset": asset, "part_path": part_path, "response": dict(metadata)}
        staged_by_candidate.setdefault(str(asset["candidate_id"]), []).append(item)

    validation_by_candidate: dict[str, Sequence[str]] = {}
    for candidate_id, staged in staged_by_candidate.items():
        try:
            validation_by_candidate[candidate_id] = staged_validator(candidate_id, staged)
        except (OSError, Tier3ValidationError, ValueError) as error:
            validation_by_candidate[candidate_id] = [f"candidate_validation_failed:{error}"]
            for item in staged:
                part_path = item["part_path"]
                quarantine_path = _quarantine(part_path, quarantine_root, str(error))
                quarantined += int(quarantine_path.exists())
                asset = item["asset"]
                row = _empty_manifest_row()
                row.update(
                    {
                        "run_id": run_id,
                        "generated_at_utc": generated_at_utc,
                        "started_at_utc": item["response"].get("started_at_utc", ""),
                        "completed_at_utc": utc_now(),
                        "record_type": "asset",
                        "status": "quarantined",
                        "candidate_id": candidate_id,
                        "asset_role": asset["role"],
                        "accession_version": asset["accession_version"],
                        "source_url": asset["url"],
                        "expected_bytes": asset["expected_size_bytes"],
                        "observed_bytes": quarantine_path.stat().st_size if quarantine_path.exists() else 0,
                        "cumulative_transferred_bytes": transfer_state["transferred_bytes"],
                        "retries": item["response"].get("retries", 0),
                        "validation_outcomes_json": json.dumps(validation_by_candidate[candidate_id], separators=(",", ":")),
                        "failure_code": "REFERENCE_LINKAGE_VALIDATION_FAILED",
                        "failure_source": str(part_path),
                        "failure_message": str(error),
                        "staging_path": str(part_path),
                        "quarantine_path": str(quarantine_path),
                    }
                )
                rows.append(row)

    for candidate_id, staged in staged_by_candidate.items():
        if validation_by_candidate[candidate_id] and str(validation_by_candidate[candidate_id][0]).startswith(
            "candidate_validation_failed:"
        ):
            continue
        for item in staged:
            asset = item["asset"]
            part_path = item["part_path"]
            checksum = asset.get("source_checksum")
            try:
                promoted = promote_verified_part(
                    part_path,
                    object_root,
                    expected_size=int(asset["expected_size_bytes"]),
                    expected_source_checksum=checksum,
                    before_reverify=before_reverify,
                )
            except (OSError, Tier3ValidationError) as error:
                quarantine_path = _quarantine(part_path, quarantine_root, str(error))
                quarantined += int(quarantine_path.exists())
                row = _empty_manifest_row()
                row.update(
                    {
                        "run_id": run_id,
                        "generated_at_utc": generated_at_utc,
                        "started_at_utc": item["response"].get("started_at_utc", ""),
                        "completed_at_utc": utc_now(),
                        "record_type": "asset",
                        "status": "quarantined",
                        "candidate_id": candidate_id,
                        "asset_role": asset["role"],
                        "accession_version": asset["accession_version"],
                        "source_url": asset["url"],
                        "expected_bytes": asset["expected_size_bytes"],
                        "observed_bytes": quarantine_path.stat().st_size if quarantine_path.exists() else "0",
                        "cumulative_transferred_bytes": transfer_state["transferred_bytes"],
                        "source_checksum_algorithm": checksum.get("algorithm", "") if checksum else "",
                        "source_checksum_value": checksum.get("value", "") if checksum else "",
                        "provider_md5": checksum.get("value", "") if checksum and checksum.get("algorithm") == "md5" else "",
                        "http_status": item["response"].get("http_status", ""),
                        "response_headers_json": json.dumps(item["response"].get("response_headers", {}), sort_keys=True, separators=(",", ":")),
                        "retries": item["response"].get("retries", 0),
                        "validation_outcomes_json": json.dumps(validation_by_candidate[candidate_id], separators=(",", ":")),
                        "failure_code": "STAGED_VALIDATION_FAILED",
                        "failure_source": str(part_path),
                        "failure_message": str(error),
                        "staging_path": str(part_path),
                        "quarantine_path": str(quarantine_path),
                    }
                )
                rows.append(row)
                continue
            checksum_verified = promoted["source_checksum_verified"]
            row = _empty_manifest_row()
            row.update(
                {
                    "run_id": run_id,
                    "generated_at_utc": generated_at_utc,
                    "started_at_utc": item["response"].get("started_at_utc", ""),
                    "completed_at_utc": utc_now(),
                    "record_type": "asset",
                    "status": "promoted",
                    "candidate_id": candidate_id,
                    "asset_role": asset["role"],
                    "accession_version": asset["accession_version"],
                    "source_url": asset["url"],
                    "expected_bytes": asset["expected_size_bytes"],
                    "observed_bytes": promoted["bytes"],
                    "cumulative_transferred_bytes": transfer_state["transferred_bytes"],
                    "expected_sha256": "local_after_complete_staging",
                    "observed_sha256": promoted["local_sha256"],
                    "source_checksum_algorithm": checksum.get("algorithm", "") if checksum else "",
                    "source_checksum_value": checksum.get("value", "") if checksum else "",
                    "source_checksum_verified": "true" if checksum_verified else ("false" if checksum_verified is False else "not_available"),
                    "provider_md5": checksum.get("value", "") if checksum and checksum.get("algorithm") == "md5" else "",
                    "http_status": item["response"].get("http_status", ""),
                    "response_headers_json": json.dumps(item["response"].get("response_headers", {}), sort_keys=True, separators=(",", ":")),
                    "retries": item["response"].get("retries", 0),
                    "validation_outcomes_json": json.dumps(
                        [*validation_by_candidate[candidate_id], "size_verified", "local_sha256_computed", "local_sha256_reverified", "atomic_read_only_promotion"],
                        separators=(",", ":"),
                    ),
                    "staging_path": str(part_path),
                    "promoted_path": promoted["object_path"],
                }
            )
            rows.append(row)
            inventory.append(
                {
                    "run_id": run_id,
                    "promoted_at_utc": row["completed_at_utc"],
                    "candidate_id": candidate_id,
                    "asset_role": asset["role"],
                    "accession_version": asset["accession_version"],
                    "source_url": asset["url"],
                    "bytes": promoted["bytes"],
                    "local_sha256": promoted["local_sha256"],
                    "source_checksum_algorithm": row["source_checksum_algorithm"],
                    "source_checksum_value": row["source_checksum_value"],
                    "source_checksum_verified": row["source_checksum_verified"],
                    "object_path": promoted["object_path"],
                    "mode_octal": promoted["mode_octal"],
                    "authorization_tuple_digest": authorization["authorization_tuple_digest"],
                    "decision_sha256": authorization["decision_sha256"],
                }
            )
    return rows, inventory, int(transfer_state["transferred_bytes"]), quarantined


def run(
    *,
    gate_path: Path = DEFAULT_GATE,
    manifest_path: Path = DEFAULT_MANIFEST,
    root_config_path: Path = DEFAULT_ROOT_CONFIG,
    output_manifest_path: Path = DEFAULT_OUTPUT_MANIFEST,
    output_report_path: Path = DEFAULT_OUTPUT_REPORT,
    output_inventory_path: Path = DEFAULT_OUTPUT_INVENTORY,
    refusal_evidence_path: Path = DEFAULT_REFUSAL_EVIDENCE,
    downloader: DownloadFunction | None = None,
    before_reverify: PromoteHook | None = None,
    staged_validator: StagedValidator = default_staged_validator,
) -> dict[str, Any]:
    run_id = make_run_id()
    generated_at_utc = utc_now()
    untrusted_payload = _load_untrusted_gate(gate_path)
    try:
        payload, authorization, assets = preflight_authorize(gate_path, manifest_path, root_config_path)
    except (Tier3ValidationError, FileNotFoundError, KeyError, TypeError, ValueError, OSError) as error:
        rows = refusal_rows(
            run_id=run_id,
            generated_at_utc=generated_at_utc,
            gate_path=gate_path,
            gate_payload=untrusted_payload,
            error=error,
        )
        write_manifest(output_manifest_path, rows)
        write_inventory(output_inventory_path, [])
        evidence = _refusal_evidence(
            run_id=run_id, generated_at_utc=generated_at_utc, gate_path=gate_path, payload=untrusted_payload, error=error
        )
        _write_json(refusal_evidence_path, evidence)
        _write_text(
            output_report_path,
            render_report(
                run_id=run_id,
                generated_at_utc=generated_at_utc,
                gate_path=gate_path,
                gate_payload=untrusted_payload,
                manifest_path=manifest_path,
                root_config_path=root_config_path,
                output_manifest_path=output_manifest_path,
                output_inventory_path=output_inventory_path,
                refusal_evidence_path=refusal_evidence_path,
                status="refused_preflight",
                transferred_bytes=0,
                promoted_count=0,
                quarantined_count=0,
                error=error,
            ),
        )
        return {
            "run_id": run_id,
            "generated_at_utc": generated_at_utc,
            "status": "refused_preflight",
            "failure_code": classify_failure_code(error),
            "failure_message": str(error),
            "transferred_bytes": 0,
            "provider_requests_attempted": 0,
            "promoted_count": 0,
            "quarantined_count": 0,
            "output_manifest": str(output_manifest_path),
            "output_report": str(output_report_path),
            "output_inventory": str(output_inventory_path),
            "refusal_evidence": str(refusal_evidence_path),
            "row_count": len(rows),
        }

    refusal_evidence_path.unlink(missing_ok=True)
    rows, inventory, transferred, quarantined = _acquire_authorized(
        payload=payload,
        authorization=authorization,
        assets=assets,
        root_config_path=root_config_path,
        run_id=run_id,
        generated_at_utc=generated_at_utc,
        downloader=downloader or _default_downloader,
        before_reverify=before_reverify,
        staged_validator=staged_validator,
    )
    write_manifest(output_manifest_path, rows)
    write_inventory(output_inventory_path, inventory)
    status = "complete" if len(inventory) == len(assets) else "complete_with_quarantine"
    _write_text(
        output_report_path,
        render_report(
            run_id=run_id,
            generated_at_utc=generated_at_utc,
            gate_path=gate_path,
            gate_payload=payload,
            manifest_path=manifest_path,
            root_config_path=root_config_path,
            output_manifest_path=output_manifest_path,
            output_inventory_path=output_inventory_path,
            refusal_evidence_path=refusal_evidence_path,
            status=status,
            transferred_bytes=transferred,
            promoted_count=len(inventory),
            quarantined_count=quarantined,
            error=None,
        ),
    )
    return {
        "run_id": run_id,
        "generated_at_utc": generated_at_utc,
        "status": status,
        "transferred_bytes": transferred,
        "provider_requests_attempted": len(assets),
        "promoted_count": len(inventory),
        "quarantined_count": quarantined,
        "output_manifest": str(output_manifest_path),
        "output_report": str(output_report_path),
        "output_inventory": str(output_inventory_path),
        "row_count": len(rows),
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gate", type=Path, default=DEFAULT_GATE)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--root-config", type=Path, default=DEFAULT_ROOT_CONFIG)
    parser.add_argument("--output-manifest", type=Path, default=DEFAULT_OUTPUT_MANIFEST)
    parser.add_argument("--output-report", type=Path, default=DEFAULT_OUTPUT_REPORT)
    parser.add_argument("--output-inventory", type=Path, default=DEFAULT_OUTPUT_INVENTORY)
    parser.add_argument("--refusal-evidence", type=Path, default=DEFAULT_REFUSAL_EVIDENCE)
    args = parser.parse_args(argv)
    result = run(
        gate_path=args.gate,
        manifest_path=args.manifest,
        root_config_path=args.root_config,
        output_manifest_path=args.output_manifest,
        output_report_path=args.output_report,
        output_inventory_path=args.output_inventory,
        refusal_evidence_path=args.refusal_evidence,
    )
    print(json.dumps(result, sort_keys=True))
    return 0 if result["status"] in {"complete", "complete_with_quarantine", "refused_preflight"} else 1


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
