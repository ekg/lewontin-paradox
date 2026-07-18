#!/usr/bin/env python3
"""Run the repaired, gate-bound VGP pilot or emit zero-use refusal evidence.

The literal repaired GO decision is the first authorization boundary.  A GO
still does not authorize compute until all nested gate bindings, the live gate
reproduction, the acquisition ledger, the immutable-object inventory, local
payload hashes, exact reference/annotation roles, pair evidence, pinned Guix
environment, and strict resource caps agree.  The committed decision is
NO_GO, so the normal repository invocation deliberately never calls Slurm.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import stat
import subprocess
import sys
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Mapping, Sequence

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from analysis import acquire_vgp_pilot as acquisition
from analysis import gate_vgp_pilot as gate
from analysis.tier3_common import Tier3ValidationError, sha256_file


PROJECT_ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = PROJECT_ROOT / "analysis"
DEFAULT_GATE = ANALYSIS / "vgp_pilot_gate.json"
DEFAULT_MANIFEST = ANALYSIS / "vgp_pilot_manifest.tsv"
DEFAULT_ROOT_CONFIG = ANALYSIS / "vgp_data_root_config.json"
DEFAULT_ACQUISITION = ANALYSIS / "vgp_pilot_acquisition_manifest.tsv"
DEFAULT_INVENTORY = ANALYSIS / "vgp_pilot_immutable_object_inventory.tsv"
DEFAULT_SWEEPGA_BUILD = ANALYSIS / "sweepga_origin_main_build.json"
DEFAULT_IMPG_HANDOFF = ANALYSIS / "sweepga_impg_observed.json"
DEFAULT_WORKER = ANALYSIS / "slurm" / "run_repaired_vgp.sh"
DEFAULT_OUTPUT_RUN_MANIFEST = ANALYSIS / "vgp_pilot_run_manifest.tsv"
DEFAULT_OUTPUT_SLURM_TELEMETRY = ANALYSIS / "vgp_pilot_slurm_telemetry.tsv"
DEFAULT_OUTPUT_RESULTS = ANALYSIS / "vgp_pilot_results.tsv"
DEFAULT_OUTPUT_EXCLUSIONS = ANALYSIS / "vgp_pilot_exclusions.tsv"
DEFAULT_OUTPUT_REFUSALS = ANALYSIS / "vgp_pilot_refusals.tsv"
DEFAULT_OUTPUT_REPORT = ANALYSIS / "vgp_pilot_run_report.md"

HEX64 = re.compile(r"^[0-9a-f]{64}$")
GIB = 1024**3
HARD_CAPS = acquisition.HARD_CAP_CEILINGS
FORBIDDEN_SCOPE_KEYS = (
    "full_catalog_download_authorized",
    "raw_population_bulk_download_authorized",
    "demographic_inference_authorized",
)

RUN_MANIFEST_FIELDS = [
    "run_id", "generated_at_utc", "record_type", "status", "candidate_id",
    "gate_decision", "decision_sha256", "authorization_tuple_digest",
    "gate_file_sha256", "manifest_sha256", "acquisition_manifest_sha256",
    "immutable_inventory_sha256", "root_config_sha256", "input_bundle_digest",
    "root_contract_digest", "cap_vector_digest", "retrieval_digest",
    "pair_evidence_digest", "measurement_contract_digest", "environment_digest",
    "sweepga_build_sha256", "impg_handoff_sha256", "worker_sha256",
    "selected_species", "promoted_objects", "compressed_input_bytes",
    "slurm_jobs_submitted", "final_state", "failure_code", "failure_source",
    "failure_message", "notes",
]

SLURM_TELEMETRY_FIELDS = [
    "run_id", "generated_at_utc", "record_type", "status", "candidate_id",
    "sbatch_command", "slurm_job_id", "slurm_array_job_id", "dependency", "requested_cpus",
    "requested_memory_gib", "requested_wall_hours", "requested_scratch_gb",
    "requested_read_gb", "requested_write_gb",
    "retry_index", "final_state", "max_rss_gib", "elapsed_seconds",
    "cpu_time_seconds", "scratch_peak_gb", "io_read_gb", "io_write_gb",
    "metadata_operations", "network_bytes", "failure_code", "notes",
]

RESULT_FIELDS = [
    "run_id", "generated_at_utc", "record_type", "status", "candidate_id",
    "scientific_name", "result_scope", "metric", "numerator", "denominator",
    "target_total", "value", "measurement_method", "artifact_sha256",
    "exclusion_reason", "failure_code", "failure_source", "notes",
]

EXCLUSION_FIELDS = [
    "run_id", "generated_at_utc", "candidate_id", "result_scope", "metric",
    "status", "threshold", "observed", "failure_code", "failure_source",
    "reason", "imputed", "demographic_input_used",
]

REFUSAL_FIELDS = [
    "run_id", "generated_at_utc", "gate_decision", "status", "failure_code",
    "failure_source", "failure_message", "sbatch_commands_issued",
    "slurm_jobs_submitted", "compute_jobs_started", "core_seconds",
    "scratch_bytes", "io_read_bytes", "io_write_bytes", "network_bytes",
    "provider_requests", "full_catalog_downloaded", "population_bulk_downloaded",
    "demographic_inferences", "evidence_sha256",
]

Submitter = Callable[[Sequence[str]], str]


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def make_run_id() -> str:
    return f"vgp-pilot-run-{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}-{os.getpid()}"


def _atomic_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=path.parent, prefix=f".{path.name}.", delete=False) as handle:
        handle.write(text)
        handle.flush()
        os.fsync(handle.fileno())
        temporary = Path(handle.name)
    os.replace(temporary, path)


def atomic_write_tsv(path: Path, fields: Sequence[str], rows: Sequence[Mapping[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", newline="", dir=path.parent, prefix=f".{path.name}.", delete=False) as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
        handle.flush()
        os.fsync(handle.fileno())
        temporary = Path(handle.name)
    os.replace(temporary, path)


def read_json_object(path: Path, label: str) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as error:
        raise Tier3ValidationError(f"invalid or unavailable {label}: {error}") from error
    if not isinstance(value, dict):
        raise Tier3ValidationError(f"{label} must be a JSON object")
    return value


def read_tsv(path: Path, label: str) -> list[dict[str, str]]:
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))
    except OSError as error:
        raise Tier3ValidationError(f"invalid or unavailable {label}: {error}") from error


def safe_sha256(path: Path) -> str:
    try:
        return sha256_file(path)
    except OSError:
        return ""


def classify_failure_code(error: Exception) -> str:
    text = str(error).lower()
    checks = (
        (("gate decision is no_go",), "GATE_NO_GO"),
        (("unknown gate decision", "not the literal go token"), "GATE_DECISION_UNKNOWN"),
        (("local payload sha-256", "content address", "payload digest"), "LOCAL_PAYLOAD_DIGEST_MISMATCH"),
        (("acquisition is incomplete", "missing promoted", "duplicate promoted", "not promoted"), "ACQUISITION_INCOMPLETE"),
        (("accession/version substitution", "approved asset identity", "url substitution"), "ACCESSION_VERSION_SUBSTITUTION"),
        (("native annotation", "reference linkage"), "NATIVE_ANNOTATION_MISMATCH"),
        (("tier3a pair", "same-individual", "correctly phased"), "TIER3A_PAIR_MISMATCH"),
        (("cap vector", "strictest integrated", "storage contract", "storage headroom"), "CAP_STORAGE_CONTRACT_WEAKENED"),
        (("manifest digest",), "MANIFEST_DIGEST_MISMATCH"),
        (("digest mismatch", "authorization tuple", "bound input bundle"), "BOUND_DIGEST_MISMATCH"),
        (("sweepga", "impg", "guix"), "ENVIRONMENT_CONTRACT_MISMATCH"),
        (("unavailable", "does not exist", "no such file"), "BOUND_INPUT_MISSING"),
    )
    for needles, code in checks:
        if any(needle in text for needle in needles):
            return code
    return "PRECHECK_FAILED"


def audit_sweepga_origin_build(path: Path) -> dict[str, Any]:
    payload = read_json_object(path, "SweepGA approved-source build record")
    resolution = payload.get("resolution", {})
    binary = payload.get("binary", {})
    if payload.get("completion_gate_passed") is not True:
        raise Tier3ValidationError("SweepGA completion gate did not pass")
    if resolution.get("accepted_syntax") != "--num-mappings 1:1":
        raise Tier3ValidationError("SweepGA mapping contract is not --num-mappings 1:1")
    if not all(resolution.get(key) is True for key in (
        "all_three_whole_haplotype_mappings_completed",
        "all_three_native_multiplicity_rechecks_passed",
        "all_three_impg_handoffs_completed",
    )):
        raise Tier3ValidationError("SweepGA approved-source whole-haplotype/1:1/IMPG audit is incomplete")
    expected = "fa7f0edb9b7e275c288db254046020e136d4267dd5ee043379227ef80da0573b"
    if binary.get("sha256_build_1") != expected or binary.get("sha256_build_2") != expected:
        raise Tier3ValidationError("SweepGA reproducible binary digest differs from the approved source build")
    return payload


def audit_impg_handoff(path: Path) -> dict[str, Any]:
    payload = read_json_object(path, "IMPG handoff record")
    environment = payload.get("environment", {})
    biological = payload.get("biological", {})
    annotation = biological.get("annotation", {})
    mapping = biological.get("sweepga_mapping", {})
    query = biological.get("annotation_query_qc", {})
    if not str(environment.get("guix_profile", "")).startswith("/gnu/store/"):
        raise Tier3ValidationError("IMPG handoff lacks a pinned GNU Guix profile")
    if environment.get("sweepga_commit") != "018e4ce49d2c125820e0ac50dc5feaa02d423683":
        raise Tier3ValidationError("SweepGA approved source commit changed")
    if environment.get("impg_commit") != "101df81eb28a809c8fac97d297acd9fcfbbfa048":
        raise Tier3ValidationError("IMPG approved source commit changed")
    if annotation.get("native_vs_projected") != "native_exact_assembly_submitted_annotation":
        raise Tier3ValidationError("native annotation mismatch in IMPG contract")
    if mapping.get("max_query_overlap_depth") != 1 or mapping.get("max_target_overlap_depth") != 1:
        raise Tier3ValidationError("SweepGA output violates 1:1 mapping multiplicity")
    if biological.get("normalized_vcf_indexed") is not True or biological.get("normalized_bcf_indexed") is not True:
        raise Tier3ValidationError("IMPG VCF/BCF/index contract is incomplete")
    if biological.get("callable_denominator_bp", 0) <= 0 or query.get("queryable_gene_count", 0) <= 0:
        raise Tier3ValidationError("IMPG reference/region contract has an impossible denominator")
    return payload


def audit_pair_contract(payload: Mapping[str, Any], selected: Sequence[str]) -> None:
    row_map = {row["candidate_id"]: row for row in payload["row_audit"]["rows"]}
    pair_map = {row["candidate_id"]: row for row in payload["pair_evidence"]["rows"]}
    for candidate_id in selected:
        row = row_map[candidate_id]
        pair = pair_map.get(candidate_id, {})
        if row.get("tier3a_required"):
            h1, h2 = pair.get("h1_accession_version"), pair.get("h2_accession_version")
            affirmative = pair.get("same_individual_status") in {"affirmed", "verified", "proven"}
            phased = pair.get("phase_evidence_status") in {"affirmed", "verified", "correctly_phased"}
            if not (h1 and h2 and h1 != h2 and affirmative and phased):
                raise Tier3ValidationError(f"Tier3A pair mismatch; exact same-individual correctly phased H1/H2 evidence is absent: {candidate_id}")


def audit_acquisition(
    payload: Mapping[str, Any], expected_assets: Sequence[Mapping[str, Any]],
    acquisition_path: Path, inventory_path: Path,
) -> dict[str, Any]:
    ledger = read_tsv(acquisition_path, "acquisition manifest")
    inventory = read_tsv(inventory_path, "immutable object inventory")
    promoted = [row for row in ledger if row.get("record_type") == "asset" and row.get("status") == "promoted"]
    expected_keys = {
        (str(a["candidate_id"]), str(a["role"]), str(a["accession_version"]), str(a["url"])): a
        for a in expected_assets
    }
    promoted_keys = [(r["candidate_id"], r["asset_role"], r["accession_version"], r["source_url"]) for r in promoted]
    inventory_keys = [(r["candidate_id"], r["asset_role"], r["accession_version"], r["source_url"]) for r in inventory]
    if len(promoted_keys) != len(set(promoted_keys)) or len(inventory_keys) != len(set(inventory_keys)):
        raise Tier3ValidationError("acquisition is incomplete: duplicate promoted asset identity")
    if set(promoted_keys) != set(expected_keys) or set(inventory_keys) != set(expected_keys):
        extra = (set(promoted_keys) | set(inventory_keys)) - set(expected_keys)
        if extra:
            raise Tier3ValidationError("accession/version substitution or URL substitution in acquisition")
        raise Tier3ValidationError("acquisition is incomplete: missing promoted exact assets")
    if any(r.get("status") != "promoted" for r in ledger if r.get("record_type") == "asset"):
        raise Tier3ValidationError("acquisition is incomplete: non-promoted asset rows exist")
    inv_by_key = {key: row for key, row in zip(inventory_keys, inventory)}
    decision_sha = str(payload["decision_sha256"])
    tuple_digest = str(payload["authorization_boundary"]["authorization_tuple_digest"])
    total = 0
    candidates: set[str] = set()
    for row, key in zip(promoted, promoted_keys):
        expected = expected_keys[key]
        inv = inv_by_key[key]
        if int(row["observed_bytes"]) != int(expected["expected_size_bytes"]) or int(inv["bytes"]) != int(expected["expected_size_bytes"]):
            raise Tier3ValidationError(f"acquisition is incomplete: finite byte size mismatch for {key[0]} {key[1]}")
        digest = row.get("observed_sha256", "")
        if not HEX64.fullmatch(digest) or inv.get("local_sha256") != digest:
            raise Tier3ValidationError(f"local payload SHA-256 ledger mismatch for {key[0]} {key[1]}")
        object_path = Path(inv["object_path"])
        if not object_path.is_file() or object_path.stat().st_size != int(inv["bytes"]):
            raise Tier3ValidationError(f"acquisition is incomplete: promoted object unavailable for {key[0]} {key[1]}")
        if sha256_file(object_path) != digest:
            raise Tier3ValidationError(f"local payload SHA-256 changed after promotion for {key[0]} {key[1]}")
        if stat.S_IMODE(object_path.stat().st_mode) & 0o222:
            raise Tier3ValidationError(f"storage contract weakened: promoted object is writable for {key[0]} {key[1]}")
        if inv.get("decision_sha256") != decision_sha or inv.get("authorization_tuple_digest") != tuple_digest:
            raise Tier3ValidationError("bound digest mismatch between acquisition inventory and repaired GO")
        outcomes = row.get("validation_outcomes_json", "")
        if "local_sha256_reverified" not in outcomes or "atomic_read_only_promotion" not in outcomes:
            raise Tier3ValidationError(f"acquisition is incomplete: promotion validations absent for {key[0]} {key[1]}")
        if key[1] == "native_h1_annotation" and "h1_native_annotation_reference_linkage_revalidated" not in outcomes:
            raise Tier3ValidationError(f"native annotation reference linkage mismatch for {key[0]}")
        total += int(inv["bytes"])
        candidates.add(key[0])
    return {"promoted_objects": len(inventory), "compressed_input_bytes": total, "candidate_ids": sorted(candidates)}


def validate_denominator_packet(packet: Mapping[str, Any], contract: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Validate post-alignment denominators; never infer or impute them."""
    thresholds = contract["minimum_thresholds"]
    required = ("callable_bases", "callable_fraction", "queryable_gene_bases", "queryable_gene_count")
    exclusions: list[dict[str, Any]] = []
    for metric in required:
        value = packet.get(metric)
        minimum = thresholds[metric]
        if not isinstance(value, (int, float)) or isinstance(value, bool) or value <= 0 or value < minimum:
            exclusions.append({"metric": metric, "observed": value, "threshold": minimum, "reason": "missing, biologically impossible, or below repaired gate threshold"})
    for total in ("target_gene_total", "target_base_total"):
        value = packet.get(total)
        if not isinstance(value, int) or isinstance(value, bool) or value <= 0:
            exclusions.append({"metric": total, "observed": value, "threshold": ">0", "reason": "target total missing or biologically impossible"})
    method = packet.get("measurement_method")
    digest = packet.get("artifact_sha256")
    if not isinstance(method, str) or not method or not isinstance(digest, str) or not HEX64.fullmatch(digest):
        exclusions.append({"metric": "denominator_provenance", "observed": "missing", "threshold": "method+SHA-256", "reason": "post-alignment measurement method or artifact digest missing"})
    return exclusions


def derive_resource_plan(payload: Mapping[str, Any], selected: Sequence[str]) -> dict[str, Any]:
    dimensions = payload["cap_vector"]["dimensions"]
    acquisition.verify_cap_contract(payload)
    # Enforce every integrated dimension, not only the six task-level hard
    # ceilings.  A future GO may add stricter I/O, inode, bandwidth, wall, or
    # persistent-storage sources; the minimum finite source always wins.
    for name, dimension in dimensions.items():
        if not isinstance(dimension, Mapping):
            raise Tier3ValidationError(f"cap vector dimension is malformed: {name}")
        try:
            limit = float(dimension["limit"])
            observed_value = float(dimension["observed"])
            source_limits = [float(source["limit"]) for source in dimension.get("sources", [])]
        except (KeyError, TypeError, ValueError) as error:
            raise Tier3ValidationError(f"cap vector dimension is not finite: {name}") from error
        if source_limits and limit != min(source_limits):
            raise Tier3ValidationError(f"cap vector relaxed {name}; strictest integrated cap does not win")
        if observed_value > limit or dimension.get("passes") is not True:
            raise Tier3ValidationError(f"cap vector {name} exceeds its strictest integrated limit")
    storage = payload.get("storage_audit", {})
    allocation = storage.get("enforceable_allocation", {}) if isinstance(storage, Mapping) else {}
    if storage.get("adequate") is not True or allocation.get("status") != "known" or allocation.get("headroom_pass") is not True:
        raise Tier3ValidationError("storage contract has no adequate enforceable allocation and required headroom")
    observed = payload["cap_vector"]["proposed_metadata_ready"]
    count = len(selected)
    if count < 1 or count > int(dimensions["species"]["limit"]):
        raise Tier3ValidationError("cap vector species limit violated")
    total_core_hours = float(observed["core_hours"])
    total_scratch = float(observed["scratch_gib"])
    memory = float(observed["memory_per_job_gib"])
    cpus = int(float(observed.get("cpus_per_element", 1)))
    wall = float(observed["aggregate_wall_hours"])
    if total_core_hours > float(dimensions["core_hours"]["limit"]):
        raise Tier3ValidationError("cap vector core-hour limit violated")
    if total_scratch > float(dimensions["scratch_gib"]["limit"]):
        raise Tier3ValidationError("cap vector scratch limit violated")
    if memory > float(dimensions["memory_per_job_gib"]["limit"]):
        raise Tier3ValidationError("cap vector memory-per-job limit violated")
    if int(float(dimensions["concurrent_species"]["limit"])) > 2:
        raise Tier3ValidationError("cap vector concurrent-species contract was weakened")
    return {
        "cpus": cpus,
        "memory_gib": memory,
        "wall_hours": max(1.0, wall / count),
        "scratch_gib": max(1.0, total_scratch / count),
        "total_core_hours": total_core_hours,
        "concurrency": min(2, int(float(dimensions["concurrent_species"]["limit"]))),
        "max_retries": 2,
    }


def _default_submitter(argv: Sequence[str]) -> str:
    result = subprocess.run(list(argv), check=True, text=True, capture_output=True)
    job_id = result.stdout.strip().split(";", 1)[0]
    if not job_id.isdigit():
        raise Tier3ValidationError(f"sbatch did not return a numeric job id: {result.stdout!r}")
    return job_id


def _wait_terminal(job_ids: Sequence[str], timeout_seconds: int = 7 * 24 * 3600) -> list[dict[str, str]]:
    deadline = time.monotonic() + timeout_seconds
    terminal = {"COMPLETED", "FAILED", "CANCELLED", "TIMEOUT", "OUT_OF_MEMORY", "PREEMPTED", "NODE_FAIL", "BOOT_FAIL", "DEADLINE"}
    while time.monotonic() < deadline:
        command = ["sacct", "-n", "-P", "-j", ",".join(job_ids), "--format=JobIDRaw,State,ElapsedRaw,TotalCPU,MaxRSS,MaxDiskRead,MaxDiskWrite"]
        result = subprocess.run(command, check=True, text=True, capture_output=True)
        records = []
        for line in result.stdout.splitlines():
            fields = line.split("|")
            if len(fields) >= 7 and fields[0] in job_ids:
                records.append(dict(zip(("job_id", "state", "elapsed", "total_cpu", "max_rss", "disk_read", "disk_write"), fields[:7])))
        state_by_id = {row["job_id"]: row["state"].split("+", 1)[0] for row in records}
        if all(state_by_id.get(job_id) in terminal for job_id in job_ids):
            return records
        time.sleep(30)
    raise Tier3ValidationError("attributable Slurm jobs remained nonterminal past bounded wait")


def _base_manifest_row(paths: Mapping[str, Path], payload: Mapping[str, Any] | None, run_id: str, now: str) -> dict[str, Any]:
    boundary = payload.get("authorization_boundary", {}) if payload else {}
    decision = payload.get("decision", {}) if payload else {}
    return {
        "run_id": run_id, "generated_at_utc": now, "record_type": "run_summary",
        "gate_decision": decision.get("status", "MISSING_OR_UNREADABLE"),
        "decision_sha256": payload.get("decision_sha256", "") if payload else "",
        "authorization_tuple_digest": boundary.get("authorization_tuple_digest", ""),
        "gate_file_sha256": safe_sha256(paths["gate"]),
        "manifest_sha256": safe_sha256(paths["manifest"]),
        "acquisition_manifest_sha256": safe_sha256(paths["acquisition"]),
        "immutable_inventory_sha256": safe_sha256(paths["inventory"]),
        "root_config_sha256": safe_sha256(paths["root_config"]),
        "input_bundle_digest": boundary.get("input_bundle_digest", ""),
        "root_contract_digest": boundary.get("root_contract_digest", ""),
        "cap_vector_digest": boundary.get("cap_vector_digest", ""),
        "retrieval_digest": boundary.get("retrieval_checksum_obligations_digest", ""),
        "pair_evidence_digest": boundary.get("pair_evidence_digest", ""),
        "measurement_contract_digest": boundary.get("measurement_contract_digest", ""),
        "environment_digest": boundary.get("environment_digest", ""),
        "sweepga_build_sha256": safe_sha256(paths["sweepga"]),
        "impg_handoff_sha256": safe_sha256(paths["impg"]),
        "worker_sha256": safe_sha256(paths["worker"]),
    }


def _refusal_outputs(
    *, paths: Mapping[str, Path], outputs: Mapping[str, Path], payload: Mapping[str, Any] | None,
    run_id: str, now: str, error: Exception,
) -> dict[str, Any]:
    code = classify_failure_code(error)
    base = _base_manifest_row(paths, payload, run_id, now)
    base.update({
        "status": "refused_preflight", "selected_species": 0, "promoted_objects": 0,
        "compressed_input_bytes": 0, "slurm_jobs_submitted": 0, "final_state": "NOT_SUBMITTED",
        "failure_code": code, "failure_source": str(paths["gate"]), "failure_message": str(error),
        "notes": "authorization boundary closed before sbatch; zero compute and zero network use",
    })
    run_rows = [base]
    blockers = payload.get("blockers", []) if payload else []
    for blocker in blockers:
        run_rows.append({**base, "record_type": "gate_blocker", "failure_code": blocker.get("code", ""), "failure_source": blocker.get("source", ""), "failure_message": blocker.get("message", "")})
    telemetry = [{
        "run_id": run_id, "generated_at_utc": now, "record_type": "run_summary", "status": "refused_preflight",
        "retry_index": 0, "final_state": "NOT_SUBMITTED", "elapsed_seconds": 0, "cpu_time_seconds": 0,
        "scratch_peak_gb": 0, "io_read_gb": 0, "io_write_gb": 0, "metadata_operations": 0,
        "network_bytes": 0, "failure_code": code, "notes": "no sbatch command, job id, dependency, sacct use, or compute-node activity",
    }]
    results = [{
        "run_id": run_id, "generated_at_utc": now, "record_type": "run_summary", "status": "excluded",
        "result_scope": "pilot", "metric": "validated_species_count", "numerator": 0, "denominator": 0,
        "target_total": 0, "value": 0, "exclusion_reason": str(error), "failure_code": code,
        "failure_source": str(paths["gate"]), "notes": "no denominator or biological result was imputed",
    }]
    exclusions = [{
        "run_id": run_id, "generated_at_utc": now, "candidate_id": "", "result_scope": "all",
        "metric": "all_results", "status": "excluded", "threshold": "literal repaired GO plus complete exact acquisition",
        "observed": base["gate_decision"], "failure_code": code, "failure_source": str(paths["gate"]),
        "reason": str(error), "imputed": "false", "demographic_input_used": "false",
    }]
    for blocker in blockers:
        exclusions.append({
            "run_id": run_id, "generated_at_utc": now, "candidate_id": "", "result_scope": "authorization_boundary",
            "metric": "preflight_blocker", "status": "excluded", "threshold": "pass", "observed": "fail",
            "failure_code": blocker.get("code", ""), "failure_source": blocker.get("source", ""),
            "reason": blocker.get("message", ""), "imputed": "false", "demographic_input_used": "false",
        })
    refusal = {
        "run_id": run_id, "generated_at_utc": now, "gate_decision": base["gate_decision"], "status": "NOT_SUBMITTED",
        "failure_code": code, "failure_source": str(paths["gate"]), "failure_message": str(error),
        "sbatch_commands_issued": 0, "slurm_jobs_submitted": 0, "compute_jobs_started": 0, "core_seconds": 0,
        "scratch_bytes": 0, "io_read_bytes": 0, "io_write_bytes": 0, "network_bytes": 0, "provider_requests": 0,
        "full_catalog_downloaded": "false", "population_bulk_downloaded": "false", "demographic_inferences": 0,
    }
    refusal["evidence_sha256"] = gate.sha256_json({k: v for k, v in refusal.items() if k != "evidence_sha256"})
    atomic_write_tsv(outputs["run_manifest"], RUN_MANIFEST_FIELDS, run_rows)
    atomic_write_tsv(outputs["telemetry"], SLURM_TELEMETRY_FIELDS, telemetry)
    atomic_write_tsv(outputs["results"], RESULT_FIELDS, results)
    atomic_write_tsv(outputs["exclusions"], EXCLUSION_FIELDS, exclusions)
    atomic_write_tsv(outputs["refusals"], REFUSAL_FIELDS, [refusal])
    _atomic_write(outputs["report"], render_refusal_report(base, refusal, paths, outputs, blockers))
    return {
        "run_id": run_id, "status": "refused_preflight", "final_state": "NOT_SUBMITTED", "failure_code": code,
        "slurm_jobs_submitted": 0, "compute_jobs_started": 0, "core_seconds": 0, "network_bytes": 0,
        **{key: str(value) for key, value in outputs.items()},
    }


def render_refusal_report(base: Mapping[str, Any], refusal: Mapping[str, Any], paths: Mapping[str, Path], outputs: Mapping[str, Path], blockers: Sequence[Mapping[str, Any]]) -> str:
    blocker_lines = "\n".join(f"- `{b.get('code', '')}` — {b.get('message', '')}" for b in blockers) or "- None copied from a readable gate."
    digest_lines = "\n".join(f"- `{key}`: `{base.get(key, '')}`" for key in (
        "gate_file_sha256", "decision_sha256", "authorization_tuple_digest", "manifest_sha256",
        "acquisition_manifest_sha256", "immutable_inventory_sha256", "root_config_sha256", "input_bundle_digest",
        "root_contract_digest", "cap_vector_digest", "retrieval_digest", "pair_evidence_digest",
        "measurement_contract_digest", "environment_digest", "sweepga_build_sha256", "impg_handoff_sha256", "worker_sha256",
    ))
    output_lines = "\n".join(f"- `{key}`: `{path}` (`{safe_sha256(path)}`)" for key, path in outputs.items() if key != "report")
    return f"""# Repaired VGP pilot run report

- Run ID: `{base['run_id']}`
- Gate decision: `{base['gate_decision']}`
- Final state: `NOT_SUBMITTED`
- Failure: `{refusal['failure_code']}` — {refusal['failure_message']}

The repaired authorization boundary closed before any Slurm or provider command. Exactly zero sbatch commands, jobs, compute starts, core-seconds, scratch bytes, I/O bytes, network bytes, provider requests, full-catalog downloads, population-bulk downloads, and demographic inferences occurred. No callable/queryable denominator or biological result was imputed.

## Bound and observed digests

{digest_lines}

## Gate blockers

{blocker_lines}

## Promoted evidence tables

{output_lines}

The Slurm worker contract remains dormant. If a future independently regenerated exact GO passes every repaired binding and complete acquisition check, it requires pinned GNU Guix, node-local `$SLURM_TMPDIR`, no compute-node network, SweepGA whole-haplotype `--num-mappings 1:1`, IMPG native partition/query VCF/BCF contracts, atomic promotion, bounded retries, sentinels, dependency records, `sacct`, and scratch/I/O telemetry. VGP H1/H2 is never population or demographic input.
"""


def run(
    *, gate_path: Path = DEFAULT_GATE, manifest_path: Path = DEFAULT_MANIFEST,
    root_config_path: Path = DEFAULT_ROOT_CONFIG, acquisition_manifest_path: Path = DEFAULT_ACQUISITION,
    inventory_path: Path = DEFAULT_INVENTORY, sweepga_build_path: Path = DEFAULT_SWEEPGA_BUILD,
    impg_handoff_path: Path = DEFAULT_IMPG_HANDOFF, worker_path: Path = DEFAULT_WORKER,
    output_run_manifest_path: Path = DEFAULT_OUTPUT_RUN_MANIFEST,
    output_slurm_telemetry_path: Path = DEFAULT_OUTPUT_SLURM_TELEMETRY,
    output_results_path: Path = DEFAULT_OUTPUT_RESULTS, output_exclusions_path: Path = DEFAULT_OUTPUT_EXCLUSIONS,
    output_refusals_path: Path = DEFAULT_OUTPUT_REFUSALS, output_report_path: Path = DEFAULT_OUTPUT_REPORT,
    submitter: Submitter | None = None,
) -> dict[str, Any]:
    run_id, now = make_run_id(), utc_now()
    # Backward-compatible callers historically supplied only the three primary
    # output paths.  Keep every newly added ledger in that caller's output
    # directory instead of mutating promoted repository artifacts during a
    # temporary review rerun.
    if output_run_manifest_path != DEFAULT_OUTPUT_RUN_MANIFEST:
        if output_exclusions_path == DEFAULT_OUTPUT_EXCLUSIONS:
            output_exclusions_path = output_run_manifest_path.parent / DEFAULT_OUTPUT_EXCLUSIONS.name
        if output_refusals_path == DEFAULT_OUTPUT_REFUSALS:
            output_refusals_path = output_run_manifest_path.parent / DEFAULT_OUTPUT_REFUSALS.name
        if output_report_path == DEFAULT_OUTPUT_REPORT:
            output_report_path = output_run_manifest_path.parent / DEFAULT_OUTPUT_REPORT.name
    paths = {"gate": gate_path, "manifest": manifest_path, "root_config": root_config_path, "acquisition": acquisition_manifest_path, "inventory": inventory_path, "sweepga": sweepga_build_path, "impg": impg_handoff_path, "worker": worker_path}
    outputs = {"run_manifest": output_run_manifest_path, "telemetry": output_slurm_telemetry_path, "results": output_results_path, "exclusions": output_exclusions_path, "refusals": output_refusals_path, "report": output_report_path}
    payload: Mapping[str, Any] | None = None
    try:
        payload = read_json_object(gate_path, "repaired VGP gate")
        acquisition.verify_gate_internal_bindings(payload)
        acquisition.verify_cap_contract(payload)
        decision = payload.get("decision", {}).get("status")
        if decision != "GO":
            if decision == "NO_GO":
                raise Tier3ValidationError("gate decision is NO_GO; repaired VGP compute is not authorized")
            raise Tier3ValidationError(f"unknown gate decision {decision!r}; decision is not the literal GO token")
        authorization = gate.authorize_gate_action(gate_path, manifest_path, root_config_path, "compute")
        for key in FORBIDDEN_SCOPE_KEYS:
            if payload.get("authorization_scope", {}).get(key) is not False:
                raise Tier3ValidationError(f"cap vector scope weakened: {key} must remain false")
        expected_assets = acquisition._selected_assets(payload)
        selected = sorted(payload["selection_audit"]["independently_expected_selected_candidate_ids"])
        audit_pair_contract(payload, selected)
        audit_sweepga_origin_build(sweepga_build_path)
        audit_impg_handoff(impg_handoff_path)
        acquired = audit_acquisition(payload, expected_assets, acquisition_manifest_path, inventory_path)
        if acquired["candidate_ids"] != selected:
            raise Tier3ValidationError("acquisition is incomplete: promoted candidate set differs from exact repaired GO")
        plan = derive_resource_plan(payload, selected)
        if acquired["compressed_input_bytes"] > int(float(payload["cap_vector"]["dimensions"]["compressed_inputs_gib"]["limit"]) * GIB):
            raise Tier3ValidationError("cap vector compressed-input limit violated by promoted objects")
        if safe_sha256(worker_path) == "":
            raise Tier3ValidationError("Slurm worker contract unavailable")
    except (Tier3ValidationError, FileNotFoundError, KeyError, TypeError, ValueError, OSError) as error:
        return _refusal_outputs(paths=paths, outputs=outputs, payload=payload, run_id=run_id, now=now, error=error)

    # This branch is unreachable for the committed NO_GO.  It is intentionally
    # explicit and injectable so tests can prove that preflight failures never
    # reach a submission function.
    submit = submitter or _default_submitter
    array_spec = f"0-{len(selected) - 1}%{plan['concurrency']}"
    argv = [
        "sbatch", "--parsable", f"--array={array_spec}", f"--cpus-per-task={plan['cpus']}",
        f"--mem={int(plan['memory_gib'])}G", f"--time={max(1, int(plan['wall_hours'] * 60))}",
        f"--tmp={max(1, int(plan['scratch_gib']))}G", "--export=NONE", str(worker_path),
        "--run-id", run_id, "--gate-sha256", safe_sha256(gate_path), "--authorization-tuple", authorization["authorization_tuple_digest"],
        "--acquisition-manifest", str(acquisition_manifest_path), "--inventory", str(inventory_path),
    ]
    job_id = submit(argv)
    if submitter is not None:
        # An injected submitter is a planning/test boundary and must not query
        # the live cluster.  Production waits for sacct below.
        return {"run_id": run_id, "status": "submitted_test_boundary", "slurm_jobs_submitted": 1, "job_id": job_id, "sbatch_argv": argv}
    records = _wait_terminal([job_id])
    return {"run_id": run_id, "status": "terminal", "slurm_jobs_submitted": 1, "job_id": job_id, "sacct": records}


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gate", type=Path, default=DEFAULT_GATE)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--root-config", type=Path, default=DEFAULT_ROOT_CONFIG)
    parser.add_argument("--acquisition-manifest", type=Path, default=DEFAULT_ACQUISITION)
    parser.add_argument("--inventory", type=Path, default=DEFAULT_INVENTORY)
    parser.add_argument("--sweepga-build", type=Path, default=DEFAULT_SWEEPGA_BUILD)
    parser.add_argument("--impg-handoff", type=Path, default=DEFAULT_IMPG_HANDOFF)
    parser.add_argument("--worker", type=Path, default=DEFAULT_WORKER)
    parser.add_argument("--run-manifest-out", type=Path, default=DEFAULT_OUTPUT_RUN_MANIFEST)
    parser.add_argument("--slurm-telemetry-out", type=Path, default=DEFAULT_OUTPUT_SLURM_TELEMETRY)
    parser.add_argument("--results-out", type=Path, default=DEFAULT_OUTPUT_RESULTS)
    parser.add_argument("--exclusions-out", type=Path, default=DEFAULT_OUTPUT_EXCLUSIONS)
    parser.add_argument("--refusals-out", type=Path, default=DEFAULT_OUTPUT_REFUSALS)
    parser.add_argument("--report-out", type=Path, default=DEFAULT_OUTPUT_REPORT)
    args = parser.parse_args(argv)
    result = run(
        gate_path=args.gate, manifest_path=args.manifest, root_config_path=args.root_config,
        acquisition_manifest_path=args.acquisition_manifest, inventory_path=args.inventory,
        sweepga_build_path=args.sweepga_build, impg_handoff_path=args.impg_handoff, worker_path=args.worker,
        output_run_manifest_path=args.run_manifest_out, output_slurm_telemetry_path=args.slurm_telemetry_out,
        output_results_path=args.results_out, output_exclusions_path=args.exclusions_out,
        output_refusals_path=args.refusals_out, output_report_path=args.report_out,
    )
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
