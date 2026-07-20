#!/usr/bin/env python3
"""Build and preflight the executable VGP ten-pair authorization packet.

This module supersedes (without rewriting) historical refusal artifacts.  It
never launches biological computation.  Generation reconciles the immutable
roster and CAS inventory; preflight rehashes and gzip-validates the twenty
assembly inputs, authenticates the pinned tool capture, checks storage, and
asks Slurm to validate the exact canary packets with ``sbatch --test-only``.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from analysis.vgp_10_pilot import (
    PilotError,
    canonical_json,
    parse_fasta,
    sequence_dictionary,
    sha256_file,
    verify_environment_capture,
    write_bed,
    Interval,
)


ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = ROOT / "analysis"
PRIMARY = ANALYSIS / "vgp_10_pair_manifest.tsv"
ACQUISITION = ANALYSIS / "vgp_10_pilot_acquisition_manifest.tsv"
INVENTORY = ANALYSIS / "vgp_10_pilot_object_inventory.tsv"
CAPTURE = ANALYSIS / "guix/vgp_10_pilot/realization.json"
AUTHORIZATION = ANALYSIS / "vgp_pilot_authorization_v2.json"
PREFLIGHT = ANALYSIS / "vgp_pilot_authorization_preflight_v2.json"
PACKET_DIR = ANALYSIS / "slurm/vgp_10_pilot/authorized/v2.0.0"
DATA_ROOT_CONFIG = ANALYSIS / "vgp_data_root_config.json"
DATA_ROOT_MIGRATION = ANALYSIS / "vgp_data_root_migration_v1.json"
AUTHORIZATION_ID = "vgp10-auth-20260718-v2"
SCHEMA_VERSION = "vgp-pilot-execution-authorization-v2.0.0"
GENERATED_AT = "2026-07-18T17:00:00Z"
CANARY = "P07"
GIB = 1024**3
CANONICAL_DATA_ROOT = Path("/moosefs/erikg/vgp")

HARD_PRE_EXECUTION_GATES = (
    "exact_same_individual_h1_h2_pair_and_accession_provenance",
    "accepted_cas_input_sha256_and_byte_size",
    "readable_mutually_comparable_assemblies",
    "accepted_sweepga_and_impg_executable_digests",
    "writable_storage_with_required_headroom",
    "valid_slurm_controller_partition_and_packet",
)
HARD_RESULT_GATES = (
    "sweepga_alignment_success",
    "retained_query_multiplicity_lte_1",
    "retained_target_multiplicity_lte_1",
    "every_vcf_ref_matches_h1",
    "manifest_bound_h2_reconstruction",
    "callable_plus_ordered_mask_reasons_equals_h1_universe_exactly",
    "callable_bp_gte_100000000",
    "callable_fraction_gte_0.60",
)
CONFIDENCE_COVARIATES = (
    "final_sequence_qv",
    "busco_completeness_missingness_and_duplication",
    "kmer_and_copy_number_audit",
    "standalone_repeat_report",
    "exact_raw_read_chemistry",
    "raw_read_validation",
    "annotation",
    "independent_ne_evidence",
)
PSMC_SENSITIVITY_MUTATION_RATES = (5e-9, 1e-8, 2e-8)
PSMC_SENSITIVITY_GENERATION_YEARS = (1.0, 2.0, 4.0)
PSMC_SENSITIVITY_SOURCE = "predeclared_generic_sensitivity_grid_not_species_calibration"


class AuthorizationError(RuntimeError):
    """The current packet cannot be authorized or preflighted."""


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_psmc_scaling_scenarios(path: Path) -> None:
    """Freeze a labeled factorial sensitivity grid without implying calibration."""

    fields = (
        "scenario_id", "mutation_rate_per_generation", "generation_time_years",
        "mutation_rate_source", "generation_time_source",
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_suffix(path.suffix + ".partial")
    with partial.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for mutation_rate in PSMC_SENSITIVITY_MUTATION_RATES:
            for generation_years in PSMC_SENSITIVITY_GENERATION_YEARS:
                writer.writerow({
                    "scenario_id": (
                        f"SENS_MU{mutation_rate:.1E}_G{generation_years:g}Y"
                        .replace("+", "")
                    ),
                    "mutation_rate_per_generation": f"{mutation_rate:.8g}",
                    "generation_time_years": f"{generation_years:g}",
                    "mutation_rate_source": PSMC_SENSITIVITY_SOURCE,
                    "generation_time_source": PSMC_SENSITIVITY_SOURCE,
                })
    partial.replace(path)


def _one(rows: Iterable[dict[str, str]], **keys: str) -> dict[str, str]:
    matches = [row for row in rows if all(row.get(key) == value for key, value in keys.items())]
    if len(matches) != 1:
        raise AuthorizationError(f"row must resolve exactly once: {keys}")
    return matches[0]


def _resource_packet(row: Mapping[str, str], h1: Mapping[str, str], h2: Mapping[str, str]) -> dict[str, Any]:
    compressed = int(h1["observed_bytes"]) + int(h2["observed_bytes"])
    assembly_bp = int(row["h1_length_bp"]) + int(row["h2_length_bp"])
    contigs = int(row["h1_contigs"]) + int(row["h2_contigs"])
    scratch = max(64 * GIB, compressed * 16)
    return {
        "measurement_basis": {
            "compressed_input_bytes": compressed,
            "assembly_sequence_bp": assembly_bp,
            "assembly_contigs": contigs,
            "source": "immutable CAS object sizes plus frozen NCBI assembly metadata",
        },
        "stages": {
            "preflight": {"cpus_per_task": 4, "slurm_time": "01:00:00", "slurm_mem": "16G", "scratch_bytes_high": compressed * 3},
            "mapping": {"cpus_per_task": 32, "slurm_time": "2-00:00:00", "slurm_mem": "128G", "scratch_bytes_high": scratch},
            "impg": {"cpus_per_task": 32, "slurm_time": "2-00:00:00", "slurm_mem": "128G", "scratch_bytes_high": scratch},
            "variants": {"cpus_per_task": 16, "slurm_time": "12:00:00", "slurm_mem": "64G", "scratch_bytes_high": scratch},
            "consensus": {"cpus_per_task": 16, "slurm_time": "1-00:00:00", "slurm_mem": "128G", "scratch_bytes_high": scratch},
            "psmc": {"cpus_per_task": 4, "slurm_time": "08:00:00", "slurm_mem": "16G", "scratch_bytes_high": max(16 * GIB, compressed * 2)},
            "psmc_finalize": {"cpus_per_task": 4, "slurm_time": "02:00:00", "slurm_mem": "16G", "scratch_bytes_high": 16 * GIB},
            "annotation": {"cpus_per_task": 8, "slurm_time": "08:00:00", "slurm_mem": "32G", "scratch_bytes_high": max(16 * GIB, compressed * 2)},
        },
        "whole_canary_packet": {
            "cpus_per_task": 32,
            "slurm_time": "3-00:00:00",
            "initial_memory_gib": 128,
            "node_local_scratch_bytes_high": scratch,
            "oom_retry_memory_gib": [256, 512],
        },
    }


def build_authorization(
    primary_path: Path = PRIMARY,
    acquisition_path: Path = ACQUISITION,
    inventory_path: Path = INVENTORY,
    capture_path: Path = CAPTURE,
) -> dict[str, Any]:
    primaries = read_tsv(primary_path)
    acquisitions = read_tsv(acquisition_path)
    inventory = read_tsv(inventory_path)
    if [row["selection_id"] for row in primaries] != [f"P{i:02d}" for i in range(1, 11)]:
        raise AuthorizationError("authorization requires exactly the frozen P01..P10 primary roster")
    capture = json.loads(capture_path.read_text(encoding="utf-8"))
    tools = {row["name"]: row for row in capture["executables"]}
    for required in ("sweepga", "impg"):
        if required not in tools:
            raise AuthorizationError(f"captured executable absent: {required}")
    pairs: list[dict[str, Any]] = []
    for design in primaries:
        sid = design["selection_id"]
        acquired = _one(acquisitions, selection_id=sid, roster_type="primary")
        if acquired["activation_status"] != "active_primary":
            raise AuthorizationError(f"primary is not approved: {sid}")
        if acquired["core_identity_status"] != "pass":
            raise AuthorizationError(f"exact-pair identity is not verified: {sid}")
        if acquired["core_acquisition_status"] != "verified":
            raise AuthorizationError(f"core CAS acquisition is incomplete: {sid}")
        for field in ("biosample", "individual_or_isolate", "h1_accession_version", "h2_accession_version"):
            if acquired[field] != design[field]:
                raise AuthorizationError(f"design/acquisition provenance drift for {sid}.{field}")
        h1 = _one(inventory, selection_id=sid, roster_type="primary", side="h1", object_role="genome_fasta")
        h2 = _one(inventory, selection_id=sid, roster_type="primary", side="h2", object_role="genome_fasta")
        for side, obj in (("h1", h1), ("h2", h2)):
            if obj["status"] not in {"reused", "verified"}:
                raise AuthorizationError(f"accepted {side} assembly absent from CAS: {sid}")
            if obj["local_sha256"] != obj["expected_local_sha256"] or len(obj["local_sha256"]) != 64:
                raise AuthorizationError(f"accepted {side} digest mismatch: {sid}")
            if int(obj["observed_bytes"]) <= 0 or obj["observed_bytes"] != obj["expected_bytes"]:
                raise AuthorizationError(f"accepted {side} size mismatch: {sid}")
        ratio = int(design["h1_length_bp"]) / int(design["h2_length_bp"])
        if not 0.8 <= ratio <= 1.25:
            raise AuthorizationError(f"assembly spans are not mutually comparable: {sid}")
        missing_covariates = [
            "standalone_repeat_report",
            "final_sequence_qv",
            "raw_read_validation",
            "kmer_and_copy_number_audit",
            "independent_ne_evidence",
        ]
        if acquired["completeness_status"] == "missing":
            missing_covariates.append("busco_completeness_missingness_and_duplication")
        if sid == "P09":
            missing_covariates.append("exact_raw_read_chemistry")
        annotation_status = acquired["annotation_branch_status"]
        pair = {
            "selection_id": sid,
            "species": design["species"],
            "biosample": design["biosample"],
            "individual_or_isolate": design["individual_or_isolate"],
            "h1_accession_version": design["h1_accession_version"],
            "h2_accession_version": design["h2_accession_version"],
            "h1": _object_record(h1),
            "h2": _object_record(h2),
            "assembly_span_ratio_h1_over_h2": ratio,
            "core_authorization": "AUTHORIZED",
            "authorized_outputs": ["core_diversity", "unscaled_psmc"],
            "annotation_authorization": "CONDITIONAL_EXACT_BINDING" if "available" in annotation_status else "NOT_AUTHORIZED_ANNOTATION_ONLY",
            "annotation_absence_blocks_core": False,
            "confidence_tier_at_authorization": "C",
            "missing_confidence_covariates": sorted(set(missing_covariates)),
            "missing_covariates_block_core": False,
            "sequence_derived_repeat_or_low_complexity_mask": "GENERATE_DURING_MAPPING",
        }
        pair["resources"] = _resource_packet(design, h1, h2)
        pairs.append(pair)
    canary = min(pairs, key=lambda row: (row["resources"]["measurement_basis"]["compressed_input_bytes"], row["selection_id"]))
    if canary["selection_id"] != CANARY:
        raise AuthorizationError(f"immutable size ordering changed canary: {canary['selection_id']}")
    return {
        "schema_version": SCHEMA_VERSION,
        "authorization_id": AUTHORIZATION_ID,
        "generated_at_utc": GENERATED_AT,
        "supersedes_policy_not_history": "vgp10-review-20260718-v1 bounded-repair refusal policy",
        "historical_artifacts_mutated": False,
        "authorization_scope": "ten frozen primary same-individual H1/H2 pairs",
        "canonical_data_root": str(CANONICAL_DATA_ROOT),
        "legacy_data_root_is_migration_input_only": True,
        "hard_pre_execution_gates": list(HARD_PRE_EXECUTION_GATES),
        "hard_per_pair_result_gates": list(HARD_RESULT_GATES),
        "confidence_covariates_not_authorization_gates": list(CONFIDENCE_COVARIATES),
        "missing_optional_qc_global_job_count_effect": 0,
        "authorized_pair_count": len(pairs),
        "authorized_selection_ids": [row["selection_id"] for row in pairs],
        "pairs": pairs,
        "canary": {
            "selection_id": canary["selection_id"],
            "species": canary["species"],
            "selection_rule": "minimum immutable summed CAS compressed assembly bytes; selection_id lexical tie-break",
            "measured_compressed_input_bytes": canary["resources"]["measurement_basis"]["compressed_input_bytes"],
            "job_packets": {
                "initial_128_gib": "analysis/slurm/vgp_10_pilot/authorized/v2.0.0/P07.canary.sbatch",
                "oom_retry_256_gib": "analysis/slurm/vgp_10_pilot/authorized/v2.0.0/P07.canary.mem256.sbatch",
                "oom_retry_512_gib": "analysis/slurm/vgp_10_pilot/authorized/v2.0.0/P07.canary.mem512.sbatch",
            },
            "worker": "analysis/slurm/vgp_10_pilot/authorized/v2.0.0/run_canary.sh",
        },
        "environment": {
            "guix_channel_commit": capture["channel_commit"],
            "capture": str(capture_path.relative_to(ROOT)),
            "capture_sha256": sha256_file(capture_path),
            "sweepga": tools["sweepga"],
            "impg": tools["impg"],
        },
        "source_artifacts": {
            str(path.relative_to(ROOT)): sha256_file(path)
            for path in (
                primary_path, acquisition_path, inventory_path, capture_path,
                DATA_ROOT_CONFIG, DATA_ROOT_MIGRATION,
            )
        },
        "execution_requires_fresh_preflight": str(PREFLIGHT.relative_to(ROOT)),
        "biological_jobs_submitted_by_authorization": 0,
    }


def _object_record(row: Mapping[str, str]) -> dict[str, Any]:
    digest = row["local_sha256"]
    return {
        "accession_version": row["accession_version"],
        "cas_path": str(CANONICAL_DATA_ROOT / "objects/sha256" / digest[:2] / digest[2:4] / digest),
        "compressed_bytes": int(row["observed_bytes"]),
        "sha256": row["local_sha256"],
        "compression": "gzip",
    }


def validate_authorization(value: Mapping[str, Any], verify_source_digests: bool = True) -> dict[str, Any]:
    if value.get("schema_version") != SCHEMA_VERSION or value.get("authorization_id") != AUTHORIZATION_ID:
        raise AuthorizationError("unexpected authorization version")
    if value.get("authorized_selection_ids") != [f"P{i:02d}" for i in range(1, 11)]:
        raise AuthorizationError("authorization is not exactly P01..P10")
    if value.get("authorized_pair_count") != 10:
        raise AuthorizationError("authorization must contain ten primaries")
    pairs = value.get("pairs", [])
    if len(pairs) != 10 or any(row.get("core_authorization") != "AUTHORIZED" for row in pairs):
        raise AuthorizationError("every primary must carry explicit core authorization")
    if any(row.get("missing_covariates_block_core") is not False for row in pairs):
        raise AuthorizationError("optional confidence covariates cannot veto core")
    if value.get("missing_optional_qc_global_job_count_effect") != 0:
        raise AuthorizationError("optional QC changed the global authorized job count")
    if value.get("canary", {}).get("selection_id") != CANARY:
        raise AuthorizationError("canary is not immutable minimum-size P07")
    if value.get("biological_jobs_submitted_by_authorization") != 0:
        raise AuthorizationError("authorization generator must not submit biology")
    if value.get("canonical_data_root") != str(CANONICAL_DATA_ROOT):
        raise AuthorizationError("active authorization does not use the canonical shared VGP root")
    if value.get("legacy_data_root_is_migration_input_only") is not True:
        raise AuthorizationError("legacy root role is not constrained to migration input")
    if verify_source_digests:
        for relative, expected in value["source_artifacts"].items():
            path = ROOT / relative
            if not path.is_file() or sha256_file(path) != expected:
                raise AuthorizationError(f"authorization source drift: {relative}")
    return {"authorized_pairs": 10, "canary": CANARY, "optional_qc_vetoes": 0}


def _gzip_fasta_readable(path: Path) -> dict[str, Any]:
    gzip_bin = shutil.which("gzip")
    if gzip_bin is not None:
        result = subprocess.run([gzip_bin, "-t", str(path)], check=False, capture_output=True)
        if result.returncode:
            raise AuthorizationError(f"gzip stream readability gate failed: {path}")
    else:
        # The pinned analysis profile deliberately need not contain GNU gzip.
        # Reading to EOF makes Python's gzip layer verify the full DEFLATE
        # stream, CRC, and uncompressed size without retaining its contents.
        try:
            with gzip.open(path, "rb") as compressed:
                while compressed.read(16 * 1024 * 1024):
                    pass
        except (EOFError, OSError) as exc:
            raise AuthorizationError(f"gzip stream readability gate failed: {path}") from exc
    with gzip.open(path, "rt", encoding="ascii") as handle:
        header = handle.readline().rstrip("\r\n")
        first_sequence = handle.readline().strip().upper()
    if not header.startswith(">") or not first_sequence:
        raise AuthorizationError(f"gzip payload is not a nonempty FASTA: {path}")
    if set(first_sequence) - set("ACGTNRYKMSWBDHV"):
        raise AuthorizationError(f"invalid FASTA alphabet in first sequence row: {path}")
    return {"gzip_stream_complete": True, "first_header": header[1:].split()[0]}


def preflight(value: Mapping[str, Any], storage_root: Path, run_scheduler_checks: bool = True) -> dict[str, Any]:
    validate_authorization(value)
    pair_checks = []
    for pair in value["pairs"]:
        checked: dict[str, Any] = {"selection_id": pair["selection_id"]}
        for side in ("h1", "h2"):
            obj = pair[side]
            path = Path(obj["cas_path"])
            if not path.is_file() or path.stat().st_size != obj["compressed_bytes"]:
                raise AuthorizationError(f"CAS file/size gate failed: {pair['selection_id']} {side}")
            if sha256_file(path) != obj["sha256"]:
                raise AuthorizationError(f"CAS digest gate failed: {pair['selection_id']} {side}")
            checked[side] = _gzip_fasta_readable(path)
        observed_ratio = float(pair["assembly_span_ratio_h1_over_h2"])
        if not 0.8 <= observed_ratio <= 1.25:
            raise AuthorizationError(f"observed assemblies are not mutually comparable: {pair['selection_id']}")
        checked["observed_h1_h2_span_ratio"] = observed_ratio
        checked["authorized"] = True
        pair_checks.append(checked)
    capture = verify_environment_capture(ROOT / value["environment"]["capture"])
    required_bytes = sum(
        int(pair["resources"]["whole_canary_packet"]["node_local_scratch_bytes_high"])
        for pair in value["pairs"]
    )
    storage_root = storage_root.resolve()
    usage = shutil.disk_usage(storage_root)
    storage = {
        "path": str(storage_root),
        "writable": os.access(storage_root, os.W_OK),
        "free_bytes": usage.free,
        "ten_pair_conservative_required_bytes": required_bytes,
        "headroom_pass": usage.free >= required_bytes,
    }
    if not storage["writable"] or not storage["headroom_pass"]:
        raise AuthorizationError("writable storage/headroom gate failed")
    scheduler: dict[str, Any] = {"checked": run_scheduler_checks}
    if run_scheduler_checks:
        for name in ("sbatch", "scontrol", "sinfo", "sacct"):
            if shutil.which(name) is None:
                raise AuthorizationError(f"Slurm client absent: {name}")
        ping = subprocess.run(["scontrol", "ping"], check=False, text=True, capture_output=True)
        if ping.returncode or "UP" not in ping.stdout:
            raise AuthorizationError(f"Slurm controller gate failed: {ping.stdout}{ping.stderr}")
        partitions = subprocess.run(
            ["sinfo", "-h", "-o", "%P|%a|%l|%c|%m"], check=True, text=True, capture_output=True
        ).stdout.splitlines()
        if not any("|up|" in row and int(row.rsplit("|", 1)[-1].rstrip("+")) >= 131072 for row in partitions):
            raise AuthorizationError("no up Slurm partition advertises at least 128 GiB/node")
        packet_results = []
        for relative in value["canary"]["job_packets"].values():
            result = subprocess.run(["sbatch", "--test-only", str(ROOT / relative)], check=False, text=True, capture_output=True)
            if result.returncode:
                raise AuthorizationError(f"sbatch --test-only rejected {relative}: {result.stderr or result.stdout}")
            packet_results.append({"packet": relative, "result": (result.stdout or result.stderr).strip()})
        scheduler.update({"controller": ping.stdout.strip(), "partitions": partitions, "packets": packet_results, "pass": True})
    return {
        "schema_version": "vgp-pilot-authorization-preflight-v2.0.0",
        "authorization_id": value["authorization_id"],
        "authorization_payload_sha256": hashlib.sha256(
            canonical_json(value).encode("utf-8")
        ).hexdigest(),
        "performed_at_utc": GENERATED_AT,
        "pair_checks": pair_checks,
        "authorized_selection_ids": [row["selection_id"] for row in pair_checks],
        "authorized_pair_count": len(pair_checks),
        "environment": {
            "profile": capture["profile"],
            "sweepga_sha256": value["environment"]["sweepga"]["sha256"],
            "impg_sha256": value["environment"]["impg"]["sha256"],
            "pass": True,
        },
        "storage": storage,
        "scheduler": scheduler,
        "biological_execution": False,
        "slurm_jobs_submitted": 0,
        "global_authorization_result": "GO_10_OF_10",
    }


def materialize_input(value: Mapping[str, Any], selection_id: str, data_root: Path) -> dict[str, Any]:
    """Stage one authorized pair into node-local scratch for the worker."""
    validate_authorization(value)
    data_root = data_root.resolve()
    matches = [row for row in value["pairs"] if row["selection_id"] == selection_id]
    if len(matches) != 1:
        raise AuthorizationError(f"authorized selection does not resolve once: {selection_id}")
    pair = matches[0]
    input_dir = data_root / "pilot/inputs" / selection_id
    input_dir.mkdir(parents=True, exist_ok=True)
    assets: dict[str, Any] = {}
    sequences: dict[str, dict[str, str]] = {}
    for side in ("h1", "h2"):
        source = Path(pair[side]["cas_path"])
        if sha256_file(source) != pair[side]["sha256"]:
            raise AuthorizationError(f"CAS digest drift while staging {selection_id} {side}")
        destination = input_dir / f"{side}.fa"
        with gzip.open(source, "rb") as incoming, destination.open("wb") as outgoing:
            shutil.copyfileobj(incoming, outgoing, length=16 * 1024 * 1024)
        sequences[side] = parse_fasta(destination)
        assets[f"{side}_fasta"] = {
            "path": str(destination),
            "sha256": sha256_file(destination),
            "size_bytes": destination.stat().st_size,
            "source_cas_sha256": pair[side]["sha256"],
            "sequence_dictionary": sequence_dictionary(sequences[side]),
        }
    h1_universe = [Interval(name, 0, len(seq)) for name, seq in sequences["h1"].items()]
    query_universe = [Interval(name, 0, len(seq)) for name, seq in sequences["h2"].items()]
    write_bed(input_dir / "h1_universe.bed", h1_universe)
    write_bed(input_dir / "eligible_query_regions.bed", h1_universe + query_universe)
    (input_dir / "exclusions").mkdir(exist_ok=True)
    manifest = {
        "canonical_vgp_root": str(data_root),
        "selection_id": selection_id,
        "biosample": pair["biosample"],
        "individual_or_isolate": pair["individual_or_isolate"],
        "h1_accession_version": pair["h1_accession_version"],
        "h2_accession_version": pair["h2_accession_version"],
        "orientation": "H1_reference_H2_query",
        "assets": assets,
        "confidence_covariates": {name: None for name in CONFIDENCE_COVARIATES},
        "selective_validation": {},
        "annotation": None,
        "result_gates": {"minimum_callable_bp": 100_000_000, "minimum_callable_fraction": 0.60},
        "authorization_id": value["authorization_id"],
    }
    (input_dir / "input-manifest.json").write_text(canonical_json(manifest), encoding="utf-8")
    (input_dir / "resources.json").write_text(canonical_json(pair["resources"]), encoding="utf-8")
    write_psmc_scaling_scenarios(input_dir / "psmc_scaling_scenarios.tsv")
    return {"selection_id": selection_id, "input_dir": str(input_dir), "assets": assets}


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    generate = sub.add_parser("generate")
    generate.add_argument("--output", type=Path, default=AUTHORIZATION)
    validate = sub.add_parser("validate")
    validate.add_argument("--authorization", type=Path, default=AUTHORIZATION)
    check = sub.add_parser("preflight")
    check.add_argument("--authorization", type=Path, default=AUTHORIZATION)
    check.add_argument("--storage-root", type=Path, default=CANONICAL_DATA_ROOT)
    check.add_argument("--skip-scheduler", action="store_true")
    check.add_argument("--output", type=Path, default=PREFLIGHT)
    materialize = sub.add_parser("materialize-input")
    materialize.add_argument("--authorization", type=Path, default=AUTHORIZATION)
    materialize.add_argument("--selection-id", required=True)
    materialize.add_argument("--data-root", required=True, type=Path)
    args = parser.parse_args(argv)
    try:
        if args.command == "generate":
            result = build_authorization()
            args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        else:
            value = json.loads(args.authorization.read_text(encoding="utf-8"))
            if args.command == "validate":
                result = validate_authorization(value)
            elif args.command == "preflight":
                result = preflight(value, args.storage_root, not args.skip_scheduler)
                args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
            else:
                result = materialize_input(value, args.selection_id, args.data_root)
        print(json.dumps(result, sort_keys=True))
        return 0
    except (AuthorizationError, PilotError, OSError, ValueError, json.JSONDecodeError) as error:
        print(f"authorization error: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
