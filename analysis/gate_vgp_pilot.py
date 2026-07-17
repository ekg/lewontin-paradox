#!/usr/bin/env python3
"""Recompute and enforce the bounded VGP pilot authorization gate."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping, MutableMapping, Sequence

from analysis import freeze_vgp_manifest as freeze
from analysis.tier3_common import Tier3ValidationError, sha256_file


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = PROJECT_ROOT / "analysis" / "vgp_pilot_manifest.tsv"
DEFAULT_SIZE_BUDGET = PROJECT_ROOT / "analysis" / "vgp_pilot_size_budget.tsv"
DEFAULT_FREEZE_PROVENANCE = PROJECT_ROOT / "analysis" / "vgp_phase1_freeze_provenance.json"
DEFAULT_ROOT_CONFIG = PROJECT_ROOT / "analysis" / "vgp_data_root_config.json"
DEFAULT_ROOT_VALIDATION = PROJECT_ROOT / "analysis" / "vgp_data_root_validation.json"
DEFAULT_DECISIONS = PROJECT_ROOT / "analysis" / "vertebrate_scaleout_decisions.tsv"
DEFAULT_EXECUTION_PLAN = PROJECT_ROOT / "analysis" / "vertebrate_scaleout_execution_plan.md"
DEFAULT_RESOURCE_BUDGET = PROJECT_ROOT / "analysis" / "vertebrate_scaleout_resource_budget.tsv"
DEFAULT_GATE_JSON = PROJECT_ROOT / "analysis" / "vgp_pilot_gate.json"
DEFAULT_GATE_REVIEW = PROJECT_ROOT / "analysis" / "vgp_pilot_gate_review.md"

RELEVANT_BUDGET_ROWS = (
    "stratified_pilot",
    "stage_assets_per_species",
    "exact_preflight_per_species",
    "sweepga_mapping_per_species",
    "impg_index_partition_query_per_species",
    "denominator_summary_per_species",
    "tier3c_composition_per_species",
)

PAIR_MODALITY = freeze.SEED_MODALITY_TIER3A
COMPOSITION_MODALITY = freeze.SEED_MODALITY_TIER3C

EXPECTED_SOURCE_COUNTS = {
    "unique_species": 714,
    "completed": 223,
    "completed_annotated": 120,
    "triple_eligible": 40,
    "triple_eligible_fish": 13,
    "completed_refseq_fish": 46,
}

HEX64 = re.compile(r"^[0-9a-f]{64}$")
HEX32 = re.compile(r"^[0-9a-f]{32}$")
VERSIONED_ACCESSION = re.compile(r"^(GC[AF]_\d+\.\d+)$")
HTTPS_URL = re.compile(r"^https://\S+$")


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def canonical_json_bytes(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("utf-8")


def sha256_json(value: Any) -> str:
    import hashlib

    return hashlib.sha256(canonical_json_bytes(value)).hexdigest()


def write_json(path: Path, value: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def parse_int(value: Any) -> int | None:
    if value is None:
        return None
    text = str(value).strip().replace(",", "")
    if not text or text.upper() in {"NA", "N/A"} or "UNAVAILABLE" in text or "UNRESOLVED" in text:
        return None
    try:
        return int(text)
    except ValueError:
        return None


def parse_float(value: Any) -> float | None:
    if value is None:
        return None
    text = str(value).strip().replace(",", "")
    if not text or text.upper() in {"NA", "N/A"} or "UNAVAILABLE" in text or "UNRESOLVED" in text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def parse_compact_number(text: str) -> float:
    value = text.strip().lower().replace(",", "")
    multiplier = 1.0
    if value.endswith("k"):
        value = value[:-1]
        multiplier = 1000.0
    return float(value) * multiplier


def parse_numeric_phrase(text: str) -> float:
    lowered = text.strip().lower()
    word_map = {
        "zero": 0.0,
        "one": 1.0,
        "two": 2.0,
        "three": 3.0,
        "four": 4.0,
        "five": 5.0,
        "six": 6.0,
        "seven": 7.0,
        "eight": 8.0,
        "nine": 9.0,
        "ten": 10.0,
    }
    if lowered in word_map:
        return word_map[lowered]
    return parse_compact_number(lowered)


def split_modalities(value: str) -> set[str]:
    return {item.strip() for item in value.split(";") if item.strip()}


def is_https_url(value: str) -> bool:
    return bool(value and HTTPS_URL.match(value.strip()))


def has_sha256(value: str) -> bool:
    return bool(value and HEX64.match(value.strip().lower()))


def has_md5(value: str) -> bool:
    return bool(value and HEX32.match(value.strip().lower()))


def versioned_accession(value: str) -> bool:
    return bool(value and VERSIONED_ACCESSION.match(value.strip()))


def blocker(code: str, message: str, source: str) -> dict[str, str]:
    return {"code": code, "message": message, "source": source}


def read_source_catalog(provenance: Mapping[str, Any]) -> dict[str, Any]:
    source_path = Path(provenance["source_catalog"]["path"])
    text = source_path.read_text(encoding="utf-8")
    rows = freeze.parse_vgp_rows(text)
    counts = freeze.compute_counts(rows)
    observed = counts["observed"]
    discrepancies = []
    for key, expected in EXPECTED_SOURCE_COUNTS.items():
        observed_value = observed[key]
        if observed_value != expected:
            discrepancies.append(
                {
                    "metric": key,
                    "expected": expected,
                    "observed": observed_value,
                    "delta": observed_value - expected,
                }
            )
    return {
        "path": str(source_path),
        "sha256": sha256_file(source_path),
        "line_count": text.count("\n"),
        "counts": counts,
        "discrepancies": discrepancies,
    }


def audit_manifest_row(row: Mapping[str, str]) -> dict[str, Any]:
    issues: list[dict[str, str]] = []
    candidate_id = row["candidate_id"]
    source_ref = f"analysis/vgp_pilot_manifest.tsv:{candidate_id}"
    modalities = split_modalities(row.get("seed_modalities", ""))
    pair_required = PAIR_MODALITY in modalities

    def add(code: str, message: str) -> None:
        issues.append({"code": code, "message": message})

    h1 = row.get("h1_accession_version", "")
    if not versioned_accession(h1):
        add("H1_ACCESSION_INVALID", "missing or unversioned H1 accession")
    if not is_https_url(row.get("h1_fasta_url", "")):
        add("H1_URL_MISSING", "missing exact H1 FASTA URL")
    if not has_md5(row.get("h1_provider_md5", "")):
        add("H1_MD5_MISSING", "missing advertised H1 MD5")
    if not has_sha256(row.get("h1_fasta_sha256", "")):
        add("H1_SHA256_MISSING", "missing H1 SHA-256")
    if parse_int(row.get("h1_fasta_compressed_bytes")) is None:
        add("H1_SIZE_MISSING", "missing H1 compressed size")

    if pair_required:
        if not versioned_accession(row.get("h2_accession_version", "")):
            add("H2_ACCESSION_INVALID", "missing or unversioned H2 accession for paired row")
        if not is_https_url(row.get("h2_fasta_url", "")):
            add("H2_URL_MISSING", "missing exact H2 FASTA URL for paired row")
        if not has_md5(row.get("h2_provider_md5", "")):
            add("H2_MD5_MISSING", "missing advertised H2 MD5 for paired row")
        if not has_sha256(row.get("h2_fasta_sha256", "")):
            add("H2_SHA256_MISSING", "missing H2 SHA-256 for paired row")
        if parse_int(row.get("h2_fasta_compressed_bytes")) is None:
            add("H2_SIZE_MISSING", "missing H2 compressed size for paired row")
        if row.get("same_individual_status") != "yes":
            add("PAIR_NOT_SAME_INDIVIDUAL", "paired row lacks exact same-individual evidence")
        if not is_https_url(row.get("pair_evidence_url", "")):
            add("PAIR_EVIDENCE_URL_MISSING", "paired row lacks exact pair evidence URL")

    annotation_native = row.get("annotation_native_status", "").strip().lower()
    if not is_https_url(row.get("annotation_gff_url", "")):
        add("ANNOTATION_URL_MISSING", "missing annotation GFF URL")
    if not has_md5(row.get("annotation_provider_md5", "")):
        add("ANNOTATION_MD5_MISSING", "missing annotation MD5")
    if not has_sha256(row.get("annotation_gff_sha256", "")):
        add("ANNOTATION_SHA256_MISSING", "missing annotation SHA-256")
    if parse_int(row.get("annotation_gff_compressed_bytes")) is None:
        add("ANNOTATION_SIZE_MISSING", "missing annotation compressed size")
    if (
        not annotation_native
        or "unresolved" in annotation_native
        or "missing" in annotation_native
        or "project" in annotation_native
        or "lift" in annotation_native
    ):
        add("ANNOTATION_NOT_NATIVE", "annotation is not proven native to the exact H1")
    if row.get("annotation_reference_accession_version") != h1:
        add("ANNOTATION_REFERENCE_MISMATCH", "annotation reference does not exactly match H1 accession")

    taxid = parse_int(row.get("ncbi_taxid"))
    if taxid is None or taxid <= 0:
        add("TAXID_INVALID", "missing or invalid NCBI taxonomy identifier")

    callable_bases = parse_int(row.get("callable_bases"))
    queryable_gene_count = parse_int(row.get("queryable_gene_count"))
    queryable_gene_bases = parse_int(row.get("queryable_gene_bases"))
    if callable_bases is None or callable_bases <= 0:
        add("CALLABLE_BASES_UNRESOLVED", "callable-base denominator is unresolved or non-positive")
    if queryable_gene_count is None or queryable_gene_count <= 0:
        add("QUERYABLE_GENE_COUNT_UNRESOLVED", "queryable gene-count denominator is unresolved or non-positive")
    if queryable_gene_bases is None or queryable_gene_bases <= 0:
        add("QUERYABLE_GENE_BASES_UNRESOLVED", "queryable CDS-base denominator is unresolved or non-positive")
    if row.get("callability_reference_accession_version") != h1:
        add("CALLABILITY_REFERENCE_MISMATCH", "callability reference does not exactly match H1 accession")
    variant_reference = row.get("variant_reference_accession_version", "")
    if variant_reference and variant_reference != h1:
        add("VARIANT_REFERENCE_MISMATCH", "variant reference does not exactly match H1 accession")

    size_sum = sum(
        parsed
        for parsed in (
            parse_int(row.get("h1_fasta_compressed_bytes")),
            parse_int(row.get("h2_fasta_compressed_bytes")) if pair_required else 0,
            parse_int(row.get("annotation_gff_compressed_bytes")),
        )
        if parsed is not None
    )
    declared_download = parse_int(row.get("predicted_download_bytes_exact"))
    if declared_download is None:
        add("DECLARED_DOWNLOAD_SIZE_MISSING", "predicted exact download bytes are missing")
    elif declared_download != size_sum:
        add("DECLARED_DOWNLOAD_SIZE_MISMATCH", "declared exact download bytes do not equal advertised asset sizes")

    composition_ready = not any(
        issue["code"]
        in {
            "H1_ACCESSION_INVALID",
            "H1_URL_MISSING",
            "H1_MD5_MISSING",
            "H1_SHA256_MISSING",
            "H1_SIZE_MISSING",
            "ANNOTATION_URL_MISSING",
            "ANNOTATION_MD5_MISSING",
            "ANNOTATION_SHA256_MISSING",
            "ANNOTATION_SIZE_MISSING",
            "ANNOTATION_NOT_NATIVE",
            "ANNOTATION_REFERENCE_MISMATCH",
            "TAXID_INVALID",
            "CALLABLE_BASES_UNRESOLVED",
            "QUERYABLE_GENE_COUNT_UNRESOLVED",
            "QUERYABLE_GENE_BASES_UNRESOLVED",
            "CALLABILITY_REFERENCE_MISMATCH",
            "VARIANT_REFERENCE_MISMATCH",
            "DECLARED_DOWNLOAD_SIZE_MISSING",
            "DECLARED_DOWNLOAD_SIZE_MISMATCH",
        }
        for issue in issues
    )
    diversity_ready = composition_ready and pair_required and not any(
        issue["code"]
        in {
            "H2_ACCESSION_INVALID",
            "H2_URL_MISSING",
            "H2_MD5_MISSING",
            "H2_SHA256_MISSING",
            "H2_SIZE_MISSING",
            "PAIR_NOT_SAME_INDIVIDUAL",
            "PAIR_EVIDENCE_URL_MISSING",
        }
        for issue in issues
    )
    return {
        "candidate_id": candidate_id,
        "pair_required": pair_required,
        "selected": row.get("pilot_selected") == "yes",
        "composition_ready": composition_ready,
        "diversity_ready": diversity_ready,
        "issues": issues,
        "source": source_ref,
    }


def manifest_budget_consistency(
    manifest_rows: Sequence[Mapping[str, str]],
    size_budget_rows: Sequence[Mapping[str, str]],
) -> tuple[list[dict[str, str]], dict[str, Any]]:
    blockers: list[dict[str, str]] = []
    candidate_budget = {row["candidate_id"]: row for row in size_budget_rows if row["row_type"] == "candidate"}
    aggregate_selected = next((row for row in size_budget_rows if row["row_type"] == "aggregate_selected"), None)
    if aggregate_selected is None:
        blockers.append(
            blocker(
                "SIZE_BUDGET_AGGREGATE_SELECTED_MISSING",
                "analysis/vgp_pilot_size_budget.tsv lacks the aggregate_selected row",
                "analysis/vgp_pilot_size_budget.tsv",
            )
        )
        aggregate_selected = {
            "download_bytes_exact": "",
            "persistent_storage_bytes_exact": "",
            "core_hours_high": "",
            "scratch_gb_high": "",
            "inode_count_high": "",
            "moosefs_read_gb_high": "",
            "moosefs_write_gb_high": "",
            "metadata_operations_high": "",
        }

    selected_rows = [row for row in manifest_rows if row.get("pilot_selected") == "yes"]
    selected_download = 0
    selected_persistent = 0
    for row in manifest_rows:
        budget_row = candidate_budget.get(row["candidate_id"])
        if budget_row is None:
            blockers.append(
                blocker(
                    "SIZE_BUDGET_ROW_MISSING",
                    f"analysis/vgp_pilot_size_budget.tsv lacks candidate row {row['candidate_id']}",
                    "analysis/vgp_pilot_size_budget.tsv",
                )
            )
            continue
        exact_download = sum(
            parsed
            for parsed in (
                parse_int(row.get("h1_fasta_compressed_bytes")),
                parse_int(row.get("h2_fasta_compressed_bytes")),
                parse_int(row.get("annotation_gff_compressed_bytes")),
            )
            if parsed is not None
        )
        budget_download = parse_int(budget_row.get("download_bytes_exact"))
        budget_persistent = parse_int(budget_row.get("persistent_storage_bytes_exact"))
        if budget_download is None:
            blockers.append(
                blocker(
                    "SIZE_BUDGET_MISSING_DOWNLOAD_BYTES",
                    f"candidate {row['candidate_id']} is missing exact download bytes in the size budget",
                    "analysis/vgp_pilot_size_budget.tsv",
                )
            )
        elif budget_download != exact_download:
            blockers.append(
                blocker(
                    "SIZE_BUDGET_DOWNLOAD_MISMATCH",
                    f"candidate {row['candidate_id']} exact download bytes disagree with advertised asset sizes",
                    "analysis/vgp_pilot_size_budget.tsv",
                )
            )
        if budget_persistent is None:
            blockers.append(
                blocker(
                    "SIZE_BUDGET_MISSING_PERSISTENT_BYTES",
                    f"candidate {row['candidate_id']} is missing exact persistent bytes in the size budget",
                    "analysis/vgp_pilot_size_budget.tsv",
                )
            )
        if row.get("pilot_selected") == "yes":
            selected_download += budget_download or 0
            selected_persistent += budget_persistent or 0

    aggregate_download = parse_int(aggregate_selected.get("download_bytes_exact"))
    if aggregate_download is None:
        blockers.append(
            blocker(
                "SIZE_BUDGET_AGGREGATE_DOWNLOAD_MISSING",
                "aggregate_selected exact download bytes are missing",
                "analysis/vgp_pilot_size_budget.tsv",
            )
        )
    elif aggregate_download != selected_download:
        blockers.append(
            blocker(
                "SIZE_BUDGET_AGGREGATE_DOWNLOAD_MISMATCH",
                "aggregate_selected exact download bytes do not match the sum of selected rows",
                "analysis/vgp_pilot_size_budget.tsv",
            )
        )

    aggregate_persistent = parse_int(aggregate_selected.get("persistent_storage_bytes_exact"))
    if aggregate_persistent is None:
        blockers.append(
            blocker(
                "SIZE_BUDGET_AGGREGATE_PERSISTENT_MISSING",
                "aggregate_selected exact persistent bytes are missing",
                "analysis/vgp_pilot_size_budget.tsv",
            )
        )
    elif aggregate_persistent != selected_persistent:
        blockers.append(
            blocker(
                "SIZE_BUDGET_AGGREGATE_PERSISTENT_MISMATCH",
                "aggregate_selected exact persistent bytes do not match the sum of selected rows",
                "analysis/vgp_pilot_size_budget.tsv",
            )
        )
    return blockers, {
        "selected_row_count": len(selected_rows),
        "selected_download_bytes_exact": selected_download,
        "selected_persistent_bytes_exact": selected_persistent,
        "aggregate_selected_download_bytes_exact": aggregate_download,
        "aggregate_selected_persistent_bytes_exact": aggregate_persistent,
        "candidate_budget_row_count": len(candidate_budget),
    }


def parse_root_validation(root_config: Mapping[str, Any], evidence: Mapping[str, Any]) -> tuple[list[dict[str, str]], dict[str, Any]]:
    blockers: list[dict[str, str]] = []
    if evidence.get("root") != root_config.get("root"):
        blockers.append(
            blocker(
                "ROOT_VALIDATION_ROOT_MISMATCH",
                "root validation evidence does not describe the configured VGP data root",
                "analysis/vgp_data_root_validation.json",
            )
        )
    system = evidence.get("system_evidence", {})
    quota_state = system.get("quota_state")
    inode_state = system.get("inode_state")
    if quota_state is None:
        blockers.append(
            blocker(
                "QUOTA_METADATA_MISSING",
                "root validation evidence does not contain quota metadata",
                "analysis/vgp_data_root_validation.json",
            )
        )
        quota_state = {"status": "missing", "available_interfaces": []}
    if inode_state is None:
        blockers.append(
            blocker(
                "INODE_METADATA_MISSING",
                "root validation evidence does not contain inode metadata",
                "analysis/vgp_data_root_validation.json",
            )
        )
        inode_state = {"status": "missing", "output": ""}

    if quota_state.get("status") != "reported":
        blockers.append(
            blocker(
                "QUOTA_UNAVAILABLE",
                quota_state.get("message", "No exact quota interface was reported for the configured VGP root"),
                "analysis/vgp_data_root_validation.json",
            )
        )
    if inode_state.get("status") != "reported":
        blockers.append(
            blocker(
                "INODE_EVIDENCE_UNAVAILABLE",
                "inode availability is not reported for the configured VGP root",
                "analysis/vgp_data_root_validation.json",
            )
        )

    free_bytes = None
    free_inodes = None
    for command in system.get("commands", []):
        text = command.get("stdout", "")
        if "block_size=" in text and "blocks_available=" in text:
            match = re.search(r"block_size=(\d+).*blocks_available=(\d+).*files_free=(\d+)", text)
            if match:
                block_size = int(match.group(1))
                blocks_available = int(match.group(2))
                free_bytes = block_size * blocks_available
                free_inodes = int(match.group(3))
                break
    if free_inodes is None:
        output = inode_state.get("output", "")
        for raw in output.splitlines()[1:]:
            parts = raw.split()
            if len(parts) >= 4:
                try:
                    free_inodes = int(parts[3])
                    break
                except ValueError:
                    continue
    safe_free_bytes = math.floor(free_bytes * 0.75) if free_bytes is not None else None
    safe_free_inodes = math.floor(free_inodes * 0.75) if free_inodes is not None else None
    return blockers, {
        "root": evidence.get("root"),
        "quota_status": quota_state.get("status"),
        "quota_interfaces": quota_state.get("available_interfaces", []),
        "quota_message": quota_state.get("message"),
        "inode_status": inode_state.get("status"),
        "free_bytes": free_bytes,
        "safe_free_bytes": safe_free_bytes,
        "free_inodes": free_inodes,
        "safe_free_inodes": safe_free_inodes,
    }


def parse_small_cap_from_execution_plan(text: str) -> dict[str, float]:
    match = re.search(r"\*\*small cap:\*\*\s*([^|]+)", text)
    if not match:
        raise Tier3ValidationError("failed to locate the small cap clause in the execution plan")
    body = match.group(1)
    patterns = {
        "slots": r"(\d+)\s+slots",
        "cpus_per_element": r"(\d+)\s+CPU/element",
        "core_hours": r"([0-9.]+)\s+core-h",
        "memory_gib_per_job": r"([0-9.]+)\s+GiB/job",
        "wall_hours": r"([0-9.]+)\s+h",
        "scratch_gb": r"([0-9.]+)\s+GB scratch",
        "input_gb": r"([0-9.]+)\s+GB input/download",
        "output_gb": r"\+\s*([0-9.]+)\s+GB output",
        "file_inodes": r"([0-9.kK]+)\s+files",
        "read_gb": r"([0-9.]+)\s+GB read",
        "write_gb": r"/([0-9.]+)\s+GB write",
        "metadata_operations": r"([0-9.kK]+)\s+ops",
        "bandwidth_mib_s": r"([0-9.]+)\s+MiB/s",
        "tier3a_concurrency": r"Tier3A\s+(\d+)",
        "tier3c_concurrency": r"Tier3C\s+(\d+)",
        "transfer_concurrency": r"transfers\s+(\d+)",
    }
    values: dict[str, float] = {}
    for key, pattern in patterns.items():
        found = re.search(pattern, body)
        if not found:
            raise Tier3ValidationError(f"failed to parse {key} from the execution-plan small cap clause")
        text_value = found.group(1)
        values[key] = parse_compact_number(text_value) if any(ch in text_value.lower() for ch in "k") else float(text_value)
    return values


def parse_decision_caps(decision_rows: Sequence[Mapping[str, str]]) -> dict[str, float]:
    by_id = {row["decision_id"]: row for row in decision_rows}
    d011 = by_id["D011"]["resolution_hard_gate"]
    d012 = by_id["D012"]["resolution_hard_gate"]
    d018 = by_id["D018"]["resolution_hard_gate"]
    slots_match = re.search(r"exact\s+([a-z0-9.]+)\s+rows", d011, flags=re.IGNORECASE)
    if not slots_match:
        raise Tier3ValidationError("failed to parse the D011 exact-row cap")
    return {
        "slots": parse_numeric_phrase(slots_match.group(1)),
        "input_gb": float(re.search(r"<=([0-9.]+)\s+GB download/input", d011).group(1)),
        "core_hours": float(re.search(r"<=([0-9.]+)\s+core-h", d012).group(1)),
        "memory_gib_per_job": float(re.search(r"([0-9.]+)\s+GiB/job", d012).group(1)),
        "wall_hours": float(re.search(r"([0-9.]+)\s+h", d012).group(1)),
        "scratch_gb": float(re.search(r"([0-9.]+)\s+GB scratch", d012).group(1)),
        "output_gb": float(re.search(r"([0-9.]+)\s+GB output", d012).group(1)),
        "transfer_concurrency": float(re.search(r"transfers <=(\d+)", d018).group(1)),
        "pause_fraction": float(re.search(r"pause at (\d+)%", d018).group(1)) / 100.0,
        "stop_fraction": float(re.search(r"stop resumption at (\d+)%", d018).group(1)) / 100.0,
        "quota_headroom_fraction": float(re.search(r"keep (\d+)% quota headroom", d018).group(1)) / 100.0,
    }


def collect_resource_budget_rows(rows: Sequence[Mapping[str, str]]) -> dict[str, dict[str, Mapping[str, str]]]:
    by_stage: dict[str, dict[str, Mapping[str, str]]] = {}
    for row in rows:
        stage = row["stage_or_dataset"]
        if stage not in RELEVANT_BUDGET_ROWS:
            continue
        by_stage.setdefault(stage, {})[row["scenario"]] = row
    missing = [stage for stage in RELEVANT_BUDGET_ROWS if stage not in by_stage]
    if missing:
        raise Tier3ValidationError(f"missing relevant resource-budget rows: {', '.join(sorted(missing))}")
    return by_stage


def selected_stage_totals(
    manifest_rows: Sequence[Mapping[str, str]],
    budget_rows: Mapping[str, Mapping[str, Mapping[str, str]]],
) -> dict[str, float]:
    selected = [row for row in manifest_rows if row.get("pilot_selected") == "yes"]
    pair_selected = [row for row in selected if PAIR_MODALITY in split_modalities(row.get("seed_modalities", ""))]
    composition_selected = [row for row in selected if COMPOSITION_MODALITY in split_modalities(row.get("seed_modalities", ""))]

    def row_value(stage: str, field: str) -> float:
        parsed = parse_float(budget_rows[stage]["high"][field])
        if parsed is None:
            raise Tier3ValidationError(f"resource-budget field {stage}.{field} is not a finite number")
        return parsed

    stage_count = len(pair_selected)
    composition_count = len(composition_selected)
    return {
        "selected_species_count": float(len(selected)),
        "selected_pair_species_count": float(stage_count),
        "selected_composition_species_count": float(composition_count),
        "aggregate_core_hours": (
            stage_count
            * (
                row_value("stage_assets_per_species", "core_hours")
                + row_value("exact_preflight_per_species", "core_hours")
                + row_value("sweepga_mapping_per_species", "core_hours")
                + row_value("impg_index_partition_query_per_species", "core_hours")
                + row_value("denominator_summary_per_species", "core_hours")
            )
            + composition_count * row_value("tier3c_composition_per_species", "core_hours")
        ),
        "aggregate_wall_hours": (
            stage_count
            * (
                row_value("stage_assets_per_species", "wall_hours_per_element")
                + row_value("exact_preflight_per_species", "wall_hours_per_element")
                + row_value("sweepga_mapping_per_species", "wall_hours_per_element")
                + row_value("impg_index_partition_query_per_species", "wall_hours_per_element")
                + row_value("denominator_summary_per_species", "wall_hours_per_element")
            )
            + composition_count * row_value("tier3c_composition_per_species", "wall_hours_per_element")
        ),
        "peak_memory_gib_per_job": max(
            [0.0]
            + ([row_value("stage_assets_per_species", "peak_resident_or_requested_memory_gib_per_element")] if stage_count else [])
            + ([row_value("exact_preflight_per_species", "peak_resident_or_requested_memory_gib_per_element")] if stage_count else [])
            + ([row_value("sweepga_mapping_per_species", "peak_resident_or_requested_memory_gib_per_element")] if stage_count else [])
            + ([row_value("impg_index_partition_query_per_species", "peak_resident_or_requested_memory_gib_per_element")] if stage_count else [])
            + ([row_value("denominator_summary_per_species", "peak_resident_or_requested_memory_gib_per_element")] if stage_count else [])
            + ([row_value("tier3c_composition_per_species", "peak_resident_or_requested_memory_gib_per_element")] if composition_count else [])
        ),
        "peak_scratch_gb": max(
            [0.0]
            + ([row_value("stage_assets_per_species", "local_scratch_peak_gb")] if stage_count else [])
            + ([row_value("exact_preflight_per_species", "local_scratch_peak_gb")] if stage_count else [])
            + ([row_value("sweepga_mapping_per_species", "local_scratch_peak_gb")] if stage_count else [])
            + ([row_value("impg_index_partition_query_per_species", "local_scratch_peak_gb")] if stage_count else [])
            + ([row_value("denominator_summary_per_species", "local_scratch_peak_gb")] if stage_count else [])
            + ([row_value("tier3c_composition_per_species", "local_scratch_peak_gb")] if composition_count else [])
        ),
        "persistent_output_gb": (
            stage_count
            * (
                row_value("stage_assets_per_species", "persistent_output_gb")
                + row_value("exact_preflight_per_species", "persistent_output_gb")
                + row_value("sweepga_mapping_per_species", "persistent_output_gb")
                + row_value("impg_index_partition_query_per_species", "persistent_output_gb")
                + row_value("denominator_summary_per_species", "persistent_output_gb")
            )
            + composition_count * row_value("tier3c_composition_per_species", "persistent_output_gb")
        ),
        "file_inodes": (
            stage_count
            * (
                row_value("stage_assets_per_species", "file_inode_count")
                + row_value("exact_preflight_per_species", "file_inode_count")
                + row_value("sweepga_mapping_per_species", "file_inode_count")
                + row_value("impg_index_partition_query_per_species", "file_inode_count")
                + row_value("denominator_summary_per_species", "file_inode_count")
            )
            + composition_count * row_value("tier3c_composition_per_species", "file_inode_count")
        ),
        "moosefs_read_gb": (
            stage_count
            * (
                row_value("stage_assets_per_species", "moosefs_read_gb")
                + row_value("exact_preflight_per_species", "moosefs_read_gb")
                + row_value("sweepga_mapping_per_species", "moosefs_read_gb")
                + row_value("impg_index_partition_query_per_species", "moosefs_read_gb")
                + row_value("denominator_summary_per_species", "moosefs_read_gb")
            )
            + composition_count * row_value("tier3c_composition_per_species", "moosefs_read_gb")
        ),
        "moosefs_write_gb": (
            stage_count
            * (
                row_value("stage_assets_per_species", "moosefs_write_gb")
                + row_value("exact_preflight_per_species", "moosefs_write_gb")
                + row_value("sweepga_mapping_per_species", "moosefs_write_gb")
                + row_value("impg_index_partition_query_per_species", "moosefs_write_gb")
                + row_value("denominator_summary_per_species", "moosefs_write_gb")
            )
            + composition_count * row_value("tier3c_composition_per_species", "moosefs_write_gb")
        ),
        "metadata_operations": (
            stage_count
            * (
                row_value("stage_assets_per_species", "metadata_operations")
                + row_value("exact_preflight_per_species", "metadata_operations")
                + row_value("sweepga_mapping_per_species", "metadata_operations")
                + row_value("impg_index_partition_query_per_species", "metadata_operations")
                + row_value("denominator_summary_per_species", "metadata_operations")
            )
            + composition_count * row_value("tier3c_composition_per_species", "metadata_operations")
        ),
        "peak_bandwidth_mib_s": max(
            [0.0]
            + ([row_value("stage_assets_per_species", "peak_aggregate_bandwidth_mib_s")] if stage_count else [])
            + ([row_value("exact_preflight_per_species", "peak_aggregate_bandwidth_mib_s")] if stage_count else [])
            + ([row_value("sweepga_mapping_per_species", "peak_aggregate_bandwidth_mib_s")] if stage_count else [])
            + ([row_value("impg_index_partition_query_per_species", "peak_aggregate_bandwidth_mib_s")] if stage_count else [])
            + ([row_value("denominator_summary_per_species", "peak_aggregate_bandwidth_mib_s")] if stage_count else [])
            + ([row_value("tier3c_composition_per_species", "peak_aggregate_bandwidth_mib_s")] if composition_count else [])
        ),
    }


def build_cap_vector(
    plan_caps: Mapping[str, float],
    decision_caps: Mapping[str, float],
    budget_rows: Mapping[str, Mapping[str, Mapping[str, str]]],
    manifest_rows: Sequence[Mapping[str, str]],
    size_budget_summary: Mapping[str, Any],
    quota_summary: Mapping[str, Any],
) -> dict[str, Any]:
    stratified_high = budget_rows["stratified_pilot"]["high"]
    stage_totals = selected_stage_totals(manifest_rows, budget_rows)

    def budget_value(field: str) -> float:
        parsed = parse_float(stratified_high[field])
        if parsed is None:
            raise Tier3ValidationError(f"stratified_pilot.high {field} is not a finite number")
        return parsed

    exact_download_gb = (size_budget_summary["selected_download_bytes_exact"] or 0) / 1_000_000_000
    exact_persistent_gb = (size_budget_summary["selected_persistent_bytes_exact"] or 0) / 1_000_000_000
    safe_free_gb = (quota_summary["safe_free_bytes"] or 0) / 1_000_000_000
    safe_free_inodes = float(quota_summary["safe_free_inodes"] or 0)

    raw_dimensions: dict[str, dict[str, Any]] = {
        "selected_species_count": {
            "unit": "count",
            "sources": {
                "selected_manifest_count": stage_totals["selected_species_count"],
            },
        },
        "selected_pair_species_count": {
            "unit": "count",
            "sources": {
                "selected_pair_count": stage_totals["selected_pair_species_count"],
            },
        },
        "transfer_concurrency": {
            "unit": "count",
            "sources": {
                "execution_plan_small_cap": plan_caps["transfer_concurrency"],
                "decision_D018": decision_caps["transfer_concurrency"],
                "resource_budget_stage_assets": parse_float(budget_rows["stage_assets_per_species"]["high"]["concurrency_cap"]),
            },
        },
        "tier3a_concurrency": {
            "unit": "count",
            "sources": {
                "execution_plan_small_cap": plan_caps["tier3a_concurrency"],
                "resource_budget_sweepga": parse_float(budget_rows["sweepga_mapping_per_species"]["high"]["concurrency_cap"]),
            },
        },
        "tier3c_concurrency": {
            "unit": "count",
            "sources": {
                "execution_plan_small_cap": plan_caps["tier3c_concurrency"],
                "resource_budget_tier3c": parse_float(budget_rows["tier3c_composition_per_species"]["high"]["concurrency_cap"]),
            },
        },
        "cpus_per_element": {
            "unit": "count",
            "sources": {
                "execution_plan_small_cap": plan_caps["cpus_per_element"],
                "resource_budget_stratified_pilot": budget_value("cpus_per_element"),
            },
        },
        "aggregate_core_hours": {
            "unit": "core-h",
            "sources": {
                "execution_plan_small_cap": plan_caps["core_hours"],
                "decision_D012": decision_caps["core_hours"],
                "resource_budget_stratified_pilot": budget_value("core_hours"),
                "resource_budget_selected_stage_sum_high": stage_totals["aggregate_core_hours"],
            },
        },
        "peak_memory_gib_per_job": {
            "unit": "GiB",
            "sources": {
                "execution_plan_small_cap": plan_caps["memory_gib_per_job"],
                "decision_D012": decision_caps["memory_gib_per_job"],
                "resource_budget_stratified_pilot": budget_value("peak_resident_or_requested_memory_gib_per_element"),
                "resource_budget_selected_stage_max_high": stage_totals["peak_memory_gib_per_job"],
            },
        },
        "aggregate_wall_hours": {
            "unit": "h",
            "sources": {
                "execution_plan_small_cap": plan_caps["wall_hours"],
                "decision_D012": decision_caps["wall_hours"],
                "resource_budget_stratified_pilot": budget_value("catalog_or_stage_wall_hours"),
                "resource_budget_selected_stage_sum_high": stage_totals["aggregate_wall_hours"],
            },
        },
        "peak_local_scratch_gb": {
            "unit": "GB",
            "sources": {
                "execution_plan_small_cap": plan_caps["scratch_gb"],
                "decision_D012": decision_caps["scratch_gb"],
                "resource_budget_stratified_pilot": budget_value("local_scratch_peak_gb"),
                "resource_budget_selected_stage_max_high": stage_totals["peak_scratch_gb"],
                "filesystem_safe_available": safe_free_gb,
            },
        },
        "persistent_input_gb": {
            "unit": "GB",
            "sources": {
                "execution_plan_small_cap": plan_caps["input_gb"],
                "decision_D011": decision_caps["input_gb"],
                "resource_budget_stratified_pilot": budget_value("persistent_input_gb"),
                "selected_manifest_exact": exact_download_gb,
                "filesystem_safe_available": safe_free_gb,
            },
        },
        "persistent_output_gb": {
            "unit": "GB",
            "sources": {
                "execution_plan_small_cap": plan_caps["output_gb"],
                "decision_D012": decision_caps["output_gb"],
                "resource_budget_stratified_pilot": budget_value("persistent_output_gb"),
                "resource_budget_selected_stage_sum_high": stage_totals["persistent_output_gb"],
                "filesystem_safe_available": safe_free_gb,
            },
        },
        "file_inodes": {
            "unit": "count",
            "sources": {
                "execution_plan_small_cap": plan_caps["file_inodes"],
                "resource_budget_stratified_pilot": budget_value("file_inode_count"),
                "resource_budget_selected_stage_sum_high": stage_totals["file_inodes"],
                "filesystem_safe_available": safe_free_inodes,
            },
        },
        "moosefs_read_gb": {
            "unit": "GB",
            "sources": {
                "execution_plan_small_cap": plan_caps["read_gb"],
                "resource_budget_stratified_pilot": budget_value("moosefs_read_gb"),
                "resource_budget_selected_stage_sum_high": stage_totals["moosefs_read_gb"],
            },
        },
        "moosefs_write_gb": {
            "unit": "GB",
            "sources": {
                "execution_plan_small_cap": plan_caps["write_gb"],
                "resource_budget_stratified_pilot": budget_value("moosefs_write_gb"),
                "resource_budget_selected_stage_sum_high": stage_totals["moosefs_write_gb"],
            },
        },
        "metadata_operations": {
            "unit": "count",
            "sources": {
                "execution_plan_small_cap": plan_caps["metadata_operations"],
                "resource_budget_stratified_pilot": budget_value("metadata_operations"),
                "resource_budget_selected_stage_sum_high": stage_totals["metadata_operations"],
            },
        },
        "peak_bandwidth_mib_s": {
            "unit": "MiB/s",
            "sources": {
                "execution_plan_small_cap": plan_caps["bandwidth_mib_s"],
                "resource_budget_stratified_pilot": budget_value("peak_aggregate_bandwidth_mib_s"),
                "resource_budget_selected_stage_max_high": stage_totals["peak_bandwidth_mib_s"],
            },
        },
        "quota_headroom_fraction": {
            "unit": "fraction",
            "sources": {
                "decision_D018": decision_caps["quota_headroom_fraction"],
            },
        },
        "pause_fraction": {
            "unit": "fraction",
            "sources": {
                "decision_D018": decision_caps["pause_fraction"],
            },
        },
        "stop_fraction": {
            "unit": "fraction",
            "sources": {
                "decision_D018": decision_caps["stop_fraction"],
            },
        },
    }

    winners: dict[str, Any] = {}
    for dimension, spec in raw_dimensions.items():
        source_items = [
            {"source": source_name, "value": float(value)}
            for source_name, value in spec["sources"].items()
            if value is not None and math.isfinite(float(value))
        ]
        if not source_items:
            raise Tier3ValidationError(f"no finite sources were available for cap dimension {dimension}")
        winner = min(source_items, key=lambda item: item["value"])
        winners[dimension] = {
            "unit": spec["unit"],
            "value": winner["value"],
            "winner_source": winner["source"],
            "sources": source_items,
        }

    payload = {
        "dimensions": winners,
        "selected_stage_totals_high": stage_totals,
    }
    payload["sha256"] = sha256_json(payload)
    return payload


def review_markdown(gate: Mapping[str, Any]) -> str:
    lines = [
        "# VGP pilot gate review",
        "",
        f"- Decision: `{gate['decision']['status']}`",
        f"- Decision SHA-256: `{gate['decision_sha256']}`",
        f"- Manifest digest: `{gate['authorization_boundary']['manifest_digest']}`",
        f"- Root contract digest: `{gate['authorization_boundary']['root_contract_digest']}`",
        f"- Cap vector digest: `{gate['authorization_boundary']['cap_vector_digest']}`",
        "",
        "## Reproduced counts",
        "",
        f"- Raw frozen catalog lines: `{gate['reproduction']['source_catalog']['line_count']}`",
        f"- Raw frozen candidate seeds: `{gate['reproduction']['manifest_candidate_count']}`",
        f"- Selected rows in manifest: `{gate['reproduction']['selected_row_count']}`",
        f"- Independently composition-ready rows: `{gate['row_audit']['summary']['composition_ready_count']}`",
        f"- Independently diversity-ready rows: `{gate['row_audit']['summary']['diversity_ready_count']}`",
        "",
        "## Winning caps",
        "",
    ]
    for name, payload in sorted(gate["cap_vector"]["dimensions"].items()):
        lines.append(
            f"- `{name}`: `{payload['value']}` {payload['unit']} "
            f"(winner: `{payload['winner_source']}`)"
        )
    lines.extend(
        [
            "",
            "## Blockers",
            "",
        ]
    )
    if gate["blockers"]:
        for item in gate["blockers"]:
            lines.append(f"- `{item['code']}`: {item['message']} ({item['source']})")
    else:
        lines.append("- none")
    lines.extend(
        [
            "",
            "## Dominant row-level failures",
            "",
        ]
    )
    for code, count in gate["row_audit"]["summary"]["issue_counts"].items():
        lines.append(f"- `{code}`: `{count}` rows")
    lines.append("")
    return "\n".join(lines)


def build_gate(
    *,
    manifest_path: Path = DEFAULT_MANIFEST,
    size_budget_path: Path = DEFAULT_SIZE_BUDGET,
    freeze_provenance_path: Path = DEFAULT_FREEZE_PROVENANCE,
    root_config_path: Path = DEFAULT_ROOT_CONFIG,
    root_validation_path: Path = DEFAULT_ROOT_VALIDATION,
    decisions_path: Path = DEFAULT_DECISIONS,
    execution_plan_path: Path = DEFAULT_EXECUTION_PLAN,
    resource_budget_path: Path = DEFAULT_RESOURCE_BUDGET,
    gate_out: Path = DEFAULT_GATE_JSON,
    review_out: Path = DEFAULT_GATE_REVIEW,
) -> dict[str, Any]:
    manifest_rows = load_tsv(manifest_path)
    size_budget_rows = load_tsv(size_budget_path)
    provenance = load_json(freeze_provenance_path)
    root_config = load_json(root_config_path)
    root_validation = load_json(root_validation_path)
    decision_rows = load_tsv(decisions_path)
    execution_plan_text = execution_plan_path.read_text(encoding="utf-8")
    resource_budget_rows = load_tsv(resource_budget_path)

    blockers: list[dict[str, str]] = []

    source_catalog = read_source_catalog(provenance)
    if source_catalog["sha256"] != provenance["source_catalog"]["sha256"]:
        blockers.append(
            blocker(
                "SOURCE_CATALOG_DIGEST_MISMATCH",
                "the raw frozen VGP catalog digest does not match the recorded provenance digest",
                "analysis/vgp_phase1_freeze_provenance.json",
            )
        )
    if source_catalog["line_count"] != provenance["source_catalog"]["line_count"]:
        blockers.append(
            blocker(
                "SOURCE_CATALOG_LINECOUNT_MISMATCH",
                "the raw frozen VGP catalog line count does not match the recorded provenance line count",
                "analysis/vgp_phase1_freeze_provenance.json",
            )
        )
    if source_catalog["discrepancies"]:
        blockers.append(
            blocker(
                "SOURCE_COUNT_DISCREPANCY_UNRESOLVED",
                "the frozen raw VGP catalog still disagrees with the planning headline counts and no explicit signed discrepancy resolution is bundled here",
                "analysis/vgp_phase1_freeze_provenance.json",
            )
        )

    row_audits = [audit_manifest_row(row) for row in manifest_rows]
    selected_rows = [row for row in manifest_rows if row.get("pilot_selected") == "yes"]
    selected_audits = [audit for audit in row_audits if audit["selected"]]
    composition_ready = sum(1 for audit in row_audits if audit["composition_ready"])
    diversity_ready = sum(1 for audit in row_audits if audit["diversity_ready"])
    issue_counts = Counter(issue["code"] for audit in row_audits for issue in audit["issues"])

    if len(selected_rows) > 6:
        blockers.append(
            blocker(
                "SELECTED_COUNT_EXCEEDS_SIX",
                f"the manifest selects {len(selected_rows)} rows, exceeding the six-species ceiling",
                "analysis/vgp_pilot_manifest.tsv",
            )
        )
    if len(selected_rows) == 0:
        blockers.append(
            blocker(
                "NO_SELECTED_ROWS",
                "the frozen pilot manifest selects zero rows, so no bounded pilot can be authorized",
                "analysis/vgp_pilot_manifest.tsv",
            )
        )
    if composition_ready == 0:
        blockers.append(
            blocker(
                "ZERO_COMPOSITION_ELIGIBLE_ROWS",
                "no manifest row independently satisfies exact H1/native-annotation/denominator requirements",
                "analysis/vgp_pilot_manifest.tsv",
            )
        )
    if diversity_ready == 0:
        blockers.append(
            blocker(
                "ZERO_DIVERSITY_ELIGIBLE_ROWS",
                "no manifest row independently satisfies the paired same-individual diversity requirements",
                "analysis/vgp_pilot_manifest.tsv",
            )
        )
    for audit in selected_audits:
        for issue in audit["issues"]:
            blockers.append(
                blocker(
                    f"SELECTED_{issue['code']}",
                    f"selected row {audit['candidate_id']} failed: {issue['message']}",
                    audit["source"],
                )
            )

    size_budget_blockers, size_budget_summary = manifest_budget_consistency(manifest_rows, size_budget_rows)
    blockers.extend(size_budget_blockers)

    root_blockers, quota_summary = parse_root_validation(root_config, root_validation)
    blockers.extend(root_blockers)

    plan_caps = parse_small_cap_from_execution_plan(execution_plan_text)
    decision_caps = parse_decision_caps(decision_rows)
    relevant_budget_rows = collect_resource_budget_rows(resource_budget_rows)
    cap_vector = build_cap_vector(
        plan_caps,
        decision_caps,
        relevant_budget_rows,
        manifest_rows,
        size_budget_summary,
        quota_summary,
    )

    required_dimensions = (
        "selected_species_count",
        "transfer_concurrency",
        "tier3a_concurrency",
        "tier3c_concurrency",
        "cpus_per_element",
        "aggregate_core_hours",
        "peak_memory_gib_per_job",
        "aggregate_wall_hours",
        "peak_local_scratch_gb",
        "persistent_input_gb",
        "persistent_output_gb",
        "file_inodes",
        "moosefs_read_gb",
        "moosefs_write_gb",
        "metadata_operations",
        "peak_bandwidth_mib_s",
        "quota_headroom_fraction",
        "pause_fraction",
        "stop_fraction",
    )
    for dimension in required_dimensions:
        if dimension not in cap_vector["dimensions"] or not math.isfinite(cap_vector["dimensions"][dimension]["value"]):
            blockers.append(
                blocker(
                    "CAP_VECTOR_DIMENSION_MISSING",
                    f"cap dimension {dimension} is not a finite number",
                    "analysis/vgp_pilot_gate.json",
                )
            )

    projected_persistent_gb = (
        cap_vector["dimensions"]["persistent_input_gb"]["value"]
        + cap_vector["dimensions"]["persistent_output_gb"]["value"]
    )
    safe_free_gb = (quota_summary["safe_free_bytes"] or 0) / 1_000_000_000
    projected_persistent_within_free = projected_persistent_gb <= safe_free_gb if quota_summary["safe_free_bytes"] is not None else False
    if quota_summary["free_inodes"] is not None and cap_vector["dimensions"]["file_inodes"]["value"] > quota_summary["safe_free_inodes"]:
        blockers.append(
            blocker(
                "INODE_HEADROOM_INSUFFICIENT",
                "worst-case persistent inode use would exceed the 25 percent safety margin",
                "analysis/vgp_data_root_validation.json",
            )
        )
    if quota_summary["safe_free_bytes"] is not None and not projected_persistent_within_free:
        blockers.append(
            blocker(
                "FREE_SPACE_HEADROOM_INSUFFICIENT",
                "worst-case persistent byte use would exceed the 25 percent filesystem safety margin",
                "analysis/vgp_data_root_validation.json",
            )
        )

    manifest_digest = sha256_file(manifest_path)
    root_contract_digest = sha256_file(root_config_path)
    decision_status = "GO" if not blockers else "NO_GO"
    authorization_boundary = {
        "manifest_digest": manifest_digest,
        "root_contract_digest": root_contract_digest,
        "cap_vector_digest": cap_vector["sha256"],
        "authorize_commands": {
            "acquire": [
                "python3",
                "-m",
                "analysis.gate_vgp_pilot",
                "authorize",
                "--gate",
                str(gate_out),
                "--manifest",
                str(manifest_path),
                "--root-config",
                str(root_config_path),
                "--action",
                "acquire",
            ],
            "compute": [
                "python3",
                "-m",
                "analysis.gate_vgp_pilot",
                "authorize",
                "--gate",
                str(gate_out),
                "--manifest",
                str(manifest_path),
                "--root-config",
                str(root_config_path),
                "--action",
                "compute",
            ],
        },
    }

    gate: dict[str, Any] = {
        "schema_version": "1.0",
        "generated_at_utc": utc_now(),
        "task_id": "gate-vgp-pilot",
        "inputs": {
            "manifest": {"path": str(manifest_path), "sha256": manifest_digest},
            "size_budget": {"path": str(size_budget_path), "sha256": sha256_file(size_budget_path)},
            "freeze_provenance": {"path": str(freeze_provenance_path), "sha256": sha256_file(freeze_provenance_path)},
            "root_config": {"path": str(root_config_path), "sha256": root_contract_digest},
            "root_validation": {"path": str(root_validation_path), "sha256": sha256_file(root_validation_path)},
            "decisions": {"path": str(decisions_path), "sha256": sha256_file(decisions_path)},
            "execution_plan": {"path": str(execution_plan_path), "sha256": sha256_file(execution_plan_path)},
            "resource_budget": {"path": str(resource_budget_path), "sha256": sha256_file(resource_budget_path)},
        },
        "reproduction": {
            "source_catalog": source_catalog,
            "manifest_candidate_count": len(manifest_rows),
            "selected_row_count": len(selected_rows),
            "selected_candidate_ids": [row["candidate_id"] for row in selected_rows],
            "size_budget_summary": size_budget_summary,
            "freeze_provenance_candidate_summary": provenance.get("candidate_summary", {}),
        },
        "row_audit": {
            "summary": {
                "composition_ready_count": composition_ready,
                "diversity_ready_count": diversity_ready,
                "selected_ready_count": sum(1 for audit in selected_audits if audit["composition_ready"]),
                "issue_counts": dict(sorted(issue_counts.items())),
            },
            "selected_rows": [
                {
                    "candidate_id": audit["candidate_id"],
                    "composition_ready": audit["composition_ready"],
                    "diversity_ready": audit["diversity_ready"],
                    "issues": audit["issues"],
                }
                for audit in selected_audits
            ],
            "blocked_examples": [
                {
                    "candidate_id": audit["candidate_id"],
                    "issues": audit["issues"][:3],
                }
                for audit in row_audits[:8]
            ],
        },
        "quota_evidence": quota_summary,
        "cap_vector": cap_vector,
        "authorization_boundary": authorization_boundary,
        "decision": {
            "status": decision_status,
            "go_requires_exact_match": True,
            "projected_persistent_gb_against_safe_free_gb": {
                "projected": projected_persistent_gb,
                "safe_free_gb": safe_free_gb,
                "within_safe_margin": projected_persistent_within_free,
            },
        },
        "blockers": blockers,
    }
    gate["decision_sha256"] = sha256_json({key: value for key, value in gate.items() if key != "decision_sha256"})
    write_json(gate_out, gate)
    review_out.write_text(review_markdown(gate), encoding="utf-8")
    return gate


def load_gate(path: Path) -> dict[str, Any]:
    gate = load_json(path)
    expected = sha256_json({key: value for key, value in gate.items() if key != "decision_sha256"})
    if gate.get("decision_sha256") != expected:
        raise Tier3ValidationError("gate decision hash does not match the gate payload")
    expected_cap = sha256_json(
        {
            "dimensions": gate["cap_vector"]["dimensions"],
            "selected_stage_totals_high": gate["cap_vector"]["selected_stage_totals_high"],
        }
    )
    if gate["cap_vector"].get("sha256") != expected_cap:
        raise Tier3ValidationError("cap vector hash does not match the gate payload")
    return gate


def authorize_gate_action(gate_path: Path, manifest_path: Path, root_config_path: Path, action: str) -> dict[str, Any]:
    gate = load_gate(gate_path)
    if action not in {"acquire", "compute"}:
        raise Tier3ValidationError(f"unknown gate action {action!r}")
    manifest_digest = sha256_file(manifest_path)
    if manifest_digest != gate["authorization_boundary"]["manifest_digest"]:
        raise Tier3ValidationError(
            f"manifest digest mismatch: expected {gate['authorization_boundary']['manifest_digest']}, observed {manifest_digest}"
        )
    root_digest = sha256_file(root_config_path)
    if root_digest != gate["authorization_boundary"]["root_contract_digest"]:
        raise Tier3ValidationError(
            f"root contract digest mismatch: expected {gate['authorization_boundary']['root_contract_digest']}, observed {root_digest}"
        )
    if gate["decision"]["status"] != "GO":
        raise Tier3ValidationError(
            f"gate decision is {gate['decision']['status']}; {action} is not authorized for this manifest/root/cap vector"
        )
    return {
        "action": action,
        "manifest_digest": manifest_digest,
        "root_contract_digest": root_digest,
        "cap_vector_digest": gate["authorization_boundary"]["cap_vector_digest"],
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=False)

    build = subparsers.add_parser("build", help="build the pilot gate")
    build.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    build.add_argument("--size-budget", type=Path, default=DEFAULT_SIZE_BUDGET)
    build.add_argument("--freeze-provenance", type=Path, default=DEFAULT_FREEZE_PROVENANCE)
    build.add_argument("--root-config", type=Path, default=DEFAULT_ROOT_CONFIG)
    build.add_argument("--root-validation", type=Path, default=DEFAULT_ROOT_VALIDATION)
    build.add_argument("--decisions", type=Path, default=DEFAULT_DECISIONS)
    build.add_argument("--execution-plan", type=Path, default=DEFAULT_EXECUTION_PLAN)
    build.add_argument("--resource-budget", type=Path, default=DEFAULT_RESOURCE_BUDGET)
    build.add_argument("--gate-out", type=Path, default=DEFAULT_GATE_JSON)
    build.add_argument("--review-out", type=Path, default=DEFAULT_GATE_REVIEW)

    authorize = subparsers.add_parser("authorize", help="enforce the gate before acquisition or compute")
    authorize.add_argument("--gate", type=Path, default=DEFAULT_GATE_JSON)
    authorize.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    authorize.add_argument("--root-config", type=Path, default=DEFAULT_ROOT_CONFIG)
    authorize.add_argument("--action", required=True, choices=("acquire", "compute"))

    args = parser.parse_args(argv)
    command = args.command or "build"
    if command == "build":
        build_gate(
            manifest_path=args.manifest,
            size_budget_path=args.size_budget,
            freeze_provenance_path=args.freeze_provenance,
            root_config_path=args.root_config,
            root_validation_path=args.root_validation,
            decisions_path=args.decisions,
            execution_plan_path=args.execution_plan,
            resource_budget_path=args.resource_budget,
            gate_out=args.gate_out,
            review_out=args.review_out,
        )
        return 0
    authorize_gate_action(args.gate, args.manifest, args.root_config, args.action)
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
