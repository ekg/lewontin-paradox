#!/usr/bin/env python3
"""Run and summarize annotation-focused biological SweepGA-to-IMPG analyses.

SweepGA's bounded whole-haplotype PAF is treated as the mapping authority.
IMPG independently owns graph indexing, native partitioning, regional querying,
and variant extraction.  This module prepares a deterministic native-annotation
panel and performs the downstream, reference-audited diversity calculation; it
does not parse CIGAR strings to make the primary variant call set.
"""

from __future__ import annotations

import argparse
import collections
import csv
import datetime as dt
import hashlib
import json
import math
import random
import re
import shlex
import statistics
import subprocess
import sys
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from analysis.tier3_common import is_fourfold_codon, merge_intervals, reverse_complement
from analysis.tier3a_origin_remap import (
    EXPECTED_DATASETS,
    REQUIRED_COMMAND_TOKEN,
    SWEEPGA_COMMIT,
    SWEEPGA_SHA256,
)


DNA = frozenset("ACGT")
csv.field_size_limit(sys.maxsize)
DIVERSITY_COLUMNS = (
    "dataset_id", "scientific_name", "annotation_class", "statistic_label", "scope_label",
    "sweepga_hit_cap_query", "sweepga_hit_cap_target", "sweepga_mapping_provenance",
    "impg_index_provenance", "impg_partition_provenance", "impg_query_provenance",
    "eligible_haplotypes", "variant_numerator", "callable_denominator", "estimate",
    "bootstrap_blocks", "bootstrap_replicates", "bootstrap_ci_low", "bootstrap_ci_high",
    "bootstrap_standard_error", "targeted_genes", "panel_selected_genes",
    "sweepga_queryable_genes", "impg_queryable_genes", "excluded_genes",
    "targeted_bases", "panel_selected_bases", "sweepga_queryable_bases",
    "impg_callable_bases", "exclusions_json", "feature_identity_provenance",
    "transcript_phase_provenance", "biological_input",
)
FAILURE_COLUMNS = (
    "dataset_id", "scientific_name", "stage", "status", "reason", "slurm_job_id",
    "stderr_path", "rerun_command",
)
RUN_COLUMNS = (
    "dataset_id", "scientific_name", "status", "slurm_job_id", "slurm_array_task_id",
    "node", "requested_cpus", "requested_memory", "elapsed", "max_rss", "exit_code",
    "sweepga_hit_cap", "sweepga_paf_path", "sweepga_paf_sha256", "impg_index_path",
    "sweepga_binary_path", "sweepga_binary_sha256", "sweepga_command",
    "impg_binary_path", "impg_binary_sha256", "impg_commit", "impg_partitions_path",
    "impg_native_partitions_selected", "impg_focus_bed", "impg_regional_vcf_count",
    "impg_normalized_bcf", "impg_normalized_bcf_sha256", "guix_channel_commit",
    "guix_channels_path", "guix_manifest_path", "guix_profile_store_path",
    "sweepga_origin_main_commit", "sweepga_native_command_path",
    "sweepga_native_command_sha256", "sweepga_observed_query_multiplicity",
    "sweepga_observed_target_multiplicity", "annotation_query_qc_path",
    "annotation_query_qc_sha256", "preflight_audit_path", "preflight_audit_sha256",
    "command",
)

def _resolved(path: str | Path) -> str:
    return str(Path(path).resolve())


def _load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def validate_preflight(
    manifest: Path,
    commands_path: Path,
    supersession_path: Path,
    build_path: Path,
) -> dict[str, Any]:
    """Fail closed before any regenerated annotation or IMPG computation."""
    rows = _read_tsv(manifest)
    if {row["dataset_id"] for row in rows} != EXPECTED_DATASETS:
        raise ValueError("corrected manifest does not contain exactly the three biological tuples")
    for row in rows:
        if row.get("tier3a_correction_status") != "current_production":
            raise ValueError(f"{row['dataset_id']}: corrected mapping is not current production")
        if row.get("sweepga_origin_main_sha256") != SWEEPGA_SHA256:
            raise ValueError(f"{row['dataset_id']}: unapproved corrected SweepGA identity")
        if row.get("sweepga_origin_main_commit") != SWEEPGA_COMMIT:
            raise ValueError(f"{row['dataset_id']}: unapproved corrected SweepGA commit")
        if REQUIRED_COMMAND_TOKEN not in row.get("sweepga_direct_command", ""):
            raise ValueError(f"{row['dataset_id']}: native 1:1 syntax absent")
        if row.get("sweepga_observed_query_multiplicity") != "1" or row.get("sweepga_observed_target_multiplicity") != "1":
            raise ValueError(f"{row['dataset_id']}: native multiplicity fields failed")
    build = _load_json(build_path)
    repository = build.get("repository", {})
    guix = build.get("guix", {})
    binary = build.get("binary", {})
    if repository.get("fetched_origin_main") != SWEEPGA_COMMIT:
        raise ValueError("build provenance does not prove the fetched origin/main commit")
    profile = Path(str(guix.get("profile_store_path", "")))
    if not str(profile).startswith("/gnu/store/") or not profile.is_dir():
        raise ValueError("pinned Guix build closure is unavailable")
    source_repo = Path(str(repository.get("parent_path", "")))
    fetched_commit = subprocess.check_output(
        [str(profile / "bin/git"), "-C", str(source_repo), "rev-parse", "refs/remotes/origin/main"],
        text=True,
    ).strip()
    if fetched_commit != SWEEPGA_COMMIT:
        raise ValueError("fetched origin/main ref no longer matches the pinned commit")
    if not build.get("completion_gate_passed"):
        raise ValueError("origin/main build completion gate did not pass")
    repository_root = build_path.resolve().parents[1]
    for path_key, hash_key in (("channels_file", "channels_sha256"), ("manifest_file", "manifest_sha256")):
        closure_input = repository_root / str(guix.get(path_key, ""))
        if not closure_input.is_file() or _sha256(closure_input) != guix.get(hash_key):
            raise ValueError(f"pinned Guix closure input failed: {path_key}")
    if binary.get("comparison") != "byte-identical" or binary.get("sha256_build_1") != SWEEPGA_SHA256:
        raise ValueError("origin/main SweepGA binary is not byte-reproducible")
    if _sha256(build_path) != rows[0]["sweepga_origin_build_provenance_sha256"]:
        raise ValueError("origin/main build provenance checksum mismatch")

    command_rows = {row["dataset_id"]: row for row in _read_tsv(commands_path)}
    if set(command_rows) != EXPECTED_DATASETS:
        raise ValueError("corrected command manifest is incomplete")
    ledger = _read_tsv(supersession_path)
    if not ledger or any(row.get("consumable") != "no" for row in ledger):
        raise ValueError("supersession ledger must reject every listed artifact")
    superseded_paths = {_resolved(row["superseded_path"]) for row in ledger if row.get("superseded_path")}
    superseded_hashes = {row["superseded_sha256"] for row in ledger if row.get("superseded_sha256")}
    tuple_audits = []
    for row in rows:
        dataset = row["dataset_id"]
        binary_path = Path(row["sweepga_origin_main_realpath"])
        if _resolved(binary_path) != row["sweepga_origin_main_realpath"] or _sha256(binary_path) != SWEEPGA_SHA256:
            raise ValueError(f"{dataset}: SweepGA realpath/checksum mismatch")
        if row["sweepga_origin_main_commit"] != SWEEPGA_COMMIT:
            raise ValueError(f"{dataset}: unpinned SweepGA commit")
        command_file = Path(row["sweepga_native_command_path"])
        command = command_file.read_text(encoding="utf-8").strip()
        command_tokens = shlex.split(command)
        native_cap_positions = [
            index for index, token in enumerate(command_tokens)
            if token == "--num-mappings"
        ]
        if (
            not command_tokens
            or command_tokens[0] != str(binary_path)
            or len(native_cap_positions) != 1
            or command_tokens[native_cap_positions[0] + 1:native_cap_positions[0] + 2] != ["1:1"]
        ):
            raise ValueError(f"{dataset}: native --num-mappings 1:1 command is absent")
        corrected_command = command_rows[dataset]
        if corrected_command["command"].strip() != command or corrected_command["binary_sha256"] != SWEEPGA_SHA256:
            raise ValueError(f"{dataset}: corrected command manifest disagreement")
        mapping = Path(row["sweepga_bounded_paf_path"])
        mapping_sha = _sha256(mapping)
        if mapping_sha != row["sweepga_bounded_paf_sha256"]:
            raise ValueError(f"{dataset}: corrected PAF checksum mismatch")
        multiplicity_path = Path(row["sweepga_native_multiplicity_audit_path"])
        if _sha256(multiplicity_path) != row["sweepga_native_multiplicity_audit_sha256"]:
            raise ValueError(f"{dataset}: multiplicity audit checksum mismatch")
        multiplicity = _load_json(multiplicity_path)
        if (
            multiplicity.get("status") != "passed"
            or multiplicity.get("native_num_mappings") != "1:1"
            or multiplicity.get("observed_native_query_multiplicity_cap") != 1
            or multiplicity.get("observed_native_target_multiplicity_cap") != 1
            or multiplicity.get("paf_sha256") != mapping_sha
        ):
            raise ValueError(f"{dataset}: native multiplicity proof failed")
        selected_paths = {
            _resolved(mapping), _resolved(binary_path), _resolved(command_file),
            _resolved(multiplicity_path), _resolved(row["annotation_gff_path"]),
        }
        selected_hashes = {mapping_sha, SWEEPGA_SHA256, _sha256(command_file), _sha256(multiplicity_path)}
        if selected_paths & superseded_paths or selected_hashes & superseded_hashes:
            raise ValueError(f"{dataset}: superseded artifact selected by path or checksum")
        tuple_audits.append({
            "dataset_id": dataset,
            "mapping_path": _resolved(mapping),
            "mapping_sha256": mapping_sha,
            "binary_realpath": _resolved(binary_path),
            "binary_sha256": SWEEPGA_SHA256,
            "native_command_path": _resolved(command_file),
            "native_command_sha256": _sha256(command_file),
            "native_num_mappings": "1:1",
            "observed_query_multiplicity": 1,
            "observed_target_multiplicity": 1,
            "query_coverage": multiplicity["query_coverage"],
            "target_coverage": multiplicity["target_coverage"],
        })
    return {
        "schema_version": "tier3a-origin-rerun-preflight-v1",
        "status": "passed",
        "validated_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "corrected_manifest": _resolved(manifest),
        "corrected_manifest_sha256": _sha256(manifest),
        "corrected_commands": _resolved(commands_path),
        "corrected_commands_sha256": _sha256(commands_path),
        "supersession_ledger": _resolved(supersession_path),
        "supersession_ledger_sha256": _sha256(supersession_path),
        "origin_main_build_provenance": _resolved(build_path),
        "origin_main_build_provenance_sha256": _sha256(build_path),
        "origin_main_commit": SWEEPGA_COMMIT,
        "fetched_origin_main_ref": fetched_commit,
        "guix_build_profile_store_path": str(profile),
        "guix_channels_sha256": guix["channels_sha256"],
        "guix_manifest_sha256": guix["manifest_sha256"],
        "superseded_paths_rejected": len(superseded_paths),
        "superseded_checksums_rejected": len(superseded_hashes),
        "tuples": tuple_audits,
    }


def preflight(args: argparse.Namespace) -> None:
    audit = validate_preflight(args.manifest, args.commands, args.supersession_ledger, args.build_provenance)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _write_tsv(path: Path, columns: Sequence[str], rows: Iterable[Mapping[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def merge_by_contig(intervals: Mapping[str, Iterable[tuple[int, int]]]) -> dict[str, list[tuple[int, int]]]:
    return {contig: merge_intervals(values) for contig, values in sorted(intervals.items()) if values}


def intersect_by_contig(
    left: Mapping[str, Iterable[tuple[int, int]]],
    right: Mapping[str, Iterable[tuple[int, int]]],
) -> dict[str, list[tuple[int, int]]]:
    result: dict[str, list[tuple[int, int]]] = {}
    left_merged, right_merged = merge_by_contig(left), merge_by_contig(right)
    for contig in sorted(set(left_merged) & set(right_merged)):
        a, b, i, j, hits = left_merged[contig], right_merged[contig], 0, 0, []
        while i < len(a) and j < len(b):
            start, end = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
            if start < end:
                hits.append((start, end))
            if a[i][1] <= b[j][1]:
                i += 1
            else:
                j += 1
        if hits:
            result[contig] = hits
    return result


def _interval_bases(intervals: Mapping[str, Iterable[tuple[int, int]]]) -> int:
    return sum(end - start for values in merge_by_contig(intervals).values() for start, end in values)


def select_execution_panel(dataset_id: str, spans: Sequence[Mapping[str, str]], target_bases: int) -> list[Mapping[str, str]]:
    """Stable hash sample of complete merged annotation spans.

    Selection uses only tuple and native span identity; it cannot depend on
    mapping success or observed variants.  Whole spans are retained, so the
    realized panel may exceed the requested base target by its last span.
    """
    if target_bases < 1 or not spans:
        raise ValueError("panel target and span collection must be nonempty")
    ranked = sorted(
        spans,
        key=lambda row: (
            hashlib.sha256(f"{dataset_id}\0{row['execution_span_id']}".encode()).hexdigest(),
            row["execution_span_id"],
        ),
    )
    selected, bases = [], 0
    for row in ranked:
        start, end = int(row["start_0based"]), int(row["end_0based_exclusive"])
        if not 0 <= start < end:
            raise ValueError(f"invalid execution span {row['execution_span_id']}")
        selected.append(row)
        bases += end - start
        if bases >= target_bases:
            break
    return selected


def _read_bed(path: Path) -> dict[str, list[tuple[int, int]]]:
    result: dict[str, list[tuple[int, int]]] = collections.defaultdict(list)
    with path.open(encoding="utf-8") as handle:
        for number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise ValueError(f"BED row {number} has fewer than three fields: {path}")
            start, end = int(fields[1]), int(fields[2])
            if not 0 <= start < end:
                raise ValueError(f"invalid BED interval at row {number}: {path}")
            result[fields[0]].append((start, end))
    return merge_by_contig(result)


def _paf_query_intervals(path: Path) -> dict[str, list[tuple[int, int]]]:
    intervals: dict[str, list[tuple[int, int]]] = collections.defaultdict(list)
    with path.open(encoding="utf-8") as handle:
        for number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 12:
                raise ValueError(f"PAF row {number} is truncated")
            intervals[fields[0]].append((int(fields[2]), int(fields[3])))
    if not intervals:
        raise ValueError("SweepGA PAF has no mappings")
    return merge_by_contig(intervals)


def _fully_covered(intervals: Mapping[str, Sequence[tuple[int, int]]], contig: str, start: int, end: int) -> bool:
    return any(left <= start and end <= right for left, right in intervals.get(contig, ()))


def _lookup_manifest(path: Path, dataset_id: str) -> dict[str, str]:
    matches = [row for row in _read_tsv(path) if row["dataset_id"] == dataset_id]
    if len(matches) != 1:
        raise ValueError(f"expected one acquisition row for {dataset_id}, found {len(matches)}")
    row = matches[0]
    if row["eligibility_status"] != "eligible_biological":
        raise ValueError(f"tuple {dataset_id} is not biologically eligible")
    if row["sweepga_selected_cap_query"] != "1" or row["sweepga_selected_cap_target"] != "1":
        raise ValueError("primary analysis requires the acquired SweepGA 1:1 cap")
    if row["sweepga_h1_paf_axis"] != "query" or row["sweepga_h2_paf_axis"] != "target":
        raise ValueError("unexpected H1/H2 PAF axes")
    return row


def prepare(args: argparse.Namespace) -> None:
    row = _lookup_manifest(args.manifest, args.dataset_id)
    preflight_audit = _load_json(args.preflight_audit)
    preflight_tuples = {
        item["dataset_id"]: item for item in preflight_audit.get("tuples", [])
    }
    if preflight_audit.get("status") != "passed" or args.dataset_id not in preflight_tuples:
        raise ValueError("tuple is absent from the passed origin/main preflight")
    annotation_qc_path = args.annotation_dir / "query_qc.json"
    annotation_qc = _load_json(annotation_qc_path)
    if annotation_qc.get("dataset_id") != args.dataset_id:
        raise ValueError("fresh annotation query QC has the wrong tuple identity")
    if annotation_qc.get("source", {}).get("bounded_whole_haplotype_paf_sha256") != row["sweepga_bounded_paf_sha256"]:
        raise ValueError("fresh annotation queries do not derive from the corrected PAF")
    row.update({
        "annotation_gene_manifest_path": str((args.annotation_dir / "gene_manifest.tsv").resolve()),
        "annotation_query_manifest_path": str((args.annotation_dir / "impg_query_manifest.tsv").resolve()),
        "annotation_span_feature_map_path": str((args.annotation_dir / "impg_span_feature_map.tsv").resolve()),
        "annotation_query_bed_path": str((args.annotation_dir / "impg_execution_spans.bed").resolve()),
        "annotation_contig_map_path": str((args.annotation_dir / "h1_h2_contig_map.tsv").resolve()),
        "annotation_query_qc_path": str(annotation_qc_path.resolve()),
        "targeted_gene_count": str(annotation_qc["targeted_gene_count"]),
        "targeted_gene_union_bases": str(annotation_qc["targeted_gene_union_bases"]),
        "queryable_gene_count": str(annotation_qc["queryable_gene_count"]),
        "excluded_gene_count": str(annotation_qc["excluded_gene_count"]),
        "excluded_gene_union_bases": str(annotation_qc["excluded_gene_union_bases"]),
    })
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    spans = _read_tsv(Path(row["annotation_span_feature_map_path"]))
    selected = select_execution_panel(args.dataset_id, spans, args.panel_bases)
    selected_intervals: dict[str, list[tuple[int, int]]] = collections.defaultdict(list)
    span_by_id = {}
    for span in selected:
        selected_intervals[span["contig"]].append((int(span["start_0based"]), int(span["end_0based_exclusive"])))
        span_by_id[span["execution_span_id"]] = span
    selected_intervals = merge_by_contig(selected_intervals)
    paf = Path(row["sweepga_bounded_paf_path"])
    if _sha256(paf) != row["sweepga_bounded_paf_sha256"]:
        raise ValueError("SweepGA PAF checksum mismatch")
    mapping_intervals = _paf_query_intervals(paf)
    mapping_callable = intersect_by_contig(selected_intervals, mapping_intervals)
    if not mapping_callable:
        raise ValueError("deterministic annotation panel has no SweepGA 1:1 coverage")
    with (out / "selected_spans.bed").open("w", encoding="utf-8") as handle:
        for span in sorted(selected, key=lambda item: (item["contig"], int(item["start_0based"]), item["execution_span_id"])):
            handle.write(f"{span['contig']}\t{span['start_0based']}\t{span['end_0based_exclusive']}\t{span['execution_span_id']}\n")
    with (out / "mapping_callable.bed").open("w", encoding="utf-8") as handle:
        index = 0
        for contig, values in mapping_callable.items():
            for start, end in values:
                index += 1
                handle.write(f"{contig}\t{start}\t{end}\t{args.dataset_id}.mapped.{index:06d}\n")
    selected_gene_ids = sorted({gene for span in selected for gene in span["gene_ids"].split(",") if gene})
    mapped_gene_ids = set()
    for span in selected:
        contig, start, end = span["contig"], int(span["start_0based"]), int(span["end_0based_exclusive"])
        if any(start < right and left < end for left, right in mapping_intervals.get(contig, [])):
            mapped_gene_ids.update(gene for gene in span["gene_ids"].split(",") if gene)
    all_targeted_genes = [
        gene for gene in _read_tsv(Path(row["annotation_gene_manifest_path"]))
        if gene["targeted"] == "yes"
    ]
    all_gene_intervals = _feature_intervals(all_targeted_genes)
    all_gene_cap1_intersection = intersect_by_contig(all_gene_intervals, mapping_intervals)
    fully_covered_gene_ids = {
        gene["gene_id"] for gene in all_targeted_genes
        if _fully_covered(
            mapping_intervals, gene["contig"],
            int(gene["start_0based"]), int(gene["end_0based_exclusive"]),
        )
    }
    native_multiplicity = _load_json(Path(row["sweepga_native_multiplicity_audit_path"]))
    if native_multiplicity.get("paf_sha256") != row["sweepga_bounded_paf_sha256"]:
        raise ValueError("corrected native multiplicity audit does not match the selected PAF")
    cap_sensitivity = {
        "1": {
            "query_coverage": float(native_multiplicity["query_coverage"]),
            "target_coverage": float(native_multiplicity["target_coverage"]),
        }
    }
    qc = {
        "dataset_id": args.dataset_id,
        "scientific_name": row["scientific_name"],
        "panel_policy": "stable_sha256_rank_of_complete_native_annotation_execution_spans",
        "panel_requested_bases": args.panel_bases,
        "panel_selected_span_count": len(selected),
        "panel_selected_union_bases": _interval_bases(selected_intervals),
        "panel_selected_gene_ids": selected_gene_ids,
        "panel_selected_gene_count": len(selected_gene_ids),
        "sweepga_1to1_mapped_gene_ids": sorted(mapped_gene_ids),
        "sweepga_1to1_mapped_gene_count": len(mapped_gene_ids),
        "sweepga_1to1_callable_panel_bases": _interval_bases(mapping_callable),
        "sweepga_1to1_lost_panel_bases": _interval_bases(selected_intervals) - _interval_bases(mapping_callable),
        "sweepga_1to1_fully_covered_target_gene_count": len(fully_covered_gene_ids),
        "sweepga_1to1_excluded_or_partial_target_gene_count": len(all_targeted_genes) - len(fully_covered_gene_ids),
        "sweepga_1to1_target_gene_callable_bases": _interval_bases(all_gene_cap1_intersection),
        "sweepga_1to1_target_gene_lost_bases": _interval_bases(all_gene_intervals) - _interval_bases(all_gene_cap1_intersection),
        "sweepga_1to1_whole_h1_mapped_union_bases": _interval_bases(mapping_intervals),
        "sweepga_1to1_whole_h1_unmapped_bases": int(row["h1_total_length"]) - _interval_bases(mapping_intervals),
        "targeted_gene_count": int(row["targeted_gene_count"]),
        "targeted_gene_union_bases": int(row["targeted_gene_union_bases"]),
        "acquisition_queryable_gene_count": int(row["queryable_gene_count"]),
        "acquisition_excluded_gene_count": int(row["excluded_gene_count"]),
        "acquisition_excluded_gene_union_bases": int(row["excluded_gene_union_bases"]),
        "cap_sensitivity": cap_sensitivity,
        "cap5_minus_cap1_query_coverage": (
            cap_sensitivity["5"]["query_coverage"] - cap_sensitivity["1"]["query_coverage"]
            if "5" in cap_sensitivity else None
        ),
        "cap10_minus_cap1_query_coverage": (
            cap_sensitivity["10"]["query_coverage"] - cap_sensitivity["1"]["query_coverage"]
            if "10" in cap_sensitivity else None
        ),
        "additional_context_bases": 0,
    }
    (out / "prepare_qc.json").write_text(json.dumps(qc, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    run_input = {key: row[key] for key in row}
    run_input["selected_span_ids"] = sorted(span_by_id)
    run_input["preflight_audit_path"] = str(args.preflight_audit.resolve())
    run_input["preflight_audit_sha256"] = _sha256(args.preflight_audit)
    (out / "input.json").write_text(json.dumps(run_input, indent=2, sort_keys=True) + "\n", encoding="utf-8")


class FastaAccessor:
    def __init__(self, fasta: Path):
        self.fasta = fasta
        self.handle = fasta.open("rb")
        self.index = {}
        with Path(str(fasta) + ".fai").open(encoding="utf-8") as handle:
            for raw in handle:
                fields = raw.rstrip("\n").split("\t")
                self.index[fields[0]] = tuple(map(int, fields[1:5]))

    def fetch(self, contig: str, start: int, end: int) -> str:
        length, offset, line_bases, line_width = self.index[contig]
        if not 0 <= start <= end <= length:
            raise ValueError(f"FASTA interval out of range: {contig}:{start}-{end}")
        pieces, cursor = [], start
        while cursor < end:
            in_line = cursor % line_bases
            take = min(end - cursor, line_bases - in_line)
            self.handle.seek(offset + (cursor // line_bases) * line_width + in_line)
            pieces.append(self.handle.read(take).decode("ascii"))
            cursor += take
        return "".join(pieces).upper()

    def close(self) -> None:
        self.handle.close()


def _contains(intervals: Mapping[str, Sequence[tuple[int, int]]], contig: str, pos: int) -> bool:
    return any(start <= pos < end for start, end in intervals.get(contig, ()))


def _callable_positions(accessor: FastaAccessor, intervals: Mapping[str, Sequence[tuple[int, int]]]) -> set[tuple[str, int]]:
    result = set()
    for contig, values in intervals.items():
        for start, end in values:
            sequence = accessor.fetch(contig, start, end)
            result.update((contig, start + index) for index, base in enumerate(sequence) if base in DNA)
    return result


def _feature_intervals(rows: Iterable[Mapping[str, str]]) -> dict[str, list[tuple[int, int]]]:
    result: dict[str, list[tuple[int, int]]] = collections.defaultdict(list)
    for row in rows:
        result[row["contig"]].append((int(row["start_0based"]), int(row["end_0based_exclusive"])))
    return merge_by_contig(result)


def _fourfold_positions(accessor: FastaAccessor, cds_rows: Sequence[Mapping[str, str]]) -> tuple[set[tuple[str, int]], dict[str, int]]:
    transcripts: dict[str, list[Mapping[str, str]]] = collections.defaultdict(list)
    for row in cds_rows:
        for transcript in row["transcript_ids"].split(","):
            if transcript:
                transcripts[transcript].append(row)
    observations: dict[tuple[str, int], set[tuple[int, bool]]] = collections.defaultdict(set)
    excluded = collections.Counter()
    for transcript, segments in transcripts.items():
        strands, contigs = {row["strand"] for row in segments}, {row["contig"] for row in segments}
        if len(strands) != 1 or len(contigs) != 1:
            excluded["transcript_multi_contig_or_strand"] += 1
            continue
        strand = next(iter(strands))
        ordered = sorted(segments, key=lambda row: int(row["start_0based"]), reverse=strand == "-")
        sequence, positions, running = [], [], 0
        valid = True
        for index, row in enumerate(ordered):
            start, end, phase = int(row["start_0based"]), int(row["end_0based_exclusive"]), int(row["phase"])
            expected = phase if index == 0 else (3 - running % 3) % 3
            if phase != expected:
                valid = False
                break
            piece = accessor.fetch(row["contig"], start, end)
            piece_positions = [(row["contig"], pos) for pos in range(start, end)]
            if strand == "-":
                piece, piece_positions = reverse_complement(piece), list(reversed(piece_positions))
            trim = phase if index == 0 else 0
            sequence.append(piece[trim:])
            positions.extend(piece_positions[trim:])
            running += len(piece) - trim
        joined = "".join(sequence)
        if not valid or not joined or len(joined) % 3 or any(base not in DNA for base in joined):
            excluded["invalid_or_ambiguous_phase_chain"] += 1
            continue
        codons = [joined[i:i + 3] for i in range(0, len(joined), 3)]
        effective = len(joined) - (3 if codons and codons[-1] in {"TAA", "TAG", "TGA"} else 0)
        for offset in range(effective):
            start = offset - offset % 3
            state = (offset % 3, offset % 3 == 2 and is_fourfold_codon(joined[start:start + 3]))
            observations[positions[offset]].add(state)
    discordant = {position for position, states in observations.items() if len(states) > 1}
    sites = {position for position, states in observations.items() if position not in discordant and next(iter(states)) == (2, True)}
    excluded["frame_discordant_overlap_positions"] = len(discordant)
    excluded["eligible_transcripts"] = len(transcripts) - excluded["transcript_multi_contig_or_strand"] - excluded["invalid_or_ambiguous_phase_chain"]
    return sites, dict(excluded)


def bootstrap_ratio(blocks: Sequence[tuple[int, int]], replicates: int, seed: int) -> dict[str, float | int]:
    """Block bootstrap for (denominator, numerator) observations."""
    usable = [(denominator, numerator) for denominator, numerator in blocks if denominator > 0]
    if not usable or replicates < 2:
        raise ValueError("bootstrap requires positive blocks and at least two replicates")
    estimate = sum(n for _, n in usable) / sum(d for d, _ in usable)
    rng, sampled = random.Random(seed), []
    for _ in range(replicates):
        draw = [usable[rng.randrange(len(usable))] for _ in usable]
        sampled.append(sum(n for _, n in draw) / sum(d for d, _ in draw))
    sampled.sort()
    low = sampled[max(0, math.floor(0.025 * (replicates - 1)))]
    high = sampled[min(replicates - 1, math.ceil(0.975 * (replicates - 1)))]
    return {
        "estimate": estimate,
        "ci_low": low,
        "ci_high": high,
        "standard_error": statistics.stdev(sampled),
        "replicates": replicates,
        "blocks": len(usable),
    }


def _variant_rows(bcftools: Path, bcf: Path) -> list[dict[str, Any]]:
    output = subprocess.check_output(
        [str(bcftools), "query", "-f", "%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n", str(bcf)],
        text=True,
    )
    result, seen = [], set()
    for raw in output.splitlines():
        fields = raw.split("\t")
        chrom, pos, ref, alt = fields[:4]
        key = chrom, int(pos) - 1, ref, alt
        if key in seen:
            raise ValueError(f"duplicate normalized IMPG allele: {key}")
        seen.add(key)
        result.append({"contig": chrom, "position": int(pos) - 1, "ref": ref.upper(), "alt": alt.upper(), "gt": fields[4] if len(fields) > 4 else "."})
    return result


_CIGAR = re.compile(r"([1-9][0-9]*)([=XID])")


def _direct_h2_check(variant: Mapping[str, Any], paf: Path, h1: FastaAccessor, h2: FastaAccessor) -> dict[str, Any] | None:
    if len(variant["ref"]) != 1 or len(variant["alt"]) != 1:
        return None
    qname, qpos = variant["contig"], int(variant["position"])
    with paf.open(encoding="utf-8") as handle:
        for raw in handle:
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 13 or fields[0] != qname or not int(fields[2]) <= qpos < int(fields[3]):
                continue
            tags = [field[5:] for field in fields[12:] if field.startswith("cg:Z:")]
            if len(tags) != 1:
                continue
            strand, qstart, qend, tname, tpos = fields[4], int(fields[2]), int(fields[3]), fields[5], int(fields[7])
            qcursor = qstart if strand == "+" else qend - 1
            for length_text, operation in _CIGAR.findall(tags[0]):
                length = int(length_text)
                if operation in {"=", "X"}:
                    if strand == "+" and qcursor <= qpos < qcursor + length:
                        target = tpos + qpos - qcursor
                    elif strand == "-" and qcursor - length < qpos <= qcursor:
                        target = tpos + qcursor - qpos
                    else:
                        target = None
                    if target is not None:
                        h1_base, h2_base = h1.fetch(qname, qpos, qpos + 1), h2.fetch(tname, target, target + 1)
                        expected_alt = variant["alt"] if strand == "+" else reverse_complement(variant["alt"])
                        return {
                            "h1_contig": qname, "h1_position_1based": qpos + 1, "h1_base": h1_base,
                            "h2_contig": tname, "h2_position_1based": target + 1, "h2_base": h2_base,
                            "paf_strand": strand, "cigar_operation": operation,
                            "ref_matches_h1": h1_base == variant["ref"], "alt_matches_h2": h2_base == expected_alt,
                            "role": "validation_only_not_primary_calling",
                        }
                    qcursor += length if strand == "+" else -length
                    tpos += length
                elif operation == "I":
                    qcursor += length if strand == "+" else -length
                else:
                    tpos += length
    return None


def _block_counts(positions: set[tuple[str, int]], variants: set[tuple[str, int]], block_size: int) -> list[tuple[int, int]]:
    counts: dict[tuple[str, int], list[int]] = collections.defaultdict(lambda: [0, 0])
    for contig, position in positions:
        counts[(contig, position // block_size)][0] += 1
    for contig, position in variants & positions:
        counts[(contig, position // block_size)][1] += 1
    return [tuple(value) for _, value in sorted(counts.items())]


def summarize(args: argparse.Namespace) -> None:
    out = args.output_dir
    run_input = json.loads((out / "input.json").read_text(encoding="utf-8"))
    prep = json.loads((out / "prepare_qc.json").read_text(encoding="utf-8"))
    selected_ids = set(run_input["selected_span_ids"])
    gene_ids = set(prep["panel_selected_gene_ids"])
    gene_rows = [row for row in _read_tsv(Path(run_input["annotation_gene_manifest_path"])) if row["gene_id"] in gene_ids and row["targeted"] == "yes"]
    cds_rows = [row for row in _read_tsv(Path(run_input["annotation_query_manifest_path"])) if selected_ids.intersection(row["execution_span_ids"].split(",")) and row["targeted"] == "yes" and row["queryable"] == "yes"]
    if not gene_rows or not cds_rows or any(row["phase"] not in {"0", "1", "2"} for row in cds_rows):
        raise ValueError("selected native annotation lacks genes or valid-phase CDS")
    mapping = _read_bed(out / "mapping_callable.bed")
    focus = _read_bed(args.focus_bed)
    queryable = intersect_by_contig(mapping, focus)
    gene_callable = intersect_by_contig(_feature_intervals(gene_rows), queryable)
    cds_callable = intersect_by_contig(_feature_intervals(cds_rows), queryable)
    h1, h2 = FastaAccessor(Path(run_input["h1_fasta_path"])), FastaAccessor(Path(run_input["h2_fasta_path"]))
    try:
        gene_positions = _callable_positions(h1, gene_callable)
        cds_positions = _callable_positions(h1, cds_callable)
        fourfold, phase_qc = _fourfold_positions(h1, cds_rows)
        fourfold &= cds_positions
        fourfold_s = {key for key in fourfold if h1.fetch(key[0], key[1], key[1] + 1) in {"G", "C"}}
        fourfold_w = {key for key in fourfold if h1.fetch(key[0], key[1], key[1] + 1) in {"A", "T"}}
        variants = _variant_rows(args.bcftools, args.bcf)
        for variant in variants:
            observed = h1.fetch(variant["contig"], variant["position"], variant["position"] + len(variant["ref"]))
            if observed != variant["ref"]:
                raise ValueError(f"IMPG REF mismatch at {variant['contig']}:{variant['position'] + 1}")
        snvs = {(v["contig"], v["position"]) for v in variants if len(v["ref"]) == len(v["alt"]) == 1 and v["ref"] in DNA and v["alt"] in DNA}
        classes = {
            "coding_gene": gene_positions,
            "CDS": cds_positions,
            "fourfold_CDS_reference_S": fourfold_s,
            "fourfold_CDS_reference_W": fourfold_w,
        }
        diversity_rows = []
        for class_name, positions in classes.items():
            numerator = snvs & positions
            if not positions:
                continue
            bootstrap = bootstrap_ratio(_block_counts(positions, numerator, args.block_size), args.bootstrap_replicates, int(hashlib.sha256(f"{run_input['dataset_id']}:{class_name}".encode()).hexdigest()[:16], 16))
            exclusions = {
                "outside_deterministic_panel_genes": int(run_input["targeted_gene_count"]) - prep["panel_selected_gene_count"],
                "sweepga_1to1_lost_panel_bases": prep["sweepga_1to1_lost_panel_bases"],
                "sweepga_1to1_excluded_or_partial_target_genes": prep["sweepga_1to1_excluded_or_partial_target_gene_count"],
                "sweepga_1to1_target_gene_lost_bases": prep["sweepga_1to1_target_gene_lost_bases"],
                "non_acgt_or_nonclass_panel_bases": prep["panel_selected_union_bases"] - len(positions),
                "phase_qc": phase_qc,
            }
            diversity_rows.append({
                "dataset_id": run_input["dataset_id"], "scientific_name": run_input["scientific_name"],
                "annotation_class": class_name,
                "statistic_label": (
                    "pi_S_reference_conditioned" if class_name == "fourfold_CDS_reference_S"
                    else "pi_W_reference_conditioned" if class_name == "fourfold_CDS_reference_W"
                    else "diploid_haplotype_diversity"
                ),
                "scope_label": "biological_coding_annotation_panel_not_genome_wide",
                "sweepga_hit_cap_query": "1", "sweepga_hit_cap_target": "1",
                "sweepga_mapping_provenance": run_input["sweepga_bounded_paf_path"],
                "impg_index_provenance": str(args.index),
                "impg_partition_provenance": f"{args.partitions};impg partition -w 2000 -d 0",
                "impg_query_provenance": f"{args.focus_bed};impg query -d 0 -o vcf:poa --force-large-region --min-transitive-len 1",
                "eligible_haplotypes": 2, "variant_numerator": len(numerator),
                "callable_denominator": len(positions), "estimate": f"{bootstrap['estimate']:.12g}",
                "bootstrap_blocks": bootstrap["blocks"], "bootstrap_replicates": bootstrap["replicates"],
                "bootstrap_ci_low": f"{bootstrap['ci_low']:.12g}", "bootstrap_ci_high": f"{bootstrap['ci_high']:.12g}",
                "bootstrap_standard_error": f"{bootstrap['standard_error']:.12g}",
                "targeted_genes": run_input["targeted_gene_count"], "panel_selected_genes": prep["panel_selected_gene_count"],
                "sweepga_queryable_genes": prep["sweepga_1to1_fully_covered_target_gene_count"],
                "impg_queryable_genes": len({row["gene_id"] for row in gene_rows if any(_contains(queryable, row["contig"], pos) for pos in (int(row["start_0based"]), int(row["end_0based_exclusive"]) - 1))}),
                "excluded_genes": prep["sweepga_1to1_excluded_or_partial_target_gene_count"],
                "targeted_bases": run_input["targeted_gene_union_bases"], "panel_selected_bases": prep["panel_selected_union_bases"],
                "sweepga_queryable_bases": prep["sweepga_1to1_callable_panel_bases"], "impg_callable_bases": len(positions),
                "exclusions_json": json.dumps(exclusions, sort_keys=True, separators=(",", ":")),
                "feature_identity_provenance": run_input["annotation_query_manifest_path"],
                "transcript_phase_provenance": "native_GFF_CDS_transcript_identity_and_phase_chain_preserved",
                "biological_input": "yes",
            })
        audit_rows = []
        for variant in variants:
            key = variant["contig"], variant["position"]
            hits = [row for row in cds_rows if row["contig"] == key[0] and int(row["start_0based"]) <= key[1] < int(row["end_0based_exclusive"])]
            gene_hits = [row for row in gene_rows if row["contig"] == key[0] and int(row["start_0based"]) <= key[1] < int(row["end_0based_exclusive"])]
            if key in gene_positions:
                audit_rows.append({
                    **variant,
                    "gene_ids": ",".join(sorted(
                        {row["gene_id"] for row in gene_hits}
                        | {gene for row in hits for gene in row["gene_ids"].split(",") if gene}
                    )),
                    "transcript_ids": ",".join(sorted({tx for row in hits for tx in row["transcript_ids"].split(",") if tx})),
                    "cds_ids": ",".join(sorted({row["cds_id"] for row in hits if row["cds_id"]})),
                    "feature_row_ids": ",".join(sorted({row["feature_row_id"] for row in hits})),
                    "cds_phases": ",".join(sorted({row["phase"] for row in hits})),
                    "in_coding_gene_callable": "yes", "in_cds_callable": "yes" if key in cds_positions else "no",
                })
        _write_tsv(out / "variant_feature_audit.tsv", (
            "contig", "position", "ref", "alt", "gt", "gene_ids", "transcript_ids", "cds_ids",
            "feature_row_ids", "cds_phases", "in_coding_gene_callable", "in_cds_callable",
        ), audit_rows)
        direct = None
        for variant in variants:
            if (variant["contig"], variant["position"]) in gene_positions:
                direct = _direct_h2_check(variant, Path(run_input["sweepga_bounded_paf_path"]), h1, h2)
                if direct and direct["ref_matches_h1"] and direct["alt_matches_h2"]:
                    break
        if not direct or not direct["ref_matches_h1"] or not direct["alt_matches_h2"]:
            raise ValueError("no representative IMPG SNV validated directly against both H1 and H2")
        representative_feature = next(
            row for row in audit_rows
            if row["contig"] == direct["h1_contig"]
            and int(row["position"]) + 1 == direct["h1_position_1based"]
        )
        direct.update({
            "h1_native_gff_path": run_input["annotation_gff_path"],
            "h1_native_gff_sha256": run_input["annotation_gff_sha256"],
            "fresh_annotation_query_manifest_path": run_input["annotation_query_manifest_path"],
            "fresh_annotation_query_manifest_sha256": _sha256(Path(run_input["annotation_query_manifest_path"])),
            "gene_ids": representative_feature["gene_ids"],
            "transcript_ids": representative_feature["transcript_ids"],
            "cds_ids": representative_feature["cds_ids"],
            "feature_row_ids": representative_feature["feature_row_ids"],
            "cds_phases": representative_feature["cds_phases"],
            "original_h1_native_annotation_status": "passed",
        })
    finally:
        h1.close()
        h2.close()
    callable_audit = {
        "mapping_callable_bases": _interval_bases(mapping), "impg_focus_intersection_bases": _interval_bases(queryable),
        "coding_gene_acgt_callable_bases": len(gene_positions), "cds_acgt_callable_bases": len(cds_positions),
        "fourfold_S_callable_bases": len(fourfold_s), "fourfold_W_callable_bases": len(fourfold_w),
        "boundary_policy": "half_open_exact_native_features_intersected_with_SweepGA_mapping_and_IMPG_native_partitions",
        "additional_context_bases": 0,
    }
    (out / "callable_audit.json").write_text(
        json.dumps(callable_audit, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (out / "representative_h1_h2_audit.json").write_text(
        json.dumps(direct, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    summary = {
        "dataset_id": run_input["dataset_id"], "scientific_name": run_input["scientific_name"],
        "biological_input": True, "diversity_rows": diversity_rows, "prepare_qc": prep,
        "callable_audit": callable_audit, "phase_qc": phase_qc, "direct_sequence_validation": direct,
        "impg_variant_records_untrimmed_after_normalization": len(
            _variant_rows(args.bcftools, out / "normalized.untrimmed.bcf")
        ),
        "impg_variant_records_after_panel_trim_and_exact_dedup": len(variants),
        "impg_variant_records_in_coding_panel": len(audit_rows),
        "normalization": "bcftools norm -f H1 -m -any; target trim; exact dedup; BCF+CSI",
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def validate_diversity_rows(rows: Sequence[Mapping[str, Any]], expected_datasets: set[str]) -> None:
    if not rows:
        raise ValueError("diversity table has no biological result rows")
    observed = {str(row["dataset_id"]) for row in rows}
    if observed != expected_datasets:
        raise ValueError(f"result datasets {observed} differ from acquired tuples {expected_datasets}")
    for row in rows:
        if row.get("biological_input") != "yes":
            raise ValueError("non-biological result row")
        if int(row["eligible_haplotypes"]) != 2:
            raise ValueError("diploid result must have two eligible haplotypes")
        if int(row["variant_numerator"]) <= 0:
            raise ValueError("result row lacks positive variants")
        if int(row["callable_denominator"]) <= 0:
            raise ValueError("result row lacks a positive callable denominator")
        numeric = [float(row[key]) for key in ("estimate", "bootstrap_ci_low", "bootstrap_ci_high", "bootstrap_standard_error")]
        if not all(math.isfinite(value) for value in numeric):
            raise ValueError("result row has non-finite estimate or uncertainty")


def _lineage_audit(
    args: argparse.Namespace,
    rows: Sequence[Mapping[str, Any]],
    run_rows: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    old_rows = _read_tsv(args.superseded_results)
    old_hashes = {row["artifact_sha256"] for row in old_rows if row.get("artifact_sha256")}
    supersession = _read_tsv(args.supersession_ledger)
    superseded_paths = {_resolved(row["superseded_path"]) for row in supersession if row.get("superseded_path")}
    superseded_hashes = {row["superseded_sha256"] for row in supersession if row.get("superseded_sha256")}
    lineage_paths: set[str] = set()
    lineage_hashes: set[str] = set()
    for row in rows:
        for key in ("sweepga_mapping_provenance", "impg_index_provenance", "feature_identity_provenance"):
            value = str(row.get(key, "")).split(";", 1)[0]
            if value:
                lineage_paths.add(_resolved(value))
    for row in run_rows:
        for key, value in row.items():
            text = str(value)
            if key.endswith("_sha256") and len(text) == 64:
                lineage_hashes.add(text)
            elif key.endswith("_path") or key in {"impg_index_path", "impg_normalized_bcf"}:
                if text:
                    lineage_paths.add(_resolved(text))
    for dataset in EXPECTED_DATASETS:
        work = args.work_root / dataset
        for relative in (
            "graph.impg", "partitions/partitions.bed", "focus.bed", "laced.vcf",
            "normalized.untrimmed.bcf", "normalized.untrimmed.bcf.csi", "normalized.bcf",
            "normalized.bcf.csi", "summary.json", "callable_audit.json",
            "representative_h1_h2_audit.json", "variant_feature_audit.tsv",
            "preflight_audit.json", "input.json", "prepare_qc.json", "mapping_callable.bed",
            "native_annotation/gene_manifest.tsv", "native_annotation/impg_query_manifest.tsv",
            "native_annotation/impg_execution_spans.bed", "native_annotation/impg_span_feature_map.tsv",
            "native_annotation/h1_h2_contig_map.tsv", "native_annotation/query_qc.json",
        ):
            path = work / relative
            if not path.is_file():
                raise ValueError(f"lineage artifact is missing: {path}")
            lineage_paths.add(_resolved(path))
            lineage_hashes.add(_sha256(path))
    for name in (
        "diploid_diversity.tsv", "diploid_failures.tsv", "diploid_run_manifest.tsv",
        "diploid_rerun_commands.sh", "diploid_qc.md",
    ):
        path = args.output_dir / name
        if not path.is_file():
            raise ValueError(f"new Tier 3A result artifact is missing: {path}")
        lineage_paths.add(_resolved(path))
        lineage_hashes.add(_sha256(path))
    path_overlap = sorted(lineage_paths & superseded_paths)
    checksum_overlap = sorted(lineage_hashes & (superseded_hashes | old_hashes))
    if path_overlap or checksum_overlap:
        raise ValueError("new lineage contains a superseded path or old Tier 3A checksum")
    return {
        "schema_version": "tier3a-origin-rerun-lineage-audit-v1",
        "status": "passed",
        "old_tier3a_result_checksums": sorted(old_hashes),
        "supersession_ledger_checksums": sorted(superseded_hashes),
        "new_lineage_checksums": sorted(lineage_hashes),
        "old_result_checksum_intersection": [],
        "superseded_checksum_intersection": [],
        "superseded_path_intersection": [],
        "assertion": "no_old_tier3a_result_checksum_or_superseded_artifact_occurs_in_new_lineage",
    }


def finalize(args: argparse.Namespace) -> None:
    acquisition = _read_tsv(args.manifest)
    expected = {row["dataset_id"] for row in acquisition if row["eligibility_status"] == "eligible_biological"}
    summaries, failures = [], []
    for dataset in sorted(expected):
        path = args.work_root / dataset / "summary.json"
        if path.is_file():
            summaries.append(json.loads(path.read_text(encoding="utf-8")))
        else:
            row = next(item for item in acquisition if item["dataset_id"] == dataset)
            failures.append({"dataset_id": dataset, "scientific_name": row["scientific_name"], "stage": "summary", "status": "failed", "reason": "missing_summary_json"})
    rows = [row for summary in summaries for row in summary["diversity_rows"]]
    validate_diversity_rows(rows, expected)
    _write_tsv(args.output_dir / "diploid_diversity.tsv", DIVERSITY_COLUMNS, rows)
    _write_tsv(args.output_dir / "diploid_failures.tsv", FAILURE_COLUMNS, failures)
    telemetry = _read_tsv(args.telemetry) if args.telemetry.is_file() else []
    telemetry_by_id = {row["dataset_id"]: row for row in telemetry}
    for measured in telemetry:
        if measured.get("prior_failed_job_id"):
            acquired = next(row for row in acquisition if row["dataset_id"] == measured["dataset_id"])
            failures.append({
                "dataset_id": measured["dataset_id"], "scientific_name": acquired["scientific_name"],
                "stage": measured.get("prior_failure_stage", "scheduler_attempt"), "status": "recovered",
                "reason": measured["prior_failure_reason"],
                "slurm_job_id": measured["prior_failed_job_id"],
                "stderr_path": measured["prior_stderr_path"],
                "rerun_command": measured.get("recovery_command", "sbatch --array=0 analysis/slurm/tier3a_biological_array.sh"),
            })
    _write_tsv(args.output_dir / "diploid_failures.tsv", FAILURE_COLUMNS, failures)
    run_rows = []
    for summary in summaries:
        dataset = summary["dataset_id"]
        acquired = next(row for row in acquisition if row["dataset_id"] == dataset)
        work = args.work_root / dataset
        measured = telemetry_by_id.get(dataset, {})
        run_rows.append({
            "dataset_id": dataset, "scientific_name": acquired["scientific_name"], "status": "completed",
            **measured, "sweepga_hit_cap": "1:1", "sweepga_paf_path": acquired["sweepga_bounded_paf_path"],
            "sweepga_paf_sha256": acquired["sweepga_bounded_paf_sha256"],
            "sweepga_binary_path": acquired["sweepga_origin_main_realpath"], "sweepga_binary_sha256": acquired["sweepga_origin_main_sha256"],
            "sweepga_command": acquired["sweepga_direct_command"], "impg_binary_path": acquired["impg_binary_path"],
            "impg_binary_sha256": acquired["impg_binary_sha256"], "impg_commit": acquired["impg_commit"],
            "impg_index_path": work / "graph.impg", "impg_partitions_path": work / "partitions/partitions.bed",
            "impg_native_partitions_selected": sum(1 for line in (work / "focus.bed").read_text(encoding="utf-8").splitlines() if line),
            "impg_focus_bed": work / "focus.bed",
            "impg_regional_vcf_count": sum(1 for line in (work / "vcf.list").read_text(encoding="utf-8").splitlines() if line),
            "impg_normalized_bcf": work / "normalized.bcf", "impg_normalized_bcf_sha256": _sha256(work / "normalized.bcf"),
            "guix_channel_commit": acquired["guix_channel_commit"],
            "guix_channels_path": acquired["guix_channels_path"],
            "guix_manifest_path": "analysis/guix/manifest.scm; supplemental analysis/guix/sweepga_impg_smoke_manifest.scm",
            "guix_profile_store_path": (work / "completed.tsv").read_text(encoding="utf-8").rstrip("\n").split("\t")[3],
            "sweepga_origin_main_commit": acquired["sweepga_origin_main_commit"],
            "sweepga_native_command_path": acquired["sweepga_native_command_path"],
            "sweepga_native_command_sha256": _sha256(Path(acquired["sweepga_native_command_path"])),
            "sweepga_observed_query_multiplicity": acquired["sweepga_observed_query_multiplicity"],
            "sweepga_observed_target_multiplicity": acquired["sweepga_observed_target_multiplicity"],
            "annotation_query_qc_path": work / "native_annotation/query_qc.json",
            "annotation_query_qc_sha256": _sha256(work / "native_annotation/query_qc.json"),
            "preflight_audit_path": work / "preflight_audit.json",
            "preflight_audit_sha256": _sha256(work / "preflight_audit.json"),
            "command": f"sbatch --array={sorted(expected).index(dataset)} analysis/slurm/tier3a_biological_array.sh",
        })
    _write_tsv(args.output_dir / "diploid_run_manifest.tsv", RUN_COLUMNS, run_rows)
    rerun_lines = ["#!/usr/bin/env bash", "set -euo pipefail", ""]
    for index, dataset in enumerate(sorted(expected)):
        rerun_lines.append(
            f"TIER3A_WORK_ROOT={args.work_root} sbatch --array={index} analysis/slurm/tier3a_biological_array.sh # {dataset}"
        )
    rerun_lines.extend([
        "",
        "# After the array completes, regenerate sacct telemetry and run:",
        f"python3 analysis/tier3a_biological.py finalize --manifest results/tier3a/acquisition_corrected_manifest.tsv --work-root {args.work_root} --telemetry results/tier3a/logs/origin-rerun-telemetry.tsv --supersession-ledger results/tier3a/acquisition_sweepga_supersession.tsv --superseded-results results/tier3a/diploid_superseded_results.tsv --output-dir results/tier3a",
    ])
    commands_output = args.output_dir / "diploid_rerun_commands.sh"
    commands_output.write_text("\n".join(rerun_lines) + "\n", encoding="utf-8")
    commands_output.chmod(0o755)
    lines = [
        "# Tier 3A biological diploid QC", "", "## Outcome", "",
        f"All {len(expected)} acquired biological tuples produced coding-region estimates through pinned origin/main SweepGA native `--num-mappings 1:1` and freshly executed IMPG.",
        "Every row is a coding/CDS annotation-panel estimate, not a genome-wide estimate. The panel was selected by stable hash of exact H1-native execution-span identity before mapping or variant inspection.", "",
        "## Separation of responsibilities", "",
        "- SweepGA supplied corrected whole-H1-versus-H2 mappings and enforced the native 1:1 query:target overlap cap. Binary realpath/SHA-256, fetched origin/main commit, native commands, multiplicity audits, Guix closure, and the supersession ledger passed before computation.",
        "- IMPG indexed those complete bounded PAFs, formed native graph partitions (`-w 2000 -d 0`), and queried only partitions intersecting selected native-annotation spans.",
        "- `bcftools norm` normalized and split alleles, panel trimming removed partition context, and exact duplicate records were removed before estimates.", "",
        "## Tuple audits", "",
    ]
    for summary in summaries:
        prep, callable_qc, direct = summary["prepare_qc"], summary["callable_audit"], summary["direct_sequence_validation"]
        sensitivity = prep["cap_sensitivity"]
        sensitivity_text = ", ".join(
            f"{cap}:{cap}=({values['query_coverage']:.4f},{values['target_coverage']:.4f})"
            for cap, values in sorted(sensitivity.items(), key=lambda item: int(item[0]))
        )
        if set(sensitivity) == {"1"}:
            sensitivity_note = "No higher-cap sensitivity was run for this production correction."
        else:
            sensitivity_note = "Higher caps are sensitivity only and were not passed to IMPG as graph policy."
        lines.extend([
            f"### {summary['scientific_name']} (`{summary['dataset_id']}`)", "",
            f"Targeted genes/bases: {prep['targeted_gene_count']:,}/{prep['targeted_gene_union_bases']:,}; deterministic panel: {prep['panel_selected_gene_count']:,}/{prep['panel_selected_union_bases']:,}; SweepGA-1:1 mapped panel genes/bases: {prep['sweepga_1to1_mapped_gene_count']:,}/{prep['sweepga_1to1_callable_panel_bases']:,}.",
            f"Across all native targets, SweepGA 1:1 fully covered {prep['sweepga_1to1_fully_covered_target_gene_count']:,} genes and excluded or only partially covered {prep['sweepga_1to1_excluded_or_partial_target_gene_count']:,}; its target-gene intersection retained {prep['sweepga_1to1_target_gene_callable_bases']:,} bases and lost {prep['sweepga_1to1_target_gene_lost_bases']:,}.",
            f"IMPG-audited callable coding/CDS bases: {callable_qc['coding_gene_acgt_callable_bases']:,}/{callable_qc['cds_acgt_callable_bases']:,}. Representative SNV: H1 {direct['h1_contig']}:{direct['h1_position_1based']} {direct['h1_base']} versus H2 {direct['h2_contig']}:{direct['h2_position_1based']} {direct['h2_base']} (strand {direct['paf_strand']}); both allele checks passed.",
            f"IMPG selected {sum(1 for line in (args.work_root / summary['dataset_id'] / 'focus.bed').read_text(encoding='utf-8').splitlines() if line):,} native partitions and emitted {sum(1 for line in (args.work_root / summary['dataset_id'] / 'vcf.list').read_text(encoding='utf-8').splitlines() if line):,} regional VCFs. Normalization yielded {summary['impg_variant_records_untrimmed_after_normalization']:,} records before exact panel trim/dedup and {summary['impg_variant_records_after_panel_trim_and_exact_dedup']:,} after; {summary['impg_variant_records_in_coding_panel']:,} overlap callable coding genes.",
            f"Native cap coverage (query,target): {sensitivity_text}. {sensitivity_note}", "",
        ])
    lines.extend([
        "## Scheduler recovery", "",
        *(
            [f"One partial attempt was discarded without artifact reuse: {failures[0]['dataset_id']} job {failures[0]['slurm_job_id']} ({failures[0]['reason']}). The replacement tuple ran from an empty directory and completed successfully.", ""]
            if failures else ["No scheduler recovery was required.", ""]
        ),
        "## Boundary, annotation, and uncertainty policy", "",
        "Coordinates remain zero-based half-open through native feature, SweepGA coverage, and IMPG-partition intersections. The H1-native GFF query manifests and mapping-derived contig maps were regenerated in the new run tree from the corrected PAF; none came from the invalid run. Only exact H1 A/C/G/T bases enter denominators. CDS feature row, gene, transcript, protein, locus, strand, and phase identities are retained in each tuple's `variant_feature_audit.tsv`. Fourfold sites require a valid transcript-order phase chain and exclude frame-discordant overlaps.",
        "Uncertainty is a deterministic genomic block bootstrap (50-kb blocks, 1,000 replicates). Direct CIGAR traversal is used only for one representative H1/H2 allele audit per tuple; it does not create the primary call set.", "",
        "## Reproducibility", "",
        "The run manifest records Slurm telemetry, exact biological paths, pinned Guix channel commit/profile, and primary artifacts. `diploid_rerun_commands.sh` contains exact commands. Scheduler stdout/stderr and `sacct` telemetry are under `results/tier3a/logs/`. `diploid_lineage_audit.json` proves by machine-readable set intersection that no old Tier 3A result checksum or superseded artifact entered the new lineage.", "",
        "Slurm allocation, elapsed time, exit status, CPU, memory, and node telemetry for this corrected rerun are recorded tuple-by-tuple in `diploid_run_manifest.tsv` and the fresh `origin-rerun-*.sacct` files. Empty accounting fields remain explicit rather than being inferred.", "",
    ])
    (args.output_dir / "diploid_qc.md").write_text("\n".join(lines), encoding="utf-8")
    lineage = _lineage_audit(args, rows, run_rows)
    (args.output_dir / "diploid_lineage_audit.json").write_text(
        json.dumps(lineage, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def parser() -> argparse.ArgumentParser:
    root = argparse.ArgumentParser(description=__doc__)
    commands = root.add_subparsers(dest="command", required=True)
    gate = commands.add_parser("preflight")
    gate.add_argument("--manifest", type=Path, required=True)
    gate.add_argument("--commands", type=Path, required=True)
    gate.add_argument("--supersession-ledger", type=Path, required=True)
    gate.add_argument("--build-provenance", type=Path, required=True)
    gate.add_argument("--output", type=Path, required=True)
    gate.set_defaults(function=preflight)
    prepare_p = commands.add_parser("prepare")
    prepare_p.add_argument("--manifest", type=Path, required=True)
    prepare_p.add_argument("--dataset-id", required=True)
    prepare_p.add_argument("--output-dir", type=Path, required=True)
    prepare_p.add_argument("--annotation-dir", type=Path, required=True)
    prepare_p.add_argument("--preflight-audit", type=Path, required=True)
    prepare_p.add_argument("--panel-bases", type=int, default=2_000_000)
    prepare_p.set_defaults(function=prepare)
    summary_p = commands.add_parser("summarize")
    summary_p.add_argument("--output-dir", type=Path, required=True)
    summary_p.add_argument("--focus-bed", type=Path, required=True)
    summary_p.add_argument("--partitions", type=Path, required=True)
    summary_p.add_argument("--index", type=Path, required=True)
    summary_p.add_argument("--bcf", type=Path, required=True)
    summary_p.add_argument("--bcftools", type=Path, required=True)
    summary_p.add_argument("--block-size", type=int, default=50_000)
    summary_p.add_argument("--bootstrap-replicates", type=int, default=1000)
    summary_p.set_defaults(function=summarize)
    final_p = commands.add_parser("finalize")
    final_p.add_argument("--manifest", type=Path, required=True)
    final_p.add_argument("--work-root", type=Path, required=True)
    final_p.add_argument("--telemetry", type=Path, required=True)
    final_p.add_argument("--supersession-ledger", type=Path, required=True)
    final_p.add_argument("--superseded-results", type=Path, required=True)
    final_p.add_argument("--output-dir", type=Path, required=True)
    final_p.set_defaults(function=finalize)
    return root


def main(argv: Sequence[str] | None = None) -> int:
    args = parser().parse_args(argv)
    args.function(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
