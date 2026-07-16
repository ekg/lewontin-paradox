#!/usr/bin/env python3
"""Audit and publish the provenance-locked Tier 3A SweepGA correction.

The production selector in this module accepts only the byte-reproducible
origin/main executable and native ``--num-mappings 1:1`` mappings.  Publication
is all-or-nothing: three completed SweepGA -> IMPG result directories are
validated before corrected manifests replace their temporary siblings.
"""

from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import json
import os
import shlex
import sys
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


EXPECTED_DATASETS = {
    "menidia_menidia_fMenMen1",
    "spinachia_spinachia_SK-2024b",
    "tautogolabrus_adspersus_fTauAds1",
}
SWEEPGA_COMMIT = "018e4ce49d2c125820e0ac50dc5feaa02d423683"
SWEEPGA_SHA256 = "fa7f0edb9b7e275c288db254046020e136d4267dd5ee043379227ef80da0573b"
LEGACY_SWEEPGA_SHA256 = "1a5440529f5eff91cb7d82876a83a5282df66fb5e2c4b1a9c6caa0bdb83de7b1"
REQUIRED_COMMAND_TOKEN = "--num-mappings 1:1"
CORRECTION_COLUMNS = (
    "tier3a_correction_status",
    "tier3a_correction_schema",
    "sweepga_origin_main_realpath",
    "sweepga_origin_main_sha256",
    "sweepga_origin_main_commit",
    "sweepga_origin_build_provenance_path",
    "sweepga_origin_build_provenance_sha256",
    "sweepga_native_command_path",
    "sweepga_native_multiplicity_audit_path",
    "sweepga_native_multiplicity_audit_sha256",
    "sweepga_observed_query_multiplicity",
    "sweepga_observed_target_multiplicity",
    "sweepga_mapping_records",
    "sweepga_query_contigs_covered",
    "sweepga_target_contigs_covered",
    "sweepga_mapping_slurm_job_id",
    "sweepga_mapping_slurm_array_task_id",
    "sweepga_mapping_slurm_state",
    "sweepga_mapping_slurm_elapsed",
    "impg_corrected_index_path",
    "impg_corrected_partitions_path",
    "impg_corrected_focus_bed_path",
    "impg_corrected_normalized_bcf_path",
    "impg_corrected_normalized_bcf_sha256",
    "impg_corrected_normalized_bcf_csi_path",
    "impg_corrected_normalized_vcf_gz_path",
    "impg_corrected_normalized_vcf_gz_sha256",
    "impg_corrected_normalized_vcf_tbi_path",
    "impg_positive_variant_records",
    "impg_callable_denominator_bases",
    "representative_variant_audit_path",
    "representative_variant_h1_h2_status",
    "slurm_job_id",
    "slurm_array_task_id",
    "slurm_node",
    "slurm_requested_cpus",
    "slurm_requested_memory",
    "slurm_elapsed_seconds_internal",
    "slurm_state",
    "slurm_elapsed",
    "slurm_total_cpu",
    "slurm_max_rss",
    "scheduler_sacct_path",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    csv.field_size_limit(sys.maxsize)
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, columns: Sequence[str], rows: Iterable[Mapping[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row.get(column, "") for column in columns})
        handle.flush()
        os.fsync(handle.fileno())
    temporary.replace(path)


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    temporary.replace(path)


def fai_lengths(path: Path) -> dict[str, int]:
    result = {}
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            fields = raw.rstrip("\n").split("\t")
            result[fields[0]] = int(fields[1])
    if not result:
        raise ValueError(f"empty FASTA index: {path}")
    return result


def merged_bases(intervals: Iterable[tuple[int, int]]) -> int:
    ordered = sorted(intervals)
    if not ordered:
        return 0
    total, start, end = 0, *ordered[0]
    for left, right in ordered[1:]:
        if left <= end:
            end = max(end, right)
        else:
            total += end - start
            start, end = left, right
    return total + end - start


def axis_audit(groups: Mapping[str, list[tuple[int, int, int]]], threshold: float) -> dict[str, int]:
    """Describe retained intervals; this is not the native cap decision rule."""
    violating_pairs = 0
    max_raw_depth = 0
    max_threshold_multiplicity = 1 if groups else 0
    for intervals in groups.values():
        active: list[tuple[int, int, int]] = []
        for start, end, record_index in sorted(intervals):
            active = [item for item in active if item[1] > start]
            max_raw_depth = max(max_raw_depth, len(active) + 1)
            conflicts = 0
            for other_start, other_end, _ in active:
                overlap = min(end, other_end) - max(start, other_start)
                denominator = min(end - start, other_end - other_start)
                if denominator > 0 and overlap / denominator > threshold:
                    violating_pairs += 1
                    conflicts += 1
            max_threshold_multiplicity = max(max_threshold_multiplicity, conflicts + 1)
            active.append((start, end, record_index))
    return {
        "raw_max_interval_depth": max_raw_depth,
        "overlap_threshold_violating_pairs": violating_pairs,
        "observed_threshold_multiplicity": max_threshold_multiplicity,
    }


def audit_mapping(args: argparse.Namespace) -> None:
    h1, h2 = fai_lengths(args.h1_fai), fai_lengths(args.h2_fai)
    query_intervals: dict[str, list[tuple[int, int, int]]] = collections.defaultdict(list)
    target_intervals: dict[str, list[tuple[int, int, int]]] = collections.defaultdict(list)
    pair_rows: dict[tuple[str, str], dict[str, Any]] = {}
    record_count = 0
    with args.paf.open(encoding="utf-8") as handle:
        for number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 12:
                raise ValueError(f"truncated PAF row {number}")
            qname, qlen, qstart, qend = fields[0], int(fields[1]), int(fields[2]), int(fields[3])
            tname, tlen, tstart, tend = fields[5], int(fields[6]), int(fields[7]), int(fields[8])
            if h1.get(qname) != qlen or h2.get(tname) != tlen:
                raise ValueError(f"PAF row {number} does not use H1 query and H2 target axes")
            if not 0 <= qstart < qend <= qlen or not 0 <= tstart < tend <= tlen:
                raise ValueError(f"invalid PAF coordinates at row {number}")
            record_count += 1
            query_intervals[qname].append((qstart, qend, record_count))
            target_intervals[tname].append((tstart, tend, record_count))
            pair = pair_rows.setdefault((qname, tname), {
                "h1_contig": qname,
                "h1_length": qlen,
                "h2_contig": tname,
                "h2_length": tlen,
                "strands": set(),
                "mapping_records": 0,
                "h1": [],
                "h2": [],
            })
            pair["strands"].add(fields[4])
            pair["mapping_records"] += 1
            pair["h1"].append((qstart, qend))
            pair["h2"].append((tstart, tend))
    if not record_count:
        raise ValueError("production PAF is empty")
    query_audit = axis_audit(query_intervals, args.overlap_threshold)
    target_audit = axis_audit(target_intervals, args.overlap_threshold)
    def paf_core(path: Path) -> list[tuple[str, ...]]:
        with path.open(encoding="utf-8") as handle:
            return [tuple(raw.rstrip("\n").split("\t")[:12]) for raw in handle if raw.strip() and not raw.startswith("#")]

    recheck_core = paf_core(args.native_recheck_paf)
    production_core = paf_core(args.paf)
    if not recheck_core or production_core != recheck_core:
        raise ValueError("native 1:1 recheck changed ordered mandatory PAF fields")
    query_bases = sum(merged_bases((start, end) for start, end, _ in values) for values in query_intervals.values())
    target_bases = sum(merged_bases((start, end) for start, end, _ in values) for values in target_intervals.values())
    mapping_sha = sha256(args.paf)
    audit = {
        "schema_version": "tier3a-origin-main-native-multiplicity-audit-v1",
        "status": "passed",
        "audit_mode": "pinned_native_1to1_recheck_read_only_no_production_replacement",
        "paf_path": str(args.paf.resolve()),
        "paf_sha256": mapping_sha,
        "record_count": record_count,
        "h1_axis": "query",
        "h2_axis": "target",
        "overlap_threshold": args.overlap_threshold,
        "native_num_mappings": "1:1",
        "native_recheck_paf_path": str(args.native_recheck_paf.resolve()),
        "native_recheck_paf_sha256": sha256(args.native_recheck_paf),
        "native_recheck_ordered_mandatory_paf_fields_identical": True,
        "observed_native_query_multiplicity_cap": 1,
        "observed_native_target_multiplicity_cap": 1,
        "query_axis": query_audit,
        "target_axis": target_audit,
        "query_contigs_covered": len(query_intervals),
        "target_contigs_covered": len(target_intervals),
        "query_union_bases": query_bases,
        "target_union_bases": target_bases,
        "query_total_bases": sum(h1.values()),
        "target_total_bases": sum(h2.values()),
        "query_coverage": query_bases / sum(h1.values()),
        "target_coverage": target_bases / sum(h2.values()),
    }
    write_json(args.output, audit)
    contig_rows = []
    for index, pair in enumerate(sorted(pair_rows.values(), key=lambda item: (item["h1_contig"], item["h2_contig"])), 1):
        contig_rows.append({
            "relation_id": f"contig_relation_{index:06d}",
            "h1_contig": pair["h1_contig"],
            "h1_length": pair["h1_length"],
            "h2_contig": pair["h2_contig"],
            "h2_length": pair["h2_length"],
            "strands": ",".join(sorted(pair["strands"])),
            "mapping_records": pair["mapping_records"],
            "h1_aligned_union_bases": merged_bases(pair["h1"]),
            "h2_aligned_union_bases": merged_bases(pair["h2"]),
            "status": "origin_main_native_1to1_mapping",
        })
    write_tsv(args.contig_map, (
        "relation_id", "h1_contig", "h1_length", "h2_contig", "h2_length", "strands",
        "mapping_records", "h1_aligned_union_bases", "h2_aligned_union_bases", "status",
    ), contig_rows)


def lookup_base(path: Path, dataset: str) -> tuple[list[str], dict[str, str]]:
    rows = read_tsv(path)
    matches = [row for row in rows if row["dataset_id"] == dataset]
    if len(matches) != 1:
        raise ValueError(f"expected one base-manifest row for {dataset}")
    return list(rows[0]), matches[0]


def corrected_row(base: dict[str, str], work: Path, stage: bool) -> dict[str, Any]:
    mapping = work / "mapping/production.1to1.paf"
    command_path = work / "mapping/command.txt"
    audit_path = work / "mapping/native_multiplicity_audit.json"
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    command = command_path.read_text(encoding="utf-8").strip()
    if REQUIRED_COMMAND_TOKEN not in command:
        raise ValueError(f"missing exact native option in command for {base['dataset_id']}")
    if base.get("sweepga_binary_sha256") == SWEEPGA_SHA256:
        raise ValueError("base manifest unexpectedly already identifies the correction")
    query_dir = work / "annotation_query"
    query_qc_path = query_dir / "query_qc.json"
    query_qc = json.loads(query_qc_path.read_text(encoding="utf-8"))
    artifacts = query_qc["artifacts"]
    row: dict[str, Any] = dict(base)
    repo_root = Path(__file__).resolve().parents[1]
    origin_manifest = repo_root / "analysis/guix/sweepga_origin_main_manifest.scm"
    origin_build_provenance = repo_root / "analysis/sweepga_origin_main_build.json"
    for field in (
        "raw_wfmash_paf_path", "raw_wfmash_paf_sha256", "callable_bed_path", "callable_bed_sha256",
        "callable_bases", "heterozygous_snvs", "normalized_bcf_path", "normalized_bcf_sha256",
        "normalized_bcf_csi_path", "normalized_bcf_csi_sha256", "callable_construction",
        "callable_exclusions_json", "sweepga_direct_cap_recheck_path", "sweepga_direct_cap_recheck_sha256",
        "sweepga_cap_sensitivity_path", "sweepga_cap_sensitivity_sha256",
        "sweepga_cap5_query_coverage", "sweepga_cap5_target_coverage",
        "sweepga_cap10_query_coverage", "sweepga_cap10_target_coverage",
        "wfmash_sensitivity_bcf_path", "wfmash_sensitivity_bcf_sha256",
        "smoke_qc_path", "smoke_qc_sha256", "smoke_callable_bases", "smoke_parser_status",
        "handoff_validation_path", "handoff_validation_sha256", "impg_truth_regions_path",
        "impg_truth_regions_sha256", "impg_laced_vcf_path", "impg_laced_vcf_sha256",
        "impg_normalized_bcf_path", "impg_normalized_bcf_sha256",
        "impg_normalized_bcf_csi_path", "impg_normalized_bcf_csi_sha256",
        "guix_tool_versions_path", "guix_tool_versions_sha256",
        "staged_object_inventory_path", "staged_object_inventory_sha256",
    ):
        row[field] = ""
    row.update({
        "modality": "corrected_origin_main_SweepGA_to_IMPG_native_partition_query",
        "primary_mapping_engine": "SweepGA unmodified origin/main reproducible Guix build",
        "primary_mapping_policy": "whole_H1_and_H2;native_--num-mappings_1:1;no_posthoc_mapping_replacement",
        "sweepga_selected_cap_query": "1",
        "sweepga_selected_cap_target": "1",
        "sweepga_bounded_paf_path": str(mapping.resolve()),
        "sweepga_bounded_paf_sha256": audit["paf_sha256"],
        "sweepga_direct_command": command,
        "sweepga_direct_threads": os.environ.get("SLURM_CPUS_PER_TASK", "8"),
        "sweepga_direct_query_coverage": f"{audit['query_coverage']:.12f}",
        "sweepga_direct_target_coverage": f"{audit['target_coverage']:.12f}",
        "sweepga_direct_raw_max_query_interval_depth": audit["query_axis"]["raw_max_interval_depth"],
        "sweepga_direct_raw_max_target_interval_depth": audit["target_axis"]["raw_max_interval_depth"],
        "sweepga_direct_cap_fixed_point_status": "not_used_native_output_read_only_multiplicity_audit_passed",
        "sweepga_h1_paf_axis": "query",
        "sweepga_h2_paf_axis": "target",
        "sweepga_cap1_query_coverage": f"{audit['query_coverage']:.12f}",
        "sweepga_cap1_target_coverage": f"{audit['target_coverage']:.12f}",
        "sweepga_commit": SWEEPGA_COMMIT,
        "sweepga_binary_path": str(Path(os.environ.get("SWEEPGA_PINNED_REALPATH", "/moosefs/erikg/tier3scratch/sweepga-origin-main-018e4ce/bin-1/sweepga")).resolve()),
        "sweepga_binary_sha256": SWEEPGA_SHA256,
        "guix_manifest_path": "analysis/guix/sweepga_origin_main_manifest.scm",
        "guix_manifest_sha256": sha256(origin_manifest),
        "guix_profile_store_path": "/gnu/store/yfffyhdm3a9bsah4gzw9dzri623af3f6-profile",
        "guix_toolchain_manifest_path": "analysis/guix/sweepga_origin_main_manifest.scm",
        "guix_toolchain_manifest_sha256": sha256(origin_manifest),
        "production_variant_status": "staged" if stage else "completed_corrected_origin_main_IMPG",
        "annotation_gene_manifest_path": artifacts["gene_manifest.tsv"]["path"],
        "annotation_gene_manifest_sha256": artifacts["gene_manifest.tsv"]["sha256"],
        "annotation_query_manifest_path": artifacts["impg_query_manifest.tsv"]["path"],
        "annotation_query_manifest_sha256": artifacts["impg_query_manifest.tsv"]["sha256"],
        "annotation_execution_spans_path": artifacts["impg_execution_spans.bed"]["path"],
        "annotation_execution_spans_sha256": artifacts["impg_execution_spans.bed"]["sha256"],
        "annotation_span_feature_map_path": artifacts["impg_span_feature_map.tsv"]["path"],
        "annotation_span_feature_map_sha256": artifacts["impg_span_feature_map.tsv"]["sha256"],
        "haplotype_contig_map_path": artifacts["h1_h2_contig_map.tsv"]["path"],
        "haplotype_contig_map_sha256": artifacts["h1_h2_contig_map.tsv"]["sha256"],
        "haplotype_contig_map_rows": sum(1 for line in Path(artifacts["h1_h2_contig_map.tsv"]["path"]).open() if line.strip()) - 1,
        "annotation_query_qc_path": str(query_qc_path.resolve()),
        "annotation_query_qc_sha256": sha256(query_qc_path),
        "targeted_gene_count": query_qc["targeted_gene_count"],
        "targeted_gene_union_bases": query_qc["targeted_gene_union_bases"],
        "queryable_gene_count": query_qc["queryable_gene_count"],
        "queryable_gene_union_bases": query_qc["queryable_gene_union_bases"],
        "excluded_gene_count": query_qc["excluded_gene_count"],
        "excluded_gene_union_bases": query_qc["excluded_gene_union_bases"],
        "queryable_CDS_rows": query_qc["queryable_CDS_rows"],
        "excluded_CDS_rows": query_qc["excluded_CDS_rows"],
        "queryable_CDS_phase_counts_json": json.dumps(query_qc["queryable_CDS_phase_counts"], sort_keys=True, separators=(",", ":")),
        "annotation_execution_span_count": query_qc["execution_span_count"],
        "annotation_execution_span_union_bases": query_qc["execution_span_union_bases"],
    })
    if not stage:
        summary_path = work / "summary.json"
        success_path = work / "success.json"
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        success = json.loads(success_path.read_text(encoding="utf-8"))
        direct = summary["direct_sequence_validation"]
        bcf = work / "normalized.bcf"
        vcf = work / "normalized.vcf.gz"
        denominators = [int(item["callable_denominator"]) for item in summary["diversity_rows"]]
        variants = [int(item["variant_numerator"]) for item in summary["diversity_rows"]]
        row.update({
            "tier3a_correction_status": "current_production",
            "tier3a_correction_schema": "tier3a-origin-main-correction-v1",
            "sweepga_origin_main_realpath": str(Path(row["sweepga_binary_path"]).resolve()),
            "sweepga_origin_main_sha256": SWEEPGA_SHA256,
            "sweepga_origin_main_commit": SWEEPGA_COMMIT,
            "sweepga_origin_build_provenance_path": "analysis/sweepga_origin_main_build.json",
            "sweepga_origin_build_provenance_sha256": sha256(origin_build_provenance),
            "sweepga_native_command_path": str(command_path.resolve()),
            "sweepga_native_multiplicity_audit_path": str(audit_path.resolve()),
            "sweepga_native_multiplicity_audit_sha256": sha256(audit_path),
            "sweepga_observed_query_multiplicity": audit["observed_native_query_multiplicity_cap"],
            "sweepga_observed_target_multiplicity": audit["observed_native_target_multiplicity_cap"],
            "sweepga_mapping_records": audit["record_count"],
            "sweepga_query_contigs_covered": audit["query_contigs_covered"],
            "sweepga_target_contigs_covered": audit["target_contigs_covered"],
            "impg_corrected_index_path": str((work / "graph.impg").resolve()),
            "impg_corrected_partitions_path": str((work / "partitions/partitions.bed").resolve()),
            "impg_corrected_focus_bed_path": str((work / "focus.bed").resolve()),
            "impg_corrected_normalized_bcf_path": str(bcf.resolve()),
            "impg_corrected_normalized_bcf_sha256": sha256(bcf),
            "impg_corrected_normalized_bcf_csi_path": str(Path(str(bcf) + ".csi").resolve()),
            "impg_corrected_normalized_vcf_gz_path": str(vcf.resolve()),
            "impg_corrected_normalized_vcf_gz_sha256": sha256(vcf),
            "impg_corrected_normalized_vcf_tbi_path": str(Path(str(vcf) + ".tbi").resolve()),
            "impg_positive_variant_records": min(variants),
            "impg_callable_denominator_bases": min(denominators),
            "representative_variant_audit_path": str(summary_path.resolve()),
            "representative_variant_h1_h2_status": "passed" if direct["ref_matches_h1"] and direct["alt_matches_h2"] else "failed",
            "slurm_job_id": success["slurm_job_id"],
            "slurm_array_task_id": success["slurm_array_task_id"],
            "slurm_node": success["node"],
            "slurm_requested_cpus": success["requested_cpus"],
            "slurm_requested_memory": success["requested_memory"],
            "slurm_elapsed_seconds_internal": success["elapsed_seconds"],
        })
    return row


def stage_manifest(args: argparse.Namespace) -> None:
    base_columns, base = lookup_base(args.base_manifest, args.dataset_id)
    row = corrected_row(base, args.work_dir, stage=True)
    write_tsv(args.output, base_columns, [row])


def validate_current_rows(rows: Sequence[Mapping[str, str]]) -> None:
    if {row["dataset_id"] for row in rows} != EXPECTED_DATASETS:
        raise ValueError("corrected manifest does not contain exactly the three Tier 3A tuples")
    for row in rows:
        if row.get("tier3a_correction_status") != "current_production":
            raise ValueError(f"tuple is not current production: {row['dataset_id']}")
        if row.get("sweepga_origin_main_sha256") != SWEEPGA_SHA256 or row.get("sweepga_binary_sha256") != SWEEPGA_SHA256:
            raise ValueError(f"unapproved SweepGA binary selected: {row['dataset_id']}")
        if row.get("sweepga_origin_main_commit") != SWEEPGA_COMMIT:
            raise ValueError(f"unapproved SweepGA commit selected: {row['dataset_id']}")
        if REQUIRED_COMMAND_TOKEN not in row.get("sweepga_direct_command", ""):
            raise ValueError(f"native 1:1 syntax absent: {row['dataset_id']}")
        if row.get("sweepga_observed_query_multiplicity") != "1" or row.get("sweepga_observed_target_multiplicity") != "1":
            raise ValueError(f"native multiplicity audit failed: {row['dataset_id']}")
        if LEGACY_SWEEPGA_SHA256 in "\t".join(row.values()):
            raise ValueError(f"legacy SweepGA identity leaked into corrected row: {row['dataset_id']}")
        if int(row.get("impg_positive_variant_records", "0")) <= 0:
            raise ValueError(f"no positive IMPG variants: {row['dataset_id']}")
        if int(row.get("impg_callable_denominator_bases", "0")) <= 0:
            raise ValueError(f"no callable denominator: {row['dataset_id']}")
        if row.get("representative_variant_h1_h2_status") != "passed":
            raise ValueError(f"representative variant audit failed: {row['dataset_id']}")
        if row.get("slurm_state") != "COMPLETED":
            raise ValueError(f"corrected downstream Slurm job did not complete: {row['dataset_id']}")
        for field in (
            "sweepga_bounded_paf_path", "sweepga_native_multiplicity_audit_path",
            "impg_corrected_normalized_bcf_path", "impg_corrected_normalized_bcf_csi_path",
            "impg_corrected_normalized_vcf_gz_path", "impg_corrected_normalized_vcf_tbi_path",
        ):
            path = Path(row[field])
            if not path.is_file() or not path.stat().st_size:
                raise ValueError(f"missing corrected artifact {field}: {path}")


def finalize(args: argparse.Namespace) -> None:
    base_rows = read_tsv(args.base_manifest)
    base_by_id = {row["dataset_id"]: row for row in base_rows}
    if set(base_by_id) != EXPECTED_DATASETS:
        raise ValueError("base acquisition manifest tuple set changed")
    rows = [corrected_row(base_by_id[dataset], args.work_root / dataset, stage=False) for dataset in sorted(EXPECTED_DATASETS)]
    with args.sacct.open(encoding="utf-8", newline="") as handle:
        sacct_rows = list(csv.DictReader(handle, delimiter="|"))
    sacct_by_job = {
        item["JobID"]: item for item in sacct_rows
        if item.get("JobID") and not item["JobID"].endswith(".batch")
    }
    for row in rows:
        row["scheduler_sacct_path"] = str(args.sacct.resolve())
        row["sweepga_mapping_slurm_job_id"] = args.mapping_job_id
        row["sweepga_mapping_slurm_array_task_id"] = str(sorted(EXPECTED_DATASETS).index(row["dataset_id"]))
        mapping_job = sacct_by_job[f"{args.mapping_job_id}_{row['sweepga_mapping_slurm_array_task_id']}"]
        downstream_job = sacct_by_job[f"{row['slurm_job_id']}_{row['slurm_array_task_id']}"]
        row.update({
            "sweepga_mapping_slurm_state": mapping_job["State"],
            "sweepga_mapping_slurm_elapsed": mapping_job["Elapsed"],
            "slurm_node": downstream_job["NodeList"],
            "slurm_requested_cpus": downstream_job["AllocCPUS"],
            "slurm_requested_memory": downstream_job["ReqMem"],
            "slurm_state": downstream_job["State"],
            "slurm_elapsed": downstream_job["Elapsed"],
            "slurm_total_cpu": downstream_job["TotalCPU"],
            "slurm_max_rss": downstream_job["MaxRSS"],
        })
    columns = list(base_rows[0]) + [column for column in CORRECTION_COLUMNS if column not in base_rows[0]]
    validate_current_rows([{column: str(row.get(column, "")) for column in columns} for row in rows])
    prior_run_manifest = read_tsv(args.base_manifest.parent / "diploid_run_manifest.tsv")
    prior_run_by_id = {row["dataset_id"]: row for row in prior_run_manifest}
    supersession_rows = []
    for old, new in zip(sorted(base_rows, key=lambda item: item["dataset_id"]), rows):
        supersession_rows.extend([
            {
                "dataset_id": old["dataset_id"], "artifact_class": "SweepGA_production_mapping",
                "superseded_path": old["sweepga_bounded_paf_path"], "superseded_sha256": old["sweepga_bounded_paf_sha256"],
                "supersession_reason": "produced_without_accepted_reproducible_origin_main_binary_identity",
                "replacement_path": new["sweepga_bounded_paf_path"], "replacement_sha256": new["sweepga_bounded_paf_sha256"],
                "consumable": "no",
            },
            {
                "dataset_id": old["dataset_id"], "artifact_class": "mapping_derived_downstream_tree",
                "superseded_path": old.get("annotation_query_manifest_path", ""), "superseded_sha256": old.get("annotation_query_manifest_sha256", ""),
                "supersession_reason": "derived_from_superseded_mapping",
                "replacement_path": new["annotation_query_manifest_path"], "replacement_sha256": new["annotation_query_manifest_sha256"],
                "consumable": "no",
            },
            {
                "dataset_id": old["dataset_id"], "artifact_class": "mapping_derived_contig_map",
                "superseded_path": old.get("haplotype_contig_map_path", ""), "superseded_sha256": old.get("haplotype_contig_map_sha256", ""),
                "supersession_reason": "derived_from_superseded_mapping",
                "replacement_path": new["haplotype_contig_map_path"], "replacement_sha256": new["haplotype_contig_map_sha256"],
                "consumable": "no",
            },
            {
                "dataset_id": old["dataset_id"], "artifact_class": "IMPG_implicit_graph_index",
                "superseded_path": prior_run_by_id[old["dataset_id"]]["impg_index_path"], "superseded_sha256": "",
                "supersession_reason": "IMPG_index_consumed_superseded_mapping",
                "replacement_path": new["impg_corrected_index_path"], "replacement_sha256": "",
                "consumable": "no",
            },
            {
                "dataset_id": old["dataset_id"], "artifact_class": "IMPG_normalized_biological_BCF",
                "superseded_path": prior_run_by_id[old["dataset_id"]]["impg_normalized_bcf"], "superseded_sha256": prior_run_by_id[old["dataset_id"]]["impg_normalized_bcf_sha256"],
                "supersession_reason": "IMPG_calls_consumed_superseded_mapping",
                "replacement_path": new["impg_corrected_normalized_bcf_path"], "replacement_sha256": new["impg_corrected_normalized_bcf_sha256"],
                "consumable": "no",
            },
        ])
    write_tsv(args.output_manifest, columns, rows)
    write_tsv(args.supersession_ledger, (
        "dataset_id", "artifact_class", "superseded_path", "superseded_sha256",
        "supersession_reason", "replacement_path", "replacement_sha256", "consumable",
    ), supersession_rows)
    command_rows = []
    qc_lines = [
        "# Tier 3A origin/main SweepGA correction QC", "", "All three tuples passed the atomic publication gate.", "",
        f"Pinned unmodified origin/main commit: `{SWEEPGA_COMMIT}`.",
        f"Pinned byte-reproducible GNU Guix binary SHA-256: `{SWEEPGA_SHA256}`.",
        "Production used native `--num-mappings 1:1`; no mapping was replaced by post-hoc filtering.",
        "Higher-cap sensitivities were not run in this correction.", "", "## Tuple gates", "",
    ]
    for row in rows:
        command_rows.append({
            "dataset_id": row["dataset_id"], "binary_realpath": row["sweepga_origin_main_realpath"],
            "binary_sha256": row["sweepga_origin_main_sha256"], "command": row["sweepga_direct_command"],
            "impg_index_command": f"{row['impg_binary_path']} index -a {row['sweepga_bounded_paf_path']} -i {row['impg_corrected_index_path']} -t 8",
            "impg_partition_command": f"{row['impg_binary_path']} partition -a {row['sweepga_bounded_paf_path']} -i {row['impg_corrected_index_path']} -w 2000 -d 0 --min-missing-size 1 --min-boundary-distance 0 -o bed",
            "impg_query_command": f"{row['impg_binary_path']} query -a {row['sweepga_bounded_paf_path']} -i {row['impg_corrected_index_path']} -b {row['impg_corrected_focus_bed_path']} -d 0 -o vcf:poa --force-large-region --min-transitive-len 1",
            "slurm_job_id": row["slurm_job_id"], "slurm_array_task_id": row["slurm_array_task_id"],
        })
        summary = json.loads(Path(row["representative_variant_audit_path"]).read_text(encoding="utf-8"))
        direct = summary["direct_sequence_validation"]
        qc_lines.extend([
            f"### {row['scientific_name']} (`{row['dataset_id']}`)", "",
            f"Mapping records: {int(row['sweepga_mapping_records']):,}; query/target coverage: {float(row['sweepga_direct_query_coverage']):.4f}/{float(row['sweepga_direct_target_coverage']):.4f}; observed threshold multiplicity: {row['sweepga_observed_query_multiplicity']}:{row['sweepga_observed_target_multiplicity']}.",
            f"IMPG biological minimum numerator/denominator across reported annotation classes: {row['impg_positive_variant_records']}/{row['impg_callable_denominator_bases']}.",
            f"Representative allele: H1 {direct['h1_contig']}:{direct['h1_position_1based']} {direct['h1_base']} to H2 {direct['h2_contig']}:{direct['h2_position_1based']} {direct['h2_base']}; direct REF/ALT checks passed.",
            f"SweepGA mapping ran in Slurm array {row['sweepga_mapping_slurm_job_id']}_{row['sweepga_mapping_slurm_array_task_id']} for {row['sweepga_mapping_slurm_elapsed']} (the allocation later exited `{row['sweepga_mapping_slurm_state']}` at the superseded initial audit). Corrected downstream gate {row['slurm_job_id']}_{row['slurm_array_task_id']} completed on {row['slurm_node']} in {row['slurm_elapsed']} with {row['slurm_requested_cpus']} CPUs and {row['slurm_requested_memory']} requested; accounting reported TotalCPU `{row['slurm_total_cpu'] or 'unavailable'}` and MaxRSS `{row['slurm_max_rss'] or 'unavailable'}`. Full initial/retry sacct is recorded at `{args.sacct}`.", "",
        ])
    write_tsv(args.commands, (
        "dataset_id", "binary_realpath", "binary_sha256", "command",
        "impg_index_command", "impg_partition_command", "impg_query_command",
        "slurm_job_id", "slurm_array_task_id",
    ), command_rows)
    temporary_qc = args.qc.with_name(args.qc.name + ".tmp")
    temporary_qc.write_text("\n".join(qc_lines).rstrip() + "\n", encoding="utf-8")
    temporary_qc.replace(args.qc)
    validate_current_rows(read_tsv(args.output_manifest))


def validate_manifest(args: argparse.Namespace) -> None:
    validate_current_rows(read_tsv(args.manifest))


def parser() -> argparse.ArgumentParser:
    root = argparse.ArgumentParser(description=__doc__)
    commands = root.add_subparsers(dest="command", required=True)
    audit = commands.add_parser("audit-mapping")
    audit.add_argument("--paf", type=Path, required=True)
    audit.add_argument("--h1-fai", type=Path, required=True)
    audit.add_argument("--h2-fai", type=Path, required=True)
    audit.add_argument("--output", type=Path, required=True)
    audit.add_argument("--contig-map", type=Path, required=True)
    audit.add_argument("--native-recheck-paf", type=Path, required=True)
    audit.add_argument("--overlap-threshold", type=float, default=0.95)
    audit.set_defaults(function=audit_mapping)
    stage = commands.add_parser("stage-manifest")
    stage.add_argument("--base-manifest", type=Path, required=True)
    stage.add_argument("--dataset-id", required=True)
    stage.add_argument("--work-dir", type=Path, required=True)
    stage.add_argument("--output", type=Path, required=True)
    stage.set_defaults(function=stage_manifest)
    final = commands.add_parser("finalize")
    final.add_argument("--base-manifest", type=Path, required=True)
    final.add_argument("--work-root", type=Path, required=True)
    final.add_argument("--sacct", type=Path, required=True)
    final.add_argument("--mapping-job-id", required=True)
    final.add_argument("--output-manifest", type=Path, required=True)
    final.add_argument("--supersession-ledger", type=Path, required=True)
    final.add_argument("--commands", type=Path, required=True)
    final.add_argument("--qc", type=Path, required=True)
    final.set_defaults(function=finalize)
    validate = commands.add_parser("validate-manifest")
    validate.add_argument("--manifest", type=Path, required=True)
    validate.set_defaults(function=validate_manifest)
    return root


def main(argv: Sequence[str] | None = None) -> int:
    args = parser().parse_args(argv)
    args.function(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
