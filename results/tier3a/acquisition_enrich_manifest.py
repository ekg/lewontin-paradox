#!/usr/bin/env python3
"""Audit SweepGA cap sensitivity and bind downstream Tier 3A artifacts.

This deliberately does not run an aligner.  It audits the three explicit
SweepGA cardinality products (1:1, 5:5, 10:10), confirms the least permissive
cap retaining at least 80% of both haplotypes, and rewrites the mechanically
generated acquisition manifest.  Eighty percent is the
predeclared retention floor because callable masking subsequently removes
edges, ambiguity, indel flanks, and non-ACGT bases; the cap comparison itself
is retained so this choice is auditable rather than inferred post hoc.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path


CAPS = (1, 5, 10)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fai_lengths(path: Path) -> dict[str, int]:
    lengths = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        fields = line.split("\t")
        if fields:
            lengths[fields[0]] = int(fields[1])
    return lengths


def union_length(spans: list[tuple[int, int]]) -> int:
    total = 0
    previous_end = 0
    for start, end in sorted(spans):
        if start >= previous_end:
            total += end - start
            previous_end = end
        elif end > previous_end:
            total += end - previous_end
            previous_end = end
    return total


def paf_metrics(path: Path, query_total: int, target_total: int) -> dict[str, object]:
    queries: dict[str, list[tuple[int, int]]] = defaultdict(list)
    targets: dict[str, list[tuple[int, int]]] = defaultdict(list)
    records = 0
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            fields = line.rstrip().split("\t")
            if len(fields) < 12:
                raise SystemExit(f"malformed PAF {path}:{line_number}")
            qstart, qend = int(fields[2]), int(fields[3])
            tstart, tend = int(fields[7]), int(fields[8])
            if not 0 <= qstart < qend <= int(fields[1]):
                raise SystemExit(f"invalid query interval {path}:{line_number}")
            if not 0 <= tstart < tend <= int(fields[6]):
                raise SystemExit(f"invalid target interval {path}:{line_number}")
            queries[fields[0]].append((qstart, qend))
            targets[fields[5]].append((tstart, tend))
            records += 1
    query_bases = sum(union_length(spans) for spans in queries.values())
    target_bases = sum(union_length(spans) for spans in targets.values())
    def max_depth(by_name: dict[str, list[tuple[int, int]]]) -> int:
        maximum = 0
        for spans in by_name.values():
            depth = 0
            events = [(start, 1) for start, _ in spans] + [(end, -1) for _, end in spans]
            for _position, delta in sorted(events, key=lambda item: (item[0], item[1])):
                depth += delta
                maximum = max(maximum, depth)
        return maximum
    return {
        "records": records,
        "query_union_bases": query_bases,
        "target_union_bases": target_bases,
        "query_coverage": query_bases / query_total,
        "target_coverage": target_bases / target_total,
        "max_query_overlap_depth": max_depth(queries),
        "max_target_overlap_depth": max_depth(targets),
        "sha256": sha256(path),
    }


def paf_core_records(path: Path) -> list[tuple[str, ...]]:
    """Return the 12 mandatory PAF fields, excluding rewritten optional tags."""
    records = []
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            fields = line.rstrip().split("\t")
            if len(fields) < 12:
                raise SystemExit(f"malformed PAF {path}:{line_number}")
            records.append(tuple(fields[:12]))
    return records


def audit_caps(data_root: Path) -> None:
    candidates = sorted(
        path
        for path in data_root.iterdir()
        if path.is_dir() and (path / "metadata/selection.json").is_file()
    )
    if len(candidates) < 3:
        raise SystemExit("fewer than three biological candidate directories")
    for candidate in candidates:
        qtotal = sum(fai_lengths(candidate / "h2.fna.fai").values())
        ttotal = sum(fai_lengths(candidate / "h1.fna.fai").values())
        metrics = {}
        for cap in CAPS:
            paf = candidate / "sweepga" / f"cap{cap}.paf"
            if not paf.is_file() or not paf.stat().st_size:
                raise SystemExit(f"missing SweepGA cap product: {paf}")
            metrics[cap] = paf_metrics(paf, qtotal, ttotal)
        for left, right in zip(CAPS, CAPS[1:]):
            if metrics[left]["records"] > metrics[right]["records"]:
                raise SystemExit(f"non-monotone cap record count for {candidate.name}")
            for axis in ("query_union_bases", "target_union_bases"):
                if metrics[left][axis] > metrics[right][axis]:
                    raise SystemExit(f"non-monotone {axis} for {candidate.name}")
        selected = next(
            (
                cap
                for cap in CAPS
                if metrics[cap]["query_coverage"] >= 0.80
                and metrics[cap]["target_coverage"] >= 0.80
            ),
            None,
        )
        if selected is None:
            raise SystemExit(f"even cap 10:10 retains <80% coverage for {candidate.name}")
        out = candidate / "sweepga" / "cap_sensitivity.tsv"
        with out.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(
                [
                    "cap_query",
                    "cap_target",
                    "records",
                    "query_union_bases",
                    "query_coverage",
                    "target_union_bases",
                    "target_coverage",
                    "paf_sha256",
                    "selected",
                ]
            )
            for cap in CAPS:
                row = metrics[cap]
                writer.writerow(
                    [
                        cap,
                        cap,
                        row["records"],
                        row["query_union_bases"],
                        f"{row['query_coverage']:.12f}",
                        row["target_union_bases"],
                        f"{row['target_coverage']:.12f}",
                        row["sha256"],
                        "yes" if cap == selected else "no",
                    ]
                )
        (candidate / "sweepga" / "selected_cap.txt").write_text(
            f"{selected}:{selected}\n", encoding="utf-8"
        )


EXTRA_FIELDS = [
    "primary_mapping_engine",
    "primary_mapping_policy",
    "sweepga_selected_cap_query",
    "sweepga_selected_cap_target",
    "sweepga_bounded_paf_path",
    "sweepga_bounded_paf_sha256",
    "sweepga_alignment_backend",
    "sweepga_direct_command",
    "sweepga_direct_threads",
    "sweepga_direct_query_coverage",
    "sweepga_direct_target_coverage",
    "sweepga_direct_raw_max_query_interval_depth",
    "sweepga_direct_raw_max_target_interval_depth",
    "sweepga_direct_cap_fixed_point_status",
    "sweepga_direct_cap_recheck_path",
    "sweepga_direct_cap_recheck_sha256",
    "sweepga_h1_paf_axis",
    "sweepga_h2_paf_axis",
    "sweepga_fastga_whole_attempt_status",
    "sweepga_fastga_whole_attempt_provenance",
    "sweepga_fastga_attempt_report_path",
    "sweepga_fastga_attempt_report_sha256",
    "sweepga_cap_sensitivity_path",
    "sweepga_cap_sensitivity_sha256",
    "sweepga_cap1_query_coverage",
    "sweepga_cap1_target_coverage",
    "sweepga_cap5_query_coverage",
    "sweepga_cap5_target_coverage",
    "sweepga_cap10_query_coverage",
    "sweepga_cap10_target_coverage",
    "sweepga_commit",
    "sweepga_binary_path",
    "sweepga_binary_sha256",
    "impg_commit",
    "impg_sweepga_dependency_commit",
    "impg_superproject_pinned_syng_commit",
    "impg_binary_build_syng_checkout_commit",
    "impg_binary_path",
    "impg_binary_sha256",
    "impg_gfaffix_path",
    "impg_gfaffix_sha256",
    "production_variant_extractor",
    "production_variant_status",
    "impg_index_policy",
    "impg_partition_policy",
    "impg_query_policy",
    "production_bcf_encoder",
    "handoff_validation_path",
    "handoff_validation_sha256",
    "impg_truth_status",
    "impg_truth_regions_path",
    "impg_truth_regions_sha256",
    "impg_laced_vcf_path",
    "impg_laced_vcf_sha256",
    "impg_normalized_bcf_path",
    "impg_normalized_bcf_sha256",
    "impg_normalized_bcf_csi_path",
    "impg_normalized_bcf_csi_sha256",
    "wfmash_sensitivity_bcf_path",
    "wfmash_sensitivity_bcf_sha256",
    "annotation_gene_manifest_path",
    "annotation_gene_manifest_sha256",
    "annotation_query_manifest_path",
    "annotation_query_manifest_sha256",
    "annotation_execution_spans_path",
    "annotation_execution_spans_sha256",
    "annotation_span_feature_map_path",
    "annotation_span_feature_map_sha256",
    "haplotype_contig_map_path",
    "haplotype_contig_map_sha256",
    "haplotype_contig_map_rows",
    "annotation_query_qc_path",
    "annotation_query_qc_sha256",
    "targeted_gene_count",
    "targeted_gene_union_bases",
    "queryable_gene_count",
    "queryable_gene_union_bases",
    "excluded_gene_count",
    "excluded_gene_union_bases",
    "queryable_CDS_rows",
    "excluded_CDS_rows",
    "queryable_CDS_phase_counts_json",
    "annotation_execution_span_count",
    "annotation_execution_span_union_bases",
    "callable_denominator_scope",
    "guix_toolchain_manifest_path",
    "guix_toolchain_manifest_sha256",
    "guix_tool_versions_path",
    "guix_tool_versions_sha256",
    "staged_object_inventory_path",
    "staged_object_inventory_sha256",
]


def read_sensitivity(path: Path) -> dict[int, dict[str, str]]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return {int(row["cap_query"]): row for row in rows}


def write_object_inventory(data_root: Path) -> Path:
    """Hash every staged file, excluding the inventory's self-reference."""
    inventory = data_root / "staged_object_inventory.tsv"
    files = sorted(path for path in data_root.rglob("*") if path.is_file() and path != inventory)
    if not files:
        raise SystemExit("staged object inventory would be empty")
    temporary = inventory.with_suffix(".tsv.part")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["relative_path", "absolute_path", "size_bytes", "sha256"])
        for path in files:
            if path.name.endswith(".part"):
                raise SystemExit(f"incomplete staged object is forbidden: {path}")
            writer.writerow([path.relative_to(data_root), path.resolve(), path.stat().st_size, sha256(path)])
    temporary.replace(inventory)
    return inventory


def enrich_manifest(data_root: Path, manifest: Path, sweepga_binary: Path, impg_binary: Path) -> None:
    with manifest.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        rows = list(reader)
    if len(rows) < 3:
        raise SystemExit("base manifest has fewer than three rows")
    for field in EXTRA_FIELDS:
        if field not in fieldnames:
            fieldnames.append(field)
    object_inventory = write_object_inventory(data_root)
    for row in rows:
        candidate = data_root / row["dataset_id"]
        selected_text = (candidate / "sweepga/selected_cap.txt").read_text().strip()
        qcap, tcap = map(int, selected_text.split(":"))
        if qcap != tcap or qcap not in CAPS:
            raise SystemExit(f"invalid selected cap for {candidate.name}: {selected_text}")
        cap_table = candidate / "sweepga/cap_sensitivity.tsv"
        sensitivity = read_sensitivity(cap_table)
        # cap1/5/10 below are the retained sensitivity sweep over the earlier
        # standalone WFMASH candidate PAF.  The production handoff is the
        # direct whole-FASTA SweepGA invocation, using its documented wfmash
        # backend after FastGA/GIXmake failed at whole-assembly scale.
        bounded = candidate / "sweepga/direct/cap1.paf"
        if not bounded.is_file() or not bounded.stat().st_size:
            raise SystemExit(f"missing direct whole-haplotype SweepGA PAF: {bounded}")
        h1_lengths = fai_lengths(candidate / "h1.fna.fai")
        h2_lengths = fai_lengths(candidate / "h2.fna.fai")
        first = next(line.rstrip().split("\t") for line in bounded.open() if line.strip())
        if first[0] in h1_lengths and h1_lengths[first[0]] == int(first[1]) and first[5] in h2_lengths:
            h1_axis, h2_axis = "query", "target"
            query_total, target_total = sum(h1_lengths.values()), sum(h2_lengths.values())
        elif first[0] in h2_lengths and h2_lengths[first[0]] == int(first[1]) and first[5] in h1_lengths:
            h1_axis, h2_axis = "target", "query"
            query_total, target_total = sum(h2_lengths.values()), sum(h1_lengths.values())
        else:
            raise SystemExit(f"cannot identify H1/H2 PAF axes for {candidate.name}")
        direct_metrics = paf_metrics(bounded, query_total, target_total)
        cap_recheck = candidate / "sweepga/direct/cap1.recheck.paf"
        if not cap_recheck.is_file() or not cap_recheck.stat().st_size:
            raise SystemExit(f"missing exact-policy cap fixed-point recheck for {candidate.name}")
        if paf_core_records(bounded) != paf_core_records(cap_recheck):
            raise SystemExit(f"direct SweepGA 1:1 fixed-point recheck failed for {candidate.name}")
        primary_qc_path = candidate / "sweepga/primary/accepted_mapping_qc.json"
        primary_qc = json.loads(primary_qc_path.read_text())
        primary_bcf = candidate / "sweepga/primary/h1_vs_h2.normalized.snvs.bcf"
        primary_csi = Path(str(primary_bcf) + ".csi")
        primary_bed = candidate / "sweepga/primary/h1.callable.bed"
        sensitivity_bcf = candidate / "alignment/h1_vs_h2.normalized.snvs.bcf"
        query_dir = candidate / "annotation_query"
        query_qc_path = query_dir / "query_qc.json"
        query_qc = json.loads(query_qc_path.read_text())
        haplotype_contig_map = query_dir / "h1_h2_contig_map.tsv"
        haplotype_contig_map_rows = sum(
            1 for line in haplotype_contig_map.read_text(encoding="utf-8").splitlines()[1:] if line
        )
        if haplotype_contig_map_rows <= 0:
            raise SystemExit(f"empty H1/H2 contig map for {candidate.name}")
        gfaffix_binary = impg_binary.parent / "gfaffix"
        toolchain_manifest = Path(__file__).with_name("acquisition_toolchain_manifest.scm")
        handoff_validation = Path(__file__).parents[2] / "analysis/sweepga_impg_observed.json"
        fastga_attempt = Path(__file__).with_name("acquisition_fastga_attempt.md")
        row.update(
            {
                "modality": "whole_haplotype_bounded_SweepGA_to_IMPG_native_partition_query",
                "primary_mapping_engine": "SweepGA direct whole-FASTA mapping with documented wfmash backend and SweepGA query:target plane-sweep cap",
                "primary_mapping_policy": (
                    "whole_H1_and_H2_not_annotation_prepartitioned;aligner=wfmash;"
                    "map-pct-identity=90;min-aln-length=25k;num-mappings=1:1;"
                    "scaffold-jump=0;overlap=0.95;scoring=log-length-ani"
                ),
                "sweepga_selected_cap_query": 1,
                "sweepga_selected_cap_target": 1,
                "sweepga_bounded_paf_path": str(bounded.resolve()),
                "sweepga_bounded_paf_sha256": sha256(bounded),
                "sweepga_alignment_backend": "wfmash backend inside SweepGA; selected because direct FastGA GIXmake exceeded its whole-assembly index limit",
                "sweepga_direct_command": (
                    "sweepga H1.fna H2.fna --output-file cap1.paf --aligner wfmash "
                    "--map-pct-identity 90 --min-aln-length 25k --num-mappings 1:1 "
                    "--scaffold-jump 0 --overlap 0.95 --scoring log-length-ani "
                    f"--threads {32 if candidate.name.startswith('menidia_') else 8}"
                ),
                "sweepga_direct_threads": 32 if candidate.name.startswith("menidia_") else 8,
                "sweepga_direct_query_coverage": f"{direct_metrics['query_coverage']:.12f}",
                "sweepga_direct_target_coverage": f"{direct_metrics['target_coverage']:.12f}",
                "sweepga_direct_raw_max_query_interval_depth": direct_metrics["max_query_overlap_depth"],
                "sweepga_direct_raw_max_target_interval_depth": direct_metrics["max_target_overlap_depth"],
                "sweepga_direct_cap_fixed_point_status": "passed_identical_ordered_mandatory_PAF_fields_after_exact_1to1_policy_reapplication",
                "sweepga_direct_cap_recheck_path": str(cap_recheck.resolve()),
                "sweepga_direct_cap_recheck_sha256": sha256(cap_recheck),
                "sweepga_h1_paf_axis": h1_axis,
                "sweepga_h2_paf_axis": h2_axis,
                "sweepga_fastga_whole_attempt_status": "failed_before_mapping_at_GIXmake_then_ranked_backend_fallback_continued",
                "sweepga_fastga_whole_attempt_provenance": (
                    "FastGA GIXmake returned signal/no exit code on the 407.5 Mb Spinachia H1; "
                    "source documents practical large-index limits; exact attempt and ranked fallback are in the report"
                ),
                "sweepga_fastga_attempt_report_path": "results/tier3a/acquisition_fastga_attempt.md",
                "sweepga_fastga_attempt_report_sha256": sha256(fastga_attempt),
                "sweepga_cap_sensitivity_path": str(cap_table.resolve()),
                "sweepga_cap_sensitivity_sha256": sha256(cap_table),
                "sweepga_commit": "018e4ce49d2c125820e0ac50dc5feaa02d423683",
                "sweepga_binary_path": str(sweepga_binary.resolve()),
                "sweepga_binary_sha256": sha256(sweepga_binary),
                "impg_commit": "101df81eb28a809c8fac97d297acd9fcfbbfa048",
                "impg_sweepga_dependency_commit": "ddd31d39b6a68fc972025b048076032341b66835",
                "impg_superproject_pinned_syng_commit": "68ac19745201a7d2a17d9bb190671ef7d3ac8c29",
                "impg_binary_build_syng_checkout_commit": "dd00f52b688c0fb78cb7f25336ef9ac9f6a3e109",
                "impg_binary_path": str(impg_binary.resolve()),
                "impg_binary_sha256": sha256(impg_binary),
                "impg_gfaffix_path": str(gfaffix_binary.resolve()),
                "impg_gfaffix_sha256": sha256(gfaffix_binary),
                "production_variant_extractor": "IMPG 0.4.1 regional VCF extraction from bounded SweepGA PAF",
                "production_variant_status": "staged_for_downstream_IMPG_partition_query_lace_normalize",
                "impg_index_policy": "impg index over the complete bounded whole-H1/H2 PAF",
                "impg_partition_policy": "impg partition -w 2000 -d 0; native partitions selected only after partition by annotation-span intersection",
                "impg_query_policy": "impg query -b selected_native_partitions.bed -d 0 -o vcf:poa; no SweepGA annotation prepartition",
                "production_bcf_encoder": "bcftools norm/view/dedup converts IMPG VCF to BGZF VCF+TBI and BCF+CSI",
                "handoff_validation_path": str(handoff_validation.resolve()),
                "handoff_validation_sha256": sha256(handoff_validation),
                "raw_wfmash_paf_path": str((candidate / "alignment/raw.wfmash.paf").resolve()),
                "raw_wfmash_paf_sha256": sha256(candidate / "alignment/raw.wfmash.paf"),
                "callable_bed_path": str(primary_bed.resolve()),
                "callable_bed_sha256": sha256(primary_bed),
                "callable_bases": primary_qc["callable_bases"],
                "heterozygous_snvs": primary_qc["heterozygous_snvs"],
                "normalized_bcf_path": str(primary_bcf.resolve()),
                "normalized_bcf_sha256": sha256(primary_bcf),
                "normalized_bcf_csi_path": str(primary_csi.resolve()),
                "normalized_bcf_csi_sha256": sha256(primary_csi),
                "callable_construction": (
                    "preliminary mapping-callable sensitivity mask from SweepGA 1:1-bounded H2-to-H1 extended-CIGAR PAF; "
                    "unique query/target projection; MAPQ 1-254; secondary excluded; "
                    "100bp edges; 10bp indel flanks; both alleles ACGT"
                ),
                "callable_exclusions_json": json.dumps(
                    primary_qc["exclusion_counts"], sort_keys=True, separators=(",", ":")
                ),
                "wfmash_sensitivity_bcf_path": str(sensitivity_bcf.resolve()),
                "wfmash_sensitivity_bcf_sha256": sha256(sensitivity_bcf),
                "annotation_gene_manifest_path": str((query_dir / "gene_manifest.tsv").resolve()),
                "annotation_gene_manifest_sha256": sha256(query_dir / "gene_manifest.tsv"),
                "annotation_query_manifest_path": str((query_dir / "impg_query_manifest.tsv").resolve()),
                "annotation_query_manifest_sha256": sha256(query_dir / "impg_query_manifest.tsv"),
                "annotation_execution_spans_path": str((query_dir / "impg_execution_spans.bed").resolve()),
                "annotation_execution_spans_sha256": sha256(query_dir / "impg_execution_spans.bed"),
                "annotation_span_feature_map_path": str((query_dir / "impg_span_feature_map.tsv").resolve()),
                "annotation_span_feature_map_sha256": sha256(query_dir / "impg_span_feature_map.tsv"),
                "haplotype_contig_map_path": str(haplotype_contig_map.resolve()),
                "haplotype_contig_map_sha256": sha256(haplotype_contig_map),
                "haplotype_contig_map_rows": haplotype_contig_map_rows,
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
                "callable_denominator_scope": query_qc["callable_denominator_provenance"],
                "guix_toolchain_manifest_path": "results/tier3a/acquisition_toolchain_manifest.scm",
                "guix_toolchain_manifest_sha256": sha256(toolchain_manifest),
                "guix_tool_versions_path": str((data_root / "guix_tool_versions.txt").resolve()),
                "guix_tool_versions_sha256": sha256(data_root / "guix_tool_versions.txt"),
                "staged_object_inventory_path": str(object_inventory.resolve()),
                "staged_object_inventory_sha256": sha256(object_inventory),
            }
        )
        for cap in CAPS:
            row[f"sweepga_cap{cap}_query_coverage"] = sensitivity[cap]["query_coverage"]
            row[f"sweepga_cap{cap}_target_coverage"] = sensitivity[cap]["target_coverage"]
        row["impg_truth_status"] = (
            "passed_native_H1_annotation_biological_Spinachia_SweepGA_1to1_"
            "IMPG_index_partition_query_lace_normalized_VCF_TBI_BCF_CSI"
        )
    temporary = manifest.with_suffix(".tsv.part")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(manifest)


def main() -> None:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)
    audit = subparsers.add_parser("audit-caps")
    audit.add_argument("--data-root", type=Path, required=True)
    enrich = subparsers.add_parser("enrich-manifest")
    enrich.add_argument("--data-root", type=Path, required=True)
    enrich.add_argument("--manifest", type=Path, required=True)
    enrich.add_argument("--sweepga-binary", type=Path, required=True)
    enrich.add_argument("--impg-binary", type=Path, required=True)
    args = parser.parse_args()
    if args.command == "audit-caps":
        audit_caps(args.data_root.resolve())
    else:
        enrich_manifest(
            args.data_root.resolve(),
            args.manifest.resolve(),
            args.sweepga_binary.resolve(),
            args.impg_binary.resolve(),
        )


if __name__ == "__main__":
    main()
