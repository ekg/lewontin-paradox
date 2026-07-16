#!/usr/bin/env python3
"""Compile exact-H1 GFF3 coding features into deterministic IMPG query inputs.

This program does not partition a mapping.  It derives annotation targets from
the provider's native GFF3, verifies CDS phase chains, and screens whole gene
spans against the H1-axis union of an already bounded whole-H1/H2 SweepGA PAF
(query or target is detected and recorded). The emitted BED contains only
overlap/touch-merged annotation spans. A
consumer must intersect those spans with *IMPG's own* ``partition`` BED before
calling ``impg query``.
"""

from __future__ import annotations

import argparse
import bisect
import collections
import csv
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import unquote


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_attributes(text: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for item in text.split(";"):
        if not item:
            continue
        key, separator, value = item.partition("=")
        if separator:
            result[unquote(key)] = unquote(value)
    return result


def read_fai(path: Path) -> dict[str, int]:
    result: dict[str, int] = {}
    with path.open(encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 2 or fields[0] in result:
                raise SystemExit(f"invalid FAI row {line_number}: {path}")
            result[fields[0]] = int(fields[1])
    if not result:
        raise SystemExit(f"empty FAI: {path}")
    return result


def merge(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    result: list[list[int]] = []
    for start, end in sorted(intervals):
        if start >= end:
            continue
        if result and start <= result[-1][1]:
            result[-1][1] = max(result[-1][1], end)
        else:
            result.append([start, end])
    return [(start, end) for start, end in result]


def union_length(by_contig: dict[str, list[tuple[int, int]]]) -> int:
    return sum(end - start for values in by_contig.values() for start, end in merge(values))


def read_h1_coverage(path: Path, fai: dict[str, int]) -> tuple[dict[str, list[tuple[int, int]]], str]:
    intervals: dict[str, list[tuple[int, int]]] = collections.defaultdict(list)
    observed_axis: str | None = None
    with path.open(encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 12:
                raise SystemExit(f"PAF row {line_number} has fewer than 12 columns")
            query_match = fields[0] in fai and fai[fields[0]] == int(fields[1])
            target_match = fields[5] in fai and fai[fields[5]] == int(fields[6])
            if query_match == target_match:
                raise SystemExit(f"PAF row {line_number} does not identify exactly one H1 axis")
            axis = "query" if query_match else "target"
            if observed_axis is not None and axis != observed_axis:
                raise SystemExit("H1 switches PAF axes between mapping records")
            observed_axis = axis
            if axis == "query":
                contig, declared_length = fields[0], int(fields[1])
                start, end = int(fields[2]), int(fields[3])
            else:
                contig, declared_length = fields[5], int(fields[6])
                start, end = int(fields[7]), int(fields[8])
            if not 0 <= start < end <= declared_length:
                raise SystemExit(f"invalid PAF H1 interval at row {line_number}")
            intervals[contig].append((start, end))
    if not intervals:
        raise SystemExit(f"bounded SweepGA PAF is empty: {path}")
    assert observed_axis is not None
    return {contig: merge(values) for contig, values in intervals.items()}, observed_axis


def write_haplotype_contig_map(
    paf: Path, h1_fai: dict[str, int], h2_fai: dict[str, int], output: Path
) -> None:
    pairs: dict[tuple[str, str], dict[str, object]] = {}
    mapped_h1: set[str] = set()
    mapped_h2: set[str] = set()
    with paf.open(encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            qname, qlen, qstart, qend = fields[0], int(fields[1]), int(fields[2]), int(fields[3])
            tname, tlen, tstart, tend = fields[5], int(fields[6]), int(fields[7]), int(fields[8])
            if qname in h1_fai and tname in h2_fai:
                if h1_fai[qname] != qlen or h2_fai[tname] != tlen:
                    raise SystemExit(f"PAF H1/H2 length mismatch at row {line_number}")
                h1_name, h1_interval, h2_name, h2_interval = qname, (qstart, qend), tname, (tstart, tend)
            elif qname in h2_fai and tname in h1_fai:
                if h2_fai[qname] != qlen or h1_fai[tname] != tlen:
                    raise SystemExit(f"PAF H1/H2 length mismatch at row {line_number}")
                h1_name, h1_interval, h2_name, h2_interval = tname, (tstart, tend), qname, (qstart, qend)
            else:
                raise SystemExit(f"PAF H1/H2 dictionary mismatch at row {line_number}")
            item = pairs.setdefault((h1_name, h2_name), {"h1": [], "h2": [], "strands": set(), "records": 0})
            item["h1"].append(h1_interval)
            item["h2"].append(h2_interval)
            item["strands"].add(fields[4])
            item["records"] += 1
            mapped_h1.add(h1_name)
            mapped_h2.add(h2_name)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "relation_id", "h1_contig", "h1_length", "h2_contig", "h2_length", "strands",
            "mapping_records", "h1_aligned_union_bases", "h2_aligned_union_bases", "status",
        ])
        relation = 0
        for (h1_name, h2_name), item in sorted(pairs.items()):
            relation += 1
            writer.writerow([
                f"contig_relation_{relation:06d}", h1_name, h1_fai[h1_name], h2_name,
                h2_fai[h2_name], ",".join(sorted(item["strands"])), item["records"],
                sum(end - start for start, end in merge(item["h1"])),
                sum(end - start for start, end in merge(item["h2"])), "selected_bounded_mapping",
            ])
        for h1_name in sorted(set(h1_fai) - mapped_h1):
            relation += 1
            writer.writerow([f"contig_relation_{relation:06d}", h1_name, h1_fai[h1_name], "", "", "", 0, 0, 0, "H1_no_selected_mapping"])
        for h2_name in sorted(set(h2_fai) - mapped_h2):
            relation += 1
            writer.writerow([f"contig_relation_{relation:06d}", "", "", h2_name, h2_fai[h2_name], "", 0, 0, 0, "H2_no_selected_mapping"])


def fully_covered(intervals: list[tuple[int, int]], start: int, end: int) -> bool:
    starts = [value[0] for value in intervals]
    index = bisect.bisect_right(starts, start) - 1
    return index >= 0 and intervals[index][0] <= start and intervals[index][1] >= end


@dataclass(frozen=True)
class Gene:
    line: int
    seqid: str
    start: int
    end: int
    strand: str
    gene_id: str
    gene_name: str
    locus_tag: str
    biotype: str
    pseudo: bool
    coordinate_valid: bool


@dataclass(frozen=True)
class Transcript:
    line: int
    transcript_id: str
    parent_gene_ids: tuple[str, ...]
    seqid: str
    start: int
    end: int
    strand: str
    pseudo: bool


@dataclass(frozen=True)
class CDS:
    line: int
    seqid: str
    start: int
    end: int
    strand: str
    phase: str
    cds_id: str
    parent_transcript_ids: tuple[str, ...]
    protein_id: str
    locus_tag: str


def load_annotation(gff: Path, fai: dict[str, int]) -> tuple[dict[str, Gene], dict[str, Transcript], list[CDS]]:
    genes: dict[str, Gene] = {}
    transcripts: dict[str, Transcript] = {}
    cds_rows: list[CDS] = []
    with gff.open(encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise SystemExit(f"GFF3 row {line_number} does not have 9 columns")
            seqid, _source, kind, start_text, end_text, _score, strand, phase, attributes = fields
            if kind not in {"gene", "mRNA", "transcript", "CDS"}:
                continue
            try:
                start, end = int(start_text), int(end_text)
            except ValueError as error:
                raise SystemExit(f"noninteger GFF3 coordinates at row {line_number}") from error
            attrs = parse_attributes(attributes)
            start0 = start - 1
            parents = tuple(value for value in attrs.get("Parent", "").split(",") if value)
            pseudo = attrs.get("pseudo", "").lower() == "true" or "pseudo" in attrs
            if kind == "gene":
                gene_id = attrs.get("ID", "")
                if not gene_id or gene_id in genes:
                    raise SystemExit(f"missing/duplicate gene ID at GFF3 row {line_number}")
                coordinate_valid = (
                    seqid in fai and 0 <= start0 < end <= fai[seqid] and strand in {"+", "-"}
                )
                genes[gene_id] = Gene(
                    line_number, seqid, start0, end, strand, gene_id,
                    attrs.get("gene", attrs.get("Name", gene_id)), attrs.get("locus_tag", ""),
                    attrs.get("gene_biotype", ""), pseudo, coordinate_valid,
                )
            elif kind in {"mRNA", "transcript"}:
                transcript_id = attrs.get("ID", "")
                if transcript_id and transcript_id not in transcripts:
                    transcripts[transcript_id] = Transcript(
                        line_number, transcript_id, parents, seqid, start0, end, strand, pseudo
                    )
            else:
                cds_rows.append(
                    CDS(
                        line_number, seqid, start0, end, strand, phase, attrs.get("ID", ""),
                        parents, attrs.get("protein_id", attrs.get("Name", "")),
                        attrs.get("locus_tag", ""),
                    )
                )
    if not genes or not transcripts or not cds_rows:
        raise SystemExit("GFF3 lacks genes, transcripts, or CDS features")
    return genes, transcripts, cds_rows


def transcript_phase_qc(
    transcript: Transcript, rows: list[CDS], fai: dict[str, int]
) -> tuple[bool, str]:
    if not rows:
        return False, "no_CDS"
    if transcript.pseudo:
        return False, "pseudo_transcript"
    if transcript.strand not in {"+", "-"}:
        return False, "invalid_transcript_strand"
    ordered = sorted(rows, key=lambda row: (row.start, row.end, row.line), reverse=transcript.strand == "-")
    coding_modulo: int | None = None
    for index, row in enumerate(ordered):
        if row.seqid != transcript.seqid or row.strand != transcript.strand:
            return False, "CDS_transcript_contig_or_strand_mismatch"
        if row.seqid not in fai or not 0 <= row.start < row.end <= fai[row.seqid]:
            return False, "invalid_CDS_coordinate"
        if not row.cds_id:
            return False, "missing_CDS_ID"
        if row.phase not in {"0", "1", "2"}:
            return False, "invalid_CDS_phase"
        phase = int(row.phase)
        if index == 0:
            coding_modulo = (row.end - row.start - phase) % 3
        else:
            assert coding_modulo is not None
            expected = (3 - coding_modulo) % 3
            if phase != expected:
                return False, f"phase_chain_mismatch_expected_{expected}_observed_{phase}"
            coding_modulo = (coding_modulo + row.end - row.start) % 3
    return True, "passed"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset-id", required=True)
    parser.add_argument("--gff", type=Path, required=True)
    parser.add_argument("--fai", type=Path, required=True)
    parser.add_argument("--h2-fai", type=Path, required=True)
    parser.add_argument("--bounded-paf", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--annotation-release", required=True)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    fai = read_fai(args.fai)
    h2_fai = read_fai(args.h2_fai)
    coverage, h1_paf_axis = read_h1_coverage(args.bounded_paf, fai)
    genes, transcripts, cds_rows = load_annotation(args.gff, fai)

    cds_by_transcript: dict[str, list[CDS]] = collections.defaultdict(list)
    orphan_rows: set[int] = set()
    for row in cds_rows:
        if not row.parent_transcript_ids:
            orphan_rows.add(row.line)
        for transcript_id in row.parent_transcript_ids:
            cds_by_transcript[transcript_id].append(row)

    transcript_qc: dict[str, tuple[bool, str]] = {}
    for transcript_id, transcript in transcripts.items():
        transcript_qc[transcript_id] = transcript_phase_qc(
            transcript, cds_by_transcript.get(transcript_id, []), fai
        )

    gene_transcripts: dict[str, list[str]] = collections.defaultdict(list)
    for transcript_id, transcript in transcripts.items():
        for gene_id in transcript.parent_gene_ids:
            gene_transcripts[gene_id].append(transcript_id)

    gene_status: dict[str, tuple[bool, bool, str]] = {}
    for gene_id, gene in genes.items():
        transcript_ids = sorted(set(gene_transcripts.get(gene_id, [])))
        has_cds = any(cds_by_transcript.get(item) for item in transcript_ids)
        targeted = gene.biotype == "protein_coding" and has_cds
        reasons: list[str] = []
        if not targeted:
            reasons.append("not_protein_coding_or_no_CDS")
        if gene.pseudo:
            reasons.append("pseudo_gene")
        if not gene.coordinate_valid:
            reasons.append("invalid_or_non_H1_gene_coordinate")
        valid_transcripts = [item for item in transcript_ids if transcript_qc.get(item) == (True, "passed")]
        if targeted and not valid_transcripts:
            reasons.append("no_phase_valid_CDS_transcript")
        mapped = gene.coordinate_valid and fully_covered(
            coverage.get(gene.seqid, []), gene.start, gene.end
        )
        if targeted and gene.coordinate_valid and not mapped:
            reasons.append("not_fully_covered_by_selected_bounded_mapping")
        queryable = targeted and not gene.pseudo and gene.coordinate_valid and bool(valid_transcripts) and mapped
        gene_status[gene_id] = (targeted, queryable, "passed" if queryable else ";".join(reasons))

    queryable_genes = [genes[item] for item, status in gene_status.items() if status[1]]
    spans_by_contig: dict[str, list[tuple[int, int]]] = collections.defaultdict(list)
    for gene in queryable_genes:
        spans_by_contig[gene.seqid].append((gene.start, gene.end))

    span_records: list[tuple[str, int, int, str]] = []
    for contig in fai:
        for number, (start, end) in enumerate(merge(spans_by_contig.get(contig, [])), 1):
            span_records.append((contig, start, end, f"{args.dataset_id}.target.{contig}.{number:06d}"))
    if not span_records:
        raise SystemExit("annotation query manifest has no queryable execution spans")

    spans_for_lookup: dict[str, list[tuple[int, int, str]]] = collections.defaultdict(list)
    for record in span_records:
        spans_for_lookup[record[0]].append((record[1], record[2], record[3]))

    gene_to_span: dict[str, str] = {}
    span_gene_ids: dict[str, list[str]] = collections.defaultdict(list)
    for gene in queryable_genes:
        candidates = spans_for_lookup[gene.seqid]
        starts = [value[0] for value in candidates]
        index = bisect.bisect_right(starts, gene.start) - 1
        if index < 0 or candidates[index][1] < gene.end:
            raise SystemExit(f"queryable gene lacks merged execution span: {gene.gene_id}")
        span_id = candidates[index][2]
        gene_to_span[gene.gene_id] = span_id
        span_gene_ids[span_id].append(gene.gene_id)

    gene_manifest = args.output_dir / "gene_manifest.tsv"
    with gene_manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "dataset_id", "annotation_release", "gene_id", "gene_name", "locus_tag", "biotype",
            "contig", "start_1based", "end_1based", "start_0based", "end_0based_exclusive",
            "strand", "targeted", "queryable", "execution_span_id", "exclusion_reason",
            "transcript_ids", "source_gff_line",
        ])
        for gene_id, gene in sorted(genes.items(), key=lambda item: (item[1].seqid, item[1].start, item[0])):
            targeted, queryable, reason = gene_status[gene_id]
            writer.writerow([
                args.dataset_id, args.annotation_release, gene_id, gene.gene_name, gene.locus_tag,
                gene.biotype, gene.seqid, gene.start + 1, gene.end, gene.start, gene.end, gene.strand,
                "yes" if targeted else "no", "yes" if queryable else "no",
                gene_to_span[gene_id] if queryable else "", "" if queryable else reason,
                ",".join(sorted(set(gene_transcripts.get(gene_id, [])))), gene.line,
            ])

    query_manifest = args.output_dir / "impg_query_manifest.tsv"
    queryable_cds_rows = 0
    excluded_cds_rows = 0
    phase_counts: collections.Counter[str] = collections.Counter()
    queryable_feature_ids: dict[str, list[str]] = collections.defaultdict(list)
    with query_manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "feature_row_id", "dataset_id", "annotation_release", "gene_ids", "transcript_ids",
            "cds_id", "protein_id", "locus_tag", "contig", "start_1based", "end_1based",
            "start_0based", "end_0based_exclusive", "strand", "phase", "targeted", "queryable",
            "execution_span_ids", "exclusion_reason", "source_gff_line",
        ])
        for row in sorted(cds_rows, key=lambda item: item.line):
            parent_ids = row.parent_transcript_ids
            resolved_gene_ids = sorted({
                gene_id
                for transcript_id in parent_ids
                for gene_id in transcripts.get(transcript_id, Transcript(0, "", (), "", 0, 0, "", False)).parent_gene_ids
                if gene_id in genes
            })
            reasons: list[str] = []
            if not parent_ids:
                reasons.append("missing_transcript_parent")
            missing_transcripts = [item for item in parent_ids if item not in transcripts]
            if missing_transcripts:
                reasons.append("unresolved_transcript_parent")
            invalid_transcripts = [
                transcript_qc[item][1] for item in parent_ids
                if item in transcript_qc and not transcript_qc[item][0]
            ]
            if invalid_transcripts:
                reasons.extend(invalid_transcripts)
            if not resolved_gene_ids:
                reasons.append("unresolved_gene_parent")
            selected_genes = [gene_id for gene_id in resolved_gene_ids if gene_status[gene_id][1]]
            targeted = any(gene_status[gene_id][0] for gene_id in resolved_gene_ids)
            queryable = bool(selected_genes) and all(
                transcript_qc.get(item, (False, ""))[0] for item in parent_ids
            )
            if targeted and not selected_genes:
                reasons.extend(gene_status[gene_id][2] for gene_id in resolved_gene_ids if not gene_status[gene_id][1])
            span_ids = sorted({gene_to_span[item] for item in selected_genes})
            if queryable and not span_ids:
                raise SystemExit(f"queryable CDS at GFF row {row.line} lacks execution span")
            feature_row_id = f"{args.dataset_id}.gff.CDS.{row.line:09d}"
            if queryable:
                queryable_cds_rows += 1
                phase_counts[row.phase] += 1
                for span_id in span_ids:
                    queryable_feature_ids[span_id].append(feature_row_id)
            else:
                excluded_cds_rows += 1
            writer.writerow([
                feature_row_id, args.dataset_id, args.annotation_release, ",".join(resolved_gene_ids),
                ",".join(parent_ids), row.cds_id, row.protein_id, row.locus_tag, row.seqid,
                row.start + 1, row.end, row.start, row.end, row.strand, row.phase,
                "yes" if targeted else "no", "yes" if queryable else "no", ",".join(span_ids),
                "" if queryable else ";".join(sorted(set(reasons))), row.line,
            ])
    if queryable_cds_rows == 0 or set(phase_counts) != {"0", "1", "2"}:
        raise SystemExit("query manifest lacks queryable CDS rows or all three verified phases")

    execution_bed = args.output_dir / "impg_execution_spans.bed"
    with execution_bed.open("w", encoding="utf-8") as handle:
        for contig, start, end, span_id in span_records:
            handle.write(f"{contig}\t{start}\t{end}\t{span_id}\n")

    feature_map = args.output_dir / "impg_span_feature_map.tsv"
    with feature_map.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "execution_span_id", "contig", "start_0based", "end_0based_exclusive",
            "gene_ids", "transcript_ids", "cds_ids", "feature_row_ids",
        ])
        for contig, start, end, span_id in span_records:
            member_gene_ids = sorted(span_gene_ids[span_id])
            member_transcripts = sorted({item for gene_id in member_gene_ids for item in gene_transcripts[gene_id] if transcript_qc.get(item, (False, ""))[0]})
            member_cds = sorted({row.cds_id for item in member_transcripts for row in cds_by_transcript[item]})
            writer.writerow([
                span_id, contig, start, end, ",".join(member_gene_ids), ",".join(member_transcripts),
                ",".join(member_cds), ",".join(queryable_feature_ids[span_id]),
            ])

    haplotype_contig_map = args.output_dir / "h1_h2_contig_map.tsv"
    write_haplotype_contig_map(args.bounded_paf, fai, h2_fai, haplotype_contig_map)

    targeted_genes = [genes[item] for item, status in gene_status.items() if status[0]]
    excluded_genes = [genes[item] for item, status in gene_status.items() if status[0] and not status[1]]
    targeted_intervals: dict[str, list[tuple[int, int]]] = collections.defaultdict(list)
    queryable_intervals: dict[str, list[tuple[int, int]]] = collections.defaultdict(list)
    excluded_intervals: dict[str, list[tuple[int, int]]] = collections.defaultdict(list)
    reason_counts: collections.Counter[str] = collections.Counter()
    for gene in targeted_genes:
        if gene.coordinate_valid:
            targeted_intervals[gene.seqid].append((gene.start, gene.end))
    for gene in queryable_genes:
        queryable_intervals[gene.seqid].append((gene.start, gene.end))
    for gene in excluded_genes:
        if gene.coordinate_valid:
            excluded_intervals[gene.seqid].append((gene.start, gene.end))
        for reason in gene_status[gene.gene_id][2].split(";"):
            reason_counts[reason] += 1

    query_qc = args.output_dir / "query_qc.json"
    qc = {
        "schema_version": "tier3a-annotation-query-v1",
        "dataset_id": args.dataset_id,
        "annotation_release": args.annotation_release,
        "annotation_status": "native_exact_H1",
        "genetic_code": 1,
        "coordinate_system": {"GFF3": "1-based closed", "BED_PAF": "0-based half-open"},
        "source": {
            "gff_path": str(args.gff.resolve()), "gff_sha256": sha256(args.gff),
            "fai_path": str(args.fai.resolve()), "fai_sha256": sha256(args.fai),
            "bounded_whole_haplotype_paf_path": str(args.bounded_paf.resolve()),
            "bounded_whole_haplotype_paf_sha256": sha256(args.bounded_paf),
            "h1_paf_axis": h1_paf_axis,
        },
        "targeted_gene_count": len(targeted_genes),
        "targeted_gene_union_bases": union_length(targeted_intervals),
        "queryable_gene_count": len(queryable_genes),
        "queryable_gene_union_bases": union_length(queryable_intervals),
        "excluded_gene_count": len(excluded_genes),
        "excluded_gene_union_bases": union_length(excluded_intervals),
        "excluded_gene_reason_counts": dict(sorted(reason_counts.items())),
        "total_CDS_rows": len(cds_rows),
        "queryable_CDS_rows": queryable_cds_rows,
        "excluded_CDS_rows": excluded_cds_rows,
        "queryable_CDS_phase_counts": dict(sorted(phase_counts.items())),
        "phase_validation": "transcript-order GFF3 phase-chain continuity passed for every queryable transcript",
        "execution_span_count": len(span_records),
        "execution_span_union_bases": sum(end - start for _, start, end, _ in span_records),
        "execution_merge_rule": "merge only overlapping or touching queryable gene spans on the same H1 contig",
        "impg_partition_rule": "run impg partition on the whole bounded PAF, then select IMPG-native partition rows intersecting execution spans",
        "sweepga_prepartitioned_by_annotation": False,
        "callable_denominator_provenance": (
            "queryable_gene_union_bases is the predeclared annotation target denominator; final analysis "
            "denominator is its intersection with selected bounded mapping, successfully queried IMPG-native "
            "partitions, exact H1 REF validation, and explicit ambiguity/exclusion filters"
        ),
        "artifacts": {},
    }
    for path in (gene_manifest, query_manifest, execution_bed, feature_map, haplotype_contig_map):
        qc["artifacts"][path.name] = {
            "path": str(path.resolve()), "sha256": sha256(path), "size_bytes": path.stat().st_size
        }
    query_qc.write_text(json.dumps(qc, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(qc, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
