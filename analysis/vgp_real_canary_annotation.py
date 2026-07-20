#!/usr/bin/env python3
"""Run exact-native annotation partitions on promoted VGP pilot outputs."""

from __future__ import annotations

import argparse
import bisect
import json
import os
from pathlib import Path
from typing import Iterable, Mapping, Sequence

from analysis.tier3_common import (
    collect_fourfold_sites,
    fasta_dictionary,
    parse_gff,
    read_fasta,
    sha256_file,
)


Interval = tuple[str, int, int]
Variant = tuple[str, int, str, str]  # contig, zero-based position, REF, ALT


class AnnotationCanaryError(RuntimeError):
    """The optional exact-native annotation branch failed a hard binding."""


def parse_bed(path: Path) -> list[Interval]:
    rows: list[Interval] = []
    for number, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if not raw or raw.startswith("#"):
            continue
        fields = raw.split("\t")
        if len(fields) < 3:
            raise AnnotationCanaryError(f"{path}:{number}: truncated BED")
        start, end = int(fields[1]), int(fields[2])
        if start < 0 or end <= start:
            raise AnnotationCanaryError(f"{path}:{number}: invalid BED interval")
        rows.append((fields[0], start, end))
    return rows


def parse_vcf(path: Path) -> list[Variant]:
    rows: list[Variant] = []
    for number, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if not raw or raw.startswith("#"):
            continue
        fields = raw.split("\t")
        if len(fields) < 5:
            raise AnnotationCanaryError(f"{path}:{number}: truncated VCF")
        for alt in fields[4].split(","):
            rows.append((fields[0], int(fields[1]) - 1, fields[3].upper(), alt.upper()))
    return rows


def merge_intervals(rows: Iterable[Interval]) -> list[Interval]:
    result: list[list[object]] = []
    for contig, start, end in sorted(rows, key=lambda row: (row[0], row[1], row[2])):
        if result and result[-1][0] == contig and start <= int(result[-1][2]):
            result[-1][2] = max(int(result[-1][2]), end)
        else:
            result.append([contig, start, end])
    return [(str(contig), int(start), int(end)) for contig, start, end in result]


class IntervalIndex:
    def __init__(self, rows: Iterable[Interval]):
        self.rows: dict[str, list[tuple[int, int]]] = {}
        self.starts: dict[str, list[int]] = {}
        for contig, start, end in merge_intervals(rows):
            self.rows.setdefault(contig, []).append((start, end))
        self.starts = {contig: [start for start, _ in values] for contig, values in self.rows.items()}

    def contains(self, contig: str, start: int, end: int) -> bool:
        starts = self.starts.get(contig, [])
        position = bisect.bisect_right(starts, start) - 1
        return position >= 0 and self.rows[contig][position][0] <= start and end <= self.rows[contig][position][1]

    def intersect_bp(self, other: "IntervalIndex") -> int:
        total = 0
        for contig, left_rows in self.rows.items():
            right_rows = other.rows.get(contig, [])
            left_index = right_index = 0
            while left_index < len(left_rows) and right_index < len(right_rows):
                left_start, left_end = left_rows[left_index]
                right_start, right_end = right_rows[right_index]
                total += max(0, min(left_end, right_end) - max(left_start, right_start))
                if left_end <= right_end:
                    left_index += 1
                else:
                    right_index += 1
        return total


def annotation_partition_counts(
    callable_rows: Sequence[Interval], cds_rows: Sequence[Interval],
    fourfold: set[tuple[str, int]], h1: Mapping[str, str], variants: Sequence[Variant],
) -> list[dict[str, object]]:
    callable_index, cds_index = IntervalIndex(callable_rows), IntervalIndex(cds_rows)
    cds_callable = callable_index.intersect_bp(cds_index)
    cds_variants = [row for row in variants if callable_index.contains(row[0], row[1], row[1] + len(row[2]))
                    and cds_index.contains(row[0], row[1], row[1] + len(row[2]))]
    callable_fourfold = {
        key for key in fourfold if callable_index.contains(key[0], key[1], key[1] + 1)
    }
    fourfold_w = {key for key in callable_fourfold if h1[key[0]][key[1]] in "AT"}
    fourfold_s = {key for key in callable_fourfold if h1[key[0]][key[1]] in "GC"}
    snps = [row for row in variants if len(row[2]) == len(row[3]) == 1]
    snp_keys = {(row[0], row[1]): row for row in snps}
    fourfold_hits = set(snp_keys) & callable_fourfold
    w_hits, s_hits = set(snp_keys) & fourfold_w, set(snp_keys) & fourfold_s
    ws = {key for key in w_hits if snp_keys[key][3] in "GC"}
    sw = {key for key in s_hits if snp_keys[key][3] in "AT"}
    specs = (
        ("CDS", cds_callable, len(cds_variants), "normalized heterozygous allele records / callable CDS bp"),
        ("fourfold", len(callable_fourfold), len(fourfold_hits), "heterozygous SNPs / callable fourfold bp"),
        ("fourfold_W", len(fourfold_w), len(w_hits), "heterozygous SNPs / callable fourfold AT bp"),
        ("fourfold_S", len(fourfold_s), len(s_hits), "heterozygous SNPs / callable fourfold GC bp"),
        ("WS", len(fourfold_w), len(ws), "AT-to-GC SNPs / callable fourfold AT bp"),
        ("SW", len(fourfold_s), len(sw), "GC-to-AT SNPs / callable fourfold GC bp"),
    )
    return [{
        "partition": name,
        "callable_bp": denominator,
        "heterozygous_variants": numerator,
        "estimate": numerator / denominator if denominator else None,
        "estimator": estimator,
    } for name, denominator, numerator, estimator in specs]


def run(args: argparse.Namespace) -> dict[str, object]:
    h1 = read_fasta(args.h1_fasta)
    annotation = parse_gff(args.annotation_gff)
    observed_dictionary = fasta_dictionary(h1)
    if annotation.sequence_regions != observed_dictionary:
        raise AnnotationCanaryError("exact-native GFF sequence dictionary does not equal pilot H1")
    fourfold, frame_discordant = collect_fourfold_sites(h1, annotation, genetic_code=1)
    cds_rows = merge_intervals(
        (segment.contig, segment.start, segment.end)
        for transcript_id in annotation.canonical_transcripts.values()
        for segment in annotation.transcripts[transcript_id].segments
    )
    partitions = annotation_partition_counts(
        parse_bed(args.callable_bed), cds_rows, fourfold, h1, parse_vcf(args.normalized_vcf)
    )
    if not any(row["partition"] == "fourfold" and int(row["callable_bp"]) > 0 for row in partitions):
        raise AnnotationCanaryError("exact annotation produced no callable fourfold partition")
    value = {
        "schema_version": args.schema_version,
        "task_id": args.task_id,
        "authorization_id": args.authorization_id,
        "selection_id": args.selection_id,
        "canonical_vgp_root": args.canonical_root,
        "annotation_status": "exact_native",
        "assembly_accession_version": args.assembly_accession_version,
        "annotation_accession_version": args.annotation_accession_version,
        "annotation_gff": {
            "canonical_source_path": args.annotation_source_path,
            "sha256": sha256_file(args.annotation_gff),
            "size_bytes": args.annotation_gff.stat().st_size,
        },
        "sequence_dictionary_equal": True,
        "sequence_regions": len(annotation.sequence_regions),
        "canonical_transcripts": len(annotation.canonical_transcripts),
        "nuclear_genetic_code": 1,
        "frame_discordant_overlap_positions_excluded": len(frame_discordant),
        "partitions": partitions,
        "gc3_status": "not_reported_without_a_separate_frozen_gc3_overlap_policy",
        "atomic_verified_promotion": True,
        "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    if args.output.exists():
        raise AnnotationCanaryError(f"refusing to overwrite annotation result: {args.output}")
    partial = args.output.with_suffix(args.output.suffix + ".partial")
    partial.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    expected = sha256_file(partial)
    partial.replace(args.output)
    if sha256_file(args.output) != expected:
        raise AnnotationCanaryError("annotation output digest changed during atomic promotion")
    return value


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--h1-fasta", type=Path, required=True)
    result.add_argument("--annotation-gff", type=Path, required=True)
    result.add_argument("--annotation-source-path", required=True)
    result.add_argument("--callable-bed", type=Path, required=True)
    result.add_argument("--normalized-vcf", type=Path, required=True)
    result.add_argument("--canonical-root", default="/moosefs/erikg/vgp")
    result.add_argument("--selection-id", default="P07")
    result.add_argument("--assembly-accession-version", default="GCA_048126635.1")
    result.add_argument("--annotation-accession-version", default="GCA_048126635.1-GB_2025_08_04")
    result.add_argument("--authorization-id", default="vgp10-auth-20260718-v2")
    result.add_argument("--task-id", default="run-vgp-real-canary")
    result.add_argument("--schema-version", default="vgp-real-canary-exact-annotation-v1")
    result.add_argument("--output", type=Path, required=True)
    return result


def main(argv: Sequence[str] | None = None) -> int:
    args = parser().parse_args(argv)
    try:
        value = run(args)
    except (AnnotationCanaryError, OSError, ValueError) as error:
        print(f"ERROR: {error}")
        return 2
    print(json.dumps({"output": str(args.output), "partitions": value["partitions"]}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
