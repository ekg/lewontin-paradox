#!/usr/bin/env python3
"""Build a conservative whole-assembly H1 callable mask and SNV BCF from WFMASH PAF.

The parser is interval-based so chromosome assemblies do not require a Python
object per base. It implements the frozen Tier 3A policy: one-to-one primary
PAF records, extended CIGAR, 100-bp alignment edges, 10-bp indel flanks,
unique target projection, unique query-base projection, and A/C/G/T in both
haplotypes. Query and target overlaps are excluded exactly by interval sweeps.
"""

from __future__ import annotations

import argparse
import bisect
import collections
import hashlib
import json
import re
import subprocess
from pathlib import Path
from typing import Iterable, Iterator, Sequence

import pysam

CIGAR = re.compile(r"(\d+)([=XID])")
DNA = frozenset("ACGT")
EDGE = 100
INDEL_FLANK = 10


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_cigar(value: str) -> list[tuple[int, str]]:
    result = [(int(length), operation) for length, operation in CIGAR.findall(value)]
    if not result or "".join(f"{n}{op}" for n, op in result) != value:
        raise ValueError(f"invalid/non-extended CIGAR: {value[:100]!r}")
    return result


def merge(intervals: Iterable[tuple[int, int]]) -> list[tuple[int, int]]:
    result: list[list[int]] = []
    for start, end in sorted(intervals):
        if start >= end:
            continue
        if result and start <= result[-1][1]:
            result[-1][1] = max(result[-1][1], end)
        else:
            result.append([start, end])
    return [(start, end) for start, end in result]


def subtract(
    intervals: Iterable[tuple[int, int]], exclusions: Sequence[tuple[int, int]]
) -> Iterator[tuple[int, int]]:
    excluded = merge(exclusions)
    for start, end in intervals:
        pieces = [(start, end)]
        for left, right in excluded:
            if right <= start:
                continue
            if left >= end:
                break
            next_pieces = []
            for piece_start, piece_end in pieces:
                if right <= piece_start or left >= piece_end:
                    next_pieces.append((piece_start, piece_end))
                else:
                    if piece_start < left:
                        next_pieces.append((piece_start, left))
                    if right < piece_end:
                        next_pieces.append((right, piece_end))
            pieces = next_pieces
        yield from pieces


def uniquely_covered(intervals: Iterable[tuple[int, int]]) -> tuple[list[tuple[int, int]], int]:
    events: dict[int, int] = collections.defaultdict(int)
    for start, end in intervals:
        events[start] += 1
        events[end] -= 1
    coverage = previous = 0
    unique: list[tuple[int, int]] = []
    multiple = 0
    for position in sorted(events):
        if previous < position:
            if coverage == 1:
                unique.append((previous, position))
            elif coverage > 1:
                multiple += position - previous
        coverage += events[position]
        previous = position
    return merge(unique), multiple


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def run(command: list[str]) -> None:
    completed = subprocess.run(command, text=True, capture_output=True, check=False)
    if completed.returncode:
        raise SystemExit(
            f"command failed ({completed.returncode}): {' '.join(command)}\n{completed.stderr}"
        )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h1", type=Path, required=True)
    parser.add_argument("--h2", type=Path, required=True)
    parser.add_argument("--paf", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--bcftools", type=Path, required=True)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    h1 = pysam.FastaFile(str(args.h1))
    h2 = pysam.FastaFile(str(args.h2))
    h1_lengths = dict(zip(h1.references, h1.lengths))
    h2_lengths = dict(zip(h2.references, h2.lengths))
    records: list[dict] = []
    operations: collections.Counter[str] = collections.Counter()
    ambiguous_mapping_records = 0
    ambiguous_mapping_target_bases = 0
    ambiguous_target_exclusions: dict[str, list[tuple[int, int]]] = collections.defaultdict(list)

    with args.paf.open(encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip().split("\t")
            if len(fields) < 13:
                raise SystemExit(f"PAF line {line_number} has fewer than 13 fields")
            qname, qlen, qstart, qend, strand = fields[:5]
            tname, tlen, tstart, tend = fields[5:9]
            qlen, qstart, qend = map(int, (qlen, qstart, qend))
            tlen, tstart, tend = map(int, (tlen, tstart, tend))
            mapq = int(fields[11])
            if qname not in h2_lengths or h2_lengths[qname] != qlen:
                raise SystemExit(f"PAF query dictionary mismatch at line {line_number}")
            if tname not in h1_lengths or h1_lengths[tname] != tlen:
                raise SystemExit(f"PAF target dictionary mismatch at line {line_number}")
            if strand not in {"+", "-"}:
                raise SystemExit(f"invalid PAF strand at line {line_number}")
            if "tp:A:S" in fields[12:]:
                raise SystemExit(f"secondary PAF line {line_number}")
            tags = [field[5:] for field in fields[12:] if field.startswith("cg:Z:")]
            if len(tags) != 1:
                raise SystemExit(f"PAF line {line_number} lacks exactly one cg:Z tag")
            cigar = parse_cigar(tags[0])
            q_consumed = sum(n for n, op in cigar if op in "=XI")
            t_consumed = sum(n for n, op in cigar if op in "=XD")
            if q_consumed != qend - qstart or t_consumed != tend - tstart:
                raise SystemExit(f"CIGAR span mismatch at PAF line {line_number}")
            if mapq <= 0 or mapq == 255:
                ambiguous_mapping_records += 1
                ambiguous_mapping_target_bases += tend - tstart
                ambiguous_target_exclusions[tname].append((tstart, tend))
                continue
            for length, operation in cigar:
                operations[operation] += length
            records.append({
                "id": len(records), "qname": qname, "qstart": qstart, "qend": qend,
                "tname": tname, "tstart": tstart, "tend": tend, "strand": strand,
                "cigar": cigar,
            })
    if not records:
        raise SystemExit("PAF has no records")

    # Project exact overlaps between mappings of the same H2 bases back to H1
    # and exclude those H1 spans. Non-overlapping portions of each record remain.
    query_exclusions: dict[int, list[tuple[int, int]]] = collections.defaultdict(list)
    by_query: dict[str, list[dict]] = collections.defaultdict(list)
    for record in records:
        by_query[record["qname"]].append(record)
    for values in by_query.values():
        active: list[dict] = []
        for record in sorted(values, key=lambda x: (x["qstart"], x["qend"])):
            active = [old for old in active if old["qend"] > record["qstart"]]
            for old in active:
                overlap = (max(old["qstart"], record["qstart"]), min(old["qend"], record["qend"]))
                query_exclusions[old["id"]].append(overlap)
                query_exclusions[record["id"]].append(overlap)
            active.append(record)

    candidates: dict[str, list[tuple[int, int]]] = collections.defaultdict(list)
    record_candidates: dict[int, list[tuple[int, int]]] = {}
    query_overlap_target_bases = 0
    edge_bases = indel_flank_bases = 0
    for record in records:
        lower, upper = record["tstart"] + EDGE, record["tend"] - EDGE
        if lower >= upper:
            edge_bases += record["tend"] - record["tstart"]
            continue
        edge_bases += 2 * EDGE
        target_position = record["tstart"]
        query_offset = 0
        aligned: list[tuple[int, int]] = []
        indels: list[tuple[int, int]] = []
        nonunique_query_targets: list[tuple[int, int]] = []
        for length, operation in record["cigar"]:
            if operation in "=X":
                aligned.append((max(target_position, lower), min(target_position + length, upper)))
                if record["strand"] == "+":
                    query_left = record["qstart"] + query_offset
                    query_right = query_left + length
                else:
                    query_right = record["qend"] - query_offset
                    query_left = query_right - length
                for exclude_start, exclude_end in query_exclusions.get(record["id"], []):
                    overlap_start = max(query_left, exclude_start)
                    overlap_end = min(query_right, exclude_end)
                    if overlap_start >= overlap_end:
                        continue
                    if record["strand"] == "+":
                        target_left = target_position + overlap_start - query_left
                        target_right = target_position + overlap_end - query_left
                    else:
                        target_left = target_position + query_right - overlap_end
                        target_right = target_position + query_right - overlap_start
                    nonunique_query_targets.append((target_left, target_right))
                query_offset += length
                target_position += length
            elif operation == "D":
                indels.append((max(lower, target_position - INDEL_FLANK), min(upper, target_position + length + INDEL_FLANK)))
                target_position += length
            else:
                # Match range(anchor-flank, anchor+flank+1) in the frozen parser.
                indels.append((max(lower, target_position - INDEL_FLANK), min(upper, target_position + INDEL_FLANK + 1)))
                query_offset += length
        merged_indels = merge(indels)
        nonunique_query_targets = merge(nonunique_query_targets)
        query_overlap_target_bases += sum(end - start for start, end in nonunique_query_targets)
        indel_flank_bases += sum(end - start for start, end in merged_indels)
        retained_for_record = list(
            subtract(
                aligned,
                merged_indels + nonunique_query_targets + ambiguous_target_exclusions.get(record["tname"], []),
            )
        )
        record_candidates[record["id"]] = merge(retained_for_record)
        candidates[record["tname"]].extend(retained_for_record)

    unique: dict[str, list[tuple[int, int]]] = {}
    target_multiple_bases = 0
    for contig, values in candidates.items():
        unique[contig], excluded = uniquely_covered(values)
        target_multiple_bases += excluded

    # Require unambiguous H1 bases. For = operations H2 is identical; for X,
    # H2 ambiguity is checked again while producing the VCF below.
    callable_intervals: dict[str, list[tuple[int, int]]] = {}
    ambiguous_h1_bases = 0
    for contig, values in unique.items():
        retained: list[tuple[int, int]] = []
        for start, end in values:
            sequence = h1.fetch(contig, start, end).upper()
            for match in re.finditer(r"[ACGT]+", sequence):
                retained.append((start + match.start(), start + match.end()))
            ambiguous_h1_bases += len(sequence) - sum(
                match.end() - match.start() for match in re.finditer(r"[ACGT]+", sequence)
            )
        callable_intervals[contig] = merge(retained)

    callable_bed = args.output_dir / "h1.callable.bed"
    callable_bases = 0
    with callable_bed.open("w", encoding="utf-8") as output:
        for contig in h1.references:
            for start, end in callable_intervals.get(contig, []):
                output.write(f"{contig}\t{start}\t{end}\n")
                callable_bases += end - start
    if not callable_bases:
        raise SystemExit("whole-assembly callable denominator is zero")

    starts = {contig: [start for start, _ in values] for contig, values in callable_intervals.items()}
    def is_callable(contig: str, position: int) -> bool:
        values = callable_intervals.get(contig, [])
        index = bisect.bisect_right(starts.get(contig, []), position) - 1
        return index >= 0 and position < values[index][1]

    record_starts = {
        record_id: [start for start, _ in values]
        for record_id, values in record_candidates.items()
    }
    def is_record_candidate(record_id: int, position: int) -> bool:
        values = record_candidates.get(record_id, [])
        index = bisect.bisect_right(record_starts.get(record_id, []), position) - 1
        return index >= 0 and position < values[index][1]

    raw_vcf = args.output_dir / "h1_vs_h2.raw.snvs.vcf"
    snvs = ambiguous_h2 = 0
    ambiguous_h2_positions: dict[str, list[tuple[int, int]]] = collections.defaultdict(list)
    with raw_vcf.open("w", encoding="utf-8") as output:
        output.write("##fileformat=VCFv4.2\n")
        output.write("##source=wfmash-e040aa10-tier3a-interval-parser-v1\n")
        for contig in h1.references:
            output.write(f"##contig=<ID={contig},length={h1_lengths[contig]}>\n")
        output.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        output.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{args.sample}\n")
        for record in records:
            query = h2.fetch(record["qname"], record["qstart"], record["qend"])
            if record["strand"] == "-":
                query = reverse_complement(query)
            query_position = 0
            target_position = record["tstart"]
            for length, operation in record["cigar"]:
                if operation in "=X":
                    if operation == "X":
                        target = h1.fetch(record["tname"], target_position, target_position + length).upper()
                        alternate = query[query_position:query_position + length].upper()
                        for offset, (ref, alt) in enumerate(zip(target, alternate)):
                            position = target_position + offset
                            if not is_callable(record["tname"], position) or not is_record_candidate(record["id"], position):
                                continue
                            if ref not in DNA or alt not in DNA:
                                ambiguous_h2 += 1
                                ambiguous_h2_positions[record["tname"]].append((position, position + 1))
                                continue
                            if ref == alt:
                                raise SystemExit("WFMASH X operation contains an equal base")
                            output.write(
                                f"{record['tname']}\t{position + 1}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT\t0/1\n"
                            )
                            snvs += 1
                    query_position += length
                    target_position += length
                elif operation == "I":
                    query_position += length
                else:
                    target_position += length

    # X columns with an ambiguous H2 allele were initially present in the H1
    # ACGT mask; remove them before freezing the denominator.
    if ambiguous_h2_positions:
        callable_bases = 0
        with callable_bed.open("w", encoding="utf-8") as output:
            for contig in h1.references:
                retained = list(
                    subtract(
                        callable_intervals.get(contig, []),
                        ambiguous_h2_positions.get(contig, []),
                    )
                )
                callable_intervals[contig] = retained
                for start, end in retained:
                    output.write(f"{contig}\t{start}\t{end}\n")
                    callable_bases += end - start

    sorted_bcf = args.output_dir / "h1_vs_h2.sorted.bcf"
    final_bcf = args.output_dir / "h1_vs_h2.normalized.snvs.bcf"
    run([str(args.bcftools), "sort", "-Ob", "-o", str(sorted_bcf), str(raw_vcf)])
    run([
        str(args.bcftools), "norm", "-f", str(args.h1), "-d", "exact", "-Ob",
        "-o", str(final_bcf), str(sorted_bcf),
    ])
    run([str(args.bcftools), "index", "--csi", "-f", str(final_bcf)])
    snvs_before_exact_dedup = snvs
    with pysam.VariantFile(str(final_bcf)) as variants:
        snvs = sum(1 for _ in variants)
    raw_vcf.unlink()
    sorted_bcf.unlink()

    qc = {
        "status": "eligible",
        "policy": "tier3-decisions-v1",
        "parser": "acquisition_build_mask.py-v1",
        "h1_role": "target_reference",
        "h2_role": "query_phased_haplotype",
        "wfmash_commit": "e040aa10e87cab44ed5a4db005e784be62b0bd21",
        "wfmash_arguments": ["-p", "90", "-w", "5k", "-l", "25k", "-o", "-4"],
        "mapping_records": len(records),
        "ambiguous_mapping_records_excluded": ambiguous_mapping_records,
        "records_with_overlapping_query_bases": len(query_exclusions),
        "callable_bases": callable_bases,
        "heterozygous_snvs": snvs,
        "heterozygous_snvs_before_exact_dedup": snvs_before_exact_dedup,
        "operation_counts": dict(sorted(operations.items())),
        "exclusion_counts": {
            "alignment_edge_bases": edge_bases,
            "indel_flank_bases_before_overlap_union": indel_flank_bases,
            "overlapping_query_record_target_bases": query_overlap_target_bases,
            "multiple_target_projection_bases": target_multiple_bases,
            "ambiguous_h1_bases": ambiguous_h1_bases,
            "ambiguous_h2_X_bases": ambiguous_h2,
            "ambiguous_mapping_target_bases_before_union": ambiguous_mapping_target_bases,
            "exact_duplicate_variant_records": snvs_before_exact_dedup - snvs,
        },
        "artifacts": {
            "h1_fasta": {"path": str(args.h1.resolve()), "sha256": sha256(args.h1)},
            "h2_fasta": {"path": str(args.h2.resolve()), "sha256": sha256(args.h2)},
            "raw_paf": {"path": str(args.paf.resolve()), "sha256": sha256(args.paf)},
            "callable_bed": {"path": str(callable_bed.resolve()), "sha256": sha256(callable_bed)},
            "normalized_bcf": {"path": str(final_bcf.resolve()), "sha256": sha256(final_bcf)},
            "normalized_bcf_csi": {"path": str(Path(str(final_bcf) + '.csi').resolve()), "sha256": sha256(Path(str(final_bcf) + '.csi'))},
        },
    }
    qc_path = args.output_dir / "accepted_mapping_qc.json"
    qc_path.write_text(json.dumps(qc, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(qc, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
