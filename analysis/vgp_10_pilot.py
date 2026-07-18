#!/usr/bin/env python3
"""Fail-closed implementation helpers for the ten-pair VGP pilot.

This module deliberately keeps scientific accounting in small, deterministic
Python functions while the production entry points invoke SweepGA, IMPG,
bcftools, and PSMC from the pinned Guix profile.  Coordinates are zero-based,
half-open everywhere except while parsing VCF records.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import random
import re
import shlex
import shutil
import sys
import tempfile
from collections import defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable, Mapping, Sequence


ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = ROOT / "analysis"
PRIMARY_MANIFEST = ANALYSIS / "vgp_10_pair_manifest.tsv"
DESIGN_MANIFEST = ANALYSIS / "vgp_analysis_manifest.json"
ENVIRONMENT_LOCK = ANALYSIS / "guix/vgp_10_pilot/environment-lock.json"

REASON_ORDER = (
    "not_eligible_contig",
    "organellar",
    "sex_linked_primary_exclusion",
    "unplaced_or_unlocalized",
    "h1_gap_or_N",
    "h2_gap_or_N",
    "not_1to1",
    "mapping_breakpoint",
    "low_base_accuracy",
    "repeat_or_low_complexity_primary",
    "duplication_or_collapse",
    "phase_uncertain",
    "other_predeclared",
)

IUPAC = {
    frozenset(("A", "G")): "R",
    frozenset(("C", "T")): "Y",
    frozenset(("G", "C")): "S",
    frozenset(("A", "T")): "W",
    frozenset(("G", "T")): "K",
    frozenset(("A", "C")): "M",
}


class PilotError(ValueError):
    """A hard, reason-preserving pilot contract failure."""


@dataclass(frozen=True, order=True)
class Interval:
    contig: str
    start: int
    end: int

    def __post_init__(self) -> None:
        if not self.contig or self.start < 0 or self.end <= self.start:
            raise PilotError(f"invalid interval: {self}")

    @property
    def length(self) -> int:
        return self.end - self.start


@dataclass(frozen=True)
class PafRecord:
    query: str
    query_length: int
    query_start: int
    query_end: int
    strand: str
    target: str
    target_length: int
    target_start: int
    target_end: int
    matches: int
    block_length: int
    mapq: int
    tags: tuple[str, ...]


@dataclass(frozen=True)
class Variant:
    contig: str
    pos0: int
    ref: str
    alt: str

    def __post_init__(self) -> None:
        if self.pos0 < 0 or not self.ref or not self.alt:
            raise PilotError(f"invalid variant: {self}")
        allowed = set("ACGTN")
        if set(self.ref.upper()) - allowed or set(self.alt.upper()) - allowed:
            raise PilotError(f"non-sequence or symbolic allele: {self}")


def sha256_file(path: Path | str) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def canonical_json(value: object) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n"


def read_bed(path: Path | str) -> list[Interval]:
    rows: list[Interval] = []
    with Path(path).open(encoding="utf-8") as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise PilotError(f"{path}:{number}: BED requires three columns")
            rows.append(Interval(fields[0], int(fields[1]), int(fields[2])))
    return rows


def write_bed(path: Path | str, intervals: Iterable[Interval]) -> None:
    Path(path).write_text(
        "".join(f"{item.contig}\t{item.start}\t{item.end}\n" for item in intervals),
        encoding="utf-8",
    )


def merge_intervals(intervals: Iterable[Interval]) -> list[Interval]:
    merged: list[Interval] = []
    for current in sorted(intervals):
        if merged and merged[-1].contig == current.contig and current.start <= merged[-1].end:
            previous = merged[-1]
            merged[-1] = Interval(previous.contig, previous.start, max(previous.end, current.end))
        else:
            merged.append(current)
    return merged


def intersect_intervals(left: Iterable[Interval], right: Iterable[Interval]) -> list[Interval]:
    by_contig: dict[str, list[Interval]] = defaultdict(list)
    for item in merge_intervals(right):
        by_contig[item.contig].append(item)
    result: list[Interval] = []
    for a in merge_intervals(left):
        for b in by_contig[a.contig]:
            if b.start >= a.end:
                break
            start, end = max(a.start, b.start), min(a.end, b.end)
            if start < end:
                result.append(Interval(a.contig, start, end))
    return result


def subtract_intervals(universe: Iterable[Interval], remove: Iterable[Interval]) -> list[Interval]:
    by_contig: dict[str, list[Interval]] = defaultdict(list)
    for item in merge_intervals(remove):
        by_contig[item.contig].append(item)
    result: list[Interval] = []
    for item in merge_intervals(universe):
        cursor = item.start
        for excluded in by_contig[item.contig]:
            if excluded.end <= cursor:
                continue
            if excluded.start >= item.end:
                break
            if cursor < excluded.start:
                result.append(Interval(item.contig, cursor, min(excluded.start, item.end)))
            cursor = max(cursor, excluded.end)
            if cursor >= item.end:
                break
        if cursor < item.end:
            result.append(Interval(item.contig, cursor, item.end))
    return result


def interval_bp(intervals: Iterable[Interval]) -> int:
    return sum(item.length for item in merge_intervals(intervals))


def parse_fasta(path: Path | str) -> dict[str, str]:
    result: dict[str, list[str]] = {}
    current: str | None = None
    with Path(path).open(encoding="utf-8") as handle:
        for number, line in enumerate(handle, 1):
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split()[0]
                if not current or current in result:
                    raise PilotError(f"{path}:{number}: duplicate/empty FASTA name")
                result[current] = []
            elif current is None:
                raise PilotError(f"{path}:{number}: sequence before FASTA header")
            else:
                result[current].append(line.upper())
    if not result:
        raise PilotError(f"empty FASTA: {path}")
    return {name: "".join(parts) for name, parts in result.items()}


def sequence_dictionary(sequences: Mapping[str, str]) -> list[dict[str, object]]:
    return [
        {"name": name, "length": len(sequence), "md5": hashlib.md5(sequence.encode()).hexdigest()}
        for name, sequence in sequences.items()
    ]


def parse_paf(path: Path | str) -> list[PafRecord]:
    records: list[PafRecord] = []
    with Path(path).open(encoding="utf-8") as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                raise PilotError(f"{path}:{number}: PAF has fewer than 12 fields")
            record = PafRecord(
                fields[0], int(fields[1]), int(fields[2]), int(fields[3]), fields[4],
                fields[5], int(fields[6]), int(fields[7]), int(fields[8]),
                int(fields[9]), int(fields[10]), int(fields[11]), tuple(fields[12:]),
            )
            if record.strand not in {"+", "-"}:
                raise PilotError(f"{path}:{number}: invalid strand")
            if not (0 <= record.query_start < record.query_end <= record.query_length):
                raise PilotError(f"{path}:{number}: invalid query coordinates")
            if not (0 <= record.target_start < record.target_end <= record.target_length):
                raise PilotError(f"{path}:{number}: invalid target coordinates")
            if sum(tag.startswith("cg:Z:") for tag in record.tags) != 1:
                raise PilotError(f"{path}:{number}: exactly one cg:Z tag is required")
            records.append(record)
    if not records:
        raise PilotError("SweepGA emitted an empty alignment")
    return records


def paf_h1_intervals(records: Sequence[PafRecord]) -> list[Interval]:
    return merge_intervals(Interval(row.target, row.target_start, row.target_end) for row in records)


def non_acgt_intervals(sequences: Mapping[str, str]) -> list[Interval]:
    return [
        Interval(contig, match.start(), match.end())
        for contig, sequence in sequences.items()
        for match in re.finditer(r"[^ACGT]+", sequence)
    ]


def low_complexity_intervals(
    sequences: Mapping[str, str], homopolymer_bases: int = 10,
    dinucleotide_copies: int = 10, trinucleotide_copies: int = 7,
) -> list[Interval]:
    """Generate a conservative sequence-derived simple-repeat mask.

    This is deliberately available when no standalone repeat report exists.
    It masks long homopolymers and exact di-/trinucleotide tandem runs; richer
    repeat annotations remain optional confidence evidence and may be added as
    separately reason-coded sensitivity inputs.
    """
    if min(homopolymer_bases, dinucleotide_copies, trinucleotide_copies) < 2:
        raise PilotError("low-complexity thresholds must be at least two")
    rows: list[Interval] = []
    patterns = (
        re.compile(rf"([ACGT])\1{{{homopolymer_bases - 1},}}"),
        re.compile(rf"([ACGT]{{2}})\1{{{dinucleotide_copies - 1},}}"),
        re.compile(rf"([ACGT]{{3}})\1{{{trinucleotide_copies - 1},}}"),
    )
    for contig, sequence in sequences.items():
        for pattern in patterns:
            rows.extend(Interval(contig, match.start(), match.end()) for match in pattern.finditer(sequence))
    return merge_intervals(rows)


def project_h2_non_acgt_to_h1(
    records: Sequence[PafRecord], h2_sequences: Mapping[str, str]
) -> list[Interval]:
    """Project query non-ACGT and alignment deletions through exact cg:Z PAF."""
    query_bad: dict[str, list[Interval]] = defaultdict(list)
    for interval in non_acgt_intervals(h2_sequences):
        query_bad[interval.contig].append(interval)
    projected: list[Interval] = []
    for row in records:
        if row.query not in h2_sequences or len(h2_sequences[row.query]) != row.query_length:
            raise PilotError(f"PAF query dictionary differs from H2 FASTA: {row.query}")
        cigar = next(tag[5:] for tag in row.tags if tag.startswith("cg:Z:"))
        operations = [(int(length), op) for length, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar)]
        if "".join(f"{length}{op}" for length, op in operations) != cigar:
            raise PilotError("unsupported or malformed cg:Z CIGAR")
        query_cursor = row.query_start if row.strand == "+" else row.query_end
        target_cursor = row.target_start
        for length, op in operations:
            if length <= 0:
                raise PilotError("zero-length CIGAR operation")
            consumes_query = op in "M=XI"
            consumes_target = op in "M=XD"
            if op not in "M=XID":
                raise PilotError(f"unsupported cg:Z operation for mask projection: {op}")
            if consumes_query and consumes_target:
                if row.strand == "+":
                    qstart, qend = query_cursor, query_cursor + length
                else:
                    qstart, qend = query_cursor - length, query_cursor
                for bad in query_bad[row.query]:
                    start, end = max(qstart, bad.start), min(qend, bad.end)
                    if start >= end:
                        continue
                    if row.strand == "+":
                        tstart = target_cursor + start - qstart
                        tend = target_cursor + end - qstart
                    else:
                        tstart = target_cursor + qend - end
                        tend = target_cursor + qend - start
                    projected.append(Interval(row.target, tstart, tend))
            elif op == "D":
                # Target sequence has no corresponding H2 base.
                projected.append(Interval(row.target, target_cursor, target_cursor + length))
            if consumes_query:
                query_cursor += length if row.strand == "+" else -length
            if consumes_target:
                target_cursor += length
        expected_query = row.query_end if row.strand == "+" else row.query_start
        if query_cursor != expected_query or target_cursor != row.target_end:
            raise PilotError("cg:Z consumption does not match PAF coordinate spans")
    return merge_intervals(projected)


def _maximum_depth(records: Sequence[PafRecord], axis: str) -> int:
    groups: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for record in records:
        if axis == "query":
            groups[record.query].append((record.query_start, record.query_end))
        elif axis == "target":
            groups[record.target].append((record.target_start, record.target_end))
        else:
            raise PilotError(f"unknown PAF axis: {axis}")
    maximum = 0
    for values in groups.values():
        # End events sort before start events at equal half-open coordinates.
        events = sorted([(start, 1) for start, _ in values] + [(end, -1) for _, end in values],
                        key=lambda item: (item[0], item[1]))
        depth = 0
        for _, delta in events:
            depth += delta
            maximum = max(maximum, depth)
    return maximum


def audit_sweepga_paf(path: Path | str, h1_names: set[str], h2_names: set[str]) -> dict[str, object]:
    records = parse_paf(path)
    if any(record.query not in h2_names or record.target not in h1_names for record in records):
        raise PilotError("PAF orientation is not exact H2-query to H1-target")
    query_depth = _maximum_depth(records, "query")
    target_depth = _maximum_depth(records, "target")
    if query_depth > 1 or target_depth > 1:
        raise PilotError(
            f"retained SweepGA multiplicity exceeds 1:1: query={query_depth}, target={target_depth}"
        )
    return {
        "orientation": "H1_reference_H2_query",
        "records": len(records),
        "maximum_query_overlap_depth": query_depth,
        "maximum_target_overlap_depth": target_depth,
        "paf_sha256": sha256_file(path),
    }


def sweepga_command(sweepga: str, h1: str, h2: str, output: str, threads: int) -> list[str]:
    # SweepGA accepts sequence files directly; H2 is the first/query sequence,
    # H1 the second/target sequence.  The audit above proves the emitted roles.
    return [
        sweepga, h2, h1, "--output-file", output, "--num-mappings", "1:1",
        "--scaffold-jump", "0", "--overlap", "0.95", "--scoring",
        "log-length-ani", "--threads", str(threads),
    ]


def impg_commands(
    impg: str, paf: str, index: str, partitions_dir: str, focus_bed: str,
    calls_dir: str, vcf_list: str, laced_vcf: str, h1: str, h2: str,
    threads: int,
) -> list[list[str]]:
    return [
        [impg, "index", "-a", paf, "-i", index, "-t", str(threads)],
        [impg, "partition", "-a", paf, "-i", index, "-d", "0", "-o", "bed",
         "--output-folder", partitions_dir, "-t", str(threads)],
        [impg, "query", "-a", paf, "-i", index, "-b", focus_bed, "-d", "0",
         "-o", "vcf:poa", "--sequence-files", h1, h2, "-O", calls_dir,
         "-t", str(threads)],
        [impg, "lace", "-l", vcf_list, "--format", "vcf", "-o", laced_vcf,
         "--reference", h1, "--compress", "none", "-t", str(threads)],
    ]


def bcftools_commands(bcftools: str, laced: str, h1: str, one_to_one_bed: str,
                      output_prefix: str) -> list[list[str]]:
    # Two norm passes make ownership explicit: decomposition/left-normalization,
    # exact region trim, then exact duplicate removal and final serialization.
    raw = output_prefix + ".split.vcf.gz"
    trimmed = output_prefix + ".trimmed.vcf.gz"
    vcfgz = output_prefix + ".vcf.gz"
    bcf = output_prefix + ".bcf"
    return [
        [bcftools, "norm", "-f", h1, "-m", "-any", "-Oz", "-o", raw, laced],
        [bcftools, "view", "-R", one_to_one_bed, "-Oz", "-o", trimmed, raw],
        [bcftools, "norm", "-f", h1, "-d", "exact", "-Oz", "-o", vcfgz, trimmed],
        [bcftools, "index", "-f", "-t", vcfgz],
        [bcftools, "norm", "-f", h1, "-d", "exact", "-Ob", "-o", bcf, trimmed],
        [bcftools, "index", "-f", bcf],
    ]


def assert_tool_roles(stage: str, argv: Sequence[str]) -> None:
    text = " ".join(argv).lower()
    if stage == "mapping":
        if Path(argv[0]).name != "sweepga" or "--num-mappings 1:1" not in text:
            raise PilotError("mapping must be owned by SweepGA with native --num-mappings 1:1")
        if any(token in text for token in ("vcf", "variant", "call")):
            raise PilotError("SweepGA is a whole-haplotype mapper, never a variant caller")
    elif stage in {"index", "partition", "query", "lace"}:
        if Path(argv[0]).name != "impg" or stage not in argv:
            raise PilotError(f"IMPG must own the {stage} stage")
        if stage == "query" and not {"--sequence-files", "-b"}.issubset(set(argv)):
            raise PilotError("IMPG query requires both sequence files and selected native partitions")
    elif stage == "normalize":
        if Path(argv[0]).name != "bcftools" or "norm" not in argv:
            raise PilotError("bcftools must own normalization and deduplication")
    else:
        raise PilotError(f"unknown tool role: {stage}")


def select_native_partitions(
    partitions_path: Path | str, eligible_regions: Iterable[Interval], output_path: Path | str
) -> list[tuple[Interval, str]]:
    native: list[tuple[Interval, str]] = []
    with Path(partitions_path).open(encoding="utf-8") as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                raise PilotError(f"native partition row {number} lacks IMPG partition identifier")
            native.append((Interval(fields[0], int(fields[1]), int(fields[2])), fields[3]))
    if not native:
        raise PilotError("IMPG native partitions are empty")
    eligible = merge_intervals(eligible_regions)
    selected: list[tuple[Interval, str]] = []
    for interval, name in native:
        if intersect_intervals([interval], eligible):
            selected.append((interval, name))
    if not selected:
        raise PilotError("no IMPG native partition intersects eligible query regions")
    Path(output_path).write_text(
        "".join(f"{item.contig}\t{item.start}\t{item.end}\t{name}\n" for item, name in selected),
        encoding="utf-8",
    )
    return selected


def construct_reason_mask(
    universe: Iterable[Interval], flags: Mapping[str, Iterable[Interval]]
) -> tuple[list[Interval], dict[str, list[Interval]], dict[str, object]]:
    unknown = set(flags) - set(REASON_ORDER)
    if unknown:
        raise PilotError(f"unregistered exclusion reasons: {sorted(unknown)}")
    universe_rows = merge_intervals(universe)
    if not universe_rows:
        raise PilotError("declared H1 coordinate universe is empty")
    ordered_flags = {reason: merge_intervals(flags.get(reason, [])) for reason in REASON_ORDER}
    breakpoints: dict[str, set[int]] = defaultdict(set)
    for interval in universe_rows:
        breakpoints[interval.contig].update((interval.start, interval.end))
    for values in ordered_flags.values():
        for flag in values:
            for universe_interval in universe_rows:
                if universe_interval.contig == flag.contig:
                    start = max(flag.start, universe_interval.start)
                    end = min(flag.end, universe_interval.end)
                    if start < end:
                        breakpoints[flag.contig].update((start, end))
    accounting: dict[str, list[Interval]] = {reason: [] for reason in REASON_ORDER}
    callable_rows: list[Interval] = []
    for universe_interval in universe_rows:
        points = sorted(point for point in breakpoints[universe_interval.contig]
                        if universe_interval.start <= point <= universe_interval.end)
        for start, end in zip(points, points[1:]):
            if start == end:
                continue
            reason = next((
                reason for reason in REASON_ORDER
                if any(flag.contig == universe_interval.contig and flag.start < end and start < flag.end
                       for flag in ordered_flags[reason])
            ), None)
            target = callable_rows if reason is None else accounting[reason]
            target.append(Interval(universe_interval.contig, start, end))
    callable_rows = merge_intervals(callable_rows)
    accounting = {reason: merge_intervals(values) for reason, values in accounting.items()}
    universe_bp = interval_bp(universe_rows)
    callable_bp = interval_bp(callable_rows)
    reason_bp = {reason: interval_bp(values) for reason, values in accounting.items()}
    if callable_bp + sum(reason_bp.values()) != universe_bp:
        raise PilotError("callable plus reason-coded complement does not reconcile to universe")
    reconciliation = {
        "coordinate_system": "H1_zero_based_half_open",
        "reason_order": list(REASON_ORDER),
        "universe_bp": universe_bp,
        "callable_bp": callable_bp,
        "callable_fraction": callable_bp / universe_bp,
        "excluded_bp_by_primary_reason": reason_bp,
        "accounting_discrepancy_bp": 0,
    }
    return callable_rows, accounting, reconciliation


def parse_vcf(path: Path | str) -> list[Variant]:
    variants: list[Variant] = []
    seen: set[Variant] = set()
    with Path(path).open(encoding="utf-8") as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                raise PilotError(f"{path}:{number}: truncated VCF")
            for alt in fields[4].split(","):
                variant = Variant(fields[0], int(fields[1]) - 1, fields[3].upper(), alt.upper())
                if variant in seen:
                    raise PilotError("normalized VCF contains an exact duplicate")
                seen.add(variant)
                variants.append(variant)
    return sorted(variants, key=lambda item: (item.contig, item.pos0, item.ref, item.alt))


def validate_ref_and_reconstruct_h2(
    h1: Mapping[str, str], h2: Mapping[str, str], variants: Sequence[Variant],
    contig_map: Mapping[str, str] | None = None,
    aligned_regions: Sequence[Mapping[str, object]] | None = None,
) -> dict[str, object]:
    contig_map = contig_map or {name: name for name in h1}
    by_contig: dict[str, list[Variant]] = defaultdict(list)
    for variant in variants:
        if variant.contig not in h1:
            raise PilotError(f"variant contig absent from H1: {variant.contig}")
        observed = h1[variant.contig][variant.pos0:variant.pos0 + len(variant.ref)]
        if observed != variant.ref:
            raise PilotError(f"H1 REF mismatch at {variant.contig}:{variant.pos0 + 1}")
        by_contig[variant.contig].append(variant)
    def apply(sequence: str, selected: Sequence[Variant], offset: int = 0) -> str:
        cursor = 0
        output: list[str] = []
        for variant in sorted(selected, key=lambda item: item.pos0):
            position = variant.pos0 - offset
            if position < cursor:
                raise PilotError(f"overlapping normalized variants on {variant.contig}")
            output.extend((sequence[cursor:position], variant.alt))
            cursor = position + len(variant.ref)
        output.append(sequence[cursor:])
        return "".join(output)

    reconstructed: dict[str, str] = {}
    assigned: set[Variant] = set()
    if aligned_regions is not None:
        for ordinal, region in enumerate(aligned_regions):
            h1_name, h2_name = str(region["h1_contig"]), str(region["h2_contig"])
            h1_start, h1_end = int(region["h1_start"]), int(region["h1_end"])
            h2_start, h2_end = int(region["h2_start"]), int(region["h2_end"])
            strand = str(region.get("strand", "+"))
            if h1_name not in h1 or h2_name not in h2 or strand not in {"+", "-"}:
                raise PilotError("invalid manifest-bound concordance region")
            selected = [v for v in by_contig[h1_name]
                        if h1_start <= v.pos0 and v.pos0 + len(v.ref) <= h1_end]
            assigned.update(selected)
            rebuilt = apply(h1[h1_name][h1_start:h1_end], selected, h1_start)
            observed = h2[h2_name][h2_start:h2_end]
            if strand == "-":
                observed = observed.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]
            if rebuilt != observed:
                raise PilotError(f"alternate-sequence reconstruction failed for region {ordinal}")
            reconstructed[f"region_{ordinal}"] = rebuilt
        if assigned != set(variants):
            raise PilotError("not every normalized variant belongs to a validated concordance region")
    else:
        for h1_name, sequence in h1.items():
            reconstructed[h1_name] = apply(sequence, by_contig[h1_name])
            h2_name = contig_map.get(h1_name)
            if h2_name is None or h2_name not in h2:
                raise PilotError(f"no manifest-bound H2 contig for {h1_name}")
            if reconstructed[h1_name] != h2[h2_name]:
                raise PilotError(f"alternate-sequence reconstruction failed for {h1_name}")
    return {
        "h1_ref_checks": len(variants),
        "h1_ref_mismatches": 0,
        "h2_contigs_reconstructed": len(reconstructed),
        "h2_reconstruction_failures": 0,
        "reconstructed_sha256": hashlib.sha256(
            canonical_json(reconstructed).encode()
        ).hexdigest(),
    }


def _position_set(intervals: Iterable[Interval]) -> dict[str, list[tuple[int, int]]]:
    result: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for item in merge_intervals(intervals):
        result[item.contig].append((item.start, item.end))
    return result


def _contained(intervals: Mapping[str, list[tuple[int, int]]], contig: str, start: int, end: int) -> bool:
    return any(left <= start and end <= right for left, right in intervals.get(contig, []))


def build_diploid_consensus(
    h1: Mapping[str, str], callable_intervals: Iterable[Interval], variants: Sequence[Variant],
    indel_flank: int = 10,
) -> tuple[dict[str, str], dict[str, object]]:
    if indel_flank < 0:
        raise PilotError("indel flank cannot be negative")
    callable_map = _position_set(callable_intervals)
    output = {name: ["N"] * len(sequence) for name, sequence in h1.items()}
    for name, regions in callable_map.items():
        if name not in h1:
            raise PilotError(f"callable contig absent from H1: {name}")
        for start, end in regions:
            if end > len(h1[name]):
                raise PilotError(f"callable interval exceeds H1 sequence: {name}")
            output[name][start:end] = list(h1[name][start:end])
    masked_indel_bp: dict[str, set[int]] = defaultdict(set)
    snps = indels = 0
    for variant in variants:
        if not _contained(callable_map, variant.contig, variant.pos0, variant.pos0 + len(variant.ref)):
            raise PilotError("variant is not wholly supported by callable H1 sequence")
        if len(variant.ref) == len(variant.alt) == 1:
            code = IUPAC.get(frozenset((variant.ref, variant.alt)))
            if code is None:
                raise PilotError("heterozygous SNP must contain two distinct canonical alleles")
            output[variant.contig][variant.pos0] = code
            snps += 1
        else:
            start = max(0, variant.pos0 - indel_flank)
            end = min(len(h1[variant.contig]), variant.pos0 + len(variant.ref) + indel_flank)
            for position in range(start, end):
                output[variant.contig][position] = "N"
                masked_indel_bp[variant.contig].add(position)
            indels += 1
    consensus = {name: "".join(sequence) for name, sequence in output.items()}
    expected_callable = interval_bp(callable_intervals) - sum(len(values) for values in masked_indel_bp.values())
    observed_callable = sum(base != "N" for sequence in consensus.values() for base in sequence)
    if observed_callable != expected_callable:
        raise PilotError("diploid consensus callable positions do not reconcile to mask and indel policy")
    return consensus, {
        "masked_base_encoding": "N",
        "heterozygous_snp_encoding": "IUPAC",
        "heterozygous_snps": snps,
        "heterozygous_indels": indels,
        "indel_flank_bp": indel_flank,
        "indel_masked_h1_bp": sum(len(values) for values in masked_indel_bp.values()),
        "consensus_callable_bp": observed_callable,
        "non_callable_bases_encoded_homozygous_reference": 0,
    }


def consensus_to_psmcfa(consensus: Mapping[str, str], bin_size: int = 100) -> str:
    if bin_size <= 0:
        raise PilotError("PSMC bin size must be positive")
    records: list[str] = []
    for name, sequence in consensus.items():
        encoded: list[str] = []
        for start in range(0, len(sequence), bin_size):
            window = sequence[start:start + bin_size]
            if len(window) < bin_size or "N" in window:
                encoded.append("N")
            elif any(base in set(IUPAC.values()) for base in window):
                encoded.append("K")  # canonical fq2psmcfa heterozygous bin
            else:
                encoded.append("T")  # canonical fq2psmcfa homozygous callable bin
        records.append(f">{name}\n")
        records.extend("".join(encoded)[offset:offset + 60] + "\n" for offset in range(0, len(encoded), 60))
    return "".join(records)


def bootstrap_psmcfa(
    consensus: Mapping[str, str], units: Sequence[Interval], sampled_unit_indices: Sequence[int],
    bin_size: int = 100,
) -> str:
    """Encode sampled blocks as distinct records so no PSMC bin spans a boundary."""
    sampled: dict[str, str] = {}
    for ordinal, index in enumerate(sampled_unit_indices, 1):
        if index < 0 or index >= len(units):
            raise PilotError("bootstrap unit index is out of range")
        unit = units[index]
        if unit.contig not in consensus or unit.end > len(consensus[unit.contig]):
            raise PilotError("bootstrap unit exceeds consensus sequence")
        sampled[f"block_{ordinal:06d}_{unit.contig}_{unit.start}_{unit.end}"] = \
            consensus[unit.contig][unit.start:unit.end]
    return consensus_to_psmcfa(sampled, bin_size=bin_size)


def write_fasta(path: Path | str, sequences: Mapping[str, str], width: int = 60) -> None:
    lines: list[str] = []
    for name, sequence in sequences.items():
        lines.append(f">{name}\n")
        lines.extend(sequence[start:start + width] + "\n" for start in range(0, len(sequence), width))
    Path(path).write_text("".join(lines), encoding="utf-8")


def freeze_bootstrap_units(callable_intervals: Iterable[Interval], block_bp: int) -> list[Interval]:
    if block_bp <= 0:
        raise PilotError("bootstrap block length must be positive")
    units: list[Interval] = []
    for interval in merge_intervals(callable_intervals):
        for start in range(interval.start, interval.end, block_bp):
            units.append(Interval(interval.contig, start, min(start + block_bp, interval.end)))
    if not units:
        raise PilotError("no bootstrap units can be frozen from an empty callable mask")
    return units


def bootstrap_manifest(
    callable_intervals: Iterable[Interval], selection_id: str, design_sha256: str,
    attempts: int = 200, block_bp: int = 5_000_000,
) -> list[dict[str, object]]:
    if attempts < 100:
        raise PilotError("PSMC requires at least 100 block bootstrap attempts")
    units = freeze_bootstrap_units(callable_intervals, block_bp)
    seed_material = f"{design_sha256}:{selection_id}:psmc:{block_bp}"
    seed = int(hashlib.sha256(seed_material.encode()).hexdigest()[:16], 16)
    result: list[dict[str, object]] = []
    for replicate in range(1, attempts + 1):
        rng = random.Random(seed + replicate)
        sampled = [rng.randrange(len(units)) for _ in units]
        result.append({
            "replicate": replicate,
            "block_bp": block_bp,
            "seed": seed + replicate,
            "unit_count": len(units),
            "sampled_unit_indices": sampled,
        })
    return result


def scale_unscaled_trajectory(
    rows: Sequence[Mapping[str, float]], scenarios: Sequence[Mapping[str, object]]
) -> tuple[list[dict[str, float]], list[dict[str, object]]]:
    unscaled = [
        {"interval": int(row["interval"]), "time_2N0": float(row["time_2N0"]),
         "lambda": float(row["lambda"])} for row in rows
    ]
    scaled: list[dict[str, object]] = []
    for scenario in scenarios:
        mutation_rate = float(scenario["mutation_rate_per_generation"])
        generation_time = float(scenario["generation_time_years"])
        n0 = float(scenario["theta_0"]) / (4.0 * mutation_rate)
        for row in unscaled:
            scaled.append({
                "scenario_id": str(scenario["scenario_id"]),
                "interval": row["interval"],
                "time_years": row["time_2N0"] * 2.0 * n0 * generation_time,
                "effective_size": row["lambda"] * n0,
                "mutation_rate_per_generation": mutation_rate,
                "generation_time_years": generation_time,
                "mutation_rate_source": str(scenario["mutation_rate_source"]),
                "generation_time_source": str(scenario["generation_time_source"]),
            })
    return unscaled, scaled


def parse_psmc_unscaled(path: Path | str) -> tuple[list[dict[str, float]], float]:
    """Read the final PSMC iteration without applying mutation/time scaling."""
    by_iteration: dict[int, list[dict[str, float]]] = defaultdict(list)
    theta: dict[int, float] = {}
    with Path(path).open(encoding="utf-8") as handle:
        for number, line in enumerate(handle, 1):
            fields = line.split()
            if not fields:
                continue
            try:
                if fields[0] == "RS" and len(fields) >= 5:
                    by_iteration[int(fields[1])].append({
                        "interval": int(fields[2]), "time_2N0": float(fields[3]),
                        "lambda": float(fields[4]),
                    })
                elif fields[0] == "TR" and len(fields) >= 3:
                    theta[int(fields[1])] = float(fields[2])
            except ValueError as error:
                raise PilotError(f"{path}:{number}: malformed PSMC numeric record") from error
    if not by_iteration:
        raise PilotError("PSMC output contains no RS trajectory records")
    iteration = max(by_iteration)
    if iteration not in theta or not (theta[iteration] > 0):
        raise PilotError("final PSMC iteration lacks a positive unscaled theta_0")
    rows = sorted(by_iteration[iteration], key=lambda row: row["interval"])
    if len({row["interval"] for row in rows}) != len(rows):
        raise PilotError("final PSMC trajectory repeats an interval")
    return rows, theta[iteration]


def validate_annotation_binding(
    h1_accession_version: str, h1_dictionary: Sequence[Mapping[str, object]],
    annotation: Mapping[str, object] | None,
) -> dict[str, object]:
    if annotation is None:
        return {"status": "not_available", "core_eligible": True, "annotation_outputs_allowed": False}
    if not annotation.get("annotation_accession_version"):
        raise PilotError("annotation accession/version is required")
    if len(str(annotation.get("gff_sha256", ""))) != 64:
        raise PilotError("annotation GFF SHA256 is required")
    if annotation.get("assembly_accession_version") == h1_accession_version:
        if annotation.get("sequence_dictionary") != list(h1_dictionary):
            raise PilotError("exact native annotation sequence dictionary does not equal H1 dictionary")
        return {"status": "exact_native", "core_eligible": True, "annotation_outputs_allowed": True,
                "annotation_accession_version": annotation["annotation_accession_version"]}
    liftover = annotation.get("validated_liftover")
    if not isinstance(liftover, Mapping):
        raise PilotError("annotation substitution is forbidden without a manifest-bound validated liftover")
    required = {"source_accession_version", "target_accession_version", "chain_sha256",
                "validation_sha256", "manifest_sha256", "passed"}
    if set(liftover) < required or liftover.get("target_accession_version") != h1_accession_version:
        raise PilotError("validated liftover is incomplete or targets a different H1 accession")
    if liftover.get("passed") is not True or any(
        len(str(liftover[key])) != 64 for key in ("chain_sha256", "validation_sha256", "manifest_sha256")
    ):
        raise PilotError("manifest-bound liftover validation failed")
    return {"status": "validated_liftover", "core_eligible": True, "annotation_outputs_allowed": True,
            "annotation_accession_version": annotation["annotation_accession_version"]}


def summarize_annotation_partitions(
    callable_intervals: Iterable[Interval], variants: Sequence[Variant],
    feature_intervals: Mapping[str, Iterable[Interval]],
) -> list[dict[str, object]]:
    """Intersect exact-native/validated-liftover feature partitions with core data."""
    permitted = {"CDS", "fourfold", "nonsynonymous", "synonymous", "WS", "SW", "GC3"}
    unknown = set(feature_intervals) - permitted
    if unknown:
        raise PilotError(f"unregistered annotation partitions: {sorted(unknown)}")
    callable_rows = merge_intervals(callable_intervals)
    result: list[dict[str, object]] = []
    weak, strong = set("AT"), set("GC")
    for name in sorted(permitted):
        rows = intersect_intervals(callable_rows, feature_intervals.get(name, []))
        denominator = interval_bp(rows)
        selected = [variant for variant in variants
                    if _contained(_position_set(rows), variant.contig, variant.pos0,
                                  variant.pos0 + len(variant.ref))]
        if name == "WS":
            selected = [v for v in selected if len(v.ref) == len(v.alt) == 1
                        and v.ref in weak and v.alt in strong]
        elif name == "SW":
            selected = [v for v in selected if len(v.ref) == len(v.alt) == 1
                        and v.ref in strong and v.alt in weak]
        result.append({
            "partition": name,
            "callable_h1_bp": denominator,
            "variant_records": len(selected),
            "diversity_per_callable_h1_bp": len(selected) / denominator if denominator else None,
        })
    return result


def confidence_tier(evidence: Mapping[str, object]) -> str:
    """Classify confidence without turning validation covariates into vetoes.

    False values for scientific execution/result invariants make a result
    unusable (``X``).  Missing result evidence is pending/low confidence
    (``C``).  Assembly and read-validation fields distinguish confidence
    tiers only; they never authorize or refuse core execution.
    """
    hard = ("exact_pair", "accepted_input_digests", "mutually_comparable_assemblies",
            "mapping_1to1_pass", "ref_reconstruction_pass", "h2_reconstruction_pass",
            "mask_accounting_pass", "callability_pass", "consensus_pass",
            "reproducibility_pass")
    if any(evidence.get(key) is False for key in hard):
        return "X"
    if any(evidence.get(key) is None for key in hard):
        return "C"
    covariates = ("qv_pass", "completeness_pass", "collapse_pass", "repeat_audit",
                  "exact_read_chemistry", "raw_read_validation", "kmer_validation",
                  "copy_number_validation", "published_estimate_validation",
                  "long_range_switch_validation", "independent_ne_validation")
    observed = [evidence.get(key) for key in covariates]
    if all(value is True for value in observed):
        return "A"
    if any(value is None for value in observed):
        return "C"
    return "B"


def estimate_resources(
    h1_bytes: int, h2_bytes: int, h1_bp: int, h2_bp: int, contigs: int,
    native_partitions: int = 0, calibration: Mapping[str, float] | None = None,
    threads: int = 8,
) -> dict[str, object]:
    values = (h1_bytes, h2_bytes, h1_bp, h2_bp, contigs, native_partitions)
    if any(value < 0 for value in values) or h1_bp == 0 or h2_bp == 0 or threads <= 0:
        raise PilotError("resource estimation requires measured positive sequence lengths")
    c = {
        "map_cpu_hours_per_gbp": 4.0, "map_rss_gib_per_gbp": 5.0,
        "scratch_bytes_per_input_byte": 4.0, "partition_cpu_hours_per_1000": 0.5,
        "read_bytes_per_input_byte": 5.0, "write_bytes_per_input_byte": 2.0,
        "metadata_operations_per_contig": 8.0,
    }
    if calibration:
        unknown = set(calibration) - set(c)
        if unknown:
            raise PilotError(f"unknown calibration coefficients: {sorted(unknown)}")
        c.update({key: float(value) for key, value in calibration.items()})
    gbp = (h1_bp + h2_bp) / 1e9
    input_bytes = h1_bytes + h2_bytes
    map_cpu = gbp * c["map_cpu_hours_per_gbp"]
    partition_cpu = native_partitions / 1000 * c["partition_cpu_hours_per_1000"]
    return {
        "basis": {"h1_bytes": h1_bytes, "h2_bytes": h2_bytes, "h1_bp": h1_bp,
                  "h2_bp": h2_bp, "contigs": contigs, "native_partitions": native_partitions},
        "coefficients": c,
        "threads": threads,
        "map_cpu_hours_estimate": map_cpu,
        "map_wall_hours_estimate": map_cpu / threads,
        "map_peak_rss_gib_estimate": max(2.0, gbp * c["map_rss_gib_per_gbp"]),
        "scratch_bytes_estimate": int(input_bytes * c["scratch_bytes_per_input_byte"]),
        "partition_cpu_hours_estimate": partition_cpu,
        "partition_wall_hours_estimate": partition_cpu / threads,
        "read_bytes_estimate": int(input_bytes * c["read_bytes_per_input_byte"]),
        "write_bytes_estimate": int(input_bytes * c["write_bytes_per_input_byte"]),
        "metadata_operations_estimate": int(contigs * c["metadata_operations_per_contig"] + native_partitions),
        "scheduler_limits": "set per job from approved high estimates; no global memory or byte ceiling",
    }


def resource_prediction_errors(
    predicted: Mapping[str, float], observed: Mapping[str, float]
) -> dict[str, float]:
    dimensions = ("wall_time", "cpu_hours", "peak_rss", "scratch_high_water", "read_bytes", "write_bytes")
    result: dict[str, float] = {}
    for dimension in dimensions:
        expected = float(predicted.get(dimension, 0))
        actual = float(observed.get(dimension, 0))
        if expected <= 0 or actual <= 0:
            raise PilotError(f"resource calibration requires positive {dimension}")
        result[dimension] = abs(expected - actual) / actual
    return result


def load_primary_pair(selection_id: str, manifest: Path | str = PRIMARY_MANIFEST) -> dict[str, str]:
    with Path(manifest).open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    matches = [row for row in rows if row["selection_id"] == selection_id]
    if len(matches) != 1:
        raise PilotError(f"selection_id must resolve exactly once: {selection_id}")
    return matches[0]


def validate_pair_input_manifest(
    value: Mapping[str, object], manifest: Path | str = PRIMARY_MANIFEST,
    verify_files: bool = True,
) -> dict[str, object]:
    """Bind acquired files to one immutable design row.

    The authorization boundary is intentionally narrow: exact pair provenance,
    accepted digests, and readable mutually comparable assemblies.  QV, BUSCO,
    repeat reports, copy-number/k-mer audits, chemistry, raw-read validation,
    annotation, and independent Ne evidence are retained as confidence
    covariates.  They are not universal pre-execution gates.
    """
    selection_id = str(value.get("selection_id", ""))
    pair = load_primary_pair(selection_id, manifest)
    bindings = (
        ("biosample", "biosample"), ("individual_or_isolate", "individual_or_isolate"),
        ("h1_accession_version", "h1_accession_version"),
        ("h2_accession_version", "h2_accession_version"),
    )
    for input_key, design_key in bindings:
        if value.get(input_key) != pair[design_key]:
            raise PilotError(f"pair provenance mismatch for {input_key}")
    if value.get("orientation") != "H1_reference_H2_query":
        raise PilotError("input manifest must orient exact H1 as reference and H2 as query")
    assets = value.get("assets")
    if not isinstance(assets, Mapping) or set(assets) < {"h1_fasta", "h2_fasta"}:
        raise PilotError("input manifest requires exact H1 and H2 FASTA assets")
    dictionaries: dict[str, list[dict[str, object]]] = {}
    for role in ("h1_fasta", "h2_fasta"):
        asset = assets[role]
        if not isinstance(asset, Mapping):
            raise PilotError(f"invalid {role} asset record")
        if len(str(asset.get("sha256", ""))) != 64 or int(asset.get("size_bytes", 0)) <= 0:
            raise PilotError(f"{role} lacks frozen checksum and positive byte size")
        path = Path(str(asset.get("path", "")))
        if verify_files:
            if not path.is_file() or path.stat().st_size != int(asset["size_bytes"]):
                raise PilotError(f"{role} file/size mismatch")
            if sha256_file(path) != asset["sha256"]:
                raise PilotError(f"{role} SHA256 mismatch")
            dictionaries[role] = sequence_dictionary(parse_fasta(path))
            if asset.get("sequence_dictionary") != dictionaries[role]:
                raise PilotError(f"{role} sequence dictionary mismatch")
    if verify_files:
        h1_bp = sum(int(row["length"]) for row in dictionaries["h1_fasta"])
        h2_bp = sum(int(row["length"]) for row in dictionaries["h2_fasta"])
        ratio = h1_bp / h2_bp
        if not (0.5 <= ratio <= 2.0):
            raise PilotError("H1/H2 assemblies are not mutually comparable in span")
    else:
        h1_bp = h2_bp = None
        ratio = None
    covariates = value.get("confidence_covariates", value.get("core_qc", {}))
    if covariates is None:
        covariates = {}
    if not isinstance(covariates, Mapping):
        raise PilotError("confidence covariates must be an object when supplied")
    selective = value.get("selective_validation", {}) or {}
    if not isinstance(selective, Mapping):
        raise PilotError("selective validation evidence must be an object when supplied")
    return {
        "selection_id": selection_id,
        "pair_design_sha256": sha256_file(manifest),
        "exact_pair": True,
        "orientation": "H1_reference_H2_query",
        "accepted_input_digests": True,
        "mutually_comparable_assemblies": True,
        "h1_sequence_bp": h1_bp,
        "h2_sequence_bp": h2_bp,
        "h1_h2_span_ratio": ratio,
        "confidence_covariates_are_authorization_gates": False,
        "confidence_covariates": dict(covariates),
        "selective_validation_is_universal_core_gate": False,
        "tracked_selective_validation": dict(selective),
        "dictionaries": dictionaries,
    }


def verify_environment_lock(lock_path: Path | str = ENVIRONMENT_LOCK, require_realized: bool = True) -> dict[str, object]:
    lock = json.loads(Path(lock_path).read_text(encoding="utf-8"))
    for relative, expected in lock["input_sha256"].items():
        observed = sha256_file(ROOT / relative)
        if observed != expected:
            raise PilotError(f"environment input digest drift: {relative}")
    if require_realized:
        required = ("profile", "derivation", "closure_sha256", "executables")
        if any(not lock["realization"].get(field) for field in required):
            raise PilotError("Guix environment has not been captured; ambient tools are forbidden")
        for executable in lock["realization"]["executables"]:
            path = Path(executable["path"])
            if not path.is_file() or sha256_file(path) != executable["sha256"]:
                raise PilotError(f"Guix executable identity mismatch: {path}")
    return lock


def verify_environment_capture(
    capture_path: Path | str, lock_path: Path | str = ENVIRONMENT_LOCK
) -> dict[str, object]:
    lock = verify_environment_lock(lock_path, require_realized=False)
    capture = json.loads(Path(capture_path).read_text(encoding="utf-8"))
    if capture.get("channel_commit") != lock["channel_commit"]:
        raise PilotError("captured Guix channel commit differs from lock")
    expected_manifest = lock["input_sha256"]["analysis/guix/vgp_10_pilot/manifest.scm"]
    if capture.get("manifest_sha256") != expected_manifest:
        raise PilotError("captured Guix manifest digest differs from lock")
    if capture.get("source_identities") != lock["source_identities"]:
        raise PilotError("captured source/submodule/companion identities differ from lock")
    for field in ("profile", "derivation", "closure_sha256", "executables", "reproducibility"):
        if not capture.get(field):
            raise PilotError(f"environment capture lacks {field}")
    if not str(capture["profile"]).startswith("/gnu/store/") or not Path(capture["profile"]).exists():
        raise PilotError("captured Guix profile is unavailable")
    if not str(capture["derivation"]).startswith("/gnu/store/") or not Path(capture["derivation"]).exists():
        raise PilotError("captured Guix profile derivation is unavailable")
    for executable in capture["executables"]:
        path = Path(executable["path"])
        if not path.is_file() or sha256_file(path) != executable["sha256"]:
            raise PilotError(f"captured executable identity mismatch: {path}")
    if capture["reproducibility"].get("psmc_guix_check") != "passed":
        raise PilotError("PSMC Guix rebuild check did not pass")
    return capture


def atomic_promote(partial: Path, final: Path, sentinel_payload: Mapping[str, object]) -> None:
    if final.exists():
        marker = final / ".complete.json"
        if not marker.is_file():
            raise PilotError(f"final directory exists without completion sentinel: {final}")
        return
    if not partial.is_dir():
        raise PilotError(f"partial stage directory is absent: {partial}")
    payload = dict(sentinel_payload)
    payload["files"] = {
        str(path.relative_to(partial)): sha256_file(path)
        for path in sorted(partial.rglob("*")) if path.is_file()
    }
    (partial / ".complete.json").write_text(canonical_json(payload), encoding="utf-8")
    final.parent.mkdir(parents=True, exist_ok=True)
    os.replace(partial, final)


def materialize_mask_consensus_psmc(
    h1_fasta: Path | str, h2_fasta: Path | str, normalized_vcf: Path | str,
    universe_bed: Path | str, exclusion_beds: Mapping[str, Path | str], output_dir: Path | str,
    contig_map: Mapping[str, str] | None = None, selection_id: str = "fixture",
    design_sha256: str | None = None, attempts: int = 200,
    aligned_regions: Sequence[Mapping[str, object]] | None = None,
) -> dict[str, object]:
    """Materialize the deterministic mask/consensus/bootstrap join packet."""
    output = Path(output_dir)
    if output.exists() and any((output / name).exists() for name in ("masks", "consensus")):
        raise PilotError(f"join output already contains a mask or consensus: {output}")
    output.mkdir(parents=True, exist_ok=True)
    h1, h2 = parse_fasta(h1_fasta), parse_fasta(h2_fasta)
    variants = parse_vcf(normalized_vcf)
    concordance = validate_ref_and_reconstruct_h2(
        h1, h2, variants, contig_map, aligned_regions=aligned_regions
    )
    flags = {reason: read_bed(path) for reason, path in exclusion_beds.items()}
    callable_rows, exclusions, reconciliation = construct_reason_mask(read_bed(universe_bed), flags)
    masks = output / "masks"
    masks.mkdir()
    write_bed(masks / "callable.bed", callable_rows)
    for reason in REASON_ORDER:
        write_bed(masks / f"exclusions.{reason}.bed", exclusions[reason])
    (masks / "mask_reconciliation.json").write_text(canonical_json(reconciliation), encoding="utf-8")
    consensus, consensus_qc = build_diploid_consensus(h1, callable_rows, variants, indel_flank=10)
    consensus_dir = output / "consensus"
    consensus_dir.mkdir()
    write_fasta(consensus_dir / "consensus.fa", consensus)
    (consensus_dir / "input.psmcfa").write_text(consensus_to_psmcfa(consensus), encoding="utf-8")
    bootstrap_dir = consensus_dir / "bootstrap"
    bootstrap_dir.mkdir()
    design_sha256 = design_sha256 or sha256_file(DESIGN_MANIFEST)
    units = freeze_bootstrap_units(callable_rows, 5_000_000)
    manifest = bootstrap_manifest(callable_rows, selection_id, design_sha256, attempts=attempts)
    write_bed(consensus_dir / "bootstrap_units.5mb.bed", units)
    with (consensus_dir / "bootstrap_manifest.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("replicate", "block_bp", "seed", "unit_count", "sampled_unit_indices"))
        for row in manifest:
            writer.writerow((row["replicate"], row["block_bp"], row["seed"], row["unit_count"],
                             ",".join(map(str, row["sampled_unit_indices"]))))
    (bootstrap_dir / "README.txt").write_text(
        "Replicate PSMCFA is generated atomically in job-local scratch from the frozen unit BED "
        "and sampled indices; it is not duplicated in durable storage.\n", encoding="utf-8"
    )
    # Freeze both block-length sensitivities, but defer their expensive runs.
    for block_bp in (1_000_000, 10_000_000):
        write_bed(consensus_dir / f"bootstrap_units.{block_bp // 1_000_000}mb.bed",
                  freeze_bootstrap_units(callable_rows, block_bp))
    qc = {
        "concordance": concordance,
        "mask": reconciliation,
        "consensus": consensus_qc,
        "bootstrap_attempts": attempts,
        "primary_block_bp": 5_000_000,
        "sensitivity_block_bp": [1_000_000, 10_000_000],
        "blocks_cross_contigs": False,
        "blocks_cross_mask_discontinuities": False,
    }
    (output / "join_qc.json").write_text(canonical_json(qc), encoding="utf-8")
    return qc


def _print_commands(selection_id: str, work_dir: str, threads: int) -> dict[str, object]:
    pair = load_primary_pair(selection_id)
    base = Path(work_dir)
    h1 = str(base / "inputs/h1.fa")
    h2 = str(base / "inputs/h2.fa")
    paf = str(base / "mapping/h2_to_h1.1to1.paf")
    impg_dir = base / "impg"
    stages = {
        "mapping": sweepga_command("sweepga", h1, h2, paf, threads),
        "impg": impg_commands(
            "impg", paf, str(impg_dir / "h1_h2.impg"), str(impg_dir / "partitions"),
            str(impg_dir / "focus.native.bed"), str(impg_dir / "calls"),
            str(impg_dir / "vcf.list"), str(impg_dir / "laced.vcf"), h1, h2, threads,
        ),
        "normalize": bcftools_commands(
            "bcftools", str(impg_dir / "laced.vcf"), h1,
            str(base / "mapping/h1.1to1.bed"), str(base / "variants/normalized"),
        ),
    }
    return {
        "selection_id": selection_id,
        "h1_accession_version": pair["h1_accession_version"],
        "h2_accession_version": pair["h2_accession_version"],
        "orientation": "H1_reference_H2_query",
        "commands": stages,
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    plan = sub.add_parser("plan", help="render the exact external-tool stage order")
    plan.add_argument("selection_id")
    plan.add_argument("work_dir")
    plan.add_argument("--threads", type=int, default=8)
    audit = sub.add_parser("audit-paf", help="audit H2-query/H1-target SweepGA PAF multiplicity")
    audit.add_argument("paf")
    audit.add_argument("h1_fasta")
    audit.add_argument("h2_fasta")
    env = sub.add_parser("verify-environment", help="fail closed on uncaptured Guix identity")
    env.add_argument("--allow-unrealized", action="store_true")
    capture = sub.add_parser("verify-capture", help="verify a realized Guix capture against the lock")
    capture.add_argument("capture_json")
    resource = sub.add_parser("estimate-resources", help="estimate from exact input sizes and calibration")
    resource.add_argument("h1_fasta")
    resource.add_argument("h2_fasta")
    resource.add_argument("--native-partitions", type=int, default=0)
    resource.add_argument("--threads", type=int, default=8)
    resource.add_argument("--calibration-json")
    emit = sub.add_parser("emit-bootstrap", help="materialize one frozen bootstrap in scratch")
    emit.add_argument("consensus_fasta")
    emit.add_argument("units_bed")
    emit.add_argument("manifest_tsv")
    emit.add_argument("replicate", type=int)
    emit.add_argument("output")
    args = parser.parse_args(argv)
    try:
        if args.command == "plan":
            result = _print_commands(args.selection_id, args.work_dir, args.threads)
        elif args.command == "audit-paf":
            result = audit_sweepga_paf(
                args.paf, set(parse_fasta(args.h1_fasta)), set(parse_fasta(args.h2_fasta))
            )
        elif args.command == "verify-environment":
            result = verify_environment_lock(require_realized=not args.allow_unrealized)
        elif args.command == "verify-capture":
            result = verify_environment_capture(args.capture_json)
        elif args.command == "estimate-resources":
            h1, h2 = parse_fasta(args.h1_fasta), parse_fasta(args.h2_fasta)
            calibration = (json.loads(Path(args.calibration_json).read_text())
                           if args.calibration_json else None)
            result = estimate_resources(
                Path(args.h1_fasta).stat().st_size, Path(args.h2_fasta).stat().st_size,
                sum(map(len, h1.values())), sum(map(len, h2.values())), len(h1) + len(h2),
                args.native_partitions, calibration, args.threads,
            )
        else:
            sequences = parse_fasta(args.consensus_fasta)
            units = read_bed(args.units_bed)
            with Path(args.manifest_tsv).open(newline="", encoding="utf-8") as handle:
                matches = [row for row in csv.DictReader(handle, delimiter="\t")
                           if int(row["replicate"]) == args.replicate]
            if len(matches) != 1:
                raise PilotError("bootstrap replicate does not resolve exactly once")
            sampled = [int(value) for value in matches[0]["sampled_unit_indices"].split(",")]
            Path(args.output).write_text(bootstrap_psmcfa(sequences, units, sampled), encoding="utf-8")
            result = {"replicate": args.replicate, "output": args.output,
                      "sha256": sha256_file(args.output), "boundary_aware": True}
    except (PilotError, OSError, json.JSONDecodeError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
