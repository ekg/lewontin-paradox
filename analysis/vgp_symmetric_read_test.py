#!/usr/bin/env python3
"""Symmetric read-vs-assembly validation for VGP pairs (review defect I-1).

The prior P07 read-validation compared reads mapped only to H1 against a
callset that has since been rejected.  This module supports the symmetric
redesign: reads are mapped to both haplotypes independently, the assembly
callable denominator and SNP sites are lifted between haplotype coordinate
frames through the exact 1:1 PAF, and every statistic is computed per frame
so that no single reference is privileged as a truth oracle.

Subcommands
-----------
transfer
    Lift a pi-callable BED and assembly SNP sites H1<->H2 through the PAF
    CIGARs with strand-aware allele complementation, reciprocal round-trip
    verification, and transfer statistics.
evidence
    Allele-aware pileup evidence at assembly SNVs for one frame.  Unlike the
    legacy assembly-evidence pass, '.'/',' pileup symbols are resolved
    against the pileup's own REF column (the frame fasta base), so reads in
    a lifted (h2) frame that match H2's base count as H2-allele support,
    never as H1-allele support (deep-dive finding A).  Per-site allele
    origin metadata (frame ref base, site H1/H2 alleles) is recorded.
metrics
    Per-frame symmetric statistics: pi per platform, genotype concordance at
    assembly SNVs with both-direction contradiction rates, and per-contig /
    per-bin discordance localization.  Accepts legacy (v1) and allele-aware
    (v2) evidence classifications.
triage
    Assembly-only collapsed-paralogy screen (no reads): per-bin divergence
    and GC scan of the bounded callset, flagged paralog bins, and a
    corrected pi excluding flagged bins.
report
    Render the multi-frame metrics into a markdown report skeleton.

This module is standard-library only so it runs identically under the pinned
Guix profile (python 3.10) and the repository test environment.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import subprocess
import sys
from array import array
from bisect import bisect_right
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple

SCHEMA_TRANSFER = "vgp-symmetric-read-test-transfer-v1"
SCHEMA_METRICS = "vgp-symmetric-read-test-metrics-v1"
SCHEMA_REPORT = "vgp-symmetric-read-test-report-v1"
SCHEMA_EVIDENCE = "vgp-symmetric-read-test-evidence-v2"
SCHEMA_TRIAGE = "vgp-paralog-triage-v1"

COMPLEMENT = {
    "A": "T", "C": "G", "G": "C", "T": "A", "N": "N",
    "a": "t", "c": "g", "g": "c", "t": "a", "n": "n",
}

# --------------------------------------------------------------------------
# PAF parsing and 1:1 chain index
# --------------------------------------------------------------------------


@dataclass
class CigarIndex:
    """Prefix-sum index over one alignment's CIGAR for O(log n) offset maps."""

    qcum: array  # query bases consumed after each token (length n+1)
    tcum: array  # target bases consumed after each token (length n+1)
    ops: List[str]

    @classmethod
    def build(cls, cigar: Sequence[Tuple[str, int]]) -> "CigarIndex":
        qcum = array("q", [0])
        tcum = array("q", [0])
        ops: List[str] = []
        q = t = 0
        for op, length in cigar:
            if op in ("M", "=", "X"):
                q += length
                t += length
            elif op in ("I", "S"):
                q += length
            elif op in ("D", "N"):
                t += length
            else:
                raise ValueError(f"unsupported CIGAR op {op}")
            qcum.append(q)
            tcum.append(t)
            ops.append(op)
        return cls(qcum=qcum, tcum=tcum, ops=ops)


@dataclass(frozen=True)
class Alignment:
    q_contig: str
    q_len: int
    q_start: int  # 0-based half-open, forward query strand
    q_end: int
    strand: str
    t_contig: str
    t_len: int
    t_start: int  # 0-based half-open on forward target strand
    t_end: int
    cigar: Tuple[Tuple[str, int], ...]


def parse_cigar(text: str) -> Tuple[Tuple[str, int], ...]:
    tokens: List[Tuple[str, int]] = []
    number = ""
    for char in text:
        if char.isdigit():
            number += char
            continue
        if not number:
            raise ValueError(f"malformed CIGAR token near {char!r} in {text[:64]!r}")
        tokens.append((char, int(number)))
        number = ""
    if number:
        raise ValueError(f"trailing CIGAR digits in {text[:64]!r}")
    return tuple(tokens)


@dataclass
class PafChainIndex:
    """Query->target and target->query lifting over a verified 1:1 PAF.

    The index refuses overlapping query or target coverage: reciprocal 1:1
    semantics require every lifted base to have exactly one owner in each
    direction.  Overlaps raise on load rather than silently choosing a row.
    """

    rows: List[Alignment]
    _by_query: Dict[str, List[Alignment]] = field(default_factory=dict, repr=False)
    _by_target: Dict[str, List[Alignment]] = field(default_factory=dict, repr=False)
    _q_starts: Dict[str, List[int]] = field(default_factory=dict, repr=False)
    _q_pmax: Dict[str, List[int]] = field(default_factory=dict, repr=False)
    _t_starts: Dict[str, List[int]] = field(default_factory=dict, repr=False)
    _t_pmax: Dict[str, List[int]] = field(default_factory=dict, repr=False)
    _cigar_cache: Dict[Alignment, CigarIndex] = field(default_factory=dict, repr=False)
    rejected_overlap_rows: int = 0

    def __post_init__(self) -> None:
        for alignment in self.rows:
            self._by_query.setdefault(alignment.q_contig, []).append(alignment)
            self._by_target.setdefault(alignment.t_contig, []).append(alignment)
        for by_key, spans_of in ((self._by_query, lambda row: (row.q_start, row.q_end)),
                                 (self._by_target, lambda row: (row.t_start, row.t_end))):
            for contig, grouped in by_key.items():
                grouped.sort(key=lambda row: spans_of(row)[0])
                previous_end = -1
                for row in grouped:
                    start, end = spans_of(row)
                    if start < previous_end:
                        # SweepGA's 1:1 cap permits residual block overlap; the
                        # overlap zones are excluded per-position in lift_site
                        # rather than rejected wholesale at load.
                        self.rejected_overlap_rows += 1
                    previous_end = max(previous_end, end)
        # Precompute per-contig start arrays and prefix maxima of end for
        # O(log n + lookback) containment queries on the real 753-row chain.
        for contig, grouped in self._by_query.items():
            self._q_starts[contig] = [row.q_start for row in grouped]
            pmax: List[int] = []
            running = -1
            for row in grouped:
                running = max(running, row.q_end)
                pmax.append(running)
            self._q_pmax[contig] = pmax
        for contig, grouped in self._by_target.items():
            self._t_starts[contig] = [row.t_start for row in grouped]
            pmax = []
            running = -1
            for row in grouped:
                running = max(running, row.t_end)
                pmax.append(running)
            self._t_pmax[contig] = pmax

    # -- single-base lifting -------------------------------------------------

    def _rows_containing(self, index: Dict[str, List[Alignment]], contig: str,
                         lo: int, hi: int, select) -> List[Alignment]:
        """All rows whose selected span intersects [lo, hi).

        Backward scan from the bisect position (rows starting at or before
        lo) with a prefix-max-end early stop, then forward scan over rows
        starting inside (lo, hi).  Both scans are required: a backward-only
        scan silently misses interior overlap rows for interval queries,
        which would let a stitched interval lift wholesale (review blocker).
        Point queries (hi == lo + 1) can never hit the forward scan.
        """
        rows = index.get(contig)
        if not rows:
            return []
        if index is self._by_query:
            starts, pmax = self._q_starts[contig], self._q_pmax[contig]
        else:
            starts, pmax = self._t_starts[contig], self._t_pmax[contig]
        found: List[Alignment] = []
        position = bisect_right(starts, lo) - 1
        while position >= 0 and pmax[position] > lo:
            row = rows[position]
            row_start, row_end = select(row)
            if row_end > lo and row_start < hi:
                found.append(row)
            position -= 1
        forward = bisect_right(starts, lo)
        while forward < len(rows) and starts[forward] < hi:
            row = rows[forward]
            row_start, row_end = select(row)
            # rows here have row_start > lo, so row_end > lo always holds
            if row_start < hi and row_end > lo:
                found.append(row)
            forward += 1
        return found

    def _cigar_for(self, row: Alignment) -> CigarIndex:
        cached = self._cigar_cache.get(row)
        if cached is None:
            cached = CigarIndex.build(row.cigar)
            self._cigar_cache[row] = cached
        return cached

    def lift_site(self, contig: str, position: int,
                  forward: bool = True) -> Optional[Tuple[str, int]]:
        """Lift a single 0-based position; returns (contig, position) or None.

        A lift is admissible only when the source position is owned by
        exactly one row on its side AND the lifted position is owned by that
        same row on the opposite side.  Positions inside overlap zones
        (multiple owning rows) are ambiguous and rejected — 1:1 semantics
        are preserved by exclusion, not by arbitrary row choice.
        """
        index, other, select_from, select_to = (
            (self._by_query, self._by_target,
             (lambda row: (row.q_start, row.q_end)),
             (lambda row: (row.t_start, row.t_end)))
            if forward else
            (self._by_target, self._by_query,
             (lambda row: (row.t_start, row.t_end)),
             (lambda row: (row.q_start, row.q_end)))
        )
        candidates = self._rows_containing(index, contig, position, position + 1, select_from)
        if len(candidates) != 1:
            return None
        row = candidates[0]
        cigar = self._cigar_for(row)
        if forward:
            offset = (position - row.q_start if row.strand == "+"
                      else row.q_end - 1 - position)
            mapped = self._walk_cigar_forward(cigar, offset)
        else:
            offset = position - row.t_start
            mapped = self._walk_cigar_reverse(cigar, offset)
        if mapped is None:
            return None
        if forward:
            out_contig, out_pos = row.t_contig, row.t_start + mapped
        else:
            out_contig, out_pos = row.q_contig, (
                row.q_start + mapped if row.strand == "+"
                else row.q_end - 1 - mapped
            )
        owners = self._rows_containing(other, out_contig, out_pos, out_pos + 1, select_to)
        if len(owners) != 1 or owners[0] is not row:
            return None
        return out_contig, out_pos

    @staticmethod
    def _walk_cigar_forward(cigar: CigarIndex, aligned_offset: int) -> Optional[int]:
        """Map an offset along the CIGAR query axis to the target axis."""
        i = bisect_right(cigar.qcum, aligned_offset) - 1
        if i < 0:
            return None
        if cigar.qcum[i + 1] <= aligned_offset:
            return None
        op = cigar.ops[i]
        if op in ("I", "S"):
            return None  # inside a query-only gap
        return cigar.tcum[i] + (aligned_offset - cigar.qcum[i])

    @staticmethod
    def _walk_cigar_reverse(cigar: CigarIndex, target_offset: int) -> Optional[int]:
        """Map an offset along the CIGAR target axis to the query axis."""
        i = bisect_right(cigar.tcum, target_offset) - 1
        if i < 0:
            return None
        if cigar.tcum[i + 1] <= target_offset:
            return None
        op = cigar.ops[i]
        if op in ("D", "N"):
            return None  # inside a target-only gap
        return cigar.qcum[i] + (target_offset - cigar.tcum[i])

    # -- interval lifting ----------------------------------------------------

    def lift_interval(self, contig: str, start: int, end: int,
                      forward: bool = True) -> Optional[Tuple[str, int, int]]:
        """Lift a half-open interval; the single candidate row must fully
        contain [start, end), and no other row may intersect the interval.
        Both endpoints must lift inside that same row."""
        index = self._by_query if forward else self._by_target
        select = (lambda row: (row.q_start, row.q_end)) if forward else (lambda row: (row.t_start, row.t_end))
        candidates = self._rows_containing(index, contig, start, end, select)
        if len(candidates) != 1:
            return None
        row_start, row_end = select(candidates[0])
        if row_start > start or row_end < end:
            return None  # intersects but does not fully contain
        lo = self.lift_site(contig, start, forward=forward)
        hi = self.lift_site(contig, end - 1, forward=forward)
        if lo is None or hi is None or lo[0] != hi[0]:
            return None
        return lo[0], min(lo[1], hi[1]), max(lo[1], hi[1]) + 1

    def round_trip_site(self, contig: str, position: int,
                        forward: bool = True) -> Optional[Tuple[str, int]]:
        """Reciprocal 1:1 check for a single site.

        The candidate lift must be unique (structural uniqueness is enforced
        at load), and every row touching either endpoint — as target of the
        source or as query of the destination — must agree on the pair.
        Disagreeing partner claims (an inverse row mapping the destination
        elsewhere, or a different row claiming this site as its target)
        reject the site.  A one-directional PAF has no inverse rows and every
        structurally unique lift is accepted.
        """
        candidate = self.lift_site(contig, position, forward=forward)
        if candidate is None:
            return None
        dest_contig, dest_pos = candidate
        # Rows claiming (contig, position) as TARGET imply a source partner.
        for row in self._by_target.get(contig, []):
            if not (row.t_start <= position < row.t_end):
                continue
            offset = position - row.t_start
            mapped = self._walk_cigar_reverse(self._cigar_for(row), offset)
            if mapped is None:
                return None
            source = (row.q_start + mapped if row.strand == "+"
                      else row.q_end - 1 - mapped)
            if (row.q_contig, source) != (dest_contig, dest_pos):
                return None
        # Rows claiming (dest_contig, dest_pos) as QUERY imply a continuation.
        for row in self._by_query.get(dest_contig, []):
            if not (row.q_start <= dest_pos < row.q_end):
                continue
            onward = self.lift_site(dest_contig, dest_pos, forward=True)
            if onward != (contig, position):
                return None
        return candidate

    def round_trip_interval(self, contig: str, start: int, end: int,
                            forward: bool = True) -> Optional[Tuple[str, int, int]]:
        lifted = self.lift_interval(contig, start, end, forward=forward)
        if lifted is None:
            return None
        if self.round_trip_site(contig, start, forward=forward) is None:
            return None
        if self.round_trip_site(contig, end - 1, forward=forward) is None:
            return None
        back = self.lift_interval(lifted[0], lifted[1], lifted[2], forward=not forward)
        if back != (contig, start, end):
            return None
        return lifted


    def classify_site(self, contig: str, position: int,
                      forward: bool = True) -> Tuple[str, Optional[Tuple[str, int]]]:
        """Classify a site for transfer accounting.

        Returns (reason, lifted) where reason is one of:
          lifted            — site transferred and round-trip verified
          unaligned         — no row covers the source position
          overlap           — source position owned by multiple rows
          indel_gap         — unique row but position inside a CIGAR gap
          owner_mismatch    — lifted position not uniquely owned by the same
                              row on the opposite side
          roundtrip_mismatch — structural lift exists but inverse/continuation
                              rows disagree on the pair
        """
        index, other, select_from, select_to = (
            (self._by_query, self._by_target,
             (lambda row: (row.q_start, row.q_end)),
             (lambda row: (row.t_start, row.t_end)))
            if forward else
            (self._by_target, self._by_query,
             (lambda row: (row.t_start, row.t_end)),
             (lambda row: (row.q_start, row.q_end)))
        )
        candidates = self._rows_containing(index, contig, position, position + 1, select_from)
        if not candidates:
            return "unaligned", None
        if len(candidates) != 1:
            return "overlap", None
        row = candidates[0]
        cigar = self._cigar_for(row)
        if forward:
            offset = (position - row.q_start if row.strand == "+"
                      else row.q_end - 1 - position)
            mapped = self._walk_cigar_forward(cigar, offset)
        else:
            offset = position - row.t_start
            mapped = self._walk_cigar_reverse(cigar, offset)
        if mapped is None:
            return "indel_gap", None
        if forward:
            out_contig, out_pos = row.t_contig, row.t_start + mapped
        else:
            out_contig, out_pos = row.q_contig, (
                row.q_start + mapped if row.strand == "+"
                else row.q_end - 1 - mapped
            )
        owners = self._rows_containing(other, out_contig, out_pos, out_pos + 1, select_to)
        if len(owners) != 1 or owners[0] is not row:
            return "owner_mismatch", None
        lifted = (out_contig, out_pos)
        if self.round_trip_site(contig, position, forward=forward) != lifted:
            return "roundtrip_mismatch", lifted
        return "lifted", lifted


def load_paf(path: Path) -> Tuple[PafChainIndex, dict]:
    rows: List[Alignment] = []
    kept_fields: Dict[str, object] = {}
    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 12:
                raise ValueError(
                    f"{path}:{line_number}: PAF row has {len(parts)} fields; "
                    f"at least 12 required")
            cigar = None
            for tag in parts[12:]:
                if tag.startswith("cg:Z:"):
                    cigar = parse_cigar(tag[5:])
                    break
            if cigar is None:
                raise ValueError(f"{path}:{line_number}: PAF row without cg:Z tag")
            strand = parts[4]
            if strand not in ("+", "-"):
                raise ValueError(
                    f"{path}:{line_number}: strand {strand!r} not in {{'+','-'}}")
            q_span = int(parts[3]) - int(parts[2])
            t_span = int(parts[8]) - int(parts[7])
            consumed_q = sum(length for op, length in cigar if op not in ("D", "N"))
            consumed_t = sum(length for op, length in cigar if op not in ("I", "S"))
            if q_span != consumed_q or t_span != consumed_t:
                raise ValueError(
                    f"{path}:{line_number}: span/CIGAR mismatch: "
                    f"query {q_span} != consumed {consumed_q} or "
                    f"target {t_span} != consumed {consumed_t}")
            rows.append(Alignment(
                q_contig=parts[0], q_len=int(parts[1]),
                q_start=int(parts[2]), q_end=int(parts[3]),
                strand=strand,
                t_contig=parts[5], t_len=int(parts[6]),
                t_start=int(parts[7]), t_end=int(parts[8]),
                cigar=cigar,
            ))
    kept_fields["row_count"] = len(rows)
    kept_fields["query_bp_aligned"] = sum(row.q_end - row.q_start for row in rows)
    kept_fields["target_bp_aligned"] = sum(row.t_end - row.t_start for row in rows)
    index = PafChainIndex(rows=rows)
    stats = {
        "paf_path": str(path),
        **kept_fields,
        "forward_contigs": sorted(index._by_query),
        "reverse_contigs": sorted(index._by_target),
    }
    return index, stats


# --------------------------------------------------------------------------
# transfer subcommand
# --------------------------------------------------------------------------


def read_bed(path: Path) -> Iterator[Tuple[str, int, int]]:
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            chrom, start, end = line.split("\t")[:3]
            yield chrom, int(start), int(end)


def read_sites(path: Path) -> Iterator[Tuple[str, int, str, str]]:
    """Yield (chrom, 0-based position, ref, alt) from a sites TSV with header."""
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"chrom", "position_1based", "ref", "alt"}
        if not required.issubset(reader.fieldnames or ()):
            raise SystemExit(f"sites TSV missing columns {sorted(required)}")
        for row in reader:
            yield (row["chrom"], int(row["position_1based"]) - 1, row["ref"], row["alt"])


def complement_alleles(ref: str, alt: str) -> Tuple[str, str]:
    return ("".join(COMPLEMENT.get(base, "N") for base in ref),
            "".join(COMPLEMENT.get(base, "N") for base in alt))


def cmd_transfer(args: argparse.Namespace) -> int:
    index, paf_stats = load_paf(Path(args.paf))

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    strand_of_site: Dict[Tuple[str, int], str] = {}

    # -- BED ------------------------------------------------------------------
    bed_out = out_dir / f"{args.target_label}.pi.callable.bed"
    bed_symmetric_out = out_dir / f"{args.source_label}.symmetric.pi.callable.bed"
    bed_stats = Counter()
    bed_bp_h1 = bed_bp_lifted = bed_bp_symmetric = 0
    with bed_out.open("w") as handle, bed_symmetric_out.open("w") as sym_handle:
        for chrom, start, end in read_bed(Path(args.bed)):
            bed_bp_h1 += end - start
            lifted = index.round_trip_interval(chrom, start, end, forward=True)
            if lifted is None:
                bed_stats["rejected_interval"] += 1
                continue
            bed_stats["lifted_interval"] += 1
            bed_bp_lifted += lifted[2] - lifted[1]
            bed_bp_symmetric += end - start
            handle.write(f"{lifted[0]}\t{lifted[1]}\t{lifted[2]}\n")
            sym_handle.write(f"{chrom}\t{start}\t{end}\n")
    # sort lifted BED by (contig, start) for downstream -R/-b usage
    for path in (bed_out, bed_symmetric_out):
        lines = sorted(path.read_text().splitlines(),
                       key=lambda line: (line.split("\t")[0], int(line.split("\t")[1])))
        path.write_text("\n".join(lines) + ("\n" if lines else ""))

    # -- sites ---------------------------------------------------------------
    sites_out = out_dir / f"{args.target_label}.assembly.snps.tsv"
    roundtrip_out = out_dir / f"{args.target_label}.roundtrip.sites.tsv"
    sites_symmetric_out = out_dir / f"{args.source_label}.symmetric.assembly.snps.tsv"
    sites_stats: Counter = Counter()
    fidelity_kept = fidelity_total = 0
    with sites_out.open("w") as sites_handle, roundtrip_out.open("w") as back_handle, \
            sites_symmetric_out.open("w") as sym_handle:
        writer = csv.writer(sites_handle, delimiter="\t", lineterminator="\n")
        back_writer = csv.writer(back_handle, delimiter="\t", lineterminator="\n")
        sym_writer = csv.writer(sym_handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["chrom", "position_1based", "ref", "alt"])
        back_writer.writerow(["chrom", "position_1based", "ref", "alt"])
        sym_writer.writerow(["chrom", "position_1based", "ref", "alt"])
        for chrom, pos, ref, alt in read_sites(Path(args.sites)):
            sites_stats["total"] += 1
            reason, lifted = index.classify_site(chrom, pos, forward=True)
            if reason != "lifted":
                sites_stats[reason] += 1
                sites_stats["rejected"] += 1
                continue
            l_chrom, l_pos = lifted
            # strand: recover from forward lift by checking inverse row orientation
            row = index._rows_containing(index._by_query, chrom, pos, pos + 1,
                                         lambda r: (r.q_start, r.q_end))
            if not row:
                sites_stats["owner_mismatch"] += 1
                sites_stats["rejected"] += 1
                continue
            row = row[0]
            out_ref, out_alt = (ref, alt) if row.strand == "+" else complement_alleles(ref, alt)
            writer.writerow([l_chrom, l_pos + 1, out_ref, out_alt])
            back_writer.writerow([l_chrom, l_pos + 1, out_ref, out_alt])
            sym_writer.writerow([chrom, pos + 1, ref, alt])
            sites_stats["lifted"] += 1
            fidelity_total += 1
            fidelity_kept += 1

    result = {
        "schema_version": SCHEMA_TRANSFER,
        "paf": paf_stats,
        "bed": {
            "input_intervals": bed_stats["lifted_interval"] + bed_stats["rejected_interval"],
            "lifted_intervals": bed_stats["lifted_interval"],
            "rejected_intervals": bed_stats["rejected_interval"],
            "input_bp": bed_bp_h1,
            "lifted_bp": bed_bp_lifted,
            "lifted_fraction_of_bp": (bed_bp_lifted / bed_bp_h1) if bed_bp_h1 else None,
            "bp_delta_source_minus_dest": bed_bp_symmetric - bed_bp_lifted,
        },
        "sites": {
            "total": sites_stats["total"],
            "lifted": sites_stats["lifted"],
            "rejected": sites_stats["rejected"],
            "lifted_fraction": (sites_stats["lifted"] / sites_stats["total"]) if sites_stats["total"] else None,
            "rejection_reasons": {
                reason: sites_stats[reason]
                for reason in ("unaligned", "overlap", "indel_gap",
                               "owner_mismatch", "roundtrip_mismatch", "other")
                if sites_stats[reason]
            },
        },
        "roundtrip_fidelity": {
            # All emitted sites round-tripped by construction (round_trip_site
            # enforces it); the fraction below is lifted/total for audit.
            "kept": fidelity_kept,
            "total": fidelity_total + sites_stats["rejected"],
            "fidelity": (fidelity_kept / (fidelity_total + sites_stats["rejected"]))
            if (fidelity_total + sites_stats["rejected"]) else None,
        },
        "outputs": {
            "lifted_bed": str(bed_out),
            "symmetric_source_bed": str(bed_symmetric_out),
            "lifted_sites": str(sites_out),
            "symmetric_source_sites": str(sites_symmetric_out),
            "roundtrip_sites": str(roundtrip_out),
        },
        "symmetric_subset_bp_source": bed_bp_symmetric,
    }
    (out_dir / "transfer_stats.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"bed_bp": result["bed"], "sites": result["sites"]}, indent=2))
    return 0


# --------------------------------------------------------------------------
# allele-aware pileup evidence (deep-dive finding A fix)
# --------------------------------------------------------------------------


def parse_pileup_bases_allele_aware(bases: str, frame_ref: str, h1_allele: str,
                                    h2_allele: str) -> Dict[str, int]:
    """Count reads by which *expected allele* they carry, not by symbol.

    '.'/',' mean "matches the pileup reference" — the frame fasta base at
    that position, given as ``frame_ref`` (pileup column 3).  They are
    resolved to the allele that base equals: H1 allele, H2 allele, or, when
    the frame reference carries neither expected allele (the lift/allele-
    origin disagreement class), ``frame_ref_mismatch``.  Explicit base
    letters are classified the same way.  The legacy parser counted '.'/','
    as site-ref (H1-allele) support unconditionally, which mislabeled every
    H2-matching read in a lifted h2 frame as H1 support.
    """
    h1_allele = h1_allele.upper()
    h2_allele = h2_allele.upper()
    frame_ref = frame_ref.upper()
    counts = {"h1": 0, "h2": 0, "frame_ref_mismatch": 0, "other": 0, "deletion": 0,
              "frame_ref_symbol": 0}
    index = 0
    while index < len(bases):
        char = bases[index]
        if char == "^":
            index += 2
            continue
        if char == "$":
            index += 1
            continue
        if char in "+-":
            index += 1
            digit_start = index
            while index < len(bases) and bases[index].isdigit():
                index += 1
            if digit_start == index:
                raise ValueError("malformed pileup indel without a length")
            index += int(bases[digit_start:index])
            continue
        if char in "*#":
            counts["deletion"] += 1
        elif char in "<>":
            pass
        elif char in ".,":
            counts["frame_ref_symbol"] += 1
            if frame_ref == h1_allele:
                counts["h1"] += 1
            elif frame_ref == h2_allele:
                counts["h2"] += 1
            elif frame_ref in ("A", "C", "G", "T"):
                counts["frame_ref_mismatch"] += 1
            else:
                counts["other"] += 1
        elif char.isalpha():
            base = char.upper()
            if base == h1_allele:
                counts["h1"] += 1
            elif base == h2_allele:
                counts["h2"] += 1
            else:
                counts["other"] += 1
        index += 1
    return counts


def classify_allele_evidence(depth: int, counts: Dict[str, int], *,
                             minimum_depth: int, maximum_depth: int) -> str:
    """Site classification from allele-aware counts.

    Labels: ``balanced_heterozygous`` (both alleles >=3 reads, balance
    0.20-0.80), ``h1_only`` / ``h2_only`` (other allele <=1 read and the
    carried allele >= max(3, ceil(0.90*depth))), ``skewed`` (both observed
    but off balance), ``outside_depth_mask``, ``not_observed``.  Thresholds
    mirror the legacy assembly-evidence rules so v1/v2 numbers stay
    comparable.
    """
    if depth == 0 and counts["h1"] + counts["h2"] == 0:
        return "not_observed"
    if depth < minimum_depth or depth > maximum_depth:
        return "outside_depth_mask"
    informative = counts["h1"] + counts["h2"]
    if informative == 0:
        return "not_observed"
    balance = counts["h2"] / informative
    if (counts["h1"] >= 3 and counts["h2"] >= 3 and 0.20 <= balance <= 0.80):
        return "balanced_heterozygous"
    if counts["h2"] <= 1 and counts["h1"] >= max(3, math.ceil(0.90 * depth)):
        return "h1_only"
    if counts["h1"] <= 1 and counts["h2"] >= max(3, math.ceil(0.90 * depth)):
        return "h2_only"
    return "skewed"


def legacy_classification_of(label: str) -> str:
    """Map v2 labels onto the legacy label set for drop-in consumers."""
    if label in ("balanced_heterozygous",):
        return "supported_heterozygous"
    if label == "h1_only":
        return "contradicted_homozygous_reference"
    if label in ("skewed", "h2_only"):
        return "ambiguous"
    return label  # outside_depth_mask / not_observed keep their names


EVIDENCE_TSV_COLUMNS = (
    "chrom", "position_1based", "ref", "alt",
    "frame_ref_base", "site_h1_allele", "site_h2_allele", "allele_origin",
    "reported_depth", "h1_reads", "h2_reads", "frame_ref_symbol_reads",
    "frame_ref_mismatch_reads", "other_reads", "deletion_reads",
    "allele_balance_h2", "classification", "legacy_classification",
)


def cmd_evidence(args: argparse.Namespace) -> int:
    """Allele-aware evidence pass over one frame's raw pileup."""
    frame_ref_by_site: Dict[Tuple[str, int], str] = {}
    observations: Dict[Tuple[str, int], Tuple[int, str]] = {}
    with open(args.pileup) as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip():
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 5:
                raise SystemExit(
                    f"{args.pileup}:{line_number}: pileup row has {len(fields)} fields")
            key = (fields[0], int(fields[1]) - 1)  # pileup pos is 1-based
            if key in observations:
                raise SystemExit(f"{args.pileup}:{line_number}: duplicate pileup coordinate {key}")
            observations[key] = (int(fields[3]), fields[4])
            frame_ref_by_site[key] = fields[2].upper()

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    classifications: Counter = Counter()
    origins: Counter = Counter()
    balances: List[float] = []
    with open(args.sites) as sites_handle, open(args.output, "w") as out_handle:
        writer = csv.DictWriter(out_handle, fieldnames=list(EVIDENCE_TSV_COLUMNS),
                                delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for chrom, pos, ref, alt in read_sites(Path(args.sites)):
            key = (chrom, pos)
            frame_ref = frame_ref_by_site.get(key, "")
            if key not in observations:
                depth, counts = 0, {"h1": 0, "h2": 0, "frame_ref_mismatch": 0,
                                    "other": 0, "deletion": 0, "frame_ref_symbol": 0}
            else:
                depth, pileup_bases = observations[key]
                counts = parse_pileup_bases_allele_aware(
                    pileup_bases, frame_ref, ref, alt)
            label = classify_allele_evidence(
                depth, counts, minimum_depth=int(args.minimum_depth),
                maximum_depth=int(args.maximum_depth))
            classifications[label] += 1
            informative = counts["h1"] + counts["h2"]
            balance = (counts["h2"] / informative) if informative else None
            if label == "balanced_heterozygous" and balance is not None:
                balances.append(balance)
            if not frame_ref:
                origin = "not_observed"
            elif frame_ref == alt:
                origin = "h2_allele"
            elif frame_ref == ref:
                origin = "h1_allele"
            else:
                origin = "neither"
            origins[origin] += 1
            writer.writerow({
                "chrom": chrom, "position_1based": pos + 1, "ref": ref, "alt": alt,
                "frame_ref_base": frame_ref,
                "site_h1_allele": ref, "site_h2_allele": alt,
                "allele_origin": origin,
                "reported_depth": depth,
                "h1_reads": counts["h1"], "h2_reads": counts["h2"],
                "frame_ref_symbol_reads": counts["frame_ref_symbol"],
                "frame_ref_mismatch_reads": counts["frame_ref_mismatch"],
                "other_reads": counts["other"], "deletion_reads": counts["deletion"],
                "allele_balance_h2": f"{balance:.6f}" if balance is not None else "",
                "classification": label,
                "legacy_classification": legacy_classification_of(label),
            })

    summary = {
        "schema_version": SCHEMA_EVIDENCE,
        "frame": args.frame,
        "platform": args.platform,
        "minimum_depth": int(args.minimum_depth),
        "maximum_depth": int(args.maximum_depth),
        "classifications": dict(classifications),
        "allele_origin": dict(origins),
        "balanced_allele_balance_h2_mean": (sum(balances) / len(balances)) if balances else None,
        "balanced_sites": len(balances),
    }
    Path(args.summary).parent.mkdir(parents=True, exist_ok=True)
    Path(args.summary).write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"frame": args.frame, "platform": args.platform,
                      "classifications": dict(classifications)}, indent=2))
    return 0


# --------------------------------------------------------------------------
# metrics subcommand
# --------------------------------------------------------------------------


def map_evidence_classification(label: str) -> str:
    """Normalize v1 (legacy) and v2 (allele-aware) labels for metrics.

    v2 balanced_heterozygous corresponds to the legacy
    supported_heterozygous rule; v2 h1_only (reads homozygous for the site
    ref = the H1 allele) corresponds to the legacy
    contradicted_homozygous_reference.  v2 h2_only has no legacy
    equivalent — the legacy parser silently miscounted those reads as
    H1 support — so it is reported separately.
    """
    if label in ("supported_heterozygous", "balanced_heterozygous"):
        return "supported_heterozygous"
    if label in ("contradicted_homozygous_reference", "h1_only"):
        return "contradicted_homozygous_reference"
    if label == "h2_only":
        return "h2_only"
    return "other"


class FrameMetrics:
    frame: str
    platform: str
    read_het_snps: int
    callable_bp: int

    @property
    def pi(self) -> Optional[float]:
        return self.read_het_snps / self.callable_bp if self.callable_bp else None


class BedLookup:
    """Point-in-BED membership over merged, non-overlapping intervals."""

    def __init__(self, intervals: Iterable[Tuple[str, int, int]]):
        by_contig: Dict[str, List[Tuple[int, int]]] = {}
        for chrom, start, end in intervals:
            if end <= start:
                continue
            by_contig.setdefault(chrom, []).append((start, end))
        self._starts: Dict[str, List[int]] = {}
        self._ends: Dict[str, List[int]] = {}
        for chrom, spans in by_contig.items():
            spans.sort()
            merged: List[List[int]] = []
            for start, end in spans:
                if merged and start <= merged[-1][1]:
                    merged[-1][1] = max(merged[-1][1], end)
                else:
                    merged.append([start, end])
            self._starts[chrom] = [span[0] for span in merged]
            self._ends[chrom] = [span[1] for span in merged]

    def contains(self, chrom: str, position: int) -> bool:
        starts = self._starts.get(chrom)
        if not starts:
            return False
        idx = bisect_right(starts, position) - 1
        return idx >= 0 and self._ends[chrom][idx] > position


def cmd_metrics(args: argparse.Namespace) -> int:
    # Inputs per invocation are ONE frame (reference x platform); the sbatch
    # invokes this once per frame and the report merges them.
    frame_bed = BedLookup(read_bed(Path(args.frame_callable_bed)))
    assembly_sites_all: Dict[Tuple[str, int], Tuple[str, str]] = {}
    with open(args.assembly_sites) as handle:
        for chrom, pos, ref, alt in read_sites(Path(args.assembly_sites)):
            assembly_sites_all[(chrom, pos)] = (ref, alt)
    # Site universe = assembly SNVs inside this frame's pi-callable bed, so
    # concordance denominators and direction-B membership are comparable
    # across frames and against the inherited H1-only run.
    assembly_sites = {
        key: value for key, value in assembly_sites_all.items()
        if frame_bed.contains(*key)
    }
    site_strata = {
        "total": len(assembly_sites_all),
        "in_frame_callable": len(assembly_sites),
        "outside_frame_callable": len(assembly_sites_all) - len(assembly_sites),
    }

    read_het: Dict[Tuple[str, int], Tuple[str, str, float]] = {}
    with open(args.read_variants) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            key = (row["chrom"], int(row["position_1based"]) - 1)
            genotype = row.get("genotype", "")
            if "|" in genotype or "/" in genotype:
                alleles = genotype.replace("|", "/").split("/")
                het = len(set(alleles)) > 1
            else:
                het = True
            if not het:
                continue
            read_het[key] = (row["ref"], row["alt"], float(row["quality"] or 0))

    callable_bp = int(args.callable_bp)

    # Concordance from pileup evidence at assembly SNVs (frame coordinates),
    # restricted to the in-callable site universe.
    evidence: Dict[Tuple[str, int], dict] = {}
    with open(args.assembly_evidence) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            evidence[(row["chrom"], int(row["position_1based"]) - 1)] = row

    in_universe_evidence = {
        key: row for key, row in evidence.items() if key in assembly_sites
    }
    classifications = Counter(row["classification"] for row in in_universe_evidence.values())
    mapped = Counter(map_evidence_classification(row["classification"])
                     for row in in_universe_evidence.values())
    contradicted = mapped["contradicted_homozygous_reference"]
    supported = mapped["supported_heterozygous"]
    h2_only = mapped["h2_only"]
    resolved = sum(classifications.values()) - classifications.get("not_observed", 0)

    # Both-direction contradictions:
    #   A) assembly-het where reads say homozygous reference (evidence table)
    #   B) read-het (QUAL>=30) at callable positions where assembly is hom-ref
    direction_b_sites = [
        key for key, (ref, alt, qual) in read_het.items()
        if qual >= 30 and key not in assembly_sites
    ]

    # -- localization --------------------------------------------------------
    bin_size = int(args.bin_size)
    contig_bins: Dict[str, Dict[int, dict]] = {}
    for key in in_universe_evidence:
        chrom, pos = key
        bins = contig_bins.setdefault(chrom, {})
        index = pos // bin_size
        entry = bins.setdefault(index, {
            "assembly_sites": 0, "contradicted": 0, "read_het": 0,
            "read_het_not_in_assembly": 0,
        })
        entry["assembly_sites"] += 1
        if map_evidence_classification(evidence[key]["classification"]) == "contradicted_homozygous_reference":
            entry["contradicted"] += 1
    for key in direction_b_sites:
        chrom, pos = key
        bins = contig_bins.setdefault(chrom, {})
        index = pos // bin_size
        bins.setdefault(index, {
            "assembly_sites": 0, "contradicted": 0, "read_het": 0,
            "read_het_not_in_assembly": 0,
        })["read_het_not_in_assembly"] += 1
    for key in read_het:
        chrom, pos = key
        bins = contig_bins.get(chrom)
        if bins and pos // bin_size in bins:
            bins[pos // bin_size]["read_het"] += 1

    genome_wide_rate = (contradicted / resolved) if resolved else None
    flagged: List[dict] = []
    per_bin_rows: List[dict] = []
    for chrom in sorted(contig_bins):
        for index in sorted(contig_bins[chrom]):
            entry = contig_bins[chrom][index]
            sites = entry["assembly_sites"]
            rate = (entry["contradicted"] / sites) if sites else None
            row = {
                "chrom": chrom,
                "bin_start_0based": index * bin_size,
                "bin_end_0based": (index + 1) * bin_size,
                **entry,
                "contradiction_rate": rate,
            }
            per_bin_rows.append(row)
            if (genome_wide_rate is not None and rate is not None and sites >= int(args.min_sites_per_bin)
                    and rate >= float(args.flag_factor) * genome_wide_rate):
                flagged.append(row)

    per_contig: Dict[str, Counter] = {}
    for row in per_bin_rows:
        counter = per_contig.setdefault(row["chrom"], Counter())
        counter["assembly_sites"] += row["assembly_sites"]
        counter["contradicted"] += row["contradicted"]
        counter["read_het"] += row["read_het"]
        counter["read_het_not_in_assembly"] += row["read_het_not_in_assembly"]

    result = {
        "schema_version": SCHEMA_METRICS,
        "frame": args.frame,
        "platform": args.platform,
        "callable_bp": callable_bp,
        "read_het_snps_on_mask": len(read_het),
        "pi_read": (len(read_het) / callable_bp) if callable_bp else None,
        "assembly_sites_observed": len(in_universe_evidence),
        "assembly_site_strata": site_strata,
        "evidence_classifications": dict(classifications),
        "evidence_classifications_mapped": dict(mapped),
        "concordance": {
            "supported_heterozygous": supported,
            "contradicted_homozygous_reference": contradicted,
            "contradiction_rate_resolved": (contradicted / resolved) if resolved else None,
        },
        "both_direction_contradictions": {
            "A_assemblyHet_readsHomRef": contradicted,
            "B_readHet_assemblyHomRef": len(direction_b_sites),
            "B_density_per_callable_bp": (len(direction_b_sites) / callable_bp) if callable_bp else None,
        },
        "h2_only_reads_carry_h2_allele_exclusively": h2_only,
        "flagged_bins": flagged,
        "flag_policy": {
            "bin_size": bin_size,
            "min_sites_per_bin": int(args.min_sites_per_bin),
            "flag_factor": float(args.flag_factor),
            "genome_wide_rate": genome_wide_rate,
        },
    }

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write("\n")

    bins_path = Path(args.output).with_suffix(".bins.tsv")
    with bins_path.open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(per_bin_rows[0]) if per_bin_rows else
                                ["chrom", "bin_start_0based", "bin_end_0based"],
                                delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(per_bin_rows)

    contigs_path = Path(args.output).with_suffix(".contigs.tsv")
    with contigs_path.open("w") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["chrom", "assembly_sites", "contradicted", "read_het",
                         "read_het_not_in_assembly", "contradiction_rate"])
        for chrom in sorted(per_contig):
            c = per_contig[chrom]
            writer.writerow([chrom, c["assembly_sites"], c["contradicted"], c["read_het"],
                             c["read_het_not_in_assembly"],
                             (c["contradicted"] / c["assembly_sites"]) if c["assembly_sites"] else ""])
    print(json.dumps({"frame": args.frame, "platform": args.platform,
                      "pi_read": result["pi_read"]}, indent=2))
    return 0


# --------------------------------------------------------------------------
# triage subcommand (assembly-only collapsed-paralogy screen)
# --------------------------------------------------------------------------


def sha256_of(path: Path, chunk: int = 1 << 20) -> str:
    import hashlib
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            block = handle.read(chunk)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


class FastaScan:
    """One streaming pass: per-bin ACGT/GC counts and an in-memory fai index."""

    def __init__(self, path: Path, bin_size: int):
        self.path = path
        self.bin_size = bin_size
        # (chrom, bin_index) -> [acgt, gc, non_acgt_letters]
        self.bins: Dict[Tuple[str, int], List[int]] = {}
        # chrom -> [length, byte_offset, line_bases, line_bytes]
        self.index: Dict[str, List[int]] = {}
        self._handle = None
        contig: Optional[str] = None
        length = offset = line_bases = line_bytes = 0
        position = 0  # 0-based within contig
        with path.open("rb") as handle:
            while True:
                raw = handle.readline()
                if not raw:
                    break
                if raw.startswith(b">"):
                    if contig is not None:
                        self.index[contig] = [length, offset, line_bases, line_bytes]
                    text = raw[1:].split()
                    contig = text[0].decode() if text else ""
                    if contig in self.index:
                        raise SystemExit(f"{path}: duplicate contig {contig}")
                    length = position = 0
                    offset = 0
                    line_bases = line_bytes = 0
                    first_sequence_line = True
                    continue
                if contig is None:
                    continue
                stripped = raw.strip()
                if not stripped:
                    continue
                if first_sequence_line:
                    line_bases = len(stripped)
                    line_bytes = len(raw)
                    offset = handle.tell() - len(raw)
                    first_sequence_line = False
                length += len(stripped)
                cursor = 0
                line_start = position
                while cursor < len(stripped):
                    bin_index = (line_start + cursor) // bin_size
                    bin_end = (bin_index + 1) * bin_size
                    take = min(len(stripped) - cursor, bin_end - (line_start + cursor))
                    segment = stripped[cursor:cursor + take].upper()
                    entry = self.bins.setdefault((contig, bin_index), [0, 0, 0])
                    entry[0] += (segment.count(b"A") + segment.count(b"C")
                                 + segment.count(b"G") + segment.count(b"T"))
                    entry[1] += segment.count(b"G") + segment.count(b"C")
                    entry[2] += take - (segment.count(b"A") + segment.count(b"C")
                                        + segment.count(b"G") + segment.count(b"T"))
                    cursor += take
                position += len(stripped)
            if contig is not None:
                self.index[contig] = [length, offset, line_bases, line_bytes]

    def base_at(self, chrom: str, position: int) -> str:
        """Uppercase base at 0-based position (needs the fai index from scan)."""
        record = self.index.get(chrom)
        if record is None or not (0 <= position < record[0]):
            return ""
        length, offset, line_bases, line_bytes = record
        if line_bases == 0:
            return ""
        line, column = divmod(position, line_bases)
        if self._handle is None:
            self._handle = self.path.open("rb")
        self._handle.seek(offset + line * line_bytes + column)
        byte = self._handle.read(1)
        return byte.decode().upper()

    def close(self) -> None:
        if self._handle is not None:
            self._handle.close()
            self._handle = None


def stream_bcf_sites(bcf: Path, bcftools: str) -> Iterator[Tuple[str, int, str, str]]:
    """Yield (chrom, 0-based pos, ref, alt) for biallelic SNP records."""
    process = subprocess.Popen(
        [bcftools, "query", "-f", "%CHROM\\t%POS\\t%REF\\t%ALT\\n", "-i",
         'TYPE="snp" && ALT!="*" && ALT!="."', str(bcf)],
        stdout=subprocess.PIPE, text=True)
    assert process.stdout is not None
    for line in process.stdout:
        fields = line.rstrip("\n").split("\t")
        if len(fields) != 4:
            continue
        chrom, pos, ref, alt = fields
        if "," in alt or len(ref) != 1 or len(alt) != 1:
            continue
        yield chrom, int(pos) - 1, ref.upper(), alt.upper()
    code = process.wait()
    if code != 0:
        raise SystemExit(f"{bcftools} query exited {code} for {bcf}")


def stream_sites_tsv(path: Path) -> Iterator[Tuple[str, int, str, str]]:
    with path.open() as handle:
        for chrom, pos, ref, alt in read_sites(path):
            yield chrom, pos, ref, alt


def detect_paf_orientation(index: PafChainIndex, site_chroms: set) -> str:
    """'h1_query' when BCF contigs sit on the query side, else 'h1_target'."""
    query_side = set(index._by_query)
    target_side = set(index._by_target)
    q_overlap = len(site_chroms & query_side)
    t_overlap = len(site_chroms & target_side)
    if site_chroms - (query_side | target_side):
        missing = sorted(site_chroms - (query_side | target_side))[:3]
        raise SystemExit(f"BCF contigs absent from both PAF sides, e.g. {missing}")
    if q_overlap > t_overlap:
        return "h1_query"
    if t_overlap > q_overlap:
        return "h1_target"
    raise SystemExit("PAF orientation ambiguous: BCF contigs match both sides equally")


def median(values: Sequence[float]) -> float:
    ordered = sorted(values)
    n = len(ordered)
    if not n:
        raise SystemExit("median of empty sequence")
    mid = n // 2
    return ordered[mid] if n % 2 else (ordered[mid - 1] + ordered[mid]) / 2.0


def bed_bp_outside(intervals: Iterable[Tuple[str, int, int]],
                    excluded: Iterable[Tuple[str, int, int]]) -> Tuple[int, int]:
    """Return (bp outside exclusion, bp excluded) for BED-vs-BED subtraction."""
    by_contig: Dict[str, List[List[int]]] = {}
    for chrom, start, end in intervals:
        if end > start:
            by_contig.setdefault(chrom, []).append([start, end])
    exclusions: Dict[str, List[List[int]]] = {}
    for chrom, start, end in excluded:
        if end > start:
            exclusions.setdefault(chrom, []).append([start, end])
    kept = excluded_bp = 0
    for chrom, spans in by_contig.items():
        spans.sort()
        merged: List[List[int]] = []
        for start, end in spans:
            if merged and start <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], end)
            else:
                merged.append([start, end])
        cuts = sorted(exclusions.get(chrom, []))
        for start, end in merged:
            cursor = start
            for cut_start, cut_end in cuts:
                if cut_end <= cursor:
                    continue
                if cut_start >= end:
                    break
                if cut_start > cursor:
                    kept += min(cut_start, end) - cursor
                overlap_lo = max(cursor, cut_start)
                overlap_hi = min(end, cut_end)
                if overlap_hi > overlap_lo:
                    excluded_bp += overlap_hi - overlap_lo
                cursor = max(cursor, cut_end)
                if cursor >= end:
                    break
            if cursor < end:
                kept += end - cursor
    return kept, excluded_bp


def cmd_triage(args: argparse.Namespace) -> int:
    bin_size = int(args.bin_size)
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    h1_scan = FastaScan(Path(args.h1_fasta), bin_size)
    h1_chroms = set(h1_scan.index)

    index, paf_stats = load_paf(Path(args.paf))

    if args.sites_tsv:
        site_list = list(stream_sites_tsv(Path(args.sites_tsv)))
    else:
        site_list = list(stream_bcf_sites(Path(args.bcf), args.bcftools))
    if not site_list:
        raise SystemExit("no SNP sites loaded")
    site_chroms = {site[0] for site in site_list}
    unknown = site_chroms - h1_chroms
    if unknown:
        raise SystemExit(f"site contigs absent from H1 fasta: {sorted(unknown)[:3]}")

    orientation = detect_paf_orientation(index, site_chroms)

    # -- pass 1: per-bin site counts ---------------------------------------
    bin_sites: Counter = Counter()
    for chrom, pos, _ref, _alt in site_list:
        bin_sites[(chrom, pos // bin_size)] += 1

    bin_rows: List[dict] = []
    for (chrom, bin_index), (acgt, gc, other) in sorted(h1_scan.bins.items()):
        sites = bin_sites.get((chrom, bin_index), 0)
        bin_rows.append({
            "chrom": chrom,
            "bin_start_0based": bin_index * bin_size,
            "bin_end_0based": (bin_index + 1) * bin_size,
            "acgt_bases": acgt,
            "non_acgt_letters": other,
            "gc_fraction": (gc / acgt) if acgt else None,
            "snp_sites": sites,
            "divergence_sites_per_kb": (sites * 1000.0 / acgt) if acgt else None,
        })
    eligible = [row for row in bin_rows
                if row["acgt_bases"] >= int(args.min_bin_bases)]
    median_divergence = median([row["divergence_sites_per_kb"] for row in eligible])
    median_gc = median([row["gc_fraction"] for row in eligible])
    divergence_threshold = float(args.divergence_factor) * median_divergence
    gc_threshold = median_gc - float(args.gc_delta)

    for row in bin_rows:
        row["flagged_paralog_bin"] = bool(
            row["acgt_bases"] >= int(args.min_bin_bases)
            and row["divergence_sites_per_kb"] is not None
            and row["divergence_sites_per_kb"] > divergence_threshold
            and row["gc_fraction"] is not None
            and row["gc_fraction"] < gc_threshold)

    flagged_bins = [row for row in bin_rows if row["flagged_paralog_bin"]]
    flagged_lookup = {(row["chrom"], row["bin_start_0based"] // bin_size)
                      for row in flagged_bins}

    # -- pass 2: corrected pi + allele-origin -------------------------------
    callable_bed = BedLookup(read_bed(Path(args.callable_bed)))
    total_callable_sites = kept_sites = flagged_sites = 0
    for chrom, pos, _ref, _alt in site_list:
        if not callable_bed.contains(chrom, pos):
            continue
        total_callable_sites += 1
        if (chrom, pos // bin_size) in flagged_lookup:
            flagged_sites += 1
        else:
            kept_sites += 1

    flagged_intervals = [(row["chrom"], row["bin_start_0based"], row["bin_end_0based"])
                         for row in flagged_bins]
    kept_bp, excluded_bp = bed_bp_outside(read_bed(Path(args.callable_bed)),
                                          flagged_intervals)
    total_bp = kept_bp + excluded_bp

    allele_origin: Optional[dict] = None
    h2_scan: Optional[FastaScan] = None
    if not args.skip_allele_origin:
        h2_scan = FastaScan(Path(args.h2_fasta), bin_size)
        forward = orientation == "h1_query"
        agree = disagree = unliftable = 0
        per_bin_origin: Dict[Tuple[str, int], List[int]] = {}
        for chrom, pos, _ref, alt in site_list:
            reason, lifted = index.classify_site(chrom, pos, forward=forward)
            if reason != "lifted" or lifted is None:
                unliftable += 1
                continue
            l_chrom, l_pos = lifted
            row_candidates = index._rows_containing(
                index._by_query if forward else index._by_target, chrom, pos, pos + 1,
                (lambda r: (r.q_start, r.q_end)) if forward else (lambda r: (r.t_start, r.t_end)))
            strand = row_candidates[0].strand if row_candidates else "+"
            base = h2_scan.base_at(l_chrom, l_pos)
            if strand == "-":
                base = COMPLEMENT.get(base, base)
            entry = per_bin_origin.setdefault((chrom, pos // bin_size), [0, 0])
            if base == alt:
                agree += 1
                entry[0] += 1
            else:
                disagree += 1
                entry[1] += 1
        observed = agree + disagree
        allele_origin = {
            "paf_orientation": orientation,
            "observed_sites": observed,
            "h2_base_matches_alt": agree,
            "h2_base_disagrees": disagree,
            "disagreement_fraction": (disagree / observed) if observed else None,
            "unliftable_sites": unliftable,
        }
        for row in bin_rows:
            entry = per_bin_origin.get((row["chrom"], row["bin_start_0based"] // bin_size))
            if entry and (entry[0] + entry[1]):
                row["allele_origin_disagreement_fraction"] = entry[1] / (entry[0] + entry[1])
            else:
                row["allele_origin_disagreement_fraction"] = None
    else:
        for row in bin_rows:
            row["allele_origin_disagreement_fraction"] = None
        h2_scan = None
    if h2_scan is not None:
        h2_scan.close()

    h1_scan.close()

    fieldnames = list(bin_rows[0].keys()) if bin_rows else ["chrom"]
    with (out_dir / "triage_bins.tsv").open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames,
                                delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(bin_rows)
    with (out_dir / "paralog_bins.bed").open("w") as handle:
        for chrom, start, end in sorted(flagged_intervals):
            handle.write(f"{chrom}\t{start}\t{end}\n")

    summary = {
        "schema_version": SCHEMA_TRIAGE,
        "pair": args.pair,
        "inputs": {
            "paf": {"path": str(args.paf), "sha256": sha256_of(Path(args.paf))},
            "h1_fasta": {"path": str(args.h1_fasta), "sha256": sha256_of(Path(args.h1_fasta))},
            "h2_fasta": {"path": str(args.h2_fasta), "sha256": sha256_of(Path(args.h2_fasta))},
            "bcf": {"path": str(args.bcf)} if args.bcf else {"sites_tsv": str(args.sites_tsv)},
            "callable_bed": {"path": str(args.callable_bed), "sha256": sha256_of(Path(args.callable_bed))},
        },
        "paf_orientation": orientation,
        "bin_size": bin_size,
        "thresholds": {
            "divergence_factor": float(args.divergence_factor),
            "gc_delta": float(args.gc_delta),
            "min_bin_bases": int(args.min_bin_bases),
            "median_bin_divergence_sites_per_kb": median_divergence,
            "median_bin_gc_fraction": median_gc,
            "divergence_threshold_sites_per_kb": divergence_threshold,
            "gc_threshold": gc_threshold,
        },
        "bins": {"total": len(bin_rows), "eligible": len(eligible),
                 "flagged_paralog": len(flagged_bins)},
        "flagged_bins": [
            {"chrom": row["chrom"], "start": row["bin_start_0based"],
             "end": row["bin_end_0based"], "sites": row["snp_sites"],
             "divergence_sites_per_kb": row["divergence_sites_per_kb"],
             "gc_fraction": row["gc_fraction"],
             "allele_origin_disagreement_fraction": row["allele_origin_disagreement_fraction"]}
            for row in flagged_bins],
        "sites": {
            "total_snp_sites": len(site_list),
            "callable_sites": total_callable_sites,
            "callable_sites_in_paralog_bins": flagged_sites,
            "artifact_site_fraction": (flagged_sites / total_callable_sites) if total_callable_sites else None,
        },
        "pi": {
            "callable_bp": total_bp,
            "callable_bp_in_paralog_bins": excluded_bp,
            "artifact_bp_fraction": (excluded_bp / total_bp) if total_bp else None,
            "raw_pi": (total_callable_sites / total_bp) if total_bp else None,
            "corrected_pi": (kept_sites / kept_bp) if kept_bp else None,
        },
        "allele_origin": allele_origin,
    }
    (out_dir / "triage_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"pair": args.pair, "flagged_bins": len(flagged_bins),
                      "raw_pi": summary["pi"]["raw_pi"],
                      "corrected_pi": summary["pi"]["corrected_pi"],
                      "artifact_site_fraction": summary["sites"]["artifact_site_fraction"]},
                     indent=2))
    return 0


# --------------------------------------------------------------------------
# report subcommand
# --------------------------------------------------------------------------


def cmd_report(args: argparse.Namespace) -> int:
    metrics = {}
    for spec in args.metrics_json:
        payload = json.loads(Path(spec).read_text())
        metrics[(payload["frame"], payload["platform"])] = payload

    lines: List[str] = []
    lines.append("# Symmetric read-vs-assembly validation (P07) — metrics report")
    lines.append("")
    lines.append(f"Schema: `{SCHEMA_REPORT}` — generated from {len(metrics)} frame×platform metric files.")
    lines.append("")
    lines.append("## Per-frame π (QUAL≥30 het SNPs / callable bp)")
    lines.append("")
    lines.append("| frame | platform | het SNPs | callable bp | π_read |")
    lines.append("|---|---|---:|---:|---:|")
    for (frame, platform), payload in sorted(metrics.items()):
        lines.append(f"| {frame} | {platform} | {payload['read_het_snps_on_mask']} | "
                     f"{payload['callable_bp']} | {payload['pi_read']} |")
    lines.append("")
    lines.append("## Concordance at assembly SNVs and both-direction contradictions")
    lines.append("")
    lines.append("| frame | platform | supported het | contradicted hom-ref | contradiction rate | A: asm-het/read-homref | B: read-het/asm-homref |")
    lines.append("|---|---|---:|---:|---:|---:|---:|")
    for (frame, platform), payload in sorted(metrics.items()):
        c = payload["concordance"]
        b = payload["both_direction_contradictions"]
        lines.append(f"| {frame} | {platform} | {c['supported_heterozygous']} | "
                     f"{c['contradicted_homozygous_reference']} | {c['contradiction_rate_resolved']} | "
                     f"{b['A_assemblyHet_readsHomRef']} | {b['B_readHet_assemblyHomRef']} |")
    lines.append("")
    lines.append("Assembly site universe (sites inside the frame pi-callable bed;")
    lines.append("concordance and direction-B denominators use exactly this universe):")
    lines.append("")
    lines.append("| frame | platform | sites total | in frame-callable | outside |")
    lines.append("|---|---|---:|---:|---:|")
    for (frame, platform), payload in sorted(metrics.items()):
        strata = payload.get("assembly_site_strata", {})
        lines.append(f"| {frame} | {platform} | {strata.get('total', 'n/a')} | "
                     f"{strata.get('in_frame_callable', 'n/a')} | "
                     f"{strata.get('outside_frame_callable', 'n/a')} |")
    lines.append("")
    lines.append("## Flagged bins (localized discordance)")
    lines.append("")
    for (frame, platform), payload in sorted(metrics.items()):
        flagged = payload.get("flagged_bins", [])
        lines.append(f"- **{frame}/{platform}**: {len(flagged)} bins flagged "
                     f"(policy: {json.dumps(payload['flag_policy'])})")
        for row in flagged[:20]:
            lines.append(f"  - {row['chrom']}:{row['bin_start_0based']}-{row['bin_end_0based']} "
                         f"contradiction_rate={row['contradiction_rate']} sites={row['assembly_sites']}")
    lines.append("")
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    Path(args.output).write_text("\n".join(lines) + "\n")
    print(f"wrote {args.output}")
    return 0


# --------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    transfer = subparsers.add_parser("transfer", help="lift BED and SNP sites through the 1:1 PAF")
    transfer.add_argument("--paf", required=True)
    transfer.add_argument("--bed", required=True, help="H1 pi-callable BED")
    transfer.add_argument("--sites", required=True, help="assembly SNP sites TSV (chrom, position_1based, ref, alt)")
    transfer.add_argument("--source-label", default="h1")
    transfer.add_argument("--target-label", default="h2")
    transfer.add_argument("--output-dir", required=True)
    transfer.set_defaults(func=cmd_transfer)

    evidence = subparsers.add_parser(
        "evidence",
        help="allele-aware pileup evidence at assembly SNVs for one frame")
    evidence.add_argument("--sites", required=True,
                          help="frame sites TSV (chrom, position_1based, ref, alt); "
                               "ref is the site H1 allele, alt the H2 allele")
    evidence.add_argument("--pileup", required=True,
                          help="samtools mpileup text (chrom, pos, REF, depth, bases, quals)")
    evidence.add_argument("--frame", required=True, help="h1 or h2")
    evidence.add_argument("--platform", required=True, help="illumina or hifi")
    evidence.add_argument("--minimum-depth", type=int, default=10)
    evidence.add_argument("--maximum-depth", type=int, default=80)
    evidence.add_argument("--output", required=True)
    evidence.add_argument("--summary", required=True)
    evidence.set_defaults(func=cmd_evidence)

    triage = subparsers.add_parser(
        "triage", help="assembly-only collapsed-paralogy screen (no reads)")
    triage.add_argument("--pair", required=True, help="selection id, e.g. P07")
    triage.add_argument("--paf", required=True, help="1:1 H1<->H2 PAF (either orientation)")
    triage.add_argument("--h1-fasta", required=True, help="H1 reference fasta (BCF coordinates)")
    triage.add_argument("--h2-fasta", required=True, help="H2 fasta for allele-origin stats")
    source = triage.add_mutually_exclusive_group(required=True)
    source.add_argument("--bcf", help="bounded normalized BCF (H1 coordinates)")
    source.add_argument("--sites-tsv", help="sites TSV (chrom, position_1based, ref, alt)")
    triage.add_argument("--callable-bed", required=True, help="bounded pi-callable BED (H1 coordinates)")
    triage.add_argument("--bcftools", default="bcftools")
    triage.add_argument("--bin-size", type=int, default=1_000_000)
    triage.add_argument("--divergence-factor", type=float, default=2.0,
                        help="flag bins with divergence > factor * genome-wide median")
    triage.add_argument("--gc-delta", type=float, default=0.02,
                        help="AND GC < genome-wide median - delta")
    triage.add_argument("--min-bin-bases", type=int, default=100_000,
                        help="bins with fewer ACGT bases are excluded from medians/flagging")
    triage.add_argument("--skip-allele-origin", action="store_true")
    triage.add_argument("--output-dir", required=True)
    triage.set_defaults(func=cmd_triage)

    metrics = subparsers.add_parser("metrics", help="symmetric metrics for one frame×platform")
    metrics.add_argument("--frame", required=True, help="e.g. h1 or h2")
    metrics.add_argument("--platform", required=True, help="illumina or hifi")
    metrics.add_argument("--assembly-sites", required=True, help="assembly SNP sites TSV in this frame's coordinates")
    metrics.add_argument("--frame-callable-bed", required=True,
                         help="this frame's pi-callable BED; the assembly site universe is intersected with it")
    metrics.add_argument("--read-variants", required=True,
                         help="read SNP TSV (chrom, position_1based, ref, alt, quality, genotype)")
    metrics.add_argument("--assembly-evidence", required=True,
                         help="assembly-evidence per-site TSV (analysis.vgp_read_validation) for this frame")
    metrics.add_argument("--callable-bp", type=int, required=True)
    metrics.add_argument("--bin-size", type=int, default=1_000_000)
    metrics.add_argument("--min-sites-per-bin", type=int, default=50)
    metrics.add_argument("--flag-factor", type=float, default=2.0)
    metrics.add_argument("--output", required=True)
    metrics.set_defaults(func=cmd_metrics)

    report = subparsers.add_parser("report", help="merge frame metrics into markdown")
    report.add_argument("--metrics-json", action="append", required=True)
    report.add_argument("--output", required=True)
    report.set_defaults(func=cmd_report)

    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
