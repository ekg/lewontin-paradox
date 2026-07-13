#!/usr/bin/env python3
"""Shared, fail-closed utilities for the Tier 3 analysis pipelines.

Coordinates exposed by this module are always zero-based and half-open unless
the function name explicitly says otherwise.  The implementation intentionally
uses only the Python standard library so that acquisition and provenance checks
can run before optional scientific Python packages are imported.
"""

from __future__ import annotations

import collections
import gzip
import hashlib
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Mapping, Optional, Sequence, Set, Tuple, Union


DNA = frozenset("ACGT")
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


class Tier3ValidationError(ValueError):
    """An input violates a frozen Tier 3 policy or provenance invariant."""


def _open_text(path: Union[str, Path]):
    path = Path(path)
    if path.suffix in {".gz", ".bgz"}:
        return gzip.open(str(path), "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def sha256_file(path: Union[str, Path], chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify_file(path: Union[str, Path], expected_sha256: str, expected_size: int) -> None:
    path = Path(path)
    if not path.is_file():
        raise Tier3ValidationError(f"file does not exist: {path}")
    observed_size = path.stat().st_size
    if observed_size != expected_size:
        raise Tier3ValidationError(
            f"size mismatch for {path}: expected {expected_size}, observed {observed_size}"
        )
    observed_sha256 = sha256_file(path)
    if observed_sha256 != expected_sha256:
        raise Tier3ValidationError(
            f"SHA-256 mismatch for {path}: expected {expected_sha256}, observed {observed_sha256}"
        )


def read_fasta(path: Union[str, Path]) -> Dict[str, str]:
    """Read FASTA without silently normalizing duplicate or empty contig IDs."""

    sequences: Dict[str, List[str]] = collections.OrderedDict()
    name: Optional[str] = None
    with _open_text(path) as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].split(None, 1)[0]
                if not name:
                    raise Tier3ValidationError(f"empty FASTA name at line {line_number}")
                if name in sequences:
                    raise Tier3ValidationError(f"duplicate FASTA contig {name!r}")
                sequences[name] = []
            else:
                if name is None:
                    raise Tier3ValidationError(f"FASTA sequence before header at line {line_number}")
                if re.search(r"\s", line):
                    raise Tier3ValidationError(f"whitespace inside FASTA sequence at line {line_number}")
                sequences[name].append(line.upper())
    if not sequences:
        raise Tier3ValidationError("FASTA has no records")
    return {key: "".join(parts) for key, parts in sequences.items()}


def fasta_dictionary(fasta: Mapping[str, str]) -> Dict[str, int]:
    return {name: len(sequence) for name, sequence in fasta.items()}


def resolve_contig_aliases(
    fasta_contigs: Mapping[str, int],
    annotation_contigs: Mapping[str, int],
    declared_aliases: Mapping[str, str],
) -> Dict[str, str]:
    """Validate an explicit annotation-to-FASTA, length-preserving bijection."""

    resolved: Dict[str, str] = {}
    targets: Set[str] = set()
    for annotation_name, annotation_length in annotation_contigs.items():
        if annotation_name in fasta_contigs:
            target = annotation_name
            if annotation_name in declared_aliases and declared_aliases[annotation_name] != target:
                raise Tier3ValidationError(f"declared alias for exact contig {annotation_name!r} is inconsistent")
        elif annotation_name in declared_aliases:
            target = declared_aliases[annotation_name]
        else:
            raise Tier3ValidationError(f"undeclared annotation contig alias {annotation_name!r}")
        if target not in fasta_contigs:
            raise Tier3ValidationError(f"annotation contig {annotation_name!r} maps to absent FASTA contig {target!r}")
        if target in targets:
            raise Tier3ValidationError(f"contig mapping is not one-to-one: multiple names map to {target!r}")
        if fasta_contigs[target] != annotation_length:
            raise Tier3ValidationError(
                f"contig length mismatch for {annotation_name!r}->{target!r}: "
                f"annotation {annotation_length}, FASTA {fasta_contigs[target]}"
            )
        targets.add(target)
        resolved[annotation_name] = target
    extra_aliases = set(declared_aliases) - set(annotation_contigs)
    if extra_aliases:
        raise Tier3ValidationError(f"aliases declared for absent annotation contigs: {sorted(extra_aliases)!r}")
    return resolved


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


def parse_attributes(text: str) -> Dict[str, str]:
    attributes: Dict[str, str] = {}
    for item in text.strip().strip(";").split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
        elif " " in item:  # GTF key "value"
            key, value = item.split(None, 1)
            value = value.strip().strip('"')
        else:
            raise Tier3ValidationError(f"malformed GFF/GTF attribute {item!r}")
        attributes[key.strip()] = value.strip()
    return attributes


@dataclass(frozen=True)
class CDSSegment:
    contig: str
    start: int
    end: int
    strand: str
    phase: int
    line_number: int


@dataclass
class Transcript:
    transcript_id: str
    gene_id: str
    segments: List[CDSSegment] = field(default_factory=list)
    provider_canonical: bool = False


@dataclass
class GFFAnnotation:
    sequence_regions: Dict[str, int]
    transcripts: Dict[str, Transcript]
    canonical_transcripts: Dict[str, str]


def parse_gff(path: Union[str, Path]) -> GFFAnnotation:
    sequence_regions: Dict[str, int] = {}
    transcripts: Dict[str, Transcript] = {}
    transcript_gene: Dict[str, str] = {}
    canonical_ids: Set[str] = set()

    with _open_text(path) as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.rstrip("\n")
            if line.startswith("##sequence-region"):
                fields = line.split()
                if len(fields) != 4:
                    raise Tier3ValidationError(f"malformed ##sequence-region at line {line_number}")
                name, first, last = fields[1], int(fields[2]), int(fields[3])
                if first != 1 or last < first:
                    raise Tier3ValidationError(f"invalid ##sequence-region at line {line_number}")
                if name in sequence_regions and sequence_regions[name] != last:
                    raise Tier3ValidationError(f"conflicting ##sequence-region for {name!r}")
                sequence_regions[name] = last
                continue
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != 9:
                raise Tier3ValidationError(f"GFF/GTF line {line_number} does not have nine columns")
            contig, _source, feature, start_text, end_text, _score, strand, phase_text, attrs_text = fields
            attrs = parse_attributes(attrs_text)
            start, end = int(start_text) - 1, int(end_text)
            if start < 0 or end <= start:
                raise Tier3ValidationError(f"invalid coordinates at GFF line {line_number}")
            if feature.lower() in {"mrna", "transcript"}:
                transcript_id = attrs.get("ID") or attrs.get("transcript_id")
                gene_id = attrs.get("Parent") or attrs.get("gene_id")
                if not transcript_id or not gene_id:
                    raise Tier3ValidationError(f"transcript lacks stable ID/parent at line {line_number}")
                transcript_gene[transcript_id] = gene_id.split(",", 1)[0]
                tags = {tag.strip().lower() for tag in attrs.get("tag", "").split(",")}
                if tags.intersection({"canonical", "mane_select", "appris_principal_1"}):
                    canonical_ids.add(transcript_id)
                transcripts.setdefault(transcript_id, Transcript(transcript_id, transcript_gene[transcript_id]))
                continue
            if feature != "CDS":
                continue
            if strand not in {"+", "-"}:
                raise Tier3ValidationError(f"CDS has invalid strand at line {line_number}")
            if phase_text not in {"0", "1", "2"}:
                raise Tier3ValidationError(f"CDS has invalid phase at line {line_number}")
            parent_text = attrs.get("Parent") or attrs.get("transcript_id")
            if not parent_text:
                raise Tier3ValidationError(f"CDS lacks transcript parent at line {line_number}")
            for transcript_id in parent_text.split(","):
                transcript_id = transcript_id.strip()
                gene_id = transcript_gene.get(transcript_id, attrs.get("gene_id", transcript_id))
                transcript = transcripts.setdefault(transcript_id, Transcript(transcript_id, gene_id))
                segment = CDSSegment(contig, start, end, strand, int(phase_text), line_number)
                if segment in transcript.segments:
                    raise Tier3ValidationError(f"duplicate CDS segment at line {line_number}")
                transcript.segments.append(segment)

    if not sequence_regions:
        raise Tier3ValidationError("annotation has no ##sequence-region declarations")
    if not transcripts or not any(transcript.segments for transcript in transcripts.values()):
        raise Tier3ValidationError("annotation has no transcript-associated CDS")

    genes: Dict[str, List[Transcript]] = collections.defaultdict(list)
    for transcript in transcripts.values():
        if not transcript.segments:
            continue
        transcript.provider_canonical = transcript.transcript_id in canonical_ids
        genes[transcript.gene_id].append(transcript)
    canonical: Dict[str, str] = {}
    for gene_id, candidates in genes.items():
        flagged = [candidate for candidate in candidates if candidate.provider_canonical]
        if len(flagged) > 1:
            raise Tier3ValidationError(f"gene {gene_id!r} has multiple provider-canonical transcripts")
        if flagged:
            chosen = flagged[0]
        else:
            chosen = sorted(
                candidates,
                key=lambda candidate: (-_phase_trimmed_length(candidate), candidate.transcript_id.encode("utf-8")),
            )[0]
        canonical[gene_id] = chosen.transcript_id
    return GFFAnnotation(sequence_regions, transcripts, canonical)


def _ordered_segments(transcript: Transcript) -> List[CDSSegment]:
    if not transcript.segments:
        raise Tier3ValidationError(f"transcript {transcript.transcript_id!r} has no CDS")
    contigs = {segment.contig for segment in transcript.segments}
    strands = {segment.strand for segment in transcript.segments}
    if len(contigs) != 1 or len(strands) != 1:
        raise Tier3ValidationError(f"transcript {transcript.transcript_id!r} spans contigs or strands")
    reverse = next(iter(strands)) == "-"
    segments = sorted(transcript.segments, key=lambda segment: segment.start, reverse=reverse)
    for previous, current in zip(sorted(transcript.segments, key=lambda item: item.start), sorted(transcript.segments, key=lambda item: item.start)[1:]):
        if current.start < previous.end:
            raise Tier3ValidationError(f"transcript {transcript.transcript_id!r} has overlapping CDS segments")
    return segments


def _phase_trimmed_length(transcript: Transcript) -> int:
    return sum(segment.end - segment.start - segment.phase for segment in _ordered_segments(transcript))


def reconstruct_cds_with_positions(
    fasta: Mapping[str, str], transcript: Transcript
) -> Tuple[str, List[Tuple[str, int]]]:
    sequence_parts: List[str] = []
    positions: List[Tuple[str, int]] = []
    for segment in _ordered_segments(transcript):
        if segment.contig not in fasta:
            raise Tier3ValidationError(f"CDS contig {segment.contig!r} is absent from FASTA")
        contig_sequence = fasta[segment.contig]
        if segment.end > len(contig_sequence):
            raise Tier3ValidationError(f"CDS at line {segment.line_number} exceeds FASTA contig")
        piece = contig_sequence[segment.start : segment.end].upper()
        piece_positions = [(segment.contig, position) for position in range(segment.start, segment.end)]
        if segment.strand == "-":
            piece = reverse_complement(piece)
            piece_positions.reverse()
        if segment.phase >= len(piece):
            raise Tier3ValidationError(f"phase consumes CDS segment at line {segment.line_number}")
        sequence_parts.append(piece[segment.phase :])
        positions.extend(piece_positions[segment.phase :])
    sequence = "".join(sequence_parts)
    if len(sequence) % 3:
        raise Tier3ValidationError(
            f"phase-handled CDS length for {transcript.transcript_id!r} is not divisible by three"
        )
    if not sequence or any(base not in DNA for base in sequence):
        raise Tier3ValidationError(f"CDS {transcript.transcript_id!r} is empty or ambiguous")
    return sequence, positions


def reconstruct_cds(fasta: Mapping[str, str], transcript: Transcript) -> str:
    return reconstruct_cds_with_positions(fasta, transcript)[0]


# NCBI nuclear translation table 1.  Other codes must be supplied explicitly
# by future policy versions rather than being guessed from organism names.
_GENETIC_CODE_1 = {
    codon: amino
    for amino, codons in {
        "F": ("TTT", "TTC"), "L": ("TTA", "TTG", "CTT", "CTC", "CTA", "CTG"),
        "I": ("ATT", "ATC", "ATA"), "M": ("ATG",), "V": ("GTT", "GTC", "GTA", "GTG"),
        "S": ("TCT", "TCC", "TCA", "TCG", "AGT", "AGC"), "P": ("CCT", "CCC", "CCA", "CCG"),
        "T": ("ACT", "ACC", "ACA", "ACG"), "A": ("GCT", "GCC", "GCA", "GCG"),
        "Y": ("TAT", "TAC"), "*": ("TAA", "TAG", "TGA"), "H": ("CAT", "CAC"),
        "Q": ("CAA", "CAG"), "N": ("AAT", "AAC"), "K": ("AAA", "AAG"),
        "D": ("GAT", "GAC"), "E": ("GAA", "GAG"), "C": ("TGT", "TGC"),
        "W": ("TGG",), "R": ("CGT", "CGC", "CGA", "CGG", "AGA", "AGG"),
        "G": ("GGT", "GGC", "GGA", "GGG"),
    }.items()
    for codon in codons
}


def is_fourfold_codon(codon: str, genetic_code: int = 1) -> bool:
    if genetic_code != 1:
        raise Tier3ValidationError(f"unsupported nuclear genetic code {genetic_code}; only table 1 is frozen")
    if len(codon) != 3 or any(base not in DNA for base in codon):
        return False
    amino_acids = {_GENETIC_CODE_1[codon[:2] + base] for base in "ACGT"}
    return len(amino_acids) == 1 and "*" not in amino_acids


def collect_fourfold_sites(
    fasta: Mapping[str, str], annotation: GFFAnnotation, genetic_code: int = 1
) -> Tuple[Set[Tuple[str, int]], Set[Tuple[str, int]]]:
    """Return unambiguous 4D sites and overlap/frame-discordant exclusions."""

    observations: Dict[Tuple[str, int], Set[Tuple[int, bool]]] = collections.defaultdict(set)
    for transcript_id in annotation.canonical_transcripts.values():
        transcript = annotation.transcripts[transcript_id]
        sequence, positions = reconstruct_cds_with_positions(fasta, transcript)
        amino_acids = [_GENETIC_CODE_1[sequence[offset : offset + 3]] for offset in range(0, len(sequence), 3)]
        if "*" in amino_acids[:-1]:
            raise Tier3ValidationError(f"internal stop in transcript {transcript_id!r}")
        effective_length = len(sequence) - (3 if amino_acids and amino_acids[-1] == "*" else 0)
        for offset in range(effective_length):
            codon_start = offset - offset % 3
            is_4d = offset % 3 == 2 and is_fourfold_codon(sequence[codon_start : codon_start + 3], genetic_code)
            observations[positions[offset]].add((offset % 3, is_4d))
    excluded = {position for position, states in observations.items() if len(states) > 1}
    sites = {
        position
        for position, states in observations.items()
        if position not in excluded and next(iter(states)) == (2, True)
    }
    return sites, excluded


def parse_gt(
    genotype: Union[str, Sequence[Optional[int]]],
    expected_ploidy: int,
    haploidize_heterozygous: bool = False,
) -> Tuple[Optional[int], ...]:
    if expected_ploidy not in {1, 2}:
        raise Tier3ValidationError(f"unsupported expected ploidy {expected_ploidy}")
    if isinstance(genotype, str):
        fields = re.split(r"[/|]", genotype)
        alleles: Tuple[Optional[int], ...] = tuple(None if item == "." else int(item) for item in fields)
    else:
        alleles = tuple(genotype)
    if haploidize_heterozygous and expected_ploidy == 1:
        called = {allele for allele in alleles if allele is not None}
        if len(alleles) == 2 and len(called) == 1 and None not in alleles:
            return (next(iter(called)),)
        return (None,)
    if len(alleles) != expected_ploidy:
        raise Tier3ValidationError(
            f"genotype ploidy {len(alleles)} does not match declared ploidy {expected_ploidy}"
        )
    if any(allele is not None and allele < 0 for allele in alleles):
        raise Tier3ValidationError("negative genotype allele index")
    return alleles


def allele_pairwise_diversity(
    genotypes: Iterable[Sequence[Optional[int]]], minimum_called: int
) -> Optional[float]:
    called = [allele for genotype in genotypes for allele in genotype if allele is not None]
    if len(called) < minimum_called or len(called) < 2:
        return None
    counts = collections.Counter(called)
    n = len(called)
    return 1.0 - sum(count * (count - 1) for count in counts.values()) / (n * (n - 1))


def merge_intervals(intervals: Iterable[Tuple[int, int]]) -> List[Tuple[int, int]]:
    merged: List[List[int]] = []
    for start, end in sorted(intervals):
        if start < 0 or end < start:
            raise Tier3ValidationError(f"invalid half-open interval [{start}, {end})")
        if start == end:
            continue
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return [(start, end) for start, end in merged]


def intersect_intervals(
    left: Iterable[Tuple[int, int]], right: Iterable[Tuple[int, int]]
) -> List[Tuple[int, int]]:
    left_merged, right_merged = merge_intervals(left), merge_intervals(right)
    result: List[Tuple[int, int]] = []
    i = j = 0
    while i < len(left_merged) and j < len(right_merged):
        start = max(left_merged[i][0], right_merged[j][0])
        end = min(left_merged[i][1], right_merged[j][1])
        if start < end:
            result.append((start, end))
        if left_merged[i][1] <= right_merged[j][1]:
            i += 1
        else:
            j += 1
    return result


_CIGAR_TOKEN = re.compile(r"([1-9][0-9]*)([=XIDM])")


def parse_extended_cigar(text: str) -> List[Tuple[int, str]]:
    tokens = [(int(length), operation) for length, operation in _CIGAR_TOKEN.findall(text)]
    if not tokens or "".join(f"{length}{operation}" for length, operation in tokens) != text:
        raise Tier3ValidationError(f"invalid CIGAR {text!r}")
    if any(operation == "M" for _, operation in tokens):
        raise Tier3ValidationError("PAF requires extended =/X/I/D CIGAR; M is forbidden")
    return tokens


@dataclass
class PAFTraversalResult:
    callable_positions: Set[Tuple[str, int]]
    snv_positions: Set[Tuple[str, int]]
    operation_counts: Dict[str, int]
    exclusion_counts: Dict[str, int]


def traverse_paf(
    paf_lines: Iterable[str],
    target_fasta: Mapping[str, str],
    query_fasta: Mapping[str, str],
    edge_exclusion_bp: int = 100,
    indel_flank_bp: int = 10,
) -> PAFTraversalResult:
    """Traverse accepted PAF alignments and construct an H1 callable mask.

    Every target position must have exactly one projection.  Edge, gap-flank,
    ambiguous, and multiple-projection exclusions are reason-counted.
    """

    if edge_exclusion_bp < 0 or indel_flank_bp < 0:
        raise Tier3ValidationError("PAF exclusion widths must be non-negative")
    projections: Dict[Tuple[str, int], List[Tuple[str, str]]] = collections.defaultdict(list)
    excluded_reasons: Dict[Tuple[str, int], Set[str]] = collections.defaultdict(set)
    operation_counts = {operation: 0 for operation in "=XID"}

    for line_number, raw in enumerate(paf_lines, 1):
        if not raw.strip() or raw.startswith("#"):
            continue
        fields = raw.rstrip("\n").split("\t")
        if len(fields) < 12:
            raise Tier3ValidationError(f"PAF line {line_number} has fewer than 12 columns")
        qname, qlen_text, qstart_text, qend_text, strand = fields[:5]
        tname, tlen_text, tstart_text, tend_text = fields[5:9]
        qlen, qstart, qend = int(qlen_text), int(qstart_text), int(qend_text)
        tlen, tstart, tend = int(tlen_text), int(tstart_text), int(tend_text)
        if strand not in {"+", "-"}:
            raise Tier3ValidationError(f"invalid PAF strand at line {line_number}")
        if qname not in query_fasta or tname not in target_fasta:
            raise Tier3ValidationError(f"PAF line {line_number} names absent query/target sequence")
        if qlen != len(query_fasta[qname]) or tlen != len(target_fasta[tname]):
            raise Tier3ValidationError(f"PAF sequence length mismatch at line {line_number}")
        if not (0 <= qstart <= qend <= qlen and 0 <= tstart <= tend <= tlen):
            raise Tier3ValidationError(f"PAF coordinate bounds invalid at line {line_number}")
        cigar_tags = [field[5:] for field in fields[12:] if field.startswith("cg:Z:")]
        if len(cigar_tags) != 1:
            raise Tier3ValidationError(f"PAF line {line_number} must have exactly one cg:Z tag")
        cigar = parse_extended_cigar(cigar_tags[0])
        query_piece = query_fasta[qname][qstart:qend].upper()
        if strand == "-":
            query_piece = reverse_complement(query_piece)
        query_offset = 0
        target_position = tstart
        alignment_indel_exclusions: Set[int] = set()

        for length, operation in cigar:
            operation_counts[operation] += length
            if operation in {"=", "X"}:
                for delta in range(length):
                    if query_offset + delta >= len(query_piece) or target_position + delta >= tend:
                        raise Tier3ValidationError(f"CIGAR exceeds PAF span at line {line_number}")
                    target_base = target_fasta[tname][target_position + delta].upper()
                    query_base = query_piece[query_offset + delta].upper()
                    if target_base in DNA and query_base in DNA:
                        observed_equal = target_base == query_base
                        if observed_equal != (operation == "="):
                            raise Tier3ValidationError(
                                f"CIGAR {operation} contradicts FASTA bases at {tname}:{target_position + delta + 1}"
                            )
                    projections[(tname, target_position + delta)].append((operation, query_base))
                query_offset += length
                target_position += length
            elif operation == "I":
                anchor = target_position
                for position in range(max(tstart, anchor - indel_flank_bp), min(tend, anchor + indel_flank_bp + 1)):
                    alignment_indel_exclusions.add(position)
                query_offset += length
            elif operation == "D":
                for position in range(
                    max(tstart, target_position - indel_flank_bp),
                    min(tend, target_position + length + indel_flank_bp),
                ):
                    alignment_indel_exclusions.add(position)
                target_position += length
        if query_offset != qend - qstart or target_position != tend:
            raise Tier3ValidationError(
                f"CIGAR consumption disagrees with PAF spans at line {line_number}: "
                f"query {query_offset}/{qend-qstart}, target {target_position-tstart}/{tend-tstart}"
            )
        for position in range(tstart, min(tend, tstart + edge_exclusion_bp)):
            excluded_reasons[(tname, position)].add("alignment_edge")
        for position in range(max(tstart, tend - edge_exclusion_bp), tend):
            excluded_reasons[(tname, position)].add("alignment_edge")
        for position in alignment_indel_exclusions:
            excluded_reasons[(tname, position)].add("indel_flank")

    callable_positions: Set[Tuple[str, int]] = set()
    snv_positions: Set[Tuple[str, int]] = set()
    for key, values in projections.items():
        if len(values) != 1:
            excluded_reasons[key].add("multiple_projection")
            continue
        operation, query_base = values[0]
        target_base = target_fasta[key[0]][key[1]].upper()
        if target_base not in DNA or query_base not in DNA:
            excluded_reasons[key].add("ambiguous_base")
            continue
        if excluded_reasons.get(key):
            continue
        callable_positions.add(key)
        if operation == "X":
            snv_positions.add(key)
    exclusion_counts = collections.Counter(
        reason for reasons in excluded_reasons.values() for reason in reasons
    )
    return PAFTraversalResult(
        callable_positions,
        snv_positions,
        operation_counts,
        dict(sorted(exclusion_counts.items())),
    )
