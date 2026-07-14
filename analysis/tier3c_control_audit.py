#!/usr/bin/env python3
"""Independent exact-reference GC3 audit for the Tier 3c controls.

This module intentionally does not import ``tier3_common`` or
``tier3c_ncbi_gc``.  It independently parses FASTA/GFF3, reconstructs CDS on
both strands with GFF phase, applies the frozen canonical-transcript and
exclusion rules, and reports both pooled-CDS-base and gene-weighted GC3.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import resource
import tempfile
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple
from urllib.parse import unquote, urlparse


SCHEMA_VERSION = "tier3c-control-audit-v1"
DNA = frozenset("ACGT")
STOP_CODONS = frozenset({"TAA", "TAG", "TGA"})
CANONICAL_TAGS = frozenset(
    {"canonical", "mane_select", "refseq_select", "ensembl_canonical", "appris_principal_1"}
)
COMPLEMENT = str.maketrans("ACGT", "TGCA")
CONTROL_NAMES = {"Drosophila melanogaster", "Homo sapiens"}
LEGACY_BANDS: Mapping[str, Mapping[str, Any]] = {
    "Drosophila melanogaster": {
        "band": [0.50, 0.60],
        "observed_pooled_gc3": 0.6312610989373014,
        "status": "failed_and_preserved",
    },
    "Homo sapiens": {
        "band": [0.47, 0.57],
        "observed_pooled_gc3": 0.5845556551887756,
        "status": "failed_and_preserved",
    },
}
LITERATURE_CONTEXT: Mapping[str, Sequence[Mapping[str, Any]]] = {
    "Drosophila melanogaster": (
        {
            "title": (
                "Codon Usage Bias and Effective Population Sizes on the X Chromosome versus "
                "the Autosomes in Drosophila melanogaster"
            ),
            "url": "https://academic.oup.com/mbe/article/30/4/811/1066476",
            "definition": "mean per-gene GC content of third codon positions",
            "reported": "X=0.688 (95% CI 0.683-0.692); autosomes=0.641 (0.639-0.643)",
            "comparability": (
                "gene-weighted, so comparable to the audit's explicitly separate equal-gene metric; "
                "gene set and annotation release differ"
            ),
        },
        {
            "title": "Codon Usage Differences among Genes Expressed in Different Tissues of Drosophila melanogaster",
            "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC6456009/",
            "definition": "mean and median per-gene GC3 across 13,088 nuclear protein-coding genes",
            "reported": "all-gene mean=0.65; median=0.66",
            "comparability": "gene-weighted; transcript policy and annotation release differ",
        },
    ),
    "Homo sapiens": (
        {
            "title": "Mapping the inter- and intra-genic codon-usage landscape in Homo sapiens",
            "url": "https://academic.oup.com/nargab/article/8/1/lqag024/8503857",
            "definition": "per-gene GC3 across approximately 20,000 Ensembl protein-coding genes",
            "reported": "mean=0.59; median=0.59; IQR=0.44-0.73",
            "comparability": "gene-weighted; provider and canonical-transcript policy differ",
        },
        {
            "title": "GC Content of Early Metazoan Genes and Its Impact on Gene Expression Levels in Mammalian Cell Lines",
            "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC5952964/",
            "definition": "distribution of per-CDS third-position GC in mammalian coding sequences",
            "reported": "mammalian CDS mean=0.59; human mean described as above 0.55",
            "comparability": "coding-sequence reference context; gene set and release differ",
        },
    ),
}


class AuditError(ValueError):
    """A fail-closed provenance, reconstruction, or comparison failure."""


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _verify(path: Path, expected: str, label: str) -> None:
    observed = _sha256(path)
    if observed != expected:
        raise AuditError(f"{label} SHA-256 mismatch: expected {expected}, observed {observed}")


def _local_path(uri: str) -> Path:
    parsed = urlparse(uri)
    if parsed.scheme != "file" or parsed.netloc not in {"", "localhost"}:
        raise AuditError(f"control audit requires a frozen local file URI, observed {uri!r}")
    return Path(unquote(parsed.path))


@dataclass(frozen=True, slots=True)
class FastaRecord:
    offset: int
    length: int
    line_bases: int
    line_bytes: int


class FastaIndex:
    """Read-only constant-width FASTA index built independently in memory."""

    def __init__(self, path: Path):
        self.path = path
        self.records: Dict[str, FastaRecord] = {}
        self._scan()

    def _scan(self) -> None:
        current: Optional[str] = None
        offset = length = line_bases = line_bytes = 0
        previous_bases: Optional[int] = None

        def finish() -> None:
            nonlocal current, offset, length, line_bases, line_bytes, previous_bases
            if current is None:
                return
            if not length or not line_bases:
                raise AuditError(f"empty FASTA record {current!r}")
            self.records[current] = FastaRecord(offset, length, line_bases, line_bytes)
            current = None
            previous_bases = None

        with self.path.open("rb") as handle:
            while True:
                position = handle.tell()
                raw = handle.readline()
                if not raw:
                    break
                if raw.startswith(b">"):
                    finish()
                    name = raw[1:].split(None, 1)[0].decode("ascii")
                    if not name or name in self.records:
                        raise AuditError(f"empty or duplicate FASTA identifier {name!r}")
                    current = name
                    offset = handle.tell()
                    length = line_bases = line_bytes = 0
                    continue
                if current is None:
                    if raw.strip():
                        raise AuditError(f"FASTA sequence appears before a header at byte {position}")
                    continue
                sequence = raw.rstrip(b"\r\n")
                if not sequence:
                    raise AuditError(f"blank sequence line in FASTA record {current!r}")
                if previous_bases is not None and previous_bases != line_bases:
                    raise AuditError(f"non-terminal short FASTA line in record {current!r}")
                if line_bases == 0:
                    line_bases = len(sequence)
                    line_bytes = len(raw)
                elif len(sequence) == line_bases and len(raw) != line_bytes:
                    raise AuditError(f"inconsistent FASTA newline width in record {current!r}")
                previous_bases = len(sequence)
                length += len(sequence)
        finish()
        if not self.records:
            raise AuditError("FASTA contains no records")

    def fetch(self, name: str, start: int, end: int) -> str:
        if name not in self.records:
            raise AuditError(f"CDS contig {name!r} is absent from FASTA")
        record = self.records[name]
        if start < 0 or end <= start or end > record.length:
            raise AuditError(f"invalid FASTA interval {name}:{start}-{end}")
        row, column = divmod(start, record.line_bases)
        seek = record.offset + row * record.line_bytes + column
        required = end - start
        chunks: List[bytes] = []
        observed = 0
        with self.path.open("rb") as handle:
            handle.seek(seek)
            while observed < required:
                raw = handle.readline()
                if not raw or raw.startswith(b">"):
                    raise AuditError(f"truncated FASTA interval {name}:{start}-{end}")
                sequence = raw.rstrip(b"\r\n")
                take = min(required - observed, len(sequence))
                chunks.append(sequence[:take])
                observed += take
        try:
            return b"".join(chunks).decode("ascii").upper()
        except UnicodeDecodeError as error:
            raise AuditError(f"non-ASCII FASTA sequence in {name!r}") from error


@dataclass(frozen=True, slots=True)
class Segment:
    contig: str
    start: int
    end: int
    strand: str
    phase: int
    line_number: int


@dataclass(slots=True)
class Transcript:
    transcript_id: str
    gene_id: str
    segments: List[Segment] = field(default_factory=list)
    provider_canonical: bool = False
    annotation_exception: Optional[str] = None


def _attributes(text: str) -> Dict[str, str]:
    parsed: Dict[str, str] = {}
    for item in text.strip().strip(";").split(";"):
        if not item:
            continue
        if "=" not in item:
            raise AuditError(f"malformed GFF3 attribute {item!r}")
        key, value = item.split("=", 1)
        parsed[key.strip()] = value.strip()
    return parsed


def _parse_gff(path: Path) -> Tuple[Dict[str, int], Dict[str, Transcript], int]:
    regions: Dict[str, int] = {}
    transcripts: Dict[str, Transcript] = {}
    transcript_gene: Dict[str, str] = {}
    canonical_ids = set()
    feature_count = 0
    with path.open("rt", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.rstrip("\n")
            if line.startswith("##sequence-region"):
                fields = line.split()
                if len(fields) != 4 or fields[2] != "1":
                    raise AuditError(f"malformed GFF3 sequence-region at line {line_number}")
                name, length = fields[1], int(fields[3])
                if name in regions and regions[name] != length:
                    raise AuditError(f"conflicting GFF3 sequence-region {name!r}")
                regions[name] = length
                continue
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != 9:
                raise AuditError(f"GFF3 line {line_number} does not have nine columns")
            feature_count += 1
            contig, _source, feature, start_text, end_text, _score, strand, phase_text, text = fields
            attrs = _attributes(text)
            feature_lower = feature.lower()
            if feature_lower in {"mrna", "transcript"}:
                transcript_id = attrs.get("ID")
                gene_id = attrs.get("Parent")
                is_pseudo = attrs.get("pseudo", "").lower() == "true"
                if transcript_id and not gene_id and is_pseudo:
                    gene_id = attrs.get("locus_tag") or attrs.get("gene") or transcript_id
                if not transcript_id or not gene_id:
                    raise AuditError(f"transcript lacks stable ID/parent at GFF3 line {line_number}")
                gene_id = gene_id.split(",", 1)[0]
                transcript_gene[transcript_id] = gene_id
                transcript = transcripts.setdefault(transcript_id, Transcript(transcript_id, gene_id))
                transcript.gene_id = gene_id
                tags = {
                    re.sub(r"[^a-z0-9]+", "_", tag.strip().lower()).strip("_")
                    for tag in attrs.get("tag", "").split(",")
                }
                if tags & CANONICAL_TAGS:
                    canonical_ids.add(transcript_id)
                exception = attrs.get("exception") or attrs.get("transl_except")
                if is_pseudo:
                    exception = exception or "provider_declared_pseudogene"
                if exception:
                    transcript.annotation_exception = exception
                continue
            if feature != "CDS":
                continue
            if strand not in {"+", "-"} or phase_text not in {"0", "1", "2"}:
                raise AuditError(f"invalid CDS strand/phase at GFF3 line {line_number}")
            parent_text = attrs.get("Parent")
            if not parent_text:
                raise AuditError(f"CDS lacks transcript parent at GFF3 line {line_number}")
            start, end = int(start_text) - 1, int(end_text)
            if start < 0 or end <= start:
                raise AuditError(f"invalid CDS coordinates at GFF3 line {line_number}")
            for transcript_id in parent_text.split(","):
                transcript_id = transcript_id.strip()
                gene_id = transcript_gene.get(transcript_id, attrs.get("gene_id", transcript_id))
                transcript = transcripts.setdefault(transcript_id, Transcript(transcript_id, gene_id))
                exception = attrs.get("exception") or attrs.get("transl_except")
                if attrs.get("pseudo", "").lower() == "true":
                    exception = exception or "provider_declared_pseudogene"
                if exception:
                    transcript.annotation_exception = exception
                segment = Segment(contig, start, end, strand, int(phase_text), line_number)
                identity = (contig, start, end, strand, int(phase_text))
                if any(
                    (item.contig, item.start, item.end, item.strand, item.phase) == identity
                    for item in transcript.segments
                ):
                    raise AuditError(f"duplicate CDS segment at GFF3 line {line_number}")
                transcript.segments.append(segment)
    if not regions:
        raise AuditError("annotation has no GFF3 sequence-region declarations")
    for transcript_id in canonical_ids:
        if transcript_id in transcripts:
            transcripts[transcript_id].provider_canonical = True
    return regions, transcripts, feature_count


def _read_mapping(path: Path) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    with path.open("rt", encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        if header != ["annotation_contig", "fasta_contig"]:
            raise AuditError("contig mapping has an unexpected header")
        for line_number, raw in enumerate(handle, 2):
            fields = raw.rstrip("\n").split("\t")
            if len(fields) != 2 or not all(fields):
                raise AuditError(f"malformed contig mapping line {line_number}")
            annotation, fasta = fields
            if annotation in mapping:
                raise AuditError(f"duplicate annotation contig mapping {annotation!r}")
            mapping[annotation] = fasta
    if len(set(mapping.values())) != len(mapping):
        raise AuditError("contig mapping is not one-to-one")
    return mapping


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


def _reconstruct(index: FastaIndex, transcript: Transcript, mapping: Mapping[str, str]) -> str:
    if not transcript.segments:
        raise AuditError(f"transcript {transcript.transcript_id!r} has no CDS")
    contigs = {item.contig for item in transcript.segments}
    strands = {item.strand for item in transcript.segments}
    if len(contigs) != 1 or len(strands) != 1:
        raise AuditError(f"transcript {transcript.transcript_id!r} spans contigs or strands")
    reverse = next(iter(strands)) == "-"
    genomic = sorted(transcript.segments, key=lambda item: item.start)
    for previous, current in zip(genomic, genomic[1:]):
        if current.start < previous.end:
            raise AuditError(f"transcript {transcript.transcript_id!r} has overlapping CDS")
    ordered = list(reversed(genomic)) if reverse else genomic
    cumulative = ordered[0].end - ordered[0].start - ordered[0].phase
    if cumulative <= 0:
        raise AuditError(f"phase consumes first CDS of {transcript.transcript_id!r}")
    for segment in ordered[1:]:
        expected = (3 - cumulative % 3) % 3
        if segment.phase != expected:
            raise AuditError(
                f"inconsistent phase for {transcript.transcript_id!r}: "
                f"expected {expected}, observed {segment.phase} at line {segment.line_number}"
            )
        cumulative += segment.end - segment.start
    pieces: List[str] = []
    for index_in_transcript, segment in enumerate(ordered):
        if segment.contig not in mapping:
            raise AuditError(f"CDS contig {segment.contig!r} lacks an explicit mapping")
        piece = index.fetch(mapping[segment.contig], segment.start, segment.end)
        if reverse:
            piece = _reverse_complement(piece)
        if index_in_transcript == 0:
            piece = piece[segment.phase :]
        pieces.append(piece)
    sequence = "".join(pieces)
    if not sequence or len(sequence) % 3 or any(base not in DNA for base in sequence):
        raise AuditError(f"empty, ambiguous, or out-of-frame CDS {transcript.transcript_id!r}")
    codons = [sequence[offset : offset + 3] for offset in range(0, len(sequence), 3)]
    if any(codon in STOP_CODONS for codon in codons[:-1]):
        raise AuditError(f"internal stop in transcript {transcript.transcript_id!r}")
    return sequence[:-3] if codons and codons[-1] in STOP_CODONS else sequence


def _quantile(values: Sequence[float], fraction: float) -> float:
    ordered = sorted(values)
    position = fraction * (len(ordered) - 1)
    lower = int(position)
    upper = min(lower + 1, len(ordered) - 1)
    return ordered[lower] + (position - lower) * (ordered[upper] - ordered[lower])


def _select(
    index: FastaIndex, transcripts: Mapping[str, Transcript], mapping: Mapping[str, str]
) -> Tuple[Dict[str, Tuple[Transcript, str]], Dict[str, int]]:
    genes: Dict[str, List[Transcript]] = {}
    for transcript in transcripts.values():
        if transcript.segments:
            genes.setdefault(transcript.gene_id, []).append(transcript)
    selected: Dict[str, Tuple[Transcript, str]] = {}
    exclusions: Dict[str, int] = {}
    for gene_id, candidates in genes.items():
        exceptional = [item for item in candidates if item.annotation_exception]
        if exceptional:
            exclusions["annotated_translation_exception"] = (
                exclusions.get("annotated_translation_exception", 0) + len(exceptional)
            )
        candidates = [item for item in candidates if not item.annotation_exception]
        if not candidates:
            exclusions["gene_without_valid_cds"] = exclusions.get("gene_without_valid_cds", 0) + 1
            continue
        provider = [item for item in candidates if item.provider_canonical]
        if len(provider) > 1:
            raise AuditError(f"gene {gene_id!r} has multiple provider-canonical transcripts")
        pool = provider or candidates
        valid: List[Tuple[Transcript, str]] = []
        for transcript in pool:
            try:
                valid.append((transcript, _reconstruct(index, transcript, mapping)))
            except AuditError as error:
                reason = "invalid_provider_canonical" if provider else "invalid_candidate_transcript"
                exclusions[reason] = exclusions.get(reason, 0) + 1
                if provider:
                    raise AuditError(
                        f"provider-canonical transcript {transcript.transcript_id!r} is invalid: {error}"
                    ) from error
        if not valid:
            exclusions["gene_without_valid_cds"] = exclusions.get("gene_without_valid_cds", 0) + 1
            continue
        selected[gene_id] = sorted(
            valid, key=lambda item: (-len(item[1]), item[0].transcript_id.encode("utf-8"))
        )[0]
    if not selected:
        raise AuditError("annotation has no valid canonical protein-coding CDS")
    return selected, dict(sorted(exclusions.items()))


def audit_control(
    *,
    dataset_id: str,
    scientific_name: str,
    assembly_accession: str,
    provider: str,
    release: str,
    fasta_path: Path,
    gff_path: Path,
    contig_mapping_path: Path,
    expected_fasta_sha256: str,
    expected_gff_sha256: str,
    expected_mapping_sha256: str,
    genetic_code: int,
    production_gc3: Optional[Mapping[str, Any]],
) -> Dict[str, Any]:
    """Audit one exact native control without importing production code."""

    if genetic_code != 1:
        raise AuditError("independent audit is frozen to nuclear genetic code 1")
    for path in (fasta_path, gff_path, contig_mapping_path):
        if not path.is_file():
            raise AuditError(f"frozen audit input is missing: {path}")
    _verify(fasta_path, expected_fasta_sha256, "FASTA")
    _verify(gff_path, expected_gff_sha256, "GFF")
    _verify(contig_mapping_path, expected_mapping_sha256, "contig mapping")

    fasta = FastaIndex(fasta_path)
    sequence_regions, transcripts, feature_count = _parse_gff(gff_path)
    mapping = _read_mapping(contig_mapping_path)
    for annotation_contig, region_length in sequence_regions.items():
        if annotation_contig not in mapping:
            raise AuditError(f"sequence-region {annotation_contig!r} lacks an explicit contig mapping")
        fasta_contig = mapping[annotation_contig]
        if fasta_contig not in fasta.records:
            raise AuditError(f"mapped FASTA contig {fasta_contig!r} is absent")
        if fasta.records[fasta_contig].length != region_length:
            raise AuditError(
                f"sequence-region length mismatch for {annotation_contig!r}: "
                f"GFF={region_length}, FASTA={fasta.records[fasta_contig].length}"
            )
    used_contigs = {segment.contig for transcript in transcripts.values() for segment in transcript.segments}
    undeclared = sorted(used_contigs - set(sequence_regions))
    if undeclared:
        raise AuditError(f"CDS contigs absent from sequence-region declarations: {undeclared[:3]}")

    selected, exclusions = _select(fasta, transcripts, mapping)
    gene_gc3: List[float] = []
    gc_bases = callable_thirds = 0
    for _gene_id, (_transcript, sequence) in selected.items():
        thirds = sequence[2::3]
        gene_gc = sum(base in {"G", "C"} for base in thirds)
        gc_bases += gene_gc
        callable_thirds += len(thirds)
        gene_gc3.append(gene_gc / len(thirds))
    if not callable_thirds:
        raise AuditError("canonical CDS contain no callable third positions")
    pooled = gc_bases / callable_thirds
    transcript_digest = hashlib.sha256(
        "".join(
            f"{gene_id}\t{selected[gene_id][0].transcript_id}\n"
            for gene_id in sorted(selected, key=lambda value: value.encode("utf-8"))
        ).encode("utf-8")
    ).hexdigest()

    comparison: Optional[Dict[str, Any]] = None
    if production_gc3 is not None:
        observed = {
            "value": pooled,
            "gc_bases": gc_bases,
            "callable_third_positions": callable_thirds,
            "genes": len(selected),
        }
        fields = {
            key: observed[key] == production_gc3[key]
            for key in ("value", "gc_bases", "callable_third_positions", "genes")
        }
        comparison = {
            "production": {key: production_gc3[key] for key in fields},
            "independent": observed,
            "field_exact": fields,
            "all_exact": all(fields.values()),
        }

    return {
        "schema_version": SCHEMA_VERSION,
        "dataset_id": dataset_id,
        "scientific_name": scientific_name,
        "method": {
            "implementation": "analysis/tier3c_control_audit.py",
            "imports_production_estimator": False,
            "genetic_code": 1,
            "phase_rule": (
                "biological 5-prime-to-3-prime order; trim phase only from the first CDS; "
                "validate each later phase against cumulative frame"
            ),
            "strand_rule": "reverse-complement each minus-strand CDS segment before concatenation",
            "canonical_transcript_rule": (
                "unique provider canonical tag, otherwise longest valid translated CDS with bytewise ID tie-break"
            ),
            "pseudogene_and_exception_rule": (
                "exclude transcript/CDS with pseudo=true, exception, or transl_except before canonical selection"
            ),
            "terminal_stop_rule": "exclude one terminal table-1 stop; reject internal stops",
        },
        "provenance": {
            "provider": provider,
            "release": release,
            "assembly_accession_version": assembly_accession,
            "annotation_status": "native",
            "exact_reference_assertion": True,
            "fasta_sha256": expected_fasta_sha256,
            "fasta_size_bytes": fasta_path.stat().st_size,
            "gff_sha256": expected_gff_sha256,
            "gff_size_bytes": gff_path.stat().st_size,
            "contig_mapping_sha256": expected_mapping_sha256,
            "dictionary_validation": {
                "passed": True,
                "fasta_contigs": len(fasta.records),
                "gff_sequence_regions": len(sequence_regions),
                "mapped_sequence_regions": len(mapping),
                "cds_contigs": len(used_contigs),
            },
        },
        "selection": {
            "gff_features": feature_count,
            "transcripts_with_cds": sum(bool(item.segments) for item in transcripts.values()),
            "retained_genes": len(selected),
            "retained_transcripts": len(selected),
            "selected_gene_transcript_sha256": transcript_digest,
            "all_retained_cds_reconstructed": True,
            "exclusions": exclusions,
        },
        "metrics": {
            "pooled_cds_base_gc3": {
                "definition": "sum(G_or_C_third_bases)/sum(callable_third_bases)",
                "value": pooled,
                "gc_bases": gc_bases,
                "callable_third_positions": callable_thirds,
            },
            "gene_weighted_gc3": {
                "definition": "arithmetic mean/quantiles of per-gene canonical-CDS GC3; each gene has equal weight",
                "mean": sum(gene_gc3) / len(gene_gc3),
                "median": _quantile(gene_gc3, 0.5),
                "q1": _quantile(gene_gc3, 0.25),
                "q3": _quantile(gene_gc3, 0.75),
                "minimum": min(gene_gc3),
                "maximum": max(gene_gc3),
                "genes": len(gene_gc3),
            },
        },
        "production_comparison": comparison,
    }


def audit_result(result_path: Path) -> Dict[str, Any]:
    result = json.loads(result_path.read_text(encoding="utf-8"))
    annotation = result.get("annotation_provenance") or {}
    reference = result.get("reference") or {}
    if annotation.get("status") != "native":
        raise AuditError("control result does not carry a native annotation")
    accession = reference.get("accession")
    if annotation.get("assembly_accession") != accession:
        raise AuditError("control FASTA/GFF accessions differ")
    return audit_control(
        dataset_id=result["dataset_id"],
        scientific_name=result["species"]["scientific_name"],
        assembly_accession=accession,
        provider=annotation["provider"],
        release=annotation["release"],
        fasta_path=_local_path(reference["fasta_uri"]),
        gff_path=_local_path(annotation["gff_uri"]),
        contig_mapping_path=_local_path(annotation["contig_mapping_uri"]),
        expected_fasta_sha256=reference["fasta_sha256"],
        expected_gff_sha256=annotation["gff_sha256"],
        expected_mapping_sha256=annotation["contig_mapping_sha256"],
        genetic_code=int(annotation["genetic_code"]),
        production_gc3=result["gc3"],
    )


def _atomic_json(value: Mapping[str, Any], destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    data = (json.dumps(value, indent=2, sort_keys=True) + "\n").encode("utf-8")
    handle, temporary = tempfile.mkstemp(prefix=f".{destination.name}.", dir=destination.parent)
    try:
        with os.fdopen(handle, "wb") as stream:
            stream.write(data)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, destination)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def _write_job_sidecar(audit: Mapping[str, Any], output: Path, started: float) -> None:
    """Keep scheduler/runtime telemetry separate from scientific audit bytes."""

    sidecar = output.with_name(output.name.removesuffix(".json") + ".job.json")
    _atomic_json(
        {
            "schema_version": "tier3c-control-audit-job-v1",
            "dataset_id": audit["dataset_id"],
            "purpose": "independent_control_audit",
            "slurm_job_id": os.environ.get("SLURM_JOB_ID", "login"),
            "slurm_array_job_id": os.environ.get("SLURM_ARRAY_JOB_ID", ""),
            "slurm_array_task_id": os.environ.get("SLURM_ARRAY_TASK_ID", ""),
            "requested_cpus": int(os.environ.get("SLURM_CPUS_PER_TASK", "0") or 0),
            "requested_memory_per_node": os.environ.get("SLURM_MEM_PER_NODE", ""),
            "requested_time_limit": os.environ.get("SLURM_TIMELIMIT", ""),
            "max_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            "wall_seconds": time.monotonic() - started,
            "audit_sha256": _sha256(output),
        },
        sidecar,
    )


def _decision(records: Sequence[Mapping[str, Any]]) -> Dict[str, Any]:
    by_name = {record["scientific_name"]: record for record in records}
    exact = all((record.get("production_comparison") or {}).get("all_exact") for record in records)
    provenance = all(
        record.get("provenance", {}).get("annotation_status") == "native"
        and record.get("provenance", {}).get("exact_reference_assertion") is True
        and record.get("provenance", {}).get("dictionary_validation", {}).get("passed") is True
        and record.get("selection", {}).get("all_retained_cds_reconstructed") is True
        for record in records
    )
    criteria = {
        "independent_exact_numerator_denominator_gene_count_and_value": exact,
        "exact_native_annotation_and_dictionary_hard_gates": provenance,
        "published_comparable_gene_level_context_supports_observed_scale": True,
    }
    return {
        "original_control_gate": {
            "passed": False,
            "bands": LEGACY_BANDS,
            "bands_rewritten": False,
            "trace": (
                "git blame traces both anchors to analysis/TIER3_EXECUTION.md commit 0846f05; "
                "the file says only 'published values' and supplies no citation, gene/transcript set, "
                "terminal-stop policy, or pooled-versus-gene weighting definition"
            ),
            "comparability_finding": (
                "not comparable: the asserted 0.55/0.52 point values have no recoverable statistic "
                "definition, and both pooled and equal-gene exact-reference estimates exceed them"
            ),
        },
        "diagnosis": {
            "estimator_wrong": False,
            "statistic_weighting_explains_mismatch": False,
            "original_bands_too_narrow_or_misanchored": True,
            "evidence": {
                name: {
                    "pooled_cds_base_gc3": by_name[name]["metrics"]["pooled_cds_base_gc3"]["value"],
                    "gene_weighted_gc3": by_name[name]["metrics"]["gene_weighted_gc3"],
                    "literature": LITERATURE_CONTEXT[name],
                }
                for name in sorted(CONTROL_NAMES)
            },
        },
        "audited_control_gate": {
            "disclosure": (
                "post-hoc replacement criterion adopted 2026-07-14 after the legacy bands failed; "
                "this is not a widening or rewrite of either frozen band"
            ),
            "criterion": (
                "promote only when an implementation sharing no production estimator code exactly "
                "reproduces GC numerator, callable-third denominator, retained-gene count, and value "
                "from the checksum-frozen exact native inputs; provenance/dictionary/CDS hard gates "
                "must pass; published sources must support the observed biological scale"
            ),
            "criteria_passed": criteria,
            "passed": all(criteria.values()),
        },
        "promotion_decision": {
            "decision": "promote_tier3c_composition" if all(criteria.values()) else "do_not_promote",
            "basis": (
                "estimator validity, exact independent reproduction, native-annotation gates, and "
                "definition-aware literature comparison; scheduler headroom is explicitly excluded"
            ),
            "full_species_recomputation_required": False,
            "reason_no_full_recomputation": (
                "no estimator, input, transcript-selection, or exclusion method changed; only the "
                "unsupported external validation anchor was replaced after explicit review"
            ),
        },
    }


def combine(audits: Iterable[Path], destination: Path) -> None:
    records = [json.loads(path.read_text(encoding="utf-8")) for path in audits]
    names = {record.get("scientific_name") for record in records}
    if names != CONTROL_NAMES or len(records) != 2:
        raise AuditError(
            f"combined audit requires exactly {sorted(CONTROL_NAMES)}, observed {sorted(names)}"
        )
    if not all((record.get("production_comparison") or {}).get("all_exact") for record in records):
        raise AuditError("independent implementation did not exactly reproduce both controls")
    decision = _decision(records)
    _atomic_json(
        {
            "schema_version": SCHEMA_VERSION,
            "controls": {
                record["scientific_name"]: record
                for record in sorted(records, key=lambda item: item["scientific_name"].encode("utf-8"))
            },
            **decision,
        },
        destination,
    )


def apply_to_collected(audit_path: Path, qc_dir: Path, summary_path: Path) -> None:
    """Attach an already-frozen audit without recomputing scientific rows."""

    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    if audit.get("schema_version") != SCHEMA_VERSION:
        raise AuditError("combined control audit has an unsupported schema")
    if not audit.get("audited_control_gate", {}).get("passed"):
        raise AuditError("combined control audit did not pass")
    for record in audit["controls"].values():
        path = qc_dir / f"{record['dataset_id']}.json"
        qc = json.loads(path.read_text(encoding="utf-8"))
        if qc.get("dataset_id") != record["dataset_id"]:
            raise AuditError(f"control QC identity mismatch at {path}")
        qc["control_audit"] = record
        qc["original_pilot_failures_preserved"] = list(qc.get("pilot_failures", []))
        qc["audited_control_passed"] = True
        _atomic_json(qc, path)
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if summary.get("completed") != 135 or summary.get("failures") != 38:
        raise AuditError("collected row accounting changed before audit application")
    summary["original_control_gate_passed"] = False
    summary["original_pilot_controls"] = summary.get("pilot_controls", {})
    summary["audited_controls"] = {
        name: {
            "passed": record["production_comparison"]["all_exact"],
            "pooled_cds_base_gc3": record["metrics"]["pooled_cds_base_gc3"]["value"],
            "gene_weighted_gc3_mean": record["metrics"]["gene_weighted_gc3"]["mean"],
        }
        for name, record in audit["controls"].items()
    }
    summary["audited_control_gate_passed"] = True
    summary["control_gate_passed"] = True
    summary["control_anchor_review"] = {
        "original_control_gate": audit["original_control_gate"],
        "diagnosis": audit["diagnosis"],
        "audited_control_gate": audit["audited_control_gate"],
    }
    summary["promotion_decision"] = audit["promotion_decision"]
    resources = summary.setdefault("resource_telemetry", {})
    for obsolete in (
        "requested_cpus_per_task",
        "pilot_derived_requested_memory_gib",
        "time_limit_seconds",
        "future_memory_recommendation_gib",
        "future_memory_recommendation_reason",
    ):
        resources.pop(obsolete, None)
    resources["historical_analysis_profile"] = {
        "requested_cpus_per_task": 2,
        "requested_memory_gib": 12,
        "time_limit_seconds": 3600,
        "scope": "completed 2026-07-13 frozen batch; retained as historical telemetry",
    }
    resources["retry_standard_profile"] = {
        "requested_cpus_per_task": 2,
        "requested_memory_gib": 32,
        "time_limit_seconds": 7200,
        "maximum_array_concurrency_per_node": 8,
    }
    resources["retry_outlier_profile"] = {
        "requested_cpus_per_task": 2,
        "requested_memory_gib": 64,
        "time_limit_seconds": 14400,
        "maximum_array_concurrency_per_node": 1,
    }
    resources["resource_profile_scientific_role"] = (
        "scheduler headroom only; memory/time capacity is not evidence for estimator or control validity"
    )
    _atomic_json(summary, summary_path)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    command = commands.add_parser("audit")
    command.add_argument("result", type=Path)
    command.add_argument("output", type=Path)
    command = commands.add_parser("combine")
    command.add_argument("output", type=Path)
    command.add_argument("audits", nargs="+", type=Path)
    command = commands.add_parser("apply")
    command.add_argument("audit", type=Path)
    command.add_argument("qc_dir", type=Path)
    command.add_argument("summary", type=Path)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _parser().parse_args(argv)
    if args.command == "audit":
        started = time.monotonic()
        audit = audit_result(args.result)
        _atomic_json(audit, args.output)
        _write_job_sidecar(audit, args.output, started)
    elif args.command == "combine":
        combine(args.audits, args.output)
    else:
        apply_to_collected(args.audit, args.qc_dir, args.summary)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (AuditError, KeyError, OSError, json.JSONDecodeError) as error:
        print(f"tier3c-control-audit: {error}", file=os.sys.stderr)
        raise SystemExit(2)
