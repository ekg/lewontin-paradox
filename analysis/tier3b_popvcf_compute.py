#!/usr/bin/env python3
"""Compute Tier 3b population diversity from pre-called VCF/BCF resources.

The central safety property is that a sparse VCF never supplies its own
invariant denominator.  Callers must declare an all-sites/gVCF source or pass
an exact-selected-cohort callable BED.  Coordinates are zero-based internally.
"""

from __future__ import annotations

import argparse
import bisect
import collections
import gzip
import hashlib
import json
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, Iterator, List, Mapping, Optional, Sequence, Set, Tuple

if __package__ in (None, ""):  # permit ``python analysis/<script>.py``
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from analysis.tier3_common import (
    DNA,
    Tier3ValidationError,
    allele_pairwise_diversity,
    collect_fourfold_sites,
    fasta_dictionary,
    merge_intervals,
    parse_gff,
    read_fasta,
    reconstruct_cds,
    resolve_contig_aliases,
    sha256_file,
)


POLICY_ID = "tier3-decisions-v1"
BLOCK_SIZE_BP = 1_000_000
BOOTSTRAP_REPLICATES = 10_000
MINIMUM_BOOTSTRAP_BLOCKS = 20
MINIMUM_4D_CLASS_SITES = 10_000
DENOMINATOR_KINDS = frozenset(("all_sites_vcf", "gvcf", "cohort_callable_mask"))


@dataclass(frozen=True)
class VariantRecord:
    contig: str
    position: int
    ref: str
    alts: Tuple[str, ...]
    filters: Tuple[str, ...]
    info: Mapping[str, str]
    genotypes: Mapping[str, Tuple[Optional[int], ...]]


def read_vcf_header(path: Path) -> Tuple[Dict[str, int], Tuple[str, ...]]:
    """Read only contig/sample headers, including from very large BCFs."""

    path = Path(path)
    if path.suffix == ".bcf":
        try:
            import pysam
        except ImportError as error:
            raise Tier3ValidationError("pysam is required to read BCF in the pinned Guix environment") from error
        with pysam.VariantFile(str(path)) as handle:
            contigs = {
                name: int(record.length)
                for name, record in handle.header.contigs.items()
                if record.length is not None
            }
            return contigs, tuple(handle.header.samples)
    contigs: Dict[str, int] = {}
    with _open_text(path) as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith("##contig="):
                parsed = _parse_contig_header(line)
                assert parsed is not None
                if parsed[0] in contigs:
                    raise Tier3ValidationError("duplicate VCF contig declaration {!r}".format(parsed[0]))
                contigs[parsed[0]] = parsed[1]
            elif line.startswith("#CHROM\t"):
                samples = tuple(line.split("\t")[9:])
                if len(samples) != len(set(samples)):
                    raise Tier3ValidationError("duplicate VCF sample IDs")
                return contigs, samples
    raise Tier3ValidationError("VCF has no #CHROM header")


def iter_vcf_records(path: Path, samples: Optional[Sequence[str]] = None) -> Iterator[VariantRecord]:
    """Stream VCF/BCF records without retaining a population file in memory."""

    path = Path(path)
    if samples is None:
        _contigs, samples = read_vcf_header(path)
    expected_samples = tuple(samples)
    if path.suffix == ".bcf":
        try:
            import pysam
        except ImportError as error:
            raise Tier3ValidationError("pysam is required to read BCF in the pinned Guix environment") from error
        with pysam.VariantFile(str(path)) as handle:
            if tuple(handle.header.samples) != expected_samples:
                raise Tier3ValidationError("BCF samples changed between header and record pass")
            for record in handle:
                info = {key: str(value) for key, value in record.info.items()}
                alts = tuple(record.alts or ())
                if any(alt.startswith("<") for alt in alts) and record.stop > record.start + len(record.ref):
                    info.setdefault("END", str(record.stop))
                yield VariantRecord(
                    record.contig,
                    record.start,
                    record.ref,
                    alts,
                    tuple(key for key in record.filter.keys() if key not in ("PASS", ".")),
                    info,
                    {
                        sample: tuple(record.samples[sample].get("GT") or (None,))
                        for sample in expected_samples
                    },
                )
        return

    observed_samples: Optional[Tuple[str, ...]] = None
    with _open_text(path) as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.rstrip("\n")
            if line.startswith("##") or not line:
                continue
            if line.startswith("#CHROM\t"):
                observed_samples = tuple(line.split("\t")[9:])
                if observed_samples != expected_samples:
                    raise Tier3ValidationError("VCF samples changed between header and record pass")
                continue
            if line.startswith("#"):
                continue
            if observed_samples is None:
                raise Tier3ValidationError("VCF data precede #CHROM header")
            fields = line.split("\t")
            if len(fields) < 8:
                raise Tier3ValidationError("VCF line {} has fewer than eight fields".format(line_number))
            contig, pos_text, _identifier, ref, alt_text, _qual, filter_text, info_text = fields[:8]
            position = int(pos_text) - 1
            if position < 0:
                raise Tier3ValidationError("invalid VCF POS at line {}".format(line_number))
            alts = () if alt_text in ("", ".") else tuple(alt_text.split(","))
            filters = () if filter_text in ("", ".", "PASS") else tuple(filter_text.split(";"))
            info: Dict[str, str] = {}
            if info_text not in ("", "."):
                for item in info_text.split(";"):
                    key, separator, value = item.partition("=")
                    info[key] = value if separator else "true"
            genotypes: Dict[str, Tuple[Optional[int], ...]] = {}
            if expected_samples:
                if len(fields) != 9 + len(expected_samples):
                    raise Tier3ValidationError("VCF sample column count mismatch at line {}".format(line_number))
                formats = fields[8].split(":")
                if "GT" not in formats:
                    raise Tier3ValidationError("VCF record lacks GT at line {}".format(line_number))
                gt_index = formats.index("GT")
                for sample, sample_text in zip(expected_samples, fields[9:]):
                    values = sample_text.split(":")
                    genotypes[sample] = _parse_gt_text(values[gt_index] if gt_index < len(values) else ".")
            yield VariantRecord(contig, position, ref, alts, filters, info, genotypes)


def _open_text(path: Path):
    if str(path).endswith((".gz", ".bgz")):
        return gzip.open(str(path), "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def _parse_contig_header(line: str) -> Optional[Tuple[str, int]]:
    match = re.match(r"##contig=<(.+)>$", line)
    if not match:
        return None
    values: Dict[str, str] = {}
    for field in re.split(r",(?=[A-Za-z][A-Za-z0-9_]*=)", match.group(1)):
        if "=" in field:
            key, value = field.split("=", 1)
            values[key] = value.strip('"')
    if "ID" not in values or "length" not in values:
        raise Tier3ValidationError("every VCF ##contig declaration must contain ID and length")
    return values["ID"], int(values["length"])


def _parse_gt_text(text: str) -> Tuple[Optional[int], ...]:
    if text in ("", "."):
        return (None,)
    fields = re.split(r"[/|]", text)
    try:
        return tuple(None if value == "." else int(value) for value in fields)
    except ValueError as error:
        raise Tier3ValidationError("malformed GT value {!r}".format(text)) from error


def read_vcf(path: Path) -> Tuple[Dict[str, int], Tuple[str, ...], List[VariantRecord]]:
    """Materialize records for small callers; production compute streams them."""

    contigs, samples = read_vcf_header(path)
    return contigs, samples, list(iter_vcf_records(path, samples))


def read_callable_bed(path: Path, contig_lengths: Mapping[str, int]) -> Dict[str, List[Tuple[int, int]]]:
    intervals: Dict[str, List[Tuple[int, int]]] = collections.defaultdict(list)
    with _open_text(Path(path)) as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise Tier3ValidationError("callable BED line {} has fewer than three fields".format(line_number))
            contig = fields[0]
            if contig not in contig_lengths:
                raise Tier3ValidationError("callable BED contig {!r} is absent from FASTA".format(contig))
            start, end = int(fields[1]), int(fields[2])
            if not 0 <= start < end <= contig_lengths[contig]:
                raise Tier3ValidationError("callable BED line {} is outside the FASTA".format(line_number))
            intervals[contig].append((start, end))
    if not intervals:
        raise Tier3ValidationError("callable BED has no intervals")
    return {contig: merge_intervals(values) for contig, values in intervals.items()}


def _position_in_intervals(position: int, intervals: Sequence[Tuple[int, int]]) -> bool:
    # BED intervals are typically few compared with callable bases.  Binary
    # search avoids materializing whole-genome position sets.
    low, high = 0, len(intervals)
    while low < high:
        middle = (low + high) // 2
        start, end = intervals[middle]
        if position < start:
            high = middle
        elif position >= end:
            low = middle + 1
        else:
            return True
    return False


def _coerce_genotype(genotype: Sequence[Optional[int]], design: str) -> Tuple[Optional[int], ...]:
    genotype = tuple(genotype)
    if design == "wild_diploid":
        if len(genotype) != 2:
            raise Tier3ValidationError("wild-diploid genotype does not have ploidy two")
        return genotype
    if design == "inbred_lines_haploidized":
        if len(genotype) == 1:
            return genotype
        if len(genotype) == 2 and None not in genotype and genotype[0] == genotype[1]:
            return (genotype[0],)
        # Heterozygous and partially missing inbred calls are uncertain, hence missing.
        return (None,)
    if design == "haploid":
        if len(genotype) != 1:
            raise Tier3ValidationError("haploid genotype does not have ploidy one")
        return genotype
    raise Tier3ValidationError("unsupported population design {!r}".format(design))


def _validate_records(
    records: Sequence[VariantRecord],
    vcf_contigs: Mapping[str, int],
    fasta: Mapping[str, str],
    samples: Sequence[str],
) -> None:
    if set(vcf_contigs) != set(fasta):
        raise Tier3ValidationError("VCF and FASTA contig dictionaries differ")
    for contig, length in vcf_contigs.items():
        if length != len(fasta[contig]):
            raise Tier3ValidationError("VCF/FASTA contig length mismatch for {!r}".format(contig))
    seen: Set[Tuple[str, int]] = set()
    for record in records:
        if record.contig not in fasta:
            raise Tier3ValidationError("VCF contig {!r} is absent from FASTA".format(record.contig))
        if (record.contig, record.position) in seen:
            raise Tier3ValidationError(
                "multiple normalized records at {}:{}; retain one multiallelic site record".format(
                    record.contig, record.position + 1
                )
            )
        seen.add((record.contig, record.position))
        observed = fasta[record.contig][record.position : record.position + len(record.ref)].upper()
        if observed != record.ref.upper():
            raise Tier3ValidationError(
                "VCF REF mismatch at {}:{} ({} != {})".format(
                    record.contig, record.position + 1, record.ref, observed
                )
            )
        for sample in samples:
            if sample not in record.genotypes:
                raise Tier3ValidationError("selected sample {!r} is missing from a VCF record".format(sample))
            for allele in record.genotypes[sample]:
                if allele is not None and not 0 <= allele <= len(record.alts):
                    raise Tier3ValidationError("genotype allele index exceeds ALT count")


def _load_contig_mapping(path: Optional[Path]) -> Dict[str, str]:
    if path is None:
        return {}
    mapping: Dict[str, str] = {}
    with Path(path).open("r", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if fields[:2] in (["annotation_contig", "fasta_contig"], ["gff_contig", "fasta_contig"]):
                continue
            if len(fields) != 2 or not fields[0] or not fields[1]:
                raise Tier3ValidationError("invalid contig mapping line {}".format(line_number))
            if fields[0] in mapping:
                raise Tier3ValidationError("duplicate annotation contig mapping {!r}".format(fields[0]))
            mapping[fields[0]] = fields[1]
    return mapping


def audit_native_annotation(
    dataset_id: str,
    fasta_path: Path,
    gff_path: Path,
    annotation_metadata: Mapping[str, Any],
    contig_mapping_path: Optional[Path] = None,
    provider_cds_path: Optional[Path] = None,
) -> Tuple[Set[Tuple[str, int]], Dict[str, Any]]:
    """Validate the hard native/exact-reference gate and return FASTA 4D sites."""

    required = ("provider", "release", "assembly_accession", "fasta_assembly_accession", "status", "genetic_code")
    missing = [key for key in required if key not in annotation_metadata]
    if missing:
        raise Tier3ValidationError("annotation provenance lacks {}".format(", ".join(missing)))
    if annotation_metadata["status"] != "native":
        raise Tier3ValidationError("primary 4D analysis requires native annotation")
    if annotation_metadata["assembly_accession"] != annotation_metadata["fasta_assembly_accession"]:
        raise Tier3ValidationError("FASTA/GFF assembly accession+version mismatch")
    genetic_code = int(annotation_metadata["genetic_code"])
    if genetic_code != 1:
        raise Tier3ValidationError("only NCBI nuclear genetic code 1 is frozen")

    fasta = read_fasta(fasta_path)
    annotation = parse_gff(gff_path)
    declared_mapping = _load_contig_mapping(contig_mapping_path)
    resolved = resolve_contig_aliases(fasta_dictionary(fasta), annotation.sequence_regions, declared_mapping)
    annotation_fasta = dict(fasta)
    for annotation_contig, fasta_contig in resolved.items():
        annotation_fasta[annotation_contig] = fasta[fasta_contig]

    # Reconstruct every retained transcript; collect_fourfold_sites additionally
    # checks frame, ambiguity, translation, stops, and overlap disagreement.
    retained_ids = sorted(set(annotation.canonical_transcripts.values()))
    reconstructed = {
        transcript_id: reconstruct_cds(annotation_fasta, annotation.transcripts[transcript_id])
        for transcript_id in retained_ids
    }
    fourfold_annotation, overlap_exclusions = collect_fourfold_sites(annotation_fasta, annotation, genetic_code)
    fourfold = {(resolved.get(contig, contig), position) for contig, position in fourfold_annotation}

    sample_ids = sorted(
        retained_ids,
        key=lambda transcript_id: hashlib.sha256(
            (dataset_id + transcript_id).encode("utf-8")
        ).digest(),
    )[:100]
    provider_comparison = "not_deposited"
    mismatches = 0
    if provider_cds_path is not None:
        provider_sequences = read_fasta(provider_cds_path)
        provider_comparison = "passed"
        for transcript_id in sample_ids:
            if provider_sequences.get(transcript_id, "").upper() != reconstructed[transcript_id]:
                mismatches += 1
        if mismatches:
            raise Tier3ValidationError("sampled provider CDS reconstruction mismatch")

    mapping_bytes = "".join(
        "{}\t{}\n".format(source, target) for source, target in sorted(resolved.items())
    ).encode("utf-8")
    audit = {
        "provider": str(annotation_metadata["provider"]),
        "release": str(annotation_metadata["release"]),
        "assembly_accession": str(annotation_metadata["assembly_accession"]),
        "fasta_assembly_accession": str(annotation_metadata["fasta_assembly_accession"]),
        "status": "native",
        "genetic_code": genetic_code,
        "fasta_sha256": sha256_file(fasta_path),
        "fasta_size_bytes": Path(fasta_path).stat().st_size,
        "gff_sha256": sha256_file(gff_path),
        "gff_size_bytes": Path(gff_path).stat().st_size,
        "provider_cds_sha256": sha256_file(provider_cds_path) if provider_cds_path else None,
        "provider_cds_size_bytes": Path(provider_cds_path).stat().st_size if provider_cds_path else None,
        "sequence_regions": dict(sorted(annotation.sequence_regions.items())),
        "fasta_contig_dictionary": fasta_dictionary(fasta),
        "contig_mapping": dict(sorted(resolved.items())),
        "contig_mapping_sha256": hashlib.sha256(mapping_bytes).hexdigest(),
        "contig_dictionary_passed": True,
        "sample_rule": "first_100_by_sha256_dataset_id_plus_transcript_id_or_all_if_fewer",
        "sampled_cds_count": len(sample_ids),
        "sampled_cds_ids_sha256": hashlib.sha256(
            "".join(identifier + "\n" for identifier in sample_ids).encode("utf-8")
        ).hexdigest(),
        "sampled_cds_provider_comparison": provider_comparison,
        "sampled_cds_mismatches": mismatches,
        "all_retained_cds_phase_translation_passed": True,
        "retained_transcripts": len(retained_ids),
        "fourfold_sites": len(fourfold),
        "overlap_frame_exclusions": len(overlap_exclusions),
    }
    return fourfold, audit


def _interval_block_counts(
    intervals: Mapping[str, Sequence[Tuple[int, int]]],
    fasta: Mapping[str, str],
) -> Tuple[Dict[Tuple[str, int], int], int]:
    counts: Dict[Tuple[str, int], int] = collections.defaultdict(int)
    total = 0
    for contig, values in intervals.items():
        sequence = fasta[contig]
        for start, end in values:
            position = start
            while position < end:
                block = position // BLOCK_SIZE_BP
                boundary = min(end, (block + 1) * BLOCK_SIZE_BP)
                segment = sequence[position:boundary].upper()
                acgt = sum(segment.count(base) for base in "ACGT")
                counts[(contig, block)] += acgt
                total += acgt
                position = boundary
    return dict(counts), total


def _explicit_denominator_intervals(
    records: Iterable[VariantRecord], denominator_kind: str
) -> Dict[str, List[Tuple[int, int]]]:
    intervals: Dict[str, List[Tuple[int, int]]] = collections.defaultdict(list)
    invariant_evidence = False
    previous_contig: Optional[str] = None
    previous_position = -1
    finished_contigs: Set[str] = set()
    for record in records:
        if record.contig == previous_contig and record.position < previous_position:
            raise Tier3ValidationError("all-sites/gVCF records must be coordinate sorted")
        if record.contig != previous_contig:
            if previous_contig is not None:
                finished_contigs.add(previous_contig)
            if record.contig in finished_contigs:
                raise Tier3ValidationError("all-sites/gVCF contig records are not contiguous")
            previous_contig = record.contig
        previous_position = record.position
        end = record.position + 1
        symbolic = any(alt.startswith("<") for alt in record.alts)
        if denominator_kind == "gvcf" and "END" in record.info:
            end = int(record.info["END"])
            invariant_evidence = True
        if not record.alts or symbolic:
            invariant_evidence = True
        contig_intervals = intervals[record.contig]
        if contig_intervals and record.position <= contig_intervals[-1][1]:
            contig_intervals[-1] = (contig_intervals[-1][0], max(contig_intervals[-1][1], end))
        else:
            contig_intervals.append((record.position, end))
    if not invariant_evidence:
        raise Tier3ValidationError(
            "declared {} contains no invariant/reference records; variant-only input has no callability".format(
                denominator_kind
            )
        )
    return dict(intervals)


def _seed_digest(dataset_id: str, statistic_name: str) -> str:
    return hashlib.sha256(
        (POLICY_ID + "\0" + dataset_id + "\0" + statistic_name).encode("utf-8")
    ).hexdigest()


def _counter_index(seed: bytes, replicate: int, contig: str, draw: int, upper: int) -> int:
    payload = (
        seed
        + replicate.to_bytes(8, "big")
        + len(contig.encode("utf-8")).to_bytes(4, "big")
        + contig.encode("utf-8")
        + draw.to_bytes(8, "big")
    )
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "big") % upper


def _percentile(values: Sequence[float], percentile: float) -> Optional[float]:
    finite = sorted(value for value in values if math.isfinite(value))
    if not finite:
        return None
    rank = (len(finite) - 1) * percentile
    lower = int(math.floor(rank))
    upper = int(math.ceil(rank))
    if lower == upper:
        return finite[lower]
    return finite[lower] * (upper - rank) + finite[upper] * (rank - lower)


def _bootstrap(
    dataset_id: str,
    statistic_name: str,
    block_stats: Mapping[Tuple[str, int], Mapping[str, float]],
    replicates: int,
    ratio: bool,
    eligible: bool,
) -> Dict[str, Any]:
    digest = _seed_digest(dataset_id, statistic_name)
    by_contig: Dict[str, List[Mapping[str, float]]] = collections.defaultdict(list)
    for (contig, _block), stats in sorted(block_stats.items()):
        if (stats["S_den"] + stats["W_den"] if ratio else stats["den"]) > 0:
            by_contig[contig].append(stats)
    block_count = sum(len(values) for values in by_contig.values())
    result: Dict[str, Any] = {
        "block_size_bp": BLOCK_SIZE_BP,
        "eligible_blocks": block_count,
        "minimum_blocks": MINIMUM_BOOTSTRAP_BLOCKS,
        "replicates": replicates,
        "seed_digest": digest,
        "seed_first_64_bits_big_endian": int(digest[:16], 16),
        "rng": "sha256-counter-v1",
        "interval": None,
        "unavailable_reason": None,
    }
    if not eligible:
        result["unavailable_reason"] = "statistic_unavailable"
        return result
    if block_count < MINIMUM_BOOTSTRAP_BLOCKS:
        result["unavailable_reason"] = "fewer_than_20_eligible_blocks"
        return result
    # The frozen seed is the first 64 digest bits; the complete digest remains
    # in provenance so independent implementations can reproduce it exactly.
    seed = bytes.fromhex(digest[:16])
    estimates: List[float] = []
    for replicate in range(replicates):
        totals: Dict[str, float] = collections.defaultdict(float)
        for contig, blocks in sorted(by_contig.items()):
            for draw in range(len(blocks)):
                selected = blocks[_counter_index(seed, replicate, contig, draw, len(blocks))]
                for key, value in selected.items():
                    totals[key] += value
        if ratio:
            pi_s = totals["S_num"] / totals["S_den"] if totals["S_den"] else float("nan")
            pi_w = totals["W_num"] / totals["W_den"] if totals["W_den"] else float("nan")
            estimates.append(pi_s / pi_w if pi_w > 0 else float("nan"))
        else:
            estimates.append(totals["num"] / totals["den"] if totals["den"] else float("nan"))
    low, high = _percentile(estimates, 0.025), _percentile(estimates, 0.975)
    if low is None or high is None:
        result["unavailable_reason"] = "no_finite_bootstrap_replicates"
    else:
        result["interval"] = [low, high]
    return result


def compute_population_pi(
    *,
    dataset_id: str,
    vcf_path: Path,
    fasta_path: Path,
    selected_samples: Sequence[str],
    design: str,
    denominator_kind: str,
    callable_bed_path: Optional[Path] = None,
    gff_path: Optional[Path] = None,
    annotation_metadata: Optional[Mapping[str, Any]] = None,
    contig_mapping_path: Optional[Path] = None,
    provider_cds_path: Optional[Path] = None,
    bootstrap_replicates: int = BOOTSTRAP_REPLICATES,
    minimum_4d_class_sites: int = MINIMUM_4D_CLASS_SITES,
    polarization_gate: Optional[Mapping[str, Any]] = None,
) -> Dict[str, Any]:
    """Compute a policy-complete Tier 3b result.

    ``minimum_4d_class_sites`` is injectable for synthetic tests; production
    CLI use is fixed to 10,000.  Polarized SFS-B remains unavailable in policy
    v1 even if a caller supplies a proposed gate.
    """

    if denominator_kind not in DENOMINATOR_KINDS:
        raise Tier3ValidationError("explicit denominator kind is required")
    if not selected_samples or len(selected_samples) != len(set(selected_samples)):
        raise Tier3ValidationError("selected sample list must be non-empty and unique")
    if design in ("wild_diploid", "inbred_lines_haploidized") and len(selected_samples) != 20:
        raise Tier3ValidationError("primary population designs require exactly 20 sampling units")
    if design == "wild_diploid":
        nominal_chromosomes, minimum_called = 2 * len(selected_samples), int(math.ceil(1.8 * len(selected_samples)))
    elif design in ("inbred_lines_haploidized", "haploid"):
        nominal_chromosomes, minimum_called = len(selected_samples), int(math.ceil(0.9 * len(selected_samples)))
    else:
        raise Tier3ValidationError("unsupported design {!r}".format(design))

    fasta = read_fasta(fasta_path)
    vcf_contigs, vcf_samples = read_vcf_header(vcf_path)
    absent_samples = sorted(set(selected_samples) - set(vcf_samples))
    if absent_samples:
        raise Tier3ValidationError("selected samples absent from VCF: {!r}".format(absent_samples))
    if set(vcf_contigs) != set(fasta):
        raise Tier3ValidationError("VCF and FASTA contig dictionaries differ")
    for contig, length in vcf_contigs.items():
        if length != len(fasta[contig]):
            raise Tier3ValidationError("VCF/FASTA contig length mismatch for {!r}".format(contig))

    if denominator_kind == "cohort_callable_mask":
        if callable_bed_path is None:
            raise Tier3ValidationError("cohort_callable_mask requires an explicit BED")
        intervals = read_callable_bed(callable_bed_path, fasta_dictionary(fasta))
    else:
        if callable_bed_path is not None:
            raise Tier3ValidationError("callable BED conflicts with {} denominator".format(denominator_kind))
        intervals = _explicit_denominator_intervals(
            iter_vcf_records(vcf_path, vcf_samples), denominator_kind
        )
        for contig, values in intervals.items():
            if contig not in fasta:
                raise Tier3ValidationError("denominator contig is absent from FASTA")
            if any(start < 0 or end <= start or end > len(fasta[contig]) for start, end in values):
                raise Tier3ValidationError("all-sites/gVCF denominator interval is outside the FASTA")

    fourfold: Set[Tuple[str, int]] = set()
    annotation_audit: Optional[Dict[str, Any]] = None
    if gff_path is not None or annotation_metadata is not None:
        if gff_path is None or annotation_metadata is None:
            raise Tier3ValidationError("4D analysis requires both GFF and annotation provenance")
        fourfold, annotation_audit = audit_native_annotation(
            dataset_id,
            fasta_path,
            gff_path,
            annotation_metadata,
            contig_mapping_path,
            provider_cds_path,
        )

    overall_blocks, initial_callable = _interval_block_counts(intervals, fasta)
    block_stats: Dict[Tuple[str, int], Dict[str, float]] = {}
    for block, count in overall_blocks.items():
        block_stats[block] = {
            "num": 0.0, "den": float(count),
            "S_num": 0.0, "S_den": 0.0, "W_num": 0.0, "W_den": 0.0,
        }
    for contig, position in fourfold:
        if _position_in_intervals(position, intervals.get(contig, ())):
            base_class = "S" if fasta[contig][position].upper() in "GC" else "W"
            block_stats[(contig, position // BLOCK_SIZE_BP)][base_class + "_den"] += 1.0
    fourfold_by_contig: Dict[str, List[int]] = collections.defaultdict(list)
    for contig, position in sorted(fourfold):
        fourfold_by_contig[contig].append(position)

    exclusion_counts: Dict[str, int] = collections.Counter()
    called_histogram: Dict[int, int] = collections.Counter()
    informative_snvs = 0
    contig_order = {contig: index for index, contig in enumerate(vcf_contigs)}
    previous_coordinate: Optional[Tuple[int, int]] = None
    # Records that fail QC must be removed from an otherwise callable mask.
    for record in iter_vcf_records(vcf_path, vcf_samples):
        if record.contig not in fasta:
            raise Tier3ValidationError("VCF contig {!r} is absent from FASTA".format(record.contig))
        coordinate = (contig_order[record.contig], record.position)
        if previous_coordinate is not None and coordinate == previous_coordinate:
            raise Tier3ValidationError(
                "multiple normalized records at {}:{}; retain one multiallelic site record".format(
                    record.contig, record.position + 1
                )
            )
        if previous_coordinate is not None and coordinate < previous_coordinate:
            raise Tier3ValidationError("VCF records must be coordinate sorted")
        previous_coordinate = coordinate
        observed_ref = fasta[record.contig][
            record.position : record.position + len(record.ref)
        ].upper()
        if observed_ref != record.ref.upper():
            raise Tier3ValidationError(
                "VCF REF mismatch at {}:{} ({} != {})".format(
                    record.contig, record.position + 1, record.ref, observed_ref
                )
            )
        for sample in selected_samples:
            if sample not in record.genotypes:
                raise Tier3ValidationError("selected sample {!r} is missing from a VCF record".format(sample))
            for allele in record.genotypes[sample]:
                if allele is not None and not 0 <= allele <= len(record.alts):
                    raise Tier3ValidationError("genotype allele index exceeds ALT count")
        if not _position_in_intervals(record.position, intervals.get(record.contig, ())):
            exclusion_counts["outside_denominator"] += 1
            continue
        key = (record.contig, record.position // BLOCK_SIZE_BP)
        base = fasta[record.contig][record.position].upper()
        selected_genotypes = [_coerce_genotype(record.genotypes[sample], design) for sample in selected_samples]
        called = sum(allele is not None for genotype in selected_genotypes for allele in genotype)
        called_histogram[called] += 1
        symbolic_reference_block = bool(record.alts) and all(alt.startswith("<") for alt in record.alts)
        snv = (
            len(record.ref) == 1
            and record.ref.upper() in DNA
            and bool(record.alts)
            and all(len(alt) == 1 and alt.upper() in DNA for alt in record.alts)
        )
        invalid_reason: Optional[str] = None
        if record.filters:
            invalid_reason = "filtered"
        elif base not in DNA:
            invalid_reason = "ambiguous_reference"
        elif called < minimum_called:
            invalid_reason = "insufficient_called_chromosomes"
        elif record.alts and not snv and not symbolic_reference_block:
            invalid_reason = "non_snv"
        if invalid_reason:
            exclusion_counts[invalid_reason] += 1
            # A rejected gVCF reference-confidence block removes its full END
            # span. Other normalized records represent one reference anchor.
            removal_end = (
                int(record.info["END"])
                if denominator_kind == "gvcf" and symbolic_reference_block and "END" in record.info
                else record.position + 1
            )
            position = record.position
            while position < removal_end:
                block = position // BLOCK_SIZE_BP
                boundary = min(removal_end, (block + 1) * BLOCK_SIZE_BP)
                segment = fasta[record.contig][position:boundary].upper()
                acgt = sum(segment.count(candidate) for candidate in "ACGT")
                block_stats[(record.contig, block)]["den"] -= float(acgt)
                position = boundary
            positions = fourfold_by_contig.get(record.contig, [])
            first = bisect.bisect_left(positions, record.position)
            last = bisect.bisect_left(positions, removal_end)
            for position in positions[first:last]:
                class_name = "S" if fasta[record.contig][position].upper() in "GC" else "W"
                block_stats[(record.contig, position // BLOCK_SIZE_BP)][class_name + "_den"] -= 1.0
            continue
        # Invariant explicit records have diversity zero. Symbolic gVCF ALT is
        # reference-confidence evidence, not a variant allele.
        if not snv:
            continue
        contribution = allele_pairwise_diversity(selected_genotypes, minimum_called)
        if contribution is None:  # defended above, retained as a fail-closed invariant
            raise Tier3ValidationError("call-count invariant failed")
        block_stats[key]["num"] += contribution
        informative_snvs += 1
        if (record.contig, record.position) in fourfold:
            class_name = "S" if base in "GC" else "W"
            block_stats[key][class_name + "_num"] += contribution

    diversity_sum = sum(stats["num"] for stats in block_stats.values())
    callable_sites = int(round(sum(stats["den"] for stats in block_stats.values())))
    if callable_sites <= 0:
        raise Tier3ValidationError("no callable sites remain after QC")
    population_pi = diversity_sum / callable_sites

    result: Dict[str, Any] = {
        "policy_id": POLICY_ID,
        "dataset_id": dataset_id,
        "source_modality": "deposited_exact_reference_variants_plus_explicit_denominator",
        "input_provenance": {
            "vcf_sha256": sha256_file(vcf_path),
            "fasta_sha256": sha256_file(fasta_path),
            "callable_mask_sha256": sha256_file(callable_bed_path) if callable_bed_path else None,
            "selected_sample_list_sha256": hashlib.sha256(
                "".join(sample + "\n" for sample in selected_samples).encode("utf-8")
            ).hexdigest(),
        },
        "sample_design": {
            "design": design,
            "selected_samples": list(selected_samples),
            "sampling_units": len(selected_samples),
            "nominal_chromosomes": nominal_chromosomes,
            "minimum_called_chromosomes": minimum_called,
        },
        "denominator": {
            "kind": denominator_kind,
            "invariant_sites_explicit": True,
            "initial_callable_acgt_sites": initial_callable,
        },
        "population_pi": {
            "diversity_sum": diversity_sum,
            "numerator": diversity_sum,
            "callable_sites": callable_sites,
            "denominator": callable_sites,
            "callable_count": callable_sites,
            "point_estimate": population_pi,
            "nominal_chromosomes": nominal_chromosomes,
            "called_chromosome_histogram_at_records": {
                str(key): value for key, value in sorted(called_histogram.items())
            },
            "informative_snv_sites": informative_snvs,
            "bootstrap": _bootstrap(
                dataset_id, "population_pi", block_stats, bootstrap_replicates, False, True
            ),
        },
        "exclusion_counts": dict(sorted(exclusion_counts.items())),
        "annotation": annotation_audit,
        "pi_S_over_pi_W": None,
        "polarization_gate": {
            "requested": polarization_gate is not None,
            "passed": False,
            "reason": "deferred_by_tier3-decisions-v1",
        },
    }

    if fourfold:
        class_results: Dict[str, Dict[str, Any]] = {}
        for class_name in ("S", "W"):
            numerator = sum(stats[class_name + "_num"] for stats in block_stats.values())
            denominator = int(round(sum(stats[class_name + "_den"] for stats in block_stats.values())))
            class_results[class_name] = {
                "diversity_sum": numerator,
                "numerator": numerator,
                "callable_sites": denominator,
                "denominator": denominator,
                "callable_count": denominator,
                "point_estimate": numerator / denominator if denominator else None,
            }
        pi_s, pi_w = class_results["S"]["point_estimate"], class_results["W"]["point_estimate"]
        unavailable_reason: Optional[str] = None
        if min(class_results["S"]["callable_sites"], class_results["W"]["callable_sites"]) < minimum_4d_class_sites:
            unavailable_reason = "fewer_than_{}_callable_sites_in_a_class".format(minimum_4d_class_sites)
        elif pi_w is None or pi_w == 0:
            unavailable_reason = "pi_W_is_zero_or_unavailable"
        ratio_value = pi_s / pi_w if unavailable_reason is None and pi_s is not None and pi_w else None
        result["pi_S_over_pi_W"] = {
            "reference_conditioned": True,
            "class_definition": "forward_reference_GC_is_S_AT_is_W",
            "S": class_results["S"],
            "W": class_results["W"],
            "point_estimate": ratio_value,
            "unavailable_reason": unavailable_reason,
            "minimum_class_callable_sites": minimum_4d_class_sites,
            "bootstrap": _bootstrap(
                dataset_id,
                "pi_S_over_pi_W",
                block_stats,
                bootstrap_replicates,
                True,
                ratio_value is not None,
            ),
        }
    return result


def _read_samples(path: Path) -> List[str]:
    samples = [line.strip() for line in Path(path).read_text(encoding="utf-8").splitlines() if line.strip()]
    if len(samples) != len(set(samples)):
        raise Tier3ValidationError("sample list contains duplicates")
    return samples


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset-id", required=True)
    parser.add_argument("--vcf", required=True, type=Path)
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--selected-samples", required=True, type=Path)
    parser.add_argument("--design", required=True, choices=("wild_diploid", "inbred_lines_haploidized", "haploid"))
    parser.add_argument("--denominator-kind", required=True, choices=sorted(DENOMINATOR_KINDS))
    parser.add_argument("--callable-bed", type=Path)
    parser.add_argument("--gff", type=Path)
    parser.add_argument("--annotation-provenance", type=Path, help="JSON with provider/release/assembly/status/code")
    parser.add_argument("--contig-mapping", type=Path)
    parser.add_argument("--provider-cds", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args(argv)
    annotation_metadata = None
    if args.annotation_provenance:
        annotation_metadata = json.loads(args.annotation_provenance.read_text(encoding="utf-8"))
    result = compute_population_pi(
        dataset_id=args.dataset_id,
        vcf_path=args.vcf,
        fasta_path=args.fasta,
        selected_samples=_read_samples(args.selected_samples),
        design=args.design,
        denominator_kind=args.denominator_kind,
        callable_bed_path=args.callable_bed,
        gff_path=args.gff,
        annotation_metadata=annotation_metadata,
        contig_mapping_path=args.contig_mapping,
        provider_cds_path=args.provider_cds,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Tier3ValidationError as error:
        raise SystemExit("tier3b population VCF rejected: {}".format(error))
