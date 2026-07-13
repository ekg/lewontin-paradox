#!/usr/bin/env python3
"""Compute Tier 3a individual heterozygosity and exact-reference composition.

The estimator in this module is always ``individual_snv_heterozygosity``:
heterozygous SNV sites divided by explicit callable/alignable A/C/G/T bases
for one diploid individual.  It is never labelled population pi.  GC3 and
reference-conditioned fourfold statistics are emitted only when the GFF is
native to, and audited against, the exact reference FASTA.
"""

from __future__ import annotations

import argparse
import bisect
import collections
import json
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple, Union

try:
    from analysis.tier3_common import (
        DNA,
        Tier3ValidationError,
        collect_fourfold_sites,
        fasta_dictionary,
        merge_intervals,
        parse_gff,
        parse_gt,
        read_fasta,
        reconstruct_cds_with_positions,
        resolve_contig_aliases,
        sha256_file,
    )
except ModuleNotFoundError:  # direct ``python analysis/<script>.py`` execution
    from tier3_common import (  # type: ignore[no-redef]
        DNA,
        Tier3ValidationError,
        collect_fourfold_sites,
        fasta_dictionary,
        merge_intervals,
        parse_gff,
        parse_gt,
        read_fasta,
        reconstruct_cds_with_positions,
        resolve_contig_aliases,
        sha256_file,
    )


STOP_CODONS_CODE_1 = frozenset({"TAA", "TAG", "TGA"})


def _read_callable_bed(
    path: Union[str, Path], reference: Mapping[str, str]
) -> Dict[str, List[Tuple[int, int]]]:
    intervals: Dict[str, List[Tuple[int, int]]] = collections.defaultdict(list)
    with Path(path).open("r", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 3 or fields[0] not in reference:
                raise Tier3ValidationError(f"invalid callable BED contig at line {line_number}")
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError as error:
                raise Tier3ValidationError(f"non-integer callable BED coordinate at line {line_number}") from error
            if not 0 <= start < end <= len(reference[fields[0]]):
                raise Tier3ValidationError(f"callable BED out of FASTA bounds at line {line_number}")
            intervals[fields[0]].append((start, end))
    if not intervals:
        raise Tier3ValidationError("callable BED has no intervals; invariant denominator is missing")
    return {contig: merge_intervals(values) for contig, values in intervals.items()}


class _IntervalLookup:
    def __init__(self, intervals: Mapping[str, Sequence[Tuple[int, int]]]):
        self.intervals = {contig: list(values) for contig, values in intervals.items()}
        self.starts = {
            contig: [start for start, _end in values] for contig, values in self.intervals.items()
        }

    def contains(self, contig: str, position: int) -> bool:
        values = self.intervals.get(contig, [])
        index = bisect.bisect_right(self.starts.get(contig, []), position) - 1
        return index >= 0 and position < values[index][1]


def _whole_genome_gc(reference: Mapping[str, str]) -> Dict[str, Any]:
    gc = acgt = ambiguous = 0
    for sequence in reference.values():
        for base in sequence.upper():
            if base in DNA:
                acgt += 1
                gc += base in {"G", "C"}
            else:
                ambiguous += 1
    return {
        "gc_bases": gc,
        "acgt_bases": acgt,
        "ambiguous_bases_excluded": ambiguous,
        "value": gc / acgt if acgt else None,
    }


def _count_callable_reference(
    reference: Mapping[str, str], intervals: Mapping[str, Sequence[Tuple[int, int]]]
) -> Tuple[int, int]:
    callable_bases = ambiguous = 0
    for contig, values in intervals.items():
        sequence = reference[contig]
        for start, end in values:
            for base in sequence[start:end].upper():
                if base in DNA:
                    callable_bases += 1
                else:
                    ambiguous += 1
    return callable_bases, ambiguous


def _annotation_unavailable(reason: str, provenance: Optional[Mapping[str, Any]]) -> Dict[str, Any]:
    return {
        "status": reason,
        "primary_eligible": False,
        "provenance": dict(provenance) if provenance else None,
    }


def _annotation_metrics(
    reference: Mapping[str, str],
    reference_fasta: Path,
    reference_accession: Optional[str],
    annotation_gff: Optional[Path],
    provenance: Optional[Mapping[str, Any]],
) -> Tuple[Dict[str, Any], Optional[Dict[str, Any]], Optional[Set[Tuple[str, int]]]]:
    """Return annotation audit, GC3 metrics, and exact-reference 4D sites.

    Missing or non-native provenance degrades annotation-derived observables
    without suppressing whole-genome GC or denominator-backed heterozygosity.
    A content contradiction in a purported native tuple fails closed.
    """

    if annotation_gff is None or provenance is None:
        return _annotation_unavailable("unavailable_missing_native_annotation", provenance), None, None
    required = {
        "provider",
        "release",
        "assembly_accession_version",
        "status",
        "genetic_code",
        "fasta_sha256",
        "gff_sha256",
        "contig_mapping",
    }
    missing = sorted(required - set(provenance))
    if missing:
        return _annotation_unavailable(
            "unavailable_incomplete_annotation_provenance:" + ",".join(missing), provenance
        ), None, None
    if provenance["status"] != "native":
        return _annotation_unavailable("unavailable_non_native_annotation", provenance), None, None
    if not reference_accession:
        return _annotation_unavailable("unavailable_missing_reference_accession", provenance), None, None
    if provenance["assembly_accession_version"] != reference_accession:
        raise Tier3ValidationError("native annotation assembly does not match the exact reference accession")
    if provenance.get("audit_passed") is False or provenance.get("sampled_cds_mismatches", 0) != 0:
        raise Tier3ValidationError("annotation CDS reconstruction audit failed")
    if provenance["fasta_sha256"] != sha256_file(reference_fasta):
        raise Tier3ValidationError("annotation provenance FASTA checksum mismatch")
    if provenance["gff_sha256"] != sha256_file(annotation_gff):
        raise Tier3ValidationError("annotation provenance GFF checksum mismatch")
    genetic_code = int(provenance["genetic_code"])
    if genetic_code != 1:
        return _annotation_unavailable("unavailable_unsupported_genetic_code", provenance), None, None

    annotation = parse_gff(annotation_gff)
    aliases = provenance["contig_mapping"]
    if not isinstance(aliases, Mapping):
        raise Tier3ValidationError("annotation contig_mapping provenance must be an object")
    resolved = resolve_contig_aliases(
        fasta_dictionary(reference), annotation.sequence_regions, aliases
    )
    annotation_fasta = dict(reference)
    for annotation_name, fasta_name in resolved.items():
        annotation_fasta[annotation_name] = reference[fasta_name]

    # This reconstructs every retained canonical CDS and enforces phase,
    # translation, ambiguity, internal-stop, and overlap/frame policy.
    fourfold_annotation_coordinates, fourfold_excluded = collect_fourfold_sites(
        annotation_fasta, annotation, genetic_code
    )
    fourfold = {
        (resolved.get(contig, contig), position)
        for contig, position in fourfold_annotation_coordinates
    }

    frame_observations: Dict[Tuple[str, int], Set[int]] = collections.defaultdict(set)
    terminal_stops = 0
    retained_transcripts = list(annotation.canonical_transcripts.values())
    for transcript_id in retained_transcripts:
        sequence, positions = reconstruct_cds_with_positions(
            annotation_fasta, annotation.transcripts[transcript_id]
        )
        effective_length = len(sequence)
        if sequence[-3:] in STOP_CODONS_CODE_1:
            terminal_stops += 1
            effective_length -= 3
        for offset in range(effective_length):
            annotation_contig, position = positions[offset]
            frame_observations[(resolved.get(annotation_contig, annotation_contig), position)].add(
                offset % 3
            )
    discordant = {key for key, frames in frame_observations.items() if len(frames) > 1}
    third_positions = {
        key for key, frames in frame_observations.items() if frames == {2} and key not in discordant
    }
    gc = sum(reference[contig][position].upper() in {"G", "C"} for contig, position in third_positions)
    audit = {
        "status": "eligible_native_exact_reference",
        "primary_eligible": True,
        "provider": provenance["provider"],
        "release": provenance["release"],
        "assembly_accession_version": provenance["assembly_accession_version"],
        "fasta_sha256": sha256_file(reference_fasta),
        "gff_sha256": sha256_file(annotation_gff),
        "sequence_regions": dict(sorted(annotation.sequence_regions.items())),
        "resolved_contig_mapping": dict(sorted(resolved.items())),
        "declared_contig_mapping": dict(sorted(aliases.items())),
        "genetic_code": genetic_code,
        "native_vs_projected_status": "native",
        "canonical_transcripts_reconstructed": len(retained_transcripts),
        "sampled_cds_reconstruction": {
            "rule": "all_retained_canonical_transcripts",
            "count": len(retained_transcripts),
            "fasta_reconstruction_passed": True,
            "provider_sequence_mismatches": int(provenance.get("sampled_cds_mismatches", 0)),
            "provider_sequence_comparison": (
                "passed" if "sampled_cds_mismatches" in provenance else "not_deposited"
            ),
        },
    }
    gc3 = {
        "status": "eligible",
        "gc_third_positions": gc,
        "third_positions": len(third_positions),
        "value": gc / len(third_positions) if third_positions else None,
        "genes": len(annotation.canonical_transcripts),
        "transcripts": len(retained_transcripts),
        "terminal_stop_codons_excluded": terminal_stops,
        "discordant_overlap_positions_excluded": len(discordant),
        "fourfold_overlap_positions_excluded": len(fourfold_excluded),
    }
    if not third_positions:
        gc3["status"] = "unavailable_no_valid_third_positions"
    return audit, gc3, fourfold


def _audit_and_collect_variants(
    normalized_bcf: Path,
    reference: Mapping[str, str],
    callable_lookup: _IntervalLookup,
    sample: Optional[str],
) -> Tuple[str, Set[Tuple[str, int]], int, int]:
    try:
        import pysam
    except ImportError as error:  # pragma: no cover - supplied by Guix
        raise Tier3ValidationError("pysam is required; enter the pinned pure Guix environment") from error

    heterozygous: Set[Tuple[str, int]] = set()
    indels = non_heterozygous_snvs = 0
    with pysam.VariantFile(str(normalized_bcf)) as variants:
        header_contigs = {
            name: variants.header.contigs[name].length for name in variants.header.contigs
        }
        if set(header_contigs) != set(reference):
            raise Tier3ValidationError("VCF/BCF contig dictionary does not exactly match reference FASTA")
        for contig, sequence in reference.items():
            if header_contigs[contig] != len(sequence):
                raise Tier3ValidationError(f"VCF/BCF contig length mismatch for {contig!r}")
        samples = list(variants.header.samples)
        if sample is None:
            if len(samples) != 1:
                raise Tier3ValidationError("single-individual BCF requires an explicit sample selection")
            selected = samples[0]
        else:
            if sample not in samples:
                raise Tier3ValidationError(f"sample {sample!r} is absent from BCF")
            selected = sample

        seen_exact: Set[Tuple[str, int, str, Tuple[str, ...]]] = set()
        for record in variants:
            if record.contig not in reference:
                raise Tier3ValidationError(f"VCF contig {record.contig!r} absent from reference FASTA")
            observed_ref = reference[record.contig][record.start : record.start + len(record.ref)].upper()
            if observed_ref != record.ref.upper():
                raise Tier3ValidationError(
                    f"VCF REF mismatch at {record.contig}:{record.pos}: {record.ref} != {observed_ref}"
                )
            alts = tuple(record.alts or ())
            exact_key = (record.contig, record.start, record.ref, alts)
            if exact_key in seen_exact:
                raise Tier3ValidationError(f"duplicate normalized allele at {record.contig}:{record.pos}")
            seen_exact.add(exact_key)
            if not callable_lookup.contains(record.contig, record.start):
                continue
            raw_gt = record.samples[selected].get("GT")
            if raw_gt is None:
                raise Tier3ValidationError(
                    f"callable variant lacks GT at {record.contig}:{record.pos}"
                )
            gt = parse_gt(raw_gt, expected_ploidy=2)
            if any(allele is None for allele in gt):
                raise Tier3ValidationError(
                    f"callable mask contradicts missing genotype at {record.contig}:{record.pos}"
                )
            if any(int(allele) > len(alts) for allele in gt if allele is not None):
                raise Tier3ValidationError(f"genotype allele index exceeds ALT count at {record.contig}:{record.pos}")
            is_snv = (
                len(record.ref) == 1
                and bool(alts)
                and all(len(alt) == 1 and alt.upper() in DNA for alt in alts)
                and record.ref.upper() in DNA
            )
            if not is_snv:
                indels += 1
                continue
            if gt[0] != gt[1]:
                heterozygous.add((record.contig, record.start))
            else:
                non_heterozygous_snvs += 1
    return selected, heterozygous, indels, non_heterozygous_snvs


def compute_tier3a(
    *,
    dataset_id: str,
    reference_fasta: Union[str, Path],
    normalized_bcf: Optional[Union[str, Path]],
    callable_bed: Optional[Union[str, Path]],
    sample: Optional[str],
    modality: str,
    reference_accession: Optional[str],
    annotation_gff: Optional[Union[str, Path]] = None,
    annotation_provenance: Optional[Mapping[str, Any]] = None,
    minimum_4d_class_sites: int = 10_000,
) -> Dict[str, Any]:
    """Compute one Tier 3a individual row from a same-reference tuple."""

    allowed_modalities = {
        "deposited_exact_reference_variants_plus_mask",
        "direct_wfmash_extended_cigar",
    }
    if modality not in allowed_modalities:
        raise Tier3ValidationError(f"unsupported Tier 3a modality {modality!r}")
    if minimum_4d_class_sites < 1:
        raise Tier3ValidationError("minimum_4d_class_sites must be positive")
    reference_path = Path(reference_fasta).resolve()
    reference = read_fasta(reference_path)
    whole_gc = _whole_genome_gc(reference)
    annotation_path = Path(annotation_gff).resolve() if annotation_gff is not None else None
    annotation, gc3, fourfold_sites = _annotation_metrics(
        reference,
        reference_path,
        reference_accession,
        annotation_path,
        annotation_provenance,
    )

    base: Dict[str, Any] = {
        "dataset_id": dataset_id,
        "modality": modality,
        "statistic_label": "individual_snv_heterozygosity",
        "interpretation": "one diploid individual; conditional on callable/alignable reference bases",
        "population_pi": None,
        "reference": {
            "assembly_accession_version": reference_accession,
            "fasta_sha256": sha256_file(reference_path),
            "contig_dictionary": fasta_dictionary(reference),
        },
        "whole_genome_gc": whole_gc,
        "annotation": annotation,
        "gc3": gc3,
        "fourfold": None,
        "pi_S": None,
        "pi_W": None,
        "pi_S_over_pi_W": None,
    }
    if normalized_bcf is None:
        base.update(
            status="unavailable_missing_variant_source",
            individual_snv_heterozygosity=None,
            heterozygous_snvs=None,
            callable_bases=None,
        )
        return base
    if callable_bed is None:
        base.update(
            status="unavailable_missing_callable_denominator",
            individual_snv_heterozygosity=None,
            heterozygous_snvs=None,
            callable_bases=None,
            denominator={
                "status": "missing",
                "variant_only_reference_length_assumption_used": False,
            },
        )
        return base

    callable_path = Path(callable_bed).resolve()
    intervals = _read_callable_bed(callable_path, reference)
    callable_lookup = _IntervalLookup(intervals)
    callable_bases, ambiguous_masked = _count_callable_reference(reference, intervals)
    if callable_bases == 0:
        raise Tier3ValidationError("callable BED contains no unambiguous A/C/G/T reference bases")
    selected, heterozygous, indels, homozygous_snvs = _audit_and_collect_variants(
        Path(normalized_bcf).resolve(), reference, callable_lookup, sample
    )
    value = len(heterozygous) / callable_bases
    base.update(
        status="eligible",
        sample=selected,
        individual_snv_heterozygosity=value,
        heterozygous_snvs=len(heterozygous),
        callable_bases=callable_bases,
        denominator={
            "kind": (
                "h1_reference_alignable_mask"
                if modality == "direct_wfmash_extended_cigar"
                else "explicit_individual_callable_mask"
            ),
            "bed_sha256": sha256_file(callable_path),
            "callable_acgt_bases": callable_bases,
            "ambiguous_masked_bases_excluded": ambiguous_masked,
            "invariant_sites_explicit": True,
            "variant_only_reference_length_assumption_used": False,
        },
        variant_qc={
            "normalized_bcf_sha256": sha256_file(normalized_bcf),
            "heterozygous_snvs": len(heterozygous),
            "non_heterozygous_snvs": homozygous_snvs,
            "indel_or_complex_records": indels,
            "ref_allele_audit_passed": True,
            "contig_dictionary_passed": True,
        },
    )

    if fourfold_sites is not None:
        classes: Dict[str, Dict[str, Any]] = {}
        for class_name, bases in (("S", {"G", "C"}), ("W", {"A", "T"})):
            callable_class = {
                key
                for key in fourfold_sites
                if reference[key[0]][key[1]].upper() in bases
                and callable_lookup.contains(key[0], key[1])
            }
            heterozygous_class = heterozygous & callable_class
            class_value = len(heterozygous_class) / len(callable_class) if callable_class else None
            classes[class_name] = {
                "reference_base_class": class_name,
                "heterozygous_snvs": len(heterozygous_class),
                "callable_bases": len(callable_class),
                "value": class_value,
                "minimum_sites_passed": len(callable_class) >= minimum_4d_class_sites,
            }
        base["fourfold"] = classes
        if all(classes[name]["minimum_sites_passed"] for name in ("S", "W")):
            base["pi_S"] = classes["S"]["value"]
            base["pi_W"] = classes["W"]["value"]
            if base["pi_W"] not in {None, 0.0}:
                base["pi_S_over_pi_W"] = base["pi_S"] / base["pi_W"]
        else:
            base["fourfold_status"] = "unavailable_minimum_class_sites"
    return base


def common_callable_concordance(
    *,
    reference_fasta: Union[str, Path],
    left_callable_bed: Union[str, Path],
    right_callable_bed: Union[str, Path],
    left_heterozygous: Iterable[Tuple[str, int]],
    right_heterozygous: Iterable[Tuple[str, int]],
    synthetic_fixture: bool = False,
) -> Dict[str, Any]:
    """Compare modalities only on their exact common callable intersection."""

    reference = read_fasta(reference_fasta)
    left = _read_callable_bed(left_callable_bed, reference)
    right = _read_callable_bed(right_callable_bed, reference)
    right_lookup = _IntervalLookup(right)
    left_bases, _ = _count_callable_reference(reference, left)
    right_bases, _ = _count_callable_reference(reference, right)
    common: Set[Tuple[str, int]] = set()
    # This audit helper favors clarity; production masks are streamed by the
    # main estimator and callers normally invoke concordance per frozen core.
    for contig, intervals in left.items():
        sequence = reference[contig]
        for start, end in intervals:
            for position in range(start, end):
                if sequence[position].upper() in DNA and right_lookup.contains(contig, position):
                    common.add((contig, position))
    left_het = set(left_heterozygous) & common
    right_het = set(right_heterozygous) & common
    shared = left_het & right_het
    precision = len(shared) / len(right_het) if right_het else (1.0 if not left_het else 0.0)
    recall = len(shared) / len(left_het) if left_het else (1.0 if not right_het else 0.0)
    left_value = len(left_het) / len(common) if common else None
    right_value = len(right_het) / len(common) if common else None
    absolute_difference = (
        abs(left_value - right_value) if left_value is not None and right_value is not None else None
    )
    minimum_common = 1 if synthetic_fixture else 10_000_000
    diversity_tolerance = (
        max(5e-5, 0.05 * ((left_value + right_value) / 2))
        if left_value is not None and right_value is not None
        else None
    )
    coverage_passed = (
        len(common) >= minimum_common
        and len(common) / left_bases >= 0.8
        and len(common) / right_bases >= 0.8
    )
    exact_threshold = 1.0 if synthetic_fixture else 0.99
    genotype_concordance = (
        (len(common) - len(left_het ^ right_het)) / len(common) if common else None
    )
    return {
        "left_callable_bases": left_bases,
        "right_callable_bases": right_bases,
        "common_callable_bases": len(common),
        "common_fraction_left": len(common) / left_bases,
        "common_fraction_right": len(common) / right_bases,
        "snv_precision": precision,
        "snv_recall": recall,
        "heterozygous_nonheterozygous_genotype_concordance": genotype_concordance,
        "left_individual_snv_heterozygosity": left_value,
        "right_individual_snv_heterozygosity": right_value,
        "absolute_difference": absolute_difference,
        "diversity_tolerance": diversity_tolerance,
        "passed": bool(
            coverage_passed
            and precision >= exact_threshold
            and recall >= exact_threshold
            and genotype_concordance is not None
            and genotype_concordance >= exact_threshold
            and absolute_difference is not None
            and absolute_difference <= diversity_tolerance
        ),
        "population_pi": None,
    }


def _cli() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset-id", required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--reference-accession", required=True)
    parser.add_argument("--normalized-bcf", type=Path)
    parser.add_argument("--callable-bed", type=Path)
    parser.add_argument("--sample")
    parser.add_argument(
        "--modality",
        choices=[
            "deposited_exact_reference_variants_plus_mask",
            "direct_wfmash_extended_cigar",
        ],
        required=True,
    )
    parser.add_argument("--annotation-gff", type=Path)
    parser.add_argument("--annotation-provenance-json", type=Path)
    parser.add_argument("--minimum-4d-class-sites", type=int, default=10_000)
    parser.add_argument("--output", type=Path)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _cli().parse_args(argv)
    provenance = (
        json.loads(args.annotation_provenance_json.read_text(encoding="utf-8"))
        if args.annotation_provenance_json
        else None
    )
    result = compute_tier3a(
        dataset_id=args.dataset_id,
        reference_fasta=args.reference,
        normalized_bcf=args.normalized_bcf,
        callable_bed=args.callable_bed,
        sample=args.sample,
        modality=args.modality,
        reference_accession=args.reference_accession,
        annotation_gff=args.annotation_gff,
        annotation_provenance=provenance,
        minimum_4d_class_sites=args.minimum_4d_class_sites,
    )
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
