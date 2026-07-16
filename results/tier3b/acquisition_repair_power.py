#!/usr/bin/env python3
"""Fail-closed Tier 3B 4D denominator and block-power preflight.

This validator deliberately computes only frozen eligibility quantities.  It
does not create population result rows, point estimates, or intervals owned by
the downstream run.  Coordinates are zero-based internally, matching the
production population estimator.
"""

from __future__ import annotations

import argparse
import bisect
import collections
import hashlib
import json
import os
import sys
from pathlib import Path
from typing import Any, Dict, List, Mapping, Sequence, Set, Tuple

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from analysis.tier3_common import (  # noqa: E402
    DNA,
    GFFAnnotation,
    Tier3ValidationError,
    collect_fourfold_sites,
    fasta_dictionary,
    parse_gff,
    read_fasta,
    resolve_contig_aliases,
    sha256_file,
)
from analysis.tier3b_popvcf_compute import (  # noqa: E402
    _coerce_genotype,
    _position_in_intervals,
    iter_vcf_records,
    read_callable_bed,
    read_vcf_header,
)


POLICY_ID = "tier3-decisions-v1"
BLOCK_SIZE_BP = 1_000_000
MINIMUM_BLOCKS = 20
MINIMUM_CLASS_SITES = 10_000
DECLARED_EXCLUSIONS = {
    "AGAP000192-RA": (
        "native AgamP4.12 CDS is empty or ambiguous; provider GFF retained byte-identical"
    )
}


def sha256_samples(samples: Sequence[str]) -> str:
    return hashlib.sha256("".join(sample + "\n" for sample in samples).encode()).hexdigest()


def exact_fourfold_sites(fasta_path: Path, gff_path: Path) -> Tuple[Set[Tuple[str, int]], Dict[str, str]]:
    """Apply the run-side audited invalid-canonical-transcript exclusions."""

    fasta = read_fasta(fasta_path)
    annotation = parse_gff(gff_path)
    resolved = resolve_contig_aliases(
        fasta_dictionary(fasta), annotation.sequence_regions, {}
    )
    annotation_fasta = dict(fasta)
    for annotation_contig, fasta_contig in resolved.items():
        annotation_fasta[annotation_contig] = fasta[fasta_contig]

    unknown = sorted(set(DECLARED_EXCLUSIONS) - set(annotation.transcripts))
    if unknown:
        raise Tier3ValidationError("declared annotation exclusions are absent: {!r}".format(unknown))
    retained = {
        gene: transcript
        for gene, transcript in annotation.canonical_transcripts.items()
        if transcript not in DECLARED_EXCLUSIONS
    }
    automatic: Dict[str, str] = {}
    for gene, transcript in sorted(retained.items()):
        one = GFFAnnotation(annotation.sequence_regions, annotation.transcripts, {gene: transcript})
        try:
            collect_fourfold_sites(annotation_fasta, one, 1)
        except Tier3ValidationError as error:
            automatic[transcript] = str(error)
    retained = {
        gene: transcript for gene, transcript in retained.items() if transcript not in automatic
    }
    if not retained:
        raise Tier3ValidationError("annotation audit excluded every canonical transcript")
    audited = GFFAnnotation(annotation.sequence_regions, annotation.transcripts, retained)
    sites, _overlap_exclusions = collect_fourfold_sites(annotation_fasta, audited, 1)
    return {(resolved.get(contig, contig), position) for contig, position in sites}, {
        **DECLARED_EXCLUSIONS,
        **automatic,
    }


def validate_tuple(
    tuple_id: str,
    vcf_path: Path,
    callable_path: Path,
    sample_path: Path,
    fasta_path: Path,
    fourfold: Set[Tuple[str, int]],
    enforce_power: bool = True,
) -> Dict[str, Any]:
    fasta = read_fasta(fasta_path)
    intervals = read_callable_bed(callable_path, fasta_dictionary(fasta))
    samples = [line for line in sample_path.read_text().splitlines() if line]
    if len(samples) != len(set(samples)) or len(samples) != 20:
        raise Tier3ValidationError("{} does not retain 20 unique samples".format(tuple_id))

    vcf_contigs, vcf_samples = read_vcf_header(vcf_path)
    if tuple(vcf_samples) != tuple(samples):
        raise Tier3ValidationError("{} VCF sample order differs from locked cohort".format(tuple_id))
    if vcf_contigs != fasta_dictionary(fasta):
        raise Tier3ValidationError("{} VCF/FASTA dictionaries differ".format(tuple_id))

    counts: Dict[Tuple[str, int], Dict[str, int]] = collections.defaultdict(
        lambda: {"S": 0, "W": 0}
    )
    fourfold_by_contig: Dict[str, List[int]] = collections.defaultdict(list)
    for contig, position in sorted(fourfold):
        fourfold_by_contig[contig].append(position)
        if _position_in_intervals(position, intervals.get(contig, ())):
            cls = "S" if fasta[contig][position].upper() in "GC" else "W"
            counts[(contig, position // BLOCK_SIZE_BP)][cls] += 1

    initial = {cls: sum(value[cls] for value in counts.values()) for cls in ("S", "W")}
    exclusions: Dict[str, int] = collections.Counter()
    previous: Tuple[int, int] | None = None
    contig_order = {contig: index for index, contig in enumerate(vcf_contigs)}
    records = 0
    multiallelic = 0
    missing_genotype_records = 0
    for record in iter_vcf_records(vcf_path, vcf_samples):
        records += 1
        coordinate = (contig_order[record.contig], record.position)
        if previous is not None and coordinate <= previous:
            reason = "duplicate coordinate" if coordinate == previous else "unsorted coordinate"
            raise Tier3ValidationError("{} has {} at {}:{}".format(tuple_id, reason, record.contig, record.position + 1))
        previous = coordinate
        observed_ref = fasta[record.contig][record.position : record.position + len(record.ref)].upper()
        if observed_ref != record.ref.upper():
            raise Tier3ValidationError("{} REF mismatch at {}:{}".format(tuple_id, record.contig, record.position + 1))
        if len(record.alts) > 1:
            multiallelic += 1
        selected_genotypes = []
        any_missing = False
        for sample in samples:
            genotype = record.genotypes.get(sample)
            if genotype is None:
                raise Tier3ValidationError("{} record lacks selected sample {}".format(tuple_id, sample))
            if any(allele is not None and not 0 <= allele <= len(record.alts) for allele in genotype):
                raise Tier3ValidationError("{} genotype allele exceeds ALT count".format(tuple_id))
            any_missing = any_missing or any(allele is None for allele in genotype)
            selected_genotypes.append(_coerce_genotype(genotype, "wild_diploid"))
        missing_genotype_records += int(any_missing)
        if not _position_in_intervals(record.position, intervals.get(record.contig, ())):
            continue
        called = sum(allele is not None for genotype in selected_genotypes for allele in genotype)
        symbolic_reference = bool(record.alts) and all(alt.startswith("<") for alt in record.alts)
        snv = (
            len(record.ref) == 1
            and record.ref.upper() in DNA
            and bool(record.alts)
            and all(len(alt) == 1 and alt.upper() in DNA for alt in record.alts)
        )
        invalid_reason = None
        base = fasta[record.contig][record.position].upper()
        if record.filters:
            invalid_reason = "filtered"
        elif base not in DNA:
            invalid_reason = "ambiguous_reference"
        elif called < 36:
            invalid_reason = "insufficient_called_chromosomes"
        elif record.alts and not snv and not symbolic_reference:
            invalid_reason = "non_snv"
        if invalid_reason is None:
            continue
        exclusions[invalid_reason] += 1
        positions = fourfold_by_contig.get(record.contig, ())
        first = bisect.bisect_left(positions, record.position)
        last = bisect.bisect_left(positions, record.position + 1)
        for position in positions[first:last]:
            if _position_in_intervals(position, intervals.get(record.contig, ())):
                cls = "S" if fasta[record.contig][position].upper() in "GC" else "W"
                counts[(record.contig, position // BLOCK_SIZE_BP)][cls] -= 1

    final = {cls: sum(value[cls] for value in counts.values()) for cls in ("S", "W")}
    eligible_blocks = sum(1 for value in counts.values() if value["S"] + value["W"] > 0)
    if enforce_power and min(final.values()) < MINIMUM_CLASS_SITES:
        raise Tier3ValidationError("{} has insufficient exact 4D counts: {}".format(tuple_id, final))
    if enforce_power and eligible_blocks < MINIMUM_BLOCKS:
        raise Tier3ValidationError("{} has only {} eligible blocks".format(tuple_id, eligible_blocks))
    return {
        "tuple_id": tuple_id,
        "policy_id": POLICY_ID,
        "sample_count": len(samples),
        "selected_sample_list_sha256": sha256_samples(samples),
        "vcf_records_streamed": records,
        "multiallelic_records": multiallelic,
        "records_with_any_missing_genotype": missing_genotype_records,
        "initial_callable_fourfold": initial,
        "exact_callable_fourfold": final,
        "eligible_ratio_blocks": eligible_blocks,
        "block_size_bp": BLOCK_SIZE_BP,
        "minimum_blocks": MINIMUM_BLOCKS,
        "minimum_class_sites": MINIMUM_CLASS_SITES,
        "record_level_denominator_exclusions": dict(sorted(exclusions.items())),
        "vcf_sha256": sha256_file(vcf_path),
        "vcf_index_sha256": sha256_file(Path(str(vcf_path) + ".tbi")),
        "callable_bed_sha256": sha256_file(callable_path),
        "status": (
            "PASS"
            if min(final.values()) >= MINIMUM_CLASS_SITES and eligible_blocks >= MINIMUM_BLOCKS
            else "EXPECTED_INSUFFICIENT_CURRENT_TUPLE"
        ),
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate", required=True, type=Path)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--mode", choices=("current", "repair"), default="repair")
    args = parser.parse_args(argv)
    fasta = args.root / "ag1000g_phase3/staged/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa"
    gff = args.root / "ag1000g_phase3/source/Anopheles-gambiae-PEST_BASEFEATURES_AgamP4.12.gff3.gz"
    fourfold, annotation_exclusions = exact_fourfold_sites(fasta, gff)
    if args.mode == "repair":
        tag = "3R_10000000_30999999"
        specifications = (
            (
                "ag1000g_phase3_ao_coluzzii",
                args.candidate / "ao" / ("ao_luanda_coluzzii." + tag + ".all_sites.vcf.gz"),
                args.candidate / "ao" / ("ao_luanda_coluzzii." + tag + ".callable.bed"),
                args.root / "ag1000g_phase3/staged/selected.samples.txt",
            ),
            (
                "ag1000g_phase3_gm_coluzzii",
                args.candidate / "gm" / ("gm_walikunda_coluzzii." + tag + ".all_sites.vcf.gz"),
                args.candidate / "gm" / ("gm_walikunda_coluzzii." + tag + ".callable.bed"),
                args.root / "ag1000g_phase3_gm_coluzzii/staged/selected.samples.txt",
            ),
        )
    else:
        tag = "3R_10000000_10999999"
        specifications = (
            (
                "ag1000g_phase3_ao_coluzzii",
                args.root / "ag1000g_phase3/staged" / ("ao_luanda_coluzzii." + tag + ".all_sites.vcf.gz"),
                args.root / "ag1000g_phase3/staged" / ("ao_luanda_coluzzii." + tag + ".callable.bed"),
                args.root / "ag1000g_phase3/staged/selected.samples.txt",
            ),
            (
                "ag1000g_phase3_gm_coluzzii",
                args.root / "ag1000g_phase3_gm_coluzzii/staged" / ("gm_walikunda_coluzzii." + tag + ".all_sites.vcf.gz"),
                args.root / "ag1000g_phase3_gm_coluzzii/staged" / ("gm_walikunda_coluzzii." + tag + ".callable.bed"),
                args.root / "ag1000g_phase3_gm_coluzzii/staged/selected.samples.txt",
            ),
        )
    tuples = []
    for specification in specifications:
        row = validate_tuple(
            *specification, fasta, fourfold, enforce_power=args.mode == "repair"
        )
        tuples.append(row)
        print(
            json.dumps(
                {
                    "progress": "tuple_validated",
                    "tuple_id": row["tuple_id"],
                    "exact_callable_fourfold": row["exact_callable_fourfold"],
                    "eligible_ratio_blocks": row["eligible_ratio_blocks"],
                },
                sort_keys=True,
            ),
            flush=True,
        )
    if args.mode == "current":
        if not all(
            row["exact_callable_fourfold"]["W"] < MINIMUM_CLASS_SITES
            and row["eligible_ratio_blocks"] < MINIMUM_BLOCKS
            for row in tuples
        ):
            raise Tier3ValidationError("current tuples did not reproduce both frozen power failures")
    result = {
        "schema_version": "tier3b-acquisition-repair-power-v1",
        "policy_id": POLICY_ID,
        "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
        "fasta_sha256": sha256_file(fasta),
        "gff_sha256": sha256_file(gff),
        "annotation_invalid_transcript_exclusions": annotation_exclusions,
        "fourfold_sites_after_annotation_audit": len(fourfold),
        "mode": args.mode,
        "tuples": tuples,
        "status": "PASS" if args.mode == "repair" else "EXPECTED_FAIL_REPRODUCED",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_name(args.output.name + ".partial-{}".format(os.getpid()))
    temporary.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    temporary.replace(args.output)
    print(json.dumps({"status": result["status"], "tuples": result["tuples"]}, sort_keys=True))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Tier3ValidationError as error:
        raise SystemExit("repair power preflight rejected: {}".format(error))
