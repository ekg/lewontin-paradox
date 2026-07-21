#!/usr/bin/env python3
"""Reproducible raw-read validation helpers for the VGP ten-pair pilot.

The module deliberately separates measurement from interpretation.  Slurm
workers stream depth, pileup and Jellyfish observations into these helpers;
the final report then compares assembly and read estimates on the *same*
reference and callable denominator.  Such comparisons are paired sensitivity
analyses, not independent biological replicates.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import IO, Any, Iterable, Mapping, Sequence


CANONICAL_VGP_ROOT = Path("/moosefs/erikg/vgp")


@dataclass(frozen=True)
class DepthMaskSpec:
    """Inclusive depth bounds for one callable-mask sensitivity."""

    name: str
    minimum: int
    maximum: int

    def __post_init__(self) -> None:
        if not self.name or any(char not in "abcdefghijklmnopqrstuvwxyz0123456789_-" for char in self.name):
            raise ValueError(f"unsafe depth-mask name: {self.name!r}")
        if self.minimum < 0 or self.maximum < self.minimum:
            raise ValueError(f"invalid depth bounds for {self.name}: {self.minimum}..{self.maximum}")


def validate_root_config(path: Path) -> dict[str, Any]:
    """Load the single repository root contract and fail closed on drift."""

    data = json.loads(path.read_text())
    observed = Path(str(data.get("root", "")))
    if observed != CANONICAL_VGP_ROOT:
        raise ValueError(
            f"canonical VGP root must be {CANONICAL_VGP_ROOT}, observed {observed or '<missing>'}"
        )
    if not isinstance(data.get("layout"), dict):
        raise ValueError("root configuration is missing its layout mapping")
    return data


def _atomic_json(path: Path, value: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_name(path.name + ".partial")
    partial.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    partial.replace(path)


def _flush_interval(handle: IO[str], interval: tuple[str, int, int] | None) -> None:
    if interval is not None:
        handle.write(f"{interval[0]}\t{interval[1]}\t{interval[2]}\n")


def stream_depth_masks(
    rows: Iterable[str], output_dir: Path, specs: Sequence[DepthMaskSpec]
) -> dict[str, Any]:
    """Convert ordered ``samtools depth`` rows into exact BED masks.

    Input positions are 1-based; BED output is 0-based, half-open.  Runs are
    split on chromosome changes and coordinate gaps, which matters when depth
    was restricted to the assembly-derived callable BED.
    """

    if not specs:
        raise ValueError("at least one depth-mask specification is required")
    if len({spec.name for spec in specs}) != len(specs):
        raise ValueError("depth-mask names must be unique")
    output_dir.mkdir(parents=True, exist_ok=True)
    handles = {spec.name: (output_dir / f"{spec.name}.bed").open("w") for spec in specs}
    active: dict[str, tuple[str, int, int] | None] = {spec.name: None for spec in specs}
    callable_bp = Counter({spec.name: 0 for spec in specs})
    histogram: Counter[int] = Counter()
    observed = 0
    previous_chrom: str | None = None
    previous_pos = -1

    try:
        for line_number, raw in enumerate(rows, 1):
            if not raw.strip():
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise ValueError(f"depth line {line_number} has fewer than three columns")
            chrom, pos_text, depth_text = fields[:3]
            pos, depth = int(pos_text), int(depth_text)
            if pos < 1 or depth < 0:
                raise ValueError(f"invalid depth row at line {line_number}: {raw.rstrip()}")
            if previous_chrom == chrom and pos <= previous_pos:
                raise ValueError(f"depth coordinates are not strictly increasing at line {line_number}")
            if previous_chrom is not None and chrom != previous_chrom:
                for spec in specs:
                    _flush_interval(handles[spec.name], active[spec.name])
                    active[spec.name] = None
            histogram[depth] += 1
            observed += 1
            start, end = pos - 1, pos
            for spec in specs:
                state = active[spec.name]
                accepted = spec.minimum <= depth <= spec.maximum
                contiguous = state is not None and state[0] == chrom and state[2] == start
                if accepted:
                    callable_bp[spec.name] += 1
                    if contiguous:
                        active[spec.name] = (chrom, state[1], end)
                    else:
                        _flush_interval(handles[spec.name], state)
                        active[spec.name] = (chrom, start, end)
                elif state is not None:
                    _flush_interval(handles[spec.name], state)
                    active[spec.name] = None
            previous_chrom, previous_pos = chrom, pos
        for spec in specs:
            _flush_interval(handles[spec.name], active[spec.name])
    finally:
        for handle in handles.values():
            handle.close()

    summary: dict[str, Any] = {
        "schema_version": "vgp-read-validation-depth-masks-v1",
        "canonical_vgp_root": str(CANONICAL_VGP_ROOT),
        "observed_positions": observed,
        "depth_histogram": {str(key): histogram[key] for key in sorted(histogram)},
        "masks": {},
    }
    positive_depths = [depth for depth, count in histogram.items() if depth > 0 and count > 0]
    modal_positive = max(positive_depths, key=lambda depth: histogram[depth]) if positive_depths else None
    if observed and modal_positive is not None:
        low_cutoff = modal_positive * 0.5
        high_1_5_cutoff = modal_positive * 1.5
        high_2_cutoff = modal_positive * 2.0
        summary["depth_structure"] = {
            "modal_positive_depth": modal_positive,
            "zero_depth_fraction": histogram[0] / observed,
            "below_half_mode_fraction": sum(
                count for depth, count in histogram.items() if depth < low_cutoff
            ) / observed,
            "above_1_5x_mode_fraction": sum(
                count for depth, count in histogram.items() if depth > high_1_5_cutoff
            ) / observed,
            "above_2x_mode_fraction": sum(
                count for depth, count in histogram.items() if depth > high_2_cutoff
            ) / observed,
            "interpretation": (
                "low-depth enrichment is a duplication/mappability/dropout diagnostic and "
                "high-depth enrichment is a collapse/repeat diagnostic; neither alone proves an assembly error"
            ),
        }
    else:
        summary["depth_structure"] = {"status": "not_estimable", "reason": "no positive depth"}
    for spec in specs:
        count = callable_bp[spec.name]
        summary["masks"][spec.name] = {
            "minimum_depth_inclusive": spec.minimum,
            "maximum_depth_inclusive": spec.maximum,
            "callable_bp": count,
            "fraction_of_assembly_callable_input": count / observed if observed else None,
            "bed": str(output_dir / f"{spec.name}.bed"),
        }
    _atomic_json(output_dir / "depth_mask_summary.json", summary)
    with (output_dir / "depth_histogram.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("depth", "bases"))
        for depth in sorted(histogram):
            writer.writerow((depth, histogram[depth]))
    return summary


def parse_pileup_bases(bases: str, reference: str, alternate: str) -> dict[str, int]:
    """Count exact REF/ALT SNP observations in a SAMtools pileup base string."""

    reference, alternate = reference.upper(), alternate.upper()
    counts = {"ref": 0, "alt": 0, "other": 0, "deletion": 0}
    index = 0
    while index < len(bases):
        char = bases[index]
        if char == "^":
            index += 2  # start-of-read marker plus mapping-quality character
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
            length = int(bases[digit_start:index])
            index += length
            continue
        if char in ".,":
            counts["ref"] += 1
        elif char in "*#":
            counts["deletion"] += 1
        elif char in "<>":
            pass
        elif char.isalpha():
            base = char.upper()
            if base == alternate:
                counts["alt"] += 1
            elif base == reference:
                counts["ref"] += 1
            else:
                counts["other"] += 1
        index += 1
    return counts


def _read_assembly_sites(handle: IO[str]) -> list[tuple[str, int, str, str]]:
    sites: list[tuple[str, int, str, str]] = []
    seen: set[tuple[str, int, str, str]] = set()
    for line_number, raw in enumerate(handle, 1):
        if not raw.strip() or raw.startswith("#") or raw.startswith("chrom\t"):
            continue
        fields = raw.rstrip("\n").split("\t")
        if len(fields) < 4:
            raise ValueError(f"assembly-site line {line_number} has fewer than four columns")
        site = (fields[0], int(fields[1]), fields[2].upper(), fields[3].upper())
        if len(site[2]) != 1 or len(site[3]) != 1 or "," in site[3]:
            raise ValueError(f"assembly evidence accepts biallelic SNPs only: {site}")
        if site in seen:
            raise ValueError(f"duplicate assembly site: {site}")
        seen.add(site)
        sites.append(site)
    return sites


def assembly_evidence(
    assembly_sites: IO[str],
    pileup_rows: IO[str],
    output: IO[str],
    *,
    minimum_depth: int,
    maximum_depth: int,
) -> dict[str, Any]:
    """Classify direct read evidence at assembly-derived heterozygous SNPs."""

    sites = _read_assembly_sites(assembly_sites)
    by_position: dict[tuple[str, int], list[tuple[str, int, str, str]]] = {}
    for site in sites:
        by_position.setdefault((site[0], site[1]), []).append(site)
    observations: dict[tuple[str, int], tuple[int, str]] = {}
    for line_number, raw in enumerate(pileup_rows, 1):
        if not raw.strip():
            continue
        fields = raw.rstrip("\n").split("\t")
        if len(fields) < 5:
            raise ValueError(f"pileup line {line_number} has fewer than five columns")
        key = (fields[0], int(fields[1]))
        if key in observations:
            raise ValueError(f"duplicate pileup coordinate: {key}")
        observations[key] = (int(fields[3]), fields[4])

    fieldnames = (
        "chrom",
        "position_1based",
        "ref",
        "alt",
        "reported_depth",
        "parsed_ref_reads",
        "parsed_alt_reads",
        "parsed_other_reads",
        "allele_balance_alt",
        "classification",
    )
    writer = csv.DictWriter(output, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    classifications: Counter[str] = Counter()
    balances: list[float] = []
    for chrom, position, ref, alt in sites:
        observation = observations.get((chrom, position))
        if observation is None:
            depth, parsed, balance, classification = 0, {"ref": 0, "alt": 0, "other": 0}, None, "not_observed"
        else:
            depth, pileup_bases = observation
            parsed = parse_pileup_bases(pileup_bases, ref, alt)
            informative = parsed["ref"] + parsed["alt"]
            balance = parsed["alt"] / informative if informative else None
            if depth < minimum_depth or depth > maximum_depth:
                classification = "outside_depth_mask"
            elif parsed["ref"] >= 3 and parsed["alt"] >= 3 and balance is not None and 0.20 <= balance <= 0.80:
                classification = "supported_heterozygous"
                balances.append(balance)
            elif parsed["alt"] <= 1 and parsed["ref"] >= max(3, math.ceil(0.90 * depth)):
                classification = "contradicted_homozygous_reference"
            else:
                classification = "ambiguous"
        classifications[classification] += 1
        writer.writerow(
            {
                "chrom": chrom,
                "position_1based": position,
                "ref": ref,
                "alt": alt,
                "reported_depth": depth,
                "parsed_ref_reads": parsed["ref"],
                "parsed_alt_reads": parsed["alt"],
                "parsed_other_reads": parsed["other"],
                "allele_balance_alt": "" if balance is None else f"{balance:.8f}",
                "classification": classification,
            }
        )
    total_sites = len(sites)
    contradicted = classifications["contradicted_homozygous_reference"]
    unresolved = (
        classifications["ambiguous"]
        + classifications["outside_depth_mask"]
        + classifications["not_observed"]
    )
    if total_sites:
        contradicted_wilson = wilson_interval(contradicted, total_sites)
    else:
        contradicted_wilson = (None, None)
    return {
        "schema_version": "vgp-read-validation-assembly-site-evidence-v1",
        "canonical_vgp_root": str(CANONICAL_VGP_ROOT),
        "assembly_sites": len(sites),
        **{key: classifications[key] for key in (
            "supported_heterozygous",
            "contradicted_homozygous_reference",
            "ambiguous",
            "outside_depth_mask",
            "not_observed",
        )},
        "supported_alt_balance_mean": sum(balances) / len(balances) if balances else None,
        "concrete_false_positive_lower_bound_sites": contradicted,
        "concrete_false_positive_lower_bound_fraction": contradicted / total_sites if total_sites else None,
        "concrete_false_positive_fraction_wilson_95": list(contradicted_wilson),
        "unresolved_sites": unresolved,
        "candidate_false_positive_upper_bound_sites": contradicted + unresolved,
        "candidate_false_positive_upper_bound_fraction": (
            (contradicted + unresolved) / total_sites if total_sites else None
        ),
        "bound_scope": (
            "lower=depth-qualified homozygous-reference contradictions; "
            "upper=all contradictions plus ambiguous/out-of-mask/unobserved sites; "
            "read mapping and caller covariance remains systematic"
        ),
    }


def wilson_interval(successes: int, trials: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if trials <= 0 or successes < 0 or successes > trials:
        raise ValueError("invalid binomial counts")
    proportion = successes / trials
    denominator = 1.0 + z * z / trials
    centre = (proportion + z * z / (2.0 * trials)) / denominator
    radius = z * math.sqrt(proportion * (1.0 - proportion) / trials + z * z / (4.0 * trials * trials)) / denominator
    return max(0.0, centre - radius), min(1.0, centre + radius)


def _kmer_fraction_to_base_error(fraction: float, k: int) -> float:
    if not 0.0 <= fraction <= 1.0 or k <= 0:
        raise ValueError("invalid k-mer error fraction or k")
    return 1.0 - (1.0 - fraction) ** (1.0 / k)


def _error_to_qv(error_rate: float) -> float | None:
    if error_rate == 0.0:
        return None
    return -10.0 * math.log10(error_rate)


def kmer_qv(*, total_kmers: int, error_kmers: int, k: int) -> dict[str, Any]:
    """Compute a Merqury-formula QV from exact assembly-k-mer containment."""

    low, high = wilson_interval(error_kmers, total_kmers)
    fraction = error_kmers / total_kmers
    rate = _kmer_fraction_to_base_error(fraction, k)
    # Larger k-mer error fraction means lower QV, hence reversed interval ends.
    return {
        "schema_version": "vgp-read-validation-kmer-qv-v1",
        "k": k,
        "assembly_kmer_occurrences": total_kmers,
        "assembly_kmer_occurrences_below_trusted_read_threshold": error_kmers,
        "error_kmer_fraction": fraction,
        "error_rate": rate,
        "qv": _error_to_qv(rate),
        "qv_lower_95": _error_to_qv(_kmer_fraction_to_base_error(high, k)),
        "qv_upper_95": _error_to_qv(_kmer_fraction_to_base_error(low, k)),
        "interval_scope": "binomial_sampling_only; systematic read-kmer and assembly dependence excluded",
    }


def _negative_binomial_pmf(depths: Any, mean: float, dispersion: float) -> Any:
    from scipy.stats import nbinom

    probability = dispersion / (dispersion + mean)
    return nbinom.pmf(depths, dispersion, probability)


def estimate_kmer_heterozygosity(
    histogram: Mapping[int, int], *, k: int, minimum_depth: int = 5
) -> dict[str, Any]:
    """Fit a transparent two-copy/repeat mixture to an Illumina k-mer spectrum.

    The first component is the allele-specific (heterozygous) peak at half the
    shared-copy peak.  Components at 1.5x and 2x absorb low-order repeats.  The
    reported estimate is intentionally labelled as this model, not GenomeScope.
    """

    if k <= 0 or minimum_depth < 1:
        raise ValueError("invalid k-mer model settings")
    clean = {int(depth): int(count) for depth, count in histogram.items() if int(depth) >= minimum_depth and int(count) >= 0}
    if not clean or sum(clean.values()) == 0:
        return {"status": "not_estimable", "reason": "empty spectrum above minimum depth"}

    import numpy as np
    from scipy.optimize import nnls

    candidate_depths = [depth for depth in clean if depth >= max(10, minimum_depth)]
    if not candidate_depths:
        return {"status": "not_estimable", "reason": "no genomic peak candidate at depth >= 10"}
    dominant = max(candidate_depths, key=lambda depth: clean[depth])
    best: tuple[float, float, float, Any, Any, Any] | None = None
    mean_grid = np.linspace(max(10.0, dominant * 0.80), dominant * 1.20, 81)
    dispersion_grid = (30.0, 60.0, 120.0, 300.0, 1000.0, 1_000_000.0)
    for homo_mean in mean_grid:
        low_depth = max(minimum_depth, int(math.floor(0.30 * homo_mean)))
        high_depth = min(max(clean), int(math.ceil(2.60 * homo_mean)))
        depths = np.arange(low_depth, high_depth + 1, dtype=float)
        observed = np.array([clean.get(int(depth), 0) for depth in depths], dtype=float)
        positive = observed > 0
        if positive.sum() < 10:
            continue
        for dispersion in dispersion_grid:
            design = np.column_stack(
                [
                    _negative_binomial_pmf(depths, homo_mean * factor, dispersion)
                    for factor in (0.5, 1.0, 1.5, 2.0)
                ]
            )
            # Pearson-like weights prevent the large homozygous peak from
            # erasing the smaller heterozygous peak while retaining count scale.
            weights = 1.0 / np.sqrt(np.maximum(observed, 1.0))
            coefficients, _ = nnls(design * weights[:, None], observed * weights)
            fitted = design @ coefficients
            residual = float(np.sum((observed - fitted) ** 2 / np.maximum(observed, 1.0)))
            if best is None or residual < best[0]:
                best = (residual, homo_mean, dispersion, coefficients, observed, fitted)
    if best is None:
        return {"status": "not_estimable", "reason": "spectrum fit had too few populated bins"}
    residual, homo_mean, dispersion, coefficients, observed, fitted = best
    hetero_area, homo_area = float(coefficients[0]), float(coefficients[1])
    if hetero_area <= 0.0 or homo_area <= 0.0:
        return {"status": "not_estimable", "reason": "fitted heterozygous or homozygous component was zero"}
    q_no_heterozygote = 2.0 * homo_area / (hetero_area + 2.0 * homo_area)
    heterozygosity = 1.0 - q_no_heterozygote ** (1.0 / k)
    total_variation = float(((observed - observed.mean()) ** 2).sum())
    r_squared = 1.0 - float(((observed - fitted) ** 2).sum()) / total_variation if total_variation else 0.0
    return {
        "schema_version": "vgp-read-validation-kmer-heterozygosity-v1",
        "status": "estimated",
        "model": "four_component_negative_binomial_distinct_kmer_spectrum",
        "k": k,
        "minimum_fitted_depth": minimum_depth,
        "heterozygous_peak_depth": homo_mean / 2.0,
        "homozygous_peak_depth": homo_mean,
        "dispersion": dispersion,
        "heterozygous_distinct_kmer_component": hetero_area,
        "homozygous_distinct_kmer_component": homo_area,
        "repeat_1_5x_component": float(coefficients[2]),
        "repeat_2x_component": float(coefficients[3]),
        "heterozygosity_per_base": heterozygosity,
        "fit_r_squared": r_squared,
        "fit_pearson_residual": residual,
        "interpretation_limit": "model-based spectrum estimate; repeats, coverage bias, and correlated k-mers are systematic uncertainty",
    }


Variant = tuple[str, int, str, str]


def summarize_site_concordance(
    assembly_sites: set[Variant], read_sites: set[Variant], *, callable_bp: int
) -> dict[str, Any]:
    if callable_bp <= 0:
        raise ValueError("common callable denominator must be positive")
    shared = assembly_sites & read_sites
    union = assembly_sites | read_sites
    assembly_only = assembly_sites - read_sites
    read_only = read_sites - assembly_sites
    assembly_pi = len(assembly_sites) / callable_bp
    read_pi = len(read_sites) / callable_bp
    return {
        "callable_bp_common_mask": callable_bp,
        "assembly_sites": len(assembly_sites),
        "read_sites": len(read_sites),
        "shared_sites": len(shared),
        "assembly_only_sites": len(assembly_only),
        "read_only_sites": len(read_only),
        "assembly_site_recall_by_reads": len(shared) / len(assembly_sites) if assembly_sites else None,
        "read_site_precision_against_assembly": len(shared) / len(read_sites) if read_sites else None,
        "jaccard": len(shared) / len(union) if union else 1.0,
        "assembly_pi_common_mask": assembly_pi,
        "read_pi_common_mask": read_pi,
        "pi_difference_read_minus_assembly": read_pi - assembly_pi,
        "pi_ratio_read_over_assembly": read_pi / assembly_pi if assembly_pi else None,
        "concordant_pi_lower_bracket": len(shared) / callable_bp,
        "union_pi_upper_bracket": len(union) / callable_bp,
        "candidate_assembly_false_positive_upper_bound_sites": len(assembly_only),
        "candidate_assembly_false_positive_upper_bound_fraction_of_assembly_calls": (
            len(assembly_only) / len(assembly_sites) if assembly_sites else None
        ),
        "candidate_assembly_false_negative_upper_bound_sites": len(read_only),
        "candidate_assembly_false_negative_upper_bound_fraction_of_read_calls": (
            len(read_only) / len(read_sites) if read_sites else None
        ),
        "error_bound_scope": (
            "observed paired-call disagreement brackets, not complete biological confidence bounds; "
            "assembly-only includes read-caller misses and read-only includes assembly-caller misses"
        ),
        "method_covariance": "paired_shared_reference_and_callable_mask",
        "independence_status": "not_independent_replication",
    }


def _load_variant_tsv(path: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"chrom", "position_1based", "ref", "alt"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"variant TSV lacks {sorted(required)}: {path}")
        for row in reader:
            parsed: dict[str, Any] = {
                "chrom": row["chrom"],
                "position_1based": int(row["position_1based"]),
                "ref": row["ref"].upper(),
                "alt": row["alt"].upper(),
            }
            for key in ("quality", "depth"):
                if key in row and row[key] not in (None, "", "."):
                    parsed[key] = float(row[key]) if key == "quality" else int(row[key])
            if "genotype" in row:
                parsed["genotype"] = row["genotype"].replace("|", "/")
            if "allelic_depths" in row and row["allelic_depths"] not in (None, "", "."):
                parsed["allelic_depths"] = [int(value) for value in row["allelic_depths"].split(",") if value != "."]
            rows.append(parsed)
    return rows


def _quantile(values: Sequence[float], probability: float) -> float | None:
    if not values:
        return None
    ordered = sorted(values)
    position = probability * (len(ordered) - 1)
    lower = int(math.floor(position))
    upper = int(math.ceil(position))
    if lower == upper:
        return ordered[lower]
    fraction = position - lower
    return ordered[lower] * (1.0 - fraction) + ordered[upper] * fraction


def summarize_mask_variants(
    assembly_rows: Sequence[Mapping[str, Any]],
    read_rows: Sequence[Mapping[str, Any]],
    *,
    callable_bp: int,
) -> dict[str, Any]:
    """Summarize a common-mask assembly/read comparison and mapping QV bound."""

    assembly_sites: set[Variant] = {
        (str(row["chrom"]), int(row["position_1based"]), str(row["ref"]), str(row["alt"]))
        for row in assembly_rows
        if len(str(row["ref"])) == 1 and len(str(row["alt"])) == 1
    }
    read_hets: set[Variant] = set()
    homo_alt: set[Variant] = set()
    balances: list[float] = []
    for row in read_rows:
        if len(str(row["ref"])) != 1 or len(str(row["alt"])) != 1:
            continue
        genotype = str(row.get("genotype", ""))
        site = (str(row["chrom"]), int(row["position_1based"]), str(row["ref"]), str(row["alt"]))
        depths = list(row.get("allelic_depths", []))
        if genotype in {"0/1", "1/0"}:
            read_hets.add(site)
            if len(depths) >= 2 and depths[0] + depths[1] > 0:
                balances.append(depths[1] / (depths[0] + depths[1]))
        elif genotype == "1/1":
            homo_alt.add(site)
    concordance = summarize_site_concordance(assembly_sites, read_hets, callable_bp=callable_bp)
    error_low, error_high = wilson_interval(len(homo_alt), callable_bp)
    observed_error = len(homo_alt) / callable_bp
    concordance.update(
        {
            "strong_homozygous_alt_snp_sites": len(homo_alt),
            "mapping_consensus_qv": _error_to_qv(observed_error),
            "mapping_consensus_qv_lower_95": _error_to_qv(error_high),
            "mapping_consensus_qv_upper_95": _error_to_qv(error_low),
            "mapping_qv_scope": "homozygous-alt SNP discrepancies only; structural and inaccessible errors excluded",
            "read_het_alt_balance": {
                "count": len(balances),
                "mean": sum(balances) / len(balances) if balances else None,
                "q05": _quantile(balances, 0.05),
                "q25": _quantile(balances, 0.25),
                "median": _quantile(balances, 0.50),
                "q75": _quantile(balances, 0.75),
                "q95": _quantile(balances, 0.95),
                "fraction_0_3_to_0_7": (
                    sum(0.30 <= value <= 0.70 for value in balances) / len(balances) if balances else None
                ),
            },
        }
    )
    return concordance


def parse_psmc(path: Path) -> tuple[float, list[dict[str, float]]]:
    """Return the final completed PSMC optimization round only."""

    theta: float | None = None
    rows: list[dict[str, float]] = []
    for raw in path.read_text().splitlines():
        fields = raw.split("\t")
        if fields[0] == "RD":
            # Each optimization round carries a complete TR/RS trajectory.
            # Reset here so intermediate iterations cannot inflate the 64-bin
            # biological curve or dominate between-method correlations.
            theta = None
            rows = []
        elif fields[0] == "TR" and len(fields) >= 2:
            theta = float(fields[1])
        elif fields[0] == "RS" and len(fields) >= 4:
            rows.append({"interval": int(fields[1]), "time_2N0": float(fields[2]), "lambda": float(fields[3])})
    if theta is None or not rows:
        raise ValueError(f"not a finite PSMC output: {path}")
    return theta, rows


def compare_psmc(assembly_path: Path, read_path: Path) -> dict[str, Any]:
    import numpy as np

    assembly_theta, assembly = parse_psmc(assembly_path)
    read_theta, reads = parse_psmc(read_path)
    if len(assembly) != len(reads) or [row["interval"] for row in assembly] != [row["interval"] for row in reads]:
        raise ValueError("PSMC interval grids differ")
    assembly_lambda = np.array([row["lambda"] for row in assembly], dtype=float)
    read_lambda = np.array([row["lambda"] for row in reads], dtype=float)
    assembly_time = np.array([row["time_2N0"] for row in assembly], dtype=float)
    read_time = np.array([row["time_2N0"] for row in reads], dtype=float)
    if np.any(assembly_lambda <= 0) or np.any(read_lambda <= 0):
        raise ValueError("PSMC lambda must be positive")
    return {
        "assembly_theta_0": assembly_theta,
        "read_theta_0": read_theta,
        "theta_ratio_read_over_assembly": read_theta / assembly_theta,
        "theta_difference_read_minus_assembly": read_theta - assembly_theta,
        "intervals": len(assembly),
        "lambda_pearson_correlation": float(np.corrcoef(assembly_lambda, read_lambda)[0, 1]),
        "log_lambda_rmse": float(np.sqrt(np.mean((np.log(read_lambda) - np.log(assembly_lambda)) ** 2))),
        "time_2N0_pearson_correlation": float(np.corrcoef(assembly_time, read_time)[0, 1]),
        "method_covariance": "same_H1_reference_PSMC_parameterization_and_overlapping_callable_sequence",
        "independence_status": "not_independent_replication",
    }


def _merged_interval_bp(intervals: Sequence[tuple[int, int]]) -> int:
    covered = 0
    active_start = active_end = None
    for start, end in sorted(intervals):
        if active_start is None:
            active_start, active_end = start, end
        elif start <= active_end:
            active_end = max(active_end, end)
        else:
            covered += active_end - active_start
            active_start, active_end = start, end
    if active_start is not None:
        covered += active_end - active_start
    return covered


def summarize_paf(
    rows: Iterable[str], *, reference_bp: int, raw_bases: int, raw_reads: int | None = None
) -> dict[str, Any]:
    """Summarize primary PAF mappings without treating low coverage as callability."""

    alignments = 0
    mapped_query_segment_bases = 0
    alignment_block_bases = 0
    matching_bases = 0
    mapq_sum = 0
    mapq_20 = 0
    target_intervals: dict[str, list[tuple[int, int]]] = {}
    query_intervals: dict[str, list[tuple[int, int]]] = {}
    for line_number, raw in enumerate(rows, 1):
        if not raw.strip():
            continue
        fields = raw.rstrip("\n").split("\t")
        if len(fields) < 12:
            raise ValueError(f"PAF line {line_number} has fewer than 12 columns")
        query, query_length = fields[0], int(fields[1])
        query_start, query_end = int(fields[2]), int(fields[3])
        target, target_start, target_end = fields[5], int(fields[7]), int(fields[8])
        matches, block, mapq = int(fields[9]), int(fields[10]), int(fields[11])
        if min(query_length, query_start, query_end, target_start, target_end, matches, block, mapq) < 0:
            raise ValueError(f"negative PAF field at line {line_number}")
        alignments += 1
        mapped_query_segment_bases += query_end - query_start
        alignment_block_bases += block
        matching_bases += matches
        mapq_sum += mapq
        mapq_20 += mapq >= 20
        target_intervals.setdefault(target, []).append((target_start, target_end))
        query_intervals.setdefault(query, []).append((query_start, query_end))
    covered = sum(_merged_interval_bp(intervals) for intervals in target_intervals.values())
    mapped_query_union = sum(_merged_interval_bp(intervals) for intervals in query_intervals.values())
    return {
        "schema_version": "vgp-read-validation-paf-summary-v2",
        "canonical_vgp_root": str(CANONICAL_VGP_ROOT),
        "non_secondary_alignment_segments": alignments,
        "mapped_queries": len(query_intervals),
        "raw_reads_metadata": raw_reads,
        "mapped_query_fraction_of_raw_reads": len(query_intervals) / raw_reads if raw_reads else None,
        "raw_bases_metadata": raw_bases,
        "reference_bp": reference_bp,
        "nominal_physical_coverage": raw_bases / reference_bp if reference_bp else None,
        "nominal_diploid_equivalent_coverage": raw_bases / (2 * reference_bp) if reference_bp else None,
        "mapped_query_bases_segment_sum": mapped_query_segment_bases,
        "mapped_query_bases_union": mapped_query_union,
        "mapped_query_fraction_of_raw_bases": mapped_query_union / raw_bases if raw_bases else None,
        "weighted_alignment_identity": matching_bases / alignment_block_bases if alignment_block_bases else None,
        "mean_segment_mapping_quality": mapq_sum / alignments if alignments else None,
        "fraction_alignments_mapq_at_least_20": mapq_20 / alignments if alignments else None,
        "reference_breadth_at_least_one_mapping_bp": covered,
        "reference_breadth_at_least_one_mapping_fraction": covered / reference_bp if reference_bp else None,
        "callability_status": "not_established_from_single_coverage_breadth",
    }


def _load_histogram(path: Path) -> dict[int, int]:
    result: dict[int, int] = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            result[int(row["depth"])] = int(row["count"] if "count" in row else row["bases"])
    return result


def _parse_mask_spec(text: str) -> DepthMaskSpec:
    try:
        name, minimum, maximum = text.split(":")
        return DepthMaskSpec(name, int(minimum), int(maximum))
    except (ValueError, TypeError) as error:
        raise argparse.ArgumentTypeError("mask must be NAME:MIN:MAX") from error


def _command_depth_masks(args: argparse.Namespace) -> None:
    stream_depth_masks(sys.stdin, args.output_dir, args.mask)


def _command_assembly_evidence(args: argparse.Namespace) -> None:
    with args.assembly_sites.open() as sites, args.pileup.open() as pileup, args.output.open("w") as output:
        summary = assembly_evidence(
            sites,
            pileup,
            output,
            minimum_depth=args.minimum_depth,
            maximum_depth=args.maximum_depth,
        )
    _atomic_json(args.summary, summary)


def _command_kmer_model(args: argparse.Namespace) -> None:
    result = estimate_kmer_heterozygosity(
        _load_histogram(args.histogram), k=args.k, minimum_depth=args.minimum_depth
    )
    result["canonical_vgp_root"] = str(CANONICAL_VGP_ROOT)
    _atomic_json(args.output, result)


def _command_kmer_qv(args: argparse.Namespace) -> None:
    total = error = 0
    query_histogram: Counter[int] = Counter()
    samples: list[str] = []
    for line_number, raw in enumerate(sys.stdin, 1):
        if not raw.strip():
            continue
        fields = raw.split()
        if len(fields) != 2:
            raise ValueError(f"Jellyfish query line {line_number} does not have k-mer and count")
        count = int(fields[1])
        total += 1
        query_histogram[count] += 1
        if count < args.trusted_minimum:
            error += 1
            if len(samples) < args.sample_limit:
                samples.append(fields[0])
    result = kmer_qv(total_kmers=total, error_kmers=error, k=args.k)
    result.update(
        {
            "canonical_vgp_root": str(CANONICAL_VGP_ROOT),
            "trusted_read_kmer_minimum_count": args.trusted_minimum,
            "query_count_histogram": {str(key): query_histogram[key] for key in sorted(query_histogram)},
            "sample_untrusted_assembly_kmers": samples,
        }
    )
    _atomic_json(args.output, result)


def _command_fasta_to_fastq(args: argparse.Namespace) -> None:
    name: str | None = None
    sequence: list[str] = []

    def emit() -> None:
        if name is None:
            return
        joined = "".join(sequence).upper()
        sys.stdout.write(f"@{name}\n{joined}\n+\n{'I' * len(joined)}\n")

    for raw in sys.stdin:
        if raw.startswith(">"):
            emit()
            name = raw[1:].strip().split()[0]
            sequence = []
        else:
            sequence.append(raw.strip())
    emit()


def _command_mask_report(args: argparse.Namespace) -> None:
    result = summarize_mask_variants(
        _load_variant_tsv(args.assembly_variants),
        _load_variant_tsv(args.read_variants),
        callable_bp=args.callable_bp,
    )
    result.update(
        {
            "schema_version": "vgp-read-validation-common-mask-v1",
            "canonical_vgp_root": str(CANONICAL_VGP_ROOT),
            "mask_id": args.mask_id,
        }
    )
    _atomic_json(args.output, result)


def _command_psmc_compare(args: argparse.Namespace) -> None:
    result = compare_psmc(args.assembly_psmc, args.read_psmc)
    result.update(
        {
            "schema_version": "vgp-read-validation-psmc-comparison-v1",
            "canonical_vgp_root": str(CANONICAL_VGP_ROOT),
            "mask_id": args.mask_id,
        }
    )
    _atomic_json(args.output, result)


def _command_paf_summary(args: argparse.Namespace) -> None:
    with args.paf.open() as handle:
        result = summarize_paf(
            handle, reference_bp=args.reference_bp, raw_bases=args.raw_bases, raw_reads=args.raw_reads
        )
    result["selection_id"] = args.selection_id
    _atomic_json(args.output, result)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    depth = subparsers.add_parser("depth-masks", help="stream samtools depth into sensitivity BED masks")
    depth.add_argument("--output-dir", type=Path, required=True)
    depth.add_argument("--mask", action="append", type=_parse_mask_spec, required=True)
    depth.set_defaults(function=_command_depth_masks)

    evidence = subparsers.add_parser("assembly-evidence", help="classify pileup support at assembly SNPs")
    evidence.add_argument("--assembly-sites", type=Path, required=True)
    evidence.add_argument("--pileup", type=Path, required=True)
    evidence.add_argument("--output", type=Path, required=True)
    evidence.add_argument("--summary", type=Path, required=True)
    evidence.add_argument("--minimum-depth", type=int, required=True)
    evidence.add_argument("--maximum-depth", type=int, required=True)
    evidence.set_defaults(function=_command_assembly_evidence)

    model = subparsers.add_parser("kmer-model", help="fit a two-copy k-mer heterozygosity model")
    model.add_argument("--histogram", type=Path, required=True)
    model.add_argument("--k", type=int, default=21)
    model.add_argument("--minimum-depth", type=int, default=5)
    model.add_argument("--output", type=Path, required=True)
    model.set_defaults(function=_command_kmer_model)

    qv = subparsers.add_parser("kmer-qv", help="summarize streamed Jellyfish assembly queries")
    qv.add_argument("--k", type=int, default=21)
    qv.add_argument("--trusted-minimum", type=int, default=5)
    qv.add_argument("--sample-limit", type=int, default=1000)
    qv.add_argument("--output", type=Path, required=True)
    qv.set_defaults(function=_command_kmer_qv)

    fastq = subparsers.add_parser("fasta-to-fastq", help="emit constant-Q40 FASTQ for fq2psmcfa")
    fastq.set_defaults(function=_command_fasta_to_fastq)

    report = subparsers.add_parser("mask-report", help="compare assembly/read variants on one exact mask")
    report.add_argument("--assembly-variants", type=Path, required=True)
    report.add_argument("--read-variants", type=Path, required=True)
    report.add_argument("--callable-bp", type=int, required=True)
    report.add_argument("--mask-id", required=True)
    report.add_argument("--output", type=Path, required=True)
    report.set_defaults(function=_command_mask_report)

    psmc = subparsers.add_parser("psmc-compare", help="compare paired unscaled PSMC trajectories")
    psmc.add_argument("--assembly-psmc", type=Path, required=True)
    psmc.add_argument("--read-psmc", type=Path, required=True)
    psmc.add_argument("--mask-id", required=True)
    psmc.add_argument("--output", type=Path, required=True)
    psmc.set_defaults(function=_command_psmc_compare)

    paf = subparsers.add_parser("paf-summary", help="summarize primary low-coverage read mappings")
    paf.add_argument("--paf", type=Path, required=True)
    paf.add_argument("--reference-bp", type=int, required=True)
    paf.add_argument("--raw-bases", type=int, required=True)
    paf.add_argument("--raw-reads", type=int)
    paf.add_argument("--selection-id", required=True)
    paf.add_argument("--output", type=Path, required=True)
    paf.set_defaults(function=_command_paf_summary)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    args.function(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
