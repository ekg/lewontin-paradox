#!/usr/bin/env python3
"""Independently recompute the completed real-VGP pilot evidence.

This reviewer deliberately does not import the production pipeline.  It reads
the canonical BED, VCF, FASTA, PSMC, and sacct products directly and emits one
machine-readable review packet.  ``VGP_ROOT`` is the single configuration
variable from which every canonical data path is resolved.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import gzip
import json
import math
import os
import statistics
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable, Mapping, Sequence


DEFAULT_VGP_ROOT = Path(os.environ.get("VGP_ROOT", "/moosefs/erikg/vgp"))
AUTHORIZATION_ID = "vgp10-auth-20260718-v2"
PAIR_IDS = tuple(f"P{number:02d}" for number in range(1, 11))
IUPAC_HET = frozenset("RYSWKM")


class ReviewError(RuntimeError):
    """The independently observed packet violates a hard review invariant."""


def read_json(path: Path) -> dict[str, object]:
    with path.open(encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ReviewError(f"expected JSON object: {path}")
    return value


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def merge_intervals(rows: Iterable[tuple[str, int, int]]) -> list[tuple[str, int, int]]:
    merged: list[tuple[str, int, int]] = []
    for contig, start, end in sorted(rows):
        if start < 0 or end <= start:
            raise ReviewError(f"invalid interval {contig}:{start}-{end}")
        if merged and merged[-1][0] == contig and start <= merged[-1][2]:
            old_contig, old_start, old_end = merged[-1]
            merged[-1] = (old_contig, old_start, max(old_end, end))
        else:
            merged.append((contig, start, end))
    return merged


def read_bed(path: Path) -> list[tuple[str, int, int]]:
    rows: list[tuple[str, int, int]] = []
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise ReviewError(f"{path}:{number}: short BED row")
            rows.append((fields[0], int(fields[1]), int(fields[2])))
    return merge_intervals(rows)


def interval_index(
    rows: Sequence[tuple[str, int, int]],
) -> dict[str, list[tuple[int, int]]]:
    result: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for contig, start, end in rows:
        result[contig].append((start, end))
    return dict(result)


def contains(
    index: Mapping[str, Sequence[tuple[int, int]]], contig: str, start: int, end: int
) -> bool:
    values = index.get(contig, ())
    position = bisect.bisect_right(values, (start, 2**63 - 1)) - 1
    return position >= 0 and values[position][0] <= start and end <= values[position][1]


def intersect_one(
    index: Mapping[str, Sequence[tuple[int, int]]], contig: str, start: int, end: int
) -> list[tuple[str, int, int]]:
    result: list[tuple[str, int, int]] = []
    values = index.get(contig, ())
    position = max(0, bisect.bisect_right(values, (start, 2**63 - 1)) - 1)
    while position < len(values) and values[position][0] < end:
        left, right = values[position]
        overlap_start, overlap_end = max(start, left), min(end, right)
        if overlap_start < overlap_end:
            result.append((contig, overlap_start, overlap_end))
        position += 1
    return result


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open(encoding="utf-8")


def vcf_records(path: Path):
    with open_text(path) as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                raise ReviewError(f"{path}:{number}: short VCF row")
            ref = fields[3].upper()
            alts = fields[4].upper().split(",")
            if len(alts) != 1 or not ref or not alts[0] or alts[0] == ".":
                raise ReviewError(f"{path}:{number}: expected one exact alternate")
            yield fields[0], int(fields[1]) - 1, ref, alts[0]


def audit_variants(
    vcf: Path, callable_rows: Sequence[tuple[str, int, int]], indel_flank: int = 10
) -> dict[str, object]:
    """Reapply the callable and indel-flank policies in two streaming passes."""

    callable_index = interval_index(callable_rows)
    normalized = Counter()
    callable_counts = Counter()
    indel_windows: list[tuple[str, int, int]] = []
    for contig, pos0, ref, alt in vcf_records(vcf):
        kind = "snp" if len(ref) == len(alt) == 1 else "indel"
        normalized[kind] += 1
        if not contains(callable_index, contig, pos0, pos0 + len(ref)):
            continue
        callable_counts[kind] += 1
        if kind == "indel":
            indel_windows.extend(
                intersect_one(
                    callable_index,
                    contig,
                    max(0, pos0 - indel_flank),
                    pos0 + len(ref) + indel_flank,
                )
            )
    merged_indel = merge_intervals(indel_windows)
    indel_index = interval_index(merged_indel)
    masked_snps = 0
    for contig, pos0, ref, alt in vcf_records(vcf):
        if (
            len(ref) == len(alt) == 1
            and contains(callable_index, contig, pos0, pos0 + 1)
            and contains(indel_index, contig, pos0, pos0 + 1)
        ):
            masked_snps += 1
    pre_callable = sum(end - start for _, start, end in callable_rows)
    indel_masked_bp = sum(end - start for _, start, end in merged_indel)
    final_callable = pre_callable - indel_masked_bp
    final_snps = callable_counts["snp"] - masked_snps
    return {
        "method": (
            "two-pass VCF recount plus BED containment and independent merged "
            f"+/-{indel_flank}-bp callable indel mask"
        ),
        "normalized_variant_records": sum(normalized.values()),
        "normalized_snp_records": normalized["snp"],
        "normalized_indel_records": normalized["indel"],
        "pre_indel_callable_bp": pre_callable,
        "callable_variant_records": sum(callable_counts.values()),
        "callable_snp_records": callable_counts["snp"],
        "callable_indel_records": callable_counts["indel"],
        "callable_snps_masked_by_indel_flank": masked_snps,
        "indel_masked_h1_bp": indel_masked_bp,
        "final_callable_bp": final_callable,
        "heterozygous_snps": final_snps,
        "pi": final_snps / final_callable,
    }


def fasta_symbols(path: Path, allowed: frozenset[str]) -> dict[str, int]:
    counts: Counter[str] = Counter()
    with open_text(path) as handle:
        for number, line in enumerate(handle, 1):
            if line.startswith(">") or not line.strip():
                continue
            sequence = line.strip().upper()
            unexpected = set(sequence) - allowed
            if unexpected:
                raise ReviewError(f"{path}:{number}: unexpected symbols {sorted(unexpected)}")
            counts.update(sequence)
    return dict(sorted(counts.items()))


def parse_psmc(path: Path) -> tuple[list[dict[str, float]], float]:
    by_iteration: dict[int, list[dict[str, float]]] = defaultdict(list)
    theta: dict[int, float] = {}
    current: int | None = None
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.split()
            if not fields:
                continue
            if fields[0] == "RD" and len(fields) >= 2:
                current = int(fields[1])
            elif fields[0] == "RS" and current is not None and len(fields) >= 4:
                by_iteration[current].append(
                    {"interval": int(fields[1]), "time_2N0": float(fields[2]), "lambda": float(fields[3])}
                )
            elif fields[0] == "TR" and current is not None and len(fields) >= 3:
                theta[current] = float(fields[1])
    if not by_iteration:
        raise ReviewError(f"PSMC has no trajectory: {path}")
    iteration = max(by_iteration)
    rows = sorted(by_iteration[iteration], key=lambda row: row["interval"])
    if iteration not in theta or theta[iteration] <= 0:
        raise ReviewError(f"PSMC has no positive final theta: {path}")
    if len({row["interval"] for row in rows}) != len(rows):
        raise ReviewError(f"PSMC repeats final intervals: {path}")
    if any(not math.isfinite(value) for row in rows for value in row.values()):
        raise ReviewError(f"PSMC has a non-finite final trajectory: {path}")
    return rows, theta[iteration]


def quantiles(values: Sequence[float]) -> dict[str, float]:
    ordered = sorted(values)
    if not ordered:
        raise ReviewError("cannot summarize an empty distribution")

    def pick(fraction: float) -> float:
        return ordered[round((len(ordered) - 1) * fraction)]

    return {
        "min": ordered[0],
        "q025": pick(0.025),
        "median": statistics.median(ordered),
        "q975": pick(0.975),
        "max": ordered[-1],
    }


def audit_psmc(root: Path) -> dict[str, object]:
    primary_rows, primary_theta = parse_psmc(root / "psmc/replicate-000/unscaled.psmc")
    replicate_paths = sorted(root.glob("psmc/replicate-[0-9][0-9][0-9]/bootstrap.unscaled.psmc"))
    thetas: list[float] = []
    interval_counts: Counter[int] = Counter()
    bootstrap_max_lambdas: list[float] = []
    bootstrap_oldest_lambdas: list[float] = []
    finite = 0
    for path in replicate_paths:
        rows, theta = parse_psmc(path)
        finite += 1
        interval_counts[len(rows)] += 1
        thetas.append(theta)
        bootstrap_max_lambdas.append(max(row["lambda"] for row in rows))
        bootstrap_oldest_lambdas.append(rows[-1]["lambda"])
    lambdas = [row["lambda"] for row in primary_rows]
    times = [row["time_2N0"] for row in primary_rows]
    scenario_rows = read_tsv(root / "psmc/finalize/scenario_scaled_trajectories.tsv")
    scenario_ids = sorted({row["scenario_id"] for row in scenario_rows})
    sources = sorted({row["mutation_rate_source"] for row in scenario_rows})
    theta_summary = quantiles(thetas)
    return {
        "method": "native final-iteration RD/TR/RS parse of primary and every bootstrap; no production helper imported",
        "primary": {
            "theta_0_per_100bp_bin": primary_theta,
            "intervals": len(primary_rows),
            "time_2N0_range": [min(times), max(times)],
            "lambda_range": [min(lambdas), max(lambdas)],
            "recent_lambda": lambdas[0],
            "oldest_lambda": lambdas[-1],
            "minimum_lambda_interval": int(primary_rows[lambdas.index(min(lambdas))]["interval"]),
            "maximum_lambda_interval": int(primary_rows[lambdas.index(max(lambdas))]["interval"]),
            "distinct_lambda_values": len(set(lambdas)),
        },
        "bootstrap_attempts": len(replicate_paths),
        "finite_bootstraps": finite,
        "bootstrap_interval_counts": {str(key): value for key, value in sorted(interval_counts.items())},
        "bootstrap_theta_0_per_100bp_bin": theta_summary,
        "bootstrap_curve_behavior": {
            "maximum_lambda_distribution": quantiles(bootstrap_max_lambdas),
            "oldest_interval_lambda_distribution": quantiles(bootstrap_oldest_lambdas),
            "primary_theta_to_bootstrap_median_ratio": primary_theta / theta_summary["median"],
            "primary_theta_outside_observed_bootstrap_range": not (
                theta_summary["min"] <= primary_theta <= theta_summary["max"]
            ),
            "validity": "STATISTICAL_FAIL",
            "reason": (
                "bootstrap inputs sample callable BED units and therefore do not preserve the primary PSMCFA's "
                "masked-bin population; numerical finiteness cannot establish centered uncertainty"
            ),
            "interpretation": (
                "lambda distributions are reported descriptively only and are not admitted as uncertainty "
                "bands because replicate theta is systematically shifted from the primary"
            ),
        },
        "scenario_count": len(scenario_ids),
        "scenario_rows": len(scenario_rows),
        "scenario_ids": scenario_ids,
        "scaling_sources": sources,
    }


def request_memory_gib(value: str) -> float | None:
    if not value:
        return None
    per_cpu = value[-1:] in {"c", "n"}
    if per_cpu:
        value = value[:-1]
    units = {"K": 1 / 1024**2, "M": 1 / 1024, "G": 1.0, "T": 1024.0}
    suffix = value[-1].upper()
    if suffix not in units:
        raise ReviewError(f"unsupported ReqMem value: {value}")
    return float(value[:-1]) * units[suffix]


def allocation_summary(rows: Sequence[Mapping[str, str]]) -> dict[str, object]:
    allocated = [row for row in rows if int(row.get("AllocCPUS") or 0) > 0]
    core_seconds = sum(int(row.get("CPUTimeRAW") or 0) for row in allocated)
    memory_gib_seconds = 0.0
    memory_rows = 0
    for row in allocated:
        memory = request_memory_gib(row.get("ReqMem", ""))
        if memory is not None:
            if row["ReqMem"].endswith("c"):
                memory *= int(row["AllocCPUS"])
            memory_gib_seconds += memory * int(row["ElapsedRaw"])
            memory_rows += 1
    return {
        "allocations": len(allocated),
        "elapsed_allocation_hours": sum(int(row["ElapsedRaw"]) for row in allocated) / 3600,
        "allocated_core_hours": core_seconds / 3600,
        "requested_memory_gib_hours": memory_gib_seconds / 3600,
        "requested_memory_rows": memory_rows,
        "maxrss_observed_rows": sum(bool(row.get("MaxRSS")) for row in allocated),
        "state_counts": dict(sorted(Counter(row["State"].split()[0] for row in rows).items())),
    }


def p04_outcome_allocations(rows: Sequence[Mapping[str, str]]) -> dict[str, object]:
    fixed = {"1782075", "1783681", "1783682", "1783683", "1783685"}
    chosen = [
        row
        for row in rows
        if row["State"].startswith("COMPLETED")
        and (row["JobIDRaw"] in fixed or (row["JobName"] == "vgp10-P04-psmc" and int(row["AllocCPUS"] or 0) == 4))
    ]
    by_stage: dict[str, list[Mapping[str, str]]] = defaultdict(list)
    for row in chosen:
        by_stage[row["JobName"].rsplit("-", 1)[-1]].append(row)
    result = {stage: allocation_summary(stage_rows) for stage, stage_rows in sorted(by_stage.items())}
    psmc_elapsed = [int(row["ElapsedRaw"]) for row in by_stage["psmc"]]
    result["psmc"]["elapsed_seconds_distribution"] = quantiles(psmc_elapsed)
    result["psmc"]["scaleout_model"] = {
        "allocations_per_pair": 201,
        "cpus_per_allocation": 4,
        "requested_memory_gib_per_allocation": 16,
        "observed_median_minutes": statistics.median(psmc_elapsed) / 60,
        "projected_makespan_hours_at_throttle_21": math.ceil(201 / 21) * statistics.median(psmc_elapsed) / 3600,
    }
    result["total"] = allocation_summary(chosen)
    return result


def audit_pair(pair: str, root: Path) -> dict[str, object]:
    if not str(root.resolve()).startswith(str(DEFAULT_VGP_ROOT.resolve()) + os.sep):
        raise ReviewError(f"{pair} root escapes canonical VGP root: {root}")
    callable_rows = read_bed(root / "consensus/masks/callable.bed")
    variants = audit_variants(root / "consensus/normalized.vcf", callable_rows)
    mask = read_json(root / "consensus/masks/mask_reconciliation.json")
    if int(mask["callable_bp"]) != variants["pre_indel_callable_bp"]:
        raise ReviewError(f"{pair} callable BED disagrees with reconciliation JSON")
    if int(mask["universe_bp"]) != int(mask["callable_bp"]) + sum(
        int(value) for value in mask["excluded_bp_by_primary_reason"].values()
    ):
        raise ReviewError(f"{pair} callable plus reasons does not equal the H1 universe")
    consensus = fasta_symbols(root / "consensus/consensus/consensus.fa", frozenset("ACGTN") | IUPAC_HET)
    psmcfa = fasta_symbols(root / "consensus/consensus/input.psmcfa", frozenset("NKT"))
    consensus_het = sum(consensus.get(symbol, 0) for symbol in IUPAC_HET)
    consensus_callable = sum(consensus.get(symbol, 0) for symbol in "ACGT") + consensus_het
    if consensus_het != variants["heterozygous_snps"] or consensus_callable != variants["final_callable_bp"]:
        raise ReviewError(f"{pair} consensus symbols disagree with independently selected variants")
    return {
        "selection_id": pair,
        "canonical_run_root": str(root),
        "universe_bp": int(mask["universe_bp"]),
        "pre_indel_callable_fraction": variants["pre_indel_callable_bp"] / int(mask["universe_bp"]),
        "final_callable_fraction": variants["final_callable_bp"] / int(mask["universe_bp"]),
        "accounting_discrepancy_bp": int(mask["accounting_discrepancy_bp"]),
        "variant_recomputation": variants,
        "consensus_recomputation": {
            "symbols": consensus,
            "callable_bases": consensus_callable,
            "heterozygous_iupac_bases": consensus_het,
            "masked_N_bases": consensus.get("N", 0),
            "psmcfa_symbols": psmcfa,
        },
        "psmc_recomputation": audit_psmc(root),
    }


def build_review(args: argparse.Namespace) -> dict[str, object]:
    global DEFAULT_VGP_ROOT
    DEFAULT_VGP_ROOT = args.canonical_root
    closed = read_json(args.closed_world)
    if closed.get("canonical_vgp_root") != str(args.canonical_root):
        raise ReviewError("closed-world packet points outside the configured canonical root")
    pairs = {row["selection_id"]: row for row in closed["pairs"]}
    if set(pairs) != set(PAIR_IDS):
        raise ReviewError("closed world is not exactly P01-P10")
    completed = {
        "P04": audit_pair("P04", args.canonical_root / "pilot/runs/vgp10-auth-20260718-v2-pilot-v1/P04"),
        "P07": audit_pair("P07", args.canonical_root / "pilot/outputs/vgp10-auth-20260718-v2/P07/core"),
    }
    pilot_sacct = read_tsv(args.pilot_sacct)
    canary_sacct = read_tsv(args.canary_sacct)
    live = read_tsv(args.live_sacct)
    # The older canary sacct table predates the canonical-root column.  Its
    # execution manifest binds the same paths; all newly emitted review rows
    # and the newer pilot/live tables must carry the root explicitly.
    for row in [*pilot_sacct, *live]:
        if row.get("canonical_vgp_root") != str(args.canonical_root):
            raise ReviewError("sacct row points outside the configured canonical root")
    pi_values = [pair["variant_recomputation"]["pi"] for pair in completed.values()]
    callable_values = [pair["variant_recomputation"]["final_callable_bp"] for pair in completed.values()]
    snp_values = [pair["variant_recomputation"]["heterozygous_snps"] for pair in completed.values()]
    hard_failures = {
        pair: pairs[pair]["terminal_failure"]
        for pair in ("P01", "P02", "P03", "P05")
    }
    live_by_pair = {row["JobName"].split("-")[1]: row for row in live}
    dispositions = {
        "P01": "hard_invalid_primary_impg_dictionary_failure",
        "P02": "hard_invalid_primary_impg_dictionary_failure",
        "P03": "hard_invalid_primary_mapping_execution_failure",
        "P04": "valid_core_complete_annotation_not_applicable",
        "P05": "hard_invalid_primary_impg_dictionary_and_mixed_site_failure",
        "P06": "no_biological_result_corrected_mapping_failed_cause_not_revalidated",
        "P07": "valid_core_complete_exact_annotation_confidence_subset",
        "P08": "no_biological_result_mapping_running_at_cutoff",
        "P09": "no_biological_result_mapping_running_at_cutoff",
        "P10": "no_biological_result_mapping_running_at_cutoff",
    }
    return {
        "schema_version": "vgp-real-pilot-independent-review-v1",
        "task_id": "review-vgp-real-pilot",
        "authorization_id": AUTHORIZATION_ID,
        "canonical_vgp_root": str(args.canonical_root),
        "review_cutoff_utc": live[0]["review_cutoff_utc"],
        "decision": "CONDITIONAL_GO",
        "decision_scope": "bounded core-only scale-out; annotation and raw-read validation have independent confidence decisions",
        "non_retrospective_core_criteria": {
            "exact_same_individual_pair_and_frozen_digests": True,
            "exact_1to1_mapping_and_reconstruction": True,
            "callable_bp_minimum": 100_000_000,
            "pre_indel_callable_fraction_minimum": 0.60,
            "mask_accounting_discrepancy_bp_required": 0,
            "minimum_finite_bootstraps": 190,
            "bootstrap_attempts": 200,
            "unscaled_primary_required": True,
            "annotation_required_for_core": False,
            "raw_read_validation_required_for_core": False,
        },
        "numbers_first": {
            "completed_pairs": 2,
            "completed_pair_ids": ["P04", "P07"],
            "final_callable_bp_total": sum(callable_values),
            "final_callable_bp_range": [min(callable_values), max(callable_values)],
            "heterozygous_snps_total": sum(snp_values),
            "heterozygous_snps_range": [min(snp_values), max(snp_values)],
            "pi_range": [min(pi_values), max(pi_values)],
            "finite_bootstraps": 400,
            "bootstrap_attempts": 400,
            "numerically_finite_bootstraps": 400,
            "bootstrap_uncertainty_sets_admitted": 0,
            "hard_invalid_pairs": 4,
            "incomplete_without_biological_result_at_cutoff": 4,
        },
        "completed_pairs": completed,
        "pair_dispositions": dispositions,
        "hard_failures_from_frozen_ledger": hard_failures,
        "post_upstream_live_delta": {
            pair: {
                "job_id": row["JobIDRaw"],
                "state": row["State"],
                "elapsed_seconds": int(row["ElapsedRaw"]),
                "exit_code": row["ExitCode"],
                "cause": (
                    "not_revalidated_from_empty_job_stdout/stderr; do not infer from prior attempts"
                    if pair == "P06"
                    else "not_applicable_job_running_at_cutoff"
                ),
            }
            for pair, row in sorted(live_by_pair.items())
            if pair in {"P06", "P08", "P09", "P10"}
        },
        "resource_model_from_sacct": {
            "measurement_boundary": "elapsed allocation capacity and requested memory; MaxRSS/TotalCPU unavailable on this cluster and never imputed",
            "pilot_all_attempts_snapshot": allocation_summary(pilot_sacct),
            "canary_all_attempts_snapshot": allocation_summary(canary_sacct),
            "P04_outcome_producing_allocations": p04_outcome_allocations(pilot_sacct),
            "P07_final_core_allocation": allocation_summary(
                [row for row in canary_sacct if row["JobIDRaw"] == "1781798"]
            ),
            "P07_successful_rescue_allocations": allocation_summary(
                [row for row in canary_sacct if row["JobIDRaw"] in {"1781772", "1781779"}]
            ),
            "corrected_scaleout_rules": [
                "budget PSMC as 201 independent 4-CPU/16-GiB allocations per completed pair, not one array allocation",
                "do not spend the 200-bootstrap allocation budget at scale until the shifted bootstrap input population is repaired and revalidated",
                "for a P04-sized 2.372-Gbp pair, use measured stage allocations: mapping 32 CPU/160 GiB for 895 s; IMPG 48 CPU/160 GiB for 23,349 s; variants 24 CPU/160 GiB for 175 s; consensus 24 CPU/160 GiB for 1,820 s",
                "retain the measured 600-GiB free-scratch gate for P04-like IMPG and size-stratified per-pair envelopes for larger genomes; sacct does not measure scratch",
                "run a bounded wave and recalibrate after each new completed size stratum; two completions do not identify a reliable genome-size scaling slope",
                "exclude failed and cancelled allocations from the success-template model but report their capacity cost separately",
            ],
        },
        "confidence_decisions": {
            "core": "PASS for P04 and P07; hard fail or no result for all other pairs",
            "psmc": "DESCRIPTIVE_PRIMARY_ONLY: both primaries are finite; bootstrap uncertainty is STATISTICAL_FAIL because the replicate input population is shifted from the primary",
            "annotation": "PASS only for P07 exact-native partitions; NOT_APPLICABLE for P04; not a core veto",
            "raw_reads": "UNVALIDATED for both completed pairs at this cutoff; confidence covariate only",
            "assembly_qv_busco_collapse_repeat_covariates": "UNVALIDATED; tier-C confidence, not hard invalidity",
            "population_inference": "NOT_AUTHORIZED: each estimate is one diploid individual's phased-haplotype divergence",
        },
        "conditions_for_scaleout": [
            "preflight graph/FASTA sequence dictionaries and mixed-site handling before expensive IMPG extraction; never omit failed regions",
            "preserve exact pair provenance, 1:1 multiplicity, reconstruction, mask accounting, callable thresholds, and >=190/200 finite PSMC bootstraps",
            "before admitting PSMC uncertainty, rebuild bootstraps from blocks of the primary PSMCFA so masked and callable bins share the primary sampling population, then require primary theta/trajectory centering diagnostics",
            "use canonical-root-only paths and fail-closed private node-local scratch guards",
            "limit the first wave to resource strata represented by completed pairs; admit larger strata only after a completed sentinel",
            "keep annotation and raw-read confidence columns separate from core admission",
        ],
        "measured_vs_unvalidated": {
            "measured": [
                "canonical BED/VCF/consensus/PSMCFA contents for P04 and P07",
                "all 402 primary/bootstrap PSMC outputs and all 1,152 sensitivity rows",
                "sacct elapsed allocation capacity, requested memory, state, node, and exit code",
                "four frozen terminal primary failures plus the post-snapshot P06 scheduler failure",
            ],
            "unvalidated": [
                "raw-read concordance and chemistry",
                "assembly QV, BUSCO, collapse/copy-number, k-mer, and standalone repeat audits",
                "species-specific mutation rates and generation times",
                "PSMC bootstrap uncertainty: outputs are finite but the resampled input population does not preserve primary masked bins",
                "population-level diversity or demographic generalization",
                "actual MaxRSS, TotalCPU, I/O, and scratch high-water from sacct",
                "biological outputs for P01/P02/P03/P05/P06/P08/P09/P10",
            ],
        },
    }


def parser() -> argparse.ArgumentParser:
    value = argparse.ArgumentParser(description=__doc__)
    value.add_argument("--canonical-root", type=Path, default=DEFAULT_VGP_ROOT)
    value.add_argument("--closed-world", type=Path, default=Path("analysis/vgp_real_pilot_closed_world_v1.json"))
    value.add_argument("--pilot-sacct", type=Path, default=Path("analysis/vgp_real_pilot_sacct_v1.tsv"))
    value.add_argument("--canary-sacct", type=Path, default=Path("analysis/vgp_real_canary_sacct_v1.tsv"))
    value.add_argument("--live-sacct", type=Path, default=Path("analysis/vgp_real_pilot_review_live_sacct_v1.tsv"))
    value.add_argument("--output", type=Path, required=True)
    return value


def main(argv: Sequence[str] | None = None) -> int:
    args = parser().parse_args(argv)
    review = build_review(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(review, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
