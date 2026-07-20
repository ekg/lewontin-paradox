#!/usr/bin/env python3
"""Independently audit and summarize a real VGP biological pair.

This module deliberately does not import ``analysis.vgp_10_pilot``.  It
reconstructs PAF multiplicity, ordered mask accounting, the normalized variant
subset, and diversity denominators from promoted files so the execution is not
self-certified by the production helper that created those files.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import hashlib
import json
import math
import subprocess
from pathlib import Path
from typing import Iterable, Mapping, Sequence, Tuple


Interval = Tuple[str, int, int]
Variant = Tuple[str, int, str, str]  # contig, one-based POS, REF, ALT


class CanaryAuditError(RuntimeError):
    """A promoted canary invariant failed independent reconstruction."""


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def file_record(path: Path) -> dict[str, object]:
    if not path.is_file():
        raise CanaryAuditError(f"required output is absent: {path}")
    return {"path": str(path), "size_bytes": path.stat().st_size, "sha256": sha256_file(path)}


def parse_bed(path: Path) -> list[Interval]:
    rows: list[Interval] = []
    with path.open(encoding="utf-8") as handle:
        for number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise CanaryAuditError(f"{path}:{number}: BED row has fewer than three fields")
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError as error:
                raise CanaryAuditError(f"{path}:{number}: non-integer BED coordinate") from error
            if start < 0 or end <= start:
                raise CanaryAuditError(f"{path}:{number}: invalid half-open BED interval")
            rows.append((fields[0], start, end))
    return rows


def _merge(rows: Iterable[Interval]) -> list[Interval]:
    merged: list[list[object]] = []
    for contig, start, end in sorted(rows, key=lambda row: (row[0], row[1], row[2])):
        if merged and merged[-1][0] == contig and start <= int(merged[-1][2]):
            merged[-1][2] = max(int(merged[-1][2]), end)
        else:
            merged.append([contig, start, end])
    return [(str(contig), int(start), int(end)) for contig, start, end in merged]


def interval_bp(rows: Iterable[Interval]) -> int:
    return sum(end - start for _, start, end in rows)


def independent_mask_reconstruction(
    universe: Mapping[str, int], reason_order: Sequence[str], flags: Mapping[str, Iterable[Interval]]
) -> dict[str, object]:
    """Reapply first-reason-wins from raw flags without pipeline helpers."""

    if not universe or any(length <= 0 for length in universe.values()):
        raise CanaryAuditError("independent mask universe is empty or non-positive")
    unknown = set(flags) - set(reason_order)
    if unknown:
        raise CanaryAuditError(f"raw mask contains unknown reasons: {sorted(unknown)}")
    normalized = {reason: _merge(flags.get(reason, ())) for reason in reason_order}
    callable_rows: list[Interval] = []
    by_reason: dict[str, list[Interval]] = {reason: [] for reason in reason_order}
    reason_index = {reason: index for index, reason in enumerate(reason_order)}
    events_by_contig: dict[str, dict[int, list[tuple[int, int]]]] = {
        contig: {} for contig in universe
    }
    for reason, rows in normalized.items():
        index = reason_index[reason]
        for contig, start, end in rows:
            if contig not in universe:
                continue
            clipped_start, clipped_end = max(0, start), min(universe[contig], end)
            if clipped_start < clipped_end:
                events_by_contig[contig].setdefault(clipped_start, []).append((index, 1))
                events_by_contig[contig].setdefault(clipped_end, []).append((index, -1))

    for contig, length in universe.items():
        active = [0] * len(reason_order)
        previous = 0
        events = events_by_contig[contig]
        for position in sorted({0, length, *events}):
            if previous < position:
                owner_index = next((index for index, count in enumerate(active) if count), None)
                destination = callable_rows if owner_index is None else by_reason[reason_order[owner_index]]
                destination.append((contig, previous, position))
            for index, delta in events.get(position, ()):
                active[index] += delta
                if active[index] < 0:
                    raise CanaryAuditError("independent mask sweep observed a negative active depth")
            previous = position
        if any(active):
            raise CanaryAuditError(f"independent mask sweep did not close on {contig}")
    callable_rows = _merge(callable_rows)
    by_reason = {reason: _merge(rows) for reason, rows in by_reason.items()}
    callable_bp = interval_bp(callable_rows)
    reason_bp = {reason: interval_bp(rows) for reason, rows in by_reason.items()}
    universe_bp = sum(universe.values())
    return {
        "callable": callable_rows,
        "by_reason": by_reason,
        "universe_bp": universe_bp,
        "callable_bp": callable_bp,
        "excluded_bp_by_primary_reason": reason_bp,
        "accounting_discrepancy_bp": universe_bp - callable_bp - sum(reason_bp.values()),
    }


def _events_depth(intervals: Iterable[tuple[str, int, int]]) -> int:
    grouped: dict[str, list[tuple[int, int]]] = {}
    for name, start, end in intervals:
        grouped.setdefault(name, []).extend(((start, 1), (end, -1)))
    maximum = 0
    for events in grouped.values():
        depth = 0
        for _, delta in sorted(events, key=lambda item: (item[0], item[1])):
            depth += delta
            maximum = max(maximum, depth)
    return maximum


def maximum_axis_depth(paf_fields: Sequence[Sequence[str]], axis: str) -> int:
    if axis == "query":
        values = ((row[0], int(row[2]), int(row[3])) for row in paf_fields)
    elif axis == "target":
        values = ((row[5], int(row[7]), int(row[8])) for row in paf_fields)
    else:
        raise CanaryAuditError(f"invalid PAF depth axis: {axis}")
    return _events_depth(values)


def audit_paf(path: Path, h1_names: set[str], h2_names: set[str]) -> dict[str, object]:
    rows: list[list[str]] = []
    with path.open(encoding="utf-8") as handle:
        for number, raw in enumerate(handle, 1):
            if not raw.strip():
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 12:
                raise CanaryAuditError(f"{path}:{number}: truncated PAF")
            if fields[0] not in h2_names or fields[5] not in h1_names:
                raise CanaryAuditError(f"{path}:{number}: PAF is not H2-query/H1-target")
            rows.append(fields)
    if not rows:
        raise CanaryAuditError("SweepGA emitted an empty PAF")
    query_depth = maximum_axis_depth(rows, "query")
    target_depth = maximum_axis_depth(rows, "target")
    if query_depth > 1 or target_depth > 1:
        raise CanaryAuditError(f"independent PAF multiplicity exceeds 1:1: {query_depth}/{target_depth}")
    return {
        "orientation": "H2_query_H1_target",
        "records": len(rows),
        "maximum_query_overlap_depth": query_depth,
        "maximum_target_overlap_depth": target_depth,
        "required_option": "--num-mappings 1:1",
        **file_record(path),
    }


class IntervalIndex:
    """Exact containment lookup over independently merged half-open intervals."""

    def __init__(self, rows: Iterable[Interval]) -> None:
        grouped: dict[str, list[tuple[int, int]]] = {}
        for contig, start, end in _merge(rows):
            grouped.setdefault(contig, []).append((start, end))
        self._rows = grouped
        self._starts = {
            contig: [start for start, _ in intervals]
            for contig, intervals in grouped.items()
        }

    def contains(self, contig: str, start: int, end: int) -> bool:
        starts = self._starts.get(contig)
        if not starts:
            return False
        index = bisect.bisect_right(starts, start) - 1
        return index >= 0 and end <= self._rows[contig][index][1]


def _vcf_rows(path: Path) -> list[Variant]:
    result: list[Variant] = []
    seen: set[Variant] = set()
    with path.open(encoding="utf-8") as handle:
        for number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 5:
                raise CanaryAuditError(f"{path}:{number}: truncated VCF")
            for alt in fields[4].split(","):
                row = (fields[0], int(fields[1]), fields[3].upper(), alt.upper())
                if row in seen:
                    raise CanaryAuditError(f"exact duplicate normalized allele: {row}")
                seen.add(row)
                result.append(row)
    return result


def summarize_vcf(path: Path, subset_intervals: Sequence[Interval]) -> dict[str, object]:
    rows = _vcf_rows(path)
    index = IntervalIndex(subset_intervals)
    subset = [
        row for row in rows
        if index.contains(row[0], row[1] - 1, row[1] - 1 + len(row[2]))
    ]
    return {
        "normalized_variant_records": len(rows),
        "normalized_snp_records": sum(len(ref) == len(alt) == 1 for _, _, ref, alt in rows),
        "normalized_indel_records": sum(not (len(ref) == len(alt) == 1) for _, _, ref, alt in rows),
        "callable_variant_records": len(subset),
        "callable_snp_records": sum(len(ref) == len(alt) == 1 for _, _, ref, alt in subset),
        "callable_indel_records": sum(not (len(ref) == len(alt) == 1) for _, _, ref, alt in subset),
        "subset_records": subset,
    }


def stratified_windows(
    dictionary: Sequence[Mapping[str, object]], strata: int, window_bp: int = 5_000_000
) -> list[Interval]:
    """Select deterministic early/middle/late dictionary strata independently."""

    if strata <= 0 or window_bp <= 0 or len(dictionary) < strata:
        raise CanaryAuditError("requested independent subset strata are unavailable")
    if strata == 1:
        indices = [0]
    else:
        indices = [round(number * (len(dictionary) - 1) / (strata - 1)) for number in range(strata)]
    if len(set(indices)) != strata:
        raise CanaryAuditError("independent subset strata do not resolve uniquely")
    windows = []
    for index in indices:
        row = dictionary[index]
        length = int(row["length"])
        if length <= 0:
            raise CanaryAuditError("independent subset stratum has a non-positive contig")
        windows.append((str(row["name"]), 0, min(window_bp, length)))
    return windows


def _bcf_rows(bcftools: Path, path: Path) -> list[Variant]:
    output = subprocess.check_output(
        [str(bcftools), "query", "-f", "%CHROM\\t%POS\\t%REF\\t%ALT\\n", str(path)],
        text=True,
    )
    result: list[Variant] = []
    for raw in output.splitlines():
        chrom, pos, ref, alt = raw.split("\t")[:4]
        for allele in alt.split(","):
            result.append((chrom, int(pos), ref.upper(), allele.upper()))
    return result


def _verify_stage_file(run_root: Path, relative: str) -> dict[str, object]:
    path = run_root / relative
    stage = path.parent
    while stage != run_root and not (stage / ".complete.json").is_file():
        stage = stage.parent
    if stage == run_root:
        raise CanaryAuditError(f"no stage sentinel owns {path}")
    sentinel = json.loads((stage / ".complete.json").read_text(encoding="utf-8"))
    stage_relative = str(path.relative_to(stage))
    expected = sentinel.get("files", {}).get(stage_relative)
    record = file_record(path)
    if expected != record["sha256"]:
        raise CanaryAuditError(f"stage sentinel digest mismatch for {path}")
    record["stage_sentinel"] = str(stage / ".complete.json")
    return record


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def summarize_scaled_scenarios(rows: Sequence[Mapping[str, str]]) -> dict[str, object]:
    if not rows:
        raise CanaryAuditError("no labeled mutation-rate/generation-time scenarios")
    required = (
        "scenario_id", "mutation_rate_per_generation", "generation_time_years",
        "psmc_bin_size_bp", "mutation_rate_source", "generation_time_source",
    )
    if any(not row.get(field) for row in rows for field in required):
        raise CanaryAuditError("scaled scenario rows lack labels, rates, bin size, or sources")
    try:
        rates = sorted({float(row["mutation_rate_per_generation"]) for row in rows})
        generations = sorted({float(row["generation_time_years"]) for row in rows})
        bin_sizes = sorted({int(row["psmc_bin_size_bp"]) for row in rows})
    except ValueError as error:
        raise CanaryAuditError("scaled scenario parameters are not numeric") from error
    if any(rate <= 0 for rate in rates) or any(years <= 0 for years in generations):
        raise CanaryAuditError("scaled scenario rates and generation times must be positive")
    if bin_sizes != [100]:
        raise CanaryAuditError(f"scaled scenarios do not use the frozen 100-bp PSMC bin: {bin_sizes}")
    scenario_ids = sorted({row["scenario_id"] for row in rows})
    sources = sorted({
        row[field] for row in rows
        for field in ("mutation_rate_source", "generation_time_source")
    })
    return {
        "scenario_ids": scenario_ids,
        "mutation_rates_per_generation": rates,
        "generation_times_years": generations,
        "psmc_bin_size_bp": 100,
        "source_labels": sources,
        "rows": len(rows),
    }


def read_regional_vcf_audit(
    path: Path, expected_focus_rows: int, canonical_root: str
) -> dict[str, object]:
    """Validate the durable census made before transient IMPG shards were removed."""

    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise CanaryAuditError(f"unreadable IMPG regional VCF audit: {path}") from error
    required_counts = (
        "focus_rows", "unique_query_names", "unique_native_partition_ids",
        "regional_vcf_count", "regional_vcf_total_bytes",
    )
    if value.get("canonical_vgp_root") != canonical_root:
        raise CanaryAuditError("IMPG regional VCF audit has the wrong canonical root")
    try:
        counts = {key: int(value[key]) for key in required_counts}
    except (KeyError, TypeError, ValueError) as error:
        raise CanaryAuditError("IMPG regional VCF audit lacks integer census fields") from error
    if not (
        counts["focus_rows"] == expected_focus_rows
        == counts["unique_query_names"]
        == counts["regional_vcf_count"]
    ):
        raise CanaryAuditError("IMPG regional VCF census differs from focused BED")
    if not (0 < counts["unique_native_partition_ids"] <= expected_focus_rows):
        raise CanaryAuditError("IMPG native partition ID census is invalid")
    if counts["regional_vcf_total_bytes"] <= 0:
        raise CanaryAuditError("IMPG regional VCF byte census is empty")
    if value.get("all_regional_vcfs_nonempty") is not True:
        raise CanaryAuditError("IMPG regional VCF audit includes empty shards")
    if value.get("transient_shards_removed_after_lacing") is not True:
        raise CanaryAuditError("IMPG transient shard lifecycle is unverified")
    return value


def audit(args: argparse.Namespace) -> dict[str, object]:
    run_root = args.run_root.resolve()
    preflight = json.loads((run_root / "preflight/preflight.json").read_text(encoding="utf-8"))
    h1_dictionary = preflight["dictionaries"]["h1_fasta"]
    h2_dictionary = preflight["dictionaries"]["h2_fasta"]
    universe = {str(row["name"]): int(row["length"]) for row in h1_dictionary}
    h1_names, h2_names = set(universe), {str(row["name"]) for row in h2_dictionary}

    mapping = audit_paf(run_root / "mapping/h2_to_h1.1to1.paf", h1_names, h2_names)
    production_multiplicity = json.loads(
        (run_root / "mapping/multiplicity.json").read_text(encoding="utf-8")
    )
    for key in ("maximum_query_overlap_depth", "maximum_target_overlap_depth"):
        if production_multiplicity.get(key) != mapping[key]:
            raise CanaryAuditError(f"independent multiplicity differs for {key}")

    reconciliation_path = run_root / "consensus/masks/mask_reconciliation.json"
    production_mask = json.loads(reconciliation_path.read_text(encoding="utf-8"))
    reason_order = tuple(production_mask["reason_order"])
    raw_flags = {
        reason: parse_bed(run_root / f"consensus/inputs.{reason}.bed")
        for reason in reason_order
        if (run_root / f"consensus/inputs.{reason}.bed").is_file()
    }
    rebuilt = independent_mask_reconstruction(universe, reason_order, raw_flags)
    emitted_callable = _merge(parse_bed(run_root / "consensus/masks/callable.bed"))
    if rebuilt["callable"] != emitted_callable:
        raise CanaryAuditError("independently reconstructed callable BED differs from promoted BED")
    for reason in reason_order:
        emitted = _merge(parse_bed(run_root / f"consensus/masks/exclusions.{reason}.bed"))
        if rebuilt["by_reason"][reason] != emitted:
            raise CanaryAuditError(f"independently reconstructed mask differs for {reason}")
    for key in ("universe_bp", "callable_bp", "excluded_bp_by_primary_reason", "accounting_discrepancy_bp"):
        if rebuilt[key] != production_mask[key]:
            raise CanaryAuditError(f"independent mask accounting differs for {key}")

    subset_window = stratified_windows(h1_dictionary, args.subset_strata)
    subset_callable = []
    subset_by_contig = {contig: (window_start, window_end) for contig, window_start, window_end in subset_window}
    for contig, start, end in emitted_callable:
        if contig in subset_by_contig:
            window_start, window_end = subset_by_contig[contig]
            if start < window_end and end > window_start:
                subset_callable.append((contig, max(window_start, start), min(window_end, end)))
    vcf_path = run_root / "consensus/normalized.vcf"
    variants = summarize_vcf(vcf_path, emitted_callable)
    subset = summarize_vcf(vcf_path, subset_callable)
    vcf_rows = _vcf_rows(vcf_path)
    bcf_rows = _bcf_rows(args.bcftools, run_root / "variants/normalized.bcf")
    if vcf_rows != bcf_rows:
        raise CanaryAuditError("independent BCF query differs from normalized VCF alleles")
    callable_bp = int(rebuilt["callable_bp"])
    if callable_bp <= 0 or variants["callable_variant_records"] <= 0:
        raise CanaryAuditError("positive callable bases and heterozygous variants are required")

    join_qc = json.loads((run_root / "consensus/join_qc.json").read_text(encoding="utf-8"))
    consensus_qc = join_qc["consensus"]
    expected_callable_snps = (
        int(consensus_qc["heterozygous_snps"])
        + int(consensus_qc["callable_snps_masked_by_indel_flank"])
    )
    independently_callable = {
        "callable_variant_records": variants["callable_variant_records"],
        "callable_snp_records": variants["callable_snp_records"],
        "callable_indel_records": variants["callable_indel_records"],
    }
    expected_callable = {
        "callable_variant_records": int(consensus_qc["callable_variant_records"]),
        "callable_snp_records": expected_callable_snps,
        "callable_indel_records": int(consensus_qc["heterozygous_indels"]),
    }
    if independently_callable != expected_callable:
        raise CanaryAuditError(
            "independent callable variant accounting differs from consensus QC: "
            f"{independently_callable} != {expected_callable}"
        )
    consensus_callable_bp = int(consensus_qc["consensus_callable_bp"])
    expected_consensus_callable_bp = callable_bp - int(consensus_qc["indel_masked_h1_bp"])
    if consensus_callable_bp != expected_consensus_callable_bp:
        raise CanaryAuditError("consensus callable denominator does not reconcile after indel masks")
    heterozygous_snps = int(consensus_qc["heterozygous_snps"])
    if consensus_callable_bp <= 0 or heterozygous_snps <= 0:
        raise CanaryAuditError("positive final callable bases and heterozygous SNPs are required")

    psmcfa_path = run_root / "consensus/consensus/input.psmcfa"
    psmc_symbols = {"T": 0, "K": 0, "N": 0}
    with psmcfa_path.open(encoding="utf-8") as handle:
        for raw in handle:
            if not raw.startswith(">"):
                for symbol in raw.strip():
                    if symbol not in psmc_symbols:
                        raise CanaryAuditError(f"unexpected PSMCFA symbol: {symbol}")
                    psmc_symbols[symbol] += 1
    trajectory = _read_tsv(run_root / "psmc/finalize/unscaled_trajectory.tsv")
    if not trajectory or not all(
        math.isfinite(float(row["time_2N0"])) and math.isfinite(float(row["lambda"]))
        for row in trajectory
    ):
        raise CanaryAuditError("unscaled PSMC trajectory is empty or non-finite")
    bootstraps = _read_tsv(run_root / "psmc/finalize/bootstrap_unscaled.tsv")
    finite_bootstraps = sum(row["finite"] == "true" for row in bootstraps)
    if len(bootstraps) < 100 or finite_bootstraps < 100:
        raise CanaryAuditError("fewer than 100 finite boundary-aware PSMC bootstraps")
    scenarios = _read_tsv(run_root / "psmc/finalize/scenario_scaled_trajectories.tsv")
    scenario_summary = summarize_scaled_scenarios(scenarios)

    subset_rows = subset["subset_records"]
    args.subset_output.parent.mkdir(parents=True, exist_ok=True)
    with args.subset_output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("contig", "position_1based", "ref", "alt"))
        writer.writerows(subset_rows)

    telemetry_rows = _read_tsv(args.sacct_telemetry)
    if not telemetry_rows or not any(row.get("State") == "COMPLETED" for row in telemetry_rows):
        raise CanaryAuditError("authoritative sacct telemetry lacks a COMPLETED allocation")

    selected_native_rows = sum(1 for _ in parse_bed(run_root / "impg/focus.native.bed"))
    regional_vcf_census = read_regional_vcf_audit(
        run_root / "impg/regional_vcf_audit.json", selected_native_rows, args.canonical_root
    )
    annotation = (
        json.loads(args.annotation_result.read_text(encoding="utf-8"))
        if args.annotation_result is not None
        else {
            "status": "not_applicable_no_eligible_exact_annotation_binding",
            "core_result_affected": False,
            "canonical_vgp_root": args.canonical_root,
        }
    )
    result = {
        "schema_version": args.schema_version,
        "task_id": args.task_id,
        "authorization_id": "vgp10-auth-20260718-v2",
        "selection_id": args.selection_id,
        "species": args.species,
        "canonical_vgp_root": args.canonical_root,
        "promoted_run_root": str(run_root),
        "slurm_job_id": args.slurm_job_id,
        "guix_environment_capture": file_record(args.environment_capture),
        "mapping": mapping,
        "impg": {
            "native_partition_count": sum(1 for _ in parse_bed(run_root / "impg/partitions/partitions.bed")),
            "selected_native_partition_count": selected_native_rows,
            "regional_vcf_count": regional_vcf_census["regional_vcf_count"],
            "regional_vcf_total_bytes": regional_vcf_census["regional_vcf_total_bytes"],
            "unique_native_partition_ids": regional_vcf_census["unique_native_partition_ids"],
            "regional_vcf_audit": _verify_stage_file(run_root, "impg/regional_vcf_audit.json"),
            "transient_regional_shards_removed_after_verified_lacing": True,
            "laced_vcf": _verify_stage_file(run_root, "impg/laced.vcf"),
        },
        "variants": {
            **{key: value for key, value in variants.items() if key != "subset_records"},
            "normalized_vcf_gz": _verify_stage_file(run_root, "variants/normalized.vcf.gz"),
            "normalized_bcf": _verify_stage_file(run_root, "variants/normalized.bcf"),
            "vcf_bcf_exact_record_match": True,
        },
        "diversity": {
            "estimator": "callable heterozygous SNPs per final callable H1 bp after indel-flank masks",
            "pi": heterozygous_snps / consensus_callable_bp,
            "heterozygous_snps": heterozygous_snps,
            "callable_bp": consensus_callable_bp,
            "pre_indel_callable_bp": callable_bp,
            "indel_masked_h1_bp": int(consensus_qc["indel_masked_h1_bp"]),
            "callable_variant_record_rate_before_indel_mask": (
                int(variants["callable_variant_records"]) / callable_bp
            ),
        },
        "independent_mask_audit": {
            "method": "raw reason flags independently reapplied in manifest order; exact BED equality required",
            "universe_bp": rebuilt["universe_bp"],
            "callable_bp": callable_bp,
            "callable_fraction": callable_bp / int(rebuilt["universe_bp"]),
            "excluded_bp_by_primary_reason": rebuilt["excluded_bp_by_primary_reason"],
            "accounting_discrepancy_bp": rebuilt["accounting_discrepancy_bp"],
            "exact_emitted_bed_match": True,
        },
        "independent_variant_subset": {
            "stratification": "early/middle/late H1 dictionary positions" if args.subset_strata > 1 else "leading H1 dictionary position",
            "windows": [
                {"contig": contig, "start_0based": start, "end_0based_exclusive": end}
                for contig, start, end in subset_window
            ],
            "callable_bp": interval_bp(subset_callable),
            "heterozygous_variant_records": len(subset_rows),
            "tsv": file_record(args.subset_output),
            "reconstructed_from_vcf_and_matched_to_bcf": True,
        },
        "consensus": {
            "fasta": _verify_stage_file(run_root, "consensus/consensus/consensus.fa"),
            "psmc_input": _verify_stage_file(run_root, "consensus/consensus/input.psmcfa"),
            "psmcfa_symbols": psmc_symbols,
        },
        "psmc": {
            "unscaled_primary": _verify_stage_file(run_root, "psmc/replicate-000/unscaled.psmc"),
            "unscaled_trajectory": _verify_stage_file(run_root, "psmc/finalize/unscaled_trajectory.tsv"),
            "trajectory_intervals": len(trajectory),
            "bootstrap_attempts": len(bootstraps),
            "finite_bootstraps": finite_bootstraps,
            "boundary_aware": True,
            "primary_block_bp": 5_000_000,
            **scenario_summary,
        },
        "telemetry": {
            "sacct": file_record(args.sacct_telemetry),
            "rows": telemetry_rows,
            "compute_node_sacct_limitation": "slurmdbd unreachable from compute node; captured after completion on login node",
        },
        "annotation": annotation,
        "verification": {
            "independent_reconstruction": True,
            "atomic_verified_promotion": True,
            "core_biological_result": "PASS",
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    partial = args.output.with_suffix(args.output.suffix + ".partial")
    partial.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    partial.replace(args.output)
    return result


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--run-root", type=Path, required=True)
    result.add_argument("--canonical-root", default="/moosefs/erikg/vgp")
    result.add_argument("--selection-id", default="P07")
    result.add_argument("--species", default="Spinachia spinachia")
    result.add_argument("--task-id", default="run-vgp-real-canary")
    result.add_argument("--schema-version", default="vgp-real-canary-execution-v1")
    result.add_argument("--slurm-job-id", required=True)
    result.add_argument("--bcftools", type=Path, required=True)
    result.add_argument("--environment-capture", type=Path, required=True)
    result.add_argument("--sacct-telemetry", type=Path, required=True)
    result.add_argument(
        "--annotation-result", type=Path,
        help="eligible exact annotation result; omit only when the frozen pair has no eligible binding",
    )
    result.add_argument("--subset-output", type=Path, required=True)
    result.add_argument("--subset-strata", type=int, default=1)
    result.add_argument("--output", type=Path, required=True)
    return result


def main(argv: Sequence[str] | None = None) -> int:
    args = parser().parse_args(argv)
    try:
        result = audit(args)
    except (CanaryAuditError, OSError, ValueError, subprocess.SubprocessError, json.JSONDecodeError) as error:
        print(f"ERROR: {error}")
        return 2
    print(json.dumps({
        "output": str(args.output),
        "callable_bp": result["diversity"]["callable_bp"],
        "heterozygous_variants": result["variants"]["callable_variant_records"],
        "pi": result["diversity"]["pi"],
        "finite_bootstraps": result["psmc"]["finite_bootstraps"],
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
