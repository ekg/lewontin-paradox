#!/usr/bin/env python3
"""Independent completion audit for the three bounded-range VGP results."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import subprocess
import tempfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable, Mapping, Sequence

from analysis.audit_vgp_three_pair_results import (
    fasta_lengths,
    fasta_non_n_counts,
    psmc_population_counts,
)
from analysis.audit_vgp_real_canary import independent_mask_reconstruction
from analysis.vgp_10_pilot import (
    Interval,
    PilotError,
    canonical_json,
    interval_bp,
    read_bed,
    sha256_file,
)
from analysis.vgp_three_pair import audit_graph_identifiers


RUN_ID = "vgp-three-pair-20260722-v1"
BOUNDED_BASE = Path("/moosefs/erikg/vgp/pilot/three-pair") / RUN_ID
MAPPING_ROOTS = {
    "P07": Path(
        "/moosefs/erikg/vgp/pilot/clean-canary/"
        "vgp-clean-canary-20260722-v1/P07/mapping"
    ),
    "P02": Path("/moosefs/erikg/vgp/pilot/runs") / RUN_ID / "P02/mapping",
    "P03": Path("/moosefs/erikg/vgp/pilot/runs") / RUN_ID / "P03/mapping",
}
H1_FASTAS = {
    "P07": Path(
        "/moosefs/erikg/vgp/derived/clean-canary-bgzf/"
        "GCA_048126635.1.fa.gz"
    ),
    "P02": Path("/moosefs/erikg/vgp/pilot/inputs/P02/h1.fa"),
    "P03": Path("/moosefs/erikg/vgp/pilot/inputs/P03/h1.fa"),
}
P07_GRAPH_LEDGER = Path(
    "/moosefs/erikg/vgp/pilot/clean-canary/"
    "vgp-clean-canary-20260722-v1/P07/vgp_clean_canary_graph_sequence_digests.tsv"
)


def audit_closed(root: Path) -> dict[str, object]:
    sentinel = json.loads((root / ".complete.json").read_text())
    files = sentinel.get("files")
    if not isinstance(files, Mapping) or not files:
        raise PilotError(f"bounded result has no closed digest ledger: {root}")
    for relative, digest in files.items():
        path = root / relative
        if not path.is_file() or sha256_file(path) != digest:
            raise PilotError(f"bounded digest mismatch: {path}")
    return {
        "sentinel_sha256": sha256_file(root / ".complete.json"),
        "verified_files": len(files),
        "verified_bytes": sum((root / relative).stat().st_size for relative in files),
    }


def _range_map(plan: Mapping[str, object]) -> dict[str, Mapping[str, object]]:
    rows = {str(row["range_id"]): row for row in plan["ranges"]}  # type: ignore[index]
    if len(rows) != int(plan["range_count"]):
        raise PilotError("bounded plan repeats a range identifier")
    return rows


def _merge_rows(rows: Iterable[Interval]) -> list[Interval]:
    merged: list[Interval] = []
    for row in sorted(rows, key=lambda value: (value.contig, value.start, value.end)):
        if row.start < 0 or row.end <= row.start:
            raise PilotError(f"invalid independently audited interval: {row}")
        if merged and merged[-1].contig == row.contig and row.start <= merged[-1].end:
            prior = merged[-1]
            merged[-1] = Interval(prior.contig, prior.start, max(prior.end, row.end))
        else:
            merged.append(row)
    return merged


def _intersection_bp(left: Iterable[Interval], right: Iterable[Interval]) -> int:
    left_by_contig: dict[str, list[Interval]] = defaultdict(list)
    right_by_contig: dict[str, list[Interval]] = defaultdict(list)
    for row in _merge_rows(left):
        left_by_contig[row.contig].append(row)
    for row in _merge_rows(right):
        right_by_contig[row.contig].append(row)
    total = 0
    for contig, left_rows in left_by_contig.items():
        right_rows = right_by_contig.get(contig, ())
        i = j = 0
        while i < len(left_rows) and j < len(right_rows):
            row, other = left_rows[i], right_rows[j]
            total += max(0, min(row.end, other.end) - max(row.start, other.start))
            if row.end <= other.end:
                i += 1
            else:
                j += 1
    return total


def audit_range_plan(root: Path, plan: Mapping[str, object]) -> dict[str, object]:
    """Recompute exhaustive ranges and native-partition ownership from promoted files."""

    ranges = list(plan["ranges"])  # type: ignore[index]
    by_contig: dict[str, list[Mapping[str, object]]] = defaultdict(list)
    for row in ranges:
        if int(row["length"]) != int(row["end"]) - int(row["start"]):
            raise PilotError("range plan contains an inconsistent length")
        by_contig[str(row["contig"])].append(row)
    h1_lengths: dict[str, int] = {}
    total = 0
    for contig, rows in by_contig.items():
        rows.sort(key=lambda row: int(row["start"]))
        cursor = 0
        for row in rows:
            if int(row["start"]) != cursor:
                raise PilotError(f"range plan is not exhaustive/disjoint on {contig}")
            cursor = int(row["end"])
            total += int(row["length"])
        h1_lengths[contig] = cursor
    if (
        total != int(plan["h1_total_bp"])
        or len(h1_lengths) != int(plan["h1_contig_count"])
        or len(ranges) != int(plan["range_count"])
    ):
        raise PilotError("independent range plan census differs from frozen totals")

    owners = Counter()
    partition_rows = 0
    partition_bp = 0
    ranges_by_contig = {
        contig: sorted(rows, key=lambda row: int(row["start"]))
        for contig, rows in by_contig.items()
    }
    for number, line in enumerate((root / "index/partitions.bed").open(), 1):
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 4:
            raise PilotError(f"partition row {number} is truncated")
        contig, start, end = fields[0], int(fields[1]), int(fields[2])
        if contig not in h1_lengths:
            continue
        matches = [
            row for row in ranges_by_contig[contig]
            if int(row["start"]) <= start and end <= int(row["end"])
        ]
        if len(matches) != 1:
            raise PilotError(f"H1 partition row {number} does not have exactly one range owner")
        owners[str(matches[0]["range_id"])] += 1
        partition_rows += 1
        partition_bp += end - start
    for row in ranges:
        observed = owners[str(row["range_id"])]
        expected = int(row["partition_count"])
        if observed != expected:
            raise PilotError(
                f"{row['range_id']}: independently owned partitions {observed} != {expected}"
            )
    if (
        partition_rows != int(plan["native_h1_partition_rows"])
        or partition_bp != int(plan["partitioned_h1_bp"])
    ):
        raise PilotError("independent native-partition census differs from frozen plan")
    return {
        "range_count": len(ranges),
        "h1_contig_count": len(h1_lengths),
        "h1_total_bp": total,
        "native_h1_partition_rows": partition_rows,
        "partitioned_h1_bp": partition_bp,
        "range_plan_disjoint": True,
        "range_plan_exhaustive": True,
        "native_partition_one_owner": True,
        "h1_lengths": h1_lengths,
    }


def audit_range_variants(root: Path, plan: Mapping[str, object]) -> dict[str, object]:
    ranges = _range_map(plan)
    seen: set[tuple[str, int, str, str]] = set()
    records = 0
    boundary_failures = 0
    for range_id, row in ranges.items():
        if not row["query_required"]:
            continue
        path = root / f"ranges/{range_id}/normalized.vcf.gz"
        with gzip.open(path, "rt") as handle:
            for number, line in enumerate(handle, 1):
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 5:
                    raise PilotError(f"{path}:{number}: malformed normalized range VCF")
                key = (fields[0], int(fields[1]) - 1, fields[3], fields[4])
                if (
                    key[0] != row["contig"]
                    or not int(row["start"]) <= key[1] < int(row["end"])
                ):
                    boundary_failures += 1
                if key in seen:
                    raise PilotError(f"exact normalized variant has two range owners: {key}")
                seen.add(key)
                records += 1
    completion = json.loads((root / "range_completion.json").read_text())
    if records != int(completion["normalized_variant_records"]) or boundary_failures:
        raise PilotError("range variant reduction or half-open ownership failed")
    return {
        "normalized_variant_records": records,
        "unique_normalized_variant_keys": len(seen),
        "boundary_ownership_failures": boundary_failures,
        "exact_duplicate_keys_between_ranges": 0,
    }


def audit_callable_ownership(root: Path, plan: Mapping[str, object]) -> dict[str, object]:
    by_contig: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
    for row in plan["ranges"]:  # type: ignore[index]
        by_contig[str(row["contig"])].append(
            (int(row["start"]), int(row["end"]), str(row["range_id"]))
        )
    for values in by_contig.values():
        values.sort()
    per_range = Counter()
    total = 0
    for interval in read_bed(root / "consensus/masks/callable.bed"):
        remaining = interval.length
        owners = 0
        for start, end, range_id in by_contig.get(interval.contig, []):
            if start >= interval.end:
                break
            overlap = max(0, min(end, interval.end) - max(start, interval.start))
            if overlap:
                per_range[range_id] += overlap
                remaining -= overlap
                owners += 1
        if remaining:
            raise PilotError(f"callable bases lack a range owner: {interval}")
        total += interval.length
    execution = json.loads((root / "execution.json").read_text())
    expected = int(execution["diversity"]["callable_bp"])
    if total != expected or sum(per_range.values()) != expected:
        raise PilotError("callable range reduction differs from biological denominator")
    return {
        "callable_bp": total,
        "summed_per_range_callable_bp": sum(per_range.values()),
        "unowned_callable_bp": 0,
        "multiply_owned_callable_bp": 0,
        "callable_one_owner_accounting": True,
    }


def audit_mask_accounting(
    root: Path, h1_lengths: Mapping[str, int],
) -> dict[str, object]:
    """Rebuild ordered prevariant masks and reconcile the final indel-flank mask."""

    production = json.loads(
        (root / "consensus/masks/mask_reconciliation.json").read_text()
    )
    reason_order = [
        str(reason) for reason in production["reason_order"]
        if reason != "variant_indel_flank"
    ]
    flags: dict[str, list[tuple[str, int, int]]] = {}
    for reason in reason_order:
        path = root / f"consensus/inputs.{reason}.bed"
        if path.is_file():
            flags[reason] = [
                (row.contig, row.start, row.end) for row in read_bed(path)
            ]
    try:
        rebuilt = independent_mask_reconstruction(h1_lengths, reason_order, flags)
    except RuntimeError as error:
        raise PilotError(f"independent mask reconstruction failed: {error}") from error
    prevariant = _merge_rows(read_bed(root / "consensus/masks/callable.prevariant.bed"))
    rebuilt_prevariant = [
        Interval(contig, start, end) for contig, start, end in rebuilt["callable"]
    ]
    if prevariant != rebuilt_prevariant:
        raise PilotError("independently reconstructed prevariant callable BED differs")
    for reason in reason_order:
        emitted = _merge_rows(
            read_bed(root / f"consensus/masks/exclusions.{reason}.bed")
        )
        independent = [
            Interval(contig, start, end)
            for contig, start, end in rebuilt["by_reason"][reason]
        ]
        if emitted != independent:
            raise PilotError(f"independent primary mask differs for {reason}")
        if (
            int(production["excluded_bp_by_primary_reason"][reason])
            != int(rebuilt["excluded_bp_by_primary_reason"][reason])
        ):
            raise PilotError(f"independent primary mask count differs for {reason}")
    final = _merge_rows(read_bed(root / "consensus/masks/callable.bed"))
    final_bp = interval_bp(final)
    if _intersection_bp(final, prevariant) != final_bp:
        raise PilotError("final callable mask is not a subset of prevariant callable")
    indel_bp = interval_bp(prevariant) - final_bp
    if (
        int(rebuilt["universe_bp"]) != int(production["universe_bp"])
        or int(rebuilt["callable_bp"]) != int(production["prevariant_callable_bp"])
        or final_bp != int(production["callable_bp"])
        or indel_bp != int(production["variant_indel_flank_excluded_bp"])
        or indel_bp
        != int(production["excluded_bp_by_primary_reason"]["variant_indel_flank"])
        or int(production["accounting_discrepancy_bp"]) != 0
    ):
        raise PilotError("independent final mask accounting differs from production")
    return {
        "universe_bp": int(rebuilt["universe_bp"]),
        "prevariant_callable_bp": int(rebuilt["callable_bp"]),
        "final_callable_bp": final_bp,
        "variant_indel_flank_excluded_bp": indel_bp,
        "excluded_bp_by_primary_reason": production["excluded_bp_by_primary_reason"],
        "accounting_discrepancy_bp": 0,
        "final_callable_subset_of_prevariant": True,
        "independently_reconstructed": True,
    }


def _query_keys(bcftools: str, path: Path) -> list[str]:
    result = subprocess.run(
        [bcftools, "query", "-f", "%CHROM\\t%POS\\t%REF\\t%ALT\\n", str(path)],
        check=True, capture_output=True, text=True,
    )
    return result.stdout.splitlines()


def _selected_fasta_sequences(path: Path, names: set[str]) -> dict[str, str]:
    opener = gzip.open if path.suffix == ".gz" else open
    parts: dict[str, list[str]] = {}
    current: str | None = None
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                name = line[1:].split(None, 1)[0]
                current = name if name in names else None
                if current is not None:
                    if current in parts:
                        raise PilotError(f"FASTA repeats selected identifier: {current}")
                    parts[current] = []
            elif current is not None:
                parts[current].append(line.strip())
    missing = names - set(parts)
    if missing:
        raise PilotError(f"selected staged-H1 FASTA identifiers are absent: {sorted(missing)}")
    return {name: "".join(value).upper() for name, value in parts.items()}


def stratified_reaudit(
    selection_id: str, root: Path, plan: Mapping[str, object], bcftools: str,
) -> list[dict[str, object]]:
    query = [row for row in plan["ranges"] if row["query_required"]]  # type: ignore[index]
    if len(query) < 3:
        raise PilotError(f"{selection_id}: fewer than three bounded validation ranges")
    chosen = [query[0], query[len(query) // 2], query[-1]]
    regions = [
        Interval(str(row["contig"]), int(row["start"]), int(row["end"])) for row in chosen
    ]
    h1_sequences = _selected_fasta_sequences(
        H1_FASTAS[selection_id], {region.contig for region in regions}
    )
    consensus = root / "consensus/consensus/consensus.fa"
    non_n = fasta_non_n_counts(consensus, regions)
    population = psmc_population_counts(root / "consensus/consensus/input.psmcfa", regions)
    callable_rows = read_bed(root / "consensus/masks/callable.bed")
    result = []
    for label, row, region in zip(("early", "middle", "late"), chosen, regions):
        range_root = root / f"ranges/{row['range_id']}"
        vcf = range_root / "normalized.vcf.gz"
        bcf = range_root / "normalized.bcf"
        vcf_keys, bcf_keys = _query_keys(bcftools, vcf), _query_keys(bcftools, bcf)
        if vcf_keys != bcf_keys:
            raise PilotError(f"{selection_id}:{label}: range VCF and BCF differ")
        reference = h1_sequences[region.contig]
        ref_mismatches = 0
        for key in vcf_keys:
            chrom, pos, ref, _ = key.split("\t")
            pos0 = int(pos) - 1
            if (
                chrom != region.contig
                or not region.start <= pos0 < region.end
                or reference[pos0:pos0 + len(ref)].upper() != ref.upper()
            ):
                ref_mismatches += 1
        if ref_mismatches:
            raise PilotError(
                f"{selection_id}:{label}: staged-H1 REF reconstruction mismatches"
            )
        callable_bp = sum(
            max(0, min(item.end, region.end) - max(item.start, region.start))
            for item in callable_rows if item.contig == region.contig
        )
        pop = population[region]
        if set(pop) - set("NKT"):
            raise PilotError(f"{selection_id}:{label}: invalid PSMC symbol")
        result.append({
            "stratum": label, "range_id": row["range_id"], "contig": region.contig,
            "start": region.start, "end": region.end, "variant_records": len(vcf_keys),
            "vcf_bcf_exact_keys_equal": True, "callable_bp": callable_bp,
            "staged_h1_ref_alleles_revalidated": len(vcf_keys),
            "staged_h1_ref_mismatches": 0,
            "consensus_non_N_bp": non_n[region],
            "psmc_population_bins": {symbol: pop[symbol] for symbol in "NKT"},
        })
    return result


def _axis_overlap_depth(rows: Iterable[tuple[str, int, int]]) -> int:
    events: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for name, start, end in rows:
        events[name].extend(((start, 1), (end, -1)))
    maximum = 0
    for values in events.values():
        depth = 0
        for _, delta in sorted(values, key=lambda value: (value[0], value[1])):
            depth += delta
            maximum = max(maximum, depth)
    return maximum


def _graph_lengths(
    selection_id: str, root: Path,
) -> tuple[dict[str, int], dict[str, int]]:
    if selection_id == "P07":
        sides: dict[str, dict[str, int]] = {"H1": {}, "H2": {}}
        with P07_GRAPH_LEDGER.open(newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if (
                    row["side"] not in sides
                    or len(row["sequence_sha256"]) != 64
                    or row["sequence_id"] in sides[row["side"]]
                ):
                    raise PilotError("P07 graph sequence digest ledger is malformed")
                sides[row["side"]][row["sequence_id"]] = int(row["length"])
        return sides["H1"], sides["H2"]
    dictionary = json.loads(
        (root / "index/staged_fasta_dictionary.json").read_text()
    )
    roles = dictionary["roles"]
    return (
        {
            str(row["name"]): int(row["length"])
            for row in roles["h1_fasta"]["records"]
        },
        {
            str(row["name"]): int(row["length"])
            for row in roles["h2_fasta"]["records"]
        },
    )


def audit_graph_ids_independent(
    selection_id: str, root: Path, production: Mapping[str, object],
) -> dict[str, object]:
    paf = MAPPING_ROOTS[selection_id] / "h2_to_h1.1to1.paf"
    partitions = root / "index/partitions.bed"
    if selection_id != "P07":
        with tempfile.TemporaryDirectory(prefix=f"vgp-{selection_id}-graph-audit-") as tmp:
            output = Path(tmp) / "audit.json"
            rebuilt = audit_graph_identifiers(
                root / "index/staged_fasta_dictionary.json", paf, partitions, output
            )
        for key in (
            "paf_rows", "partition_rows", "partition_contigs",
            "partition_rows_by_staged_role", "staged_h1_records",
            "staged_h2_records", "digest_validated_aliases_used",
        ):
            if rebuilt[key] != production[key]:
                raise PilotError(
                    f"{selection_id}: independent graph identifier audit differs for {key}"
                )
        return {**rebuilt, "independently_recomputed": True}

    h1, h2 = _graph_lengths(selection_id, root)
    used: set[str] = set()
    paf_rows = 0
    with paf.open() as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if (
                len(fields) < 12
                or fields[0] not in h2
                or fields[5] not in h1
                or int(fields[1]) != h2.get(fields[0])
                or int(fields[6]) != h1.get(fields[5])
            ):
                raise PilotError(f"P07: graph PAF identifier failure at row {number}")
            used.update((fields[0], fields[5]))
            paf_rows += 1
    partition_rows = 0
    for number, line in enumerate(partitions.open(), 1):
        if not line.strip():
            continue
        fields = line.rstrip("\n").split("\t")
        names = [values for values in (h1, h2) if fields[0] in values]
        if len(fields) < 4 or not names:
            raise PilotError(f"P07: unresolved partition identifier at row {number}")
        length = names[0][fields[0]]
        if int(fields[1]) < 0 or int(fields[2]) > length or int(fields[1]) >= int(fields[2]):
            raise PilotError(f"P07: invalid partition coordinate at row {number}")
        used.add(fields[0])
        partition_rows += 1
    if (
        paf_rows <= 0
        or len(used) != int(production["resolved_ids"])
        or sha256_file(P07_GRAPH_LEDGER) != production["source_ledger_sha256"]
    ):
        raise PilotError("P07 independent graph identifier census differs")
    return {
        **production,
        "paf_rows": paf_rows,
        "partition_rows": partition_rows,
        "ledger_sequence_records": len(h1) + len(h2),
        "resolved_ids": len(used),
        "independently_recomputed": True,
    }


def audit_mapping(selection_id: str, root: Path) -> dict[str, object]:
    path = MAPPING_ROOTS[selection_id] / "h2_to_h1.1to1.paf"
    h1, h2 = _graph_lengths(selection_id, root)
    strands = Counter()
    rows = 0
    query_intervals: list[tuple[str, int, int]] = []
    target_intervals: list[tuple[str, int, int]] = []
    reverse_transforms = 0
    for number, line in enumerate(path.open(), 1):
        if not line.strip():
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 12 or fields[4] not in {"+", "-"}:
            raise PilotError(f"{selection_id}: malformed PAF row {number}")
        qlen, qs, qe, tlen, ts, te = map(
            int, (fields[1], fields[2], fields[3], fields[6], fields[7], fields[8])
        )
        if (
            fields[0] not in h2 or fields[5] not in h1
            or qlen != h2.get(fields[0]) or tlen != h1.get(fields[5])
            or not (0 <= qs < qe <= qlen and 0 <= ts < te <= tlen)
        ):
            raise PilotError(f"{selection_id}: invalid PAF coordinate at row {number}")
        if fields[4] == "-":
            oriented_start, oriented_end = qlen - qe, qlen - qs
            if not 0 <= oriented_start < oriented_end <= qlen:
                raise PilotError(f"{selection_id}: reverse-strand transform failed at row {number}")
            reverse_transforms += 1
        query_intervals.append((fields[0], qs, qe))
        target_intervals.append((fields[5], ts, te))
        strands[fields[4]] += 1
        rows += 1
    if not rows:
        raise PilotError(f"{selection_id}: exact mapping is empty")
    query_depth = _axis_overlap_depth(query_intervals)
    target_depth = _axis_overlap_depth(target_intervals)
    if query_depth > 1 or target_depth > 1:
        raise PilotError(f"{selection_id}: mapping is not deterministic bidirectional 1:1")
    return {
        "paf_rows": rows, "invalid_coordinates": 0,
        "strand_counts": dict(strands), "orientation": "H2_query_to_H1_reference",
        "reverse_strand_transforms_checked": reverse_transforms,
        "maximum_query_overlap_depth": query_depth,
        "maximum_target_overlap_depth": target_depth,
        "bidirectional_one_to_one_verified": True,
        "paf_sha256": sha256_file(path),
    }


def audit_psmc_population(root: Path) -> dict[str, object]:
    counts = Counter()
    records = 0
    for line in (root / "consensus/consensus/input.psmcfa").open():
        if line.startswith(">"):
            records += 1
            continue
        counts.update(line.strip())
    invalid = set(counts) - set("NKT")
    if invalid or not counts or records <= 0:
        raise PilotError(f"invalid complete primary PSMCFA population: {sorted(invalid)}")
    replicate_files = 0
    for replicate in range(1, 201):
        path = root / f"psmc/replicate-{replicate:03d}/bootstrap.unscaled.psmc"
        if not path.is_file() or path.stat().st_size <= 0:
            raise PilotError(f"missing finite PSMC bootstrap output: {replicate}")
        replicate_files += 1
    return {
        "primary_psmcfa_records": records,
        "complete_primary_bin_population": {
            symbol: counts[symbol] for symbol in "NKT"
        },
        "invalid_population_symbols": [],
        "nonempty_bootstrap_outputs": replicate_files,
        "population_independently_recomputed": True,
    }


def audit_pair(
    selection: Mapping[str, object], bcftools: str,
) -> dict[str, object]:
    selection_id = str(selection["selection_id"])
    root = BOUNDED_BASE / selection_id / "bounded-production"
    execution = json.loads((root / "execution.json").read_text())
    if execution.get("actual_core_biological_result") is not True:
        raise PilotError(f"{selection_id}: core biological result absent")
    plan = json.loads((root / "index/range_plan.json").read_text())
    required = (
        plan.get("range_plan_disjoint") is True
        and plan.get("range_plan_exhaustive") is True
        and plan.get("native_partition_one_owner") is True
        and plan.get("global_partition_assignment_ledger_materialized") is False
        and plan.get("global_impg_lace_created") is False
    )
    if not required:
        raise PilotError(f"{selection_id}: bounded plan gate failed")
    independent_plan = audit_range_plan(root, plan)
    graph = json.loads((root / "index/graph_identifier_audit.json").read_text())
    if graph.get("unresolved_ids") != 0 or graph.get("silently_omitted_regions") != 0:
        raise PilotError(f"{selection_id}: graph identifier audit failed")
    independent_graph = audit_graph_ids_independent(selection_id, root, graph)
    range_audit = audit_range_variants(root, plan)
    callable_audit = audit_callable_ownership(root, plan)
    mask_audit = audit_mask_accounting(root, independent_plan["h1_lengths"])
    block_audit = json.loads(
        (root / "consensus/bounded_consensus_block_audit.json").read_text()
    )
    if (
        block_audit.get("range_count") != plan["range_count"]
        or block_audit.get("range_blocks_reconstruct_global_sequence") is not True
        or block_audit.get("primary_psmcfa_bins_cross_contigs") is not False
        or block_audit.get("bootstrap_units_cross_contigs") is not False
    ):
        raise PilotError(f"{selection_id}: bounded consensus/PSMCFA block audit failed")
    psmc = json.loads((root / "psmc/finalize/psmc_qc.json").read_text())
    if (
        psmc.get("finite_bootstraps") != 200
        or psmc.get("primary_theta_centered") is not True
        or psmc.get("passed") is not True
    ):
        raise PilotError(f"{selection_id}: PSMC finite/centering gate failed")
    scenario_path = root / "psmc/finalize/scenario_scaled_trajectories.tsv"
    with scenario_path.open(newline="") as handle:
        scenarios = {row["scenario_id"] for row in csv.DictReader(handle, delimiter="\t")}
    if len(scenarios) != 9:
        raise PilotError(f"{selection_id}: expected nine PSMC scenarios")
    psmc_population = audit_psmc_population(root)
    annotation: Mapping[str, object]
    if selection_id == "P07":
        annotation = json.loads((root / "annotation/exact_partitions.json").read_text())
        if annotation.get("annotation_status") != "exact_native":
            raise PilotError("P07 exact native annotation absent")
    else:
        annotation = {"annotation_status": "missing_nonblocking", "core_outputs_stopped": False}
    normalized = root / "variants/normalized.vcf.gz"
    bcf = root / "variants/normalized.bcf"
    if not Path(str(normalized) + ".tbi").is_file() or not Path(str(bcf) + ".csi").is_file():
        raise PilotError(f"{selection_id}: normalized index absent")
    mapping = audit_mapping(selection_id, root)
    reconstruction = json.loads(
        (root / "variants/ref_alt_coordinate_audit.json").read_text()
    )
    if (
        reconstruction.get("invalid_coordinates") != 0
        or reconstruction.get("normalized_ref_mismatches") != 0
        or reconstruction.get("ref_alt_reconstruction_failures") != 0
        or int(reconstruction.get("paf_rows", -1)) != mapping["paf_rows"]
        or reconstruction.get("strand_counts") != mapping["strand_counts"]
    ):
        raise PilotError(f"{selection_id}: coordinate/REF-ALT production audit failed")
    return {
        "selection_id": selection_id, "species": selection["species"],
        "individual_or_isolate": selection["individual_or_isolate"],
        "failure_class": selection["failure_class"],
        "actual_core_biological_result": True,
        "diversity": execution["diversity"], "closed_output_ledger": audit_closed(root),
        "bounded_range_plan": {
            key: plan[key] for key in (
                "range_count", "h1_total_bp", "partitioned_h1_bp", "nonquery_h1_bp",
                "native_h1_partition_rows", "maximum_range_bp", "range_plan_disjoint",
                "range_plan_exhaustive", "native_partition_one_owner",
                "global_partition_assignment_ledger_materialized",
                "global_impg_lace_created",
            )
        },
        "independent_range_plan_audit": {
            key: independent_plan[key] for key in (
                "range_count", "h1_contig_count", "h1_total_bp",
                "native_h1_partition_rows", "partitioned_h1_bp",
                "range_plan_disjoint", "range_plan_exhaustive",
                "native_partition_one_owner",
            )
        },
        "range_variant_audit": range_audit,
        "callable_ownership_audit": callable_audit,
        "independent_mask_accounting": mask_audit,
        "consensus_block_audit": block_audit,
        "graph_identifier_audit": independent_graph,
        "coordinate_and_strand_audit": mapping,
        "ref_alt_reconstruction": reconstruction,
        "normalized_variants": {
            "vcf_gz_sha256": sha256_file(normalized), "bcf_sha256": sha256_file(bcf),
            "vcf_index_present": True, "bcf_index_present": True,
            "construction": "bcftools concat of normalized nonoverlapping range outputs",
        },
        "consensus_sha256": sha256_file(root / "consensus/consensus/consensus.fa"),
        "mask": execution["mask"], "psmc": psmc,
        "independent_psmc_population_audit": psmc_population,
        "scenario_count": len(scenarios),
        "independent_stratified_reaudit": stratified_reaudit(
            selection_id, root, plan, bcftools
        ),
        "annotation": annotation,
    }


def pipeline_limitations() -> list[str]:
    return [
        "FastGA remains unreliable for P03 at whole-assembly scale; the pinned WFMASH fallback is infrastructure provenance, not a biological exclusion.",
        "SweepGA/FastGA and WFMASH differ in raw gap placement on the controlled overlap despite high shared target coverage.",
        "Pinned IMPG lace with one thread did not progress; each bounded range uses the tested two-thread minimum.",
        "P02 and P03 lack cataloged exact native annotations; their core range, diversity, consensus, and PSMC outputs continue.",
        "Assembly-derived same-individual haplotype diversity is not a substitute for population sampling.",
    ]


def audit_all(selection_path: Path, output_path: Path, bcftools: str) -> dict[str, object]:
    freeze = json.loads(selection_path.read_text())
    by_id = {row["selection_id"]: row for row in freeze["selections"]}
    pairs = [audit_pair(by_id[key], bcftools) for key in ("P07", "P03", "P02")]
    comparison = json.loads(Path(
        "/moosefs/erikg/vgp/pilot/three-pair/vgp-three-pair-20260722-v1/"
        "P03/backend-comparison/comparison/comparison.json"
    ).read_text())
    value = {
        "schema_version": "vgp-three-pair-bounded-execution-v1",
        "task_id": "run-vgp-three-pair", "run_id": RUN_ID,
        "selection_freeze_sha256": sha256_file(selection_path),
        "actual_core_biological_results": len(pairs),
        "completion_gate_passed": len(pairs) == 3,
        "pairs": pairs, "controlled_fastga_wfmash_comparison": comparison,
        "canceled_global_jobs": [
            "1797004", "1797005", "1797006", "1797007", "1797008",
            "1797029", "1797030", "1797031", "1797032", "1797033",
        ],
        "cancellations_are_technical_not_biological_exclusions": True,
        "remaining_pipeline_limitations": pipeline_limitations(),
        "limitations_are_not_biological_exclusions": True,
    }
    output_path.write_text(canonical_json(value))
    return value


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--selection", type=Path, default=Path("analysis/vgp_three_pair_selection_v1.json")
    )
    parser.add_argument(
        "--output", type=Path, default=Path("analysis/vgp_three_pair_execution_v2.json")
    )
    parser.add_argument("--bcftools", required=True)
    args = parser.parse_args(argv)
    try:
        value = audit_all(args.selection, args.output, args.bcftools)
    except (PilotError, OSError, KeyError, TypeError, ValueError, json.JSONDecodeError,
            subprocess.CalledProcessError) as error:
        print(f"ERROR: {error}", file=__import__("sys").stderr)
        return 2
    print(json.dumps({
        "actual_core_biological_results": value["actual_core_biological_results"],
        "completion_gate_passed": value["completion_gate_passed"],
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
