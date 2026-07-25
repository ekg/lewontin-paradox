#!/usr/bin/env python3
"""Render the bounded three-pair execution audit as a concise evidence report."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Sequence

from analysis.vgp_10_pilot import PilotError, canonical_json, sha256_file


def materialize_independent_reaudit(execution: dict, output: Path) -> None:
    value = {
        "schema_version": "vgp-three-pair-bounded-independent-reaudit-v1",
        "task_id": "run-vgp-three-pair",
        "actual_core_biological_results": execution["actual_core_biological_results"],
        "completion_gate_passed": execution["completion_gate_passed"],
        "pairs": [
            {
                key: pair[key] for key in (
                    "selection_id", "frozen_source_audit", "closed_output_ledger",
                    "independent_range_plan_audit", "range_variant_audit",
                    "callable_ownership_audit", "independent_mask_accounting",
                    "graph_identifier_audit", "coordinate_and_strand_audit",
                    "ref_alt_reconstruction", "consensus_block_audit",
                    "independent_psmc_population_audit",
                    "independent_stratified_reaudit",
                )
            }
            for pair in execution["pairs"]
        ],
    }
    output.write_text(canonical_json(value))


def render(
    execution_path: Path, transition_path: Path, sacct_path: Path, output: Path,
) -> None:
    execution = json.loads(execution_path.read_text())
    transition = json.loads(transition_path.read_text())
    if execution.get("actual_core_biological_results") != 3:
        raise PilotError("report requires three actual bounded biological results")
    with sacct_path.open(newline="") as handle:
        telemetry = list(csv.DictReader(handle, delimiter="\t"))
    canary = transition["bounded_canary"]
    lines = [
        "# Three-pair VGP bounded-range reliability pilot",
        "",
        "## Outcome",
        "",
        "Three actual same-individual H1/H2 biological results completed. Every result was "
        "derived from bounded H1-coordinate IMPG queries; no all-genome graph, exhaustive "
        "all-partition query, global IMPG lace, or global partition-assignment ledger was "
        "created. The optional aggregate VCF/BCF is a `bcftools concat` of already normalized, "
        "nonoverlapping range products.",
        "",
        "| Pair | Failure class | Heterozygous SNPs | Callable bp | Pi | Ranges | PSMC | Annotation |",
        "|---|---|---:|---:|---:|---:|---|---|",
    ]
    for pair in execution["pairs"]:
        diversity, plan, psmc = pair["diversity"], pair["bounded_range_plan"], pair["psmc"]
        annotation = pair["annotation"]["annotation_status"]
        lines.append(
            f"| {pair['selection_id']} ({pair['species']}) | {pair['failure_class']} | "
            f"{diversity['heterozygous_snps']:,} | {diversity['callable_bp']:,} | "
            f"{diversity['pi']:.12g} | {plan['range_count']:,} | "
            f"{psmc['finite_bootstraps']}/200 finite; centered={str(psmc['primary_theta_centered']).lower()} | "
            f"{annotation} |"
        )
    lines.extend([
        "",
        "## Architecture correction and canary",
        "",
        f"The prohibited global chains were canceled: {', '.join(execution['canceled_global_jobs'])}. "
        "Their cancellation is technical provenance and is not a species exclusion.",
        "",
        f"The fresh P07 canary queried `{canary['contig']}:{canary['start']}-{canary['end']}` "
        f"({canary['range_bp']:,} bp) using {canary['native_partition_count']:,} complete native "
        f"partitions. Its {canary['normalized_variant_keys']:,} normalized IMPG keys exactly "
        f"matched the like-for-like clean P07 subset (SHA-256 `{canary['exact_variant_key_sha256']}`). "
        f"Callable accounting was {canary['callable_bp']:,} bp, peak local graph state was "
        f"{canary['peak_local_graph_state_bytes']:,} bytes, and local graph temporaries were deleted.",
        "",
        "The frozen range rule is complete consecutive native H1 partition boundaries, targeting "
        "5 Mb with a hard 20 Mb ceiling. Unaligned H1 ranges remain explicit non-query ranges and "
        "enter consensus/PSMC as non-callable sequence; they are never silently omitted.",
        "",
        "P02 job 1799308 exposed a graph-allele REF/ALT orientation error in the first range. "
        "The corrected job 1804024 completed all 925 query ranges, then correctly rejected the "
        "H1-monomorphic `ALT=.` record at CM078529.1:2167395 during reduction. Its zero-copy-failure "
        "bundle retained all 7,412 result files plus the exact plan and graph index. Resume job "
        "1804979 removed only post-reconstruction monomorphic records, strict-H1-revalidated all "
        "925 VCF/BCF pairs, issued zero IMPG queries, and completed the biological result. These "
        "failures and the resume are technical provenance, not an exclusion.",
        "",
        "## Independent validation",
        "",
    ])
    for pair in execution["pairs"]:
        lines.append(f"### {pair['selection_id']}")
        lines.append("")
        lines.append(
            f"Graph IDs unresolved={pair['graph_identifier_audit']['unresolved_ids']}; "
            f"silently omitted={pair['graph_identifier_audit']['silently_omitted_regions']}. "
            f"Range duplicate keys={pair['range_variant_audit']['exact_duplicate_keys_between_ranges']}; "
            f"unowned callable bp={pair['callable_ownership_audit']['unowned_callable_bp']}; "
            f"multiply owned callable bp={pair['callable_ownership_audit']['multiply_owned_callable_bp']}."
        )
        lines.append("")
        lines.append("| Stratum | Range | Variants | Callable bp | Consensus non-N bp | PSMC N/K/T |")
        lines.append("|---|---|---:|---:|---:|---|")
        for row in pair["independent_stratified_reaudit"]:
            pop = row["psmc_population_bins"]
            lines.append(
                f"| {row['stratum']} | {row['range_id']} "
                f"({row['contig']}:{row['start']}-{row['end']}) | "
                f"{row['variant_records']:,} | {row['callable_bp']:,} | "
                f"{row['consensus_non_N_bp']:,} | {pop['N']:,}/{pop['K']:,}/{pop['T']:,} |"
            )
        lines.append("")
    p07 = next(pair for pair in execution["pairs"] if pair["selection_id"] == "P07")
    annotation_rows = p07["annotation"].get("partitions", [])
    if annotation_rows:
        lines.extend([
            "## Exact native annotation partitions",
            "",
            "P07’s GFF is exact-native and dictionary-identical to H1. Annotation ranges were "
            "queried/intersected directly in H1 coordinates; P02/P03 lack cataloged exact-native "
            "annotations, so only their annotation products are absent.",
            "",
            "| Partition | Heterozygous variants | Callable bp | Estimate | Estimator |",
            "|---|---:|---:|---:|---|",
        ])
        for row in annotation_rows:
            lines.append(
                f"| {row['partition']} | {row['heterozygous_variants']:,} | "
                f"{row['callable_bp']:,} | {row['estimate']:.12g} | "
                f"{row['estimator']} |"
            )
        lines.append("")
    comparison = execution["controlled_fastga_wfmash_comparison"]
    lines.extend([
        "## Controlled backend comparison",
        "",
        f"P03’s byte-identical staged-FASTA control has {comparison['overlapping_target_bp']:,} bp "
        f"of common target coverage, target-coverage Jaccard "
        f"{comparison['target_coverage_jaccard']:.6f}, and exact raw-variant Jaccard "
        f"{comparison['exact_variant_jaccard']:.6f}. Both reconstructions have zero coordinate/"
        "REF-ALT failures. Backend-specific gap placement is retained as a limitation rather than "
        "a biological exclusion.",
        "",
        "## Resource telemetry",
        "",
        "| Job | Name | State | Elapsed | CPUs | MaxRSS | Node |",
        "|---|---|---|---|---:|---:|---|",
    ])
    for row in telemetry:
        lines.append(
            f"| {row.get('JobIDRaw','')} | {row.get('JobName','')} | {row.get('State','')} | "
            f"{row.get('Elapsed','')} | {row.get('AllocCPUS','')} | {row.get('MaxRSS','')} | "
            f"{row.get('NodeList','')} |"
        )
    lines.extend(["", "## Remaining limitations", ""])
    lines.extend(f"- {item}" for item in execution["remaining_pipeline_limitations"])
    lines.extend([
        "",
        "All limitations above are technical or interpretive limitations, not biological exclusions.",
        "",
        "## Evidence identities",
        "",
        f"- Selection freeze SHA-256: `{execution['selection_freeze_sha256']}`",
        f"- Execution audit SHA-256: `{sha256_file(execution_path)}`",
        f"- Bounded transition record SHA-256: `{sha256_file(transition_path)}`",
        f"- Scheduler telemetry SHA-256: `{sha256_file(sacct_path)}`",
        "",
    ])
    output.write_text("\n".join(lines), encoding="utf-8")


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--execution", type=Path, default=Path("analysis/vgp_three_pair_execution_v2.json")
    )
    parser.add_argument(
        "--transition", type=Path,
        default=Path("analysis/vgp_three_pair_bounded_transition_v1.json"),
    )
    parser.add_argument(
        "--sacct", type=Path, default=Path("analysis/vgp_three_pair_bounded_sacct_v1.tsv")
    )
    parser.add_argument(
        "--output", type=Path, default=Path("analysis/vgp_three_pair_bounded_report_v1.md")
    )
    parser.add_argument(
        "--reaudit-output", type=Path,
        default=Path("analysis/vgp_three_pair_independent_reaudit_v2.json"),
    )
    args = parser.parse_args(argv)
    try:
        render(args.execution, args.transition, args.sacct, args.output)
        materialize_independent_reaudit(
            json.loads(args.execution.read_text()), args.reaudit_output
        )
    except (PilotError, OSError, KeyError, TypeError, ValueError, json.JSONDecodeError) as error:
        print(f"ERROR: {error}", file=__import__("sys").stderr)
        return 2
    print(json.dumps({
        "independent_reaudit": str(args.reaudit_output),
        "report": str(args.output),
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
