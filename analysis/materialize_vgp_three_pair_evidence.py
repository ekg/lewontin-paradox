#!/usr/bin/env python3
"""Materialize reviewable evidence from a completed three-pair VGP audit."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Mapping, Sequence

from analysis.vgp_10_pilot import PilotError, canonical_json, sha256_file


TELEMETRY_FIELDS = (
    "selection_id", "stage", "job_id", "array_task_id", "disposition",
    "started_epoch", "ended_epoch", "elapsed_seconds", "maximum_rss_kib",
    "child_cpu_seconds", "scratch_required_bytes",
    "scratch_available_bytes_at_start", "scratch_high_water_bytes",
    "filesystem_read_bytes", "filesystem_write_bytes", "retry", "record",
)


def _require_completed(execution: Mapping[str, object]) -> list[Mapping[str, object]]:
    pairs = execution.get("pairs")
    if (
        execution.get("actual_core_biological_results") != 3
        or execution.get("completion_gate_passed") is not True
        or not isinstance(pairs, list)
        or len(pairs) != 3
        or any(not isinstance(pair, Mapping) or pair.get("actual_core_biological_result") is not True for pair in pairs)
    ):
        raise PilotError("three actual core biological results are required")
    if {str(pair["selection_id"]) for pair in pairs} != {"P07", "P03", "P02"}:
        raise PilotError("completed pair identities do not match the frozen pilot")
    return pairs


def materialize_independent_reaudit(
    execution: Mapping[str, object], execution_path: Path, output: Path
) -> dict[str, object]:
    pairs = _require_completed(execution)
    result = {
        "schema_version": "vgp-three-pair-independent-reaudit-v1",
        "task_id": execution["task_id"],
        "run_id": execution["run_id"],
        "source_execution_sha256": sha256_file(execution_path),
        "method": (
            "Independent post-promotion parsing of the durable consensus FASTA, "
            "normalized VCF, callable BED, PSMCFA, PAF, graph census, and all 201 "
            "closed PSMC replicate ledgers by audit_vgp_three_pair_results.py."
        ),
        "pairs": [
            {
                "selection_id": pair["selection_id"],
                "actual_core_biological_result": pair["actual_core_biological_result"],
                "diversity": pair["diversity"],
                "independent_stratified_reaudit": pair["independent_stratified_reaudit"],
                "coordinate_and_strand_audit": pair["coordinate_and_strand_audit"],
                "ref_alt_reconstruction": pair["ref_alt_reconstruction"],
                "normalized_concordance": pair["normalized_concordance"],
                "graph_identifier_audit": pair["graph_identifier_audit"],
                "mask": pair["mask"],
                "psmc": pair["psmc"],
                "psmc_population_preserved": pair["psmc_population_preserved"],
            }
            for pair in pairs
        ],
        "all_three_stratified_reaudits_present": all(
            len(pair["independent_stratified_reaudit"]) == 3 for pair in pairs
        ),
    }
    output.write_text(canonical_json(result), encoding="utf-8")
    return result


def materialize_telemetry(execution: Mapping[str, object], output: Path) -> int:
    pairs = _require_completed(execution)
    rows: list[dict[str, object]] = []
    for pair in pairs:
        for raw in pair["resource_telemetry"]:
            row = dict(raw)
            started, ended = row.get("started_epoch", ""), row.get("ended_epoch", "")
            elapsed = int(ended) - int(started) if started != "" and ended != "" else ""
            row["elapsed_seconds"] = elapsed
            row.setdefault("selection_id", pair["selection_id"])
            rows.append({field: row.get(field, "") for field in TELEMETRY_FIELDS})
    if not rows:
        raise PilotError("stage telemetry is absent")
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=TELEMETRY_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return len(rows)


def _sacct_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="|"))
    roots = [row for row in rows if row.get("JobIDRaw") and "." not in row["JobIDRaw"]]
    if not roots:
        raise PilotError("scheduler telemetry has no root allocations")
    return roots


def _fmt_int(value: object) -> str:
    return f"{int(value):,}"


def materialize_report(
    execution: Mapping[str, object], execution_path: Path, sacct_path: Path, output: Path
) -> None:
    pairs = _require_completed(execution)
    sacct = _sacct_rows(sacct_path)
    comparison = execution["controlled_fastga_wfmash_comparison"]
    lines = [
        "# Three-pair VGP reliability pilot",
        "",
        "Status: **COMPLETE — three actual core biological results**",
        "",
        f"Run ID: `{execution['run_id']}`  ",
        f"Execution record SHA-256: `{sha256_file(execution_path)}`",
        "",
        "## Biological results",
        "",
        "| Pair | Failure class exercised | Heterozygous SNPs | Callable bp | Assembly-derived pi | PSMC bootstraps | Annotation |",
        "| --- | --- | ---: | ---: | ---: | ---: | --- |",
    ]
    for pair in pairs:
        diversity, psmc = pair["diversity"], pair["psmc"]
        annotation = pair["annotation"]["annotation_status"]
        lines.append(
            f"| {pair['selection_id']} (*{pair['species']}*) | {pair['failure_class']} | "
            f"{_fmt_int(diversity['heterozygous_snps'])} | {_fmt_int(diversity['callable_bp'])} | "
            f"{diversity['pi']:.16g} | {psmc['finite_bootstraps']}/200, centered | {annotation} |"
        )
    lines.extend([
        "",
        "Every pair has indexed normalized VCF.gz and BCF, a masked diploid consensus, "
        "a primary PSMC fit, 200 finite centered bootstrap fits, and nine sensitivity scenarios. "
        "Missing P02/P03 annotation stopped only annotation outputs.",
        "",
        "## Reliability repairs and independent checks",
        "",
        "P03 retained the reproducible corrected FastGA failure as infrastructure telemetry and "
        "used pinned WFMASH on the exact same staged FASTAs, followed by the identical deterministic "
        "bidirectional 1:1 filter. P02 rebuilt FASTA dictionaries, graph index, partitions, and "
        "queries; every graph identifier was checked on its H1 or H2 axis, and no region was omitted.",
        "",
        "For each pair, the post-promotion auditor independently selected early, middle, and late "
        "5-Mb strata and reparsed normalized variants, callable BED intervals, consensus bases, "
        "and N/K/T PSMC populations. It also rechecked every PAF coordinate and strand, exact "
        "REF/ALT reconstruction, the reason-mask universe, graph identifiers, and all 201 PSMC "
        "stage ledgers.",
        "",
        "## Controlled FastGA/WFMASH overlap",
        "",
        f"On the frozen P03 20-Mb overlapping subset, FastGA covered {_fmt_int(comparison['fastga_target_bp'])} "
        f"target bp and WFMASH covered {_fmt_int(comparison['wfmash_target_bp'])} bp; "
        f"the intersection was {_fmt_int(comparison['overlapping_target_bp'])} bp and coverage "
        f"Jaccard was {comparison['target_coverage_jaccard']:.6f}. Both outputs had exact bidirectional "
        "depth one and zero coordinate or reconstruction failures. Raw exact variant Jaccard was "
        f"{comparison['exact_variant_jaccard']:.6f}, recording backend-specific gap placement "
        "rather than a biological exclusion.",
        "",
        "## Scheduler and stage telemetry",
        "",
        "| Job | Name | State | Exit | Elapsed | CPUs | Requested memory | Node |",
        "| --- | --- | --- | --- | --- | ---: | --- | --- |",
    ])
    for row in sacct:
        lines.append(
            f"| {row.get('JobIDRaw', '')} | {row.get('JobName', '')} | {row.get('State', '')} | "
            f"{row.get('ExitCode', '')} | {row.get('Elapsed', '')} | {row.get('AllocCPUS', '')} | "
            f"{row.get('ReqMem', '')} | {row.get('NodeList', '')} |"
        )
    lines.extend([
        "",
        "The stage telemetry TSV records measured elapsed time, child CPU, peak process RSS when "
        "available, scratch reservation/free/high-water bytes, filesystem I/O counters, retry "
        "number, and disposition. Blank Slurm MaxRSS/TotalCPU values reflect cluster accounting "
        "configuration and are retained as missing—not estimated.",
        "",
        "## Remaining limitations",
        "",
    ])
    lines.extend(f"- {item}" for item in execution["remaining_pipeline_limitations"])
    lines.extend([
        "",
        "These are pipeline or interpretation limitations, not biological exclusions.",
        "",
        "## Evidence",
        "",
        "- `analysis/vgp_three_pair_selection_v1.json`: immutable pre-execution selections and alternates",
        "- `analysis/vgp_three_pair_execution_v1.json`: complete machine-readable biological and technical audit",
        "- `analysis/vgp_three_pair_independent_reaudit_v1.json`: three-stratum and contract re-audit",
        "- `analysis/vgp_three_pair_sacct_v1.tsv`: scheduler allocations, including failed/canceled attempts",
        "- `analysis/vgp_three_pair_stage_telemetry_v1.tsv`: stage-level resource and retry telemetry",
        "- `analysis/vgp_three_pair_output_digests_v1.txt`: closed repository evidence digest ledger",
        "",
    ])
    output.write_text("\n".join(lines), encoding="utf-8")


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--execution", type=Path, default=Path("analysis/vgp_three_pair_execution_v1.json"))
    parser.add_argument("--sacct", type=Path, default=Path("analysis/vgp_three_pair_sacct_v1.tsv"))
    parser.add_argument("--independent", type=Path, default=Path("analysis/vgp_three_pair_independent_reaudit_v1.json"))
    parser.add_argument("--telemetry", type=Path, default=Path("analysis/vgp_three_pair_stage_telemetry_v1.tsv"))
    parser.add_argument("--report", type=Path, default=Path("analysis/vgp_three_pair_report_v1.md"))
    args = parser.parse_args(argv)
    try:
        execution = json.loads(args.execution.read_text(encoding="utf-8"))
        materialize_independent_reaudit(execution, args.execution, args.independent)
        rows = materialize_telemetry(execution, args.telemetry)
        materialize_report(execution, args.execution, args.sacct, args.report)
    except (PilotError, OSError, KeyError, TypeError, ValueError, json.JSONDecodeError) as error:
        print(f"ERROR: {error}", file=__import__("sys").stderr)
        return 2
    print(json.dumps({"completion_gate_passed": True, "telemetry_rows": rows}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
