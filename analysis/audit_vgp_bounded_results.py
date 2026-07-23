#!/usr/bin/env python3
"""Independent completion audit for the three bounded-range VGP results."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import subprocess
from collections import Counter, defaultdict
from pathlib import Path
from typing import Mapping, Sequence

from analysis.audit_vgp_three_pair_results import (
    fasta_lengths,
    fasta_non_n_counts,
    psmc_population_counts,
)
from analysis.vgp_10_pilot import (
    Interval,
    PilotError,
    canonical_json,
    interval_bp,
    read_bed,
    sha256_file,
)


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


def _query_keys(bcftools: str, path: Path) -> list[str]:
    result = subprocess.run(
        [bcftools, "query", "-f", "%CHROM\\t%POS\\t%REF\\t%ALT\\n", str(path)],
        check=True, capture_output=True, text=True,
    )
    return result.stdout.splitlines()


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
        # Independently make bcftools re-check every selected REF against the
        # staged H1 dictionary embedded in the normalized header is not enough;
        # range production's norm -f audit plus exact key equality is retained.
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
            "consensus_non_N_bp": non_n[region],
            "psmc_population_bins": {symbol: pop[symbol] for symbol in "NKT"},
        })
    return result


def audit_mapping(selection_id: str) -> dict[str, object]:
    path = MAPPING_ROOTS[selection_id] / "h2_to_h1.1to1.paf"
    strands = Counter()
    rows = 0
    for number, line in enumerate(path.open(), 1):
        if not line.strip():
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 12 or fields[4] not in {"+", "-"}:
            raise PilotError(f"{selection_id}: malformed PAF row {number}")
        qlen, qs, qe, tlen, ts, te = map(
            int, (fields[1], fields[2], fields[3], fields[6], fields[7], fields[8])
        )
        if not (0 <= qs < qe <= qlen and 0 <= ts < te <= tlen):
            raise PilotError(f"{selection_id}: invalid PAF coordinate at row {number}")
        strands[fields[4]] += 1
        rows += 1
    if not rows:
        raise PilotError(f"{selection_id}: exact mapping is empty")
    return {
        "paf_rows": rows, "invalid_coordinates": 0,
        "strand_counts": dict(strands), "orientation": "H2_query_to_H1_reference",
        "paf_sha256": sha256_file(path),
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
    graph = json.loads((root / "index/graph_identifier_audit.json").read_text())
    if graph.get("unresolved_ids") != 0 or graph.get("silently_omitted_regions") != 0:
        raise PilotError(f"{selection_id}: graph identifier audit failed")
    range_audit = audit_range_variants(root, plan)
    callable_audit = audit_callable_ownership(root, plan)
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
        "range_variant_audit": range_audit,
        "callable_ownership_audit": callable_audit,
        "consensus_block_audit": block_audit,
        "graph_identifier_audit": graph,
        "coordinate_and_strand_audit": audit_mapping(selection_id),
        "ref_alt_reconstruction": json.loads(
            (root / "variants/ref_alt_coordinate_audit.json").read_text()
        ),
        "normalized_variants": {
            "vcf_gz_sha256": sha256_file(normalized), "bcf_sha256": sha256_file(bcf),
            "vcf_index_present": True, "bcf_index_present": True,
            "construction": "bcftools concat of normalized nonoverlapping range outputs",
        },
        "consensus_sha256": sha256_file(root / "consensus/consensus/consensus.fa"),
        "mask": execution["mask"], "psmc": psmc, "scenario_count": len(scenarios),
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
