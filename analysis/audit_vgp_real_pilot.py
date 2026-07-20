#!/usr/bin/env python3
"""Build fail-closed, closed-world accounting for the ten real VGP pairs."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Mapping, Sequence


PAIR_IDS = tuple(f"P{number:02d}" for number in range(1, 11))
REMAINING_IDS = tuple(pair for pair in PAIR_IDS if pair != "P07")
CORE_STAGES = ("mapping", "impg", "variants", "consensus", "psmc_finalize")
TERMINAL_PRIMARY_FAILURES = {"reproducible_hard_primary_execution_failure"}


class PilotAuditError(RuntimeError):
    """The closed-world pilot packet is incomplete or internally inconsistent."""


def _json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def _tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _count_tsv(path: Path) -> int:
    with path.open(encoding="utf-8") as handle:
        return max(0, sum(1 for line in handle if line.strip()) - 1)


def _annotation_expected(input_root: Path, pair: str) -> bool:
    return _json(input_root / pair / "input-manifest.json").get("annotation") is not None


def summarize_pair(
    pair: str,
    run_root: Path,
    submissions: Sequence[Mapping[str, str]],
    input_root: Path,
) -> dict[str, object]:
    pair_rows = [row for row in submissions if row.get("selection_id") == pair]
    submitted_stages = sorted({row.get("stage", "") for row in pair_rows})
    sentinels = {
        "preflight": run_root / "preflight/.complete.json",
        "mapping": run_root / "mapping/.complete.json",
        "impg": run_root / "impg/.complete.json",
        "variants": run_root / "variants/.complete.json",
        "consensus": run_root / "consensus/.complete.json",
        "psmc_finalize": run_root / "psmc/finalize/.complete.json",
        "annotation": run_root / "annotation/.complete.json",
    }
    stage_complete = {name: path.is_file() for name, path in sentinels.items()}
    replicate_sentinels = list((run_root / "psmc").glob("replicate-[0-9][0-9][0-9]/.complete.json"))
    complete = all(stage_complete[stage] for stage in CORE_STAGES)
    annotation_expected = _annotation_expected(input_root, pair)
    if annotation_expected:
        complete = complete and stage_complete["annotation"]
    result: dict[str, object] = {
        "selection_id": pair,
        "canonical_run_root": str(run_root),
        "submission_count": len(pair_rows),
        "submitted_stages": submitted_stages,
        "submitted_job_ids": [row.get("job_id") for row in pair_rows],
        "stage_complete": stage_complete,
        "psmc_replicate_sentinels": len(replicate_sentinels),
        "annotation_expected": annotation_expected,
        "status": "complete" if complete else "submitted_incomplete",
    }
    if not complete:
        return result
    mask = _json(run_root / "consensus/masks/mask_reconciliation.json")
    consensus = _json(run_root / "consensus/join_qc.json")["consensus"]
    psmc_qc = _json(run_root / "psmc/finalize/psmc_qc.json")
    callable_bp = int(consensus["consensus_callable_bp"])
    snps = int(consensus["heterozygous_snps"])
    result.update({
        "callability": {
            "universe_bp": int(mask["universe_bp"]),
            "pre_indel_callable_bp": int(mask["callable_bp"]),
            "final_callable_bp": callable_bp,
            "final_callable_fraction": callable_bp / int(mask["universe_bp"]),
        },
        "diversity": {
            "heterozygous_snps": snps,
            "pi": snps / callable_bp,
            "estimator": "heterozygous SNPs per final callable H1 bp",
        },
        "psmc": {
            "bootstrap_attempts": int(psmc_qc["bootstrap_attempts"]),
            "finite_bootstraps": int(psmc_qc["finite_bootstraps"]),
            "qc_passed": psmc_qc["passed"] is True,
            "unscaled_intervals": _count_tsv(run_root / "psmc/finalize/unscaled_trajectory.tsv"),
            "scenario_rows": _count_tsv(run_root / "psmc/finalize/scenario_scaled_trajectories.tsv"),
            "unscaled_primary_preserved": psmc_qc["unscaled_primary_preserved"] is True,
        },
    })
    if annotation_expected:
        result["annotation"] = _json(run_root / "annotation/exact_partitions.json")
    return result


def canary_summary(canary: Mapping[str, object]) -> dict[str, object]:
    psmc = canary["psmc"]
    diversity = canary["diversity"]
    return {
        "selection_id": "P07",
        "canonical_run_root": canary["promoted_run_root"],
        "status": "complete",
        "stage_complete": {stage: True for stage in (*CORE_STAGES, "preflight", "annotation")},
        "submission_count": len(canary["telemetry"]["rows"]),
        "submitted_stages": ["canary_core", "annotation"],
        "submitted_job_ids": [row["JobIDRaw"] for row in canary["telemetry"]["rows"]],
        "psmc_replicate_sentinels": int(psmc["bootstrap_attempts"]) + 1,
        "annotation_expected": True,
        "callability": {
            "final_callable_bp": int(diversity["callable_bp"]),
            "pre_indel_callable_bp": int(diversity["pre_indel_callable_bp"]),
        },
        "diversity": {
            "heterozygous_snps": int(diversity["heterozygous_snps"]),
            "pi": float(diversity["pi"]),
            "estimator": diversity["estimator"],
        },
        "psmc": {
            "bootstrap_attempts": int(psmc["bootstrap_attempts"]),
            "finite_bootstraps": int(psmc["finite_bootstraps"]),
            "qc_passed": int(psmc["finite_bootstraps"]) >= 100,
            "unscaled_intervals": int(psmc["trajectory_intervals"]),
            "scenario_rows": int(psmc["rows"]),
            "unscaled_primary_preserved": True,
        },
        "annotation": canary["annotation"],
    }


def apply_terminal_failures(
    rows: Sequence[dict[str, object]], failures: Mapping[str, object]
) -> None:
    """Attach retained terminal primary failures without hiding prior attempts."""

    failure_rows = failures.get("failures")
    if not isinstance(failure_rows, list):
        raise PilotAuditError("execution failure ledger lacks a failure list")
    for row in rows:
        if row["status"] == "complete":
            continue
        terminal = [
            failure for failure in failure_rows
            if isinstance(failure, dict)
            and failure.get("selection_id") == row["selection_id"]
            and failure.get("classification") in TERMINAL_PRIMARY_FAILURES
            and failure.get("primary_preserved") is True
            and failure.get("alternate_activated") is False
            and failure.get("terminal_primary_failure") is not False
        ]
        if terminal:
            row["status"] = "failed_primary"
            row["terminal_failure"] = terminal[-1]


def validate_closed_world(rows: Sequence[Mapping[str, object]]) -> None:
    by_pair = {str(row["selection_id"]): row for row in rows}
    if set(by_pair) != set(PAIR_IDS) or len(rows) != len(PAIR_IDS):
        raise PilotAuditError("closed world must contain exactly P01..P10 once")
    for pair in REMAINING_IDS:
        row = by_pair[pair]
        submitted = set(row["submitted_stages"])
        missing = {"mapping", "impg", "variants", "consensus", "psmc", "psmc_finalize"} - submitted
        if missing:
            raise PilotAuditError(f"{pair} lacks real biological submissions: {sorted(missing)}")
    newly_complete = [by_pair[pair] for pair in REMAINING_IDS if by_pair[pair]["status"] == "complete"]
    if not newly_complete:
        raise PilotAuditError("no newly completed biological pair")
    for row in [by_pair["P07"], *newly_complete]:
        if float(row["diversity"]["pi"]) <= 0 or int(row["callability"]["final_callable_bp"]) <= 0:
            raise PilotAuditError(f"{row['selection_id']} has a zero diversity result")
        psmc = row["psmc"]
        if int(psmc["finite_bootstraps"]) < 100 or int(psmc["unscaled_intervals"]) <= 0:
            raise PilotAuditError(f"{row['selection_id']} lacks finite PSMC support")
        if int(psmc["scenario_rows"]) <= 0 or psmc["unscaled_primary_preserved"] is not True:
            raise PilotAuditError(f"{row['selection_id']} lacks separate scaling scenarios")


def build(args: argparse.Namespace) -> dict[str, object]:
    submissions = _tsv(args.submissions)
    if any(row.get("canonical_vgp_root") != args.canonical_root for row in submissions):
        raise PilotAuditError("submission manifest contains a noncanonical VGP root")
    failures = _json(args.failures)
    if failures.get("canonical_vgp_root") != args.canonical_root:
        raise PilotAuditError("execution failure ledger contains a noncanonical VGP root")
    canary = _json(args.canary_execution)
    if canary.get("canonical_vgp_root") != args.canonical_root:
        raise PilotAuditError("canary execution contains a noncanonical VGP root")
    rows = [canary_summary(canary)]
    rows.extend(summarize_pair(
        pair, args.run_root / pair, submissions, args.input_root
    ) for pair in REMAINING_IDS)
    rows.sort(key=lambda row: str(row["selection_id"]))
    apply_terminal_failures(rows, failures)
    validate_closed_world(rows)
    result = {
        "schema_version": "vgp-real-pilot-closed-world-v1",
        "task_id": "run-vgp-real-pilot",
        "authorization_id": "vgp10-auth-20260718-v2",
        "canonical_vgp_root": args.canonical_root,
        "run_id": args.run_id,
        "pair_count": len(rows),
        "complete_pair_count": sum(row["status"] == "complete" for row in rows),
        "newly_complete_pair_count": sum(
            row["status"] == "complete" and row["selection_id"] != "P07" for row in rows
        ),
        "pairs": rows,
        "execution_failures": failures,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    partial = args.output.with_suffix(args.output.suffix + ".partial")
    partial.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    partial.replace(args.output)
    return result


def parser() -> argparse.ArgumentParser:
    value = argparse.ArgumentParser(description=__doc__)
    value.add_argument("--canonical-root", default="/moosefs/erikg/vgp")
    value.add_argument("--run-id", default="vgp10-auth-20260718-v2-pilot-v1")
    value.add_argument("--run-root", type=Path, required=True)
    value.add_argument("--input-root", type=Path, required=True)
    value.add_argument("--submissions", type=Path, required=True)
    value.add_argument("--failures", type=Path, required=True)
    value.add_argument("--canary-execution", type=Path, required=True)
    value.add_argument("--output", type=Path, required=True)
    return value


def main(argv: Sequence[str] | None = None) -> int:
    args = parser().parse_args(argv)
    try:
        result = build(args)
    except (OSError, ValueError, KeyError, json.JSONDecodeError, PilotAuditError) as error:
        print(f"ERROR: {error}")
        return 2
    print(json.dumps({
        "output": str(args.output),
        "pair_count": result["pair_count"],
        "complete_pair_count": result["complete_pair_count"],
        "newly_complete_pair_count": result["newly_complete_pair_count"],
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
