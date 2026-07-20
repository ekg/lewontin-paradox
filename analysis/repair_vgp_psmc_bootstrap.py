#!/usr/bin/env python3
"""Repair and independently diagnose real-VGP PSMC block bootstraps.

The repair population is the primary ``input.psmcfa`` itself.  Five-megabase
blocks are represented as 50,000-bin, contig-bounded slices and therefore
carry masked ``N`` and callable ``K/T`` bins together.  The centering rule
below is declared in source and in ``centering_diagnostic.predeclared.json``
before any repaired PSMC output is evaluated.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import shutil
import statistics
import subprocess
import tempfile
from collections import Counter
from pathlib import Path
from typing import Mapping, Sequence

from analysis.vgp_10_pilot import (
    Interval,
    PSMCFA_BOOTSTRAP_DESIGN_SHA256,
    PilotError,
    bootstrap_psmcfa,
    canonical_json,
    freeze_psmcfa_bootstrap_units,
    parse_psmc_unscaled,
    parse_psmcfa,
    psmcfa_bootstrap_manifest,
    sha256_file,
    write_bed,
)


CANONICAL_VGP_ROOT = Path("/moosefs/erikg/vgp")
REPAIR_ID = "vgp-psmc-bootstrap-repair-v1"
ATTEMPTS = 200
BIN_SIZE_BP = 100
BLOCK_BP = 5_000_000
BLOCK_BINS = BLOCK_BP // BIN_SIZE_BP
PAIR_RELATIVE_ROOTS = {
    "P04": Path("pilot/runs/vgp10-auth-20260718-v2-pilot-v1/P04"),
    "P07": Path("pilot/outputs/vgp10-auth-20260718-v2/P07/core"),
}
PSMC_OPTIONS = ("-b", "-N25", "-t15", "-r5", "-p", "4+25*2+4+6")

# This is deliberately a source constant, not a threshold selected after the
# repaired runs.  A repaired pair passes only when all 200 fits are finite and
# its unchanged primary theta lies in the equal-tail central 95% bootstrap
# interval (nearest-rank q=0.025 and q=0.975, matching the pilot review).
CENTERING_DIAGNOSTIC: Mapping[str, object] = {
    "name": "primary_theta_in_equal_tail_central_95pct_bootstrap_interval",
    "metric": "final_native_iteration_theta_0_per_100bp_bin",
    "lower_quantile": 0.025,
    "upper_quantile": 0.975,
    "quantile_method": "nearest index round((n - 1) * q)",
    "required_finite_outputs": ATTEMPTS,
    "required_attempts": ATTEMPTS,
    "predeclared_before_execution": True,
}


def pair_root(canonical_root: Path, pair: str) -> Path:
    if pair not in PAIR_RELATIVE_ROOTS:
        raise PilotError(f"unsupported repaired pair: {pair}")
    root = canonical_root / PAIR_RELATIVE_ROOTS[pair]
    if not str(root.resolve()).startswith(str(canonical_root.resolve()) + os.sep):
        raise PilotError(f"pair root escapes canonical VGP root: {root}")
    return root


def primary_path(canonical_root: Path, pair: str) -> Path:
    return pair_root(canonical_root, pair) / "consensus/consensus/input.psmcfa"


def primary_psmc_path(canonical_root: Path, pair: str) -> Path:
    return pair_root(canonical_root, pair) / "psmc/replicate-000/unscaled.psmc"


def symbol_counts(records: Mapping[str, str]) -> Counter[str]:
    return Counter(symbol for sequence in records.values() for symbol in sequence)


def unit_symbol_counts(records: Mapping[str, str], units: Sequence[Interval]) -> Counter[str]:
    return Counter(
        symbol
        for unit in units
        for symbol in records[unit.contig][unit.start:unit.end]
    )


def quantiles(values: Sequence[float]) -> dict[str, float]:
    ordered = sorted(values)
    if not ordered:
        raise PilotError("cannot summarize an empty diagnostic distribution")

    def pick(fraction: float) -> float:
        return ordered[round((len(ordered) - 1) * fraction)]

    return {
        "min": ordered[0],
        "q025": pick(0.025),
        "median": statistics.median(ordered),
        "q975": pick(0.975),
        "max": ordered[-1],
    }


def atomic_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", dir=path.parent, prefix=f".{path.name}.", delete=False
    ) as handle:
        partial = Path(handle.name)
        handle.write(text)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(partial, path)


def design_paths(output_root: Path, pair: str) -> tuple[Path, Path]:
    design = output_root / pair / "design"
    return design / "bootstrap_units.5mb.psmcfa_bins.tsv", design / "bootstrap_manifest.tsv"


def prepare_pair(canonical_root: Path, output_root: Path, pair: str) -> dict[str, object]:
    source = primary_path(canonical_root, pair)
    primary = parse_psmcfa(source)
    units = freeze_psmcfa_bootstrap_units(primary, BLOCK_BINS)
    counts = symbol_counts(primary)
    frozen_counts = unit_symbol_counts(primary, units)
    if counts != frozen_counts:
        raise PilotError(f"{pair} frozen units shift the primary PSMCFA population")
    if any(unit.contig not in primary or unit.end > len(primary[unit.contig]) for unit in units):
        raise PilotError(f"{pair} frozen unit crosses a contig boundary")
    manifest = psmcfa_bootstrap_manifest(
        primary, pair, PSMCFA_BOOTSTRAP_DESIGN_SHA256, attempts=ATTEMPTS,
        block_bp=BLOCK_BP, bin_size=BIN_SIZE_BP,
    )
    units_path, manifest_path = design_paths(output_root, pair)
    units_path.parent.mkdir(parents=True, exist_ok=True)
    write_bed(units_path, units)
    with tempfile.NamedTemporaryFile(
        "w", newline="", encoding="utf-8", dir=manifest_path.parent,
        prefix=f".{manifest_path.name}.", delete=False,
    ) as handle:
        partial = Path(handle.name)
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow((
            "replicate", "block_bp", "block_bins", "bin_size", "seed", "unit_count",
            "sampled_unit_indices",
        ))
        for row in manifest:
            writer.writerow((
                row["replicate"], row["block_bp"], row["block_bins"], row["bin_size"],
                row["seed"], row["unit_count"],
                ",".join(map(str, row["sampled_unit_indices"])),
            ))
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(partial, manifest_path)
    return {
        "selection_id": pair,
        "primary_psmcfa": str(source),
        "primary_psmcfa_sha256": sha256_file(source),
        "contigs": len(primary),
        "block_bp": BLOCK_BP,
        "block_bins": BLOCK_BINS,
        "unit_count": len(units),
        "blocks_cross_contig_boundaries": False,
        "primary_psmcfa_bins": sum(counts.values()),
        "frozen_unit_bins": sum(unit.length for unit in units),
        "primary_psmcfa_symbols": dict(sorted(counts.items())),
        "frozen_unit_symbols": dict(sorted(frozen_counts.items())),
        "masked_and_callable_sampling_population_preserved": True,
        "units_path": str(units_path),
        "units_sha256": sha256_file(units_path),
        "manifest_path": str(manifest_path),
        "manifest_sha256": sha256_file(manifest_path),
        "bootstrap_attempts": len(manifest),
    }


def prepare(canonical_root: Path, output_root: Path) -> dict[str, object]:
    if canonical_root != CANONICAL_VGP_ROOT:
        raise PilotError(
            f"repair diagnostics require canonical root {CANONICAL_VGP_ROOT}, got {canonical_root}"
        )
    output_root.mkdir(parents=True, exist_ok=True)
    predeclared = {
        "schema_version": "vgp-psmc-bootstrap-centering-diagnostic-v1",
        "repair_id": REPAIR_ID,
        "canonical_vgp_root": str(canonical_root),
        "diagnostic": dict(CENTERING_DIAGNOSTIC),
    }
    atomic_text(
        output_root / "centering_diagnostic.predeclared.json", canonical_json(predeclared)
    )
    result = {
        "schema_version": "vgp-psmc-bootstrap-repair-design-v1",
        "repair_id": REPAIR_ID,
        "canonical_vgp_root": str(canonical_root),
        "sampling_population": "primary_psmcfa_NKT_bins",
        "bootstrap_design_sha256": PSMCFA_BOOTSTRAP_DESIGN_SHA256,
        "centering_diagnostic": dict(CENTERING_DIAGNOSTIC),
        "pairs": {pair: prepare_pair(canonical_root, output_root, pair) for pair in PAIR_RELATIVE_ROOTS},
    }
    atomic_text(output_root / "repair_design.json", canonical_json(result))
    return result


def read_draw(manifest_path: Path, replicate: int) -> list[int]:
    with manifest_path.open(newline="", encoding="utf-8") as handle:
        matches = [
            row for row in csv.DictReader(handle, delimiter="\t")
            if int(row["replicate"]) == replicate
        ]
    if len(matches) != 1:
        raise PilotError(f"replicate {replicate} does not resolve exactly once")
    row = matches[0]
    if (
        int(row["block_bp"]) != BLOCK_BP
        or int(row["block_bins"]) != BLOCK_BINS
        or int(row["bin_size"]) != BIN_SIZE_BP
    ):
        raise PilotError("bootstrap manifest does not match the frozen repair design")
    sampled = [int(value) for value in row["sampled_unit_indices"].split(",")]
    if len(sampled) != int(row["unit_count"]):
        raise PilotError("bootstrap manifest draw length differs from its unit count")
    return sampled


def run_one(
    canonical_root: Path, output_root: Path, pair: str, replicate: int,
    psmc_binary: Path, scratch_root: Path,
) -> dict[str, object]:
    if replicate < 1 or replicate > ATTEMPTS:
        raise PilotError(f"replicate must be 1..{ATTEMPTS}")
    output = output_root / pair / f"replicate-{replicate:03d}" / "bootstrap.unscaled.psmc"
    metadata = output.with_name("complete.json")
    if output.is_file() and metadata.is_file():
        rows, theta = parse_psmc_unscaled(output)
        if rows and math.isfinite(theta) and all(
            math.isfinite(row["time_2N0"]) and math.isfinite(row["lambda"]) for row in rows
        ):
            return json.loads(metadata.read_text(encoding="utf-8"))
    if not psmc_binary.is_file():
        raise PilotError(f"pinned PSMC binary is absent: {psmc_binary}")
    primary = parse_psmcfa(primary_path(canonical_root, pair))
    units_path, manifest_path = design_paths(output_root, pair)
    units: list[Interval] = []
    with units_path.open(encoding="utf-8") as handle:
        for line in handle:
            contig, start, end = line.rstrip("\n").split("\t")[:3]
            units.append(Interval(contig, int(start), int(end)))
    sampled = read_draw(manifest_path, replicate)
    scratch_root.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
        dir=scratch_root, prefix=f"{REPAIR_ID}-{pair}-{replicate:03d}-"
    ) as temporary:
        temp = Path(temporary)
        bootstrap_input = temp / "bootstrap.psmcfa"
        bootstrap_output = temp / "bootstrap.unscaled.psmc"
        bootstrap_input.write_text(bootstrap_psmcfa(primary, units, sampled), encoding="utf-8")
        subprocess.run(
            [str(psmc_binary), *PSMC_OPTIONS, "-o", str(bootstrap_output), str(bootstrap_input)],
            check=True,
        )
        rows, theta = parse_psmc_unscaled(bootstrap_output)
        finite = math.isfinite(theta) and all(
            math.isfinite(row["time_2N0"]) and math.isfinite(row["lambda"]) for row in rows
        )
        if not rows or not finite:
            raise PilotError(f"{pair} replicate {replicate} emitted a non-finite PSMC fit")
        output.parent.mkdir(parents=True, exist_ok=True)
        staged_output = output.parent / f".{output.name}.{os.getpid()}.partial"
        shutil.copy2(bootstrap_output, staged_output)
        os.replace(staged_output, output)
        result = {
            "schema_version": "vgp-psmc-bootstrap-repair-replicate-v1",
            "repair_id": REPAIR_ID,
            "canonical_vgp_root": str(canonical_root),
            "selection_id": pair,
            "replicate": replicate,
            "sampling_population": "primary_psmcfa_NKT_bins",
            "blocks_cross_contig_boundaries": False,
            "masked_and_callable_sampling_population_preserved": True,
            "finite": True,
            "final_intervals": len(rows),
            "theta_0_per_100bp_bin": theta,
            "output": str(output),
            "output_sha256": sha256_file(output),
            "psmc_binary": str(psmc_binary),
            "psmc_binary_sha256": sha256_file(psmc_binary),
        }
        atomic_text(metadata, canonical_json(result))
        return result


def sampled_population_summary(
    primary: Mapping[str, str], units: Sequence[Interval], manifest_path: Path,
) -> dict[str, object]:
    per_unit = [Counter(primary[unit.contig][unit.start:unit.end]) for unit in units]
    proportions: dict[str, list[float]] = {symbol: [] for symbol in "NKT"}
    bin_counts: list[float] = []
    for replicate in range(1, ATTEMPTS + 1):
        sampled = read_draw(manifest_path, replicate)
        counts = sum((per_unit[index] for index in sampled), Counter())
        total = sum(counts.values())
        bin_counts.append(float(total))
        for symbol in "NKT":
            proportions[symbol].append(counts[symbol] / total)
    return {
        "replicate_bin_count": quantiles(bin_counts),
        "replicate_symbol_proportions": {
            symbol: quantiles(values) for symbol, values in proportions.items()
        },
    }


def audit_pair(canonical_root: Path, output_root: Path, pair: str) -> dict[str, object]:
    primary = parse_psmcfa(primary_path(canonical_root, pair))
    units_path, manifest_path = design_paths(output_root, pair)
    units: list[Interval] = []
    with units_path.open(encoding="utf-8") as handle:
        for line in handle:
            contig, start, end = line.rstrip("\n").split("\t")[:3]
            units.append(Interval(contig, int(start), int(end)))
    counts = symbol_counts(primary)
    frozen_counts = unit_symbol_counts(primary, units)
    if frozen_counts != counts or sum(unit.length for unit in units) != sum(counts.values()):
        raise PilotError(f"{pair} repaired design shifts the primary PSMCFA population")
    primary_rows, primary_theta = parse_psmc_unscaled(primary_psmc_path(canonical_root, pair))
    thetas: list[float] = []
    finite = 0
    interval_counts: Counter[int] = Counter()
    for replicate in range(1, ATTEMPTS + 1):
        path = output_root / pair / f"replicate-{replicate:03d}" / "bootstrap.unscaled.psmc"
        try:
            rows, theta = parse_psmc_unscaled(path)
            valid = math.isfinite(theta) and all(
                math.isfinite(row["time_2N0"]) and math.isfinite(row["lambda"])
                for row in rows
            )
        except (OSError, PilotError):
            valid = False
        if not valid:
            continue
        finite += 1
        thetas.append(theta)
        interval_counts[len(rows)] += 1
    theta_summary = quantiles(thetas)
    centered = theta_summary["q025"] <= primary_theta <= theta_summary["q975"]
    passed = finite == ATTEMPTS and len(thetas) == ATTEMPTS and centered
    return {
        "selection_id": pair,
        "canonical_run_root": str(pair_root(canonical_root, pair)),
        "canonical_primary_psmcfa": str(primary_path(canonical_root, pair)),
        "canonical_primary_psmc": str(primary_psmc_path(canonical_root, pair)),
        "sampling_population": "primary_psmcfa_NKT_bins",
        "block_bp": BLOCK_BP,
        "block_bins": BLOCK_BINS,
        "unit_count": len(units),
        "blocks_cross_contig_boundaries": False,
        "primary_psmcfa_bins": sum(counts.values()),
        "frozen_unit_bins": sum(unit.length for unit in units),
        "primary_psmcfa_symbols": dict(sorted(counts.items())),
        "frozen_unit_symbols": dict(sorted(frozen_counts.items())),
        "masked_and_callable_sampling_population_preserved": frozen_counts == counts,
        "sampled_population_diagnostic": sampled_population_summary(
            primary, units, manifest_path
        ),
        "bootstrap_attempts": ATTEMPTS,
        "finite_bootstraps": finite,
        "bootstrap_interval_counts": {
            str(key): value for key, value in sorted(interval_counts.items())
        },
        "primary": {
            "theta_0_per_100bp_bin": primary_theta,
            "intervals": len(primary_rows),
        },
        "bootstrap_theta_0_per_100bp_bin": theta_summary,
        "centering_diagnostic": {
            **CENTERING_DIAGNOSTIC,
            "observed_lower_bound": theta_summary["q025"],
            "observed_upper_bound": theta_summary["q975"],
            "primary_inside_bounds": centered,
            "passed": passed,
        },
        "passed": passed,
    }


def audit(canonical_root: Path, output_root: Path) -> dict[str, object]:
    predeclared_path = output_root / "centering_diagnostic.predeclared.json"
    predeclared = json.loads(predeclared_path.read_text(encoding="utf-8"))
    if predeclared.get("diagnostic") != dict(CENTERING_DIAGNOSTIC):
        raise PilotError("on-disk centering diagnostic differs from the source predeclaration")
    result = {
        "schema_version": "vgp-psmc-bootstrap-repair-diagnostic-v1",
        "task_id": "repair-vgp-psmc",
        "repair_id": REPAIR_ID,
        "canonical_vgp_root": str(canonical_root),
        "diagnostic_output_root": str(output_root),
        "sampling_population": "primary_psmcfa_NKT_bins",
        "bootstrap_design_sha256": PSMCFA_BOOTSTRAP_DESIGN_SHA256,
        "centering_diagnostic_predeclaration": {
            "path": str(predeclared_path),
            "sha256": sha256_file(predeclared_path),
            "diagnostic": dict(CENTERING_DIAGNOSTIC),
        },
        "pairs": {pair: audit_pair(canonical_root, output_root, pair) for pair in PAIR_RELATIVE_ROOTS},
    }
    result["passed"] = all(pair["passed"] for pair in result["pairs"].values())
    return result


def parser() -> argparse.ArgumentParser:
    value = argparse.ArgumentParser(description=__doc__)
    value.add_argument(
        "--canonical-root", type=Path, default=CANONICAL_VGP_ROOT,
        help=f"must remain the canonical root {CANONICAL_VGP_ROOT}",
    )
    value.add_argument(
        "--output-root", type=Path,
        default=CANONICAL_VGP_ROOT / "pilot/diagnostics" / REPAIR_ID,
    )
    sub = value.add_subparsers(dest="command", required=True)
    sub.add_parser("prepare")
    one = sub.add_parser("run-one")
    one.add_argument("pair", choices=sorted(PAIR_RELATIVE_ROOTS))
    one.add_argument("replicate", type=int)
    one.add_argument("--psmc-binary", type=Path, required=True)
    one.add_argument("--scratch-root", type=Path, required=True)
    check = sub.add_parser("audit")
    check.add_argument("--output", type=Path, required=True)
    return value


def main(argv: Sequence[str] | None = None) -> int:
    args = parser().parse_args(argv)
    try:
        if args.command == "prepare":
            result = prepare(args.canonical_root, args.output_root)
        elif args.command == "run-one":
            result = run_one(
                args.canonical_root, args.output_root, args.pair, args.replicate,
                args.psmc_binary, args.scratch_root,
            )
        else:
            result = audit(args.canonical_root, args.output_root)
            atomic_text(args.output, json.dumps(result, indent=2, sort_keys=True) + "\n")
    except (OSError, PilotError, subprocess.CalledProcessError, json.JSONDecodeError) as error:
        print(f"ERROR: {error}", file=os.sys.stderr)
        return 2
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
