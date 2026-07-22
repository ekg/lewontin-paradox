#!/usr/bin/env python3
"""Materialize and independently audit the clean P07 VGP canary."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import subprocess
from collections import Counter, OrderedDict
from pathlib import Path
from typing import Any, Iterable

from analysis.vgp_10_pilot import (
    Interval,
    canonical_json,
    parse_psmc_unscaled,
    parse_fasta,
    read_bed,
    sequence_dictionary,
    sha256_file,
    write_bed,
)
from analysis.vgp_pilot_authorization import write_psmc_scaling_scenarios


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SELECTION = ROOT / "analysis/vgp_clean_canary_selection_v1.json"
AUTHORIZATION = ROOT / "analysis/vgp_pilot_authorization_v2.json"


class CanaryError(RuntimeError):
    """A clean-canary invariant failed."""


def _load(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _pair() -> dict[str, Any]:
    rows = [row for row in _load(AUTHORIZATION)["pairs"] if row["selection_id"] == "P07"]
    if len(rows) != 1:
        raise CanaryError("P07 does not resolve exactly once in the authorization")
    return rows[0]


def materialize(selection_path: Path, data_root: Path) -> dict[str, Any]:
    """Expand only the predeclared BGZF sources into the private run tree."""
    selection = _load(selection_path)
    pair = _pair()
    if selection["selection_id"] != "P07" or selection["run_id"] == "":
        raise CanaryError("unexpected clean-canary selection")
    input_dir = data_root.resolve() / "pilot/inputs/P07"
    input_dir.mkdir(parents=True, exist_ok=False)
    assets: dict[str, Any] = {}
    sequences: dict[str, OrderedDict[str, str]] = {}
    for side in ("h1", "h2"):
        record = selection["immutable_bgzf_inputs"][side]
        source = Path(record["derived_bgzf_path"])
        if not source.is_file():
            raise CanaryError(f"immutable BGZF source absent: {source}")
        if source.read_bytes()[:4] != b"\x1f\x8b\x08\x04":
            raise CanaryError(f"source is not BGZF (missing gzip FEXTRA): {source}")
        destination = input_dir / f"{side}.fa"
        with gzip.open(source, "rb") as incoming, destination.open("wb") as outgoing:
            while True:
                block = incoming.read(16 * 1024 * 1024)
                if not block:
                    break
                outgoing.write(block)
        sequences[side] = parse_fasta(destination)
        observed_dictionary = sequence_dictionary(sequences[side])
        assets[f"{side}_fasta"] = {
            "path": str(destination),
            "sha256": sha256_file(destination),
            "size_bytes": destination.stat().st_size,
            "source_bgzf_path": str(source),
            "source_bgzf_sha256": sha256_file(source),
            "source_cas_sha256": record["source_cas_sha256"],
            "sequence_dictionary": observed_dictionary,
        }
    old = Path(
        "/moosefs/erikg/vgp/pilot/outputs/vgp10-auth-20260718-v2/P07/core/preflight/preflight.json"
    )
    if old.is_file():
        expected = _load(old)["dictionaries"]
        for side in ("h1", "h2"):
            if assets[f"{side}_fasta"]["sequence_dictionary"] != expected[f"{side}_fasta"]:
                raise CanaryError(f"{side} BGZF dictionary differs from the frozen P07 result")
    h1_universe = [Interval(name, 0, len(seq)) for name, seq in sequences["h1"].items()]
    query_universe = [Interval(name, 0, len(seq)) for name, seq in sequences["h2"].items()]
    write_bed(input_dir / "h1_universe.bed", h1_universe)
    write_bed(input_dir / "eligible_query_regions.bed", h1_universe + query_universe)
    (input_dir / "exclusions").mkdir()
    annotation = dict(selection["annotation"])
    annotation.update(
        {
            "canonical_vgp_root": "/moosefs/erikg/vgp",
            "binding_status": "exact_native_sequence_dictionary_equal",
            "sequence_dictionary": assets["h1_fasta"]["sequence_dictionary"],
        }
    )
    manifest = {
        "canonical_vgp_root": str(data_root.resolve()),
        "selection_id": "P07",
        "biosample": pair["biosample"],
        "individual_or_isolate": pair["individual_or_isolate"],
        "h1_accession_version": pair["h1_accession_version"],
        "h2_accession_version": pair["h2_accession_version"],
        "orientation": "H1_reference_H2_query",
        "assets": assets,
        "confidence_covariates": {
            name: None for name in pair["missing_confidence_covariates"]
        },
        "selective_validation": {},
        "annotation": annotation,
        "result_gates": {"minimum_callable_bp": 100_000_000, "minimum_callable_fraction": 0.60},
        "authorization_id": "vgp10-auth-20260718-v2",
        "clean_run_id": selection["run_id"],
        "prior_intermediates_reused": False,
    }
    (input_dir / "input-manifest.json").write_text(canonical_json(manifest), encoding="utf-8")
    (input_dir / "resources.json").write_text(canonical_json(pair["resources"]), encoding="utf-8")
    write_psmc_scaling_scenarios(input_dir / "psmc_scaling_scenarios.tsv")
    source_manifest = {
        "schema_version": "vgp-clean-canary-staged-input-v1",
        "selection_id": "P07",
        "run_id": selection["run_id"],
        "assets": assets,
        "dictionaries_match_frozen_P07": True,
        "materialized_inside_private_scratch": True,
    }
    (input_dir / "clean_source_manifest.json").write_text(canonical_json(source_manifest), encoding="utf-8")
    return source_manifest


def _interval_overlap(rows: Iterable[Interval], contig: str, start: int, end: int) -> int:
    return sum(max(0, min(end, row.end) - max(start, row.start)) for row in rows if row.contig == contig)


def _fasta_records(path: Path) -> OrderedDict[str, str]:
    return parse_fasta(path)


def _psmcfa_records(path: Path) -> OrderedDict[str, str]:
    records: OrderedDict[str, str] = OrderedDict()
    current = None
    chunks: list[str] = []
    with path.open(encoding="ascii") as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if current is not None:
                    records[current] = "".join(chunks)
                current, chunks = line[1:].split()[0], []
            elif line:
                chunks.append(line)
    if current is not None:
        records[current] = "".join(chunks)
    return records


def _strata(dictionary: list[dict[str, Any]], width: int = 5_000_000) -> list[dict[str, Any]]:
    total = sum(int(row["length"]) for row in dictionary)
    points = [("early", 0), ("middle", total // 2), ("late", max(0, total - width))]
    result = []
    for label, offset in points:
        cumulative = 0
        for row in dictionary:
            length = int(row["length"])
            if cumulative + length > offset:
                start = offset - cumulative
                end = min(length, start + width)
                result.append({"stratum": label, "contig": row["name"], "start": start, "end": end})
                break
            cumulative += length
    return result


def _query_variant_lines(bcftools: Path, path: Path, region: str) -> list[str]:
    command = [str(bcftools), "query", "-r", region, "-f", "%CHROM\\t%POS\\t%REF\\t%ALT\\t[%GT]\\n", str(path)]
    return subprocess.run(command, check=True, text=True, capture_output=True).stdout.splitlines()


def _annotation_counts(gff: Path, strata: list[dict[str, Any]]) -> dict[str, int]:
    result = {row["stratum"]: 0 for row in strata}
    opener = gzip.open if gff.read_bytes()[:2] == b"\x1f\x8b" else open
    with opener(gff, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            start, end = int(fields[3]) - 1, int(fields[4])
            for row in strata:
                if fields[0] == row["contig"] and start < row["end"] and end > row["start"]:
                    result[row["stratum"]] += 1
    return result


def _nearest_quantile(values: list[float], fraction: float) -> float:
    """Return the predeclared PSMC diagnostic's nearest-index quantile."""
    if not values:
        raise CanaryError("cannot summarize an empty bootstrap distribution")
    ordered = sorted(values)
    return ordered[round((len(ordered) - 1) * fraction)]


def _bootstrap_audit(run_root: Path) -> dict[str, Any]:
    """Independently verify source population, finite fits, and centering."""
    primary_psmcfa = _psmcfa_records(run_root / "consensus/consensus/input.psmcfa")
    primary_symbols = Counter("".join(primary_psmcfa.values()))
    if set(primary_symbols) - {"N", "K", "T"}:
        raise CanaryError(f"unexpected primary PSMCFA symbols: {sorted(primary_symbols)}")

    units: list[Interval] = []
    with (run_root / "consensus/consensus/bootstrap_units.5mb.psmcfa_bins.tsv").open(
        encoding="utf-8"
    ) as handle:
        for line in handle:
            contig, start, end = line.rstrip("\n").split("\t")[:3]
            units.append(Interval(contig, int(start), int(end)))
    frozen_symbols: Counter[str] = Counter()
    for unit in units:
        frozen_symbols.update(primary_psmcfa[unit.contig][unit.start:unit.end])
    if frozen_symbols != primary_symbols:
        raise CanaryError("5 Mb bootstrap units do not preserve the primary N/K/T population")

    manifest_path = run_root / "consensus/consensus/bootstrap_manifest.tsv"
    with manifest_path.open(newline="", encoding="utf-8") as handle:
        manifest = list(csv.DictReader(handle, delimiter="\t"))
    if [int(row["replicate"]) for row in manifest] != list(range(1, 201)):
        raise CanaryError("bootstrap manifest is not the exact ordered set 1..200")
    for row in manifest:
        sampled = [int(value) for value in row["sampled_unit_indices"].split(",")]
        if int(row["unit_count"]) != len(units) or len(sampled) != len(units):
            raise CanaryError("bootstrap draw does not preserve the frozen unit count")
        if any(index < 0 or index >= len(units) for index in sampled):
            raise CanaryError("bootstrap draw contains an out-of-range unit")

    primary_rows, primary_theta = parse_psmc_unscaled(
        run_root / "psmc/replicate-000/unscaled.psmc"
    )
    thetas: list[float] = []
    interval_counts: Counter[int] = Counter()
    for replicate in range(1, 201):
        rows, theta = parse_psmc_unscaled(
            run_root / f"psmc/replicate-{replicate:03d}/bootstrap.unscaled.psmc"
        )
        if not math.isfinite(theta) or not all(
            math.isfinite(row["time_2N0"]) and math.isfinite(row["lambda"])
            for row in rows
        ):
            raise CanaryError(f"PSMC replicate {replicate} is not finite")
        thetas.append(theta)
        interval_counts[len(rows)] += 1
    lower = _nearest_quantile(thetas, 0.025)
    upper = _nearest_quantile(thetas, 0.975)
    centered = lower <= primary_theta <= upper
    if not centered:
        raise CanaryError(
            f"primary theta {primary_theta} is outside predeclared bootstrap interval {lower}..{upper}"
        )
    return {
        "sampling_population": "primary_psmcfa_NKT_bins",
        "block_bp": 5_000_000,
        "block_bins": 50_000,
        "unit_count": len(units),
        "blocks_cross_contig_boundaries": False,
        "primary_psmcfa_bins": sum(primary_symbols.values()),
        "frozen_unit_bins": sum(frozen_symbols.values()),
        "primary_psmcfa_symbols": dict(sorted(primary_symbols.items())),
        "frozen_unit_symbols": dict(sorted(frozen_symbols.items())),
        "masked_and_callable_sampling_population_preserved": True,
        "manifest_sha256": sha256_file(manifest_path),
        "bootstrap_attempts": 200,
        "finite_bootstraps": len(thetas),
        "bootstrap_interval_counts": {
            str(key): value for key, value in sorted(interval_counts.items())
        },
        "primary": {"theta_0_per_100bp_bin": primary_theta, "intervals": len(primary_rows)},
        "bootstrap_theta_0_per_100bp_bin": {
            "minimum": min(thetas),
            "q025": lower,
            "median": _nearest_quantile(thetas, 0.5),
            "q975": upper,
            "maximum": max(thetas),
        },
        "centering_diagnostic": {
            "name": "primary_theta_in_equal_tail_central_95pct_bootstrap_interval",
            "metric": "final_native_iteration_theta_0_per_100bp_bin",
            "lower_quantile": 0.025,
            "upper_quantile": 0.975,
            "quantile_method": "nearest index round((n - 1) * q)",
            "required_finite_outputs": 200,
            "required_attempts": 200,
            "predeclared_before_execution": True,
            "observed_lower_bound": lower,
            "observed_upper_bound": upper,
            "primary_inside_bounds": True,
            "passed": True,
        },
        "passed": True,
    }


def audit(selection_path: Path, run_root: Path, bcftools: Path, output: Path) -> dict[str, Any]:
    selection = _load(selection_path)
    required = ["preflight", "mapping", "impg", "variants", "consensus", "annotation", "psmc/finalize"]
    for stage in required:
        if not (run_root / stage / ".complete.json").is_file():
            raise CanaryError(f"stage sentinel absent: {stage}")
    for replicate in range(201):
        if not (run_root / f"psmc/replicate-{replicate:03d}/.complete.json").is_file():
            raise CanaryError(f"PSMC replicate sentinel absent: {replicate}")
    inputs = run_root.parent.parent.parent / "inputs/P07"
    if not inputs.is_dir():
        # The worker passes the local run root; its input tree is sibling to runs.
        inputs = run_root.parents[3] / "inputs/P07"
    manifest = _load(inputs / "input-manifest.json")
    join = _load(run_root / "consensus/join_qc.json")
    callable_bp = int(join["consensus"]["consensus_callable_bp"])
    snps = int(join["consensus"]["heterozygous_snps"])
    pi = snps / callable_bp
    tolerance = selection["comparison_target"]
    comparison = {
        "observed_pi": pi,
        "target_pi": tolerance["pi"],
        "pi_absolute_difference": abs(pi - tolerance["pi"]),
        "pi_within_tolerance": abs(pi - tolerance["pi"]) <= tolerance["pi_absolute_tolerance"],
        "observed_callable_bp": callable_bp,
        "target_callable_bp": tolerance["callable_bp"],
        "callable_bp_difference": callable_bp - tolerance["callable_bp"],
        "callable_within_tolerance": abs(callable_bp - tolerance["callable_bp"]) <= tolerance["callable_bp_absolute_tolerance"],
    }
    if not comparison["pi_within_tolerance"] or not comparison["callable_within_tolerance"]:
        raise CanaryError(f"clean result is outside predeclared tolerance: {comparison}")

    h1 = _fasta_records(inputs / "h1.fa")
    h2 = _fasta_records(inputs / "h2.fa")
    digest_rows = []
    lookup: dict[str, str] = {}
    for side, records in (("H1", h1), ("H2", h2)):
        for name, sequence in records.items():
            digest = hashlib.sha256(sequence.encode("ascii")).hexdigest()
            lookup[name] = digest
            digest_rows.append({"side": side, "sequence_id": name, "length": len(sequence), "sequence_sha256": digest})
    observed_ids: set[str] = set()
    with (run_root / "mapping/h2_to_h1.1to1.paf").open(encoding="ascii") as handle:
        for line in handle:
            fields = line.split("\t")
            observed_ids.update((fields[0], fields[5]))
    with (run_root / "impg/focus.native.bed").open(encoding="ascii") as handle:
        for line in handle:
            if line.strip() and not line.startswith("#"):
                observed_ids.add(line.split("\t", 1)[0])
    unresolved = sorted(observed_ids - lookup.keys())
    if unresolved:
        raise CanaryError(f"graph sequence identifiers lack sequence digests: {unresolved[:5]}")
    ledger = output.parent / "vgp_clean_canary_graph_sequence_digests.tsv"
    with ledger.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=("side", "sequence_id", "length", "sequence_sha256"), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(digest_rows)

    strata = _strata(manifest["assets"]["h1_fasta"]["sequence_dictionary"])
    callable_rows = read_bed(run_root / "consensus/masks/callable.bed")
    consensus = _fasta_records(run_root / "consensus/consensus/consensus.fa")
    psmcfa = _psmcfa_records(run_root / "consensus/consensus/input.psmcfa")
    annotation_counts = _annotation_counts(Path(manifest["annotation"]["gff_path"]), strata)
    validations = []
    vcf = run_root / "variants/normalized.vcf.gz"
    bcf = run_root / "variants/normalized.bcf"
    for row in strata:
        region = f"{row['contig']}:{row['start'] + 1}-{row['end']}"
        vcf_lines = _query_variant_lines(bcftools, vcf, region)
        bcf_lines = _query_variant_lines(bcftools, bcf, region)
        if vcf_lines != bcf_lines:
            raise CanaryError(f"VCF/BCF mismatch in {row['stratum']}")
        callable_in_region = _interval_overlap(callable_rows, row["contig"], row["start"], row["end"])
        sequence = consensus[row["contig"]][row["start"]:row["end"]]
        psmc_slice = psmcfa[row["contig"]][row["start"] // 100:math.ceil(row["end"] / 100)]
        validations.append({
            **row,
            "region": region,
            "vcf_bcf_equal": True,
            "variant_records": len(vcf_lines),
            "callable_bp": callable_in_region,
            "consensus_non_N_bp": sum(base != "N" for base in sequence),
            "annotation_feature_overlaps": annotation_counts[row["stratum"]],
            "psmc_K_bins": psmc_slice.count("K"),
            "psmc_T_bins": psmc_slice.count("T"),
            "psmc_N_bins": psmc_slice.count("N"),
        })
    psmc_qc = _load(run_root / "psmc/finalize/psmc_qc.json")
    if psmc_qc["finite_bootstraps"] != 200:
        raise CanaryError("clean canary did not produce 200 finite bootstraps")
    bootstrap_audit = _bootstrap_audit(run_root)
    annotation = _load(run_root / "annotation/exact_partitions.json")
    result = {
        "schema_version": "vgp-clean-canary-execution-v1",
        "task_id": "run-vgp-clean-canary",
        "selection_id": "P07",
        "run_id": selection["run_id"],
        "prior_intermediates_reused": False,
        "comparison": comparison,
        "diversity": {"heterozygous_snps": snps, "callable_bp": callable_bp, "pi": pi},
        "mapping": {
            "required_option": "--num-mappings 1:1",
            "paf_sha256": sha256_file(run_root / "mapping/h2_to_h1.1to1.paf"),
            "fastga_scratch_contract": _load(run_root / "mapping/fastga_scratch_contract.json"),
        },
        "graph_sequence_digest_ledger": {"path": str(ledger), "sha256": sha256_file(ledger), "resolved_ids": len(observed_ids), "unresolved_ids": 0},
        "independent_early_middle_late_validation": validations,
        "psmc": {**psmc_qc, "independent_bootstrap_audit": bootstrap_audit},
        "annotation": annotation,
        "source_manifest": _load(inputs / "clean_source_manifest.json"),
    }
    output.write_text(canonical_json(result), encoding="utf-8")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    material = sub.add_parser("materialize")
    material.add_argument("--selection", type=Path, default=DEFAULT_SELECTION)
    material.add_argument("--data-root", type=Path, required=True)
    check = sub.add_parser("audit")
    check.add_argument("--selection", type=Path, default=DEFAULT_SELECTION)
    check.add_argument("--run-root", type=Path, required=True)
    check.add_argument("--bcftools", type=Path, required=True)
    check.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.command == "materialize":
        print(canonical_json(materialize(args.selection, args.data_root)), end="")
    else:
        print(canonical_json(audit(args.selection, args.run_root, args.bcftools, args.output)), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
