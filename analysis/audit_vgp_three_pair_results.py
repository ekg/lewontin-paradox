#!/usr/bin/env python3
"""Independently audit the three completed VGP reliability-pilot pairs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path
from typing import Iterable, Mapping, Sequence

from analysis.vgp_10_pilot import (
    Interval,
    PilotError,
    canonical_json,
    interval_bp,
    intersect_intervals,
    parse_psmc_unscaled,
    read_bed,
    sha256_file,
)


def fasta_lengths(path: Path | str) -> list[tuple[str, int]]:
    rows: list[tuple[str, int]] = []
    name: str | None = None
    length = 0
    with Path(path).open(encoding="utf-8") as handle:
        for number, line in enumerate(handle, 1):
            if line.startswith(">"):
                if name is not None:
                    rows.append((name, length))
                name, length = line[1:].split()[0], 0
                if not name:
                    raise PilotError(f"empty FASTA identifier at line {number}")
            elif name is None:
                raise PilotError("sequence precedes FASTA header")
            else:
                length += len(line.strip())
    if name is not None:
        rows.append((name, length))
    if not rows or len({name for name, _ in rows}) != len(rows):
        raise PilotError("empty or duplicate FASTA dictionary")
    return rows


def select_strata(lengths: Sequence[tuple[str, int]], width: int = 5_000_000) -> list[Interval]:
    eligible = [(name, length) for name, length in lengths if length >= width]
    if len(eligible) < 3:
        raise PilotError("fewer than three sequences support stratified re-audit")
    first, middle, last = eligible[0], eligible[len(eligible) // 2], eligible[-1]
    return [
        Interval(first[0], 0, width),
        Interval(middle[0], (middle[1] - width) // 2, (middle[1] + width) // 2),
        Interval(last[0], last[1] - width, last[1]),
    ]


def _region_index(regions: Iterable[Interval]) -> dict[str, list[Interval]]:
    result: dict[str, list[Interval]] = {}
    for region in regions:
        result.setdefault(region.contig, []).append(region)
    return result


def fasta_non_n_counts(path: Path | str, regions: Sequence[Interval]) -> dict[Interval, int]:
    wanted = _region_index(regions)
    counts = {region: 0 for region in regions}
    name: str | None = None
    position = 0
    with Path(path).open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                name, position = line[1:].split()[0], 0
                continue
            sequence = line.strip().upper()
            for region in wanted.get(name or "", []):
                start, end = max(position, region.start), min(position + len(sequence), region.end)
                if start < end:
                    piece = sequence[start - position:end - position]
                    counts[region] += sum(base != "N" for base in piece)
            position += len(sequence)
    return counts


def psmc_population_counts(path: Path | str, regions: Sequence[Interval]) -> dict[Interval, Counter]:
    wanted = _region_index(regions)
    counts = {region: Counter() for region in regions}
    name: str | None = None
    bin_position = 0
    with Path(path).open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                name, bin_position = line[1:].split()[0], 0
                continue
            symbols = line.strip().upper()
            for region in wanted.get(name or "", []):
                bin_start, bin_end = region.start // 100, (region.end + 99) // 100
                start, end = max(bin_position, bin_start), min(bin_position + len(symbols), bin_end)
                if start < end:
                    counts[region].update(symbols[start - bin_position:end - bin_position])
            bin_position += len(symbols)
    return counts


def variant_counts(path: Path | str, regions: Sequence[Interval]) -> dict[Interval, dict[str, int]]:
    wanted = _region_index(regions)
    result = {region: {"records": 0, "snps": 0} for region in regions}
    with Path(path).open(encoding="utf-8") as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                raise PilotError(f"VCF row {number} is malformed")
            contig, pos0, ref, alt = fields[0], int(fields[1]) - 1, fields[3], fields[4]
            for region in wanted.get(contig, []):
                if region.start <= pos0 < region.end:
                    result[region]["records"] += 1
                    result[region]["snps"] += int(
                        len(ref) == len(alt) == 1 and ref in "ACGT" and alt in "ACGT"
                    )
    return result


def audit_closed_stage(path: Path | str) -> dict[str, object]:
    root = Path(path)
    sentinel = json.loads((root / ".complete.json").read_text(encoding="utf-8"))
    files = sentinel.get("files")
    if not isinstance(files, Mapping) or not files:
        raise PilotError(f"stage sentinel has no file ledger: {root}")
    for relative, expected in files.items():
        candidate = root / relative
        if not candidate.is_file() or sha256_file(candidate) != expected:
            raise PilotError(f"closed stage digest mismatch: {candidate}")
    return {"path": str(root), "verified_files": len(files), "sentinel_sha256": sha256_file(root / ".complete.json")}


def _nearest_interval(values: Sequence[float]) -> tuple[float, float]:
    ordered = sorted(values)
    return ordered[round(0.025 * (len(ordered) - 1))], ordered[round(0.975 * (len(ordered) - 1))]


def _bootstrap_audit(root: Path) -> dict[str, object]:
    primary_rows, primary_theta = parse_psmc_unscaled(root / "psmc/replicate-000/unscaled.psmc")
    theta: list[float] = []
    finite = 0
    for replicate in range(1, 201):
        rows, value = parse_psmc_unscaled(
            root / f"psmc/replicate-{replicate:03d}/bootstrap.unscaled.psmc"
        )
        finite += int(bool(rows))
        theta.append(value)
    low, high = _nearest_interval(theta)
    if finite != 200 or not low <= primary_theta <= high:
        raise PilotError("PSMC bootstrap population is not 200/200 finite and centered")
    qc = json.loads((root / "psmc/finalize/psmc_qc.json").read_text())
    if qc["bootstrap_attempts"] != 200 or qc["finite_bootstraps"] != 200:
        raise PilotError("PSMC finalize census disagrees with independent audit")
    return {
        "primary_theta": primary_theta,
        "nearest_index_central_95pct": [low, high],
        "primary_theta_centered": True,
        "finite_bootstraps": finite,
        "bootstrap_attempts": 200,
        "primary_intervals": len(primary_rows),
    }


def _telemetry(root: Path) -> list[dict[str, object]]:
    values = []
    for path in sorted((root / "telemetry").glob("*.json")):
        value = json.loads(path.read_text())
        value["record"] = path.name
        values.append(value)
    return values


def audit_pair(selection: Mapping[str, object], root: Path, p07_execution: Mapping[str, object] | None = None) -> dict[str, object]:
    selection_id = str(selection["selection_id"])
    required = ["mapping", "impg", "variants", "consensus", "psmc/finalize"]
    closed = [audit_closed_stage(root / stage) for stage in required]
    replicate_ledgers = [audit_closed_stage(root / f"psmc/replicate-{number:03d}") for number in range(201)]

    join = json.loads((root / "consensus/join_qc.json").read_text())
    reconstruction = json.loads((root / "variants/paf_variant_audit.json").read_text())
    if reconstruction["reconstruction_failures"] != 0:
        raise PilotError(f"{selection_id}: REF/ALT reconstruction failure")
    if join["concordance"]["h2_reconstruction_failures"] != 0:
        raise PilotError(f"{selection_id}: normalized reconstruction failure")
    mask = join["mask"]
    if mask["accounting_discrepancy_bp"] != 0:
        raise PilotError(f"{selection_id}: mask accounting discrepancy")
    consensus = join["consensus"]
    callable_bp, heterozygous_snps = consensus["consensus_callable_bp"], consensus["heterozygous_snps"]
    if callable_bp <= 0 or heterozygous_snps <= 0:
        raise PilotError(f"{selection_id}: absent core biological result")

    fasta = root / "consensus/consensus/consensus.fa"
    regions = select_strata(fasta_lengths(fasta))
    variants = variant_counts(root / "consensus/normalized.vcf", regions)
    non_n = fasta_non_n_counts(fasta, regions)
    populations = psmc_population_counts(root / "consensus/consensus/input.psmcfa", regions)
    callable_rows = read_bed(root / "consensus/masks/callable.bed")
    strata = []
    for label, region in zip(("early", "middle", "late"), regions):
        pop = populations[region]
        if set(pop) - set("NKT"):
            raise PilotError(f"{selection_id}: invalid PSMC population symbol")
        strata.append({
            "stratum": label, "contig": region.contig, "start": region.start, "end": region.end,
            "callable_bp": interval_bp(intersect_intervals([region], callable_rows)),
            "consensus_non_N_bp": non_n[region], **variants[region],
            "psmc_population_bins": {symbol: pop[symbol] for symbol in "NKT"},
        })

    scenario_path = root / "psmc/finalize/scenario_scaled_trajectories.tsv"
    with scenario_path.open(newline="", encoding="utf-8") as handle:
        scenarios = {row["scenario_id"] for row in csv.DictReader(handle, delimiter="\t")}
    graph: Mapping[str, object]
    if selection_id == "P07":
        if p07_execution is None:
            raise PilotError("P07 execution record is required")
        graph = p07_execution["graph_sequence_digest_ledger"]  # type: ignore[index]
        annotation = p07_execution["annotation"]  # type: ignore[index]
        if annotation["annotation_status"] != "exact_native":  # type: ignore[index]
            raise PilotError("P07 exact annotation is absent")
    else:
        graph = json.loads((root / "impg/graph_identifier_audit.json").read_text())
        if graph["unresolved_ids"] != 0 or graph["silently_omitted_regions"] != 0:
            raise PilotError(f"{selection_id}: graph identifier audit failed")
        annotation = {"annotation_status": "missing_nonblocking", "core_outputs_stopped": False}

    mapping_records = 0
    strands = Counter()
    with (root / "mapping/h2_to_h1.1to1.paf").open() as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12 or fields[4] not in {"+", "-"}:
                raise PilotError(f"{selection_id}: malformed PAF row {number}")
            qlen, qs, qe, tlen, ts, te = map(int, (fields[1], fields[2], fields[3], fields[6], fields[7], fields[8]))
            if not (0 <= qs < qe <= qlen and 0 <= ts < te <= tlen):
                raise PilotError(f"{selection_id}: PAF coordinate failure at row {number}")
            strands[fields[4]] += 1
            mapping_records += 1
    if mapping_records == 0:
        raise PilotError(f"{selection_id}: no exact mapping rows")

    normalized = root / "variants/normalized.vcf.gz"
    bcf = root / "variants/normalized.bcf"
    result = {
        "selection_id": selection_id, "species": selection["species"],
        "individual_or_isolate": selection["individual_or_isolate"],
        "failure_class": selection["failure_class"], "actual_core_biological_result": True,
        "diversity": {"heterozygous_snps": heterozygous_snps, "callable_bp": callable_bp,
                      "pi": heterozygous_snps / callable_bp},
        "normalized_variants": {"vcf_gz_sha256": sha256_file(normalized), "bcf_sha256": sha256_file(bcf),
                                "vcf_index_present": (Path(str(normalized) + ".tbi")).is_file(),
                                "bcf_index_present": (Path(str(bcf) + ".csi")).is_file()},
        "consensus_sha256": sha256_file(fasta), "mask": mask,
        "coordinate_and_strand_audit": {"records": mapping_records, "strand_counts": dict(strands),
                                        "invalid_coordinates": 0},
        "ref_alt_reconstruction": reconstruction,
        "normalized_concordance": join["concordance"], "graph_identifier_audit": graph,
        "psmc": _bootstrap_audit(root), "scenario_count": len(scenarios),
        "psmc_population_preserved": join["masked_and_callable_population_preserved"],
        "independent_stratified_reaudit": strata, "annotation": annotation,
        "closed_stage_ledgers": closed, "closed_psmc_replicate_ledgers": replicate_ledgers,
        "resource_telemetry": _telemetry(root),
    }
    if not all(result["normalized_variants"][key] for key in ("vcf_index_present", "bcf_index_present")):
        raise PilotError(f"{selection_id}: normalized index absent")
    return result


def audit_results(selection_path: Path, output_path: Path) -> dict[str, object]:
    freeze = json.loads(selection_path.read_text())
    by_id = {row["selection_id"]: row for row in freeze["selections"]}
    p07_execution = json.loads(Path("analysis/vgp_clean_canary_execution_v1.json").read_text())
    roots = {
        "P07": Path("/moosefs/erikg/vgp/pilot/clean-canary/vgp-clean-canary-20260722-v1/P07"),
        "P02": Path("/moosefs/erikg/vgp/pilot/runs/vgp-three-pair-20260722-v1/P02"),
        "P03": Path("/moosefs/erikg/vgp/pilot/runs/vgp-three-pair-20260722-v1/P03"),
    }
    pairs = [audit_pair(by_id[key], roots[key], p07_execution if key == "P07" else None)
             for key in ("P07", "P03", "P02")]
    comparison_path = Path(
        "/moosefs/erikg/vgp/pilot/three-pair/vgp-three-pair-20260722-v1/"
        "P03/backend-comparison/comparison/comparison.json"
    )
    comparison = json.loads(comparison_path.read_text())
    if comparison["overlapping_target_bp"] <= 0:
        raise PilotError("controlled backend comparison is absent")
    result = {
        "schema_version": "vgp-three-pair-execution-v1", "task_id": "run-vgp-three-pair",
        "run_id": freeze["run_id"], "selection_freeze_sha256": sha256_file(selection_path),
        "actual_core_biological_results": len(pairs), "completion_gate_passed": len(pairs) == 3,
        "pairs": pairs, "controlled_fastga_wfmash_comparison": comparison,
        "remaining_pipeline_limitations": [
            "FastGA remains unreliable for P03 at whole-assembly scale; pinned WFMASH is the recorded infrastructure fallback.",
            "Backend-specific gap placement lowers exact raw-variant concordance despite high shared target coverage; normalized biological estimates retain backend provenance.",
            "P02 and P03 have no cataloged exact native annotation, so only their annotation outputs are absent.",
            "Assembly-derived diversity is a same-individual haplotype comparison and is not a population sample-size substitute.",
        ],
        "limitations_are_not_biological_exclusions": True,
    }
    output_path.write_text(canonical_json(result), encoding="utf-8")
    return result


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selection", type=Path, default=Path("analysis/vgp_three_pair_selection_v1.json"))
    parser.add_argument("--output", type=Path, default=Path("analysis/vgp_three_pair_execution_v1.json"))
    args = parser.parse_args(argv)
    try:
        result = audit_results(args.selection, args.output)
    except (PilotError, OSError, KeyError, ValueError, json.JSONDecodeError) as error:
        print(f"ERROR: {error}", file=__import__("sys").stderr)
        return 2
    print(json.dumps({"actual_core_biological_results": result["actual_core_biological_results"],
                      "completion_gate_passed": result["completion_gate_passed"]}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
