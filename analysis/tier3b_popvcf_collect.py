#!/usr/bin/env python3
"""Select and normalize approved pre-called Tier 3b population VCF pilots.

This collector deliberately has no raw-read mode.  It accepts staged,
checksum-addressable VCF/BCF resources, selects one population and 20 sampling
units deterministically, subsets/normalizes with bcftools, and records enough
provenance to make reruns idempotent.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

if __package__ in (None, ""):  # permit ``python analysis/<script>.py``
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from analysis.tier3_common import Tier3ValidationError, read_fasta, sha256_file
from analysis.tier3b_popvcf_compute import DENOMINATOR_KINDS, read_callable_bed, read_vcf_header


POLICY_ID = "tier3-decisions-v1"
TARGET_UNITS = 20
PILOT_DESIGNS = {
    "dgrp": "inbred_lines_haploidized",
    "ag1000g": "wild_diploid",
}


def _truth(value: Any, default: bool = False) -> bool:
    if value is None or value == "":
        return default
    if isinstance(value, bool):
        return value
    normalized = str(value).strip().lower()
    if normalized in ("1", "true", "yes", "pass", "passed"):
        return True
    if normalized in ("0", "false", "no", "fail", "failed"):
        return False
    raise Tier3ValidationError("unrecognized boolean metadata value {!r}".format(value))


def read_sample_metadata(path: Path) -> List[Dict[str, str]]:
    """Read TSV/CSV sample metadata with stable public identifiers."""

    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        preview = handle.read(4096)
        handle.seek(0)
        try:
            dialect = csv.Sniffer().sniff(preview, delimiters="\t,")
        except csv.Error:
            dialect = csv.excel_tab
        reader = csv.DictReader(handle, dialect=dialect)
        required = {"sample_id", "population_id"}
        if reader.fieldnames is None or not required <= set(reader.fieldnames):
            raise Tier3ValidationError("sample metadata requires sample_id and population_id columns")
        rows = [dict(row) for row in reader]
    sample_ids = [row["sample_id"] for row in rows]
    if any(not value for value in sample_ids) or len(sample_ids) != len(set(sample_ids)):
        raise Tier3ValidationError("sample metadata contains empty or duplicate sample IDs")
    if any(not row["population_id"] for row in rows):
        raise Tier3ValidationError("sample metadata contains an empty population ID")
    return rows


def _qualifies(row: Mapping[str, Any]) -> bool:
    return (
        _truth(row.get("qc_pass"), True)
        and _truth(row.get("unrelated"), True)
        and not _truth(row.get("cross_or_progeny"), False)
        and not _truth(row.get("duplicate_or_related"), False)
        and not _truth(row.get("laboratory_control"), False)
        and not _truth(row.get("contaminated"), False)
    )


def deterministic_sample_selection(
    dataset_id: str,
    rows: Sequence[Mapping[str, Any]],
    population_id: Optional[str] = None,
    target_units: int = TARGET_UNITS,
) -> Dict[str, Any]:
    """Apply the frozen single-population and SHA-256 rank rules."""

    if target_units != TARGET_UNITS:
        raise Tier3ValidationError("primary Tier 3b selection is frozen at 20 sampling units")
    eligible_by_population: Dict[str, List[str]] = {}
    exclusions: List[Tuple[str, str]] = []
    for row in rows:
        sample = str(row["sample_id"])
        population = str(row["population_id"])
        if _qualifies(row):
            eligible_by_population.setdefault(population, []).append(sample)
        else:
            exclusions.append((sample, "provider_or_design_qc"))
    if population_id is None:
        qualifying = [
            (len(samples), population)
            for population, samples in eligible_by_population.items()
            if len(samples) >= target_units
        ]
        if not qualifying:
            raise Tier3ValidationError("no single population has 20 unrelated QC-passing units")
        largest = max(count for count, _population in qualifying)
        population_id = sorted(
            population.encode("utf-8")
            for count, population in qualifying
            if count == largest
        )[0].decode("utf-8")
    if len(eligible_by_population.get(population_id, ())) < target_units:
        raise Tier3ValidationError("requested population lacks 20 unrelated QC-passing units")

    ranking: List[Tuple[str, str]] = []
    for sample in eligible_by_population[population_id]:
        digest = hashlib.sha256(
            (dataset_id + "\0" + population_id + "\0" + sample).encode("utf-8")
        ).hexdigest()
        ranking.append((digest, sample))
    ranking.sort(key=lambda item: (item[0], item[1].encode("utf-8")))
    selected = [sample for _digest, sample in ranking[:target_units]]
    selected_bytes = "".join(sample + "\n" for sample in selected).encode("utf-8")
    excluded_bytes = "".join(
        "{}\t{}\n".format(sample, reason)
        for sample, reason in sorted(exclusions, key=lambda item: item[0].encode("utf-8"))
    ).encode("utf-8")
    return {
        "population_id": population_id,
        "population_eligible_units": len(eligible_by_population[population_id]),
        "eligible_population_counts": {
            key: len(value) for key, value in sorted(eligible_by_population.items())
        },
        "selected_samples": selected,
        "selected_sample_list_sha256": hashlib.sha256(selected_bytes).hexdigest(),
        "exclusion_list_sha256": hashlib.sha256(excluded_bytes).hexdigest(),
        "downsampling_rule": "sha256_dataset_population_sample_take_20",
    }


# Backwards-friendly concise name for callers and evaluator fixtures.
select_samples = deterministic_sample_selection


def _artifact(path: Path) -> Dict[str, Any]:
    path = Path(path).resolve()
    return {
        "path": str(path),
        "sha256": sha256_file(path),
        "size_bytes": path.stat().st_size,
    }


def _run_normalization(
    input_vcf: Path,
    fasta_path: Path,
    selected_path: Path,
    output_bcf: Path,
) -> List[List[str]]:
    bcftools = shutil.which("bcftools")
    if bcftools is None:
        raise Tier3ValidationError("bcftools is required from the pinned pure Guix environment")
    output_bcf.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(str(output_bcf) + ".partial")
    temporary_index = Path(str(temporary) + ".csi")
    commands = [
        [bcftools, "view", "--samples-file", str(selected_path), "--output-type", "u", str(input_vcf)],
        [
            bcftools,
            "norm",
            "--fasta-ref",
            str(fasta_path),
            "--check-ref",
            "e",
            "--multiallelics",
            "+any",
            "--output-type",
            "b",
            "--output",
            str(temporary),
        ],
        [bcftools, "index", "--csi", "--force", str(temporary)],
    ]
    try:
        first = subprocess.Popen(commands[0], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        assert first.stdout is not None
        second = subprocess.Popen(commands[1], stdin=first.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        first.stdout.close()
        _second_stdout, second_stderr = second.communicate()
        first_stderr = first.stderr.read() if first.stderr is not None else b""
        first_code = first.wait()
        if first_code or second.returncode:
            raise Tier3ValidationError(
                "bcftools subset/normalization failed: {} {}".format(
                    first_stderr.decode("utf-8", "replace"), second_stderr.decode("utf-8", "replace")
                ).strip()
            )
        indexed = subprocess.run(commands[2], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if indexed.returncode:
            raise Tier3ValidationError(
                "bcftools indexing failed: {}".format(indexed.stderr.decode("utf-8", "replace").strip())
            )
        os.replace(str(temporary), str(output_bcf))
        os.replace(str(temporary_index), str(Path(str(output_bcf) + ".csi")))
    finally:
        for partial in (temporary, temporary_index):
            if partial.exists():
                partial.unlink()
    return commands


def _existing_idempotent(provenance_path: Path, input_fingerprint: Mapping[str, Any]) -> Optional[Dict[str, Any]]:
    if not provenance_path.is_file():
        return None
    try:
        provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    if provenance.get("input_fingerprint") != input_fingerprint:
        return None
    for key in ("normalized_bcf", "normalized_bcf_index", "selected_sample_list"):
        artifact = provenance.get("outputs", {}).get(key)
        if not artifact:
            return None
        path = Path(artifact["path"])
        if not path.is_file() or _artifact(path) != artifact:
            return None
    provenance["idempotent_reuse"] = True
    return provenance


def collect_population_vcf(
    *,
    dataset_id: str,
    input_vcf: Path,
    fasta_path: Path,
    sample_metadata_path: Path,
    output_dir: Path,
    design: str,
    denominator_kind: str,
    callable_bed_path: Optional[Path] = None,
    population_id: Optional[str] = None,
    assembly_accession: str,
    pilot: Optional[str] = None,
) -> Dict[str, Any]:
    """Collect one approved pre-called pilot, with no raw-read fallback."""

    if pilot is not None:
        if pilot not in PILOT_DESIGNS:
            raise Tier3ValidationError("pilot must be dgrp or ag1000g")
        if design != PILOT_DESIGNS[pilot]:
            raise Tier3ValidationError("pilot/design mismatch")
    if design not in ("wild_diploid", "inbred_lines_haploidized"):
        raise Tier3ValidationError("primary pilot design must be wild diploid or established inbred lines")
    if denominator_kind not in DENOMINATOR_KINDS:
        raise Tier3ValidationError("pre-called resource requires an explicit denominator kind")
    if denominator_kind == "cohort_callable_mask" and callable_bed_path is None:
        raise Tier3ValidationError("cohort_callable_mask requires a staged exact-cohort BED")
    if denominator_kind != "cohort_callable_mask" and callable_bed_path is not None:
        raise Tier3ValidationError("callable BED conflicts with all-sites/gVCF denominator")

    fasta = read_fasta(fasta_path)
    if callable_bed_path:
        read_callable_bed(callable_bed_path, {name: len(sequence) for name, sequence in fasta.items()})
    rows = read_sample_metadata(sample_metadata_path)
    selection = deterministic_sample_selection(dataset_id, rows, population_id)
    # This fast audit rejects wrong dictionaries/REF before invoking bcftools.
    vcf_contigs, vcf_samples = read_vcf_header(input_vcf)
    if vcf_contigs != {name: len(sequence) for name, sequence in fasta.items()}:
        raise Tier3ValidationError("input VCF and exact-reference FASTA dictionaries differ")
    absent = sorted(set(selection["selected_samples"]) - set(vcf_samples))
    if absent:
        raise Tier3ValidationError("selected samples absent from input VCF: {!r}".format(absent))

    input_fingerprint = {
        "policy_id": POLICY_ID,
        "dataset_id": dataset_id,
        "assembly_accession": assembly_accession,
        "design": design,
        "denominator_kind": denominator_kind,
        "input_vcf": _artifact(input_vcf),
        "fasta": _artifact(fasta_path),
        "sample_metadata": _artifact(sample_metadata_path),
        "callable_bed": _artifact(callable_bed_path) if callable_bed_path else None,
        "selection": selection,
    }
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    provenance_path = output_dir / "collect.provenance.json"
    existing = _existing_idempotent(provenance_path, input_fingerprint)
    if existing is not None:
        return existing

    selected_path = output_dir / "selected.samples.txt"
    selected_path.write_text(
        "".join(sample + "\n" for sample in selection["selected_samples"]), encoding="utf-8"
    )
    normalized_bcf = output_dir / "normalized.bcf"
    commands = _run_normalization(input_vcf, fasta_path, selected_path, normalized_bcf)
    provenance = {
        "policy_id": POLICY_ID,
        "collector": "tier3b_popvcf_collect.py",
        "raw_read_mapping_or_joint_calling": "not_implemented_by_policy",
        "pilot": pilot,
        "input_fingerprint": input_fingerprint,
        "selection": selection,
        "denominator": {
            "kind": denominator_kind,
            "invariant_sites_explicit": True,
            "exact_selected_cohort": True,
            "sample_list_sha256": selection["selected_sample_list_sha256"],
        },
        "normalization": {
            "commands": commands,
            "reference_checked": True,
            "multiallelic_representation": "bcftools_norm_plus_any_one_site_record",
            "filtering": "PASS_or_dot_retained_by_compute; failures_excluded_from_denominator",
        },
        "outputs": {
            "selected_sample_list": _artifact(selected_path),
            "normalized_bcf": _artifact(normalized_bcf),
            "normalized_bcf_index": _artifact(Path(str(normalized_bcf) + ".csi")),
        },
        "idempotent_reuse": False,
    }
    temporary = Path(str(provenance_path) + ".partial")
    temporary.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(str(temporary), str(provenance_path))
    return provenance


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset-id", required=True)
    parser.add_argument("--vcf", required=True, type=Path, help="approved staged pre-called VCF/BCF")
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--sample-metadata", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--assembly-accession", required=True, help="assembly accession including version")
    parser.add_argument("--design", required=True, choices=("wild_diploid", "inbred_lines_haploidized"))
    parser.add_argument("--denominator-kind", required=True, choices=sorted(DENOMINATOR_KINDS))
    parser.add_argument("--callable-bed", type=Path)
    parser.add_argument("--population-id")
    parser.add_argument("--pilot", choices=sorted(PILOT_DESIGNS))
    args = parser.parse_args(argv)
    provenance = collect_population_vcf(
        dataset_id=args.dataset_id,
        input_vcf=args.vcf,
        fasta_path=args.fasta,
        sample_metadata_path=args.sample_metadata,
        output_dir=args.output_dir,
        design=args.design,
        denominator_kind=args.denominator_kind,
        callable_bed_path=args.callable_bed,
        population_id=args.population_id,
        assembly_accession=args.assembly_accession,
        pilot=args.pilot,
    )
    print(json.dumps(provenance, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Tier3ValidationError as error:
        raise SystemExit("tier3b collection rejected: {}".format(error))
