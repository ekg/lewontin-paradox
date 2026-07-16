#!/usr/bin/env python3
"""Run, collect, and validate recovered biological Tier 3B population tuples."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import re
import subprocess
import tempfile
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Sequence, Tuple

if __package__ in (None, ""):
    import sys

    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from analysis.tier3_common import Tier3ValidationError, sha256_file
from analysis.tier3b_popvcf_compute import compute_population_pi


REQUIRED_INPUTS = {
    "reference_path": "reference_sha256",
    "reference_fai_path": "reference_fai_sha256",
    "native_annotation_path": "native_annotation_sha256",
    "callset_path": "callset_sha256",
    "callset_index_path": "callset_index_sha256",
    "sample_list_path": "sample_list_sha256",
    "sample_metadata_path": "sample_metadata_sha256",
    "callable_mask_path": "callable_mask_sha256",
    "callable_source_path": "callable_source_sha256",
}
STATISTICS = ("population_pi", "pi_S", "pi_W", "pi_S_over_pi_W")
DIVERSITY_FIELDS = (
    "tuple_id",
    "biological",
    "scientific_name",
    "taxon_id",
    "population_id",
    "population_release",
    "design",
    "source_modality",
    "region",
    "reference_accession",
    "native_annotation_accession",
    "eligible_sample_size",
    "nominal_chromosomes",
    "statistic",
    "annotation_category",
    "variant_count",
    "numerator",
    "numerator_definition",
    "callable_site_denominator",
    "estimate",
    "component_s_numerator",
    "component_s_callable_sites",
    "component_w_numerator",
    "component_w_callable_sites",
    "uncertainty_method",
    "uncertainty_unit",
    "uncertainty_replicates",
    "uncertainty_standard_error",
    "interval_low",
    "interval_high",
    "interval_type",
    "frozen_genomic_bootstrap_status",
    "exclusions",
    "raw_result_path",
)


def _read_manifest(path: Path) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise Tier3ValidationError("acquisition manifest has no rows")
    approved = [
        row
        for row in rows
        if row.get("status") == "approved" and row.get("biological", "").lower() == "true"
    ]
    if len(approved) < 2:
        raise Tier3ValidationError("fewer than two approved biological acquisition tuples")
    tuple_ids = [row["tuple_id"] for row in approved]
    if len(tuple_ids) != len(set(tuple_ids)):
        raise Tier3ValidationError("duplicate approved tuple IDs")
    return approved


def _verify_inputs(row: Mapping[str, str]) -> None:
    for path_field, checksum_field in REQUIRED_INPUTS.items():
        path = Path(row[path_field])
        if not path.is_file():
            raise Tier3ValidationError("missing {} for {}: {}".format(path_field, row["tuple_id"], path))
        observed = sha256_file(path)
        if observed != row[checksum_field]:
            raise Tier3ValidationError(
                "{} checksum mismatch for {}: {} != {}".format(
                    path_field, row["tuple_id"], observed, row[checksum_field]
                )
            )
    samples = [
        line.strip()
        for line in Path(row["sample_list_path"]).read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if len(samples) != int(row["sample_count"]) or len(samples) != 20:
        raise Tier3ValidationError("approved tuple does not contain exactly 20 selected samples")
    if row["callable_coordinate_system"] != "0-based half-open BED":
        raise Tier3ValidationError("unexpected callable coordinate convention")
    if "1-based" not in row["vcf_coordinate_system"] or "1-based" not in row["annotation_coordinate_system"]:
        raise Tier3ValidationError("unexpected VCF/GFF coordinate convention")
    if not row["contig_compatibility"].startswith("PASS:"):
        raise Tier3ValidationError("acquisition contig/REF compatibility did not pass")


def _annotation_provenance(row: Mapping[str, str]) -> Dict[str, Any]:
    return {
        "provider": "VectorBase/MalariaGEN",
        "release": row["native_annotation_version"],
        "assembly_accession": row["reference_accession"] + "/" + row["reference_version"],
        "fasta_assembly_accession": row["reference_accession"] + "/" + row["reference_version"],
        "status": "native",
        "genetic_code": 1,
        "invalid_transcript_policy": "exclude_with_audit",
        "excluded_transcripts": {
            "AGAP000192-RA": (
                "native AgamP4.12 CDS is empty or ambiguous; exclusion was declared in "
                "acquisition_manifest.tsv and the provider GFF remains byte-identical"
            )
        },
    }


def _samples(row: Mapping[str, str]) -> List[str]:
    return [
        line.strip()
        for line in Path(row["sample_list_path"]).read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]


def run_one(manifest: Path, index: int, output_dir: Path) -> Path:
    rows = _read_manifest(manifest)
    if index < 0 or index >= len(rows):
        raise Tier3ValidationError("tuple index is outside approved manifest rows")
    row = rows[index]
    _verify_inputs(row)
    result = compute_population_pi(
        dataset_id=row["tuple_id"],
        vcf_path=Path(row["callset_path"]),
        fasta_path=Path(row["reference_path"]),
        selected_samples=_samples(row),
        design=row["design"],
        denominator_kind="cohort_callable_mask",
        callable_bed_path=Path(row["callable_mask_path"]),
        gff_path=Path(row["native_annotation_path"]),
        annotation_metadata=_annotation_provenance(row),
        sampling_unit_jackknife=True,
    )
    result["tuple"] = {
        key: row[key]
        for key in (
            "tuple_id",
            "biological",
            "scientific_name",
            "taxon_id",
            "population_id",
            "population_release",
            "design",
            "region",
            "reference_accession",
            "reference_version",
            "native_annotation_accession",
            "native_annotation_version",
            "callset_modality",
        )
    }
    result["acquisition_manifest_sha256"] = sha256_file(manifest)
    output_dir.mkdir(parents=True, exist_ok=True)
    output = output_dir / (row["tuple_id"] + ".json")
    temporary = output.with_suffix(".json.partial-{}".format(os.getpid()))
    temporary.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(output)
    return output


def _uncertainty_fields(value: Mapping[str, Any]) -> Dict[str, Any]:
    bootstrap = value.get("bootstrap")
    if isinstance(bootstrap, Mapping) and bootstrap.get("interval") is not None:
        interval = bootstrap["interval"]
        if not isinstance(interval, list) or len(interval) != 2:
            raise Tier3ValidationError("genomic bootstrap lacks a two-sided interval")
        return {
            "uncertainty_method": "chromosome_stratified_block_bootstrap",
            "uncertainty_unit": "1-Mb genomic block resampled within chromosome",
            "uncertainty_replicates": bootstrap["replicates"],
            "uncertainty_standard_error": "",
            "interval_low": interval[0],
            "interval_high": interval[1],
            "interval_type": "percentile_95_percent",
            "frozen_genomic_bootstrap_status": "available:{}_eligible_1Mb_blocks".format(
                bootstrap["eligible_blocks"]
            ),
        }

    uncertainty = value.get("uncertainty")
    if not isinstance(uncertainty, Mapping):
        raise Tier3ValidationError("biological estimate lacks sampling-unit uncertainty")
    interval = uncertainty.get("interval")
    if not isinstance(interval, list) or len(interval) != 2:
        raise Tier3ValidationError("biological estimate lacks a two-sided uncertainty interval")
    if isinstance(bootstrap, Mapping):
        bootstrap_status = "unavailable:" + str(bootstrap.get("unavailable_reason"))
    else:
        bootstrap_status = "not_applicable:class_component_uses_sampling_unit_jackknife"
    return {
        "uncertainty_method": uncertainty["method"],
        "uncertainty_unit": uncertainty["unit"],
        "uncertainty_replicates": uncertainty["replicates"],
        "uncertainty_standard_error": uncertainty["standard_error"],
        "interval_low": interval[0],
        "interval_high": interval[1],
        "interval_type": uncertainty["interval_type"],
        "frozen_genomic_bootstrap_status": bootstrap_status,
    }


def _base_row(result: Mapping[str, Any], raw_path: Path) -> Dict[str, Any]:
    row = result["tuple"]
    return {
        "tuple_id": row["tuple_id"],
        "biological": row["biological"],
        "scientific_name": row["scientific_name"],
        "taxon_id": row["taxon_id"],
        "population_id": row["population_id"],
        "population_release": row["population_release"],
        "design": row["design"],
        "source_modality": row["callset_modality"],
        "region": row["region"],
        "reference_accession": row["reference_accession"] + "/" + row["reference_version"],
        "native_annotation_accession": row["native_annotation_accession"] + "/" + row["native_annotation_version"],
        "eligible_sample_size": result["sample_design"]["sampling_units"],
        "nominal_chromosomes": result["sample_design"]["nominal_chromosomes"],
        "exclusions": json.dumps(result["exclusion_counts"], sort_keys=True, separators=(",", ":")),
        "raw_result_path": str(raw_path),
        "component_s_numerator": "",
        "component_s_callable_sites": "",
        "component_w_numerator": "",
        "component_w_callable_sites": "",
    }


def _diversity_rows(result: Mapping[str, Any], raw_path: Path) -> List[Dict[str, Any]]:
    base = _base_row(result, raw_path)
    population = result["population_pi"]
    rows: List[Dict[str, Any]] = []
    whole = dict(base)
    whole.update(
        {
            "statistic": "population_pi",
            "annotation_category": "callable_nuclear_region",
            "variant_count": population["polymorphic_snv_sites"],
            "numerator": population["numerator"],
            "numerator_definition": "sum_unbiased_pairwise_differences",
            "callable_site_denominator": population["callable_sites"],
            "estimate": population["point_estimate"],
        }
    )
    whole.update(_uncertainty_fields(population))
    rows.append(whole)

    ratio = result.get("pi_S_over_pi_W")
    if not isinstance(ratio, Mapping):
        raise Tier3ValidationError("{} lacks 4D class results".format(result["dataset_id"]))
    for class_name in ("S", "W"):
        value = ratio[class_name]
        class_row = dict(base)
        class_row.update(
            {
                "statistic": "pi_" + class_name,
                "annotation_category": "native_4D_forward_reference_" + ("GC" if class_name == "S" else "AT"),
                "variant_count": value["polymorphic_snv_sites"],
                "numerator": value["numerator"],
                "numerator_definition": "sum_unbiased_pairwise_differences_in_4D_" + class_name,
                "callable_site_denominator": value["callable_sites"],
                "estimate": value["point_estimate"],
            }
        )
        class_row.update(_uncertainty_fields(value))
        rows.append(class_row)

    if ratio.get("point_estimate") is None:
        return rows

    ratio_row = dict(base)
    ratio_row.update(
        {
            "statistic": "pi_S_over_pi_W",
            "annotation_category": "native_4D_reference_conditioned_ratio",
            "variant_count": ratio["S"]["polymorphic_snv_sites"] + ratio["W"]["polymorphic_snv_sites"],
            "numerator": ratio["S"]["numerator"],
            "numerator_definition": "S_diversity_sum; estimate_is_(S_num/S_sites)/(W_num/W_sites)",
            "callable_site_denominator": ratio["S"]["callable_sites"],
            "estimate": ratio["point_estimate"],
            "component_s_numerator": ratio["S"]["numerator"],
            "component_s_callable_sites": ratio["S"]["callable_sites"],
            "component_w_numerator": ratio["W"]["numerator"],
            "component_w_callable_sites": ratio["W"]["callable_sites"],
        }
    )
    ratio_row.update(_uncertainty_fields(ratio))
    rows.append(ratio_row)
    return rows


def _write_tsv(path: Path, fieldnames: Sequence[str], rows: Iterable[Mapping[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def _read_telemetry(path: Path) -> Dict[str, Dict[str, str]]:
    if not path.is_file():
        raise Tier3ValidationError("missing scheduler telemetry: {}".format(path))
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    result: Dict[str, Dict[str, str]] = {}
    for row in rows:
        tuple_id = row.get("tuple_id", "")
        if tuple_id:
            result[tuple_id] = row
    return result


def scheduler_telemetry(
    manifest: Path, sacct_path: Path, array_job_id: str, output: Path
) -> None:
    approved = _read_manifest(manifest)
    with sacct_path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="|"))
    by_job_id = {row.get("JobID", ""): row for row in rows}
    telemetry: List[Dict[str, Any]] = []
    for array_index, manifest_row in enumerate(approved):
        job_id = "{}_{}".format(array_job_id, array_index)
        task = by_job_id.get(job_id)
        batch = by_job_id.get(job_id + ".batch", {})
        if task is None:
            raise Tier3ValidationError("sacct lacks array task {}".format(job_id))
        state = task["State"].split()[0]
        if state != "COMPLETED" or task["ExitCode"] != "0:0":
            raise Tier3ValidationError("sacct task {} did not complete successfully".format(job_id))
        telemetry.append(
            {
                "tuple_id": manifest_row["tuple_id"],
                "job_id": job_id,
                "array_task_id": array_index,
                "state": state,
                "elapsed": task["Elapsed"],
                "total_cpu": task["TotalCPU"],
                "max_rss": batch.get("MaxRSS") or task.get("MaxRSS") or "not_reported_by_sacct",
                "req_mem": task["ReqMem"],
                "alloc_cpus": task["AllocCPUS"],
                "exit_code": task["ExitCode"],
            }
        )
    _write_tsv(output, tuple(telemetry[0]), telemetry)


def collect(manifest: Path, raw_dir: Path, output_dir: Path, telemetry_path: Path, environment: Path) -> None:
    approved = _read_manifest(manifest)
    environment_record = json.loads(environment.read_text(encoding="utf-8"))
    telemetry = _read_telemetry(telemetry_path)
    diversity: List[Dict[str, Any]] = []
    samples: List[Dict[str, Any]] = []
    annotation_exclusions: List[Dict[str, Any]] = []
    run_rows: List[Dict[str, Any]] = []
    raw_results: Dict[str, Dict[str, Any]] = {}
    for manifest_row in approved:
        tuple_id = manifest_row["tuple_id"]
        raw_path = raw_dir / (tuple_id + ".json")
        if not raw_path.is_file():
            raise Tier3ValidationError("missing raw biological result: {}".format(raw_path))
        result = json.loads(raw_path.read_text(encoding="utf-8"))
        if result.get("dataset_id") != tuple_id:
            raise Tier3ValidationError("raw result tuple identity mismatch")
        raw_results[tuple_id] = result
        diversity.extend(_diversity_rows(result, raw_path))
        declared = result["annotation"].get("declared_excluded_transcripts", {})
        for transcript_id, reason in sorted(result["annotation"]["excluded_transcripts"].items()):
            annotation_exclusions.append(
                {
                    "tuple_id": tuple_id,
                    "transcript_id": transcript_id,
                    "exclusion_source": "acquisition_declared" if transcript_id in declared else "run_validated",
                    "reason": reason,
                    "native_gff_sha256": result["annotation"]["gff_sha256"],
                }
            )
        for sample_id, summary in sorted(result["sample_summaries"].items()):
            samples.append(
                {
                    "tuple_id": tuple_id,
                    "population_id": manifest_row["population_id"],
                    "sample_id": sample_id,
                    **summary,
                }
            )
        scheduler = telemetry.get(tuple_id)
        if scheduler is None or scheduler.get("state") != "COMPLETED" or scheduler.get("exit_code") != "0:0":
            raise Tier3ValidationError("missing successful sacct row for {}".format(tuple_id))
        run_rows.append(
            {
                "tuple_id": tuple_id,
                "status": "accepted",
                "manifest_sha256": sha256_file(manifest),
                "reference_path": manifest_row["reference_path"],
                "reference_sha256": manifest_row["reference_sha256"],
                "annotation_path": manifest_row["native_annotation_path"],
                "annotation_sha256": manifest_row["native_annotation_sha256"],
                "callset_path": manifest_row["callset_path"],
                "callset_sha256": manifest_row["callset_sha256"],
                "callable_mask_path": manifest_row["callable_mask_path"],
                "callable_mask_sha256": manifest_row["callable_mask_sha256"],
                "sample_list_path": manifest_row["sample_list_path"],
                "sample_list_sha256": manifest_row["sample_list_sha256"],
                "raw_result_path": str(raw_path),
                "raw_result_sha256": sha256_file(raw_path),
                "guix_manifest": manifest_row["guix_manifest"],
                "guix_manifest_sha256": environment_record["manifest_sha256"],
                "guix_channels": manifest_row["guix_channels"],
                "guix_channels_sha256": environment_record["channels_sha256"],
                "guix_channel_commit": environment_record["channel_commit"],
                "guix_profile": environment_record["profile_store_path"],
                "tool_versions": json.dumps(environment_record["tool_versions"], sort_keys=True, separators=(",", ":")),
                "slurm_job_id": scheduler["job_id"],
                "slurm_array_task_id": scheduler["array_task_id"],
                "slurm_state": scheduler["state"],
                "slurm_elapsed": scheduler["elapsed"],
                "slurm_total_cpu": scheduler["total_cpu"],
                "slurm_max_rss": scheduler["max_rss"],
                "slurm_req_mem": scheduler["req_mem"],
                "slurm_alloc_cpus": scheduler["alloc_cpus"],
                "slurm_exit_code": scheduler["exit_code"],
            }
        )

    _write_tsv(output_dir / "population_diversity.tsv", DIVERSITY_FIELDS, diversity)
    _write_tsv(
        output_dir / "population_sample_summary.tsv",
        (
            "tuple_id",
            "population_id",
            "sample_id",
            "records_in_callable_mask",
            "called_chromosome_observations",
            "missing_chromosome_observations",
            "heterozygous_snv_sites",
            "nonreference_genotype_sites",
        ),
        samples,
    )
    _write_tsv(
        output_dir / "population_annotation_exclusions.tsv",
        ("tuple_id", "transcript_id", "exclusion_source", "reason", "native_gff_sha256"),
        annotation_exclusions,
    )
    _write_tsv(
        output_dir / "population_run_manifest.tsv",
        tuple(run_rows[0]),
        run_rows,
    )
    failure_fields = ("tuple_id", "population_id", "stage", "status", "job_id", "reason")
    failed_attempts_path = output_dir / "run_logs" / "population_failed_attempts.tsv"
    failed_attempts: List[Dict[str, str]] = []
    if failed_attempts_path.is_file():
        with failed_attempts_path.open("r", encoding="utf-8", newline="") as handle:
            failed_attempts = list(csv.DictReader(handle, delimiter="\t"))
    for tuple_id, result in sorted(raw_results.items()):
        ratio = result.get("pi_S_over_pi_W") or {}
        if ratio.get("point_estimate") is None:
            failed_attempts.append(
                {
                    "tuple_id": tuple_id,
                    "population_id": result["tuple"]["population_id"],
                    "stage": "frozen_pi_S_over_pi_W_power_gate",
                    "status": "blocked_on_repair-tier-3b",
                    "job_id": next(
                        row["slurm_job_id"] for row in run_rows if row["tuple_id"] == tuple_id
                    ),
                    "reason": str(ratio.get("unavailable_reason")),
                }
            )
    if all((result.get("pi_S_over_pi_W") or {}).get("point_estimate") is not None for result in raw_results.values()):
        failed_attempts.append(
            {
                "tuple_id": "none",
                "population_id": "none",
                "stage": "final_full_run",
                "status": "no_unresolved_failures",
                "job_id": "none",
                "reason": "all approved biological acquisition tuples completed and were retained",
            }
        )
    _write_tsv(output_dir / "population_failures.tsv", failure_fields, failed_attempts)
    _write_qc(output_dir / "population_qc.md", approved, raw_results, environment_record)


def _write_qc(
    path: Path,
    approved: Sequence[Mapping[str, str]],
    results: Mapping[str, Mapping[str, Any]],
    environment: Mapping[str, Any],
) -> None:
    ratios_complete = all(
        (result.get("pi_S_over_pi_W") or {}).get("point_estimate") is not None
        for result in results.values()
    )
    status = (
        "**PASS — every approved non-synthetic biological tuple and frozen ratio completed.**"
        if ratios_complete
        else "**BLOCKED — biological pi/pi_S/pi_W completed, but frozen pi_S/pi_W power requires `repair-tier-3b`.**"
    )
    lines = [
        "# Tier 3B recovered population analysis QC",
        "",
        "Status: " + status,
        "",
        "## Biological estimates",
        "",
        "| tuple | population | n | polymorphic SNVs | callable sites | population pi | pi_S/pi_W | secondary sample-jackknife SE (pi) |",
        "|---|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in approved:
        result = results[row["tuple_id"]]
        population = result["population_pi"]
        ratio = result["pi_S_over_pi_W"]
        ratio_text = (
            "{:.12g}".format(ratio["point_estimate"])
            if ratio["point_estimate"] is not None
            else "unavailable:" + str(ratio["unavailable_reason"])
        )
        lines.append(
            "| `{}` | `{}` | {} | {} | {} | {:.12g} | {} | {:.12g} |".format(
                row["tuple_id"],
                row["population_id"],
                result["sample_design"]["sampling_units"],
                population["polymorphic_snv_sites"],
                population["callable_sites"],
                population["point_estimate"],
                ratio_text,
                population["uncertainty"]["standard_error"],
            )
        )
    lines.extend(
        [
            "",
            "## Coordinate, reference, genotype, and annotation gates",
            "",
            "- VCF positions were converted from 1-based POS to zero-based internal coordinates; callable BED remained 0-based half-open; GFF3 was parsed as 1-based closed and converted once.",
            "- The complete eight-contig VCF and AgamP4 FASTA dictionaries matched. Every streamed REF allele was checked against the exact FASTA; duplicates and coordinate disorder were fatal.",
            "- Multiallelic A/C/G/T SNVs were counted once with all allele counts. Indels, symbolic/non-SNV records, filters, ambiguous reference bases, and sites below 36 called chromosomes were removed from numerator and denominator.",
            "- Callable denominators came only from each exact 20-sample cohort BED. No absent sparse-VCF row was treated as callable reference.",
            "- Native VectorBase AgamP4.12 annotation was consumed byte-for-byte. `AGAP000192-RA` was excluded exactly as declared upstream; all additional canonical CDS failures were enumerated with exact reasons in `population_annotation_exclusions.tsv` and the raw audits rather than editing the GFF.",
            "- 4D sites use nuclear code 1 and forward-reference S=G/C, W=A/T. The frozen 10,000-callable-site gate is enforced without relaxation; any underpowered ratio is explicit in the table and failure ledger.",
            "",
            "## Uncertainty",
            "",
            "The powered 21-Mb tuples meet the frozen minimum of 20 eligible genomic blocks. Population pi and pi_S/pi_W use the deterministic 10,000-replicate, chromosome-stratified 1-Mb block bootstrap with a percentile 95% interval as their reported uncertainty. The pi_S and pi_W component rows use a 20-replicate delete-one-selected-individual jackknife with a normal-approximation 95% interval and SE, conditional on the exact cohort callable mask. Each TSV row labels its method and resampling unit explicitly; raw JSON retains both uncertainty calculations where both apply.",
            "",
            "## Reproducible environment and scheduler",
            "",
            "- Guix channel commit: `{}`".format(environment["channel_commit"]),
            "- Guix profile: `{}`".format(environment["profile_store_path"]),
            "- Manifest SHA-256: `{}`".format(environment["manifest_sha256"]),
            "- Channels SHA-256: `{}`".format(environment["channels_sha256"]),
            "- Tool versions: `{}`".format(json.dumps(environment["tool_versions"], sort_keys=True)),
            "- Heavy tuple analyses and the independent subset reconciliation ran through Slurm. Exact job IDs and `sacct` resource records are in `population_run_manifest.tsv` and `run_logs/`.",
            "- Array concurrency was throttled to `%1` to avoid two simultaneous full scans of shared MooseFS inputs.",
            "",
            "## Independent reconciliation",
            "",
            "`population_independent_check.tsv` is produced by a separate standard-library VCF/GT parser on the first 100 kb. The production side consumes an exact tabix-indexed `bcftools view --regions` slice; the manual side reads the original BGZF stream and stops at the same coordinate boundary. It independently clips the BED denominator, checks REF, missingness, filters, multiallelic SNVs, and recomputes pairwise diversity, then requires exact denominator/variant-count agreement and numerical agreement within 1e-12 with the production subset calculation.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _parse_region(region: str, subset_bp: int = 100_000) -> Tuple[str, int, int]:
    match = re.fullmatch(r"([^:]+):(\d+)-(\d+)", region)
    if not match:
        raise Tier3ValidationError("malformed acquisition region")
    contig, start_text, end_text = match.groups()
    start = int(start_text) - 1
    end = min(int(end_text), start + subset_bp)
    return contig, start, end


def _clip_bed(path: Path, contig: str, start: int, end: int) -> List[Tuple[int, int]]:
    intervals: List[Tuple[int, int]] = []
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            fields = raw.rstrip("\n").split("\t")
            if fields[0] != contig:
                continue
            clipped_start, clipped_end = max(start, int(fields[1])), min(end, int(fields[2]))
            if clipped_start < clipped_end:
                intervals.append((clipped_start, clipped_end))
    return intervals


def _inside(position: int, intervals: Sequence[Tuple[int, int]]) -> bool:
    return any(start <= position < end for start, end in intervals)


def _manual_subset(row: Mapping[str, str], intervals: Sequence[Tuple[int, int]], contig: str, start: int, end: int) -> Dict[str, Any]:
    reference = subprocess.run(
        ["samtools", "faidx", row["reference_path"], "{}:{}-{}".format(contig, start + 1, end)],
        check=True,
        text=True,
        stdout=subprocess.PIPE,
    ).stdout
    sequence = "".join(line.strip() for line in reference.splitlines() if not line.startswith(">"))
    callable_sites = sum(
        sum(base in "ACGT" for base in sequence[left - start : right - start].upper())
        for left, right in intervals
    )
    samples = _samples(row)
    diversity_sum = 0.0
    polymorphic = 0
    exclusions: Dict[str, int] = {}
    with gzip.open(row["callset_path"], "rt", encoding="utf-8") as handle:
        sample_columns: List[int] = []
        for raw in handle:
            if raw.startswith("##"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if raw.startswith("#CHROM"):
                header_samples = fields[9:]
                sample_columns = [header_samples.index(sample) + 9 for sample in samples]
                continue
            if raw.startswith("#"):
                continue
            if fields[0] != contig:
                continue
            position = int(fields[1]) - 1
            if position >= end:
                break
            if position < start or not _inside(position, intervals):
                continue
            ref, alt, filter_text = fields[3].upper(), fields[4].upper(), fields[6]
            format_fields = fields[8].split(":")
            gt_index = format_fields.index("GT")
            genotypes: List[List[int | None]] = []
            for column in sample_columns:
                gt_text = fields[column].split(":")[gt_index]
                genotypes.append(
                    [None if allele == "." else int(allele) for allele in re.split(r"[/|]", gt_text)]
                )
            called = [allele for genotype in genotypes for allele in genotype if allele is not None]
            alts = [] if alt in ("", ".") else alt.split(",")
            snv = len(ref) == 1 and ref in "ACGT" and bool(alts) and all(len(value) == 1 and value in "ACGT" for value in alts)
            reason = None
            expected_ref = sequence[position - start : position - start + len(ref)].upper()
            if expected_ref != ref:
                raise Tier3ValidationError("independent subset REF mismatch")
            if filter_text not in ("PASS", ".", ""):
                reason = "filtered"
            elif ref not in "ACGT":
                reason = "ambiguous_reference"
            elif len(called) < 36:
                reason = "insufficient_called_chromosomes"
            elif alts and not snv:
                reason = "non_snv"
            if reason:
                callable_sites -= 1
                exclusions[reason] = exclusions.get(reason, 0) + 1
                continue
            if not snv:
                continue
            counts: Dict[int, int] = {}
            for allele in called:
                counts[allele] = counts.get(allele, 0) + 1
            n = len(called)
            contribution = 1.0 - sum(count * (count - 1) for count in counts.values()) / (n * (n - 1))
            diversity_sum += contribution
            if contribution > 0:
                polymorphic += 1
    if callable_sites <= 0 or polymorphic <= 0:
        raise Tier3ValidationError("independent subset has no usable biological signal")
    return {
        "callable_sites": callable_sites,
        "diversity_sum": diversity_sum,
        "polymorphic_snv_sites": polymorphic,
        "estimate": diversity_sum / callable_sites,
        "exclusions": exclusions,
    }


def _extract_subset_vcf(
    source: Path, destination: Path, contig: str, start: int, end: int
) -> None:
    region = "{}:{}-{}".format(contig, start + 1, end)
    command = [
        "bcftools",
        "view",
        "--no-version",
        "--regions",
        region,
        "--output-type",
        "z",
        "--output",
        str(destination),
        str(source),
    ]
    print("independent production subset command: " + " ".join(command), flush=True)
    subprocess.run(command, check=True)
    if not destination.is_file() or destination.stat().st_size == 0:
        raise Tier3ValidationError("bcftools produced no independent subset VCF")


def independent_check(manifest: Path, output_dir: Path) -> None:
    rows: List[Dict[str, Any]] = []
    output_dir.mkdir(parents=True, exist_ok=True)
    scratch_root = Path(os.environ.get("TIER3_SCRATCH_ROOT", str(output_dir)))
    scratch_root.mkdir(parents=True, exist_ok=True)
    for row in _read_manifest(manifest):
        _verify_inputs(row)
        contig, start, end = _parse_region(row["region"])
        intervals = _clip_bed(Path(row["callable_mask_path"]), contig, start, end)
        clipped_bed = output_dir / ("population_subset_" + row["tuple_id"] + ".bed")
        clipped_bed.write_text(
            "".join("{}\t{}\t{}\n".format(contig, left, right) for left, right in intervals),
            encoding="utf-8",
        )
        with tempfile.TemporaryDirectory(
            prefix="{}-subset-".format(row["tuple_id"]), dir=str(scratch_root)
        ) as temporary_directory:
            subset_vcf = Path(temporary_directory) / "subset.vcf.gz"
            _extract_subset_vcf(
                Path(row["callset_path"]), subset_vcf, contig, start, end
            )
            production = compute_population_pi(
                dataset_id=row["tuple_id"] + "_independent_subset",
                vcf_path=subset_vcf,
                fasta_path=Path(row["reference_path"]),
                selected_samples=_samples(row),
                design=row["design"],
                denominator_kind="cohort_callable_mask",
                callable_bed_path=clipped_bed,
                bootstrap_replicates=1,
            )["population_pi"]
        manual = _manual_subset(row, intervals, contig, start, end)
        checks = {
            "callable_sites": manual["callable_sites"] == production["callable_sites"],
            "polymorphic_snv_sites": manual["polymorphic_snv_sites"] == production["polymorphic_snv_sites"],
            "diversity_sum": math.isclose(manual["diversity_sum"], production["diversity_sum"], rel_tol=0, abs_tol=1e-12),
            "estimate": math.isclose(manual["estimate"], production["point_estimate"], rel_tol=0, abs_tol=1e-12),
        }
        if not all(checks.values()):
            raise Tier3ValidationError("independent subset reconciliation failed: {}".format(checks))
        rows.append(
            {
                "tuple_id": row["tuple_id"],
                "population_id": row["population_id"],
                "subset": "{}:{}-{}".format(contig, start + 1, end),
                "eligible_sample_size": 20,
                "manual_variant_count": manual["polymorphic_snv_sites"],
                "production_variant_count": production["polymorphic_snv_sites"],
                "manual_numerator": manual["diversity_sum"],
                "production_numerator": production["diversity_sum"],
                "manual_callable_sites": manual["callable_sites"],
                "production_callable_sites": production["callable_sites"],
                "manual_estimate": manual["estimate"],
                "production_estimate": production["point_estimate"],
                "absolute_difference": abs(manual["estimate"] - production["point_estimate"]),
                "status": "PASS",
            }
        )
    _write_tsv(output_dir.parent / "population_independent_check.tsv", tuple(rows[0]), rows)


def validate(diversity_path: Path, manifest: Path, independent_path: Path) -> None:
    approved = _read_manifest(manifest)
    with diversity_path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise Tier3ValidationError("population diversity table has no data rows")
    biological = {row["tuple_id"] for row in rows if row["biological"].lower() == "true"}
    approved_ids = {row["tuple_id"] for row in approved}
    if biological != approved_ids or len(biological) < 2:
        raise Tier3ValidationError("diversity table does not contain every approved biological tuple")
    for tuple_id in approved_ids:
        tuple_rows = [row for row in rows if row["tuple_id"] == tuple_id]
        if {row["statistic"] for row in tuple_rows} != set(STATISTICS):
            raise Tier3ValidationError("{} lacks required whole/annotation statistics".format(tuple_id))
    required_text = (
        "uncertainty_method",
        "uncertainty_unit",
        "uncertainty_replicates",
        "interval_low",
        "interval_high",
        "interval_type",
    )
    for row in rows:
        for field in ("variant_count", "numerator", "callable_site_denominator"):
            if not math.isfinite(float(row[field])) or float(row[field]) <= 0:
                raise Tier3ValidationError("non-positive {} in accepted row".format(field))
        if not math.isfinite(float(row["estimate"])):
            raise Tier3ValidationError("non-finite accepted estimate")
        if any(not row[field] for field in required_text):
            raise Tier3ValidationError("accepted estimate has incomplete uncertainty")
        method = row["uncertainty_method"]
        replicates = int(row["uncertainty_replicates"])
        if method == "chromosome_stratified_block_bootstrap":
            if replicates != 10000 or not row["frozen_genomic_bootstrap_status"].startswith("available:"):
                raise Tier3ValidationError("incomplete frozen genomic bootstrap evidence")
        elif method == "delete_one_sampling_unit_jackknife":
            if replicates != 20 or not row["uncertainty_standard_error"]:
                raise Tier3ValidationError("incomplete sampling-unit jackknife evidence")
        else:
            raise Tier3ValidationError("unexpected uncertainty method: {}".format(method))
    for tuple_id in approved_ids:
        by_statistic = {row["statistic"]: row for row in rows if row["tuple_id"] == tuple_id}
        for statistic in ("population_pi", "pi_S_over_pi_W"):
            if by_statistic[statistic]["uncertainty_method"] != "chromosome_stratified_block_bootstrap":
                raise Tier3ValidationError("{} {} lacks powered genomic bootstrap".format(tuple_id, statistic))
    with independent_path.open("r", encoding="utf-8", newline="") as handle:
        independent = list(csv.DictReader(handle, delimiter="\t"))
    if {row["tuple_id"] for row in independent if row["status"] == "PASS"} != approved_ids:
        raise Tier3ValidationError("independent checks do not pass for every approved tuple")
    print(
        "PASS: {} biological tuples, {} estimates, positive counts/denominators, finite estimates, complete uncertainty, independent reconciliation".format(
            len(approved_ids), len(rows)
        )
    )


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    run = subparsers.add_parser("run-one")
    run.add_argument("--manifest", type=Path, required=True)
    run.add_argument("--index", type=int, required=True)
    run.add_argument("--output-dir", type=Path, required=True)
    gather = subparsers.add_parser("collect")
    gather.add_argument("--manifest", type=Path, required=True)
    gather.add_argument("--raw-dir", type=Path, required=True)
    gather.add_argument("--output-dir", type=Path, required=True)
    gather.add_argument("--telemetry", type=Path, required=True)
    gather.add_argument("--environment", type=Path, required=True)
    check = subparsers.add_parser("independent-check")
    check.add_argument("--manifest", type=Path, required=True)
    check.add_argument("--output-dir", type=Path, required=True)
    assertion = subparsers.add_parser("validate")
    assertion.add_argument("--manifest", type=Path, required=True)
    assertion.add_argument("--diversity", type=Path, required=True)
    assertion.add_argument("--independent", type=Path, required=True)
    telemetry = subparsers.add_parser("telemetry")
    telemetry.add_argument("--manifest", type=Path, required=True)
    telemetry.add_argument("--sacct", type=Path, required=True)
    telemetry.add_argument("--array-job", required=True)
    telemetry.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    if args.command == "run-one":
        print(run_one(args.manifest, args.index, args.output_dir))
    elif args.command == "collect":
        collect(args.manifest, args.raw_dir, args.output_dir, args.telemetry, args.environment)
    elif args.command == "independent-check":
        independent_check(args.manifest, args.output_dir)
    elif args.command == "validate":
        validate(args.diversity, args.manifest, args.independent)
    else:
        scheduler_telemetry(args.manifest, args.sacct, args.array_job, args.output)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Tier3ValidationError as error:
        raise SystemExit("tier3b recovery rejected: {}".format(error))
