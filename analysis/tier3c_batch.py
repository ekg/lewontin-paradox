#!/usr/bin/env python3
"""Freeze, stage, run, and collect the production Tier 3c composition batch.

Discovery and analysis are deliberately separate.  ``discover`` freezes one
exact, same-species NCBI assembly accession for every Buffalo core row.
``stage-one`` then acquires the exact NCBI payloads, verifies the provider MD5,
removes assembly-report-labelled non-nuclear sequences, and records SHA-256
checksums.  Only ``freeze`` output is accepted by ``run-one``.

Raw FASTA/GFF files and scheduler telemetry live outside git.  ``collect``
emits the small result TSV, per-species provenance/QC JSON, and a complete
failure ledger intended for version control.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import resource
import re
import shlex
import shutil
import socket
import statistics
import subprocess
import sys
import tempfile
import time
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

from analysis.tier3_common import Tier3ValidationError, sha256_file, verify_file
from analysis.tier3_manifest import buffalo_core_rows
from analysis.tier3c_ncbi_gc import (
    USER_AGENT,
    analyze_dataset,
    atomic_write_json,
    discover_assemblies,
    ncbi_artifact_urls,
    select_assembly,
    validate_pilot,
)


BATCH_SCHEMA = "tier3c-batch-v1"
DISCOVERY_SCHEMA = "tier3c-discovery-v1"
STAGE_SCHEMA = "tier3c-stage-v1"
MINIMUM_NUCLEAR_REFERENCE_BASES = 10_000_000
RESULT_COLUMNS = (
    "dataset_id",
    "scientific_name",
    "taxon_id",
    "buffalo_diversity",
    "buffalo_pred_log10_N",
    "assembly_accession_version",
    "whole_genome_gc",
    "whole_genome_callable_bases",
    "gc3",
    "gc3_callable_third_positions",
    "callable_genes",
    "annotation_status",
)
FAILURE_COLUMNS = (
    "scientific_name",
    "phase",
    "reason_code",
    "detail",
    "reproducible_command",
)
MANIFEST_COLUMNS = (
    "dataset_id",
    "scientific_name",
    "taxon_id",
    "assembly_accession_version",
    "provider",
    "release",
    "fasta_source_uri",
    "fasta_upstream_md5",
    "fasta_sha256",
    "fasta_size_bytes",
    "gff_source_uri",
    "gff_upstream_md5",
    "gff_sha256",
    "gff_size_bytes",
    "annotation_native_exact_reference",
    "genetic_code",
)


def _slug(name: str) -> str:
    return ".".join(part.lower().replace("-", "_") for part in name.split())


def _read_json(path: Path) -> Mapping[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise Tier3ValidationError(f"cannot read {path}: {error}") from error
    if not isinstance(value, Mapping):
        raise Tier3ValidationError(f"{path} is not a JSON object")
    return value


def _atomic_tsv(path: Path, columns: Sequence[str], rows: Iterable[Mapping[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=str(path.parent), text=True)
    temporary = Path(temporary_name)
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for row in rows:
                writer.writerow({column: row.get(column, "") for column in columns})
            handle.flush()
            os.fsync(handle.fileno())
        if path.is_file() and path.read_bytes() == temporary.read_bytes():
            temporary.unlink()
            return
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def _read_tsv_exact(path: Path, columns: Sequence[str]) -> List[Dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != tuple(columns):
            raise Tier3ValidationError(
                f"{path} schema mismatch: expected {tuple(columns)!r}, observed {tuple(reader.fieldnames or ())!r}"
            )
        return [dict(row) for row in reader]


def validate_collected(output_dir: Path) -> Mapping[str, Any]:
    manifest = _read_tsv_exact(output_dir / "tier3c_manifest.tsv", MANIFEST_COLUMNS)
    data = _read_tsv_exact(output_dir / "tier3c_data.tsv", RESULT_COLUMNS)
    failures = _read_tsv_exact(output_dir / "tier3c_failure_ledger.tsv", FAILURE_COLUMNS)
    manifest_ids = [row["dataset_id"] for row in manifest]
    data_ids = [row["dataset_id"] for row in data]
    if len(manifest_ids) != len(set(manifest_ids)) or len(data_ids) != len(set(data_ids)):
        raise Tier3ValidationError("collected Tier 3c tables contain duplicate dataset IDs")
    if not set(data_ids) <= set(manifest_ids):
        raise Tier3ValidationError("result table contains a dataset absent from the frozen manifest TSV")
    for row in manifest:
        if not re.fullmatch(r"[0-9a-f]{64}", row["fasta_sha256"]):
            raise Tier3ValidationError(f"invalid FASTA SHA-256 for {row['dataset_id']}")
        if row["gff_sha256"] and not re.fullmatch(r"[0-9a-f]{64}", row["gff_sha256"]):
            raise Tier3ValidationError(f"invalid GFF SHA-256 for {row['dataset_id']}")
        for field in ("fasta_upstream_md5", "gff_upstream_md5"):
            if row[field] and not re.fullmatch(r"[0-9a-f]{32}", row[field]):
                raise Tier3ValidationError(f"invalid provider MD5 for {row['dataset_id']}:{field}")
    for row in data:
        for field in ("whole_genome_gc", "gc3"):
            if row[field] and not 0 <= float(row[field]) <= 1:
                raise Tier3ValidationError(f"{field} outside [0,1] for {row['dataset_id']}")
        for field in ("whole_genome_callable_bases", "gc3_callable_third_positions", "callable_genes"):
            if row[field] and int(row[field]) <= 0:
                raise Tier3ValidationError(f"non-positive {field} for {row['dataset_id']}")
        if not (output_dir / "tier3c_qc" / f"{row['dataset_id']}.json").is_file():
            raise Tier3ValidationError(f"missing per-species QC for {row['dataset_id']}")
    summary = _read_json(output_dir / "tier3c_qc_summary.json")
    if summary.get("completed") != len(data) or summary.get("failures") != len(failures):
        raise Tier3ValidationError("QC summary counts disagree with collected TSV schemas")
    return {
        "manifest_rows": len(manifest),
        "result_rows": len(data),
        "failure_rows": len(failures),
        "schema_and_checksum_fields_valid": True,
    }


def _retry_discovery(name: str, attempts: int = 4) -> List[Any]:
    for attempt in range(attempts):
        try:
            return discover_assemblies(name)
        except Tier3ValidationError:
            if attempt + 1 == attempts:
                raise
            time.sleep(2**attempt)
    raise AssertionError("unreachable")


def discover(buffalo: Path, output: Path, *, delay_seconds: float = 0.36) -> Mapping[str, Any]:
    rows = buffalo_core_rows(buffalo)
    datasets: List[Dict[str, Any]] = []
    for index, row in enumerate(rows):
        name = row["species"].strip()
        record: Dict[str, Any] = {
            "dataset_id": f"{_slug(name)}.tier3c",
            "scientific_name": name,
            "buffalo_diversity": row["diversity"],
            "buffalo_pred_log10_N": row["pred_log10_N"],
        }
        try:
            candidates = _retry_discovery(name)
            selected = select_assembly(candidates) if candidates else None
            if selected is None:
                record.update(status="ineligible", reason_code="missing_exact_same_species_ncbi_assembly")
            else:
                record.update(
                    status="selected",
                    scientific_name=name,
                    taxon_id=selected.species_taxid,
                    assembly_accession_version=selected.accession,
                    assembly_name=selected.assembly_name,
                    assembly_level=selected.assembly_level,
                    refseq_category=selected.refseq_category,
                    release=selected.release_date,
                    provider="NCBI RefSeq" if selected.accession.startswith("GCF_") else "NCBI GenBank",
                    artifact_urls=dict(ncbi_artifact_urls(selected)),
                    selection_policy="same_species_best_refseq_category_then_gcf_then_level_then_release",
                )
        except Tier3ValidationError as error:
            record.update(status="discovery_failed", reason_code="ncbi_discovery_error", detail=str(error))
        datasets.append(record)
        if index + 1 != len(rows):
            time.sleep(delay_seconds)
    value = {
        "schema_version": DISCOVERY_SCHEMA,
        "decision_version": "tier3-decisions-v1",
        "buffalo_sha256": sha256_file(buffalo),
        "congener_substitution_allowed": False,
        "datasets": datasets,
    }
    atomic_write_json(value, output)
    return value


def _download_unfrozen(url: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary_name = tempfile.mkstemp(prefix=f".{destination.name}.", dir=str(destination.parent))
    temporary = Path(temporary_name)
    try:
        request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
        with urllib.request.urlopen(request, timeout=600) as response, os.fdopen(fd, "wb") as handle:
            fd = -1
            shutil.copyfileobj(response, handle, length=1024 * 1024)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, destination)
    finally:
        if fd >= 0:
            os.close(fd)
        if temporary.exists():
            temporary.unlink()


def _provider_md5s(path: Path) -> Dict[str, str]:
    result: Dict[str, str] = {}
    for raw in path.read_text(encoding="utf-8").splitlines():
        fields = raw.split(None, 1)
        if len(fields) == 2:
            result[Path(fields[1].removeprefix("./")).name] = fields[0]
    return result


def _md5(path: Path) -> str:
    digest = hashlib.md5(usedforsecurity=False)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _excluded_non_nuclear(report: Path) -> Tuple[set[str], Dict[str, int]]:
    excluded: set[str] = set()
    counts: Dict[str, int] = {}
    for raw in report.read_text(encoding="utf-8").splitlines():
        if not raw or raw.startswith("#"):
            continue
        fields = raw.split("\t")
        if len(fields) < 9:
            raise Tier3ValidationError("NCBI assembly report has fewer than nine columns")
        assigned_molecule, location, unit = fields[2].lower(), fields[3].lower(), fields[7].lower()
        reason = None
        if unit == "non-nuclear":
            reason = "assembly_unit_non_nuclear"
        elif location in {"mitochondrion", "chloroplast", "plastid", "apicoplast"}:
            reason = f"location_{location}"
        elif any(token in assigned_molecule for token in ("mitochond", "chloroplast", "plastid", "apicoplast")):
            reason = "assigned_molecule_non_nuclear"
        if reason:
            counts[reason] = counts.get(reason, 0) + 1
            excluded.update(value for value in (fields[0], fields[2], fields[4], fields[6]) if value not in {"na", ""})
            if len(fields) > 9 and fields[9] not in {"na", ""}:
                excluded.add(fields[9])
    return excluded, counts


class _StagingMissingness(Tier3ValidationError):
    def __init__(self, reason_code: str, detail: str):
        super().__init__(detail)
        self.reason_code = reason_code


def _validate_nuclear_reference_size(fai: Path) -> int:
    total = sum(int(line.split("\t")[1]) for line in fai.read_text(encoding="utf-8").splitlines())
    if total < MINIMUM_NUCLEAR_REFERENCE_BASES:
        raise _StagingMissingness(
            "no_valid_nuclear_genome_assembly",
            f"retained nuclear reference has {total} bases; minimum is {MINIMUM_NUCLEAR_REFERENCE_BASES}",
        )
    return total


def _filter_fasta(source: Path, destination: Path, excluded: set[str]) -> int:
    kept = 0
    keep = False
    with gzip.open(source, "rt", encoding="utf-8") as input_handle, destination.open("w", encoding="utf-8") as output:
        for raw in input_handle:
            if raw.startswith(">"):
                name = raw[1:].split(None, 1)[0]
                keep = name not in excluded
                kept += int(keep)
            if keep:
                output.write(raw)
    if kept == 0:
        raise _StagingMissingness(
            "no_valid_nuclear_genome_assembly",
            "nuclear FASTA filtering retained no sequences",
        )
    return kept


def _filter_gff(source: Path, destination: Path, excluded: set[str]) -> Tuple[List[Tuple[str, int]], int]:
    regions: List[Tuple[str, int]] = []
    features = 0
    with gzip.open(source, "rt", encoding="utf-8") as input_handle, destination.open("w", encoding="utf-8") as output:
        for raw in input_handle:
            if raw.startswith("##FASTA"):
                break
            if raw.startswith("##sequence-region"):
                fields = raw.split()
                if len(fields) != 4:
                    raise Tier3ValidationError("malformed NCBI GFF sequence-region directive")
                if fields[1] in excluded:
                    continue
                regions.append((fields[1], int(fields[3])))
                output.write(raw)
            elif raw.startswith("#") or not raw.strip():
                output.write(raw)
            else:
                contig = raw.split("\t", 1)[0]
                if contig not in excluded:
                    features += 1
                    output.write(raw)
    return regions, features


def _bgzip_replace(path: Path) -> Path:
    """Replace a text payload with deterministic BGZF for bounded scratch use."""

    destination = Path(str(path) + ".gz")
    temporary = destination.with_name("." + destination.name + ".partial")
    try:
        with temporary.open("wb") as output:
            subprocess.run(["bgzip", "-@", "2", "-c", str(path)], stdout=output, check=True)
            output.flush()
            os.fsync(output.fileno())
        os.replace(temporary, destination)
        path.unlink()
        return destination
    finally:
        if temporary.exists():
            temporary.unlink()


def _artifact(path: Path, *, source_uri: Optional[str] = None, upstream_md5: Optional[str] = None) -> Dict[str, Any]:
    result: Dict[str, Any] = {
        "logical_name": path.name,
        "uri": path.resolve().as_uri(),
        "sha256": sha256_file(path),
        "size_bytes": path.stat().st_size,
    }
    if source_uri:
        result["source_uri"] = source_uri
    if upstream_md5:
        result["upstream_md5"] = upstream_md5
    return result


def stage_one(
    discovery_path: Path, index: int, stage_root: Path, environment_record: Path
) -> Mapping[str, Any]:
    discovery_value = _read_json(discovery_path)
    if discovery_value.get("schema_version") != DISCOVERY_SCHEMA:
        raise Tier3ValidationError("stage-one requires a Tier 3c discovery manifest")
    selected = sorted(
        (row for row in discovery_value["datasets"] if row["status"] == "selected"),
        key=lambda row: row["dataset_id"].encode("utf-8"),
    )
    if not 0 <= index < len(selected):
        raise Tier3ValidationError(f"array index {index} outside selected dataset count {len(selected)}")
    row = selected[index]
    dataset_dir = stage_root / row["dataset_id"]
    dataset_dir.mkdir(parents=True, exist_ok=True)
    sidecar = dataset_dir / "stage.json"
    if sidecar.is_file():
        previous = _read_json(sidecar)
        if (
            previous.get("status") == "staged"
            and previous.get("reference", {}).get("assembly_accession")
            == row["assembly_accession_version"]
        ):
            required = [previous["reference"][key] for key in ("fasta", "fai", "contig_dictionary")]
            if previous.get("annotation"):
                required.extend(previous["annotation"][key] for key in ("file", "contig_mapping"))
            if all(
                Path(str(artifact["uri"])[7:]).is_file()
                and sha256_file(Path(str(artifact["uri"])[7:])) == artifact["sha256"]
                for artifact in required
            ):
                return previous
    urls = row["artifact_urls"]
    started = time.monotonic()
    try:
        checksums = dataset_dir / "md5checksums.txt"
        report = dataset_dir / "assembly_report.txt"
        _download_unfrozen(urls["provider_checksums"], checksums)
        _download_unfrozen(urls["assembly_report"], report)
        md5s = _provider_md5s(checksums)
        expected_report_md5 = md5s.get(Path(urls["assembly_report"]).name)
        if not expected_report_md5 or _md5(report) != expected_report_md5:
            raise Tier3ValidationError("NCBI assembly-report provider MD5 mismatch or absent checksum")
        raw_fasta = dataset_dir / Path(urls["fasta"]).name
        _download_unfrozen(urls["fasta"], raw_fasta)
        expected_fasta_md5 = md5s.get(raw_fasta.name)
        if not expected_fasta_md5 or _md5(raw_fasta) != expected_fasta_md5:
            raise Tier3ValidationError("NCBI FASTA provider MD5 mismatch or absent checksum")
        excluded, exclusion_counts = _excluded_non_nuclear(report)
        fasta = dataset_dir / "nuclear.fna"
        retained_contigs = _filter_fasta(raw_fasta, fasta, excluded)
        fasta = _bgzip_replace(fasta)
        subprocess.run(["samtools", "faidx", str(fasta)], check=True)
        fai = Path(str(fasta) + ".fai")
        nuclear_reference_bases = _validate_nuclear_reference_size(fai)
        dictionary = dataset_dir / "nuclear.contigs.tsv"
        with dictionary.open("w", encoding="utf-8") as handle:
            handle.write("contig\tlength\n")
            for raw in fai.read_text(encoding="utf-8").splitlines():
                fields = raw.split("\t")
                handle.write(f"{fields[0]}\t{fields[1]}\n")

        annotation = None
        raw_gff = dataset_dir / Path(urls["gff"]).name
        try:
            _download_unfrozen(urls["gff"], raw_gff)
            expected_gff_md5 = md5s.get(raw_gff.name)
            if not expected_gff_md5 or _md5(raw_gff) != expected_gff_md5:
                raise Tier3ValidationError("NCBI GFF provider MD5 mismatch or absent checksum")
            gff = dataset_dir / "nuclear.gff3"
            regions, feature_count = _filter_gff(raw_gff, gff, excluded)
            if not regions or feature_count == 0:
                raise Tier3ValidationError("native NCBI GFF has no retained nuclear annotation")
            fasta_lengths = {
                fields[0]: int(fields[1])
                for fields in (line.split("\t") for line in fai.read_text(encoding="utf-8").splitlines())
            }
            if any(fasta_lengths.get(name) != length for name, length in regions):
                raise Tier3ValidationError("native GFF sequence-region dictionary disagrees with nuclear FASTA")
            mapping = dataset_dir / "annotation-contig-map.tsv"
            with mapping.open("w", encoding="utf-8") as handle:
                handle.write("annotation_contig\tfasta_contig\n")
                for name, _length in regions:
                    handle.write(f"{name}\t{name}\n")
            gff = _bgzip_replace(gff)
            annotation = {
                "provider": row["provider"],
                "release": row["release"],
                "assembly_accession": row["assembly_accession_version"],
                "status": "native",
                "native_vs_projected": "native",
                "format": "gff3",
                "file": _artifact(gff, source_uri=urls["gff"], upstream_md5=expected_gff_md5),
                "contig_mapping": _artifact(mapping),
                "genetic_code": 1,
                "exact_reference_assertion": True,
                "sequence_region_count": len(regions),
                "feature_count": feature_count,
            }
        except urllib.error.HTTPError as error:
            if error.code != 404:
                raise

        # The checksum-locked, nuclear-filtered derivatives are the analysis
        # inputs.  Provider archives are no longer needed once their MD5 and
        # source URI have been recorded, so remove them to bound scratch use.
        raw_fasta.unlink(missing_ok=True)
        raw_gff.unlink(missing_ok=True)
        value = {
            "schema_version": STAGE_SCHEMA,
            "status": "staged",
            "dataset_id": row["dataset_id"],
            "species": {
                "scientific_name": row["scientific_name"],
                "taxon_id": row["taxon_id"],
                "is_congener_substitution": False,
            },
            "buffalo_diversity": row["buffalo_diversity"],
            "buffalo_pred_log10_N": row["buffalo_pred_log10_N"],
            "reference": {
                "assembly_accession": row["assembly_accession_version"],
                "provider": row["provider"],
                "release": row["release"],
                "fasta": _artifact(fasta, source_uri=urls["fasta"], upstream_md5=expected_fasta_md5),
                "fai": _artifact(fai),
                "contig_dictionary": _artifact(dictionary),
                "nuclear_contigs_only": True,
                "assembly_report": _artifact(
                    report, source_uri=urls["assembly_report"], upstream_md5=expected_report_md5
                ),
                "provider_checksums": _artifact(checksums, source_uri=urls["provider_checksums"]),
                "non_nuclear_exclusion_counts": exclusion_counts,
                "retained_contigs": retained_contigs,
                "retained_nuclear_bases": nuclear_reference_bases,
            },
            "annotation": annotation,
            "staging": {
                "host": socket.getfqdn(),
                "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
                "slurm_array_task_id": os.environ.get("SLURM_ARRAY_TASK_ID"),
                "wall_seconds": round(time.monotonic() - started, 6),
                "max_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                "staged_bytes": sum(path.stat().st_size for path in dataset_dir.iterdir() if path.is_file()),
                "environment": _environment_from_record(environment_record),
            },
        }
        atomic_write_json(value, sidecar)
        return value
    except Exception as error:
        failure = {
            "schema_version": STAGE_SCHEMA,
            "status": "failed",
            "dataset_id": row["dataset_id"],
            "scientific_name": row["scientific_name"],
            "reason_code": getattr(error, "reason_code", "staging_or_native_annotation_audit_failed"),
            "detail": f"{type(error).__name__}: {error}",
            "wall_seconds": round(time.monotonic() - started, 6),
            "max_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
            "slurm_array_task_id": os.environ.get("SLURM_ARRAY_TASK_ID"),
        }
        atomic_write_json(failure, sidecar)
        raise


def _environment_from_record(path: Path) -> Dict[str, Any]:
    record = _read_json(path)
    tools = {
        name: {
            "version": record["tool_versions"][name],
            "executable": executable,
            "store_path": executable.split("/bin/", 1)[0],
        }
        for name, executable in record["tool_store_paths"].items()
    }
    return {
        "manager": "gnu-guix",
        "guix_environment": record["profile_store_path"],
        "channel_commit": record["channel_commit"],
        "resolved_channels_sha256": record["resolved_channels_sha256"],
        "manifest_sha256": record["manifest_sha256"],
        "profile_store_path": record["profile_store_path"],
        "store_paths": record["store_paths"],
        "tools": tools,
    }


def freeze(discovery_path: Path, stage_root: Path, environment_record: Path, output: Path) -> Mapping[str, Any]:
    discovery_value = _read_json(discovery_path)
    datasets: List[Mapping[str, Any]] = []
    failures: List[Dict[str, Any]] = []
    selected_indices = {
        row["dataset_id"]: index
        for index, row in enumerate(sorted(
            (item for item in discovery_value["datasets"] if item["status"] == "selected"),
            key=lambda item: item["dataset_id"].encode("utf-8"),
        ))
    }
    for row in discovery_value["datasets"]:
        if row["status"] != "selected":
            failures.append({
                "scientific_name": row["scientific_name"], "phase": "discovery",
                "reason_code": row.get("reason_code", row["status"]), "detail": row.get("detail", ""),
                "reproducible_command": (
                    "python3 analysis/tier3c_ncbi_gc.py discover "
                    + shlex.quote(row["scientific_name"])
                ),
            })
            continue
        sidecar = stage_root / row["dataset_id"] / "stage.json"
        if not sidecar.is_file():
            failures.append({
                "scientific_name": row["scientific_name"], "phase": "staging",
                "reason_code": "missing_stage_record", "detail": str(sidecar),
                "reproducible_command": " ".join(map(shlex.quote, (
                    "python3", "analysis/tier3c_batch.py", "stage-one", str(discovery_path),
                    str(selected_indices[row["dataset_id"]]), str(stage_root), str(environment_record),
                ))),
            })
            continue
        staged = _read_json(sidecar)
        if staged.get("status") != "staged":
            failures.append({
                "scientific_name": row["scientific_name"], "phase": "staging",
                "reason_code": staged.get("reason_code", "stage_failed"), "detail": staged.get("detail", ""),
                "reproducible_command": " ".join(map(shlex.quote, (
                    "python3", "analysis/tier3c_batch.py", "stage-one", str(discovery_path),
                    str(selected_indices[row["dataset_id"]]), str(stage_root), str(environment_record),
                ))),
            })
            continue
        for group, keys in (
            (staged["reference"], ("fasta", "fai", "contig_dictionary", "assembly_report", "provider_checksums")),
            (staged.get("annotation") or {}, ("file", "contig_mapping")),
        ):
            for key in keys:
                artifact = group.get(key)
                if artifact:
                    verify_file(Path(str(artifact["uri"])[7:]), artifact["sha256"], artifact["size_bytes"])
        datasets.append(staged)
    value = {
        "schema_version": BATCH_SCHEMA,
        "decision_version": "tier3-decisions-v1",
        "discovery_sha256": sha256_file(discovery_path),
        "environment": _environment_from_record(environment_record),
        "datasets": sorted(datasets, key=lambda row: row["dataset_id"].encode("utf-8")),
        "pre_analysis_failures": sorted(failures, key=lambda row: row["scientific_name"].encode("utf-8")),
    }
    atomic_write_json(value, output)
    return value


def run_one(batch_path: Path, index: int, output_root: Path) -> Mapping[str, Any]:
    batch = _read_json(batch_path)
    if batch.get("schema_version") != BATCH_SCHEMA:
        raise Tier3ValidationError("run-one requires a frozen Tier 3c batch manifest")
    datasets = batch["datasets"]
    if not 0 <= index < len(datasets):
        raise Tier3ValidationError(f"array index {index} outside frozen dataset count {len(datasets)}")
    dataset = datasets[index]
    output = output_root / f"{dataset['dataset_id']}.json"
    started = time.monotonic()
    failure_path = output_root / f"{dataset['dataset_id']}.failure.json"
    try:
        result = analyze_dataset(dataset, output, environment=batch["environment"])
    except Exception as error:
        # A failed rerun must never leave an older successful result available
        # for collection under the same frozen dataset ID.
        output.unlink(missing_ok=True)
        failure = {
            "schema_version": "tier3c-analysis-failure-v1",
            "dataset_id": dataset["dataset_id"],
            "scientific_name": dataset["species"]["scientific_name"],
            "reason_code": "analysis_validation_failed",
            "detail": f"{type(error).__name__}: {error}",
            "host": socket.getfqdn(),
            "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
            "slurm_array_task_id": os.environ.get("SLURM_ARRAY_TASK_ID"),
            "wall_seconds": round(time.monotonic() - started, 6),
            "max_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            "environment": batch["environment"],
        }
        atomic_write_json(failure, failure_path)
        raise
    failure_path.unlink(missing_ok=True)
    job_record = {
        "schema_version": "tier3c-job-v1",
        "dataset_id": dataset["dataset_id"],
        "result_sha256": sha256_file(output),
        "host": socket.getfqdn(),
        "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
        "slurm_array_task_id": os.environ.get("SLURM_ARRAY_TASK_ID"),
        "wall_seconds": round(time.monotonic() - started, 6),
        "max_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        "input_bytes": sum(
            artifact["size_bytes"]
            for artifact in (
                dataset["reference"]["fasta"],
                *((dataset["annotation"][key] for key in ("file", "contig_mapping")) if dataset.get("annotation") else ()),
            )
        ),
        "environment": batch["environment"],
    }
    atomic_write_json(job_record, output_root / f"{dataset['dataset_id']}.job.json")
    return result


def collect(batch_path: Path, result_root: Path, output_dir: Path) -> None:
    batch = _read_json(batch_path)
    result_rows: List[Dict[str, Any]] = []
    failures = list(batch.get("pre_analysis_failures", []))
    qc_records: List[Mapping[str, Any]] = []
    provenance_dir = output_dir / "tier3c_qc"
    provenance_dir.mkdir(parents=True, exist_ok=True)
    manifest_rows: List[Dict[str, Any]] = []
    for index, dataset in enumerate(batch["datasets"]):
        reference, annotation = dataset["reference"], dataset.get("annotation")
        manifest_rows.append({
            "dataset_id": dataset["dataset_id"],
            "scientific_name": dataset["species"]["scientific_name"],
            "taxon_id": dataset["species"]["taxon_id"],
            "assembly_accession_version": reference["assembly_accession"],
            "provider": reference["provider"],
            "release": reference["release"],
            "fasta_source_uri": reference["fasta"].get("source_uri", ""),
            "fasta_upstream_md5": reference["fasta"].get("upstream_md5", ""),
            "fasta_sha256": reference["fasta"]["sha256"],
            "fasta_size_bytes": reference["fasta"]["size_bytes"],
            "gff_source_uri": annotation["file"].get("source_uri", "") if annotation else "",
            "gff_upstream_md5": annotation["file"].get("upstream_md5", "") if annotation else "",
            "gff_sha256": annotation["file"]["sha256"] if annotation else "",
            "gff_size_bytes": annotation["file"]["size_bytes"] if annotation else "",
            "annotation_native_exact_reference": bool(
                annotation and annotation.get("status") == "native"
                and annotation.get("assembly_accession") == reference["assembly_accession"]
            ),
            "genetic_code": annotation.get("genetic_code", "") if annotation else "",
        })
        result_path = result_root / f"{dataset['dataset_id']}.json"
        command = f"python3 analysis/tier3c_batch.py run-one {batch_path} {index} {result_root}"
        if not result_path.is_file():
            failure_path = result_root / f"{dataset['dataset_id']}.failure.json"
            failure_record = _read_json(failure_path) if failure_path.is_file() else {}
            failures.append({
                "scientific_name": dataset["species"]["scientific_name"], "phase": "analysis",
                "reason_code": failure_record.get("reason_code", "missing_result"),
                "detail": failure_record.get("detail", str(result_path)),
                "reproducible_command": command,
            })
            continue
        result = _read_json(result_path)
        gc3, genome = result["gc3"], result["whole_genome_gc"]
        result_rows.append({
            "dataset_id": result["dataset_id"],
            "scientific_name": result["species"]["scientific_name"],
            "taxon_id": result["species"]["taxon_id"],
            "buffalo_diversity": dataset["buffalo_diversity"],
            "buffalo_pred_log10_N": dataset["buffalo_pred_log10_N"],
            "assembly_accession_version": result["reference"]["accession"],
            "whole_genome_gc": genome.get("value", ""),
            "whole_genome_callable_bases": genome.get("callable_bases", ""),
            "gc3": gc3.get("value", ""),
            "gc3_callable_third_positions": gc3.get("callable_third_positions", ""),
            "callable_genes": gc3.get("genes", ""),
            "annotation_status": (result.get("annotation_provenance") or {}).get("status", "unavailable"),
        })
        qc = {
            "dataset_id": result["dataset_id"],
            "species": result["species"],
            "reference": result["reference"],
            "annotation_provenance": result.get("annotation_provenance"),
            "whole_genome_gc": genome,
            "gc3": gc3,
            "environment": result["environment"],
            "result_sha256": sha256_file(result_path),
            "job": _read_json(result_root / f"{dataset['dataset_id']}.job.json")
            if (result_root / f"{dataset['dataset_id']}.job.json").is_file() else None,
            "pilot_failures": validate_pilot(result["species"]["scientific_name"], result)
            if result["species"]["scientific_name"] in {"Drosophila melanogaster", "Homo sapiens"} else [],
        }
        atomic_write_json(qc, provenance_dir / f"{result['dataset_id']}.json")
        qc_records.append(qc)
    for failure in failures:
        failure.setdefault("reproducible_command", "python3 analysis/tier3c_batch.py discover <buffalo.tsv> <discovery.json>")
    _atomic_tsv(output_dir / "tier3c_data.tsv", RESULT_COLUMNS, sorted(result_rows, key=lambda row: row["dataset_id"]))
    _atomic_tsv(
        output_dir / "tier3c_manifest.tsv", MANIFEST_COLUMNS,
        sorted(manifest_rows, key=lambda row: row["dataset_id"]),
    )
    _atomic_tsv(output_dir / "tier3c_failure_ledger.tsv", FAILURE_COLUMNS, sorted(failures, key=lambda row: row["scientific_name"]))
    values = [row["gc3"]["value"] for row in qc_records if row["gc3"].get("status") == "available"]
    genes = [row["gc3"]["genes"] for row in qc_records if row["gc3"].get("status") == "available"]
    controls = {
        row["species"]["scientific_name"]: {
            "passed": not row["pilot_failures"],
            "failures": row["pilot_failures"],
        }
        for row in qc_records
        if row["species"]["scientific_name"] in {"Drosophila melanogaster", "Homo sapiens"}
    }
    exclusion_totals: Dict[str, int] = {}
    exclusion_by_species: List[Dict[str, Any]] = []
    unavailable_reasons: Dict[str, int] = {}
    for row in qc_records:
        exclusions = ((row.get("annotation_provenance") or {}).get("cds_audit") or {}).get(
            "exclusions", {}
        )
        for reason, count in exclusions.items():
            exclusion_totals[reason] = exclusion_totals.get(reason, 0) + int(count)
        exclusion_by_species.append({
            "scientific_name": row["species"]["scientific_name"],
            "total_exclusions": sum(int(count) for count in exclusions.values()),
            "exclusions": dict(sorted(exclusions.items())),
        })
        if row["gc3"].get("status") != "available":
            reason = row["gc3"].get("reason", "unspecified")
            unavailable_reasons[reason] = unavailable_reasons.get(reason, 0) + 1
    def outliers(field: str) -> List[str]:
        observed = sorted(
            (math.log10(float(row[field])), row["scientific_name"])
            for row in result_rows if row[field] not in {"", 0}
        )
        if len(observed) < 4:
            return []
        values_only = [value for value, _name in observed]
        def quantile(fraction: float) -> float:
            position = fraction * (len(values_only) - 1)
            lower = int(position)
            upper = min(lower + 1, len(values_only) - 1)
            return values_only[lower] + (position - lower) * (values_only[upper] - values_only[lower])
        q1, q3 = quantile(0.25), quantile(0.75)
        lower, upper = q1 - 1.5 * (q3 - q1), q3 + 1.5 * (q3 - q1)
        return sorted(name for value, name in observed if value < lower or value > upper)

    gene_outlier_names = outliers("callable_genes")
    result_by_species = {row["scientific_name"]: row for row in result_rows}
    exclusion_record_by_species = {row["scientific_name"]: row for row in exclusion_by_species}
    review_notes = {
        "Malurus melanocephalus": (
            "retain with annotation-quality sensitivity flag: low callable-gene count co-occurs "
            "with extensive invalid/pseudogene CDS exclusions"
        ),
        "Oncorhynchus tshawytscha": (
            "retain flagged native result: high callable-gene count is biologically consistent "
            "with salmonid whole-genome duplication"
        ),
    }
    gene_outlier_review = []
    for name in gene_outlier_names:
        result_row = result_by_species[name]
        exclusion_record = exclusion_record_by_species[name]
        gene_outlier_review.append({
            "scientific_name": name,
            "callable_genes": result_row["callable_genes"],
            "gc3_callable_third_positions": result_row["gc3_callable_third_positions"],
            "total_annotation_exclusions": exclusion_record["total_exclusions"],
            "exclusions": exclusion_record["exclusions"],
            "review_disposition": review_notes.get(
                name, "retain native result with explicit outlier flag for synthesis sensitivity analysis"
            ),
        })

    job_telemetry = [row["job"] for row in qc_records if row.get("job")]
    resource_telemetry: Dict[str, Any] = {"jobs_recorded": len(job_telemetry)}
    if job_telemetry:
        maximum_rss_job = max(job_telemetry, key=lambda row: row["max_rss_kib"])
        maximum_wall_job = max(job_telemetry, key=lambda row: row["wall_seconds"])
        resource_telemetry.update({
            "historical_analysis_profile": {
                "requested_cpus_per_task": 2,
                "requested_memory_gib": 12,
                "time_limit_seconds": 3600,
                "scope": "completed 2026-07-13 frozen batch; retained as historical telemetry",
            },
            "retry_standard_profile": {
                "requested_cpus_per_task": 2,
                "requested_memory_gib": 32,
                "time_limit_seconds": 7200,
                "maximum_array_concurrency_per_node": 8,
            },
            "retry_outlier_profile": {
                "requested_cpus_per_task": 2,
                "requested_memory_gib": 64,
                "time_limit_seconds": 14400,
                "maximum_array_concurrency_per_node": 1,
            },
            "input_bytes_total": sum(row["input_bytes"] for row in job_telemetry),
            "max_rss_kib_median": int(statistics.median(row["max_rss_kib"] for row in job_telemetry)),
            "max_rss_kib_maximum": maximum_rss_job["max_rss_kib"],
            "max_rss_dataset": maximum_rss_job["dataset_id"],
            "wall_seconds_median": statistics.median(row["wall_seconds"] for row in job_telemetry),
            "wall_seconds_maximum": maximum_wall_job["wall_seconds"],
            "walltime_dataset": maximum_wall_job["dataset_id"],
            "resource_profile_scientific_role": (
                "scheduler headroom only; memory/time capacity is not evidence for estimator or control validity"
            ),
        })

    summary = {
        "schema_version": "tier3c-qc-summary-v1",
        "batch_manifest_sha256": sha256_file(batch_path),
        "eligible_staged": len(batch["datasets"]),
        "completed": len(result_rows),
        "failures": len(failures),
        "native_gc3_available": len(values),
        "native_gc3_unavailable": len(result_rows) - len(values),
        "gc3_range": [min(values), max(values)] if values else None,
        "callable_gene_range": [min(genes), max(genes)] if genes else None,
        "callable_base_outliers_reviewed": outliers("whole_genome_callable_bases"),
        "callable_gene_outliers_reviewed": gene_outlier_names,
        "callable_gene_outlier_review": gene_outlier_review,
        "annotation_exclusion_totals": dict(sorted(exclusion_totals.items())),
        "annotation_quality_distribution_review": {
            "gc3_unavailable_reasons": dict(sorted(unavailable_reasons.items())),
            "highest_exclusion_species": sorted(
                exclusion_by_species,
                key=lambda row: (-row["total_exclusions"], row["scientific_name"].encode("utf-8")),
            )[:10],
            "review_disposition": (
                "retain exact-native results; carry exclusion counts and missing-native status into "
                "synthesis sensitivity/weighting, with no projected annotation substituted"
            ),
        },
        "resource_telemetry": resource_telemetry,
        "pilot_controls": controls,
        "control_gate_passed": set(controls) == {"Drosophila melanogaster", "Homo sapiens"}
        and all(record["passed"] for record in controls.values()),
        "outlier_review_rule": "flag Tukey 1.5-IQR on log10 callable bases and callable genes during synthesis",
    }
    atomic_write_json(summary, output_dir / "tier3c_qc_summary.json")
    validate_collected(output_dir)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    command = commands.add_parser("discover")
    command.add_argument("buffalo", type=Path)
    command.add_argument("output", type=Path)
    command.add_argument("--delay-seconds", type=float, default=0.36)
    command = commands.add_parser("stage-one")
    command.add_argument("discovery", type=Path)
    command.add_argument("index", type=int)
    command.add_argument("stage_root", type=Path)
    command.add_argument("environment_record", type=Path)
    command = commands.add_parser("freeze")
    command.add_argument("discovery", type=Path)
    command.add_argument("stage_root", type=Path)
    command.add_argument("environment_record", type=Path)
    command.add_argument("output", type=Path)
    command = commands.add_parser("run-one")
    command.add_argument("batch", type=Path)
    command.add_argument("index", type=int)
    command.add_argument("output_root", type=Path)
    command = commands.add_parser("collect")
    command.add_argument("batch", type=Path)
    command.add_argument("result_root", type=Path)
    command.add_argument("output_dir", type=Path)
    command = commands.add_parser("validate-collected")
    command.add_argument("output_dir", type=Path)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _parser().parse_args(argv)
    if args.command == "discover":
        discover(args.buffalo, args.output, delay_seconds=args.delay_seconds)
    elif args.command == "stage-one":
        stage_one(args.discovery, args.index, args.stage_root, args.environment_record)
    elif args.command == "freeze":
        freeze(args.discovery, args.stage_root, args.environment_record, args.output)
    elif args.command == "run-one":
        run_one(args.batch, args.index, args.output_root)
    elif args.command == "collect":
        collect(args.batch, args.result_root, args.output_dir)
    elif args.command == "validate-collected":
        print(json.dumps(validate_collected(args.output_dir), sort_keys=True))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (Tier3ValidationError, subprocess.CalledProcessError, OSError) as error:
        print(f"tier3c-batch: {error}", file=sys.stderr)
        raise SystemExit(2)
