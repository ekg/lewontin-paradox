#!/usr/bin/env python3
"""Run fail-closed, atomic Tier 3 datasets inside the pinned Guix closure.

The scheduler submits one dataset per array element.  A dataset is a short,
declarative list of argv vectors; commands are never evaluated by a shell.
Each stage writes into a private directory, is validated, receives a checksum
marker, and is renamed atomically.  Resubmission verifies and reuses completed
stages while deleting abandoned ``.partial-*`` directories.
"""

from __future__ import annotations

import argparse
import fcntl
import hashlib
import json
import os
import random
import re
import shutil
import socket
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, MutableMapping, Optional, Sequence, Tuple

from analysis.tier3_common import (
    Tier3ValidationError,
    collect_fourfold_sites,
    fasta_dictionary,
    parse_gff,
    read_fasta,
    resolve_contig_aliases,
    sha256_file,
    traverse_paf,
)


ROOT = Path(__file__).resolve().parents[1]
CHANNELS = ROOT / "analysis/guix/channels.scm"
MANIFEST = ROOT / "analysis/guix/manifest.scm"
GUIX_LOAD_PATH = ROOT / "analysis/guix"
DECISION_VERSION = "tier3-decisions-v1"
CHANNEL_COMMIT = "44bbfc24e4bcc48d0e3343cd3d83452721af8c36"
STORE_PATH = re.compile(r"^/gnu/store/[0-9a-z]{32}-[^/]+$")
STORE_DERIVATION = re.compile(r"^/gnu/store/[0-9a-z]{32}-[^/]+\.drv$")
DATASET_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$")
VARIABLE = re.compile(r"\$\{([A-Z][A-Z0-9_]*)\}")
TOKEN = re.compile(r"\{(stage_dir|dataset_dir|repo_root|workflow_dir|input:[A-Za-z0-9_.-]+|stage:[A-Za-z0-9_.-]+)\}")

MODE_STAGES: Mapping[str, Tuple[str, ...]] = {
    "composition": ("annotation_4d", "qc"),
    "population_vcf": ("mapping", "normalized_bcf", "annotation_4d", "qc"),
    "deposited_vcf": ("mapping", "normalized_bcf", "annotation_4d", "qc"),
    "direct_wfmash": ("alignment", "mapping", "normalized_bcf", "annotation_4d", "qc"),
}
TIER_MODES = {
    "3a": {"deposited_vcf", "direct_wfmash"},
    "3b": {"population_vcf"},
    "3c": {"composition"},
}
FORBIDDEN_ARGUMENTS = re.compile(
    r"(?i)(?:^|[/ _-])(conda|micromamba|mamba)(?:$|[/ _-])|source\s+activate|pip\s+--user"
)
SHELL_PROGRAMS = frozenset(("sh", "bash", "dash", "zsh", "fish", "csh", "tcsh"))
PRESERVED_EXACT = frozenset(("SCRATCH", "TMPDIR"))
PRESERVED_PREFIXES = ("SLURM_", "TIER3_")


class WorkflowError(Tier3ValidationError):
    """A workflow is unsafe, incomplete, inconsistent, or non-resumable."""


def _utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def _canonical(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("utf-8")


def _digest(value: Any) -> str:
    return hashlib.sha256(_canonical(value)).hexdigest()


def _atomic_json(path: Path, value: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".partial-{os.getpid()}")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def _expand_declared_variables(value: str, environment: Optional[Mapping[str, str]] = None) -> str:
    environment = os.environ if environment is None else environment

    def replace(match: re.Match[str]) -> str:
        name = match.group(1)
        if not (name.startswith(PRESERVED_PREFIXES) or name in PRESERVED_EXACT):
            raise WorkflowError(f"workflow references undeclared environment variable {name}")
        if name not in environment:
            raise WorkflowError(f"workflow requires unset environment variable {name}")
        return environment[name]

    rendered = VARIABLE.sub(replace, value)
    if "$" in rendered:
        raise WorkflowError(f"unexpanded or non-declarative variable syntax in {value!r}")
    return rendered


def _safe_path(value: Any, label: str) -> Path:
    if not isinstance(value, str) or not value:
        raise WorkflowError(f"{label} must be a non-empty path string")
    return Path(_expand_declared_variables(value)).expanduser().resolve()


def _validate_command(stage: Mapping[str, Any], dataset_id: str) -> None:
    argv = stage.get("argv")
    if not isinstance(argv, list) or not argv or not all(isinstance(item, str) and item for item in argv):
        raise WorkflowError(f"dataset {dataset_id}: stage {stage.get('name')!r} requires a non-empty argv list")
    program = Path(argv[0]).name.lower()
    if program in SHELL_PROGRAMS or any(item in {"-c", "--command"} for item in argv[1:2] if program in SHELL_PROGRAMS):
        raise WorkflowError(f"dataset {dataset_id}: shell execution is forbidden; provide an argv vector")
    rendered = " ".join(argv)
    if FORBIDDEN_ARGUMENTS.search(rendered):
        raise WorkflowError(f"dataset {dataset_id}: Conda or another forbidden environment manager appears in argv")
    outputs = stage.get("outputs")
    if not isinstance(outputs, list) or not outputs or not all(isinstance(item, str) and item for item in outputs):
        raise WorkflowError(f"dataset {dataset_id}: stage {stage.get('name')!r} requires declared outputs")
    for output in outputs:
        candidate = Path(output)
        if candidate.is_absolute() or ".." in candidate.parts or output.startswith("."):
            raise WorkflowError(f"dataset {dataset_id}: stage output must be a safe relative path: {output!r}")


def _validate_workflow_value(value: Any, path: Path) -> Dict[str, Any]:
    if not isinstance(value, dict):
        raise WorkflowError("workflow root must be an object")
    if value.get("schema_version") != "1.0":
        raise WorkflowError("workflow schema_version must be 1.0")
    if value.get("decision_version") != DECISION_VERSION:
        raise WorkflowError(f"workflow decision_version must be {DECISION_VERSION}")
    _safe_path(value.get("scratch_root"), "scratch_root")
    _safe_path(value.get("output_root"), "output_root")
    datasets = value.get("datasets")
    if not isinstance(datasets, list) or not datasets:
        raise WorkflowError("workflow requires at least one dataset")
    observed: set[str] = set()
    for dataset in datasets:
        if not isinstance(dataset, dict):
            raise WorkflowError("each workflow dataset must be an object")
        dataset_id = dataset.get("dataset_id")
        if not isinstance(dataset_id, str) or not DATASET_ID.fullmatch(dataset_id):
            raise WorkflowError(f"unsafe dataset_id {dataset_id!r}")
        if dataset_id in observed:
            raise WorkflowError(f"duplicate dataset_id {dataset_id!r}")
        observed.add(dataset_id)
        tier, mode = dataset.get("tier"), dataset.get("mode")
        if tier not in TIER_MODES or mode not in TIER_MODES[tier]:
            raise WorkflowError(f"dataset {dataset_id}: mode {mode!r} is invalid for tier {tier!r}")
        stages = dataset.get("stages")
        if not isinstance(stages, list):
            raise WorkflowError(f"dataset {dataset_id}: stages must be a list")
        names = tuple(item.get("name") if isinstance(item, dict) else None for item in stages)
        if names != MODE_STAGES[mode]:
            raise WorkflowError(
                f"dataset {dataset_id}: stage order must be {list(MODE_STAGES[mode])!r}, observed {list(names)!r}"
            )
        for stage in stages:
            _validate_command(stage, dataset_id)
    result = dict(value)
    result["_workflow_path"] = str(path.resolve())
    return result


def load_workflow(path: Path | str) -> Dict[str, Any]:
    workflow_path = Path(path).expanduser().resolve()
    try:
        value = json.loads(workflow_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise WorkflowError(f"cannot read workflow {workflow_path}: {error}") from error
    return _validate_workflow_value(value, workflow_path)


def sanitized_stage_environment(stage_tmp: Path) -> Dict[str, str]:
    """Return a deterministic process environment with a narrow allow-list.

    PATH is deliberately copied because the caller is already inside the pure
    manifest shell.  HOME and locale are fresh fixed values, not preserved
    login state.  TMPDIR is always stage-local even if the submit environment
    supplied another value.
    """

    result = {
        "PATH": os.environ.get("PATH", os.defpath),
        "HOME": str((Path(stage_tmp) / ".home").resolve()),
        "LANG": "C.UTF-8",
        "LC_ALL": "C.UTF-8",
        "TZ": "UTC",
    }
    # These are derived from the already-audited profile by guix_job.sh.  They
    # are runtime mechanics, not preserved login/user environment state.
    for name in ("GUIX_PROFILE", "GUIX_PYTHONPATH"):
        if name in os.environ:
            result[name] = os.environ[name]
    for name, value in os.environ.items():
        if name.startswith(PRESERVED_PREFIXES) or name in PRESERVED_EXACT:
            result[name] = value
    result["TMPDIR"] = str(Path(stage_tmp).resolve())
    Path(result["HOME"]).mkdir(parents=True, exist_ok=True)
    return result


def _read_callable_mask(path: Path, reference: Mapping[str, str]) -> int:
    total = 0
    seen = False
    previous: Dict[str, Tuple[int, int]] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 3 or fields[0] not in reference:
                raise WorkflowError(f"callable mask has invalid contig at line {line_number}")
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError as error:
                raise WorkflowError(f"callable mask has non-integer coordinates at line {line_number}") from error
            if not 0 <= start < end <= len(reference[fields[0]]):
                raise WorkflowError(f"callable mask is outside the exact reference at line {line_number}")
            if fields[0] in previous and start < previous[fields[0]][1]:
                raise WorkflowError(f"callable mask overlaps or is unsorted at line {line_number}")
            previous[fields[0]] = (start, end)
            total += sum(base in "ACGT" for base in reference[fields[0]][start:end].upper())
            seen = True
    if not seen or total == 0:
        raise WorkflowError("callable mask has no unambiguous exact-reference denominator")
    return total


def _annotation_accession(provenance: Mapping[str, Any]) -> Optional[str]:
    return provenance.get("assembly_accession_version") or provenance.get("assembly_accession")


def _audit_native_annotation(
    dataset: Mapping[str, Any], preflight: Mapping[str, Any], reference_path: Path, reference: Mapping[str, str]
) -> Dict[str, Any]:
    provenance = preflight.get("annotation_provenance")
    gff_value = preflight.get("annotation_gff")
    if not isinstance(provenance, dict) or not gff_value:
        raise WorkflowError("primary annotation/4D stage requires native annotation provenance and GFF")
    required = ("provider", "release", "status", "genetic_code", "fasta_sha256", "gff_sha256", "contig_mapping")
    missing = [field for field in required if field not in provenance]
    if missing:
        raise WorkflowError(f"annotation provenance lacks fields: {missing!r}")
    if provenance["status"] != "native" or provenance.get("native_vs_projected", "native") != "native":
        raise WorkflowError("primary annotation-derived results require native, not projected, annotation")
    accession = _annotation_accession(provenance)
    if not accession or accession != dataset.get("reference_accession"):
        raise WorkflowError("annotation assembly accession/version does not match the exact reference assembly")
    fasta_accession = provenance.get("fasta_assembly_accession", accession)
    if fasta_accession != accession:
        raise WorkflowError("annotation and FASTA assembly accession/version disagree")
    if sha256_file(reference_path) != provenance["fasta_sha256"]:
        raise WorkflowError("annotation provenance FASTA checksum does not match exact reference")
    gff_path = _safe_path(gff_value, "annotation_gff")
    if not gff_path.is_file() or sha256_file(gff_path) != provenance["gff_sha256"]:
        raise WorkflowError("annotation provenance GFF checksum does not match native annotation")
    annotation = parse_gff(gff_path)
    aliases = provenance["contig_mapping"]
    if not isinstance(aliases, dict) or not all(isinstance(k, str) and isinstance(v, str) for k, v in aliases.items()):
        raise WorkflowError("annotation contig_mapping must be an explicit object")
    resolved = resolve_contig_aliases(fasta_dictionary(reference), annotation.sequence_regions, aliases)
    annotation_fasta = dict(reference)
    for annotation_contig, fasta_contig in resolved.items():
        annotation_fasta[annotation_contig] = reference[fasta_contig]
    sites, exclusions = collect_fourfold_sites(annotation_fasta, annotation, int(provenance["genetic_code"]))
    return {
        "provider": provenance["provider"],
        "release": provenance["release"],
        "assembly_accession_version": accession,
        "fasta_sha256": provenance["fasta_sha256"],
        "gff_sha256": provenance["gff_sha256"],
        "sequence_regions": dict(sorted(annotation.sequence_regions.items())),
        "contig_mapping": dict(sorted(resolved.items())),
        "genetic_code": int(provenance["genetic_code"]),
        "native_vs_projected": "native",
        "retained_fourfold_sites": len(sites),
        "overlap_exclusions": len(exclusions),
        "retained_cds_reconstruction_passed": True,
    }


def _audit_paf(path: Path, target: Mapping[str, str], query: Mapping[str, str]) -> Dict[str, Any]:
    lines = path.read_text(encoding="utf-8").splitlines()
    query_intervals: Dict[str, List[Tuple[int, int, str]]] = {}
    mapping_count = 0
    for line_number, line in enumerate(lines, 1):
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 12:
            raise WorkflowError(f"PAF line {line_number} has fewer than 12 columns")
        mapping_count += 1
        qname, qstart, qend, target_name = fields[0], int(fields[2]), int(fields[3]), fields[5]
        query_intervals.setdefault(qname, []).append((qstart, qend, target_name))
        try:
            mapq = int(fields[11])
        except ValueError as error:
            raise WorkflowError(f"PAF line {line_number} has invalid mapping quality") from error
        if mapq <= 0 or mapq == 255:
            raise WorkflowError(f"PAF line {line_number} is ambiguous rather than uniquely mapped")
    if not mapping_count:
        raise WorkflowError("PAF has no mappings")
    for qname, intervals in query_intervals.items():
        intervals.sort()
        for left, right in zip(intervals, intervals[1:]):
            if right[0] < left[1]:
                raise WorkflowError(f"PAF query {qname!r} has multiple/overlapping mappings")
    try:
        result = traverse_paf(lines, target, query, edge_exclusion_bp=0, indel_flank_bp=0)
    except Tier3ValidationError as error:
        raise WorkflowError(str(error)) from error
    if result.exclusion_counts.get("multiple_projection", 0):
        raise WorkflowError("PAF target bases have multiple mappings; uniqueness gate failed")
    if not result.callable_positions:
        raise WorkflowError("PAF has no uniquely mapped callable A/C/G/T bases")
    return {
        "mapping_records": mapping_count,
        "callable_bases_before_policy_edges": len(result.callable_positions),
        "operation_counts": result.operation_counts,
        "exclusion_counts": result.exclusion_counts,
        "extended_cigar_passed": True,
        "unique_mapping_passed": True,
    }


def _audit_variant_reference(path: Path, reference: Mapping[str, str]) -> Dict[str, Any]:
    try:
        import pysam
    except ImportError as error:
        raise WorkflowError("pysam is required from the pinned Guix manifest") from error
    records = 0
    try:
        with pysam.VariantFile(str(path)) as variants:
            dictionary = {name: variants.header.contigs[name].length for name in variants.header.contigs}
            if set(dictionary) != set(reference):
                raise WorkflowError("VCF/BCF contig dictionary does not exactly match reference FASTA")
            for contig, sequence in reference.items():
                if dictionary[contig] != len(sequence):
                    raise WorkflowError(f"VCF/BCF contig length mismatch for {contig!r}")
            for record in variants:
                observed = reference[record.contig][record.start : record.start + len(record.ref)].upper()
                if observed != record.ref.upper():
                    raise WorkflowError(
                        f"VCF/BCF REF mismatch at {record.contig}:{record.pos}: {record.ref} != {observed}"
                    )
                records += 1
    except (OSError, ValueError) as error:
        raise WorkflowError(f"cannot audit VCF/BCF {path}: {error}") from error
    return {"records": records, "contig_dictionary_passed": True, "ref_allele_audit_passed": True}


def validate_dataset_preflight(dataset: Mapping[str, Any]) -> Dict[str, Any]:
    dataset_id = dataset.get("dataset_id", "<unknown>")
    mode = dataset.get("mode")
    preflight = dataset.get("preflight")
    if not isinstance(preflight, dict):
        raise WorkflowError(f"dataset {dataset_id}: preflight must be an object")
    reference_value = preflight.get("reference_fasta")
    if not reference_value:
        raise WorkflowError(f"dataset {dataset_id}: exact reference FASTA is required")
    reference_path = _safe_path(reference_value, "reference_fasta")
    if not reference_path.is_file():
        raise WorkflowError(f"dataset {dataset_id}: reference FASTA is missing: {reference_path}")
    reference = read_fasta(reference_path)
    audit: Dict[str, Any] = {
        "dataset_id": dataset_id,
        "reference_fasta": str(reference_path),
        "reference_fasta_sha256": sha256_file(reference_path),
        "reference_dictionary": fasta_dictionary(reference),
    }

    annotation_disposition = dataset.get("annotation_derived", "primary_native")
    if annotation_disposition not in {"primary_native", "unavailable_missing_native"}:
        raise WorkflowError(
            f"dataset {dataset_id}: annotation_derived must be primary_native or unavailable_missing_native"
        )
    # Missing native annotation is an explicit unavailability, not permission
    # to use a lifted, congener, or predicted annotation as primary.
    if annotation_disposition == "primary_native":
        audit["annotation"] = _audit_native_annotation(dataset, preflight, reference_path, reference)
    else:
        if preflight.get("annotation_gff") or preflight.get("annotation_provenance"):
            raise WorkflowError(
                f"dataset {dataset_id}: unavailable primary annotation may not silently carry projected/native inputs"
            )
        audit["annotation"] = {
            "status": "unavailable_missing_native_exact_reference_annotation",
            "gc3": "unavailable",
            "fourfold": "unavailable",
            "whole_genome_gc": "eligible",
        }

    if mode in {"population_vcf", "deposited_vcf"}:
        mask_value = preflight.get("callable_mask")
        if not mask_value:
            raise WorkflowError(f"dataset {dataset_id}: an explicit callable mask is required")
        mask_path = _safe_path(mask_value, "callable_mask")
        if not mask_path.is_file():
            raise WorkflowError(f"dataset {dataset_id}: callable mask is missing: {mask_path}")
        audit["callable_mask"] = {
            "path": str(mask_path),
            "sha256": sha256_file(mask_path),
            "callable_acgt_bases": _read_callable_mask(mask_path, reference),
        }
    if mode == "population_vcf":
        samples_value = preflight.get("selected_samples")
        if not samples_value:
            raise WorkflowError(f"dataset {dataset_id}: population workflow requires an exact selected-sample list")
        samples_path = _safe_path(samples_value, "selected_samples")
        if not samples_path.is_file():
            raise WorkflowError(f"dataset {dataset_id}: selected-sample list is missing: {samples_path}")
        samples = [line.strip() for line in samples_path.read_text(encoding="utf-8").splitlines() if line.strip()]
        if not samples or len(samples) != len(set(samples)):
            raise WorkflowError(f"dataset {dataset_id}: selected-sample list is empty or contains duplicates")
        audit["selected_samples"] = {
            "path": str(samples_path),
            "sha256": sha256_file(samples_path),
            "sampling_units": len(samples),
        }

    variant_value = preflight.get("normalized_bcf") or preflight.get("variants")
    if variant_value:
        variant_path = _safe_path(variant_value, "variants")
        if not variant_path.is_file():
            raise WorkflowError(f"dataset {dataset_id}: declared VCF/BCF is missing: {variant_path}")
        audit["variants"] = {"path": str(variant_path), "sha256": sha256_file(variant_path), **_audit_variant_reference(variant_path, reference)}

    if mode == "direct_wfmash":
        if preflight.get("phase_identity_audit_passed") is not True:
            raise WorkflowError(f"dataset {dataset_id}: H1/H2 phase identity audit did not pass")
        if preflight.get("collapse_qc_passed") is not True:
            raise WorkflowError(f"dataset {dataset_id}: H1/H2 collapse/duplication QC did not pass")
        query_value, paf_value = preflight.get("query_fasta"), preflight.get("paf")
        if not query_value or not paf_value:
            raise WorkflowError(f"dataset {dataset_id}: direct WFMASH preflight requires H2 query FASTA and PAF")
        query_path, paf_path = _safe_path(query_value, "query_fasta"), _safe_path(paf_value, "paf")
        if not query_path.is_file() or not paf_path.is_file():
            raise WorkflowError(f"dataset {dataset_id}: H2 FASTA or PAF is missing")
        query = read_fasta(query_path)
        audit["query_fasta"] = {"path": str(query_path), "sha256": sha256_file(query_path)}
        audit["paf"] = {"path": str(paf_path), "sha256": sha256_file(paf_path), **_audit_paf(paf_path, reference, query)}
    return audit


def _render_argument(
    value: str,
    *,
    stage_tmp: Path,
    dataset_dir: Path,
    workflow: Mapping[str, Any],
    dataset: Mapping[str, Any],
    completed: Mapping[str, Path],
) -> str:
    preflight = dataset["preflight"]

    def replace(match: re.Match[str]) -> str:
        token = match.group(1)
        if token == "stage_dir":
            return str(stage_tmp)
        if token == "dataset_dir":
            return str(dataset_dir)
        if token == "repo_root":
            return str(ROOT)
        if token == "workflow_dir":
            return str(Path(workflow["_workflow_path"]).parent)
        kind, name = token.split(":", 1)
        if kind == "input":
            if name not in preflight or isinstance(preflight[name], (dict, list)):
                raise WorkflowError(f"command references absent/non-scalar preflight input {name!r}")
            return _expand_declared_variables(str(preflight[name]))
        if name not in completed:
            raise WorkflowError(f"command references incomplete stage {name!r}")
        return str(completed[name])

    rendered = TOKEN.sub(replace, value)
    if "{" in rendered or "}" in rendered:
        raise WorkflowError(f"unknown command placeholder in {value!r}")
    return _expand_declared_variables(rendered)


def _output_artifacts(stage_dir: Path, outputs: Sequence[str]) -> Dict[str, Dict[str, Any]]:
    result: Dict[str, Dict[str, Any]] = {}
    for relative in outputs:
        path = stage_dir / relative
        if not path.is_file() or path.stat().st_size == 0:
            raise WorkflowError(f"stage did not produce non-empty declared output {relative!r}")
        result[relative] = {"sha256": sha256_file(path), "size_bytes": path.stat().st_size}
    return result


def _validate_stage_outputs(
    stage: Mapping[str, Any], stage_dir: Path, dataset: Mapping[str, Any]
) -> Dict[str, Dict[str, Any]]:
    outputs: List[str] = stage["outputs"]
    artifacts = _output_artifacts(stage_dir, outputs)
    reference_path = _safe_path(dataset["preflight"]["reference_fasta"], "reference_fasta")
    reference = read_fasta(reference_path)
    for relative in outputs:
        path = stage_dir / relative
        if path.suffix == ".json":
            try:
                json.loads(path.read_text(encoding="utf-8"))
            except (UnicodeDecodeError, json.JSONDecodeError) as error:
                raise WorkflowError(f"declared JSON output is invalid: {relative}: {error}") from error
        if path.suffix == ".bed":
            _read_callable_mask(path, reference)
        if path.suffix == ".paf":
            query_value = dataset["preflight"].get("query_fasta")
            if not query_value:
                raise WorkflowError("cannot validate PAF output without query_fasta")
            _audit_paf(path, reference, read_fasta(_safe_path(query_value, "query_fasta")))
        if path.suffix == ".bcf":
            index_name = relative + ".csi"
            if index_name not in outputs:
                raise WorkflowError(f"normalized BCF output {relative!r} lacks a declared CSI index")
            _audit_variant_reference(path, reference)
        if path.name.endswith((".fa", ".fasta", ".fna")) and relative + ".fai" not in outputs:
            raise WorkflowError(f"FASTA output {relative!r} lacks a declared FAI index")
    return artifacts


def _marker_valid(stage_dir: Path, fingerprint: str) -> Optional[Dict[str, Any]]:
    marker_path = stage_dir / ".complete.json"
    if not marker_path.is_file():
        return None
    try:
        marker = json.loads(marker_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    if marker.get("fingerprint") != fingerprint or marker.get("status") != "complete":
        return None
    outputs = marker.get("outputs")
    if not isinstance(outputs, dict):
        return None
    for relative, expected in outputs.items():
        path = stage_dir / relative
        if not path.is_file() or path.stat().st_size != expected.get("size_bytes"):
            return None
        if sha256_file(path) != expected.get("sha256"):
            return None
    return marker


def _clean_partial_directories(parent: Path, prefix: str) -> None:
    for path in parent.glob(prefix + ".partial-*"):
        if path.is_dir():
            shutil.rmtree(path)
        else:
            path.unlink()


def _run_one_stage(
    workflow: Mapping[str, Any],
    dataset: Mapping[str, Any],
    stage: Mapping[str, Any],
    index: int,
    dataset_dir: Path,
    preflight_audit: Mapping[str, Any],
    completed: MutableMapping[str, Path],
    previous_markers: List[str],
) -> Dict[str, Any]:
    stage_name = stage["name"]
    final = dataset_dir / "stages" / f"{index:02d}-{stage_name}"
    final.parent.mkdir(parents=True, exist_ok=True)
    _clean_partial_directories(final.parent, final.name)
    fingerprint = _digest(
        {
            "decision_version": DECISION_VERSION,
            "dataset_id": dataset["dataset_id"],
            "stage": stage,
            "preflight": preflight_audit,
            "predecessor_markers": previous_markers,
        }
    )
    marker = _marker_valid(final, fingerprint)
    if marker is not None:
        completed[stage_name] = final
        previous_markers.append(sha256_file(final / ".complete.json"))
        return {"name": stage_name, "status": "complete", "resumed": True, "directory": str(final)}
    if final.exists():
        shutil.rmtree(final)
    prefix = final.name
    temporary = final.parent / f"{prefix}.partial-{os.getpid()}"
    temporary.mkdir(mode=0o700)
    argv = [
        _render_argument(
            item,
            stage_tmp=temporary,
            dataset_dir=dataset_dir,
            workflow=workflow,
            dataset=dataset,
            completed=completed,
        )
        for item in stage["argv"]
    ]
    environment = sanitized_stage_environment(temporary)
    completed_process = subprocess.run(
        argv,
        cwd=str(ROOT),
        env=environment,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )
    (temporary / ".stdout.log").write_text(completed_process.stdout, encoding="utf-8")
    (temporary / ".stderr.log").write_text(completed_process.stderr, encoding="utf-8")
    if completed_process.returncode:
        raise WorkflowError(
            f"dataset {dataset['dataset_id']} stage {stage_name} failed with exit {completed_process.returncode}; "
            f"partial directory retained at {temporary}"
        )
    outputs = _validate_stage_outputs(stage, temporary, dataset)
    marker = {
        "status": "complete",
        "decision_version": DECISION_VERSION,
        "dataset_id": dataset["dataset_id"],
        "stage": stage_name,
        "completed_at": _utc_now(),
        "fingerprint": fingerprint,
        "argv": argv,
        "outputs": outputs,
        "scheduler": {name: value for name, value in os.environ.items() if name.startswith("SLURM_")},
    }
    _atomic_json(temporary / ".complete.json", marker)
    os.replace(temporary, final)
    completed[stage_name] = final
    previous_markers.append(sha256_file(final / ".complete.json"))
    return {"name": stage_name, "status": "complete", "resumed": False, "directory": str(final)}


def _select_dataset(workflow: Mapping[str, Any], tier: str, array_index: int) -> Mapping[str, Any]:
    datasets = sorted(
        (item for item in workflow["datasets"] if item["tier"] == tier),
        key=lambda item: item["dataset_id"].encode("utf-8"),
    )
    if not 0 <= array_index < len(datasets):
        raise WorkflowError(f"array index {array_index} is outside tier {tier} dataset count {len(datasets)}")
    return datasets[array_index]


def run_array_task(workflow_path: Path | str, *, tier: str, array_index: int) -> Dict[str, Any]:
    workflow = load_workflow(workflow_path)
    dataset = _select_dataset(workflow, tier, array_index)
    preflight = validate_dataset_preflight(dataset)
    scratch_root = _safe_path(workflow["scratch_root"], "scratch_root")
    output_root = _safe_path(workflow["output_root"], "output_root")
    dataset_dir = scratch_root / tier / dataset["dataset_id"]
    dataset_dir.mkdir(parents=True, mode=0o700, exist_ok=True)
    lock_path = dataset_dir / ".workflow.lock"
    with lock_path.open("a+", encoding="utf-8") as lock:
        try:
            fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise WorkflowError(f"dataset {dataset['dataset_id']} is already running") from error
        completed: Dict[str, Path] = {}
        marker_hashes: List[str] = []
        stage_records = [
            _run_one_stage(workflow, dataset, stage, index, dataset_dir, preflight, completed, marker_hashes)
            for index, stage in enumerate(dataset["stages"])
        ]
        record = {
            "status": "complete",
            "decision_version": DECISION_VERSION,
            "dataset_id": dataset["dataset_id"],
            "tier": tier,
            "mode": dataset["mode"],
            "pilot": dataset.get("pilot"),
            "completed_at": _utc_now(),
            "workflow_sha256": sha256_file(Path(workflow["_workflow_path"])),
            "dataset_directory": str(dataset_dir),
            "preflight": preflight,
            "stages": stage_records,
            "final_marker_sha256": marker_hashes[-1],
        }
        _atomic_json(output_root / tier / f"{dataset['dataset_id']}.run.json", record)
        return record


def audit_environment_record(
    record: Mapping[str, Any], *, require_existing_profile: bool = True
) -> Mapping[str, Any]:
    if record.get("channel_commit") != CHANNEL_COMMIT:
        raise WorkflowError("environment record channel commit does not match channels.scm")
    profile = record.get("profile_store_path")
    if not isinstance(profile, str) or not STORE_PATH.fullmatch(profile):
        raise WorkflowError("environment record lacks a resolved profile store path")
    if require_existing_profile and not Path(profile).is_dir():
        raise WorkflowError(f"recorded profile store path is not visible: {profile}")
    gc_root = record.get("profile_gc_root")
    if not isinstance(gc_root, str) or not Path(gc_root).is_absolute():
        raise WorkflowError("environment record lacks an absolute profile GC root")
    if require_existing_profile and (
        not Path(gc_root).is_symlink() or str(Path(gc_root).resolve()) != profile
    ):
        raise WorkflowError("recorded profile GC root is absent or resolves to a different store path")
    derivations = record.get("derivations")
    if not isinstance(derivations, list) or not derivations or not all(
        isinstance(item, str) and STORE_DERIVATION.fullmatch(item) for item in derivations
    ):
        raise WorkflowError("environment record lacks valid derivations/derivation paths")
    store_paths = record.get("store_paths")
    if not isinstance(store_paths, list) or not store_paths or not all(
        isinstance(item, str) and STORE_PATH.fullmatch(item) for item in store_paths
    ):
        raise WorkflowError("environment record lacks valid store paths")
    versions = record.get("tool_versions")
    required_versions = {
        "python3",
        "pytest",
        "samtools",
        "bcftools",
        "bgzip",
        "tabix",
        "bedtools",
        "vcftools",
        "wfmash",
    }
    if (
        not isinstance(versions, dict)
        or not required_versions <= set(versions)
        or not all(isinstance(versions[name], str) and versions[name] for name in required_versions)
    ):
        raise WorkflowError("environment record lacks tool versions")
    resolved_channels = record.get("resolved_channels_scm")
    resolved_channels_sha256 = record.get("resolved_channels_sha256")
    if not isinstance(resolved_channels, str) or CHANNEL_COMMIT not in resolved_channels:
        raise WorkflowError("environment record lacks resolved channels at the frozen commit")
    if hashlib.sha256(resolved_channels.encode("utf-8")).hexdigest() != resolved_channels_sha256:
        raise WorkflowError("environment record resolved channels checksum is inconsistent")
    fallback = record.get("pack_fallback")
    if not isinstance(fallback, dict) or not isinstance(fallback.get("required"), bool):
        raise WorkflowError("environment record lacks pack fallback disposition")
    if fallback["required"]:
        for field in ("manifest_sha256", "runtime_store_path", "pack_sha256", "pack_path"):
            if not fallback.get(field):
                raise WorkflowError(f"required Guix pack fallback lacks {field}")
    return record


def _command_version(
    executable: Path,
    arguments: Sequence[str] = ("--version",),
    *,
    empty_fallback: Optional[str] = None,
) -> str:
    completed = subprocess.run(
        [str(executable), *arguments], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, check=False
    )
    rendered = completed.stdout.strip()
    if completed.returncode:
        raise WorkflowError(f"cannot record version for {executable}: exit {completed.returncode}")
    if not rendered and empty_fallback is not None:
        return empty_fallback
    if not rendered:
        raise WorkflowError(f"cannot record version for {executable}: command emitted no version")
    return rendered.splitlines()[0]


def record_environment(
    profile: Path,
    derivations_file: Path,
    store_paths_file: Path,
    resolved_channels_file: Path,
    output: Path,
) -> Dict[str, Any]:
    profile_store = str(profile.resolve(strict=True))
    if not STORE_PATH.fullmatch(profile_store):
        raise WorkflowError(f"profile does not resolve into /gnu/store: {profile_store}")
    derivations = sorted(set(line.strip() for line in derivations_file.read_text().splitlines() if line.strip()))
    store_paths = sorted(set(line.strip() for line in store_paths_file.read_text().splitlines() if line.strip()))
    tools = {
        "python3": ("python3", ("--version",)),
        "pytest": ("pytest", ("--version",)),
        "samtools": ("samtools", ("--version",)),
        "bcftools": ("bcftools", ("--version",)),
        "bgzip": ("bgzip", ("--version",)),
        "tabix": ("tabix", ("--version",)),
        "bedtools": ("bedtools", ("--version",)),
        "vcftools": ("vcftools", ("--version",)),
        "wfmash": ("wfmash", ("--version",)),
    }
    versions: Dict[str, str] = {}
    executable_paths: Dict[str, str] = {}
    for label, (name, arguments) in tools.items():
        executable = profile.resolve() / "bin" / name
        if not executable.is_file():
            raise WorkflowError(f"realized profile lacks {name}")
        resolved = executable.resolve()
        if not str(resolved).startswith("/gnu/store/"):
            raise WorkflowError(f"profile tool does not resolve into Guix store: {resolved}")
        fallback = (
            "wfmash commit e040aa10e87cab44ed5a4db005e784be62b0bd21 "
            "(upstream binary emits no version string)"
            if label == "wfmash"
            else None
        )
        versions[label] = _command_version(executable, arguments, empty_fallback=fallback)
        executable_paths[label] = str(resolved)
    record: Dict[str, Any] = {
        "schema_version": "1.0",
        "decision_version": DECISION_VERSION,
        "recorded_at": _utc_now(),
        "realization_host": socket.getfqdn(),
        "channel_commit": CHANNEL_COMMIT,
        "channels_file": str(CHANNELS.relative_to(ROOT)),
        "channels_sha256": sha256_file(CHANNELS),
        "resolved_channels_scm": resolved_channels_file.read_text(encoding="utf-8"),
        "resolved_channels_sha256": sha256_file(resolved_channels_file),
        "manifest_file": str(MANIFEST.relative_to(ROOT)),
        "manifest_sha256": sha256_file(MANIFEST),
        "profile_gc_root": str(profile.absolute()),
        "profile_store_path": profile_store,
        "derivations": derivations,
        "store_paths": store_paths,
        "tool_store_paths": executable_paths,
        "tool_versions": versions,
        "execution": "guix time-machine -C analysis/guix/channels.scm -- shell -m analysis/guix/manifest.scm --pure",
        "pack_fallback": {"required": False, "reason": "shared /gnu/store execution is the primary topology"},
    }
    audit_environment_record(record)
    _atomic_json(output, record)
    return record


def _run_checked(argv: Sequence[str], *, stdout_path: Optional[Path] = None) -> str:
    if stdout_path is None:
        completed = subprocess.run(argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False)
        output = completed.stdout
    else:
        with stdout_path.open("w", encoding="utf-8") as handle:
            completed = subprocess.run(argv, stdout=handle, stderr=subprocess.PIPE, text=True, check=False)
        output = ""
    if completed.returncode:
        raise WorkflowError(f"smoke command failed ({' '.join(argv)}): {completed.stderr.strip()}")
    return output


def compute_smoke(environment_record: Path, output: Path) -> Dict[str, Any]:
    record = json.loads(environment_record.read_text(encoding="utf-8"))
    audit_environment_record(record)
    imports: Dict[str, str] = {}
    try:
        from importlib.metadata import version
        import Bio
        import numpy
        import pandas
        import pyfaidx
        import pysam
        import scipy
    except ImportError as error:
        raise WorkflowError(f"pinned Python import smoke failed: {error}") from error
    imports.update(
        {
            "biopython": Bio.__version__,
            "jsonschema": version("jsonschema"),
            "numpy": numpy.__version__,
            "pandas": pandas.__version__,
            "pyfaidx": pyfaidx.__version__,
            "pysam": pysam.__version__,
            "scipy": scipy.__version__,
        }
    )
    with tempfile.TemporaryDirectory(prefix="tier3-compute-smoke-", dir=os.environ.get("TIER3_SCRATCH_ROOT")) as raw:
        temporary = Path(raw)
        rng = random.Random(13)
        sequence = "".join(rng.choice("ACGT") for _ in range(30000))
        query = sequence[:15000] + ("A" if sequence[15000] != "A" else "C") + sequence[15001:]
        target_fasta, query_fasta = temporary / "h1.fa", temporary / "h2.fa"
        target_fasta.write_text(">h1\n" + sequence + "\n", encoding="utf-8")
        query_fasta.write_text(">h2\n" + query + "\n", encoding="utf-8")
        _run_checked(["samtools", "faidx", str(target_fasta)])
        _run_checked(["samtools", "faidx", str(query_fasta)])
        aligned = temporary / "aligned.paf"
        _run_checked(
            [
                "wfmash",
                str(target_fasta),
                str(query_fasta),
                "-p",
                "90",
                "-l",
                "25000",
                "-4",
                "-t",
                "2",
            ],
            stdout_path=aligned,
        )
        aligned_text = aligned.read_text(encoding="utf-8")
        if "cg:Z:" not in aligned_text or re.search(r"cg:Z:\S*M", aligned_text):
            raise WorkflowError("WFMASH smoke did not emit an extended =/X/I/D CIGAR")
        paf_audit = _audit_paf(aligned, read_fasta(target_fasta), read_fasta(query_fasta))

        vcf = temporary / "truth.vcf"
        vcf.write_text(
            "##fileformat=VCFv4.2\n##contig=<ID=h1,length=30000>\n"
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=Genotype>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ttruth\n"
            f"h1\t51\t.\t{sequence[50]}\t{'A' if sequence[50] != 'A' else 'C'}\t.\tPASS\t.\tGT\t0|1\n",
            encoding="utf-8",
        )
        bcf = temporary / "truth.normalized.bcf"
        _run_checked(["bcftools", "norm", "-f", str(target_fasta), "-Ob", "-o", str(bcf), str(vcf)])
        _run_checked(["bcftools", "index", "--csi", str(bcf)])
        variant_audit = _audit_variant_reference(bcf, read_fasta(target_fasta))
        callable_bed = temporary / "truth.callable.bed"
        callable_bed.write_text("h1\t0\t100\n", encoding="utf-8")
        from analysis.tier3a_vgp_compute import compute_tier3a

        denominator_truth = compute_tier3a(
            dataset_id="compute.smoke.truth",
            reference_fasta=target_fasta,
            normalized_bcf=bcf,
            callable_bed=callable_bed,
            sample="truth",
            modality="deposited_exact_reference_variants_plus_mask",
            reference_accession="SMOKE_1.1",
        )
        if denominator_truth["callable_bases"] != 100 or denominator_truth["heterozygous_snvs"] != 1:
            raise WorkflowError("callable denominator smoke truth disagrees with 1/100 construction")
        if abs(denominator_truth["individual_snv_heterozygosity"] - 0.01) > 1e-12:
            raise WorkflowError("callable denominator smoke did not produce exact 0.01 heterozygosity")
        smoke = {
            "status": "passed",
            "decision_version": DECISION_VERSION,
            "run_at": _utc_now(),
            "compute_host": socket.getfqdn(),
            "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
            "login_profile_store_path": record["profile_store_path"],
            "compute_profile_store_path": str(Path(record["profile_store_path"]).resolve()),
            "store_path_identity_passed": str(Path(record["profile_store_path"]).resolve()) == record["profile_store_path"],
            "python_imports": imports,
            "fasta_index": {"path": str(target_fasta) + ".fai", "passed": Path(str(target_fasta) + ".fai").is_file()},
            "wfmash_extended_cigar": paf_audit,
            "normalized_bcf": variant_audit,
            "bcf_csi_passed": Path(str(bcf) + ".csi").is_file(),
            "callable_denominator_truth": {
                "heterozygous_snvs": 1,
                "callable_bases": 100,
                "individual_snv_heterozygosity": denominator_truth["individual_snv_heterozygosity"],
                "passed": True,
            },
        }
    if not smoke["store_path_identity_passed"]:
        raise WorkflowError("compute node did not reproduce login-node profile store path")
    _atomic_json(output, smoke)
    return smoke


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    validate = commands.add_parser("preflight", help="validate every selected dataset before submission")
    validate.add_argument("workflow", type=Path)
    validate.add_argument("--tier", choices=sorted(TIER_MODES))
    count = commands.add_parser("count", help="print array element count for one tier")
    count.add_argument("workflow", type=Path)
    count.add_argument("--tier", required=True, choices=sorted(TIER_MODES))
    run = commands.add_parser("run-array", help="run one zero-based array element atomically")
    run.add_argument("workflow", type=Path)
    run.add_argument("--tier", required=True, choices=sorted(TIER_MODES))
    run.add_argument("--index", required=True, type=int)
    environment = commands.add_parser("record-environment", help="record a realized, rooted Guix profile")
    environment.add_argument("--profile", required=True, type=Path)
    environment.add_argument("--derivations", required=True, type=Path)
    environment.add_argument("--store-paths", required=True, type=Path)
    environment.add_argument("--resolved-channels", required=True, type=Path)
    environment.add_argument("--output", required=True, type=Path)
    audit = commands.add_parser("audit-environment", help="check a recorded Guix closure")
    audit.add_argument("record", type=Path)
    smoke = commands.add_parser("compute-smoke", help="run the real compute-node pure-Guix truth smoke")
    smoke.add_argument("--environment-record", required=True, type=Path)
    smoke.add_argument("--output", required=True, type=Path)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _parser().parse_args(argv)
    if args.command == "preflight":
        workflow = load_workflow(args.workflow)
        selected = [item for item in workflow["datasets"] if args.tier is None or item["tier"] == args.tier]
        results = [validate_dataset_preflight(item) for item in selected]
        print(json.dumps({"status": "passed", "datasets": results}, indent=2, sort_keys=True))
    elif args.command == "count":
        workflow = load_workflow(args.workflow)
        print(sum(item["tier"] == args.tier for item in workflow["datasets"]))
    elif args.command == "run-array":
        print(json.dumps(run_array_task(args.workflow, tier=args.tier, array_index=args.index), indent=2, sort_keys=True))
    elif args.command == "record-environment":
        print(
            json.dumps(
                record_environment(
                    args.profile,
                    args.derivations,
                    args.store_paths,
                    args.resolved_channels,
                    args.output,
                ),
                indent=2,
                sort_keys=True,
            )
        )
    elif args.command == "audit-environment":
        record = json.loads(args.record.read_text(encoding="utf-8"))
        audit_environment_record(record)
        print(f"environment record passed: {args.record}")
    elif args.command == "compute-smoke":
        print(json.dumps(compute_smoke(args.environment_record, args.output), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (WorkflowError, Tier3ValidationError) as error:
        print(f"tier3 workflow rejected: {error}", file=sys.stderr)
        raise SystemExit(2)
