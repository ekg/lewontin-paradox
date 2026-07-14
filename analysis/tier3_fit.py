#!/usr/bin/env python3
"""Fail-closed Tier 3 synthesis, models, and dependency-free figure renderer.

The synthesis never uses a generic ``pi`` field.  Population diversity,
deposited-call individual heterozygosity, and alignment-conditioned individual
heterozygosity remain separate observable tiers through loading, fitting, and
reporting.  Annotation-derived composition is admitted only when the upstream
validated row says ``native``.  Exact scientific-name joins prohibit congener
substitution.

The pinned Buffalo input supplies census-size covariates and taxonomic labels.
It is identified by an immutable commit and SHA-256.  The committed result TSV
contains the small joined point set, model claims, negative results, and
structured missingness; ``--figure-from-results`` regenerates both figures
without raw inputs or a non-declared plotting library.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import struct
import sys
import urllib.request
import zlib
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, MutableMapping, Optional, Sequence, Tuple

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = Path(__file__).resolve().parent
BUFFALO_COMMIT = "b8f91d5c34675733db8cae8dcab625dcbb55c30a"
BUFFALO_URL = (
    "https://raw.githubusercontent.com/vsbuffalo/paradox_variation/"
    f"{BUFFALO_COMMIT}/data/combined_data.tsv"
)
BUFFALO_SHA256 = "df559451dad94b53ba8675e09811708107a57eeb6ffe8f72b944bcbbf3a1f2eb"
DECISION_VERSION = "tier3-decisions-v1"
PRIMARY_FAMILY = "primary_composition"
DEFAULT_BOOTSTRAPS = 10_000
MIN_CLADE_N = 8


class SynthesisError(ValueError):
    """An input violates a frozen synthesis invariant."""


RESULT_FIELDS = [
    "row_kind",
    "analysis_id",
    "status",
    "observable",
    "observable_tier",
    "scientific_name",
    "dataset_id",
    "assembly_accession_version",
    "clade",
    "order",
    "family",
    "genus",
    "scope",
    "model",
    "predictor",
    "predictor_value",
    "n",
    "eligible_n",
    "missing_n",
    "denominator",
    "estimate",
    "effect",
    "ci_low",
    "ci_high",
    "p_value",
    "q_value",
    "analysis_family",
    "uncertainty",
    "measurement_precision",
    "cross_species_identification",
    "sensitivity",
    "conclusion",
    "limitation",
    "annotation_status",
    "annotation_provider",
    "annotation_release",
    "annotation_assembly_accession_version",
    "reference_fasta_sha256",
    "annotation_gff_sha256",
    "annotation_contig_mapping_sha256",
    "annotation_sequence_regions_sha256",
    "annotation_genetic_code",
    "annotation_native_vs_projected",
    "annotation_contig_dictionary_validated",
    "annotation_cds_reconstruction_audit",
    "source_modality",
    "decision_version",
]


def empty_result_row() -> Dict[str, str]:
    row = {field: "" for field in RESULT_FIELDS}
    row["decision_version"] = DECISION_VERSION
    return row


def read_tsv(path: Path | str) -> List[Dict[str, str]]:
    path = Path(path)
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _finite(value: Any) -> Optional[float]:
    if value in (None, "", "NA", "NaN", "nan"):
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _unique(rows: Sequence[Mapping[str, Any]], key: str, label: str) -> Dict[str, Mapping[str, Any]]:
    index: Dict[str, Mapping[str, Any]] = {}
    for row in rows:
        value = str(row.get(key, "")).strip()
        if not value:
            raise SynthesisError(f"{label} has an empty {key}")
        if value in index:
            raise SynthesisError(f"duplicate {label} {key}: {value!r}")
        index[value] = row
    return index


def read_buffalo_rows(rows: Sequence[Mapping[str, str]]) -> List[Dict[str, str]]:
    """Validate and retain Buffalo core rows without hiding duplicates."""

    _unique(rows, "species", "Buffalo")
    core = []
    for source in rows:
        if _finite(source.get("diversity")) is None or _finite(source.get("pred_log10_N")) is None:
            continue
        row = dict(source)
        if not row.get("class") or row.get("class") == "NA":
            row["class"] = "unclassified"
        for key in ("order", "family", "genus", "phylum"):
            if not row.get(key) or row.get(key) == "NA":
                row[key] = "unavailable"
        core.append(row)
    return core


def read_buffalo(path: Path | str) -> List[Dict[str, str]]:
    return read_buffalo_rows(read_tsv(path))


def fetch_buffalo(destination: Path) -> Path:
    """Fetch the immutable Buffalo file and verify it before publication."""

    data = urllib.request.urlopen(BUFFALO_URL, timeout=120).read()
    observed = hashlib.sha256(data).hexdigest()
    if observed != BUFFALO_SHA256:
        raise SynthesisError(f"Buffalo SHA-256 mismatch: {observed}")
    destination.write_bytes(data)
    return destination


def load_composition_provenance(directory: Path | str) -> Dict[str, Dict[str, Any]]:
    """Load the small, committed per-dataset Tier 3c QC provenance."""

    records: Dict[str, Dict[str, Any]] = {}
    for path in sorted(Path(directory).glob("*.json")):
        record = json.loads(path.read_text(encoding="utf-8"))
        dataset_id = record.get("dataset_id")
        if not dataset_id:
            raise SynthesisError(f"Tier 3c QC record lacks dataset_id: {path}")
        if dataset_id in records:
            raise SynthesisError(f"duplicate Tier 3c QC dataset_id: {dataset_id!r}")
        records[str(dataset_id)] = record
    return records


def join_composition(
    composition_rows: Sequence[Mapping[str, str]],
    buffalo_rows: Sequence[Mapping[str, str]],
    provenance: Optional[Mapping[str, Mapping[str, Any]]] = None,
) -> List[Dict[str, Any]]:
    """Perform a one-to-one exact-species join and enforce provenance gates."""

    composition = _unique(composition_rows, "scientific_name", "Tier 3c")
    buffalo = _unique(buffalo_rows, "species", "Buffalo core")
    missing = sorted(set(composition) - set(buffalo))
    if missing:
        preview = ", ".join(missing[:3])
        raise SynthesisError(
            f"exact-species Buffalo join failed for {len(missing)} row(s): {preview}; "
            "congener substitution is prohibited"
        )
    joined: List[Dict[str, Any]] = []
    for name in sorted(composition):
        source = composition[name]
        metadata = buffalo[name]
        dataset_id = str(source.get("dataset_id", ""))
        provenance_fields: Dict[str, Any] = {}
        if provenance is not None:
            if dataset_id not in provenance:
                raise SynthesisError(f"missing Tier 3c QC provenance for {dataset_id!r}")
            qc = provenance[dataset_id]
            annotation = qc.get("annotation_provenance") or {}
            reference = qc.get("reference") or {}
            if reference.get("accession") != source.get("assembly_accession_version"):
                raise SynthesisError(f"QC reference accession mismatch for {dataset_id!r}")
            required_reference = ("fasta_sha256", "contig_dictionary_sha256")
            if any(not reference.get(field) for field in required_reference):
                raise SynthesisError(f"incomplete reference QC provenance for {dataset_id!r}")
            if source.get("annotation_status") == "native":
                required_annotation = (
                    "provider",
                    "release",
                    "assembly_accession",
                    "fasta_sha256",
                    "gff_sha256",
                    "contig_mapping_sha256",
                    "sequence_regions_sha256",
                    "genetic_code",
                    "native_vs_projected",
                    "contig_dictionary_validated",
                    "cds_audit",
                )
                if any(annotation.get(field) in (None, "", False) for field in required_annotation):
                    raise SynthesisError(f"incomplete native-annotation QC provenance for {dataset_id!r}")
                if annotation.get("assembly_accession") != reference.get("accession"):
                    raise SynthesisError(f"QC FASTA/GFF accession mismatch for {dataset_id!r}")
                if annotation.get("fasta_sha256") != reference.get("fasta_sha256"):
                    raise SynthesisError(f"QC FASTA/annotation hash mismatch for {dataset_id!r}")
                if annotation.get("sequence_regions_sha256") != reference.get("contig_dictionary_sha256"):
                    raise SynthesisError(f"QC sequence-region/FASTA dictionary mismatch for {dataset_id!r}")
                if annotation.get("native_vs_projected") != "native":
                    raise SynthesisError(f"non-native annotation promoted for {dataset_id!r}")
                cds_audit = annotation["cds_audit"]
                if not cds_audit.get("all_retained_cds_phase_translation_passed") or cds_audit.get(
                    "sampled_cds_mismatches"
                ) != 0:
                    raise SynthesisError(f"failed CDS reconstruction audit for {dataset_id!r}")
            provenance_fields = {
                "annotation_provider": annotation.get("provider", ""),
                "annotation_release": annotation.get("release", ""),
                "annotation_assembly_accession_version": annotation.get("assembly_accession", ""),
                "reference_fasta_sha256": reference.get("fasta_sha256", ""),
                "annotation_gff_sha256": annotation.get("gff_sha256", ""),
                "annotation_contig_mapping_sha256": annotation.get("contig_mapping_sha256", ""),
                "annotation_sequence_regions_sha256": annotation.get("sequence_regions_sha256", ""),
                "annotation_genetic_code": annotation.get("genetic_code", ""),
                "annotation_native_vs_projected": annotation.get("native_vs_projected", ""),
                "annotation_contig_dictionary_validated": annotation.get(
                    "contig_dictionary_validated", ""
                ),
                "annotation_cds_reconstruction_audit": json.dumps(
                    annotation.get("cds_audit", {}), sort_keys=True, separators=(",", ":")
                ),
            }
        embedded = _finite(source.get("buffalo_pred_log10_N"))
        pinned = _finite(metadata.get("pred_log10_N"))
        if embedded is None or pinned is None or not math.isclose(embedded, pinned, abs_tol=1e-10):
            raise SynthesisError(f"Buffalo predictor mismatch for exact species {name!r}")
        annotation_status = str(source.get("annotation_status", ""))
        raw_gc3 = _finite(source.get("gc3"))
        gc3 = raw_gc3 if annotation_status == "native" else None
        if raw_gc3 is not None and annotation_status != "native":
            gc3_reason = "annotation_not_native"
        elif raw_gc3 is None:
            gc3_reason = "native_annotation_unavailable"
        else:
            gc3_reason = ""
        whole_gc = _finite(source.get("whole_genome_gc"))
        if whole_gc is None:
            raise SynthesisError(f"Tier 3c exact FASTA whole-genome GC is missing for {name!r}")
        joined.append(
            {
                "dataset_id": dataset_id,
                "scientific_name": name,
                "taxon_id": source.get("taxon_id", ""),
                "assembly_accession_version": source.get("assembly_accession_version", ""),
                "pred_log10_N": pinned,
                "buffalo_diversity": _finite(source.get("buffalo_diversity")),
                "phylum": metadata.get("phylum", "unavailable"),
                "class": metadata.get("class", "unclassified"),
                "order": metadata.get("order", "unavailable"),
                "family": metadata.get("family", "unavailable"),
                "genus": metadata.get("genus", "unavailable"),
                "ave_rec": _finite(metadata.get("ave_rec")),
                "map_length": _finite(metadata.get("map_length")),
                "gc3": gc3,
                "gc3_denominator": int(source["gc3_callable_third_positions"])
                if gc3 is not None and source.get("gc3_callable_third_positions")
                else 0,
                "gc3_missing_reason": gc3_reason,
                "whole_genome_gc": whole_gc,
                "whole_genome_denominator": int(source.get("whole_genome_callable_bases", "0")),
                "annotation_status": annotation_status,
                **provenance_fields,
            }
        )
    return joined


@dataclass(frozen=True)
class ModelFit:
    n: int
    effect: float
    ci_low: float
    ci_high: float
    p_value: float
    intercept: float
    r_squared: float


def _design(
    rows: Sequence[Mapping[str, Any]], outcome: str, fixed_effect: Optional[str]
) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    clean = [
        row
        for row in rows
        if _finite(row.get("pred_log10_N")) is not None and _finite(row.get(outcome)) is not None
    ]
    if len(clean) < 3:
        raise SynthesisError(f"fewer than three complete rows for {outcome}")
    x = np.asarray([float(row["pred_log10_N"]) for row in clean], dtype=float)
    columns = [np.ones(len(clean)), x]
    if fixed_effect:
        levels = sorted({str(row.get(fixed_effect, "unclassified")) for row in clean})
        for level in levels[1:]:
            columns.append(np.asarray([row.get(fixed_effect) == level for row in clean], dtype=float))
    return np.column_stack(columns), np.asarray([float(row[outcome]) for row in clean]), [
        str(row.get("scientific_name", index)) for index, row in enumerate(clean)
    ]


def _ols_effect(rows: Sequence[Mapping[str, Any]], outcome: str, fixed_effect: Optional[str]) -> Tuple[float, ...]:
    design, y, _ = _design(rows, outcome, fixed_effect)
    if len(y) <= design.shape[1]:
        raise SynthesisError(f"insufficient residual degrees of freedom for {outcome}")
    coef, _residuals, rank, _singular = np.linalg.lstsq(design, y, rcond=None)
    if rank < design.shape[1]:
        raise SynthesisError(f"rank-deficient model for {outcome}")
    fitted = design @ coef
    residual = y - fitted
    rss = float(residual @ residual)
    centered = y - y.mean()
    tss = float(centered @ centered)
    r2 = 1.0 - rss / tss if tss > 0 else 1.0
    df = len(y) - design.shape[1]
    covariance = (rss / df) * np.linalg.inv(design.T @ design)
    se = math.sqrt(max(0.0, float(covariance[1, 1])))
    if se == 0:
        p_value = 0.0 if coef[1] != 0 else 1.0
    else:
        # Normal approximation is reported explicitly; bootstrap percentiles
        # carry the effect interval and are the primary uncertainty statement.
        p_value = math.erfc(abs(float(coef[1]) / se) / math.sqrt(2.0))
    return float(coef[0]), float(coef[1]), r2, p_value


def _bootstrap_seed(outcome: str, fixed_effect: Optional[str]) -> int:
    digest = hashlib.sha256(
        f"{DECISION_VERSION}\0cross-species-fit\0{outcome}\0{fixed_effect or 'none'}".encode()
    ).digest()
    return int.from_bytes(digest[:8], "big")


def fit_ols(
    rows: Sequence[Mapping[str, Any]],
    outcome: str,
    *,
    fixed_effect: Optional[str] = None,
    bootstrap_replicates: int = DEFAULT_BOOTSTRAPS,
) -> ModelFit:
    clean = [
        row
        for row in rows
        if _finite(row.get("pred_log10_N")) is not None and _finite(row.get(outcome)) is not None
    ]
    intercept, effect, r2, p_value = _ols_effect(clean, outcome, fixed_effect)
    rng = np.random.Generator(np.random.PCG64(_bootstrap_seed(outcome, fixed_effect)))
    effects: List[float] = []
    # Stratification keeps every fixed-effect level represented and separates
    # cross-species sampling uncertainty from per-assembly measurement counts.
    groups: List[List[Mapping[str, Any]]]
    if fixed_effect:
        groups = [
            [row for row in clean if str(row.get(fixed_effect)) == level]
            for level in sorted({str(row.get(fixed_effect)) for row in clean})
        ]
    else:
        groups = [clean]
    for _ in range(max(0, bootstrap_replicates)):
        sample: List[Mapping[str, Any]] = []
        for group in groups:
            sample.extend(group[int(index)] for index in rng.integers(0, len(group), len(group)))
        try:
            effects.append(_ols_effect(sample, outcome, fixed_effect)[1])
        except SynthesisError:
            continue
    if effects:
        low, high = np.percentile(np.asarray(effects), [2.5, 97.5])
    else:
        low = high = effect
    return ModelFit(len(clean), effect, float(low), float(high), p_value, intercept, r2)


def _quadratic_effect(
    rows: Sequence[Mapping[str, Any]], outcome: str, fixed_effect: Optional[str]
) -> Tuple[float, ...]:
    clean = [
        row
        for row in rows
        if _finite(row.get("pred_log10_N")) is not None and _finite(row.get(outcome)) is not None
    ]
    if len(clean) < 5:
        raise SynthesisError(f"fewer than five complete rows for quadratic {outcome}")
    x = np.asarray([float(row["pred_log10_N"]) for row in clean])
    centered = x - x.mean()
    columns = [np.ones(len(clean)), centered, centered**2]
    if fixed_effect:
        levels = sorted({str(row.get(fixed_effect, "unclassified")) for row in clean})
        for level in levels[1:]:
            columns.append(np.asarray([row.get(fixed_effect) == level for row in clean], dtype=float))
    design = np.column_stack(columns)
    y = np.asarray([float(row[outcome]) for row in clean])
    if len(y) <= design.shape[1]:
        raise SynthesisError(f"insufficient residual degrees of freedom for quadratic {outcome}")
    coef, _residuals, rank, _singular = np.linalg.lstsq(design, y, rcond=None)
    if rank < design.shape[1]:
        raise SynthesisError(f"rank-deficient quadratic model for {outcome}")
    residual = y - design @ coef
    rss = float(residual @ residual)
    tss = float((y - y.mean()) @ (y - y.mean()))
    df = len(y) - design.shape[1]
    covariance = (rss / df) * np.linalg.inv(design.T @ design)
    se = math.sqrt(max(0.0, float(covariance[2, 2])))
    p_value = math.erfc(abs(float(coef[2]) / se) / math.sqrt(2.0)) if se else (0.0 if coef[2] else 1.0)
    return float(coef[0]), float(coef[2]), (1.0 - rss / tss if tss else 1.0), p_value


def fit_quadratic(
    rows: Sequence[Mapping[str, Any]],
    outcome: str,
    *,
    fixed_effect: Optional[str] = None,
    bootstrap_replicates: int = DEFAULT_BOOTSTRAPS,
) -> ModelFit:
    """Fit the predeclared concavity sensitivity on centered log census size."""

    clean = [
        row
        for row in rows
        if _finite(row.get("pred_log10_N")) is not None and _finite(row.get(outcome)) is not None
    ]
    intercept, effect, r2, p_value = _quadratic_effect(clean, outcome, fixed_effect)
    rng = np.random.Generator(np.random.PCG64(_bootstrap_seed(outcome + ".quadratic", fixed_effect)))
    if fixed_effect:
        groups = [
            [row for row in clean if str(row.get(fixed_effect)) == level]
            for level in sorted({str(row.get(fixed_effect)) for row in clean})
        ]
    else:
        groups = [clean]
    effects = []
    for _ in range(max(0, bootstrap_replicates)):
        sample = []
        for group in groups:
            sample.extend(group[int(index)] for index in rng.integers(0, len(group), len(group)))
        try:
            effects.append(_quadratic_effect(sample, outcome, fixed_effect)[1])
        except SynthesisError:
            continue
    low, high = np.percentile(effects, [2.5, 97.5]) if effects else (effect, effect)
    return ModelFit(len(clean), effect, float(low), float(high), p_value, intercept, r2)


def split_clades(
    rows: Sequence[Mapping[str, Any]], *, min_n: int = MIN_CLADE_N
) -> Dict[str, List[Mapping[str, Any]]]:
    levels = sorted({str(row.get("class", "unclassified")) for row in rows})
    return {
        level: [row for row in rows if str(row.get("class", "unclassified")) == level]
        for level in levels
        if sum(
            _finite(row.get("gc3")) is not None
            for row in rows
            if str(row.get("class", "unclassified")) == level
        )
        >= min_n
    }


def leave_one_out(
    rows: Sequence[Mapping[str, Any]], outcome: str, *, unit: str
) -> Dict[str, float | int]:
    if unit == "species":
        levels = [str(row.get("scientific_name", index)) for index, row in enumerate(rows)]
        subsets = [rows[:index] + rows[index + 1 :] for index in range(len(rows))]
    else:
        levels = sorted({str(row.get(unit, "unclassified")) for row in rows})
        subsets = [[row for row in rows if str(row.get(unit, "unclassified")) != level] for level in levels]
    effects = []
    for subset in subsets:
        try:
            effects.append(fit_ols(subset, outcome, bootstrap_replicates=0).effect)
        except SynthesisError:
            continue
    if not effects:
        return {"n": 0, "min": math.nan, "max": math.nan}
    return {"n": len(effects), "min": min(effects), "max": max(effects)}


def _observation(
    row: Mapping[str, str], observable: str, tier: str, value: float, denominator: str, modality: str
) -> Dict[str, Any]:
    return {
        "dataset_id": row.get("dataset_id", ""),
        "scientific_name": row.get("scientific_name", ""),
        "observable": observable,
        "observable_tier": tier,
        "value": value,
        "denominator": denominator,
        "source_modality": modality,
    }


def load_diversity_observations(
    tier3a_path: Optional[Path | str], tier3b_path: Optional[Path | str]
) -> List[Dict[str, Any]]:
    observations: List[Dict[str, Any]] = []
    if tier3a_path:
        for row in read_tsv(tier3a_path):
            modality = row.get("modality", "")
            if "deposited" in modality:
                tier = "deposited_vgp_individual"
            elif "wfmash" in modality or "alignment" in modality:
                tier = "alignment_conditioned_individual"
            else:
                raise SynthesisError(f"unrecognized Tier 3a modality {modality!r}")
            value = _finite(row.get("individual_snv_heterozygosity"))
            if value is not None:
                observations.append(
                    _observation(
                        row,
                        "individual_snv_heterozygosity",
                        tier,
                        value,
                        f"{row.get('total_denominator', '')} callable/alignable bases",
                        modality,
                    )
                )
            ratio = _finite(row.get("pi_S_over_pi_W"))
            if ratio is not None:
                observations.append(
                    _observation(
                        row,
                        "pi_S_over_pi_W",
                        tier,
                        ratio,
                        f"W={row.get('fourfold_W_denominator', '')};S={row.get('fourfold_S_denominator', '')}",
                        modality,
                    )
                )
    if tier3b_path:
        for row in read_tsv(tier3b_path):
            value = _finite(row.get("population_pi"))
            if value is not None:
                observations.append(
                    _observation(
                        row,
                        "population_pi",
                        "population",
                        value,
                        f"{row.get('population_pi_denominator', '')} callable cohort sites",
                        "deposited_population_variants",
                    )
                )
            ratio = _finite(row.get("pi_S_over_pi_W"))
            if ratio is not None:
                observations.append(
                    _observation(
                        row,
                        "pi_S_over_pi_W",
                        "population",
                        ratio,
                        f"W={row.get('pi_W_denominator', '')};S={row.get('pi_S_denominator', '')}",
                        "deposited_population_variants",
                    )
                )
    return observations


def _claim(
    analysis_id: str,
    status: str,
    observable: str,
    tier: str,
    n: int,
    limitation: str,
    *,
    conclusion: str,
) -> Dict[str, str]:
    row = empty_result_row()
    row.update(
        row_kind="claim",
        analysis_id=analysis_id,
        status=status,
        observable=observable,
        observable_tier=tier,
        scope="observable_status",
        model="not_fit" if status != "deferred" else "deferred_by_frozen_plan",
        n=str(n),
        eligible_n=str(n),
        missing_n="not_applicable" if status == "deferred" else "all_predeclared_candidates",
        denominator="no eligible checksum-locked numerator/denominator tuple" if n == 0 else "see point rows",
        effect="not_estimable",
        ci_low="not_estimable",
        ci_high="not_estimable",
        uncertainty="unavailable" if n == 0 else "upstream frozen block-bootstrap intervals",
        measurement_precision="not estimable" if n == 0 else "conditional on frozen callable denominator",
        cross_species_identification="not identified" if n < 3 else "observable-specific fit only",
        conclusion=conclusion,
        limitation=limitation,
    )
    return row


def diversity_claim_rows(observations: Sequence[Mapping[str, Any]]) -> List[Dict[str, str]]:
    specifications = [
        ("population_pi", "population"),
        ("individual_snv_heterozygosity", "deposited_vgp_individual"),
        ("individual_snv_heterozygosity", "alignment_conditioned_individual"),
        ("pi_S_over_pi_W", "population"),
        ("pi_S_over_pi_W", "deposited_vgp_individual"),
        ("pi_S_over_pi_W", "alignment_conditioned_individual"),
    ]
    rows = []
    for observable, tier in specifications:
        selected = [
            item
            for item in observations
            if item["observable"] == observable and item["observable_tier"] == tier
        ]
        n = len(selected)
        rows.append(
            _claim(
                f"status.{observable}.{tier}",
                "estimated" if n else "unavailable",
                observable,
                tier,
                n,
                "No eligible upstream row passed exact-reference, native-annotation, sample, and invariant-denominator gates. Absence is not zero."
                if not n
                else "Distinct observable tier; it is not interchangeable with any other diversity modality.",
                conclusion="No Tier 3 estimate or cross-species effect is available."
                if not n
                else "Eligible estimates exist only for this named modality.",
            )
        )
    rows.append(
        _claim(
            "status.polarized_sfs_B",
            "deferred",
            "polarized_sfs_B",
            "population",
            0,
            "Frozen v1 supplies no outgroup, ancestral-error model, demography, or power threshold; pi_S_over_pi_W is not B.",
            conclusion="No SFS-B estimate or causal strength claim was made.",
        )
    )
    return rows


def _point_rows(joined: Sequence[Mapping[str, Any]]) -> List[Dict[str, str]]:
    output: List[Dict[str, str]] = []
    for source in joined:
        common = {
            "scientific_name": str(source["scientific_name"]),
            "dataset_id": str(source["dataset_id"]),
            "assembly_accession_version": str(source["assembly_accession_version"]),
            "clade": str(source["class"]),
            "order": str(source["order"]),
            "family": str(source["family"]),
            "genus": str(source["genus"]),
            "predictor": "Buffalo pred_log10_N (census-size proxy, not N_e)",
            "predictor_value": _format(source["pred_log10_N"]),
            "n": "1",
            "eligible_n": "1",
            "missing_n": "0",
            "annotation_status": str(source["annotation_status"]),
            "source_modality": "exact_single_assembly_composition",
            "annotation_provider": str(source.get("annotation_provider", "")),
            "annotation_release": str(source.get("annotation_release", "")),
            "annotation_assembly_accession_version": str(
                source.get("annotation_assembly_accession_version", "")
            ),
            "reference_fasta_sha256": str(source.get("reference_fasta_sha256", "")),
            "annotation_gff_sha256": str(source.get("annotation_gff_sha256", "")),
            "annotation_contig_mapping_sha256": str(
                source.get("annotation_contig_mapping_sha256", "")
            ),
            "annotation_sequence_regions_sha256": str(
                source.get("annotation_sequence_regions_sha256", "")
            ),
            "annotation_genetic_code": str(source.get("annotation_genetic_code", "")),
            "annotation_native_vs_projected": str(
                source.get("annotation_native_vs_projected", "")
            ),
            "annotation_contig_dictionary_validated": str(
                source.get("annotation_contig_dictionary_validated", "")
            ),
            "annotation_cds_reconstruction_audit": str(
                source.get("annotation_cds_reconstruction_audit", "")
            ),
        }
        gc3 = empty_result_row()
        gc3.update(common)
        gc3.update(
            row_kind="point",
            analysis_id=f"point.gc3.{source['dataset_id']}",
            status="estimated" if source["gc3"] is not None else "unavailable",
            observable="gc3",
            observable_tier="exact_assembly_native_annotation_composition",
            scope="species",
            model="pooled_canonical_CDS_third_positions",
            denominator=f"{source['gc3_denominator']} callable CDS third positions",
            estimate=_format(source["gc3"]),
            uncertainty="Upstream table has no genomic block interval; no binomial pseudo-replication interval is invented.",
            measurement_precision="Callable third-position count reported; bases are correlated and do not identify cross-species uncertainty.",
            cross_species_identification="not applicable to point row",
            conclusion="Exact-assembly native-annotation GC3 point estimate."
            if source["gc3"] is not None
            else "GC3 unavailable; whole-genome GC remains a separate observable.",
            limitation=str(source["gc3_missing_reason"]),
        )
        output.append(gc3)
        whole = empty_result_row()
        whole.update(common)
        whole.update(
            row_kind="point",
            analysis_id=f"point.whole_genome_gc.{source['dataset_id']}",
            status="estimated",
            observable="whole_genome_gc",
            observable_tier="exact_assembly_composition",
            scope="species",
            model="ACGT_fraction",
            denominator=f"{source['whole_genome_denominator']} callable A/C/G/T FASTA bases",
            estimate=_format(source["whole_genome_gc"]),
            uncertainty="No genomic block interval is available in the validated small table.",
            measurement_precision="Callable FASTA-base count reported; not used as cross-species regression weight.",
            cross_species_identification="not applicable to point row",
            conclusion="Whole-genome GC is available independently of annotation.",
            limitation="Weaker composition proxy; never substituted for GC3.",
        )
        output.append(whole)
    return output


def _format(value: Any) -> str:
    if value is None or value == "":
        return ""
    if isinstance(value, float) and math.isnan(value):
        return ""
    if isinstance(value, (float, np.floating)):
        return f"{float(value):.12g}"
    return str(value)


def _model_row(
    analysis_id: str,
    outcome: str,
    tier: str,
    scope: str,
    model: str,
    fit: ModelFit,
    total: int,
    denominator: str,
    *,
    clade: str = "",
    family: str = PRIMARY_FAMILY,
    sensitivity: str = "",
    limitation: str,
) -> Dict[str, str]:
    row = empty_result_row()
    row.update(
        row_kind="model",
        analysis_id=analysis_id,
        status="estimated",
        observable=outcome,
        observable_tier=tier,
        clade=clade,
        scope=scope,
        model=model,
        predictor="Buffalo pred_log10_N (census-size proxy, not N_e)",
        n=str(fit.n),
        eligible_n=str(fit.n),
        missing_n=str(total - fit.n),
        denominator=denominator,
        effect=_format(fit.effect),
        ci_low=_format(fit.ci_low),
        ci_high=_format(fit.ci_high),
        p_value=_format(fit.p_value),
        analysis_family=family,
        uncertainty=f"10,000 species-bootstrap percentile interval; p uses a normal approximation; R2={fit.r_squared:.4g}",
        measurement_precision="Per-assembly denominators are reported in point rows and are not treated as independent observations or regression weights.",
        cross_species_identification="Species bootstrap conditional on Buffalo predictor and taxonomic labels; it is not measurement or causal uncertainty.",
        sensitivity=sensitivity,
        conclusion="Effect retained regardless of sign or multiplicity-adjusted significance.",
        limitation=limitation,
    )
    return row


def apply_bh(rows: Sequence[MutableMapping[str, str]]) -> None:
    """Benjamini-Hochberg within each named, predeclared family."""

    families = sorted({row.get("analysis_family", "") for row in rows if row.get("p_value")})
    for family in families:
        indexed = [
            (index, float(row["p_value"]))
            for index, row in enumerate(rows)
            if row.get("analysis_family") == family and row.get("p_value")
        ]
        ordered = sorted(indexed, key=lambda item: item[1])
        adjusted: Dict[int, float] = {}
        running = 1.0
        m = len(ordered)
        for rank_index in range(m - 1, -1, -1):
            original_index, p_value = ordered[rank_index]
            rank = rank_index + 1
            running = min(running, p_value * m / rank)
            adjusted[original_index] = min(1.0, running)
        for index, value in adjusted.items():
            rows[index]["q_value"] = _format(value)


def synthesize(
    joined: Sequence[Mapping[str, Any]], observations: Sequence[Mapping[str, Any]]
) -> List[Dict[str, str]]:
    results = _point_rows(joined)
    gc3_rows = [row for row in joined if row.get("gc3") is not None]
    whole_rows = list(joined)
    across = fit_ols(gc3_rows, "gc3", fixed_effect="class")
    species_loo = leave_one_out(gc3_rows, "gc3", unit="species")
    clade_loo = leave_one_out(gc3_rows, "gc3", unit="class")
    sensitivity = (
        f"leave-one-species-out slopes [{_format(species_loo['min'])}, {_format(species_loo['max'])}] "
        f"over {species_loo['n']} fits; leave-one-class-out slopes "
        f"[{_format(clade_loo['min'])}, {_format(clade_loo['max'])}] over {clade_loo['n']} fits"
    )
    results.append(
        _model_row(
            "primary.gc3.across_class_fixed",
            "gc3",
            "exact_assembly_native_annotation_composition",
            "across_clades",
            "OLS with taxonomic-class fixed intercepts",
            across,
            len(joined),
            f"{len(gc3_rows)}/{len(joined)} exact-species rows have native exact-assembly annotation",
            sensitivity=sensitivity,
            limitation="Taxonomic-class fixed effects are a coarse phylogenetic control, not PGLS; no frozen branch-length tree was supplied. Predictor is census size, not N_e.",
        )
    )
    for clade, subset in split_clades(gc3_rows).items():
        model = fit_ols(subset, "gc3")
        loo = leave_one_out(subset, "gc3", unit="species")
        results.append(
            _model_row(
                f"primary.gc3.within_class.{clade.lower().replace(' ', '_')}",
                "gc3",
                "exact_assembly_native_annotation_composition",
                "within_clade",
                "within-class OLS",
                model,
                len([row for row in joined if row["class"] == clade]),
                f"native exact-assembly GC3 rows in taxonomic class {clade}",
                clade=clade,
                sensitivity=f"leave-one-species-out slopes [{_format(loo['min'])}, {_format(loo['max'])}] over {loo['n']} fits",
                limitation="Within-class taxonomy does not guarantee conserved recombination machinery; observational census-size association only.",
            )
        )
    drosophila = [row for row in gc3_rows if row.get("genus") == "Drosophila"]
    if len(drosophila) >= MIN_CLADE_N:
        model = fit_ols(drosophila, "gc3")
        loo = leave_one_out(drosophila, "gc3", unit="species")
        results.append(
            _model_row(
                "primary.gc3.within_genus.drosophila",
                "gc3",
                "exact_assembly_native_annotation_composition",
                "within_genus",
                "within-genus OLS",
                model,
                len([row for row in joined if row.get("genus") == "Drosophila"]),
                "native exact-assembly Drosophila GC3 rows",
                clade="Drosophila",
                sensitivity=f"leave-one-species-out slopes [{_format(loo['min'])}, {_format(loo['max'])}] over {loo['n']} fits",
                limitation="Genus is only a proxy for conserved machinery; no predeclared recombination-class labels exist.",
            )
        )
    whole_fit = fit_ols(whole_rows, "whole_genome_gc", fixed_effect="class")
    results.append(
        _model_row(
            "sensitivity.whole_genome_gc.across_class_fixed",
            "whole_genome_gc",
            "exact_assembly_composition",
            "across_clades",
            "OLS with taxonomic-class fixed intercepts",
            whole_fit,
            len(joined),
            f"{len(whole_rows)}/{len(joined)} exact FASTA rows",
            family="sensitivity_composition",
            sensitivity="alternate composition observable; not a replacement for GC3",
            limitation="Whole-genome GC is affected by noncoding composition and other processes; causal gBGC inference is unsupported.",
        )
    )
    quadratic = fit_quadratic(gc3_rows, "gc3", fixed_effect="class")
    quadratic_row = _model_row(
        "sensitivity.gc3.quadratic_concavity",
        "gc3",
        "exact_assembly_native_annotation_composition",
        "across_clades",
        "centered quadratic OLS with taxonomic-class fixed intercepts",
        quadratic,
        len(joined),
        f"{len(gc3_rows)}/{len(joined)} native exact-assembly GC3 rows",
        family="sensitivity_composition",
        sensitivity="negative quadratic coefficient is the predeclared concavity/saturation-shape direction",
        limitation="A quadratic is a shape diagnostic, not a mechanistic saturation curve; census size is not N_e and taxonomic fixed effects are not PGLS.",
    )
    quadratic_row["predictor"] = "centered Buffalo pred_log10_N squared"
    quadratic_row["conclusion"] = (
        "Concavity is supported only if the interval excludes zero below zero; otherwise GC3 saturation is not established."
    )
    results.append(quadratic_row)
    results.extend(diversity_claim_rows(observations))
    for analysis_id, limitation in [
        (
            "status.phylogenetic_pgls",
            "No frozen species tree with branch lengths or covariance model was supplied; class-fixed effects are reported separately and are not called PGLS.",
        ),
        (
            "status.recombination_class_model",
            "Buffalo provides sparse recombination measurements but no predeclared recombination-machinery classes; no post-hoc classes were invented.",
        ),
        (
            "sensitivity.alternate_callable_downsampling",
            "No eligible population or individual diversity tuple exists, so alternate callable-mask and deterministic downsampling sensitivities are unavailable rather than contradictory values being hidden.",
        ),
        (
            "claim.causal_gbgc",
            "Composition associations alone cannot distinguish gBGC from mutation, selection, life history, or assembly/annotation structure; diversity and SFS-B tests are unavailable/deferred.",
        ),
    ]:
        row = _claim(
            analysis_id,
            "unavailable" if analysis_id != "claim.causal_gbgc" else "not_identified",
            "gc3" if "callable" not in analysis_id else "diversity_sensitivity",
            "cross_species_identification",
            len(gc3_rows) if analysis_id in {"status.phylogenetic_pgls", "status.recombination_class_model", "claim.causal_gbgc"} else 0,
            limitation,
            conclusion="No causal gBGC effect is claimed."
            if analysis_id == "claim.causal_gbgc"
            else "Analysis was not performed; missing metadata or eligible observables are explicit.",
        )
        row["denominator"] = f"{len(gc3_rows)} native-GC3 rows; {len(joined)-len(gc3_rows)} missing" if row["n"] != "0" else row["denominator"]
        results.append(row)
    apply_bh(results)
    return results


def write_results(rows: Sequence[Mapping[str, Any]], path: Path | str) -> None:
    path = Path(path)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            RESULT_FIELDS,
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        for source in rows:
            row = empty_result_row()
            row.update({key: _format(value) for key, value in source.items() if key in row})
            writer.writerow(row)


# A compact 5x7 raster font used only to keep the committed PNG reproducible
# inside the pinned headless Guix environment.
FONT = {
    "A": ("01110", "10001", "10001", "11111", "10001", "10001", "10001"),
    "B": ("11110", "10001", "10001", "11110", "10001", "10001", "11110"),
    "C": ("01111", "10000", "10000", "10000", "10000", "10000", "01111"),
    "D": ("11110", "10001", "10001", "10001", "10001", "10001", "11110"),
    "E": ("11111", "10000", "10000", "11110", "10000", "10000", "11111"),
    "F": ("11111", "10000", "10000", "11110", "10000", "10000", "10000"),
    "G": ("01111", "10000", "10000", "10111", "10001", "10001", "01110"),
    "H": ("10001", "10001", "10001", "11111", "10001", "10001", "10001"),
    "I": ("11111", "00100", "00100", "00100", "00100", "00100", "11111"),
    "J": ("00111", "00010", "00010", "00010", "10010", "10010", "01100"),
    "K": ("10001", "10010", "10100", "11000", "10100", "10010", "10001"),
    "L": ("10000", "10000", "10000", "10000", "10000", "10000", "11111"),
    "M": ("10001", "11011", "10101", "10101", "10001", "10001", "10001"),
    "N": ("10001", "11001", "10101", "10011", "10001", "10001", "10001"),
    "O": ("01110", "10001", "10001", "10001", "10001", "10001", "01110"),
    "P": ("11110", "10001", "10001", "11110", "10000", "10000", "10000"),
    "Q": ("01110", "10001", "10001", "10001", "10101", "10010", "01101"),
    "R": ("11110", "10001", "10001", "11110", "10100", "10010", "10001"),
    "S": ("01111", "10000", "10000", "01110", "00001", "00001", "11110"),
    "T": ("11111", "00100", "00100", "00100", "00100", "00100", "00100"),
    "U": ("10001", "10001", "10001", "10001", "10001", "10001", "01110"),
    "V": ("10001", "10001", "10001", "10001", "10001", "01010", "00100"),
    "W": ("10001", "10001", "10001", "10101", "10101", "10101", "01010"),
    "X": ("10001", "10001", "01010", "00100", "01010", "10001", "10001"),
    "Y": ("10001", "10001", "01010", "00100", "00100", "00100", "00100"),
    "Z": ("11111", "00001", "00010", "00100", "01000", "10000", "11111"),
    "0": ("01110", "10001", "10011", "10101", "11001", "10001", "01110"),
    "1": ("00100", "01100", "00100", "00100", "00100", "00100", "01110"),
    "2": ("01110", "10001", "00001", "00010", "00100", "01000", "11111"),
    "3": ("11110", "00001", "00001", "01110", "00001", "00001", "11110"),
    "4": ("00010", "00110", "01010", "10010", "11111", "00010", "00010"),
    "5": ("11111", "10000", "10000", "11110", "00001", "00001", "11110"),
    "6": ("01110", "10000", "10000", "11110", "10001", "10001", "01110"),
    "7": ("11111", "00001", "00010", "00100", "01000", "01000", "01000"),
    "8": ("01110", "10001", "10001", "01110", "10001", "10001", "01110"),
    "9": ("01110", "10001", "10001", "01111", "00001", "00001", "01110"),
    ".": ("00000", "00000", "00000", "00000", "00000", "00110", "00110"),
    "-": ("00000", "00000", "00000", "11111", "00000", "00000", "00000"),
    "/": ("00001", "00010", "00010", "00100", "01000", "01000", "10000"),
    " ": ("00000",) * 7,
}


class Canvas:
    def __init__(self, width: int, height: int):
        self.width, self.height = width, height
        self.data = bytearray([255, 255, 255] * width * height)

    def pixel(self, x: int, y: int, color: Tuple[int, int, int]) -> None:
        if 0 <= x < self.width and 0 <= y < self.height:
            offset = (y * self.width + x) * 3
            self.data[offset : offset + 3] = bytes(color)

    def line(self, x0: int, y0: int, x1: int, y1: int, color=(45, 45, 45), width=1) -> None:
        dx, dy = abs(x1 - x0), -abs(y1 - y0)
        sx, sy = (1 if x0 < x1 else -1), (1 if y0 < y1 else -1)
        error = dx + dy
        while True:
            for ox in range(-(width // 2), width // 2 + 1):
                for oy in range(-(width // 2), width // 2 + 1):
                    self.pixel(x0 + ox, y0 + oy, color)
            if x0 == x1 and y0 == y1:
                break
            twice = 2 * error
            if twice >= dy:
                error += dy
                x0 += sx
            if twice <= dx:
                error += dx
                y0 += sy

    def circle(self, x: int, y: int, radius: int, color: Tuple[int, int, int]) -> None:
        for yy in range(y - radius, y + radius + 1):
            for xx in range(x - radius, x + radius + 1):
                if (xx - x) ** 2 + (yy - y) ** 2 <= radius**2:
                    self.pixel(xx, yy, color)

    def text(self, x: int, y: int, value: str, color=(25, 25, 25), scale=2) -> None:
        cursor = x
        for character in value.upper():
            glyph = FONT.get(character, FONT[" "])
            for gy, row in enumerate(glyph):
                for gx, bit in enumerate(row):
                    if bit == "1":
                        for oy in range(scale):
                            for ox in range(scale):
                                self.pixel(cursor + gx * scale + ox, y + gy * scale + oy, color)
            cursor += 6 * scale

    def png(self) -> bytes:
        raw = b"".join(b"\x00" + bytes(self.data[y * self.width * 3 : (y + 1) * self.width * 3]) for y in range(self.height))

        def chunk(kind: bytes, payload: bytes) -> bytes:
            return struct.pack(">I", len(payload)) + kind + payload + struct.pack(">I", zlib.crc32(kind + payload) & 0xFFFFFFFF)

        return b"\x89PNG\r\n\x1a\n" + chunk(b"IHDR", struct.pack(">IIBBBBB", self.width, self.height, 8, 2, 0, 0, 0)) + chunk(b"IDAT", zlib.compress(raw, 9)) + chunk(b"IEND", b"")


COLORS = [
    (35, 112, 160),
    (222, 105, 45),
    (47, 145, 90),
    (165, 75, 150),
    (190, 150, 40),
    (70, 70, 70),
    (50, 160, 175),
    (210, 80, 95),
]


def _plot_points(rows: Sequence[Mapping[str, str]], observable: str) -> List[Tuple[float, float, str]]:
    points = []
    for row in rows:
        if row.get("row_kind") == "point" and row.get("observable") == observable and row.get("status") == "estimated":
            x, y = _finite(row.get("predictor_value")), _finite(row.get("estimate"))
            if x is not None and y is not None:
                points.append((x, y, row.get("clade", "unclassified")))
    return points


def _panel(
    canvas: Canvas,
    bounds: Tuple[int, int, int, int],
    points,
    title: str,
    reported_effect: Optional[float] = None,
) -> None:
    left, top, right, bottom = bounds
    canvas.line(left, bottom, right, bottom, width=2)
    canvas.line(left, top, left, bottom, width=2)
    canvas.text(left, top - 28, title, scale=2)
    if not points:
        canvas.text(left + 60, (top + bottom) // 2, "UNAVAILABLE", color=(170, 40, 40), scale=3)
        return
    xs, ys = [p[0] for p in points], [p[1] for p in points]
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)
    if xmax == xmin:
        xmax += 1
    if ymax == ymin:
        ymax += 0.1
    levels = {level: COLORS[index % len(COLORS)] for index, level in enumerate(sorted({p[2] for p in points}))}
    for x, y, level in points:
        px = int(left + 8 + (x - xmin) / (xmax - xmin) * (right - left - 16))
        py = int(bottom - 8 - (y - ymin) / (ymax - ymin) * (bottom - top - 16))
        canvas.circle(px, py, 4, levels[level])
    if reported_effect is None:
        design = np.column_stack((np.ones(len(xs)), np.asarray(xs)))
        coefficient, *_ = np.linalg.lstsq(design, np.asarray(ys), rcond=None)
        intercept, slope = coefficient
    else:
        # Display the reported class-fixed slope through the point centroid.
        # A raw pooled line can reverse sign under deep-clade confounding and
        # would not represent the model summarized in the result table.
        slope = reported_effect
        intercept = float(np.mean(ys) - slope * np.mean(xs))
    for xa, xb in zip(np.linspace(xmin, xmax, 60)[:-1], np.linspace(xmin, xmax, 60)[1:]):
        ya, yb = intercept + slope * xa, intercept + slope * xb
        x0 = int(left + 8 + (xa - xmin) / (xmax - xmin) * (right - left - 16))
        x1 = int(left + 8 + (xb - xmin) / (xmax - xmin) * (right - left - 16))
        y0 = int(bottom - 8 - (ya - ymin) / (ymax - ymin) * (bottom - top - 16))
        y1 = int(bottom - 8 - (yb - ymin) / (ymax - ymin) * (bottom - top - 16))
        canvas.line(x0, y0, x1, y1, color=(15, 15, 15), width=2)
    canvas.text(left, bottom + 12, f"LOG10 N  N={len(points)}", scale=2)


def _pdf_bytes(rows: Sequence[Mapping[str, str]], width=1000, height=580) -> bytes:
    gc3 = _plot_points(rows, "gc3")
    whole = _plot_points(rows, "whole_genome_gc")
    commands = ["1 1 1 rg 0 0 1000 580 re f", "0.15 0.15 0.15 RG 1 w"]
    commands.append("BT /F1 18 Tf 45 545 Td (Tier 3 exact-species composition synthesis) Tj ET")
    effects = {
        "gc3": _reported_effect(rows, "primary.gc3.across_class_fixed"),
        "whole_genome_gc": _reported_effect(rows, "sensitivity.whole_genome_gc.across_class_fixed"),
    }
    for points, x0, title, observable in [
        (gc3, 55, "GC3 native annotation", "gc3"),
        (whole, 535, "Whole-genome GC sensitivity", "whole_genome_gc"),
    ]:
        commands.append(f"{x0} 95 405 390 re S")
        commands.append(f"BT /F1 13 Tf {x0} 505 Td ({title}) Tj ET")
        if points:
            xs, ys = [p[0] for p in points], [p[1] for p in points]
            xmin, xmax = min(xs), max(xs)
            ymin, ymax = min(ys), max(ys)
            xmax = xmax if xmax != xmin else xmin + 1
            ymax = ymax if ymax != ymin else ymin + 0.1
            levels = {level: index for index, level in enumerate(sorted({p[2] for p in points}))}
            for x, y, level in points:
                r, g, b = [component / 255 for component in COLORS[levels[level] % len(COLORS)]]
                px = x0 + 8 + (x - xmin) / (xmax - xmin) * 389
                py = 103 + (y - ymin) / (ymax - ymin) * 374
                commands.append(f"{r:.3f} {g:.3f} {b:.3f} rg {px-2.5:.2f} {py-2.5:.2f} 5 5 re f")
            slope = effects[observable]
            if slope is None:
                design = np.column_stack((np.ones(len(xs)), np.asarray(xs)))
                coefficient, *_ = np.linalg.lstsq(design, np.asarray(ys), rcond=None)
                intercept, slope = coefficient
            else:
                intercept = float(np.mean(ys) - slope * np.mean(xs))
            ya, yb = intercept + slope * xmin, intercept + slope * xmax
            pya = 103 + (ya - ymin) / (ymax - ymin) * 374
            pyb = 103 + (yb - ymin) / (ymax - ymin) * 374
            commands.append(f"0 0 0 RG 2 w {x0+8} {pya:.2f} m {x0+397} {pyb:.2f} l S")
        commands.append(f"0 0 0 rg BT /F1 10 Tf {x0} 70 Td (Buffalo pred_log10_N; N={len(points)}) Tj ET")
    commands.append("BT /F1 11 Tf 55 35 Td (Diversity panels unavailable: zero eligible population or individual tuples; polarized SFS-B deferred.) Tj ET")
    stream = ("\n".join(commands) + "\n").encode("ascii")
    objects = [
        b"<< /Type /Catalog /Pages 2 0 R >>",
        b"<< /Type /Pages /Kids [3 0 R] /Count 1 >>",
        b"<< /Type /Page /Parent 2 0 R /MediaBox [0 0 1000 580] /Resources << /Font << /F1 5 0 R >> >> /Contents 4 0 R >>",
        b"<< /Length " + str(len(stream)).encode() + b" >>\nstream\n" + stream + b"endstream",
        b"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica >>",
    ]
    output = bytearray(b"%PDF-1.4\n% reproducible Tier 3 figure\n")
    offsets = [0]
    for index, obj in enumerate(objects, 1):
        offsets.append(len(output))
        output.extend(f"{index} 0 obj\n".encode() + obj + b"\nendobj\n")
    xref = len(output)
    output.extend(f"xref\n0 {len(objects)+1}\n0000000000 65535 f\n".encode())
    for offset in offsets[1:]:
        output.extend(f"{offset:010d} 00000 n\n".encode())
    output.extend(f"trailer << /Size {len(objects)+1} /Root 1 0 R >>\nstartxref\n{xref}\n%%EOF\n".encode())
    return bytes(output)


def render_figure(
    rows: Sequence[Mapping[str, str]], png_path: Path | str, pdf_path: Path | str
) -> None:
    canvas = Canvas(1400, 820)
    canvas.text(55, 28, "TIER 3 EXACT SPECIES COMPOSITION SYNTHESIS", scale=3)
    _panel(
        canvas,
        (90, 125, 675, 690),
        _plot_points(rows, "gc3"),
        "GC3 NATIVE ANNOTATION",
        _reported_effect(rows, "primary.gc3.across_class_fixed"),
    )
    _panel(
        canvas,
        (790, 125, 1375, 690),
        _plot_points(rows, "whole_genome_gc"),
        "WHOLE GENOME GC SENSITIVITY",
        _reported_effect(rows, "sensitivity.whole_genome_gc.across_class_fixed"),
    )
    canvas.text(90, 750, "DIVERSITY UNAVAILABLE   SFS B DEFERRED   NO CAUSAL CLAIM", color=(145, 35, 35), scale=2)
    Path(png_path).write_bytes(canvas.png())
    Path(pdf_path).write_bytes(_pdf_bytes(rows))


def _reported_effect(rows: Sequence[Mapping[str, str]], analysis_id: str) -> Optional[float]:
    matches = [row for row in rows if row.get("analysis_id") == analysis_id]
    if len(matches) != 1:
        return None
    return _finite(matches[0].get("effect"))


def _arguments(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tier3a", type=Path, default=ANALYSIS / "tier3a_data.tsv")
    parser.add_argument("--tier3b", type=Path, default=ANALYSIS / "tier3b_data.tsv")
    parser.add_argument("--tier3c", type=Path, default=ANALYSIS / "tier3c_data.tsv")
    parser.add_argument("--buffalo", type=Path, help="verified pinned combined_data.tsv; fetched if omitted")
    parser.add_argument("--results", type=Path, default=ANALYSIS / "tier3_results.tsv")
    parser.add_argument("--png", type=Path, default=ANALYSIS / "fig_tier3.png")
    parser.add_argument("--pdf", type=Path, default=ANALYSIS / "fig_tier3.pdf")
    parser.add_argument("--figure-from-results", type=Path, help="skip joins/fits and regenerate figures from a committed result TSV")
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _arguments(argv)
    if args.figure_from_results:
        render_figure(read_tsv(args.figure_from_results), args.png, args.pdf)
        return 0
    buffalo_path = args.buffalo
    temporary = None
    if buffalo_path is None:
        import tempfile

        handle = tempfile.NamedTemporaryFile(prefix="tier3-buffalo-", suffix=".tsv", delete=False)
        handle.close()
        temporary = Path(handle.name)
        buffalo_path = fetch_buffalo(temporary)
    else:
        digest = hashlib.sha256(buffalo_path.read_bytes()).hexdigest()
        if digest != BUFFALO_SHA256:
            raise SynthesisError(f"Buffalo SHA-256 mismatch: {digest}")
    try:
        joined = join_composition(
            read_tsv(args.tier3c),
            read_buffalo(buffalo_path),
            load_composition_provenance(ANALYSIS / "tier3c_qc"),
        )
        observations = load_diversity_observations(args.tier3a, args.tier3b)
        results = synthesize(joined, observations)
        write_results(results, args.results)
        render_figure(results, args.png, args.pdf)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except SynthesisError as error:
        print(f"tier3_fit: {error}", file=sys.stderr)
        raise SystemExit(2)
