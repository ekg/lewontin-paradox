#!/usr/bin/env python3
"""Synthesize the measured VGP scale-out and raw-read validation evidence.

This module is deliberately report-only.  It reconciles committed upstream
packets, writes paper-facing tables and static SVG figures, and preserves the
covariance and confidence limits that govern interpretation.  It never submits
jobs, downloads data, or promotes canonical biological objects.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import os
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable, Mapping, Sequence
from xml.sax.saxutils import escape


ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = ROOT / "analysis"

ROOT_CONFIG = ANALYSIS / "vgp_data_root_config.json"
PILOT_WORKFLOW = ANALYSIS / "vgp_10_pilot_workflow.json"
PILOT_ENVIRONMENT_LOCK = ANALYSIS / "guix/vgp_10_pilot/environment-lock.json"
PILOT_ENVIRONMENT_REALIZATION = ANALYSIS / "guix/vgp_10_pilot/realization.json"
SCALE_SUMMARY = ANALYSIS / "vgp_real_scaleout_v1/summary.json"
CATALOG_ACCOUNTING = ANALYSIS / "vgp_real_scaleout_v1/catalog_accounting.tsv"
PAIR_ACCOUNTING = ANALYSIS / "vgp_real_scaleout_v1/pair_accounting.tsv"
SCALE_TELEMETRY = ANALYSIS / "vgp_real_scaleout_v1/scaleout_telemetry.tsv"
SCALE_SACCT = ANALYSIS / "vgp_real_scaleout_sacct_v1.tsv"
PSMC_TRAJECTORIES = ANALYSIS / "vgp_real_scaleout_v1/psmc_unscaled_trajectories.tsv"
PSMC_SCENARIOS = ANALYSIS / "vgp_real_scaleout_v1/scenario_uncertainty.tsv"
ANNOTATION_PARTITIONS = ANALYSIS / "vgp_real_scaleout_v1/exact_native_annotation_partitions.tsv"
SCRATCH_CONTRACTS = ANALYSIS / "vgp_real_scaleout_v1/fastga_scratch_contracts.tsv"
RESOURCE_PLAN = ANALYSIS / "vgp_real_scaleout_v1/per_pair_resource_plan.json"
SCALE_REPORT = ANALYSIS / "vgp_real_scaleout_v1/results.md"
READ_RESULTS = ANALYSIS / "vgp_read_validation_results_v1.json"
READ_PAIRS = ANALYSIS / "vgp_read_validation_per_pair_v1.tsv"
READ_MASKS = ANALYSIS / "vgp_read_validation_mask_sensitivity_v1.tsv"
READ_SACCT = ANALYSIS / "vgp_read_validation_sacct_v1.tsv"
READ_EVIDENCE_MANIFEST = ANALYSIS / "vgp_read_validation_evidence_manifest_v1.json"
READ_OBJECT_MANIFEST = ANALYSIS / "vgp_validation_reads_manifest_v1.json"
READ_ENVIRONMENT = ANALYSIS / "vgp_read_validation_environment_v1.json"
READ_REPORT = ANALYSIS / "vgp_read_validation_report_v1.md"
PRIOR_GENE_CONVERSION = ANALYSIS / "vgp_comprehensive_table_gene_conversion.tsv"
METHOD = ANALYSIS / "synthesize_vgp_real.py"
METHOD_TEST = ANALYSIS / "tests/test_synthesize_vgp_real.py"

INPUT_PATHS = (
    ROOT_CONFIG,
    PILOT_WORKFLOW,
    PILOT_ENVIRONMENT_LOCK,
    PILOT_ENVIRONMENT_REALIZATION,
    SCALE_SUMMARY,
    CATALOG_ACCOUNTING,
    PAIR_ACCOUNTING,
    SCALE_TELEMETRY,
    SCALE_SACCT,
    PSMC_TRAJECTORIES,
    PSMC_SCENARIOS,
    ANNOTATION_PARTITIONS,
    SCRATCH_CONTRACTS,
    RESOURCE_PLAN,
    SCALE_REPORT,
    READ_RESULTS,
    READ_PAIRS,
    READ_MASKS,
    READ_SACCT,
    READ_EVIDENCE_MANIFEST,
    READ_OBJECT_MANIFEST,
    READ_ENVIRONMENT,
    READ_REPORT,
    PRIOR_GENE_CONVERSION,
    METHOD,
    METHOD_TEST,
)

OUTPUT_DIR = ANALYSIS / "vgp_real_synthesis_v1"
REPORT = OUTPUT_DIR / "report.md"
PAPER_PAIRS = OUTPUT_DIR / "paper_pairs.tsv"
PSMC_HISTORIES = OUTPUT_DIR / "psmc_histories.tsv"
ANNOTATION_TABLE = OUTPUT_DIR / "annotation_partitions.tsv"
GENE_CONVERSION = OUTPUT_DIR / "gene_conversion_branches.tsv"
CLAIMS = OUTPUT_DIR / "claim_ledger.tsv"
JOB_LEDGER = OUTPUT_DIR / "job_ledger.tsv"
DIGEST_LEDGER = OUTPUT_DIR / "digest_ledger.tsv"
DIVERSITY_FIGURE = OUTPUT_DIR / "figure_diversity.svg"
PSMC_FIGURE = OUTPUT_DIR / "figure_psmc.svg"
MANIFEST = OUTPUT_DIR / "manifest.json"

OUTPUT_PATHS = (
    REPORT,
    PAPER_PAIRS,
    PSMC_HISTORIES,
    ANNOTATION_TABLE,
    GENE_CONVERSION,
    CLAIMS,
    JOB_LEDGER,
    DIGEST_LEDGER,
    DIVERSITY_FIGURE,
    PSMC_FIGURE,
)
OUTPUT_FILENAMES = tuple(path.name for path in OUTPUT_PATHS) + (MANIFEST.name,)

PAPER_FIELDS = (
    "canonical_vgp_root",
    "selection_id",
    "scientific_name",
    "h1_accession_version",
    "h2_accession_version",
    "upstream_confidence_tier",
    "execution_disposition",
    "hard_failure_class",
    "attempted_slurm_jobs",
    "completed_slurm_allocations",
    "callable_bp",
    "heterozygous_snps",
    "assembly_pi",
    "bootstrap_theta_q025",
    "primary_theta_0",
    "bootstrap_theta_q975",
    "raw_validation_status",
    "primary_common_callable_bp",
    "primary_common_assembly_pi",
    "primary_read_pi",
    "primary_read_over_assembly_pi",
    "kmer_heterozygosity",
    "illumina_contradiction_lower_bound",
    "hifi_contradiction_lower_bound",
    "synthesis_confidence_tier",
    "quantitative_disposition",
    "psmc_disposition",
    "annotation_disposition",
    "interpretive_scope",
)

PSMC_FIELDS = (
    "canonical_vgp_root",
    "selection_id",
    "scientific_name",
    "interval_count",
    "primary_theta_0_per_100bp_bin",
    "bootstrap_theta_q025",
    "bootstrap_theta_median",
    "bootstrap_theta_q975",
    "trajectory_time_2N0_min",
    "trajectory_time_2N0_max",
    "trajectory_lambda_start",
    "trajectory_lambda_end",
    "trajectory_lambda_min",
    "trajectory_lambda_min_interval",
    "trajectory_lambda_min_time_2N0",
    "trajectory_lambda_max",
    "trajectory_lambda_max_interval",
    "trajectory_lambda_max_time_2N0",
    "generic_scenario_count",
    "generic_time_years_min",
    "generic_time_years_max",
    "generic_effective_size_min",
    "generic_effective_size_max",
    "absolute_history_status",
    "pi_psmc_independence",
    "quantitative_disposition",
    "history_interpretation",
)

ANNOTATION_FIELDS = (
    "canonical_vgp_root",
    "selection_id",
    "assembly_accession_version",
    "annotation_accession_version",
    "annotation_gff_sha256",
    "sequence_dictionary_equal",
    "slurm_job_id",
    "partition",
    "callable_bp",
    "heterozygous_variants",
    "estimate",
    "estimator",
    "quantitative_disposition",
    "gene_conversion_disposition",
)

GENE_CONVERSION_FIELDS = (
    "canonical_vgp_root",
    "branch",
    "prior_execution_state",
    "vgp_integration_status",
    "estimate",
    "sampling_unit",
    "reason",
)

CLAIM_FIELDS = (
    "canonical_vgp_root",
    "claim_id",
    "classification",
    "conclusion",
    "estimand",
    "sampling_unit",
    "observed_scope",
    "evidence_artifacts",
    "uncertainty_covariance",
    "forbidden_inference",
)

JOB_FIELDS = (
    "canonical_vgp_root",
    "packet",
    "source_row",
    "job_id",
    "job_name",
    "partition",
    "state",
    "elapsed",
    "elapsed_seconds",
    "timelimit",
    "allocated_cpus",
    "requested_memory",
    "max_rss",
    "cpu_time_raw",
    "total_cpu",
    "exit_code",
    "node_list",
    "start",
    "end",
    "submit",
    "failure_or_nonterminal",
)

DIGEST_FIELDS = (
    "canonical_vgp_root",
    "binding_class",
    "container_path",
    "locator",
    "digest_field",
    "algorithm",
    "digest",
    "occurrences",
    "bytes",
    "verification",
)

SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
MD5_RE = re.compile(r"^[0-9a-f]{32}$")


class SynthesisError(RuntimeError):
    """A closed-world or interpretation invariant failed."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise SynthesisError(message)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_json(path: Path) -> dict[str, object]:
    value = json.loads(path.read_text(encoding="utf-8"))
    require(isinstance(value, dict), f"expected JSON object: {path}")
    return value


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.partial-{os.getpid()}")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        handle.write(text)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def atomic_tsv(path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, object]]) -> None:
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    atomic_text(path, buffer.getvalue())


def atomic_json(path: Path, value: Mapping[str, object]) -> None:
    atomic_text(path, json.dumps(value, indent=2, sort_keys=True) + "\n")


def relative(path: Path) -> str:
    return path.relative_to(ROOT).as_posix()


def fmt_float(value: object) -> str:
    return format(float(value), ".17g")


def _assert_json_roots(value: object, canonical_root: str, pointer: str = "") -> None:
    if isinstance(value, dict):
        for key, child in value.items():
            child_pointer = f"{pointer}/{key}"
            if key in {"canonical_vgp_root", "canonical_root"} and isinstance(child, str):
                require(child == canonical_root, f"canonical root drift at {child_pointer}: {child}")
            _assert_json_roots(child, canonical_root, child_pointer)
    elif isinstance(value, list):
        for index, child in enumerate(value):
            _assert_json_roots(child, canonical_root, f"{pointer}/{index}")


def _state_counts(rows: Sequence[Mapping[str, str]]) -> dict[str, int]:
    return dict(sorted(Counter(row["State"] for row in rows).items()))


def reconcile_inputs() -> dict[str, object]:
    for path in INPUT_PATHS:
        require(path.is_file(), f"missing synthesis input: {path}")

    root_config = read_json(ROOT_CONFIG)
    pilot_workflow = read_json(PILOT_WORKFLOW)
    pilot_lock = read_json(PILOT_ENVIRONMENT_LOCK)
    pilot_realization = read_json(PILOT_ENVIRONMENT_REALIZATION)
    canonical_root = str(root_config["root"])
    require(Path(canonical_root).is_absolute(), "configured VGP root must be absolute")
    require("lewontin-paradox" not in canonical_root, "canonical VGP root cannot be project-named")
    require(root_config["migration_input_only"] == "/moosefs/erikg/lewontin-paradox-data/vgp/phase1-freeze-1.0",
            "legacy root must remain migration input only")

    scale = read_json(SCALE_SUMMARY)
    reads = read_json(READ_RESULTS)
    objects = read_json(READ_OBJECT_MANIFEST)
    environment = read_json(READ_ENVIRONMENT)
    evidence_manifest = read_json(READ_EVIDENCE_MANIFEST)
    resource_plan = read_json(RESOURCE_PLAN)
    catalog_rows = read_tsv(CATALOG_ACCOUNTING)
    pair_rows = read_tsv(PAIR_ACCOUNTING)
    telemetry_rows = read_tsv(SCALE_TELEMETRY)
    scale_jobs = read_tsv(SCALE_SACCT)
    read_jobs = read_tsv(READ_SACCT)
    trajectory_rows = read_tsv(PSMC_TRAJECTORIES)
    scenario_rows = read_tsv(PSMC_SCENARIOS)
    annotation_rows = read_tsv(ANNOTATION_PARTITIONS)
    scratch_rows = read_tsv(SCRATCH_CONTRACTS)
    read_pair_rows = read_tsv(READ_PAIRS)
    read_mask_rows = read_tsv(READ_MASKS)
    prior_gc_rows = read_tsv(PRIOR_GENE_CONVERSION)

    for value in (scale, reads, objects, environment, evidence_manifest, resource_plan):
        _assert_json_roots(value, canonical_root)
    for path, rows in (
        (CATALOG_ACCOUNTING, catalog_rows),
        (PAIR_ACCOUNTING, pair_rows),
        (SCALE_TELEMETRY, telemetry_rows),
        (SCALE_SACCT, scale_jobs),
        (READ_SACCT, read_jobs),
        (PSMC_TRAJECTORIES, trajectory_rows),
        (PSMC_SCENARIOS, scenario_rows),
        (ANNOTATION_PARTITIONS, annotation_rows),
        (SCRATCH_CONTRACTS, scratch_rows),
        (READ_PAIRS, read_pair_rows),
        (READ_MASKS, read_mask_rows),
    ):
        require(rows and all(row.get("canonical_vgp_root") == canonical_root for row in rows),
                f"canonical-root reconciliation failed: {path}")

    require(len(catalog_rows) == 716, "Freeze 1 catalog accounting must contain 716 rows")
    require(len(pair_rows) == 569, "Freeze 1 link accounting must contain 569 rows")
    eligible = [row for row in pair_rows if row["selection_id"]]
    require(len(eligible) == 10 and {row["selection_id"] for row in eligible} == {f"P{i:02d}" for i in range(1, 11)},
            "audited ten-pair roster drift")
    require(scale["catalog"] == {
        "commit": "dc1b2af5a7741b97d66fb10cb2bce97f41765cdf",
        "links": 569,
        "rows": 716,
        "sha256": "9c58420484a8b76a2d6175b7c26bf709e68bdc726a67fc7541b8c2b5a2fc13a4",
    }, "frozen catalog identity drift")
    dispositions = Counter(row["execution_disposition"] for row in eligible)
    require(dispositions == {
        "HARD_EXECUTION_ERROR_NO_ESTIMATE": 3,
        "HARD_INVALID_PRIMARY": 4,
        "RUNNING_RESUMABLE_WAVE": 1,
        "VERIFIED_CORE_COMPLETE": 2,
    }, "ten-pair execution disposition drift")
    measured = [row for row in eligible if row["callable_pi"]]
    require({row["selection_id"] for row in measured} == {"P04", "P07"},
            "completed biological result roster drift")
    require(all(int(row["finite_bootstraps"]) == int(row["bootstrap_attempts"]) == 200 for row in measured),
            "PSMC bootstrap completion drift")

    require(len(scale_jobs) == 650 and len({row["JobIDRaw"] for row in scale_jobs}) == 650,
            "scale-out sacct closed world drift")
    require(_state_counts(scale_jobs) == {
        "CANCELLED by 1001": 346, "COMPLETED": 240, "FAILED": 42, "PENDING": 21, "RUNNING": 1,
    }, "scale-out sacct state drift")
    require(len(read_jobs) == 8 and len({row["JobIDRaw"] for row in read_jobs}) == 8,
            "read-validation sacct closed world drift")
    require(_state_counts(read_jobs) == {
        "CANCELLED by 1001": 2, "COMPLETED": 2, "FAILED": 4,
    }, "read-validation sacct state drift")
    require(not ({row["JobIDRaw"] for row in scale_jobs} & {row["JobIDRaw"] for row in read_jobs}),
            "scheduler packets unexpectedly overlap")
    telemetry_all = next(row for row in telemetry_rows if row["scope"] == "pilot_and_scaleout_attempts")
    require(int(telemetry_all["allocations"]) == len(scale_jobs), "scale telemetry/sacct allocation mismatch")

    require(len(trajectory_rows) == 128 and Counter(row["selection_id"] for row in trajectory_rows) == {"P04": 64, "P07": 64},
            "unscaled PSMC trajectory reconciliation drift")
    require(len(scenario_rows) == 1152 and Counter(row["selection_id"] for row in scenario_rows) == {"P04": 576, "P07": 576},
            "generic PSMC scenario reconciliation drift")
    require(all(row["mutation_rate_source"] == "predeclared_generic_sensitivity_grid_not_species_calibration" and
                row["generation_time_source"] == "predeclared_generic_sensitivity_grid_not_species_calibration"
                for row in scenario_rows), "generic scenarios were promoted to calibrations")
    require(len(annotation_rows) == 6 and {row["selection_id"] for row in annotation_rows} == {"P07"} and
            all(row["sequence_dictionary_equal"] == "true" for row in annotation_rows),
            "exact-native annotation partition drift")
    require(len(scratch_rows) == 2 and all(row.get("contract_status", row.get("status", "")).startswith("PASS")
                                           for row in scratch_rows),
            "FastGA scratch contracts are not both passing")

    require({row["selection_id"] for row in read_pair_rows} == {"P04", "P07", "P09"},
            "raw-validation pair roster drift")
    require(len(read_mask_rows) == 6 and {row["selection_id"] for row in read_mask_rows} == {"P07"},
            "raw-validation mask sensitivity drift")
    p07 = reads["pairs"]["P07"]
    require(p07["validation_assessment"]["classification"] == "concrete_haplotype_reconstruction_failure" and
            p07["validation_assessment"]["cross_technology_majority_homozygous_reference_contradiction"] is True,
            "P07 concrete reconstruction failure drift")
    require(p07["validation_assessment"]["downstream_action"] ==
            "preserve the core artifact and provenance, but do not use its pi or PSMC as validated quantitative evidence",
            "P07 downstream action drift")
    require(reads["pairs"]["P04"]["raw_validation_status"] == "not_estimable_pending_raw_reads",
            "P04 validation-pending status drift")
    require(reads["pairs"]["P09"]["raw_validation_status"] == "low_coverage_mapping_control_only",
            "P09 compatibility-control status drift")
    require(environment["frozen_channel_commit"] == "44bbfc24e4bcc48d0e3343cd3d83452721af8c36" and
            environment["closure_sha256"] == "ac0cb3601e56ef62b9ef99419de3659b2a2ba59b2aead29bc5f1928b50c83da2",
            "pinned Guix environment drift")
    require(pilot_workflow["environment"] == {
        "ambient_executables_allowed": False,
        "channels": "analysis/guix/vgp_10_pilot/channels.scm",
        "lock": "analysis/guix/vgp_10_pilot/environment-lock.json",
        "manifest": "analysis/guix/vgp_10_pilot/manifest.scm",
        "realization": "analysis/guix/vgp_10_pilot/realization.json",
    }, "core workflow environment binding drift")
    require(pilot_lock["channel_commit"] == pilot_realization["channel_commit"] ==
            environment["frozen_channel_commit"] and
            pilot_realization["closure_sha256"] ==
            "8fcdb32021f1cd8eac839509cff47ab6bdd63b656b30e243fdf78d3c4ba24f9d",
            "realized core Guix environment drift")
    require(objects["legacy_root_role"] == "migration_input_only" and
            objects["summary"]["verified"] == {"bytes": 31058137613, "objects": 4} and
            objects["summary"]["pending"] == {"bytes": 42344746693, "objects": 1},
            "raw-read acquisition accounting drift")
    require([row["branch"] for row in prior_gc_rows] == [
        "direct_pedigree_or_gamete", "population_allele_frequency_spectrum",
        "historical_phylogenetic_substitution", "non_allelic_paralog",
    ] and all(row["estimate"] == "NOT_ESTIMABLE" for row in prior_gc_rows),
            "gene-conversion branch contract changed; inspect before integration")

    return {
        "canonical_root": canonical_root,
        "scale": scale,
        "reads": reads,
        "objects": objects,
        "environment": environment,
        "pilot_environment": pilot_realization,
        "catalog_rows": catalog_rows,
        "pair_rows": pair_rows,
        "eligible": eligible,
        "measured": measured,
        "scale_jobs": scale_jobs,
        "read_jobs": read_jobs,
        "trajectory_rows": trajectory_rows,
        "scenario_rows": scenario_rows,
        "annotation_rows": annotation_rows,
        "read_pair_rows": read_pair_rows,
        "read_mask_rows": read_mask_rows,
        "prior_gc_rows": prior_gc_rows,
        "dispositions": dispositions,
    }


def build_pair_rows(data: Mapping[str, object]) -> list[dict[str, str]]:
    canonical_root = str(data["canonical_root"])
    read_by_id = {row["selection_id"]: row for row in data["read_pair_rows"]}
    reads = data["reads"]
    output: list[dict[str, str]] = []
    for row in sorted(data["eligible"], key=lambda item: item["selection_id"]):
        pair_id = row["selection_id"]
        raw = read_by_id.get(pair_id, {})
        common_bp = common_assembly = read_pi = pi_ratio = kmer_het = illumina = hifi = ""
        if pair_id == "P07":
            p07 = reads["pairs"]["P07"]
            common = p07["primary_common_mask_comparison"]
            evidence = p07["read_backed_assembly_sites"]
            # Preserve the committed paper-table decimal spellings rather than
            # round-tripping their binary JSON representations.
            common_bp = raw["primary_callable_bp"]
            common_assembly = raw["primary_assembly_pi"]
            read_pi = raw["primary_read_pi"]
            pi_ratio = raw["primary_pi_ratio"]
            kmer_het = raw["kmer_heterozygosity"]
            illumina = fmt_float(evidence["illumina"]["concrete_false_positive_lower_bound_fraction"])
            hifi = fmt_float(evidence["hifi"]["concrete_false_positive_lower_bound_fraction"])

        disposition = row["execution_disposition"]
        if pair_id == "P04":
            tier = "T1_RETAINED_RAW_VALIDATION_PENDING"
            quantitative = "RETAINED_ASSEMBLY_DERIVED_PENDING_RAW_VALIDATION"
            psmc = "RETAINED_DESCRIPTIVE_UNSCALED_RAW_PENDING"
            annotation = "NOT_AVAILABLE_NOT_CORE_VETO"
            scope = "one assembly-derived individual result; no population mean; raw validation pending"
        elif pair_id == "P07":
            tier = "T2_INVALIDATED_BY_EXACT_READS"
            quantitative = "PRESERVED_NOT_ADMITTED_QUANTITATIVELY"
            psmc = "INVALIDATED_BY_EXACT_READ_VALIDATION"
            annotation = "EXACT_NATIVE_DESCRIPTIVE_PARENT_PAIR_INVALIDATED"
            scope = "assembly artifact preserved; exact-individual Illumina and HiFi meet concrete reconstruction-failure rule"
        elif disposition == "HARD_INVALID_PRIMARY":
            tier = "X_HARD_INVALID_PRIMARY"
            quantitative = "NOT_ESTIMABLE"
            psmc = "NOT_ESTIMABLE"
            annotation = "NOT_RUN_NO_VALID_CORE"
            scope = "hard-invalid primary execution; absence is not a biological zero"
        elif disposition == "HARD_EXECUTION_ERROR_NO_ESTIMATE":
            tier = "X_EXECUTION_ERROR_NO_ESTIMATE"
            quantitative = "NOT_ESTIMABLE"
            psmc = "NOT_ESTIMABLE"
            annotation = "NOT_RUN_NO_VALID_CORE"
            scope = "concrete execution error; absence is not a biological zero"
        else:
            require(disposition == "RUNNING_RESUMABLE_WAVE", f"unexpected disposition: {disposition}")
            tier = "P_RESUMABLE_RUNNING_AT_FREEZE"
            quantitative = "NOT_ESTIMABLE_AT_FREEZE"
            psmc = "NOT_ESTIMABLE_AT_FREEZE"
            annotation = "NOT_RUN_AT_FREEZE"
            scope = "resumable running wave at closed-world cutoff; no frozen estimate"

        output.append({
            "canonical_vgp_root": canonical_root,
            "selection_id": pair_id,
            "scientific_name": row["scientific_name"],
            "h1_accession_version": row["h1_accession_version"],
            "h2_accession_version": row["h2_accession_version"],
            "upstream_confidence_tier": row["confidence_tier"],
            "execution_disposition": disposition,
            "hard_failure_class": row["hard_failure_class"],
            "attempted_slurm_jobs": row["attempted_slurm_jobs"],
            "completed_slurm_allocations": row["completed_slurm_allocations"],
            "callable_bp": row["callable_bp"],
            "heterozygous_snps": row["heterozygous_snps"],
            "assembly_pi": row["callable_pi"],
            "bootstrap_theta_q025": row["bootstrap_theta_q025"],
            "primary_theta_0": row["primary_theta_0"],
            "bootstrap_theta_q975": row["bootstrap_theta_q975"],
            "raw_validation_status": raw.get("raw_validation_status", "NOT_RUN"),
            "primary_common_callable_bp": common_bp,
            "primary_common_assembly_pi": common_assembly,
            "primary_read_pi": read_pi,
            "primary_read_over_assembly_pi": pi_ratio,
            "kmer_heterozygosity": kmer_het,
            "illumina_contradiction_lower_bound": illumina,
            "hifi_contradiction_lower_bound": hifi,
            "synthesis_confidence_tier": tier,
            "quantitative_disposition": quantitative,
            "psmc_disposition": psmc,
            "annotation_disposition": annotation,
            "interpretive_scope": scope,
        })
    return output


def build_psmc_rows(data: Mapping[str, object]) -> list[dict[str, str]]:
    by_pair: dict[str, list[dict[str, str]]] = defaultdict(list)
    scenarios: dict[str, list[dict[str, str]]] = defaultdict(list)
    eligible = {row["selection_id"]: row for row in data["eligible"]}
    for row in data["trajectory_rows"]:
        by_pair[row["selection_id"]].append(row)
    for row in data["scenario_rows"]:
        scenarios[row["selection_id"]].append(row)

    output: list[dict[str, str]] = []
    for pair_id in ("P04", "P07"):
        rows = sorted(by_pair[pair_id], key=lambda row: int(row["interval"]))
        scenario = scenarios[pair_id]
        lambdas = [float(row["lambda"]) for row in rows]
        times = [float(row["time_2N0"]) for row in rows]
        min_index = min(range(len(rows)), key=lambdas.__getitem__)
        max_index = max(range(len(rows)), key=lambdas.__getitem__)
        pair = eligible[pair_id]
        generic_ids = {row["scenario_id"] for row in scenario}
        if pair_id == "P04":
            disposition = "RETAINED_DESCRIPTIVE_UNSCALED_RAW_PENDING"
            interpretation = (
                "early lambda maximum followed by decline to a broad later low; descriptive shape only, "
                "because absolute scaling is generic and the same H1/H2 pair also supplies pi"
            )
        else:
            disposition = "INVALIDATED_BY_EXACT_READ_VALIDATION"
            interpretation = (
                "large late lambda increase is preserved as an assembly-derived artifact but excluded from "
                "biological demographic interpretation after exact-read reconstruction failure"
            )
        output.append({
            "canonical_vgp_root": str(data["canonical_root"]),
            "selection_id": pair_id,
            "scientific_name": pair["scientific_name"],
            "interval_count": str(len(rows)),
            "primary_theta_0_per_100bp_bin": rows[0]["primary_theta_0_per_100bp_bin"],
            "bootstrap_theta_q025": pair["bootstrap_theta_q025"],
            "bootstrap_theta_median": pair["bootstrap_theta_median"],
            "bootstrap_theta_q975": pair["bootstrap_theta_q975"],
            "trajectory_time_2N0_min": fmt_float(min(times)),
            "trajectory_time_2N0_max": fmt_float(max(times)),
            "trajectory_lambda_start": rows[0]["lambda"],
            "trajectory_lambda_end": rows[-1]["lambda"],
            "trajectory_lambda_min": rows[min_index]["lambda"],
            "trajectory_lambda_min_interval": rows[min_index]["interval"],
            "trajectory_lambda_min_time_2N0": rows[min_index]["time_2N0"],
            "trajectory_lambda_max": rows[max_index]["lambda"],
            "trajectory_lambda_max_interval": rows[max_index]["interval"],
            "trajectory_lambda_max_time_2N0": rows[max_index]["time_2N0"],
            "generic_scenario_count": str(len(generic_ids)),
            "generic_time_years_min": fmt_float(min(float(row["time_years"]) for row in scenario)),
            "generic_time_years_max": fmt_float(max(float(row["time_years"]) for row in scenario)),
            "generic_effective_size_min": fmt_float(min(float(row["effective_size"]) for row in scenario)),
            "generic_effective_size_max": fmt_float(max(float(row["effective_size"]) for row in scenario)),
            "absolute_history_status": "BOUNDED_9_GENERIC_SCENARIOS_NOT_SPECIES_CALIBRATION",
            "pi_psmc_independence": "SAME_PAIR_NONINDEPENDENT",
            "quantitative_disposition": disposition,
            "history_interpretation": interpretation,
        })
    return output


def build_annotation_rows(data: Mapping[str, object]) -> list[dict[str, str]]:
    output = []
    for row in data["annotation_rows"]:
        output.append({
            **{field: row[field] for field in ANNOTATION_FIELDS[:12]},
            "quantitative_disposition": "DESCRIPTIVE_ONLY_PARENT_PAIR_INVALIDATED",
            "gene_conversion_disposition": "NOT_A_CONFORMING_GENE_CONVERSION_ESTIMATE",
        })
    return output


def build_gene_conversion_rows(data: Mapping[str, object]) -> list[dict[str, str]]:
    reasons = {
        "direct_pedigree_or_gamete": "no conforming VGP meiosis or gamete estimate; H1/H2 differences are not transmission events",
        "population_allele_frequency_spectrum": "no conforming VGP population sampling or polarized allele-frequency spectrum",
        "historical_phylogenetic_substitution": "no conforming VGP branch-substitution estimate; annotation WS/SW counts are individual heterozygosity partitions",
        "non_allelic_paralog": "no conforming copy-resolved VGP paralog-tract estimate",
    }
    return [{
        "canonical_vgp_root": str(data["canonical_root"]),
        "branch": row["branch"],
        "prior_execution_state": row["execution_state"],
        "vgp_integration_status": "SEPARATE_NO_ACTUAL_CONFORMING_VGP_ESTIMATE",
        "estimate": "NOT_ESTIMABLE",
        "sampling_unit": row["sampling_unit"],
        "reason": reasons[row["branch"]],
    } for row in data["prior_gc_rows"]]


def build_claim_rows(data: Mapping[str, object]) -> list[dict[str, str]]:
    root = str(data["canonical_root"])
    rows = [
        ("EXEC-CLOSED", "supported", "All 716 catalog rows, 569 links, 10 audited pairs, and 650 scale packet scheduler allocations are accounted at the freeze.", "closed-world execution accounting", "catalog row, link, audited pair, and scheduler allocation", "716 rows; 569 links; 10 pairs; 650 allocations", "catalog_accounting.tsv;pair_accounting.tsv;vgp_real_scaleout_sacct_v1.tsv", "Scheduler states include failed, cancelled, pending, and running records; MaxRSS was unavailable and not imputed.", "Do not equate submitted allocations with independent biological observations."),
        ("DIV-ASSEMBLY", "supported", "Two completed H1/H2 assembly-derived individual estimates are 0.0021472198306856562 and 0.004604184795871289 per callable base, a 2.144254-fold observed range.", "assembly-derived callable heterozygosity", "one phased H1/H2 pair from one individual", "P04 and P07", "pair_accounting.tsv;paper_pairs.tsv", "Both estimates are method-specific; their PSMC outputs share the same pair and are non-independent.", "Do not call this a validated vertebrate distribution or a species population mean."),
        ("DIV-P04", "bounded", "P04 remains the sole retained assembly-derived quantitative result, with raw validation explicitly pending.", "callable heterozygosity", "one P04 individual", "770,780,965 callable bp; 3,548,818 SNPs; pi 0.004604184795871289", "pair_accounting.tsv;vgp_read_validation_per_pair_v1.tsv;paper_pairs.tsv", "Assembly bootstrap concerns PSMC theta rather than raw-read measurement error; exact CLR reads remain unacquired.", "Do not describe P04 as independently read-validated or as a population estimate."),
        ("DIV-P07-FAIL", "supported", "P07 meets the predeclared concrete haplotype-reconstruction-failure rule in both Illumina and HiFi and is preserved but excluded from quantitative synthesis.", "cross-technology assembly-SNP contradiction", "assembly SNP within exact-individual depth-qualified read evidence", "Illumina lower bound 0.501223; HiFi lower bound 0.531662", "vgp_read_validation_results_v1.json;vgp_read_validation_per_pair_v1.tsv", "Chemistries share the individual and reference and are informative but not fully independent.", "Do not delete the artifact, average it with read estimates, or admit its pi/PSMC as validated biology."),
        ("DIV-P07-BRACKET", "bounded", "For the P07 DP10-80 common mask, concordant and union call sets bracket pi at 0.00035313709 and 0.00089988195 under the observed paired callers.", "paired-call concordance bracket", "shared H1 coordinate within the common depth mask", "255,821,332 bp; 90,340 shared, 132,832 assembly-only, 7,037 read-only SNPs", "vgp_read_validation_mask_sensitivity_v1.tsv;vgp_read_validation_results_v1.json", "These are disagreement brackets, not complete biological confidence limits; caller failures are correlated.", "Do not label the bracket endpoints true pi or convert disagreement bounds into proven error rates."),
        ("DIV-P07-MASK", "supported", "The primary common mask retains 95.677% of inherited callable bases but 38.872% of inherited assembly differences; excluded-site density is 34.8067-fold higher.", "assembly callability sensitivity", "H1 assembly coordinate", "P07 inherited and DP10-80 masks", "vgp_read_validation_results_v1.json;vgp_read_validation_mask_sensitivity_v1.tsv", "Masking, mappability, collapse, and caller effects are entangled.", "Do not assert every excluded assembly difference is false."),
        ("DIV-P07-KMER", "suggestive", "P07 k-mer heterozygosity (0.00266527) near the inherited assembly pi but far above stringent mapped-read pi suggests strong representation/callability dependence.", "model-based k-mer heterozygosity sensitivity", "distinct read k-mer spectrum from one individual", "P07 four-component negative-binomial spectrum fit", "vgp_read_validation_results_v1.json;paper_pairs.tsv", "Repeats, coverage bias, correlated k-mers, and shared molecular material are systematic uncertainty.", "Do not select or average the k-mer value as truth or treat it as a population estimate."),
        ("PSMC-COMPUTED", "supported", "Both completed pairs have 64-interval unscaled PSMC trajectories and 200/200 finite block bootstraps centered on primary theta.", "assembly-derived PSMC trajectory", "block of the same pair-specific primary PSMCFA", "P04 and P07", "psmc_unscaled_trajectories.tsv;pair_accounting.tsv;psmc_histories.tsv", "Masked and callable bins remain in the block-bootstrap population; P07 is later invalidated by reads.", "Do not count pi and PSMC as independent evidence or retain P07 demographic biology."),
        ("PSMC-ABSOLUTE", "bounded", "Absolute time and effective-size values are bounded only by nine generic mutation-rate by generation-time sensitivity scenarios.", "scenario-scaled PSMC time and effective size", "pair-specific PSMC interval under generic scenario", "9 scenarios x 64 intervals x 2 pairs", "scenario_uncertainty.tsv;psmc_histories.tsv", "No scenario is species calibrated or preferred.", "Do not quote one absolute history as a species estimate."),
        ("PSMC-P04-SHAPE", "suggestive", "The retained P04 unscaled trajectory has an early lambda maximum (2.874536) and a later broad low (minimum 0.492194), suggesting temporal structure for follow-up.", "relative PSMC shape", "one P04 assembly-derived individual trajectory", "64 unscaled intervals", "psmc_unscaled_trajectories.tsv;psmc_histories.tsv", "Single-individual PSMC resolution, assembly sensitivity, and generic scaling limit interpretation.", "Do not name demographic events, dates, or population sizes from this shape."),
        ("ANNOT-EXACT", "supported", "Six P07 CDS/fourfold/WS/SW partitions use an exact-native annotation with an equal sequence dictionary.", "partitioned assembly heterozygosity", "callable annotated base in one P07 individual", "6 partitions; annotation job 1781559", "exact_native_annotation_partitions.tsv;annotation_partitions.tsv", "The parent P07 quantitative result is invalidated by exact reads.", "Do not use these partitions as validated diversity or gene-conversion evidence."),
        ("ANNOT-ASYM", "suggestive", "P07 fourfold W/S and WS/SW partition differences are an audit signal only.", "directional partition contrast", "annotated fourfold site in one assembly-derived pair", "P07 exact-native annotation", "exact_native_annotation_partitions.tsv;annotation_partitions.tsv", "Mutation spectrum, polarization, callability, reconstruction failure, and non-independence are unresolved.", "Do not call this gBGC, a conversion rate, or transmission bias."),
        ("POPULATION", "unidentifiable", "Species means, within-species distributions, contemporary Ne, and cross-species population relationships are not identifiable from one H1/H2 individual per completed pair.", "population diversity and effective size", "independent individual nested in a sampled population and species", "no conforming population sample", "paper_pairs.tsv;claim_ledger.tsv", "No population replication or sampling frame exists.", "Do not generalize individual assembly estimates to populations or species."),
        ("VERTEBRATE-RANGE", "unidentifiable", "A validated vertebrate diversity range is not identifiable: of two completed assembly observations, P07 is invalidated and P04 is raw-validation pending.", "validated cross-vertebrate diversity distribution", "independent validated species estimate", "0 fully raw-validated admitted quantitative pairs", "paper_pairs.tsv;vgp_read_validation_results_v1.json", "Execution attrition and validation are nonrandom; missing estimates are not zeroes.", "Do not report the 2.144254-fold assembly range as the vertebrate biological range."),
        ("LR-IMPLICATION", "unidentifiable", "This pilot exposes assembly and callability sensitivity that a Lewontin-paradox analysis must control, but it does not test compression of diversity across vertebrates or its relationship to census size.", "cross-species diversity range and diversity-census-size association", "independent population-level species estimate", "two assembly observations; zero fully validated admitted observations", "paper_pairs.tsv;claim_ledger.tsv;figure_diversity.svg", "Tiny, selected, non-population sample with one invalidated pair and one validation-pending pair.", "Do not claim support for, resolution of, or refutation of Lewontin's paradox."),
        ("GENE-CONVERSION", "unidentifiable", "No actual conforming VGP estimate exists in any of the four gene-conversion branches, which remain separate.", "direct, population, historical, and non-allelic gene conversion", "branch-specific conforming sampling unit", "four explicitly separated branches", "gene_conversion_branches.tsv;vgp_comprehensive_table_gene_conversion.tsv", "H1/H2 WS/SW annotation counts do not satisfy any branch estimand.", "Do not merge branches or treat annotation partitions as conversion estimates."),
    ]
    return [dict(zip(CLAIM_FIELDS, (root,) + row)) for row in rows]


def build_job_rows(data: Mapping[str, object]) -> list[dict[str, str]]:
    output: list[dict[str, str]] = []
    for packet, rows in (("scale_vgp_real", data["scale_jobs"]),
                         ("validate_vgp_pilot_reads", data["read_jobs"])):
        for source_row, row in enumerate(rows, 1):
            state = row["State"]
            output.append({
                "canonical_vgp_root": str(data["canonical_root"]),
                "packet": packet,
                "source_row": str(source_row),
                "job_id": row["JobIDRaw"],
                "job_name": row["JobName"],
                "partition": row.get("Partition", ""),
                "state": state,
                "elapsed": row["Elapsed"],
                "elapsed_seconds": row["ElapsedRaw"],
                "timelimit": row["Timelimit"],
                "allocated_cpus": row["AllocCPUS"],
                "requested_memory": row["ReqMem"],
                "max_rss": row["MaxRSS"],
                "cpu_time_raw": row["CPUTimeRAW"],
                "total_cpu": row["TotalCPU"],
                "exit_code": row["ExitCode"],
                "node_list": row["NodeList"],
                "start": row["Start"],
                "end": row["End"],
                "submit": row["Submit"],
                "failure_or_nonterminal": str(state != "COMPLETED").lower(),
            })
    return output


def _walk_json_digests(value: object, pointer: str = "") -> Iterable[tuple[str, str, str]]:
    if isinstance(value, dict):
        sibling_path = next((str(value[key]) for key in ("path", "source_path", "view", "local_path")
                             if isinstance(value.get(key), str) and value[key]), "")
        for key, child in value.items():
            child_pointer = f"{pointer}/{key}"
            if "sha256" in key.lower() and isinstance(child, str) and SHA256_RE.fullmatch(child):
                yield child_pointer, sibling_path or child_pointer, child
            yield from _walk_json_digests(child, child_pointer)
    elif isinstance(value, list):
        for index, child in enumerate(value):
            yield from _walk_json_digests(child, f"{pointer}/{index}")


def _walk_json_md5(value: object, pointer: str = "") -> Iterable[tuple[str, str, str]]:
    if isinstance(value, dict):
        sibling_path = next((str(value[key]) for key in ("path", "source_path", "view", "local_path", "source_url")
                             if isinstance(value.get(key), str) and value[key]), "")
        for key, child in value.items():
            child_pointer = f"{pointer}/{key}"
            explicitly_md5 = key.lower().endswith("md5")
            declared_md5 = key == "upstream_checksum" and value.get("upstream_checksum_algorithm") == "md5"
            if (explicitly_md5 or declared_md5) and isinstance(child, str) and MD5_RE.fullmatch(child):
                yield child_pointer, sibling_path or child_pointer, child
            yield from _walk_json_md5(child, child_pointer)
    elif isinstance(value, list):
        for index, child in enumerate(value):
            yield from _walk_json_md5(child, f"{pointer}/{index}")


def build_digest_rows(canonical_root: str) -> list[dict[str, str]]:
    output: list[dict[str, str]] = []
    for path in INPUT_PATHS:
        output.append({
            "canonical_vgp_root": canonical_root,
            "binding_class": "REPOSITORY_INPUT",
            "container_path": relative(path),
            "locator": relative(path),
            "digest_field": "file_sha256",
            "algorithm": "sha256",
            "digest": sha256_file(path),
            "occurrences": "1",
            "bytes": str(path.stat().st_size),
            "verification": "PASS_REHASHED",
        })

    grouped: Counter[tuple[str, str, str, str, str]] = Counter()
    for path in INPUT_PATHS:
        if path.suffix == ".json":
            for pointer, locator, digest in _walk_json_digests(read_json(path)):
                grouped[(relative(path), locator, pointer.rsplit("/", 1)[-1], "sha256", digest)] += 1
            value = read_json(path)
            for pointer, locator, digest in _walk_json_md5(value):
                grouped[(relative(path), locator, pointer.rsplit("/", 1)[-1], "md5", digest)] += 1
        elif path.suffix == ".tsv":
            rows = read_tsv(path)
            for row_number, row in enumerate(rows, 2):
                for field, value in row.items():
                    algorithm = "sha256" if "sha256" in field.lower() and SHA256_RE.fullmatch(value or "") else ""
                    if not algorithm and (field.lower().endswith("md5") or field == "upstream_checksum") and MD5_RE.fullmatch(value or ""):
                        algorithm = "md5"
                    if not algorithm:
                        continue
                    path_field = field[:-6] + "path" if field.endswith("sha256") else ""
                    locator = row.get(path_field, "") or row.get("source_path", "")
                    if not locator:
                        locator = f"{relative(path)}#{field}"
                    grouped[(relative(path), locator, field, algorithm, value)] += 1
    for (container, locator, field, algorithm, digest), occurrences in sorted(grouped.items()):
        output.append({
            "canonical_vgp_root": canonical_root,
            "binding_class": "UPSTREAM_EMBEDDED_DIGEST",
            "container_path": container,
            "locator": locator,
            "digest_field": field,
            "algorithm": algorithm,
            "digest": digest,
            "occurrences": str(occurrences),
            "bytes": "",
            "verification": "PASS_INHERITED_UPSTREAM_PACKET",
        })
    return output


def diversity_svg(data: Mapping[str, object]) -> str:
    root = escape(str(data["canonical_root"]))
    points = [
        ("P04 assembly pi (retained; raw pending)", 0.004604184795871289, "#176b46", "circle"),
        ("P07 inherited assembly pi (invalidated)", 0.0021472198306856562, "#b33b32", "cross"),
        ("P07 common-mask assembly pi", 0.0008723744742287559, "#c46a27", "cross"),
        ("P07 common-mask read pi", 0.00038064456641950404, "#365f91", "cross"),
        ("P07 k-mer model heterozygosity", 0.00266526705636827, "#78569b", "cross"),
    ]
    left, right, top = 255.0, 855.0, 112.0
    rows = []
    for index, (label, value, color, marker) in enumerate(points):
        y = top + index * 56
        x = left + value / 0.005 * (right - left)
        mark = (f'<circle cx="{x:.2f}" cy="{y}" r="7" fill="{color}"/>' if marker == "circle" else
                f'<path d="M{x-7:.2f},{y-7} L{x+7:.2f},{y+7} M{x-7:.2f},{y+7} L{x+7:.2f},{y-7}" stroke="{color}" stroke-width="3"/>')
        rows.append(f'<text x="18" y="{y+5}" font-size="15">{escape(label)}</text>{mark}'
                    f'<text x="{x+11:.2f}" y="{y+5}" font-size="13" fill="{color}">{value:.9f}</text>')
    ticks = []
    for value in (0, 0.001, 0.002, 0.003, 0.004, 0.005):
        x = left + value / 0.005 * (right - left)
        ticks.append(f'<line x1="{x:.2f}" y1="85" x2="{x:.2f}" y2="388" stroke="#dedede"/>'
                     f'<text x="{x:.2f}" y="414" text-anchor="middle" font-size="13">{value:.3f}</text>')
    return f'''<svg xmlns="http://www.w3.org/2000/svg" width="920" height="485" viewBox="0 0 920 485">
<title>Observed VGP individual diversity and P07 method sensitivity</title>
<desc>Canonical VGP root {root}. P04 is retained with raw validation pending. Every P07 point is crossed because exact reads trigger a concrete reconstruction failure. Values are per-base individual estimates, not population means.</desc>
<rect width="920" height="485" fill="white"/>
<text x="460" y="30" text-anchor="middle" font-family="sans-serif" font-size="21" font-weight="bold">VGP diversity: measured values and validation status</text>
<text x="460" y="55" text-anchor="middle" font-family="sans-serif" font-size="13">Two inherited assembly estimates span 2.144254×; only P04 is retained, with raw validation pending</text>
<g font-family="sans-serif">{''.join(ticks)}<line x1="{left}" y1="388" x2="{right}" y2="388" stroke="#222" stroke-width="1.5"/>{''.join(rows)}
<text x="555" y="442" text-anchor="middle" font-size="15">heterozygosity / callable base</text>
<text x="18" y="468" font-size="12" fill="#555">Circle = retained assembly observation; crosses = P07 sensitivity values excluded from quantitative synthesis.</text></g>
</svg>\n'''


def psmc_svg(data: Mapping[str, object]) -> str:
    root = escape(str(data["canonical_root"]))
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in data["trajectory_rows"]:
        grouped[row["selection_id"]].append(row)
    left, right, top, bottom = 92.0, 910.0, 78.0, 425.0
    max_x = math.log10(1 + 30.0)
    min_y, max_y = -0.4, 2.3

    def coordinates(rows: Sequence[Mapping[str, str]]) -> str:
        points = []
        for row in sorted(rows, key=lambda item: int(item["interval"])):
            xvalue = math.log10(1 + float(row["time_2N0"]))
            yvalue = math.log10(float(row["lambda"]))
            x = left + xvalue / max_x * (right - left)
            y = bottom - (yvalue - min_y) / (max_y - min_y) * (bottom - top)
            points.append(f"{x:.2f},{y:.2f}")
        return " ".join(points)

    x_ticks = []
    for value in (0, 0.1, 0.3, 1, 3, 10, 30):
        x = left + math.log10(1 + value) / max_x * (right - left)
        x_ticks.append(f'<line x1="{x:.2f}" y1="{top}" x2="{x:.2f}" y2="{bottom}" stroke="#e4e4e4"/>'
                       f'<text x="{x:.2f}" y="448" text-anchor="middle" font-size="12">{value:g}</text>')
    y_ticks = []
    for value in (1, 3, 10, 30, 100):
        y = bottom - (math.log10(value) - min_y) / (max_y - min_y) * (bottom - top)
        y_ticks.append(f'<line x1="{left}" y1="{y:.2f}" x2="{right}" y2="{y:.2f}" stroke="#e4e4e4"/>'
                       f'<text x="78" y="{y+4:.2f}" text-anchor="end" font-size="12">{value:g}</text>')
    return f'''<svg xmlns="http://www.w3.org/2000/svg" width="980" height="535" viewBox="0 0 980 535">
<title>Assembly-derived unscaled PSMC histories for P04 and P07</title>
<desc>Canonical VGP root {root}. P04 is a retained descriptive unscaled trajectory with raw validation pending. P07 is dashed and invalidated by exact-read reconstruction failure. Both PSMC trajectories are non-independent of same-pair pi.</desc>
<rect width="980" height="535" fill="white"/>
<text x="490" y="29" text-anchor="middle" font-family="sans-serif" font-size="21" font-weight="bold">Assembly-derived PSMC histories (unscaled)</text>
<text x="490" y="53" text-anchor="middle" font-family="sans-serif" font-size="13">Relative time and lambda only; nine absolute scenarios are generic sensitivities, not species calibrations</text>
<g font-family="sans-serif">{''.join(x_ticks)}{''.join(y_ticks)}
<line x1="{left}" y1="{bottom}" x2="{right}" y2="{bottom}" stroke="#222"/><line x1="{left}" y1="{top}" x2="{left}" y2="{bottom}" stroke="#222"/>
<polyline points="{coordinates(grouped['P04'])}" fill="none" stroke="#176b46" stroke-width="3"/>
<polyline points="{coordinates(grouped['P07'])}" fill="none" stroke="#b33b32" stroke-width="3" stroke-dasharray="8 5"/>
<rect x="638" y="75" width="280" height="54" rx="4" fill="white" fill-opacity="0.9"/>
<line x1="650" y1="91" x2="690" y2="91" stroke="#176b46" stroke-width="3"/><text x="700" y="96" font-size="13">P04 retained / raw pending</text>
<line x1="650" y1="115" x2="690" y2="115" stroke="#b33b32" stroke-width="3" stroke-dasharray="8 5"/><text x="700" y="120" font-size="13">P07 invalidated artifact</text>
<text x="500" y="478" text-anchor="middle" font-size="14">time / 2N0 (log1p axis)</text><text x="22" y="260" transform="rotate(-90 22 260)" text-anchor="middle" font-size="14">relative effective size lambda (log axis)</text>
<text x="25" y="510" font-size="12" fill="#555">PSMC and pi reuse each H1/H2 pair: the two outcomes are descriptive and non-independent.</text></g>
</svg>\n'''


def build_report(data: Mapping[str, object], pair_rows: Sequence[Mapping[str, str]],
                 psmc_rows: Sequence[Mapping[str, str]], claim_rows: Sequence[Mapping[str, str]]) -> str:
    reads = data["reads"]
    p07 = reads["pairs"]["P07"]
    primary = p07["primary_common_mask_comparison"]
    evidence = p07["read_backed_assembly_sites"]
    p04_hist = next(row for row in psmc_rows if row["selection_id"] == "P04")
    p07_hist = next(row for row in psmc_rows if row["selection_id"] == "P07")
    classes = Counter(row["classification"] for row in claim_rows)
    return f"""# Real VGP biological evidence synthesis

Generated from the closed-world upstream packets completed {reads['completed_utc']}
WG task: `synthesize-vgp-real`
Canonical shared VGP root: `{data['canonical_root']}`

## Paper-oriented result

The real VGP execution produced two assembly-derived individual diversity estimates, not zero biological estimates. P07 (*Spinachia spinachia*) has pi = 0.0021472198306856562 from 574,122 differences across 267,379,237 callable bp; P04 (*Falco naumanni*) has pi = 0.004604184795871289 from 3,548,818 differences across 770,780,965 callable bp. Their observed assembly-derived range is therefore 0.002147–0.004604 per bp, or 2.144254-fold.

That computational range is not a validated vertebrate biological range. Exact-individual Illumina and HiFi validation makes P07 a `concrete_haplotype_reconstruction_failure`, so its pi and PSMC are preserved for provenance but excluded from quantitative biological use. P04 remains the sole retained assembly-derived estimate, with its exact CLR validation run pending. Thus there are two completed assembly observations, one retained/raw-pending observation, one read-invalidated observation, and zero fully raw-validated admitted quantitative pairs.

## Execution and confidence accounting

The immutable Freeze 1 accounting covers all 716 catalog rows and all 569 links. Ten links belong to the frozen same-individual, mutually comparable audit roster. Across the pilot and scale packet, all 650 scheduler allocations are retained: 240 completed, 42 failed, 346 cancelled, 21 pending, and one running. The raw-read packet adds eight distinct allocations: two completed, four failed, and two cancelled. Nothing missing or nonterminal is converted to a biological zero.

| Synthesis tier | Pair(s) | Meaning |
|---|---|---|
| T1 retained/raw pending | P04 | Completed assembly pi and PSMC retained; exact raw validation pending |
| T2 invalidated by exact reads | P07 | Artifact preserved, but pi, PSMC, and annotation partitions excluded from quantitative biological synthesis |
| X hard-invalid primary | P01, P02, P03, P05 | Invalid primary execution; no estimate |
| X execution error | P06, P09, P10 | Concrete execution error; no estimate |
| P resumable at freeze | P08 | Running wave at cutoff; no frozen estimate |

The machine-readable pair table is `paper_pairs.tsv`; the complete 658-record scheduler reconciliation is `job_ledger.tsv`.

## Raw-read validation and assembly sensitivity

At the P07 primary DP10-80 common mask, {primary['callable_bp_common_mask']:,} bp support assembly pi = {primary['assembly_pi_common_mask']:.12g} and mapped-read pi = {primary['read_pi_common_mask']:.12g}, a read/assembly ratio of {primary['pi_ratio_read_over_assembly']:.6f}. Shared calls yield a lower concordant bracket of {primary['concordant_pi_lower_bracket']:.12g}; the union yields an upper bracket of {primary['union_pi_upper_bracket']:.12g}. These are paired-caller disagreement bounds, not complete biological confidence intervals.

The primary mask retains {p07['primary_mask_exclusion_diagnostic']['primary_callable_bp_fraction_of_core']:.3%} of inherited callable bases but {p07['primary_mask_exclusion_diagnostic']['primary_assembly_site_fraction_of_core']:.3%} of inherited assembly differences. The excluded assembly-difference density is {p07['primary_mask_exclusion_diagnostic']['excluded_to_primary_assembly_site_density_ratio']:.4f}-fold the retained density, a strong mappability/collapse/callability sensitivity signal.

More decisively, depth-qualified homozygous-reference contradictions are the majority in both technologies: Illumina {evidence['illumina']['concrete_false_positive_lower_bound_fraction']:.3%} (Wilson 95% {evidence['illumina']['concrete_false_positive_fraction_wilson_95'][0]:.3%}–{evidence['illumina']['concrete_false_positive_fraction_wilson_95'][1]:.3%}) and HiFi {evidence['hifi']['concrete_false_positive_lower_bound_fraction']:.3%} (Wilson 95% {evidence['hifi']['concrete_false_positive_fraction_wilson_95'][0]:.3%}–{evidence['hifi']['concrete_false_positive_fraction_wilson_95'][1]:.3%}). These satisfy the predeclared reconstruction-failure rule. The k-mer QV of {p07['kmer_qv']['qv']:.4f} does not meet the separate severe-sequence-error rule (QV <20), which shows why consensus quality and haplotype reconstruction must remain distinct.

P07's negative-binomial k-mer spectrum gives model-based heterozygosity {p07['kmer_heterozygosity']['heterozygosity_per_base']:.10f}, {p07['kmer_heterozygosity_to_core_pi_ratio']:.5f} times inherited assembly pi but about {p07['kmer_heterozygosity']['heterozygosity_per_base']/primary['read_pi_common_mask']:.3f} times stringent common-mask read pi. This disagreement is reported as method sensitivity; estimates are not averaged.

P09's one-cell 0.810x HiFi run is only a mapping compatibility control. P04's 42,344,746,693-byte exact CLR run remains planned. Four canonical raw objects totaling 31,058,137,613 bytes were independently rehashed and reused. One full-size corrupted P07 R2 transfer (4,431,902,981 bytes; SHA-256 `defaed9e929d8acf9d58006a3be51c26dd7c1937079f7d47e9d818420475c965`) remains quarantined; a clean retry produced canonical SHA-256 `c542f6efd9fc1d8f557c89629743fbe4a39584f24002c84185e479057cb443ac` without redownload of verified objects.

## Demographic histories

Each completed pair has a 64-interval assembly-derived PSMC trajectory and 200/200 finite block bootstraps. P04 theta0 is {p04_hist['primary_theta_0_per_100bp_bin']} per 100-bp bin (bootstrap {p04_hist['bootstrap_theta_q025']}–{p04_hist['bootstrap_theta_q975']}); its unscaled lambda rises to {p04_hist['trajectory_lambda_max']} at interval {p04_hist['trajectory_lambda_max_interval']} and falls to a broad later minimum of {p04_hist['trajectory_lambda_min']} at interval {p04_hist['trajectory_lambda_min_interval']}. This is a descriptive relative-history shape only.

P07 theta0 is {p07_hist['primary_theta_0_per_100bp_bin']} (bootstrap {p07_hist['bootstrap_theta_q025']}–{p07_hist['bootstrap_theta_q975']}) and its inherited trajectory reaches lambda {p07_hist['trajectory_lambda_max']}; exact-read PSMC has a similar theta ratio ({p07['primary_psmc_comparison']['theta_ratio_read_over_assembly']:.6f}) but strongly discordant lambda shape (Pearson r = {p07['primary_psmc_comparison']['lambda_pearson_correlation']:.6f}). The time grids and data sources overlap, so this is sensitivity, not replication, and the inherited P07 history is invalidated for biology.

For both pairs, absolute histories exist only as nine generic mutation-rate × generation-time scenarios. Across those grids, P04 spans {p04_hist['generic_effective_size_min']}–{p04_hist['generic_effective_size_max']} in scenario Ne and 0–{p04_hist['generic_time_years_max']} scenario years; P07 spans {p07_hist['generic_effective_size_min']}–{p07_hist['generic_effective_size_max']} and 0–{p07_hist['generic_time_years_max']}. No scenario is preferred or species calibrated. Crucially, PSMC and pi use the same H1/H2 pair and are explicitly non-independent.

## Annotation and gene-conversion separation

P07 has an exact-native annotation with equal assembly/annotation sequence dictionaries and six measured partitions: CDS, fourfold, fourfold-W, fourfold-S, WS, and SW. The fourfold estimate is 0.000847594, with W and S estimates 0.001166236 and 0.000680551; WS and SW normalized counts are 0.000925017 and 0.000474837. These are real assembly-derived annotation results, but the parent pair is read-invalidated. Mutation spectrum, polarization, and reconstruction error are unresolved, so the asymmetry is suggestive only.

The direct-pedigree/gamete, population AFS, historical phylogenetic, and non-allelic paralog branches remain separate in `gene_conversion_branches.tsv`. No actual conforming VGP estimate exists in any branch. H1/H2 WS/SW counts are not transmission events, population B, historical substitutions, or copy-resolved paralog tracts.

## Implication for Lewontin's paradox

The measured assembly range demonstrates that the program now contains real nonzero biological estimates. Its stronger lesson is methodological: inferred diversity can change sharply with representation, common-mask definition, and read-backed haplotype validation. Any test of Lewontin's paradox must therefore separate biological between-species variation from assembly/callability variation and must sample populations rather than single assembly individuals.

This pilot cannot itself test the paradox. Two selected assembly observations—one invalidated, one validation-pending—cannot identify a cross-vertebrate diversity distribution, a relationship with census size, or a selection/drift explanation. Reporting 2.144254-fold as the vertebrate diversity range would ignore the exact-read adjudication and nonrandom execution attrition.

## Claim classification and reproducibility

The claim ledger contains {classes['supported']} supported, {classes['bounded']} bounded, {classes['suggestive']} suggestive, and {classes['unidentifiable']} unidentifiable claims. Every claim states its sampling unit, covariance, and forbidden inference. `digest_ledger.tsv` rehashes every repository input and inventories every embedded upstream SHA-256 binding; `manifest.json` binds all emitted outputs.

Both compute environments use GNU Guix channel `{data['environment']['frozen_channel_commit']}`. The core pi/PSMC workflow profile is `{data['pilot_environment']['profile']}` (closure SHA-256 `{data['pilot_environment']['closure_sha256']}`); raw-read validation uses `{data['environment']['profile']}` (closure SHA-256 `{data['environment']['closure_sha256']}`). Both completed FastGA mappings passed live node-local `/scratch` contracts and reproduced the frozen PAFs. Canonical data, CAS objects, views, raw reads, run outputs, and promoted products resolve from the single configured VGP root shown above; the legacy Lewontin-paradox path remains migration input only.
"""


def generate(output_dir: Path = OUTPUT_DIR) -> dict[str, object]:
    data = reconcile_inputs()
    pair_rows = build_pair_rows(data)
    psmc_rows = build_psmc_rows(data)
    annotation_rows = build_annotation_rows(data)
    gc_rows = build_gene_conversion_rows(data)
    claim_rows = build_claim_rows(data)
    job_rows = build_job_rows(data)
    digest_rows = build_digest_rows(str(data["canonical_root"]))

    paths = {path.name: output_dir / path.name for path in OUTPUT_PATHS}
    atomic_tsv(paths[PAPER_PAIRS.name], PAPER_FIELDS, pair_rows)
    atomic_tsv(paths[PSMC_HISTORIES.name], PSMC_FIELDS, psmc_rows)
    atomic_tsv(paths[ANNOTATION_TABLE.name], ANNOTATION_FIELDS, annotation_rows)
    atomic_tsv(paths[GENE_CONVERSION.name], GENE_CONVERSION_FIELDS, gc_rows)
    atomic_tsv(paths[CLAIMS.name], CLAIM_FIELDS, claim_rows)
    atomic_tsv(paths[JOB_LEDGER.name], JOB_FIELDS, job_rows)
    atomic_tsv(paths[DIGEST_LEDGER.name], DIGEST_FIELDS, digest_rows)
    atomic_text(paths[DIVERSITY_FIGURE.name], diversity_svg(data))
    atomic_text(paths[PSMC_FIGURE.name], psmc_svg(data))
    atomic_text(paths[REPORT.name], build_report(data, pair_rows, psmc_rows, claim_rows))

    pis = sorted(float(row["callable_pi"]) for row in data["measured"])
    read_acquisition = data["reads"]["raw_read_acquisition"]
    manifest: dict[str, object] = {
        "schema_version": "vgp-real-biological-synthesis-v1",
        "task_id": "synthesize-vgp-real",
        "generated_at_utc": data["reads"]["completed_utc"],
        "canonical_vgp_root": data["canonical_root"],
        "root_config": relative(ROOT_CONFIG),
        "legacy_root_role": "migration_input_only",
        "closed_world": {
            "catalog_rows": 716,
            "catalog_links": 569,
            "audited_pairs": 10,
            "verified_core_results": 2,
            "hard_invalid_primary_pairs": 4,
            "execution_error_pairs": 3,
            "resumable_running_pairs_at_freeze": 1,
            "annotation_partition_rows": 6,
            "psmc_trajectory_rows": 128,
            "generic_scenario_rows": 1152,
        },
        "biological_evidence": {
            "assembly_derived_diversity": {
                "completed_pair_count": 2,
                "selection_ids": ["P04", "P07"],
                "minimum_pi": pis[0],
                "maximum_pi": pis[1],
                "max_to_min_ratio": pis[1] / pis[0],
                "validated_quantitative_pair_count": 0,
                "retained_raw_pending_pair_count": 1,
                "invalidated_by_raw_reads_pair_count": 1,
            },
            "psmc_completed_pairs": 2,
            "exact_native_annotation_pairs": 1,
            "exact_native_annotation_partitions": 6,
        },
        "jobs": {
            "scale_vgp_real": {"allocations": len(data["scale_jobs"]), "states": _state_counts(data["scale_jobs"])},
            "validate_vgp_pilot_reads": {"allocations": len(data["read_jobs"]), "states": _state_counts(data["read_jobs"])},
            "combined_distinct_scheduler_records": len(job_rows),
            "combined_failure_cancelled_or_nonterminal_records": sum(row["failure_or_nonterminal"] == "true" for row in job_rows),
        },
        "raw_reads": {
            "verified_or_reused_objects": read_acquisition["summary"]["verified"]["objects"],
            "verified_or_reused_bytes": read_acquisition["summary"]["verified"]["bytes"],
            "planned_pending_objects": read_acquisition["summary"]["pending"]["objects"],
            "planned_pending_bytes": read_acquisition["summary"]["pending"]["bytes"],
            "quarantined_failed_transfer_count": 1,
            "quarantined_failed_transfer_sha256": read_acquisition["retained_failed_R2_transfer"]["sha256"],
        },
        "environment": {
            "guix_channel_commit": data["environment"]["frozen_channel_commit"],
            "core_guix_profile": data["pilot_environment"]["profile"],
            "core_guix_closure_sha256": data["pilot_environment"]["closure_sha256"],
            "read_validation_guix_profile": data["environment"]["profile"],
            "read_validation_guix_closure_sha256": data["environment"]["closure_sha256"],
        },
        "analysis_contract": {
            "same_pair_pi_psmc_independent": False,
            "species_calibrated_absolute_psmc": False,
            "population_inference_authorized": False,
            "annotation_absence_core_veto": False,
            "gene_conversion_branches_integrated": False,
            "missing_estimates_encoded_as_zero": False,
        },
        "claim_classifications": dict(sorted(Counter(row["classification"] for row in claim_rows).items())),
        "digest_accounting": {
            "repository_inputs_rehashed": sum(row["binding_class"] == "REPOSITORY_INPUT" for row in digest_rows),
            "unique_upstream_embedded_bindings": sum(row["binding_class"] == "UPSTREAM_EMBEDDED_DIGEST" for row in digest_rows),
            "upstream_embedded_occurrences": sum(int(row["occurrences"]) for row in digest_rows if row["binding_class"] == "UPSTREAM_EMBEDDED_DIGEST"),
            "algorithms": dict(sorted(Counter(row["algorithm"] for row in digest_rows).items())),
        },
        "input_digests": {relative(path): sha256_file(path) for path in INPUT_PATHS},
    }
    manifest["output_digests"] = {
        relative(OUTPUT_DIR / path.name): sha256_file(output_dir / path.name) for path in OUTPUT_PATHS
    }
    atomic_json(output_dir / MANIFEST.name, manifest)
    errors = validate_outputs(output_dir)
    require(not errors, "generated synthesis failed validation: " + "; ".join(errors))
    return manifest


def validate_outputs(output_dir: Path = OUTPUT_DIR, verify_digests: bool = True) -> list[str]:
    errors: list[str] = []
    try:
        manifest = read_json(output_dir / MANIFEST.name)
        canonical_root = str(manifest["canonical_vgp_root"])
        pairs = read_tsv(output_dir / PAPER_PAIRS.name)
        histories = read_tsv(output_dir / PSMC_HISTORIES.name)
        annotations = read_tsv(output_dir / ANNOTATION_TABLE.name)
        branches = read_tsv(output_dir / GENE_CONVERSION.name)
        claims = read_tsv(output_dir / CLAIMS.name)
        jobs = read_tsv(output_dir / JOB_LEDGER.name)
        digests = read_tsv(output_dir / DIGEST_LEDGER.name)
        for name, rows in (("pairs", pairs), ("histories", histories), ("annotations", annotations),
                           ("branches", branches), ("claims", claims), ("jobs", jobs), ("digests", digests)):
            if not rows or any(row.get("canonical_vgp_root") != canonical_root for row in rows):
                errors.append(f"{name} canonical root mismatch")
        by_pair = {row["selection_id"]: row for row in pairs}
        if len(pairs) != 10 or set(by_pair) != {f"P{i:02d}" for i in range(1, 11)}:
            errors.append("paper pair roster mismatch")
        if by_pair.get("P07", {}).get("quantitative_disposition") != "PRESERVED_NOT_ADMITTED_QUANTITATIVELY":
            errors.append("P07 invalidated quantitative disposition changed")
        if by_pair.get("P04", {}).get("quantitative_disposition") != "RETAINED_ASSEMBLY_DERIVED_PENDING_RAW_VALIDATION":
            errors.append("P04 retained/raw-pending disposition changed")
        if any(row["pi_psmc_independence"] != "SAME_PAIR_NONINDEPENDENT" for row in histories):
            errors.append("PSMC was promoted to independent evidence")
        if any(row["gene_conversion_disposition"] != "NOT_A_CONFORMING_GENE_CONVERSION_ESTIMATE" for row in annotations):
            errors.append("annotation was promoted to gene-conversion evidence")
        if len(branches) != 4 or any(row["vgp_integration_status"] != "SEPARATE_NO_ACTUAL_CONFORMING_VGP_ESTIMATE" or
                                     row["estimate"] != "NOT_ESTIMABLE" for row in branches):
            errors.append("gene-conversion branches were merged or promoted")
        if {row["classification"] for row in claims} != {"supported", "bounded", "suggestive", "unidentifiable"}:
            errors.append("claim classification set mismatch")
        if len(jobs) != 658 or len({row["job_id"] for row in jobs}) != 658:
            errors.append("scheduler job reconciliation mismatch")
        if len([row for row in digests if row["binding_class"] == "REPOSITORY_INPUT"]) != len(INPUT_PATHS):
            errors.append("repository input digest accounting mismatch")
        for name in (DIVERSITY_FIGURE.name, PSMC_FIGURE.name):
            text = (output_dir / name).read_text(encoding="utf-8")
            if not text.startswith("<svg") or "<title>" not in text or canonical_root not in text:
                errors.append(f"invalid static SVG: {name}")
        report = (output_dir / REPORT.name).read_text(encoding="utf-8")
        if "zero fully raw-validated admitted quantitative pairs" not in report or "cannot itself test the paradox" not in report:
            errors.append("paper report interpretation guards missing")
        if verify_digests:
            for relative_path, expected in manifest["output_digests"].items():
                path = output_dir / Path(relative_path).name
                if not path.is_file() or sha256_file(path) != expected:
                    errors.append(f"output digest mismatch: {relative_path}")
    except (OSError, KeyError, TypeError, ValueError, json.JSONDecodeError, SynthesisError) as exc:
        errors.append(str(exc))
    return errors


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    parser.add_argument("--validate-only", action="store_true")
    args = parser.parse_args(argv)
    if args.validate_only:
        errors = validate_outputs(args.output_dir)
        if errors:
            for error in errors:
                print(error)
            return 1
        print(f"PASS: validated {args.output_dir}")
        return 0
    manifest = generate(args.output_dir)
    print(json.dumps({
        "canonical_vgp_root": manifest["canonical_vgp_root"],
        "completed_assembly_estimates": manifest["biological_evidence"]["assembly_derived_diversity"]["completed_pair_count"],
        "scheduler_records": manifest["jobs"]["combined_distinct_scheduler_records"],
        "output_dir": str(args.output_dir),
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
