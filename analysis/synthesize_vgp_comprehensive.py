#!/usr/bin/env python3
"""Generate the fail-closed comprehensive VGP evidence synthesis.

This is a reconciliation program, not a biological analysis launcher.  It
reads only committed, reviewed artifacts; it neither acquires sequence nor
submits compute.  In particular, zero admitted observations are represented as
NOT_ESTIMABLE, never as biological zeroes.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import os
from collections import Counter
from pathlib import Path
from typing import Iterable, Mapping, Sequence


ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = ROOT / "analysis"

CORE_MANIFEST = ANALYSIS / "vgp_core_scaleout_manifest.tsv"
CORE_QC = ANALYSIS / "vgp_core_scaleout_qc.tsv"
CORE_TELEMETRY = ANALYSIS / "vgp_core_scaleout_telemetry.tsv"
CORE_SUMMARY = ANALYSIS / "vgp_core_scaleout_summary.json"
CORE_RESULTS = ANALYSIS / "vgp_core_scaleout_results.md"
CORE_WAVES = ANALYSIS / "vgp_core_scaleout_wave_manifest.tsv"
CORE_PAPER_PAIRS = ANALYSIS / "vgp_core_scaleout_paper_pairs.tsv"
CORE_PAPER_SUMMARY = ANALYSIS / "vgp_core_scaleout_paper_summary.tsv"
CORE_SENSITIVITY = ANALYSIS / "vgp_core_scaleout_sensitivity.tsv"
CORE_SCALING = ANALYSIS / "vgp_core_scaleout_scaling_scenarios.tsv"
CORE_INDEPENDENT = ANALYSIS / "vgp_core_scaleout_independent_validation.tsv"
MIRROR_SUMMARY = ANALYSIS / "vgp_freeze1_mirror_summary.json"
MIRROR_MANIFEST = ANALYSIS / "vgp_freeze1_mirror_manifest.tsv"
REVIEW_DECISION = ANALYSIS / "vgp_10_pilot_review_decision.json"
REVIEW_GATES = ANALYSIS / "vgp_10_pilot_review_gates.tsv"
REVIEW_SCALEOUT = ANALYSIS / "vgp_10_pilot_scaleout_manifest.tsv"
REVIEW_RESOURCES = ANALYSIS / "vgp_10_pilot_scaleout_resource_manifest.tsv"
DIRECT_DATASETS = ANALYSIS / "direct_gene_conversion_dataset_manifest.tsv"
DIRECT_TRACTS = ANALYSIS / "direct_gene_conversion_tracts.tsv"
DIRECT_SUMMARY = ANALYSIS / "direct_gene_conversion_summary.tsv"
DIRECT_REPORT = ANALYSIS / "direct_gene_conversion_pilot.md"
PHYLO_CLADES = ANALYSIS / "vgp_phylogenetic_gbgc_clade_manifest.tsv"
PHYLO_QC = ANALYSIS / "vgp_phylogenetic_gbgc_qc.tsv"
PHYLO_RESULTS = ANALYSIS / "vgp_phylogenetic_gbgc_results.tsv"
PHYLO_REPORT = ANALYSIS / "vgp_phylogenetic_gbgc_pilot.md"
GC_DATASET_DESIGN = ANALYSIS / "gene_conversion_dataset_manifest.tsv"
GC_ESTIMAND_DESIGN = ANALYSIS / "gene_conversion_estimand_manifest.tsv"
GC_CLAIM_DESIGN = ANALYSIS / "gene_conversion_claim_matrix.tsv"
GUIX_CHANNELS = ANALYSIS / "guix/channels.scm"
GUIX_MANIFEST = ANALYSIS / "guix/manifest.scm"
PHYLO_GUIX_MANIFEST = ANALYSIS / "guix/vgp-phylogenetic-gbgc-preflight-manifest.scm"
VGP_DESIGN = ANALYSIS / "vgp_analysis_manifest.json"
METHOD = ANALYSIS / "synthesize_vgp_comprehensive.py"
GUIX_RUNNER = ANALYSIS / "run_vgp_comprehensive_synthesis_guix.sh"
METHOD_TEST = ANALYSIS / "tests/test_vgp_comprehensive_synthesis.py"

SYNTHESIS = ANALYSIS / "vgp_comprehensive_synthesis.md"
CLAIM_LEDGER = ANALYSIS / "vgp_comprehensive_claim_ledger.tsv"
FINAL_MANIFEST = ANALYSIS / "vgp_comprehensive_final_manifest_v2.json"
CORE_TABLE = ANALYSIS / "vgp_comprehensive_table_core.tsv"
GENE_CONVERSION_TABLE = ANALYSIS / "vgp_comprehensive_table_gene_conversion.tsv"
CLOSED_WORLD_FIGURE = ANALYSIS / "vgp_comprehensive_figure_closed_world.svg"
EVIDENCE_FIGURE = ANALYSIS / "vgp_comprehensive_figure_evidence.svg"

OUTPUT_PATHS = (
    SYNTHESIS, CLAIM_LEDGER, CORE_TABLE, GENE_CONVERSION_TABLE,
    CLOSED_WORLD_FIGURE, EVIDENCE_FIGURE,
)
OUTPUT_FILENAMES = tuple(path.name for path in OUTPUT_PATHS) + (FINAL_MANIFEST.name,)

INPUT_PATHS = (
    CORE_MANIFEST, CORE_QC, CORE_TELEMETRY, CORE_SUMMARY, CORE_RESULTS,
    CORE_WAVES, CORE_PAPER_PAIRS, CORE_PAPER_SUMMARY, CORE_SENSITIVITY,
    CORE_SCALING, CORE_INDEPENDENT, MIRROR_SUMMARY, MIRROR_MANIFEST,
    REVIEW_DECISION, REVIEW_GATES, REVIEW_SCALEOUT, REVIEW_RESOURCES,
    DIRECT_DATASETS, DIRECT_TRACTS, DIRECT_SUMMARY, DIRECT_REPORT,
    PHYLO_CLADES, PHYLO_QC, PHYLO_RESULTS, PHYLO_REPORT,
    GC_DATASET_DESIGN, GC_ESTIMAND_DESIGN, GC_CLAIM_DESIGN,
    GUIX_CHANNELS, GUIX_MANIFEST, PHYLO_GUIX_MANIFEST, VGP_DESIGN,
    METHOD, GUIX_RUNNER, METHOD_TEST,
)

CLAIM_FIELDS = (
    "claim_id", "evidence_stratum", "classification", "conclusion",
    "estimand", "sampling_unit", "observed_scope", "evidence_artifacts",
    "uncertainty_and_covariance", "data_lineage", "forbidden_extrapolation",
)
CORE_FIELDS = (
    "row_id", "evidence_layer", "estimand", "sampling_unit", "status",
    "count_or_estimate", "uncertainty", "independence", "eligibility_role",
    "phylogenetic_treatment", "lineage",
)
GC_FIELDS = (
    "row_id", "branch", "execution_state", "estimand", "sampling_unit",
    "observed_scope", "estimate", "uncertainty", "claim_classification",
    "supported_statement", "forbidden_statement", "lineage",
)


class SynthesisError(RuntimeError):
    """A reviewed-input or closed-world invariant failed."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


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


def require(condition: bool, message: str) -> None:
    if not condition:
        raise SynthesisError(message)


def reconcile_inputs() -> dict[str, object]:
    for path in INPUT_PATHS:
        require(path.is_file(), f"missing reviewed input: {path}")

    core = json.loads(CORE_SUMMARY.read_text(encoding="utf-8"))
    mirror = json.loads(MIRROR_SUMMARY.read_text(encoding="utf-8"))
    review = json.loads(REVIEW_DECISION.read_text(encoding="utf-8"))
    manifest_rows = read_tsv(CORE_MANIFEST)
    pair_rows = [row for row in manifest_rows if row["record_type"] == "linked_pair"]
    catalog_rows = [row for row in manifest_rows if row["record_type"] == "catalog_row"]
    qc_rows = read_tsv(CORE_QC)
    mirror_rows = read_tsv(MIRROR_MANIFEST)
    review_gates = read_tsv(REVIEW_GATES)
    roster = read_tsv(REVIEW_SCALEOUT)

    require(len(catalog_rows) == 716 and len(pair_rows) == len(qc_rows) == 569,
            "core closed-world row multiplicity drift")
    require(len({row["catalog_row"] for row in catalog_rows}) == 716,
            "catalog rows are not one-to-one")
    require(len({row["pair_id"] for row in pair_rows}) == 569,
            "pair identifiers are not one-to-one")
    require({row["pair_id"] for row in pair_rows} == {row["pair_id"] for row in qc_rows},
            "pair/QC identifiers differ")
    require(core["catalog"] == {
        "commit": "dc1b2af5a7741b97d66fb10cb2bce97f41765cdf",
        "sha256": "9c58420484a8b76a2d6175b7c26bf709e68bdc726a67fc7541b8c2b5a2fc13a4",
        "rows": 716,
    }, "Freeze 1 identity drift")
    require(core["eligible_pairs"] == core["completed_pairs"] == core["biological_estimates"] == 0,
            "unreviewed core biological results appeared")
    require(core["slurm_jobs_submitted"] == core["atomic_promotions"] == 0,
            "unreviewed scale-out execution appeared")
    require(core["full_biological_scaleout_authorized"] is False,
            "full scale-out authorization appeared")
    require(all(row["core_eligible"] == "false" for row in qc_rows),
            "an unreviewed eligible pair appeared")
    require(all(row["same_pair_psmc_independent_evidence"] == "false" for row in qc_rows),
            "same-pair PSMC was promoted to independent evidence")
    require(all(row["annotation_absence_core_veto"] == "false" for row in qc_rows),
            "annotation became a core eligibility gate")
    require(all(row["core_result_status"] == row["psmc_status"] == "NOT_RUN" for row in pair_rows),
            "unreviewed core or PSMC result appeared")
    require(all(row["annotation_partition_status"] == "NOT_RUN_CORE_NOT_PASSED" for row in pair_rows),
            "unreviewed exact-annotation result appeared")
    require(Counter(row["disposition"] for row in pair_rows) == {"failed": 566, "excluded": 3},
            "pair disposition count drift")
    require(Counter(row["primary_reason_code"] for row in pair_rows) == {
        "UPSTREAM_SCALEOUT_NOT_AUTHORIZED": 566, "CATALOG_SELF_LINK_NOT_PAIR": 3,
    }, "pair reason count drift")
    require(Counter(row["confidence_tier"] for row in pair_rows) == {"UNASSIGNED": 566, "X": 3},
            "confidence-tier count drift")

    names = [row["scientific_name"] for row in catalog_rows]
    require(len(set(names)) == 714, "frozen scientific-name reconciliation drift")
    require(mirror["catalog_reconciliation"]["released_rows"] == 581 and
            mirror["catalog_reconciliation"]["unreleased_rows"] == 135,
            "mirror catalog reconciliation drift")
    require(
        len(mirror_rows) == 47_870
        and sum(int(row["size_bytes"]) for row in mirror_rows if row["object_type"] == "file")
        == 3_916_877_494_936
        and {row["state"] for row in mirror_rows}
        == {"verified", "reused", "verified_upstream_conflict"}
        and all(
            row.get("canonical_vgp_root") == "/moosefs/erikg/vgp"
            and row.get("mirror_root") == "/moosefs/erikg/vgp/freeze1"
            for row in mirror_rows
        )
        and mirror["live_progress"]["remaining_files"] == 0
        and mirror["live_progress"]["remaining_bytes"] == 0,
        "mirror terminal closed-world canonical-root reconciliation drift",
    )
    require(len(roster) == 16 and Counter(row["roster_type"] for row in roster) == {
        "primary": 10, "alternate": 6,
    }, "review roster drift")
    require(review["primary_passes"] == 0 and review["primary_slots_reviewed"] == 10,
            "pilot review pass count drift")
    require(review["program_decision"] == "CONDITIONAL_GO" and
            review["full_biological_scaleout_authorized"] is False,
            "pilot review authorization drift")
    require(Counter(row["status"] for row in review_gates) == {
        "PASS": 15, "NOT_REACHED": 28, "FAIL": 11, "FAIL_NOT_ESTIMABLE": 2,
    }, "pilot review-gate count drift")

    direct_datasets = read_tsv(DIRECT_DATASETS)
    direct_tracts = read_tsv(DIRECT_TRACTS)
    direct_summary = read_tsv(DIRECT_SUMMARY)
    require(len(direct_datasets) == 2 and len(direct_tracts) == 58 and len(direct_summary) == 27,
            "direct-pilot manifest count drift")
    require(all(row["analysis_status"] == "PUBLISHED_DIRECTIONAL_CANDIDATE_NOT_RAW_REVALIDATED" and
                row["direct_rate_inclusion"] == "EXCLUDED" for row in direct_tracts),
            "published direct candidates were promoted")
    validated = next(row for row in direct_summary if row["estimand"] == "VALIDATED_ALLELIC_EVENTS")
    require(validated["status"] == "NOT_ESTIMABLE_INPUT_GATE" and validated["numerator"] == "0 admitted",
            "direct validated-event gate drift")
    direct_rate_rows = [row for row in direct_summary if row["estimand"].startswith("D_EVT_RATE")]
    require(direct_rate_rows and all(row["estimate"].startswith("NOT_ESTIMABLE") for row in direct_rate_rows),
            "direct rate was numerically populated")

    phylo_clades = read_tsv(PHYLO_CLADES)
    phylo_qc = read_tsv(PHYLO_QC)
    phylo_results = read_tsv(PHYLO_RESULTS)
    require(len(phylo_clades) == 10 and Counter(row["panel_id"] for row in phylo_clades) == {
        "H01": 5, "H02": 5,
    }, "phylogenetic panel count drift")
    require(sum(int(row["observed_sequence_bytes"]) for row in phylo_clades) == 0,
            "verified phylogenetic sequence payload appeared")
    require(len(phylo_qc) == 56 and Counter(row["status"] for row in phylo_qc) == {
        "PASS": 3, "FAIL": 23, "NOT_EVALUABLE": 2, "NOT_RUN_UPSTREAM_GATE": 28,
    }, "phylogenetic QC count drift")
    require(len(phylo_results) == 23 and all(
        (row["status"], row["estimate_value"]) in {
            ("NOT_ESTIMABLE_INPUT_GATE", "NOT_ESTIMABLE_INPUT_GATE"),
            ("NOT_MEASURED", "NOT_MEASURED"),
        } for row in phylo_results),
            "phylogenetic biological estimate appeared")

    datasets = read_tsv(GC_DATASET_DESIGN)
    estimands = read_tsv(GC_ESTIMAND_DESIGN)
    for branch in ("population", "non_allelic"):
        branch_datasets = [row for row in datasets if row["branch"] == branch]
        branch_estimands = [row for row in estimands if row["branch"] == branch]
        require(branch_datasets and branch_estimands,
                f"{branch} design rows missing")
        require(all(row["execution_state"] in {"NOT_RUN_DESIGN_ONLY", "DESIGN_ONLY_ALTERNATE"}
                    for row in branch_datasets), f"{branch} execution manifest appeared")
        require(all(row["execution_state"] == "NOT_RUN_DESIGN_ONLY" for row in branch_estimands),
                f"{branch} estimand was promoted")

    input_digests = {str(path.relative_to(ROOT)): sha256_file(path) for path in INPUT_PATHS}
    for key, digest in core["input_digests"].items():
        source = {
            "design": "analysis/vgp_analysis_manifest.json",
            "guix_channels": "analysis/guix/vgp_10_pilot/channels.scm",
            "guix_manifest": "analysis/guix/vgp_10_pilot/manifest.scm",
            "mirror_manifest": "analysis/vgp_freeze1_mirror_manifest.tsv",
            "mirror_summary": "analysis/vgp_freeze1_mirror_summary.json",
            "review_decision": "analysis/vgp_10_pilot_review_decision.json",
            "review_resources": "analysis/vgp_10_pilot_scaleout_resource_manifest.tsv",
            "review_scaleout": "analysis/vgp_10_pilot_scaleout_manifest.tsv",
        }.get(key)
        if source and Path(ROOT / source).is_file():
            observed = sha256_file(ROOT / source)
            if key == "design" and observed != digest:
                # The core summary is an immutable historical refusal packet.
                # Its design digest predates the audited relocation of the
                # active VGP root, so preserve that packet and recognize only
                # the versioned canonical-root transition here.
                current_design = json.loads((ROOT / source).read_text(encoding="utf-8"))
                require(
                    current_design.get("data_root") == "/moosefs/erikg/vgp",
                    "embedded core design drift is not the canonical-root migration",
                )
            elif key in {"mirror_manifest", "mirror_summary"} and observed != digest:
                # The core packet predates the authorized real mirror.  Accept
                # only the terminal, closed-world canonical-root transition;
                # current mirror files are independently pinned below.
                current_summary = json.loads(MIRROR_SUMMARY.read_text(encoding="utf-8"))
                require(
                    current_summary.get("canonical_vgp_root") == "/moosefs/erikg/vgp"
                    and current_summary.get("mirror_root") == "/moosefs/erikg/vgp/freeze1"
                    and current_summary.get("catalog_reconciliation", {}).get("closed_world") is True
                    and current_summary.get("inventory_totals", {}).get("full_release", {}).get("files") == 43_371
                    and current_summary.get("inventory_totals", {}).get("full_release", {}).get("bytes") == 3_916_877_494_936
                    and current_summary.get("live_progress", {}).get("remaining_files") == 0
                    and current_summary.get("live_progress", {}).get("remaining_bytes") == 0,
                    f"embedded core {key} drift is not the completed canonical mirror transition",
                )
            else:
                require(observed == digest, f"embedded core input digest drift: {key}")

    return {
        "core": core, "mirror": mirror, "review": review,
        "catalog_rows": catalog_rows, "pair_rows": pair_rows,
        "direct_tracts": direct_tracts, "direct_summary": direct_summary,
        "phylo_clades": phylo_clades, "phylo_results": phylo_results,
        "input_digests": input_digests,
    }


def core_table_rows() -> list[dict[str, str]]:
    lineage = "vgp_core_scaleout_summary.json;vgp_core_scaleout_manifest.tsv"
    return [
        {"row_id": "CORE-DIVERSITY", "evidence_layer": "cross_species_core",
         "estimand": "callable heterozygosity", "sampling_unit": "one audited diploid H1/H2 pair",
         "status": "NOT_ESTIMABLE_ZERO_ELIGIBLE_PAIRS", "count_or_estimate": "0 admitted pairs; no diversity estimate",
         "uncertainty": "not computable; within-pair windows would not create species replicates",
         "independence": "pair is the biological unit", "eligibility_role": "pre-annotation core gates only",
         "phylogenetic_treatment": "PGLS/hierarchical tree covariance preregistered; not fit with n=0", "lineage": lineage},
        {"row_id": "CORE-PSMC", "evidence_layer": "descriptive_demography",
         "estimand": "unscaled PSMC trajectory", "sampling_unit": "same individual, H1/H2 differences, mask, and consensus as diversity",
         "status": "NOT_ESTIMABLE_ZERO_ELIGIBLE_PAIRS", "count_or_estimate": "0 trajectories; 0 bootstrap attempts",
         "uncertainty": "boundary-aware block bootstrap required after a core pass",
         "independence": "SAME_PAIR_NONINDEPENDENT_DESCRIPTIVE", "eligibility_role": "cannot validate or predict its paired diversity",
         "phylogenetic_treatment": "joint/multivariate pair-cluster model required across species", "lineage": lineage},
        {"row_id": "CORE-ANNOTATION", "evidence_layer": "functional_partition",
         "estimand": "exact-annotation diversity partitions", "sampling_unit": "passing pair x exact bound partition",
         "status": "NOT_RUN_CORE_NOT_PASSED", "count_or_estimate": "0 exact-annotation subset pairs",
         "uncertainty": "not computable", "independence": "nested partitions share the pair and callable mask",
         "eligibility_role": "POST_CORE_OPTIONAL_PARTITION_ONLY", "phylogenetic_treatment": "pair-nested partition effects; no retrospective eligibility",
         "lineage": "vgp_core_scaleout_manifest.tsv;vgp_10_pilot_review_decision.json"},
        {"row_id": "CORE-UNSCALED", "evidence_layer": "psmc_scaling",
         "estimand": "unscaled interval parameters", "sampling_unit": "passing pair",
         "status": "NOT_MATERIALIZED_CORE_NOT_RUN", "count_or_estimate": "primary object absent",
         "uncertainty": "not computable", "independence": "same-pair outcome",
         "eligibility_role": "must remain primary if later executed", "phylogenetic_treatment": "retain without conversion assumptions",
         "lineage": "vgp_core_scaleout_scaling_scenarios.tsv"},
        {"row_id": "CORE-SCENARIOS", "evidence_layer": "psmc_scaling",
         "estimand": "mutation-rate/generation-time scaled trajectories", "sampling_unit": "passing pair x explicitly sourced scenario",
         "status": "NOT_MATERIALIZED_NO_APPROVED_TIER", "count_or_estimate": "no mutation or generation values selected",
         "uncertainty": "future scenario envelope; never pooled into one curve", "independence": "deterministic transforms share unscaled trajectory",
         "eligibility_role": "sensitivity only", "phylogenetic_treatment": "species-specific sources; sensitivity bounds labeled",
         "lineage": "vgp_core_scaleout_scaling_scenarios.tsv"},
    ]


def gene_conversion_rows() -> list[dict[str, str]]:
    return [
        {"row_id": "GC-DIRECT", "branch": "direct_pedigree_or_gamete", "execution_state": "EXECUTED_PREFLIGHT_INPUT_GATE",
         "estimand": "directional allelic conversion event/rate and GC resolution bias", "sampling_unit": "independent complete meiosis; event cluster nested within meiosis",
         "observed_scope": "Arabidopsis D01: 13 tetrads/52 products; 58 published candidates (44 CO, 14 NCO), all excluded; 0 admitted validated events",
         "estimate": "NOT_ESTIMABLE", "uncertainty": "raw callability, FDR/FNR, paralog masking, and power unavailable",
         "claim_classification": "not identifiable", "supported_statement": "directional published candidates were reconciled but no direct rate or bias passed its gate",
         "forbidden_statement": "do not treat candidates or H1/H2 WS/SW states as direct conversion/biased transmission; do not transfer this plant pilot across vertebrates",
         "lineage": "direct_gene_conversion_dataset_manifest.tsv;direct_gene_conversion_tracts.tsv;direct_gene_conversion_summary.tsv"},
        {"row_id": "GC-POPULATION", "branch": "population_allele_frequency_spectrum", "execution_state": "NOT_RUN_DESIGN_ONLY",
         "estimand": "polarization-aware WS/SW SFS and model-dependent B", "sampling_unit": "unrelated diploid individual within a preregistered population; linked blocks for resampling",
         "observed_scope": "design rows only; no conforming population execution manifest", "estimate": "NOT_ESTIMABLE",
         "uncertainty": "future demographic/mutation/polarization/linkage model and block uncertainty",
         "claim_classification": "design-only", "supported_statement": "requirements are specified; no population gBGC claim exists",
         "forbidden_statement": "do not infer population B from H1/H2 states or turn B into a pedigree event rate",
         "lineage": "gene_conversion_dataset_manifest.tsv;gene_conversion_estimand_manifest.tsv"},
        {"row_id": "GC-PHYLOGENETIC", "branch": "historical_phylogenetic_substitution", "execution_state": "EXECUTED_PREFLIGHT_INPUT_GATE",
         "estimand": "branch-specific WS/SW substitution asymmetry and historical gBGC-like bias", "sampling_unit": "branch x callable single-copy alignment partition; chromosome/synteny blocks",
         "observed_scope": "H01/H02, 5 taxa each; 0 verified sequence bytes, 0 callable bases, 0 substitutions",
         "estimate": "NOT_ESTIMABLE", "uncertainty": "tree, polarization, context mutation, gaps, annotation, and block uncertainty not fit",
         "claim_classification": "not identifiable", "supported_statement": "a pinned metadata preflight exists but no historical bias estimate",
         "forbidden_statement": "do not turn historical asymmetry into event counts, present transmission distortion, or population B; semi-complete genomes require explicit fragmentation sensitivity",
         "lineage": "vgp_phylogenetic_gbgc_clade_manifest.tsv;vgp_phylogenetic_gbgc_qc.tsv;vgp_phylogenetic_gbgc_results.tsv"},
        {"row_id": "GC-NONALLELIC", "branch": "non_allelic_paralog", "execution_state": "NOT_RUN_DESIGN_ONLY",
         "estimand": "copy-aware paralog tract candidates/homogenization excess", "sampling_unit": "resolved copy nested within haplotype and individual; orthogroup",
         "observed_scope": "design rows only; no conforming copy-resolved execution manifest", "estimate": "NOT_ESTIMABLE",
         "uncertainty": "future copy tree, recurrent mutation, assembly, CNV, and individual-cluster uncertainty",
         "claim_classification": "design-only", "supported_statement": "copy-aware requirements are specified; no non-allelic biological claim exists",
         "forbidden_statement": "do not call paralog homogenization allelic meiotic conversion or gBGC transmission",
         "lineage": "gene_conversion_dataset_manifest.tsv;gene_conversion_estimand_manifest.tsv"},
    ]


def claim_rows() -> list[dict[str, str]]:
    return [
        {"claim_id": "C01", "evidence_stratum": "robust_cross_species_core", "classification": "supported",
         "conclusion": "The pinned Freeze 1 closed world is completely accounted for: 716 catalog rows and 569 linked entries.",
         "estimand": "inventory/accounting", "sampling_unit": "catalog row or linked entry", "observed_scope": "716 rows; 569 links",
         "evidence_artifacts": "vgp_core_scaleout_manifest.tsv;vgp_core_scaleout_summary.json", "uncertainty_and_covariance": "exact counts; no sampling uncertainty; duplicate names retained as rows",
         "data_lineage": "catalog commit dc1b2af5a7741b97d66fb10cb2bce97f41765cdf; SHA-256 9c58420484a8b76a2d6175b7c26bf709e68bdc726a67fc7541b8c2b5a2fc13a4",
         "forbidden_extrapolation": "inventory completeness is not biological representativeness or a diversity estimate"},
        {"claim_id": "C02", "evidence_stratum": "robust_cross_species_core", "classification": "supported",
         "conclusion": "No pair was eligible or completed and no biological core job or estimate was promoted.",
         "estimand": "authorization/QC disposition", "sampling_unit": "linked pair", "observed_scope": "566 failed non-self candidates; 3 excluded self-links",
         "evidence_artifacts": "vgp_core_scaleout_qc.tsv;vgp_core_scaleout_summary.json", "uncertainty_and_covariance": "exact operational status; not a confidence interval on diversity",
         "data_lineage": "review decision and all-planned mirror digests embedded in each pair row",
         "forbidden_extrapolation": "technical non-execution is neither low diversity nor absence of vertebrate variation"},
        {"claim_id": "C03", "evidence_stratum": "robust_cross_species_core", "classification": "not identifiable",
         "conclusion": "Compression of callable diversity across vertebrates is not identifiable in the executed evidence.",
         "estimand": "cross-species range/rank of callable heterozygosity", "sampling_unit": "one diploid individual/pair per species", "observed_scope": "0 eligible pairs",
         "evidence_artifacts": "vgp_core_scaleout_paper_pairs.tsv;vgp_core_scaleout_sensitivity.tsv", "uncertainty_and_covariance": "no interval or phylogenetic model can be fit at n=0; windows are not replicates",
         "data_lineage": "closed-world pair and QC ledgers", "forbidden_extrapolation": "do not infer compression from exclusions, planning priors, or targeted-gene estimates"},
        {"claim_id": "C04", "evidence_stratum": "descriptive_demography", "classification": "not identifiable",
         "conclusion": "No unscaled or scenario-scaled PSMC history was materialized.", "estimand": "unscaled and scaled PSMC trajectory",
         "sampling_unit": "same passing pair used for diversity", "observed_scope": "0 trajectories; no mutation/generation values selected",
         "evidence_artifacts": "vgp_core_scaleout_scaling_scenarios.tsv;vgp_core_scaleout_qc.tsv", "uncertainty_and_covariance": "future PSMC and diversity must be joint outcomes clustered by pair; scenario transforms are perfectly lineage-linked",
         "data_lineage": "UNSCALED_PRIMARY and SPECIES_SCENARIOS_REQUIRED_AFTER_REVIEW rows",
         "forbidden_extrapolation": "same-pair PSMC cannot independently predict, validate, or causally explain paired diversity"},
        {"claim_id": "C05", "evidence_stratum": "functional_partition", "classification": "not identifiable",
         "conclusion": "No exact-annotation partition passed because no core pair passed; annotation remained optional.",
         "estimand": "partition-specific callable diversity", "sampling_unit": "pair x exact bound partition", "observed_scope": "0 subset pairs",
         "evidence_artifacts": "vgp_core_scaleout_manifest.tsv;vgp_10_pilot_review_decision.json", "uncertainty_and_covariance": "partitions would be nested correlated outcomes within pair",
         "data_lineage": "annotation status retained for all 569 pair rows", "forbidden_extrapolation": "do not use annotation availability retrospectively to select core pairs"},
        {"claim_id": "C06", "evidence_stratum": "bounded_association", "classification": "bounded",
         "conclusion": "Resource envelopes and threshold sensitivities are planning bounds only.", "estimand": "operational scale-out envelope",
         "sampling_unit": "planning scenario, not species", "observed_scope": "low/base/high templates; zero biological jobs",
         "evidence_artifacts": "vgp_core_scaleout_wave_manifest.tsv;vgp_core_scaleout_sensitivity.tsv", "uncertainty_and_covariance": "not fitted to biological PSMC runtime; bounds exclude unknown eligibility",
         "data_lineage": "reviewed pilot resource manifests", "forbidden_extrapolation": "do not report planning envelopes as achieved scale, statistical power, or biological uncertainty"},
        {"claim_id": "C07", "evidence_stratum": "specialized_mechanistic_direct", "classification": "not identifiable",
         "conclusion": "The direct pilot admits zero validated events, so event rate, tract distribution, CO association, and GC bias are not estimable.",
         "estimand": "direct allelic event/rate/bias", "sampling_unit": "complete meiosis; event cluster", "observed_scope": "13 tetrads; 58 published candidates excluded",
         "evidence_artifacts": "direct_gene_conversion_tracts.tsv;direct_gene_conversion_summary.tsv", "uncertainty_and_covariance": "callable opportunity and FDR/FNR absent; linked markers are clustered",
         "data_lineage": "Slurm preflight 1781172; Guix channel 44bbfc24e4bcc48d0e3343cd3d83452721af8c36", "forbidden_extrapolation": "do not call an admitted count of zero biological absence or transfer a plant rate to vertebrates"},
        {"claim_id": "C08", "evidence_stratum": "suggestive_association", "classification": "suggestive",
         "conclusion": "The published 45 S versus 37 W marker tally is an audit signal for follow-up only, not biased transmission.",
         "estimand": "published candidate-marker direction tally", "sampling_unit": "linked marker within unvalidated candidate tract", "observed_scope": "82 W/S candidate markers plus 1 ambiguous/non-SNV",
         "evidence_artifacts": "direct_gene_conversion_summary.tsv;direct_gene_conversion_tracts.tsv", "uncertainty_and_covariance": "events and linked markers are non-independent; reciprocal detection uncalibrated; power threshold unmet",
         "data_lineage": "published supplement reconciliation", "forbidden_extrapolation": "do not treat marker imbalance as direct GC bias, a rate, or cross-vertebrate evidence"},
        {"claim_id": "C09", "evidence_stratum": "specialized_historical", "classification": "not identifiable",
         "conclusion": "Historical phylogenetic WS/SW substitution bias is not estimable because the panels contain zero verified sequence bytes.",
         "estimand": "branch substitution asymmetry/historical bias", "sampling_unit": "branch x partition with block resampling", "observed_scope": "2 panels; 10 taxa; 0 callable bases",
         "evidence_artifacts": "vgp_phylogenetic_gbgc_clade_manifest.tsv;vgp_phylogenetic_gbgc_results.tsv", "uncertainty_and_covariance": "tree, polarization, mutation-context, fragmentation, and block uncertainty not fit",
         "data_lineage": "preflight Slurm job 1781129; frozen accessions and mirror states", "forbidden_extrapolation": "do not turn historical substitutions into direct events, current bias, or population B"},
        {"claim_id": "C10", "evidence_stratum": "future_population", "classification": "design-only",
         "conclusion": "Population AFS gBGC remains NOT_RUN/DESIGN_ONLY.", "estimand": "polarized WS/SW SFS and B",
         "sampling_unit": "unrelated individual within population; linked blocks", "observed_scope": "no execution manifest",
         "evidence_artifacts": "gene_conversion_dataset_manifest.tsv;gene_conversion_estimand_manifest.tsv", "uncertainty_and_covariance": "future demographic, mutation, polarization, LD, and block uncertainty",
         "data_lineage": "P_SFS/P_B design rows", "forbidden_extrapolation": "do not infer B from H1/H2 and do not convert B to a pedigree rate"},
        {"claim_id": "C11", "evidence_stratum": "future_non_allelic", "classification": "design-only",
         "conclusion": "Non-allelic paralog conversion remains NOT_RUN/DESIGN_ONLY.", "estimand": "copy-aware tracts/homogenization",
         "sampling_unit": "copy nested in haplotype/individual and orthogroup", "observed_scope": "no execution manifest",
         "evidence_artifacts": "gene_conversion_dataset_manifest.tsv;gene_conversion_estimand_manifest.tsv", "uncertainty_and_covariance": "future copy-tree, assembly/CNV, mutation, and individual-cluster uncertainty",
         "data_lineage": "N_TRACT/N_HOM design rows", "forbidden_extrapolation": "paralog similarity is not allelic conversion or biased transmission"},
        {"claim_id": "C12", "evidence_stratum": "cross_species_model", "classification": "not identifiable",
         "conclusion": "A phylogenetically aware diversity-demography comparison is specified but cannot be fit with zero eligible pairs.",
         "estimand": "joint cross-species association", "sampling_unit": "pair/species with phylogenetic covariance", "observed_scope": "n=0",
         "evidence_artifacts": "vgp_core_scaleout_summary.json;vgp_comprehensive_table_core.tsv", "uncertainty_and_covariance": "tree covariance plus shared-pair residual covariance required; bootstrap windows nested in pair",
         "data_lineage": "pre-registered analysis contract", "forbidden_extrapolation": "do not run naive OLS, treat windows as taxa, or place same-pair PSMC as independent predictor"},
    ]


def closed_world_svg() -> str:
    return """<svg xmlns="http://www.w3.org/2000/svg" width="1200" height="620" viewBox="0 0 1200 620" role="img">
<title>VGP Freeze 1 closed-world reconciliation</title>
<desc>Flow from 716 catalog rows to 569 linked entries and zero eligible pairs, estimates, or annotation subsets.</desc>
<rect width="1200" height="620" fill="#ffffff"/><text x="60" y="65" font-family="sans-serif" font-size="30" font-weight="bold">VGP Freeze 1: audited evidence, not biological zeroes</text>
<g font-family="sans-serif" text-anchor="middle"><rect x="60" y="140" width="240" height="150" rx="12" fill="#dbeafe" stroke="#1d4ed8" stroke-width="3"/><text x="180" y="195" font-size="34" font-weight="bold">716</text><text x="180" y="230" font-size="20">catalog rows</text><text x="180" y="260" font-size="16">581 released · 135 unreleased</text>
<rect x="365" y="140" width="240" height="150" rx="12" fill="#e0f2fe" stroke="#0369a1" stroke-width="3"/><text x="485" y="195" font-size="34" font-weight="bold">569</text><text x="485" y="230" font-size="20">linked entries</text><text x="485" y="260" font-size="16">566 non-self · 3 self-links</text>
<rect x="670" y="140" width="210" height="150" rx="12" fill="#fef3c7" stroke="#b45309" stroke-width="3"/><text x="775" y="195" font-size="34" font-weight="bold">0</text><text x="775" y="230" font-size="20">eligible pairs</text><text x="775" y="260" font-size="16">mirror verified objects: 0</text>
<rect x="945" y="140" width="195" height="150" rx="12" fill="#fee2e2" stroke="#b91c1c" stroke-width="3"/><text x="1042" y="195" font-size="34" font-weight="bold">0</text><text x="1042" y="230" font-size="20">core estimates</text><text x="1042" y="260" font-size="16">0 PSMC · 0 annotation</text>
<path d="M300 215 H365 M605 215 H670 M880 215 H945" stroke="#334155" stroke-width="4" marker-end="url(#a)"/></g>
<defs><marker id="a" markerWidth="10" markerHeight="10" refX="8" refY="3" orient="auto"><path d="M0,0 L0,6 L9,3 z" fill="#334155"/></marker></defs>
<rect x="90" y="380" width="1020" height="145" rx="12" fill="#f8fafc" stroke="#64748b" stroke-width="2"/><text x="600" y="425" text-anchor="middle" font-family="sans-serif" font-size="22" font-weight="bold">Interpretation boundary</text><text x="600" y="463" text-anchor="middle" font-family="sans-serif" font-size="19">0 eligible pairs means NOT ESTIMABLE—not low diversity, no diversity, or a compressed range.</text><text x="600" y="498" text-anchor="middle" font-family="sans-serif" font-size="17">Same-pair PSMC is descriptive and non-independent; annotation is optional after core acceptance.</text>
</svg>\n"""


def evidence_svg() -> str:
    rows = [
        ("Direct pedigree/gamete", "PREFLIGHT INPUT GATE", "58 published candidates excluded; rate not estimable", "#fee2e2", "#b91c1c"),
        ("Population AFS gBGC", "NOT_RUN / DESIGN_ONLY", "No population execution manifest; no B estimate", "#f1f5f9", "#475569"),
        ("Historical phylogenetic", "PREFLIGHT INPUT GATE", "0 verified sequence bytes; no substitution estimate", "#fee2e2", "#b91c1c"),
        ("Non-allelic paralogs", "NOT_RUN / DESIGN_ONLY", "No copy-resolved execution manifest; no tract claim", "#f1f5f9", "#475569"),
    ]
    blocks = []
    for index, (label, state, detail, fill, stroke) in enumerate(rows):
        y = 120 + index * 125
        blocks.append(f'<rect x="60" y="{y}" width="1080" height="95" rx="10" fill="{fill}" stroke="{stroke}" stroke-width="2"/>')
        blocks.append(f'<text x="90" y="{y + 34}" font-family="sans-serif" font-size="22" font-weight="bold">{label}</text>')
        blocks.append(f'<text x="1110" y="{y + 34}" text-anchor="end" font-family="sans-serif" font-size="18" font-weight="bold" fill="{stroke}">{state}</text>')
        blocks.append(f'<text x="90" y="{y + 70}" font-family="sans-serif" font-size="18">{detail}</text>')
    return """<svg xmlns="http://www.w3.org/2000/svg" width="1200" height="680" viewBox="0 0 1200 680" role="img">
<title>Four distinct gene-conversion evidence branches</title>
<desc>Direct and historical pilots stopped at input gates; population and non-allelic branches remain design only.</desc>
<rect width="1200" height="680" fill="#ffffff"/><text x="60" y="65" font-family="sans-serif" font-size="30" font-weight="bold">Gene-conversion evidence: estimands must not be exchanged</text>
""" + "\n".join(blocks) + """
<text x="600" y="640" text-anchor="middle" font-family="sans-serif" font-size="17">H1/H2 WS/SW states support none of these directional, population, historical, or copy-aware claims.</text>
</svg>\n"""


def synthesis_markdown(manifest: Mapping[str, object]) -> str:
    digests = manifest["input_digests"]
    return f"""# Comprehensive VGP evidence synthesis

**Evidence freeze:** 2026-07-18 UTC

**WG task:** `synthesize-vgp-program`

**Decision:** closed-world reconciliation complete; biological core and all gene-conversion effects remain fail-closed.

## Executive result

The robust result is an audit result, not a vertebrate-diversity effect. The pinned VGP Freeze 1 catalog contains exactly 716 rows (714 unique scientific names; the two duplicate-name occurrences remain separate release rows), of which 581 are released and 135 unreleased in the frozen mirror inventory. It yields 569 catalog-linked entries: 566 distinct non-self candidates and three self-links. None is eligible or complete under the reviewed gates. The mirror has 47,870 planned objects but zero verified or reused objects. Thus there are zero admitted callable-diversity estimates, zero PSMC trajectories, zero exact-annotation subsets, zero biological scale-out jobs, and no estimate of diversity-range compression.

Technical non-execution is not low diversity and is not evidence against heterozygosity. The analysis therefore classifies compression of callable diversity, demographic-shape contrasts, functional partition differences, and phylogenetically adjusted cross-species associations as **not identifiable**.

## Closed-world reconciliation

| Layer | Frozen count/state | Interpretation |
| --- | ---: | --- |
| Freeze 1 catalog | 716 rows; 714 unique names | Release-row membership is the sampling frame. |
| Frozen mirror | 581 released, 135 unreleased; 47,870 planned objects; 0 verified/reused | Inventory and transport evidence only. |
| Linked-pair ledger | 569 entries: 566 non-self, 3 self-links | Every entry has one QC/status row. |
| Pair disposition | 566 failed `UPSTREAM_SCALEOUT_NOT_AUTHORIZED`; 3 excluded `CATALOG_SELF_LINK_NOT_PAIR` | Operational dispositions, not biological classifications. |
| Confidence tier | 566 `UNASSIGNED`; 3 `X` | No Tier A/B/C biological result. |
| Ten-slot review | 0/10 primary passes; 6 alternates retained; 15 PASS, 11 FAIL, 2 FAIL_NOT_ESTIMABLE, 28 NOT_REACHED gates | `CONDITIONAL_GO` means bounded repair/re-pilot only, not scale-out. |
| Core / PSMC / exact annotation | 0 / 0 / 0 | All remain not run/not estimable. |

Figure 1 is `analysis/vgp_comprehensive_figure_closed_world.svg`; the exact paper table is `analysis/vgp_comprehensive_table_core.tsv`.

## Diversity, PSMC, functional partitions, and scaling

The biological sampling unit is one audited diploid individual represented by an exact same-individual H1/H2 pair. Windows and callable blocks quantify within-pair uncertainty; they do not manufacture independent individuals. A future cross-species analysis must use a versioned species tree and a phylogenetic mixed/PGLS or equivalent hierarchical covariance model, with pair/species as the outer unit and technical strata as covariates.

Assembly-derived diversity and PSMC are explicitly non-independent because they reuse the individual, H1/H2 differences, callable mask, and consensus. They must be represented as joint correlated outcomes or descriptive paired summaries with pair-clustered resampling. Same-pair PSMC must never be placed as an independent predictor or validator of same-pair diversity.

Annotation is not an eligibility gate. Exact accession/version and coordinate-dictionary binding may partition an already accepted core result, and partition effects remain nested within the pair. No non-annotated result may be removed retrospectively. Here no core pair passed, so the exact-annotation subset is empty for a pre-result reason.

The scenario ledger preserves both required scales. `UNSCALED_PRIMARY` is absent because core did not run. `SPECIES_SCENARIOS_REQUIRED_AFTER_REVIEW` is also absent: no mutation rate or generation time was selected. Future scenarios must cite species-specific sources, be reported separately with assumption bounds, and leave the unscaled trajectory primary.

## Four non-interchangeable gene-conversion branches

| Evidence branch | Execution state | Sampling unit | Current claim |
| --- | --- | --- | --- |
| Direct pedigree/gamete events | Executed metadata/table preflight; input gate | independent complete meiosis; events clustered within meiosis | D01 reconciles 13 tetrads/52 products and 58 published directional candidates (44 CO-associated, 14 NCO), all excluded. Zero events are admitted, so rates, tracts, CO association, and GC bias are not estimable. |
| Population allele-frequency-spectrum gBGC | **NOT_RUN/DESIGN_ONLY** | unrelated individuals within one population; LD blocks | No conforming execution manifest, frequency spectrum, or B estimate exists. |
| Historical phylogenetic substitution bias | Executed metadata preflight; input gate | branch by callable single-copy partition; chromosome/synteny blocks | H01/H02 freeze ten taxa but have zero verified sequence bytes, callable bases, or substitutions. No historical bias is estimable. |
| Non-allelic conversion among paralogs | **NOT_RUN/DESIGN_ONLY** | resolved copy nested in haplotype/individual and orthogroup | No copy-resolved execution manifest, tract, or homogenization estimate exists. |

The direct pilot's published-marker audit contains 45 S-resolved and 37 W-resolved markers plus one ambiguous/non-SNV marker. This is only a suggestive follow-up signal: linked markers and event candidates are not independent, all candidates fail raw/paralog/error gates, reciprocal detection is uncalibrated, and the registered power threshold is unmet. H1/H2 heterozygous WS/SW states cannot replace parent-of-origin or four-product direction. Historical WS/SW asymmetry, if later estimated, would be long-term substitution evidence only. Because the frozen historical panels may include semi-complete assemblies, any later model must preserve callable-gap accounting and the preregistered fragmentation sensitivity instead of treating missing sequence as absence of substitutions. Population B would remain a model-dependent frequency parameter, and paralog homogenization would remain non-allelic. Figure 2 (`analysis/vgp_comprehensive_figure_evidence.svg`) preserves these boundaries.

## Conclusion classes and claim boundary

The complete claim-by-claim ledger is `analysis/vgp_comprehensive_claim_ledger.tsv`.

- **Supported:** immutable release accounting, pair/QC reconciliation, and explicit non-execution/authorization states.
- **Bounded:** planning/resource and future scenario envelopes, never biological effects.
- **Suggestive:** the unvalidated published direct-candidate marker tally, solely as a follow-up target.
- **Design-only:** population SFS/B and non-allelic copy-aware programs.
- **Not identifiable:** callable-diversity compression, PSMC histories, functional effects, direct rates/bias, historical substitution bias, and cross-species associations.

No row supports extrapolating a direct single-species/plant or human rate across vertebrates. A future sensitivity could introduce an external rate only with provenance, taxonomic relevance, and separately labeled lower/base/upper bounds; it could not become an observed vertebrate result.

## Immutable lineage and reproduction

The final manifest binds {len(digests)} reviewed inputs, including core, review, mirror, specialized, design, and environment files, and binds every paper output by SHA-256. The final-analysis environment is GNU Guix channel commit `44bbfc24e4bcc48d0e3343cd3d83452721af8c36`; `analysis/guix/channels.scm` has SHA-256 `{digests['analysis/guix/channels.scm']}` and `analysis/guix/manifest.scm` has SHA-256 `{digests['analysis/guix/manifest.scm']}`.

Reproduce without acquisition or biological submission:

```bash
./analysis/run_vgp_comprehensive_synthesis_guix.sh
./analysis/run_vgp_comprehensive_synthesis_guix.sh --tests
./analysis/run_vgp_comprehensive_synthesis_guix.sh --all-tests
```

The generator is intentionally a local, read-only reconciliation over reviewed inputs. It refuses any promoted core result, any same-pair independence claim, any annotation veto, any admitted direct candidate, any phylogenetic estimate without inputs, or any population/non-allelic execution-state promotion.
"""


def generate(output_dir: Path = ANALYSIS) -> dict[str, object]:
    reconciled = reconcile_inputs()
    output_dir = Path(output_dir)
    core_rows = core_table_rows()
    gc_rows = gene_conversion_rows()
    claims = claim_rows()
    atomic_tsv(output_dir / CORE_TABLE.name, CORE_FIELDS, core_rows)
    atomic_tsv(output_dir / GENE_CONVERSION_TABLE.name, GC_FIELDS, gc_rows)
    atomic_tsv(output_dir / CLAIM_LEDGER.name, CLAIM_FIELDS, claims)
    atomic_text(output_dir / CLOSED_WORLD_FIGURE.name, closed_world_svg())
    atomic_text(output_dir / EVIDENCE_FIGURE.name, evidence_svg())

    pair_rows = reconciled["pair_rows"]
    catalog_rows = reconciled["catalog_rows"]
    manifest: dict[str, object] = {
        "schema_version": "vgp-comprehensive-final-manifest-v2.0.0",
        "supersedes_manifest_without_mutating_history":
            "analysis/vgp_comprehensive_final_manifest.json",
        "evidence_freeze_utc": "2026-07-18T00:00:00Z",
        "task_id": "synthesize-vgp-program",
        "decision": "CLOSED_WORLD_RECONCILED_BIOLOGICAL_EFFECTS_NOT_IDENTIFIABLE",
        "execution_boundary": "REVIEWED_RECONCILIATION_ONLY_NO_ACQUISITION_NO_BIOLOGICAL_COMPUTE",
        "catalog": {
            "commit": "dc1b2af5a7741b97d66fb10cb2bce97f41765cdf",
            "sha256": "9c58420484a8b76a2d6175b7c26bf709e68bdc726a67fc7541b8c2b5a2fc13a4",
            "rows": 716,
        },
        "closed_world": {
            "catalog_rows": 716, "catalog_unique_scientific_names": 714,
            "catalog_duplicate_name_rows": 2, "released_catalog_rows": 581,
            "unreleased_catalog_rows": 135, "mirror_inventory_objects": 47_870,
            "mirror_verified_or_reused_objects": 0, "linked_haplotype_entries": 569,
            "distinct_nonself_pair_candidates": 566, "self_links_excluded": 3,
            "eligible_pairs": 0, "completed_pairs": 0, "exact_annotation_subset_pairs": 0,
            "biological_estimates": 0, "biological_jobs_submitted": 0,
        },
        "core_pair_statuses": {
            "dispositions": dict(sorted(Counter(row["disposition"] for row in pair_rows).items())),
            "primary_reason_codes": dict(sorted(Counter(row["primary_reason_code"] for row in pair_rows).items())),
            "confidence_tiers": dict(sorted(Counter(row["confidence_tier"] for row in pair_rows).items())),
        },
        "catalog_row_statuses": {
            "dispositions": dict(sorted(Counter(row["disposition"] for row in catalog_rows).items())),
            "primary_reason_codes": dict(sorted(Counter(row["primary_reason_code"] for row in catalog_rows).items())),
        },
        "pilot_review": {
            "decision": "CONDITIONAL_GO", "authorization": "BOUNDED_REPAIR_AND_TEN_SLOT_REPILOT_ONLY",
            "primary_slots": 10, "primary_passes": 0, "alternates": 6,
            "gate_counts": {"PASS": 15, "FAIL": 11, "FAIL_NOT_ESTIMABLE": 2, "NOT_REACHED": 28},
            "full_scaleout_authorized": False,
        },
        "analysis_contract": {
            "sampling_unit": "one audited diploid individual/H1-H2 pair",
            "windows_are_independent_species_replicates": False,
            "same_pair_psmc_independent_evidence": False,
            "same_pair_model": "joint_multivariate_outcomes_with_pair_cluster_covariance",
            "annotation_absence_core_veto": False,
            "annotation_role": "post_core_optional_exact_bound_partition",
            "cross_species_phylogenetic_model": "versioned_tree_PGLS_or_hierarchical_phylogenetic_covariance",
            "cross_species_model_status": "NOT_IDENTIFIABLE_ZERO_ELIGIBLE_PAIRS",
            "unscaled_psmc_status": "NOT_MATERIALIZED_CORE_NOT_RUN",
            "mutation_generation_scenario_status": "NOT_MATERIALIZED_NO_APPROVED_TIER",
        },
        "gene_conversion": {
            "direct": {"state": "EXECUTED_PREFLIGHT_INPUT_GATE", "tetrads": 13,
                       "products": 52, "published_candidates_excluded": 58,
                       "validated_events": 0, "estimate": "NOT_ESTIMABLE"},
            "population": {"state": "NOT_RUN_DESIGN_ONLY", "estimate": "NOT_ESTIMABLE"},
            "phylogenetic": {"state": "EXECUTED_PREFLIGHT_INPUT_GATE", "panels": 2,
                              "taxa": 10, "verified_sequence_bytes": 0,
                              "estimate": "NOT_ESTIMABLE"},
            "non_allelic": {"state": "NOT_RUN_DESIGN_ONLY", "estimate": "NOT_ESTIMABLE"},
        },
        "claim_classification_counts": dict(sorted(Counter(row["classification"] for row in claims).items())),
        "input_digests": reconciled["input_digests"],
        "environment": {
            "guix_channel_commit": "44bbfc24e4bcc48d0e3343cd3d83452721af8c36",
            "channels": "analysis/guix/channels.scm",
            "manifest": "analysis/guix/manifest.scm",
            "pure_environment_required": True,
        },
        "reproducible_commands": [
            "./analysis/run_vgp_comprehensive_synthesis_guix.sh",
            "./analysis/run_vgp_comprehensive_synthesis_guix.sh --tests",
            "./analysis/run_vgp_comprehensive_synthesis_guix.sh --all-tests",
        ],
        "output_digests": {},
    }
    atomic_text(output_dir / SYNTHESIS.name, synthesis_markdown(manifest))
    manifest["output_digests"] = {
        str(Path("analysis") / path.name): sha256_file(output_dir / path.name)
        for path in OUTPUT_PATHS
    }
    atomic_json(output_dir / FINAL_MANIFEST.name, manifest)
    errors = validate_outputs(output_dir)
    if errors:
        raise SynthesisError("generated output validation failed: " + "; ".join(errors))
    return manifest


def validate_outputs(output_dir: Path = ANALYSIS, *, verify_digests: bool = True) -> list[str]:
    output_dir = Path(output_dir)
    errors: list[str] = []
    for filename in OUTPUT_FILENAMES:
        if not (output_dir / filename).is_file():
            errors.append(f"missing output: {filename}")
    if errors:
        return errors
    try:
        manifest = json.loads((output_dir / FINAL_MANIFEST.name).read_text(encoding="utf-8"))
        gc_rows = read_tsv(output_dir / GENE_CONVERSION_TABLE.name)
        claims = read_tsv(output_dir / CLAIM_LEDGER.name)
        core_rows = read_tsv(output_dir / CORE_TABLE.name)
    except Exception as exc:
        return [f"unreadable output: {exc}"]

    if len(gc_rows) != 4 or [row["branch"] for row in gc_rows] != [
        "direct_pedigree_or_gamete", "population_allele_frequency_spectrum",
        "historical_phylogenetic_substitution", "non_allelic_paralog",
    ]:
        errors.append("four gene-conversion rows are missing, reordered, or conflated")
    for branch in ("population_allele_frequency_spectrum", "non_allelic_paralog"):
        selected = [row for row in gc_rows if row["branch"] == branch]
        if len(selected) != 1 or selected[0]["execution_state"] != "NOT_RUN_DESIGN_ONLY" or selected[0]["estimate"] != "NOT_ESTIMABLE":
            errors.append(f"{branch} must remain NOT_RUN_DESIGN_ONLY with no estimate")
    direct = next((row for row in gc_rows if row["branch"] == "direct_pedigree_or_gamete"), None)
    if direct and direct["estimate"] != "NOT_ESTIMABLE":
        errors.append("direct-pilot candidate inventory was promoted to an estimate")
    phylo = next((row for row in gc_rows if row["branch"] == "historical_phylogenetic_substitution"), None)
    if phylo and phylo["estimate"] != "NOT_ESTIMABLE":
        errors.append("phylogenetic preflight was promoted to an estimate")
    if {row["classification"] for row in claims} - {
        "supported", "bounded", "suggestive", "design-only", "not identifiable",
    }:
        errors.append("claim ledger contains an unauthorized classification")
    for row in claims:
        for field in ("evidence_artifacts", "sampling_unit", "uncertainty_and_covariance", "data_lineage", "forbidden_extrapolation"):
            if not row[field]:
                errors.append(f"{row['claim_id']} lacks {field}")
    core_by_id = {row["row_id"]: row for row in core_rows}
    if core_by_id.get("CORE-PSMC", {}).get("independence") != "SAME_PAIR_NONINDEPENDENT_DESCRIPTIVE":
        errors.append("same-pair PSMC non-independence is not preserved")
    if core_by_id.get("CORE-ANNOTATION", {}).get("eligibility_role") != "POST_CORE_OPTIONAL_PARTITION_ONLY":
        errors.append("annotation optionality is not preserved")
    contract = manifest.get("analysis_contract", {})
    if contract.get("same_pair_psmc_independent_evidence") is not False:
        errors.append("manifest promotes same-pair PSMC to independent evidence")
    if contract.get("annotation_absence_core_veto") is not False:
        errors.append("manifest turns annotation into a core veto")
    if manifest.get("closed_world", {}).get("biological_estimates") != 0:
        errors.append("manifest contains a biological estimate")
    if manifest.get("gene_conversion", {}).get("population", {}).get("state") != "NOT_RUN_DESIGN_ONLY":
        errors.append("population manifest state was promoted")
    if manifest.get("gene_conversion", {}).get("non_allelic", {}).get("state") != "NOT_RUN_DESIGN_ONLY":
        errors.append("non-allelic manifest state was promoted")
    if verify_digests:
        for relative, expected in manifest.get("output_digests", {}).items():
            path = output_dir / Path(relative).name
            if not path.is_file() or sha256_file(path) != expected:
                errors.append(f"output digest mismatch: {relative}")
    return errors


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=ANALYSIS)
    parser.add_argument("--validate-only", action="store_true")
    args = parser.parse_args()
    if args.validate_only:
        errors = validate_outputs(args.output_dir)
        if errors:
            raise SynthesisError("; ".join(errors))
    else:
        generate(args.output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
