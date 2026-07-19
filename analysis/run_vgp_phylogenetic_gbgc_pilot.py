#!/usr/bin/env python3
"""Fail-closed execution ledger for the two frozen VGP phylogenetic gBGC pilots.

The approved H01/H02 analyses require ten exact assemblies, verified sequence
digests, coordinate-compatible annotations, and callable 1:1 orthology.  The
upstream mirror stopped before bulk transfer because enforceable quota evidence
was unavailable.  This program therefore emits an auditable *preflight result*,
not biological estimates.  It refuses to write repository artifacts unless it
is running in the pinned Guix environment inside Slurm, and it refuses to reuse
this blocked-output path if verified sequence inputs later appear.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import sys
from typing import Iterable


ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = ROOT / "analysis"
DATASETS = ANALYSIS / "gene_conversion_dataset_manifest.tsv"
SOURCE_INVENTORY = ANALYSIS / "vgp_freeze1_source_inventory.tsv"
MIRROR_MANIFEST = ANALYSIS / "vgp_freeze1_mirror_manifest.tsv"
MIRROR_SUMMARY = ANALYSIS / "vgp_freeze1_mirror_summary.json"

CLADE_OUTPUT = "vgp_phylogenetic_gbgc_clade_manifest.tsv"
QC_OUTPUT = "vgp_phylogenetic_gbgc_qc.tsv"
RESULT_OUTPUT = "vgp_phylogenetic_gbgc_results.tsv"
REPORT_OUTPUT = "vgp_phylogenetic_gbgc_pilot.md"

CATALOG_COMMIT = "dc1b2af5a7741b97d66fb10cb2bce97f41765cdf"
CATALOG_SHA256 = "9c58420484a8b76a2d6175b7c26bf709e68bdc726a67fc7541b8c2b5a2fc13a4"
CHANNEL_COMMIT = "44bbfc24e4bcc48d0e3343cd3d83452721af8c36"
NCBI_LICENSE = (
    "NCBI molecular sequence data carry no NCBI-imposed use restriction; "
    "submitter rights are not transferred; VGP/G10K attribution/publication status remains assembly-specific"
)
NCBI_POLICY_URL = "https://www.ncbi.nlm.nih.gov/home/about/policies/"
UNRESOLVED_CHECKSUM = "NOT_RESOLVED_NO_LOCAL_SEQUENCE_OR_PUBLISHED_DIGEST"
NOT_ESTIMABLE = "NOT_ESTIMABLE_INPUT_GATE"


PANELS = {
    "H01": [
        # species, role, accession, bp, contigs, catalog membership, annotation
        ("Spinachia_spinachia", "focal_ingroup", "GCA_048126635.1", 407541755, 278, "exact_main_haplotype", "GCA_048126635.1-GB_2025_08_04"),
        ("Pungitius_pungitius", "close_ingroup", "GCA_949316345.1", 480434100, 912, "exact_main_haplotype", "GCF_949316345.1-RS_2024_01"),
        ("Gasterosteus_aculeatus", "close_ingroup", "GCA_964276395.1", 511138727, 460, "exact_main_haplotype", "GCF_964276395.1-RS_2025_12"),
        ("Syngnathus_acus", "outgroup_1", "GCA_948146105.1", 359199823, 345, "exact_other_high_quality_haplotype", "NONE_ON_EXACT_ASSEMBLY"),
        ("Syngnathus_typhle", "outgroup_2", "GCA_048301445.1", 368283750, 229, "exact_main_haplotype", "GCA_048301445.1-GB_2025_08_16"),
    ],
    "H02": [
        ("Falco_naumanni", "focal_ingroup", "GCA_017639655.1", 1215702009, 588, "exact_main_haplotype", "GCF_017639655.2-RS_2026_06"),
        ("Falco_tinnunculus", "close_ingroup", "GCA_976974455.1", 1417225294, 1002, "not_in_freeze_catalog_close_relative", "NONE_ON_EXACT_ASSEMBLY"),
        ("Falco_peregrinus", "close_ingroup", "GCA_965282525.1", 1284654402, 728, "not_in_freeze_catalog_close_relative", "NONE_ON_EXACT_ASSEMBLY"),
        ("Falco_cherrug", "outgroup_1", "GCA_023634085.1", 1309366964, 540, "exact_main_haplotype", "GCF_023634085.1-RS_2023_04"),
        ("Falco_punctatus", "outgroup_2", "GCA_963210335.1", 1279260504, 722, "exact_main_haplotype", "GCF_963210335.1-RS_2025_09"),
    ],
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{path}: missing header")
        return reader.fieldnames, list(reader)


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    if not rows:
        raise ValueError(f"refusing empty output: {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def expected_accessions() -> set[str]:
    return {taxon[2] for taxa in PANELS.values() for taxon in taxa}


def audit_upstream() -> dict[str, object]:
    for path in (DATASETS, SOURCE_INVENTORY, MIRROR_MANIFEST, MIRROR_SUMMARY):
        if not path.is_file():
            raise ValueError(f"missing upstream artifact: {path}")

    _, design_rows = read_tsv(DATASETS)
    design = {row["dataset_id"]: row for row in design_rows}
    for panel_id, taxa in PANELS.items():
        if design.get(panel_id, {}).get("execution_state") != "AUTHORIZED_DOWNSTREAM_PILOT":
            raise ValueError(f"{panel_id}: not authorized by frozen design")
        for accession in (taxon[2] for taxon in taxa):
            if accession not in design[panel_id]["exact_accessions"]:
                raise ValueError(f"{panel_id}: accession drift for {accession}")

    summary = json.loads(MIRROR_SUMMARY.read_text(encoding="utf-8"))
    release = summary["release"]
    if release["catalog_commit"] != CATALOG_COMMIT or release["catalog_sha256"] != CATALOG_SHA256:
        raise ValueError("VGP Freeze 1 catalog identity drift")
    if release["guix_channel_commit"] != CHANNEL_COMMIT:
        raise ValueError("pinned Guix channel drift")
    launch = summary["bulk_launch"]
    storage = summary["storage"]
    if (
        summary.get("canonical_vgp_root") != "/moosefs/erikg/vgp"
        or summary.get("mirror_root") != "/moosefs/erikg/vgp/freeze1"
        or launch.get("slurm_jobs_launched") != 0
        or not str(launch.get("reason", "")).startswith(
            "capacity_write_and_inode_headroom_verified"
        )
        or storage.get("adequate") is not True
        or storage.get("quota_visibility_is_policy_gate") is not False
    ):
        raise ValueError("canonical mirror execution contract drift")

    _, source_rows = read_tsv(SOURCE_INVENTORY)
    _, mirror_rows = read_tsv(MIRROR_MANIFEST)
    wanted = expected_accessions()
    source_sequences = [
        row for row in source_rows
        if row["accession_version"] in wanted
        and row["sequence_subset"] in {"assembly_fasta", "assembly_2bit"}
    ]
    mirror_sequences = [
        row for row in mirror_rows
        if row["accession_version"] in wanted
        and row["sequence_subset"] in {"assembly_fasta", "assembly_2bit"}
    ]
    if any(
        row.get("canonical_vgp_root") != "/moosefs/erikg/vgp"
        or row.get("mirror_root") != "/moosefs/erikg/vgp/freeze1"
        for row in mirror_sequences
    ):
        raise ValueError("pilot sequence paths escaped the canonical mirror root")
    return {
        "summary": summary,
        "source_sequences": source_sequences,
        "mirror_sequences": mirror_sequences,
        "source_sha256": sha256_file(SOURCE_INVENTORY),
        "mirror_sha256": sha256_file(MIRROR_MANIFEST),
        "summary_sha256": sha256_file(MIRROR_SUMMARY),
    }


def clade_rows(audit: dict[str, object], slurm_job_id: str) -> list[dict[str, str]]:
    source_rows = audit["source_sequences"]
    assert isinstance(source_rows, list)
    by_accession: dict[str, list[dict[str, str]]] = {}
    for row in source_rows:
        by_accession.setdefault(row["accession_version"], []).append(row)
    rows: list[dict[str, str]] = []
    for panel_id, taxa in PANELS.items():
        for species, role, accession, span, contigs, catalog_status, annotation in taxa:
            inventory = by_accession.get(accession, [])
            paths = ";".join(sorted(row["source_relative_path"] for row in inventory)) or "NOT_IN_RELEASED_UCSC_INVENTORY"
            expected_bytes = sum(int(row["size_bytes"]) for row in inventory)
            annotation_status = (
                "NOT_AVAILABLE_ON_EXACT_ASSEMBLY"
                if annotation.startswith("NONE")
                else "NOT_VALIDATED_NO_SEQUENCE_DICTIONARY"
            )
            rows.append({
                "panel_id": panel_id,
                "taxon_id": species.lower(),
                "species": species.replace("_", " "),
                "polarization_role": role,
                "exact_assembly_accession": accession,
                "assembly_span_bp_design_metadata": str(span),
                "contig_count_design_metadata": str(contigs),
                "freeze_catalog_membership": catalog_status,
                "freeze_catalog_commit": CATALOG_COMMIT,
                "freeze_catalog_sha256": CATALOG_SHA256,
                "assembly_source_url": f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{accession}/dataset_report",
                "released_ucsc_inventory_status": "PRESENT" if inventory else "ABSENT",
                "released_sequence_object_count": str(len(inventory)),
                "released_sequence_expected_bytes": str(expected_bytes),
                "released_sequence_paths": paths,
                "mirror_state": "PLANNED_ZERO_BYTES" if inventory else "NOT_IN_RELEASED_UCSC_INVENTORY",
                "observed_sequence_bytes": "0",
                "upstream_sequence_checksum": UNRESOLVED_CHECKSUM,
                "local_sequence_sha256": UNRESOLVED_CHECKSUM,
                "license_or_reuse_terms": NCBI_LICENSE,
                "license_evidence_url": NCBI_POLICY_URL,
                "annotation_accession_release": annotation,
                "annotation_source_url": f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{accession}/dataset_report",
                "annotation_checksum": "NOT_RESOLVED_ANNOTATION_NOT_ACQUIRED_OR_USED",
                "annotation_license": "NOT_RESOLVED_ASSEMBLY_SPECIFIC_AUDIT_REQUIRED",
                "annotation_coordinate_evidence": annotation_status,
                "sequence_dictionary_evidence": "NOT_MEASURED_NO_LOCAL_SEQUENCE",
                "species_tree_source_and_version": "NOT_RESOLVED_NO_EXECUTABLE_TOPOLOGY_ARTIFACT",
                "species_tree_checksum": "NOT_RESOLVED_NO_EXECUTABLE_TOPOLOGY_ARTIFACT",
                "species_tree_license": "NOT_RESOLVED_NO_EXECUTABLE_TOPOLOGY_ARTIFACT",
                "orthology_source_and_version": "NOT_CREATED_INPUT_GATE",
                "orthology_checksum": "NOT_APPLICABLE_NO_ORTHOLOGY_PRODUCT",
                "orthology_license": "NOT_APPLICABLE_NO_ORTHOLOGY_PRODUCT",
                "copy_and_paralogy_control": "NOT_RUN; reciprocal_1to1_single_copy_synteny_required",
                "callable_alignment_status": "NOT_CREATED_INPUT_GATE",
                "exclusion_reason": "NO_VERIFIED_SEQUENCE_PAYLOAD;QUOTA_VISIBILITY_UNAVAILABLE;CHECKSUM_UNRESOLVED",
                "source_inventory_sha256": str(audit["source_sha256"]),
                "mirror_manifest_sha256": str(audit["mirror_sha256"]),
                "preflight_slurm_job_id": slurm_job_id,
            })
    return rows


def qc_rows(audit: dict[str, object], slurm_job_id: str) -> list[dict[str, str]]:
    source_rows = audit["source_sequences"]
    assert isinstance(source_rows, list)
    released = {row["accession_version"] for row in source_rows}
    rows: list[dict[str, str]] = []

    def add(panel: str, gate: str, category: str, requirement: str, observed: str,
            threshold: str, status: str, reason: str, evidence: str) -> None:
        rows.append({
            "panel_id": panel,
            "qc_id": gate,
            "category": category,
            "requirement": requirement,
            "observed": observed,
            "threshold": threshold,
            "status": status,
            "excluded_bases": "NOT_ESTIMABLE_NO_ALIGNMENT",
            "excluded_fraction": "NOT_ESTIMABLE_NO_ALIGNMENT",
            "reason_code": reason,
            "evidence": evidence,
            "preflight_slurm_job_id": slurm_job_id,
        })

    for panel, taxa in PANELS.items():
        panel_accessions = {taxon[2] for taxon in taxa}
        present = len(panel_accessions & released)
        add(panel, "G01_DESIGN_IDENTITY", "identity", "five frozen exact accession versions", "5 exact versions", "5", "PASS", "NONE", "gene_conversion_dataset_manifest.tsv")
        add(panel, "G02_RELEASED_INVENTORY", "acquisition", "all five sequence sources inventoried", f"{present}/5 accessions; {present * 2} FASTA/2bit objects", "5/5", "FAIL", "MISSING_RELEASED_SEQUENCE_SOURCE", "vgp_freeze1_source_inventory.tsv")
        add(panel, "G03_LOCAL_SEQUENCE", "acquisition", "five checksum-verified local assemblies", "0/5; 0 observed bytes", "5/5", "FAIL", "NO_VERIFIED_SEQUENCE_PAYLOAD", "vgp_freeze1_mirror_manifest.tsv")
        add(panel, "G04_SEQUENCE_CHECKSUM", "provenance", "published checksum plus local SHA-256 for every assembly", "0/5 resolved", "5/5", "FAIL", "CHECKSUM_UNRESOLVED", "vgp_freeze1_mirror_summary.json")
        add(panel, "G04B_LICENSE_AUDIT", "provenance", "exact assembly/annotation VGP publication status, attribution, and reuse audit", "general NCBI policy only; assembly-specific audit unresolved", "all used objects", "FAIL", "ASSEMBLY_SPECIFIC_LICENSE_AUDIT_UNRESOLVED", "clade manifest")
        add(panel, "G05_ANNOTATION_COORDINATES", "annotation", "exact native annotation or validated coordinate mapping and sequence dictionary", "0 exact dictionaries validated", "all used partitions", "FAIL", "ANNOTATION_DICTIONARY_UNAVAILABLE", "gene_conversion_dataset_manifest.tsv")
        add(panel, "G06_TAXON_POLARIZATION", "design", "at least 3 callable ingroups plus 2 callable outgroups", "0 callable ingroups; 0 callable outgroups", ">=3 ingroup; >=2 outgroup", "FAIL", "NO_CALLABLE_ALIGNMENT", "pre-result gate")
        add(panel, "G06B_TOPOLOGY_IDENTITY", "phylogeny", "versioned species topology/source/checksum/license plus alternate topology", "no executable topology artifact frozen", "all taxa and branches", "FAIL", "TOPOLOGY_PROVENANCE_UNRESOLVED", "design describes sensitivity but supplies no versioned Newick")
        add(panel, "G07_NEUTRAL_CALLABILITY", "callability", "reciprocal 1:1 neutral aligned bases callable in required taxa", "0 bases", ">=10000 bases", "FAIL", "NO_ORTHOLOGY_PRODUCT", "H_SUB preregistration")
        add(panel, "G08_H_SUB_COUNTS", "estimand", "directionally informative substitutions per branch/partition", "NOT_ESTIMABLE_NO_BRANCH_ASSIGNMENTS", ">=20 WS and >=20 SW", "FAIL", "NO_BRANCH_ASSIGNMENTS", "H_SUB preregistration")
        add(panel, "G09_H_GBGC_COUNTS", "estimand", "total directionally informative substitutions", "NOT_ESTIMABLE_NO_BRANCH_ASSIGNMENTS", ">=100", "FAIL", "NO_BRANCH_ASSIGNMENTS", "H_GBGC preregistration")
        add(panel, "G10_H_CLUSTER_CALLABILITY", "estimand", "branch-polarized WS substitutions and continuous callable sequence", "NOT_ESTIMABLE_WS; 0 callable bp", ">=50 WS and >=20000000 bp", "FAIL", "NO_BRANCH_ASSIGNMENTS", "H_CLUSTER preregistration")
        add(panel, "G11_FRAGMENTATION_MISSINGNESS", "assembly", "reason-coded aligned/callable/excluded bases by assembly and gap", "metadata contig counts only; no callable denominator", "all taxa and partitions", "NOT_EVALUABLE", "NO_ALIGNMENT_FOR_MISSINGNESS_MODEL", "clade manifest")

        controls = [
            ("C01_LABEL_SHUFFLE", "negative_control", "branch/taxon label shuffling"),
            ("C02_NULL_SIMULATION", "negative_control", "parametric nonstationary null with empirical gaps/GC/branch lengths"),
            ("C03_ALTERNATE_OUTGROUP", "polarization", "each outgroup separately and agreement analysis"),
            ("C04_ANCESTRAL_UNCERTAINTY", "polarization", "posterior integration, epsilon sensitivity, high-posterior and transversion-only"),
            ("C05_ALIGNMENT_FILTERS", "alignment", "alternate alignment-quality, repeat, duplication, and base-quality filters"),
            ("C06_MUTATION_BIAS", "mutation", "ancestral 3-mer/5-mer opportunities, GC-conservative rates, CpG sensitivity"),
            ("C07_MULTIPLE_HITS_GC_EQ", "substitution_model", "nonstationary model with multiple hits and GC equilibrium"),
            ("C08_PHYLOGENETIC_CORRECTION", "phylogeny", "species tree, local gene tree, branch length, block jackknife"),
            ("C09_ILS_INTROGRESSION", "phylogeny", "multispecies-coalescent/topology sensitivity"),
            ("C10_RECOMBINATION_CHROM_SIZE", "covariates", "map or labeled syntenic proxy, chromosome size, local GC"),
            ("C11_PARALOG_SEG_DUP", "orthology", "reciprocal 1:1 copy control and segmental-duplication exclusion"),
            ("C12_PARTITIONS", "annotation", "neutral, fourfold, coding, and noncoding partitions"),
            ("C13_MULTIPLE_TESTING", "statistics", "Holm primary panels; BH and BY exploratory sensitivity"),
            ("C14_H1_H2_SENSITIVITY", "claim_boundary", "heterozygous H1/H2 W/S states limited to QC sensitivity"),
        ]
        for gate, category, requirement in controls:
            add(panel, gate, category, requirement, "NOT_RUN", "required before interpretation", "NOT_RUN_UPSTREAM_GATE", "UPSTREAM_INPUT_GATE", "pre-registered control")

    add("GLOBAL", "E01_EXECUTION_ENVIRONMENT", "execution", "pinned Guix channel and authorized Slurm", f"Guix {CHANNEL_COMMIT}; profile {os.environ.get('GUIX_PROFILE', 'NOT_RECORDED')}; Slurm job {slurm_job_id}; metadata preflight only", "both required", "PASS", "NONE", "slurm environment plus analysis/guix/channels.scm")
    add("GLOBAL", "E02_BULK_MIRROR", "execution", "verified mirror completion", "quota_visibility_unavailable_fail_closed; 0 jobs; 0 verified objects", "complete verified inputs", "FAIL", "UPSTREAM_MIRROR_BLOCKED", "vgp_freeze1_mirror_summary.json")
    return rows


def result_rows(slurm_job_id: str) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []

    def add(panel: str, scope: str, partition: str, estimand: str, label: str,
            model: str, sensitivity: str, candidate_only: str = "false") -> None:
        rows.append({
            "panel_id": panel,
            "scope": scope,
            "partition": partition,
            "estimand_id": estimand,
            "estimand_label": label,
            "status": NOT_ESTIMABLE,
            "callable_bases": NOT_ESTIMABLE,
            "ws_opportunities": NOT_ESTIMABLE,
            "sw_opportunities": NOT_ESTIMABLE,
            "n_ws": NOT_ESTIMABLE,
            "n_sw": NOT_ESTIMABLE,
            "estimate_name": label,
            "estimate_value": NOT_ESTIMABLE,
            "interval_lower": NOT_ESTIMABLE,
            "interval_upper": NOT_ESTIMABLE,
            "uncertainty_method": "NOT_RUN; chromosome/block bootstrap plus parametric bootstrap required",
            "model": model,
            "control_or_sensitivity": sensitivity,
            "multiple_testing": "NOT_APPLIED_NO_TEST; Holm primary, BH/BY exploratory required",
            "historical_only": "true",
            "direct_estimand": "NOT_MEASURED",
            "population_estimand": "NOT_MEASURED",
            "non_allelic_estimand": "NOT_MEASURED",
            "candidate_only": candidate_only,
            "not_measured_reason": "NO_VERIFIED_SEQUENCE;NO_CALLABLE_ORTHOLOGY;MINIMUM_GATES_FAILED",
            "preflight_slurm_job_id": slurm_job_id,
        })

    for panel in PANELS:
        for partition in ("neutral", "fourfold", "coding", "noncoding"):
            add(panel, "genome_wide", partition, "H_SUB", "WS_SW_log_rate_ratio", "context-dependent nonstationary substitution model", "primary plus alternate outgroup/polarization/alignment/mutation/topology sensitivities")
            add(panel, "genome_wide", partition, "H_GBGC", "historical_gBGC_like_fixation_bias", "mutation-selection-like nonstationary model; never population B", "zero-bias parametric bootstrap plus CpG/topology/branch-length sensitivity")
        add(panel, "tract_like", "neutral_noncoding", "H_CLUSTER", "candidate_historical_WS_cluster_excess", "inhomogeneous Poisson/renewal null preserving context, gaps, chromosome, local rate", "label shuffle, multinucleotide/recurrent-mutation null, maxT scales", "true")
        add(panel, "fragmentation_sensitivity", "neutral", "H_SUB", "WS_SW_log_rate_ratio", "same primary model after reason-coded assembly/gap downsampling", "fragmentation and missingness envelope")

    for branch, estimand, reason in (
        ("DIRECT", "direct_event_or_pedigree_tract_rate", "NO_PEDIGREE_OR_GAMETE_TRANSMISSION_IN_PHYLOGENETIC_BRANCH"),
        ("POPULATION", "population_frequency_B", "NO_MULTI_INDIVIDUAL_POPULATION_FREQUENCY_DATA_IN_PHYLOGENETIC_BRANCH"),
        ("NON_ALLELIC", "paralog_homogenization_or_conversion", "PARALOGS_AND_SEGMENTAL_DUPLICATIONS_EXCLUDED;NON_ALLELIC_BRANCH_NOT_RUN"),
    ):
        row = {
            "panel_id": "GLOBAL",
            "scope": "cross_branch_boundary",
            "partition": "NOT_APPLICABLE",
            "estimand_id": branch,
            "estimand_label": estimand,
            "status": "NOT_MEASURED",
            "callable_bases": "NOT_APPLICABLE",
            "ws_opportunities": "NOT_APPLICABLE",
            "sw_opportunities": "NOT_APPLICABLE",
            "n_ws": "NOT_APPLICABLE",
            "n_sw": "NOT_APPLICABLE",
            "estimate_name": estimand,
            "estimate_value": "NOT_MEASURED",
            "interval_lower": "NOT_APPLICABLE",
            "interval_upper": "NOT_APPLICABLE",
            "uncertainty_method": "NOT_APPLICABLE",
            "model": "NOT_RUN_DESIGN_BOUNDARY",
            "control_or_sensitivity": "H1/H2 W/S states cannot be direct evidence",
            "multiple_testing": "NOT_APPLICABLE",
            "historical_only": "false",
            "direct_estimand": "NOT_MEASURED",
            "population_estimand": "NOT_MEASURED",
            "non_allelic_estimand": "NOT_MEASURED",
            "candidate_only": "false",
            "not_measured_reason": reason,
            "preflight_slurm_job_id": slurm_job_id,
        }
        rows.append(row)
    return rows


def report_text(audit: dict[str, object], slurm_job_id: str) -> str:
    return f"""# VGP phylogenetic gBGC pilot: fail-closed execution report

**Disposition:** `NOT_EXECUTED_INPUT_GATE`
**Preflight:** pinned GNU Guix channel `{CHANNEL_COMMIT}`, immutable profile `{os.environ.get('GUIX_PROFILE', 'NOT_RECORDED')}`, authorized Slurm job `{slurm_job_id}`
**Frozen VGP catalog:** commit `{CATALOG_COMMIT}`, SHA-256 `{CATALOG_SHA256}`

## Outcome

Neither H01 nor H02 produced a biological estimate. The upstream Freeze 1 mirror stopped at
`quota_visibility_unavailable_fail_closed`: zero bulk Slurm jobs launched, zero objects transferred,
and zero objects verified. H01 has released UCSC FASTA/2bit inventory entries for 3/5 exact
assemblies; H02 has entries for 1/5. All eight represented sequence objects remain `planned` with
zero observed bytes. The remaining exact close relatives or unreleased catalog entries have no
sequence object in the frozen released-source inventory. No sequence digest, sequence dictionary,
coordinate-compatible annotation, assembly-specific license audit, versioned/checksummed topology,
reciprocal 1:1 orthology product, alignment, callable base, or branch-polarized substitution was
therefore available.

Running phastBias or another nonstationary model against missing or floating inputs would violate
the design. The Slurm job performed metadata-only preflight generation and fail-closed validation;
it did not download sequence, align genomes, infer ancestors, shuffle biological labels, simulate a
fitted biological null, or estimate WS/SW asymmetry, historical bias, or clusters.

## Exact clades and polarization design

| Panel | Ingroup | Outgroup paths | Frozen identity result |
|---|---|---|---|
| H01 | *Spinachia spinachia* focal; *Pungitius pungitius*; *Gasterosteus aculeatus* | *Syngnathus acus* and *S. typhle* | 5 accessions frozen in the design/catalog; only 3 represented in released source inventory; 0 local verified |
| H02 | *Falco naumanni* focal; *F. tinnunculus*; *F. peregrinus* | *F. cherrug* and *F. punctatus* | 3 VGP catalog assemblies plus 2 exact close relatives frozen in design; only 1 represented in released source inventory; 0 local verified |

The row-level manifest records exact accession/version, catalog membership (including the H01
*S. acus* other-haplotype selection), source URL, metadata assembly span/fragmentation, annotation
release, license terms, inventory paths, checksum state, coordinate evidence, and exclusion reason.
The source inventory SHA-256 is `{audit['source_sha256']}`; the mirror manifest SHA-256 is
`{audit['mirror_sha256']}`; the mirror summary SHA-256 is `{audit['summary_sha256']}`.

## Pre-registered gates and missingness accounting

The minimum is three callable ingroup species plus two callable outgroups. H_SUB additionally needs
10,000 callable neutral aligned bases and at least 20 WS plus 20 SW substitutions per interpreted
branch/partition. H_GBGC needs at least 100 directionally informative substitutions. H_CLUSTER needs
at least 50 high-posterior branch-polarized WS substitutions and 20 Mb of continuous callable
single-copy alignment. Callable input is zero because no verified sequence was available;
biological substitution counts are not estimable, so all three gates fail before modeling.

Assembly span and contig count are retained as design metadata, but alignment missingness,
fragmentation loss, callable orthology, and excluded bases cannot be calculated without sequence.
They are reported as `NOT_ESTIMABLE_NO_ALIGNMENT`, never as zero missingness. Exclusion is accounted
for by reason: `NO_VERIFIED_SEQUENCE_PAYLOAD`, `MISSING_RELEASED_SEQUENCE_SOURCE`,
`CHECKSUM_UNRESOLVED`, `ANNOTATION_DICTIONARY_UNAVAILABLE`, `NO_ORTHOLOGY_PRODUCT`, and
`NO_CALLABLE_ALIGNMENT`.

## Models and controls not run

All planned controls remain explicit `NOT_RUN_UPSTREAM_GATE` rows: negative controls, branch/taxon
label shuffling, context/opportunity-matched null simulation, alternate outgroups, ancestral-state
posterior and polarization-error sensitivity, alignment filters, 3-mer/5-mer mutation-bias models,
CpG sensitivity, GC equilibrium and multiple hits, ILS/introgression and local-gene-tree sensitivity,
branch length and phylogenetic correction, chromosome block resampling, recombination proxies,
chromosome size, local GC, paralogy/segmental-duplication exclusion, and Holm/BH/BY correction.
Neutral, fourfold, coding, and noncoding partitions are separately ledgered. No coding partition is
treated as neutral. No native annotation or liftover is admitted until its exact sequence dictionary
and coordinate mapping validate.

## Historical results and uncertainty

Genome-wide H_SUB and H_GBGC rows, tract-like H_CLUSTER rows, and fragmentation/missingness
sensitivity rows are present for both panels. Their counts, opportunities, estimates, intervals, and
callable denominators are `NOT_ESTIMABLE_INPUT_GATE`. This is not a null result and supplies no
evidence for or against historical gBGC. Candidate clusters are not conversion tracts even when a
future run can estimate them.

This branch can only estimate long-term, model-dependent substitution signatures. It does not
observe direct conversion events or biased transmission, does not estimate a pedigree tract rate,
and does not substitute for multi-individual population-frequency evidence.
Direct, population, and non-allelic estimands are explicitly `NOT_MEASURED`.
H1/H2 heterozygous WS/SW states were not used as direct evidence and could only enter a future
QC/sensitivity analysis.

## Safe activation rule

Do not rerun this blocked-output generator as a biological analysis. First obtain enforceable quota
evidence and an authorized selective acquisition amendment covering all ten exact sequences and
needed annotations; verify provider digests plus local SHA-256; bind native annotation or a separately
validated liftover to each exact sequence dictionary; build reciprocal 1:1 single-copy syntenic
alignments; account for excluded bases by reason; then submit a separately fingerprinted Slurm
analysis under a pinned environment that includes the declared nonstationary implementation. Any
accession substitution, outgroup change, or annotation remap requires a pre-result amendment.
"""


def _rectangular_nonempty(path: Path, errors: list[str]) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file():
        errors.append(f"missing artifact: {path}")
        return [], []
    header, rows = read_tsv(path)
    if len(header) != len(set(header)):
        errors.append(f"{path}: duplicate columns")
    if not rows:
        errors.append(f"{path}: no rows")
    for line, row in enumerate(rows, start=2):
        if None in row:
            errors.append(f"{path}:{line}: extra fields")
        for field in header:
            value = row.get(field, "")
            if not value.strip():
                errors.append(f"{path}:{line}: empty {field}")
            if value != value.strip():
                errors.append(f"{path}:{line}: outer whitespace {field}")
    return header, rows


def validate_artifacts(output_dir: Path = ANALYSIS) -> list[str]:
    errors: list[str] = []
    try:
        audit_upstream()
    except (KeyError, OSError, ValueError, json.JSONDecodeError) as error:
        errors.append(f"upstream audit: {error}")
    _, clades = _rectangular_nonempty(output_dir / CLADE_OUTPUT, errors)
    _, qc = _rectangular_nonempty(output_dir / QC_OUTPUT, errors)
    _, results = _rectangular_nonempty(output_dir / RESULT_OUTPUT, errors)
    report_path = output_dir / REPORT_OUTPUT
    if not report_path.is_file():
        errors.append(f"missing artifact: {report_path}")
        report = ""
    else:
        report = report_path.read_text(encoding="utf-8")

    if {row.get("exact_assembly_accession") for row in clades} != expected_accessions():
        errors.append("clade manifest exact accessions differ from H01/H02 freeze")
    if len(clades) != 10 or {row.get("panel_id") for row in clades} != {"H01", "H02"}:
        errors.append("clade manifest must contain ten rows across H01/H02")
    for row in clades:
        if row.get("observed_sequence_bytes") != "0":
            errors.append(f"{row.get('exact_assembly_accession')}: nonzero bytes in blocked run")
        if row.get("local_sequence_sha256") != UNRESOLVED_CHECKSUM:
            errors.append(f"{row.get('exact_assembly_accession')}: invented local digest")
        if not row.get("exact_assembly_accession", "").startswith("GCA_"):
            errors.append("unversioned or invalid assembly accession")

    required_qc = {
        "G01_DESIGN_IDENTITY", "G02_RELEASED_INVENTORY", "G03_LOCAL_SEQUENCE",
        "G04_SEQUENCE_CHECKSUM", "G04B_LICENSE_AUDIT", "G05_ANNOTATION_COORDINATES",
        "G06_TAXON_POLARIZATION", "G06B_TOPOLOGY_IDENTITY",
        "G07_NEUTRAL_CALLABILITY", "G08_H_SUB_COUNTS", "G09_H_GBGC_COUNTS",
        "G10_H_CLUSTER_CALLABILITY", "G11_FRAGMENTATION_MISSINGNESS",
        *(f"C{i:02d}_{name}" for i, name in [
            (1, "LABEL_SHUFFLE"), (2, "NULL_SIMULATION"), (3, "ALTERNATE_OUTGROUP"),
            (4, "ANCESTRAL_UNCERTAINTY"), (5, "ALIGNMENT_FILTERS"), (6, "MUTATION_BIAS"),
            (7, "MULTIPLE_HITS_GC_EQ"), (8, "PHYLOGENETIC_CORRECTION"),
            (9, "ILS_INTROGRESSION"), (10, "RECOMBINATION_CHROM_SIZE"),
            (11, "PARALOG_SEG_DUP"), (12, "PARTITIONS"), (13, "MULTIPLE_TESTING"),
            (14, "H1_H2_SENSITIVITY"),
        ])
    }
    for panel in PANELS:
        ids = {row.get("qc_id") for row in qc if row.get("panel_id") == panel}
        missing = required_qc - ids
        if missing:
            errors.append(f"{panel}: missing QC/control rows {sorted(missing)}")

    historical = [row for row in results if row.get("panel_id") in PANELS]
    if any(row.get("status") != NOT_ESTIMABLE for row in historical):
        errors.append("blocked historical row promoted to a biological result")
    if {row.get("partition") for row in historical} < {"neutral", "fourfold", "coding", "noncoding", "neutral_noncoding"}:
        errors.append("required analysis partitions absent")
    for panel in PANELS:
        estimands = {row.get("estimand_id") for row in historical if row.get("panel_id") == panel}
        if estimands != {"H_SUB", "H_GBGC", "H_CLUSTER"}:
            errors.append(f"{panel}: historical estimand rows incomplete")
    boundaries = {row.get("estimand_id"): row for row in results if row.get("panel_id") == "GLOBAL"}
    if set(boundaries) != {"DIRECT", "POPULATION", "NON_ALLELIC"}:
        errors.append("cross-branch NOT_MEASURED rows incomplete")
    elif any(row.get("status") != "NOT_MEASURED" for row in boundaries.values()):
        errors.append("cross-branch estimand promoted")

    phrases = (
        "NOT_EXECUTED_INPUT_GATE", "quota_visibility_unavailable_fail_closed",
        "negative controls", "label shuffling", "alternate outgroups",
        "ancestral-state", "mutation-bias", "multiple hits", "phylogenetic correction",
        "paralogy/segmental-duplication", "fragmentation/missingness",
        "Direct, population, and non-allelic estimands are explicitly `NOT_MEASURED`",
        "H1/H2 heterozygous WS/SW states were not used as direct evidence",
    )
    for phrase in phrases:
        if phrase not in report:
            errors.append(f"report missing phrase: {phrase}")
    return errors


def write_outputs(output_dir: Path, slurm_job_id: str) -> None:
    audit = audit_upstream()
    output_dir.mkdir(parents=True, exist_ok=True)
    write_tsv(output_dir / CLADE_OUTPUT, clade_rows(audit, slurm_job_id))
    write_tsv(output_dir / QC_OUTPUT, qc_rows(audit, slurm_job_id))
    write_tsv(output_dir / RESULT_OUTPUT, result_rows(slurm_job_id))
    (output_dir / REPORT_OUTPUT).write_text(report_text(audit, slurm_job_id), encoding="utf-8")


def production_environment_errors() -> list[str]:
    errors: list[str] = []
    if os.environ.get("VGP_GBGC_PINNED_GUIX") != CHANNEL_COMMIT:
        errors.append("VGP_GBGC_PINNED_GUIX does not name the pinned channel commit")
    profile = os.environ.get("GUIX_PROFILE", "")
    if not re.fullmatch(r"/gnu/store/[0-9a-z]{32}-[^/]+", profile):
        errors.append("GUIX_PROFILE is not an immutable Guix store path")
    job_id = os.environ.get("SLURM_JOB_ID", "")
    if not re.fullmatch(r"[0-9]+", job_id):
        errors.append("SLURM_JOB_ID is absent or invalid; repository writes require Slurm")
    return errors


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true", help="write the four frozen artifacts")
    parser.add_argument("--validate", action="store_true", help="validate existing artifacts")
    parser.add_argument("--output-dir", type=Path, default=ANALYSIS)
    args = parser.parse_args(argv)
    if not args.write and not args.validate:
        parser.error("select --write and/or --validate")
    if args.write:
        environment_errors = production_environment_errors()
        if environment_errors:
            for error in environment_errors:
                print(f"ERROR: {error}", file=sys.stderr)
            return 2
        write_outputs(args.output_dir, os.environ["SLURM_JOB_ID"])
    if args.validate:
        errors = validate_artifacts(args.output_dir)
        if errors:
            for error in errors:
                print(f"ERROR: {error}", file=sys.stderr)
            return 1
        print("VGP phylogenetic gBGC preflight artifacts validate")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
