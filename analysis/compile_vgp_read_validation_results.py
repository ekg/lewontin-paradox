#!/usr/bin/env python3
"""Compile promoted VGP read-validation evidence into the review packet.

This is deliberately a pure post-processing step: it reads checksum-promoted
Slurm products, rejects canonical-root drift, and preserves the paired-method
covariance labels emitted by the biological computation.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import time
from pathlib import Path
from typing import Any, Mapping


CANONICAL_ROOT = Path("/moosefs/erikg/vgp")
MASK_IDS = ("dp5_100", "dp10_60", "dp10_80", "dp15_80", "dp20_80", "dp10_100")
PRIMARY_MASK = "dp10_80"
REPOSITORY_ROOT = Path(__file__).resolve().parents[1]


def load_json(path: Path, *, require_root: bool = True) -> dict[str, Any]:
    value = json.loads(path.read_text())
    if require_root:
        observed = value.get("canonical_vgp_root", value.get("canonical_root"))
        if observed != str(CANONICAL_ROOT):
            raise ValueError(f"canonical VGP root drift in {path}: {observed!r}")
    return value


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def evidence(path: Path) -> dict[str, Any]:
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": sha256(path)}


def atomic_text(path: Path, text: str) -> None:
    temporary = path.with_suffix(path.suffix + ".partial")
    temporary.write_text(text)
    temporary.replace(path)


def atomic_json(path: Path, value: Mapping[str, Any]) -> None:
    atomic_text(path, json.dumps(value, indent=2, sort_keys=True) + "\n")


def verify_output_manifest(run: Path) -> dict[str, Any]:
    manifest = run / "output_manifest.tsv"
    with manifest.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"empty output manifest: {manifest}")
    for row in rows:
        target = run / row["relative_path"]
        observed = sha256(target)
        if observed != row["sha256"]:
            raise ValueError(f"output digest mismatch: {target}")
    return {**evidence(manifest), "verified_objects": len(rows), "all_digests_match": True}


def coverage_summary(path: Path) -> dict[str, Any]:
    reference_bp = covered_bp = 0
    depth_bp = 0.0
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            length = int(row["endpos"]) - int(row["startpos"]) + 1
            reference_bp += length
            covered_bp += int(row["covbases"])
            depth_bp += float(row["meandepth"]) * length
    return {
        "reference_bp": reference_bp,
        "covered_bp_at_least_one": covered_bp,
        "breadth_at_least_one": covered_bp / reference_bp,
        "length_weighted_mean_depth": depth_bp / reference_bp,
        "evidence": evidence(path),
    }


def flagstat_summary(path: Path) -> dict[str, Any]:
    values: dict[str, Any] = {"evidence": evidence(path)}
    patterns = {
        "total_alignments": re.compile(r"^(\d+) \+ \d+ in total"),
        "primary_reads": re.compile(r"^(\d+) \+ \d+ primary$"),
        "primary_mapped_reads": re.compile(r"^(\d+) \+ \d+ primary mapped \(([0-9.]+)%"),
    }
    for line in path.read_text().splitlines():
        for key, pattern in patterns.items():
            match = pattern.match(line)
            if match:
                values[key] = int(match.group(1))
                if key == "primary_mapped_reads":
                    values["primary_mapped_fraction"] = float(match.group(2)) / 100.0
    return values


def f(value: Any, digits: int = 8) -> str:
    if value is None or value == "":
        return "NA"
    if isinstance(value, int):
        return f"{value:,}"
    return f"{float(value):.{digits}g}"


def pct(value: Any, digits: int = 3) -> str:
    if value is None:
        return "NA"
    return f"{100.0 * float(value):.{digits}f}%"


def p07_assessment(
    qv: Mapping[str, Any],
    primary: Mapping[str, Any],
    psmc: Mapping[str, Any],
    illumina_sites: Mapping[str, Any],
    hifi_sites: Mapping[str, Any],
) -> dict[str, Any]:
    """Apply a conservative, explicit rule to the paired validation evidence."""

    pi_ratio = float(primary["pi_ratio_read_over_assembly"])
    lambda_correlation = float(psmc["lambda_pearson_correlation"])
    qv_value = float(qv["qv"])
    illumina_contradiction = float(illumina_sites["concrete_false_positive_lower_bound_fraction"])
    hifi_contradiction = float(hifi_sites["concrete_false_positive_lower_bound_fraction"])
    severe_sequence_error = qv_value < 20.0
    cross_technology_majority_contradiction = (
        illumina_contradiction >= 0.5 and hifi_contradiction >= 0.5
    )
    material_method_discordance = (
        pi_ratio < 0.8 or pi_ratio > 1.25 or lambda_correlation < 0.5
    )
    if severe_sequence_error:
        classification = "concrete_severe_sequence_error_failure"
        downstream_action = (
            "preserve the core artifact and provenance, but do not use its pi or PSMC "
            "as validated quantitative evidence"
        )
    elif cross_technology_majority_contradiction:
        classification = "concrete_haplotype_reconstruction_failure"
        downstream_action = (
            "preserve the core artifact and provenance, but do not use its pi or PSMC "
            "as validated quantitative evidence"
        )
    elif material_method_discordance:
        classification = "material_method_discordance_without_concrete_failure"
        downstream_action = (
            "retain the core result with low validation confidence and carry the paired "
            "pi bounds and PSMC sensitivity into synthesis"
        )
    else:
        classification = "paired_validation_supports_core_with_measured_sensitivity"
        downstream_action = "retain the core result with the reported uncertainty and covariance"
    return {
        "classification": classification,
        "downstream_action": downstream_action,
        "severe_sequence_error_qv_below_20": severe_sequence_error,
        "cross_technology_majority_homozygous_reference_contradiction": (
            cross_technology_majority_contradiction
        ),
        "material_method_discordance": material_method_discordance,
        "observed": {
            "kmer_qv": qv_value,
            "primary_pi_ratio_read_over_assembly": pi_ratio,
            "primary_psmc_lambda_pearson_correlation": lambda_correlation,
            "illumina_concrete_contradiction_fraction": illumina_contradiction,
            "hifi_concrete_contradiction_fraction": hifi_contradiction,
        },
        "decision_rule": (
            "severe sequence error if k-mer QV <20; concrete reconstruction failure only "
            "if a majority of primary assembly SNPs are depth-qualified homozygous-reference "
            "contradictions in both Illumina and HiFi; otherwise pi ratio outside 0.8-1.25 "
            "or PSMC lambda r <0.5 is material method discordance, not a deletion trigger"
        ),
    }


def qv_effect(assessment: Mapping[str, Any]) -> str:
    observed = assessment["observed"]
    return (
        f"{assessment['classification']}: {assessment['downstream_action']}. "
        f"Primary common-mask read/assembly pi ratio "
        f"{observed['primary_pi_ratio_read_over_assembly']:.6g}; k-mer QV "
        f"{observed['kmer_qv']:.6g}; PSMC lambda r "
        f"{observed['primary_psmc_lambda_pearson_correlation']:.6g}; concrete "
        f"Illumina/HiFi contradiction fractions "
        f"{observed['illumina_concrete_contradiction_fraction']:.6g}/"
        f"{observed['hifi_concrete_contradiction_fraction']:.6g}."
    )


def build_packet(args: argparse.Namespace) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    p07_run = args.p07_run.resolve()
    p09_run = args.p09_run.resolve()
    p07_manifest = verify_output_manifest(p07_run)
    p09_manifest = verify_output_manifest(p09_run)
    acquisition = load_json(args.acquisition_manifest)
    environment = load_json(args.environment_manifest)
    p07_input = load_json(p07_run / "input_manifest.json")
    p07_execution = load_json(p07_run / "execution.json")
    p09_input = load_json(p09_run / "input_manifest.json")
    p09_limits = load_json(p09_run / "validation_limits.json")
    p09_mapping = load_json(args.p09_refined)
    depth = load_json(p07_run / "masks/depth_mask_summary.json")
    kmer_qv = load_json(p07_run / "kmer/H1.qv.json", require_root=False)
    kmer_het = load_json(p07_run / "kmer/heterozygosity.json", require_root=False)
    illumina_sites = load_json(p07_run / "sites/illumina.assembly_evidence.json")
    hifi_sites = load_json(p07_run / "sites/hifi.assembly_evidence.json")
    mapping_qc = {
        technology: {
            "coverage": coverage_summary(p07_run / f"{technology}.coverage.tsv"),
            "flagstat": flagstat_summary(p07_run / f"{technology}.flagstat.txt"),
            "stats_evidence": evidence(p07_run / f"{technology}.stats.txt"),
        }
        for technology in ("illumina", "hifi")
    }

    mask_rows: list[dict[str, Any]] = []
    mask_results: dict[str, Any] = {}
    for mask_id in MASK_IDS:
        comparison_path = p07_run / f"masks/{mask_id}/comparison.json"
        psmc_path = p07_run / f"masks/{mask_id}/psmc_comparison.json"
        comparison = load_json(comparison_path)
        psmc = load_json(psmc_path)
        if comparison.get("method_covariance") != "paired_shared_reference_and_callable_mask":
            raise ValueError(f"lost paired pi covariance label: {comparison_path}")
        if psmc.get("independence_status") != "not_independent_replication":
            raise ValueError(f"lost PSMC covariance label: {psmc_path}")
        row = {
            "canonical_vgp_root": str(CANONICAL_ROOT),
            "selection_id": "P07",
            "mask_id": mask_id,
            "minimum_depth": depth["masks"][mask_id]["minimum_depth_inclusive"],
            "maximum_depth": depth["masks"][mask_id]["maximum_depth_inclusive"],
            **comparison,
            **{f"psmc_{key}": value for key, value in psmc.items() if key not in {"schema_version", "canonical_vgp_root", "selection_id", "mask_id"}},
        }
        mask_rows.append(row)
        mask_results[mask_id] = {
            "comparison": comparison,
            "psmc": psmc,
            "evidence": {
                "comparison": evidence(comparison_path),
                "psmc": evidence(psmc_path),
            },
        }

    primary = mask_results[PRIMARY_MASK]["comparison"]
    primary_psmc = mask_results[PRIMARY_MASK]["psmc"]
    core_pi = 574122 / 267379237
    excluded_primary_bp = 267379237 - primary["callable_bp_common_mask"]
    excluded_primary_sites = 574122 - primary["assembly_sites"]
    excluded_primary_density = excluded_primary_sites / excluded_primary_bp
    primary_exclusion = {
        "core_callable_bp": 267379237,
        "core_assembly_sites": 574122,
        "primary_common_callable_bp": primary["callable_bp_common_mask"],
        "primary_common_assembly_sites": primary["assembly_sites"],
        "primary_callable_bp_fraction_of_core": primary["callable_bp_common_mask"] / 267379237,
        "primary_assembly_site_fraction_of_core": primary["assembly_sites"] / 574122,
        "excluded_bp": excluded_primary_bp,
        "excluded_assembly_sites": excluded_primary_sites,
        "excluded_assembly_site_density": excluded_primary_density,
        "excluded_to_primary_assembly_site_density_ratio": (
            excluded_primary_density / primary["assembly_pi_common_mask"]
        ),
        "primary_assembly_pi_to_core_pi_ratio": primary["assembly_pi_common_mask"] / core_pi,
        "primary_read_pi_to_core_pi_ratio": primary["read_pi_common_mask"] / core_pi,
        "interpretation": (
            "depth-mask exclusion is a paired mappability/collapse sensitivity diagnostic, "
            "not by itself proof that every excluded assembly difference is false"
        ),
    }
    quarantine = Path(
        "/moosefs/erikg/vgp/quarantine/vgp-validation-reads-v1/"
        "2026-07-20T183020Z/"
        "P07_SRR30200290_SRR30200290_2.fastq.gz.partial.content-mismatch-ena_md5_match"
    )
    if not quarantine.is_file():
        raise ValueError(f"missing retained quarantine evidence: {quarantine}")
    invocation_paths = [
        CANONICAL_ROOT / "staging/partials/vgp-validation-reads-v1/P07-R1-invocation-manifest.json",
        CANONICAL_ROOT / "staging/partials/vgp-validation-reads-v1/P07-R2-invocation-manifest.json",
        CANONICAL_ROOT / "staging/partials/vgp-validation-reads-v1/P07-R2-retry-manifest.json",
    ]
    for path in invocation_paths:
        if not path.is_file():
            raise ValueError(f"missing invocation evidence: {path}")

    verified_objects = [
        {
            "object_id": row["object_id"],
            "status": row["status"],
            "bytes": row["observed_bytes"],
            "sha256": row["local_sha256"],
            "view": row["accession_view_path"],
        }
        for row in acquisition["objects"]
        if row["status"] in {"verified", "reused"}
    ]
    assessment = p07_assessment(kmer_qv, primary, primary_psmc, illumina_sites, hifi_sites)
    pair_results = {
        "P04": {
            "species": "Falco naumanni",
            "individual": "bFalNau1",
            "biosample": "SAMN16870685",
            "assembly_result_status": "complete",
            "assembly_pi": 0.004604184795871289,
            "raw_validation_status": "not_estimable_pending_raw_reads",
            "raw_scope": "planned exact-individual CLR run remains unacquired",
            "qv_status": "not_estimable",
            "kmer_heterozygosity_status": "not_estimable",
            "read_pi_status": "not_estimable",
            "psmc_status": "not_estimable",
            "effect": "raw validation pending; the otherwise valid complete core result is retained, not converted to zero or deleted",
        },
        "P07": {
            "species": "Spinachia spinachia",
            "individual": p07_input["exact_individual"],
            "biosample": p07_input["biosample"],
            "assembly_result_status": "complete",
            "assembly_pi": core_pi,
            "assembly_callable_bp": 267379237,
            "assembly_heterozygous_snps": 574122,
            "raw_validation_status": "complete_exact_individual_paired_illumina_plus_hifi",
            "kmer_qv": kmer_qv,
            "kmer_heterozygosity": kmer_het,
            "kmer_heterozygosity_to_core_pi_ratio": (
                kmer_het["heterozygosity_per_base"] / core_pi
            ),
            "depth_and_collapse_diagnostics": depth,
            "primary_mask_exclusion_diagnostic": primary_exclusion,
            "mask_sensitivity": mask_results,
            "primary_mask": PRIMARY_MASK,
            "primary_common_mask_comparison": primary,
            "primary_psmc_comparison": primary_psmc,
            "read_backed_assembly_sites": {
                "illumina": illumina_sites,
                "hifi": hifi_sites,
            },
            "mapping_qc": mapping_qc,
            "validation_assessment": assessment,
            "effect": qv_effect(assessment),
            "method_covariance": {
                "pi": "paired shared H1 coordinates and exact common callable mask; call sets are not independent",
                "psmc": "same H1 reference, PSMC parameters, structural mask, and overlapping response reads; trajectory comparison is sensitivity, not replication",
                "qv": "read k-mers and assembly derive from the same individual and partially overlapping underlying molecules",
            },
            "input_provenance": p07_input,
            "input_manifest_evidence": evidence(p07_run / "input_manifest.json"),
            "execution": p07_execution,
            "promoted_output_manifest": p07_manifest,
        },
        "P09": {
            "species": "Heterodontus francisci",
            "individual": p09_input["exact_individual"],
            "biosample": p09_input["biosample"],
            "assembly_result_status": "incomplete_at_validation_freeze",
            "raw_validation_status": "low_coverage_mapping_control_only",
            "mapping": p09_mapping,
            "mapping_evidence": evidence(args.p09_refined),
            "input_provenance": p09_input,
            "input_manifest_evidence": evidence(p09_run / "input_manifest.json"),
            "limits": p09_limits,
            "qv_status": "not_estimable",
            "kmer_heterozygosity_status": "not_estimable",
        "read_pi_status": "not_estimable",
            "psmc_status": "not_estimable",
            "effect": "0.81x diploid-equivalent partial HiFi coverage establishes individual/reference compatibility only; no absent estimate is represented as zero",
            "promoted_output_manifest": p09_manifest,
        },
    }
    packet = {
        "schema_version": "vgp-read-validation-results-v1",
        "task_id": "validate-vgp-pilot-reads",
        "canonical_vgp_root": str(CANONICAL_ROOT),
        "completed_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "interpretation_rule": "validation assigns confidence and sensitivity; it changes a core result only after a concrete pair-identity, reconstruction, or severe sequence-error failure",
        "environment": {
            "manifest": evidence(args.environment_manifest),
            "frozen_channel_commit": environment["frozen_channel_commit"],
            "profile_store_path": environment["profile"],
            "closure_sha256": environment["closure_sha256"],
        },
        "raw_read_acquisition": {
            "manifest": evidence(args.acquisition_manifest),
            "summary": acquisition["summary"],
            "verified_objects": verified_objects,
            "invocation_manifests": [evidence(path) for path in invocation_paths],
            "retained_failed_R2_transfer": {
                **evidence(quarantine),
                "observed_sha256_from_invocation": "defaed9e929d8acf9d58006a3be51c26dd7c1937079f7d47e9d818420475c965",
                "classification": "full-size content mismatch; ENA MD5 failed and gzip trailer was invalid; clean retry passed",
            },
        },
        "pairs": pair_results,
        "commands_and_telemetry": {
            "p07_executed_worker": evidence(p07_run / "executed_worker.sh"),
            "p07_worker_stdout_snapshot": evidence(p07_run / "telemetry/worker.stdout.snapshot"),
            "p07_worker_stderr_snapshot": evidence(p07_run / "telemetry/worker.stderr.snapshot"),
            "p09_executed_worker": evidence(p09_run / "executed_worker.sh"),
            "slurm_accounting": evidence(args.sacct),
            "slurm_job_attempts": {
                "P09_clean_node_environment_failure": "1787121",
                "P09_split_index_cancelled_after_diagnosis": "1787122",
                "P09_corrected_single_index_complete": "1787124",
                "P07_pending_submission_cancelled_without_compute": "1787557",
                "P07_first_full_run_failed_after_biological_compute_on_psmc_cli": "1787561",
                "P07_optimization_restart_terminated_before_biological_compute": "1787663",
                "P07_reused_status_preflight_failure_before_biological_compute": "1787664",
                "P07_full_validation_complete": p07_execution["slurm_job_id"],
            },
            "repository_commands": [
                "PYTHONPATH=. python3 analysis/acquire_vgp_validation_reads.py --only P07 --delay-seconds 0",
                "PYTHONPATH=. python3 analysis/acquire_vgp_validation_reads.py --verify-only",
                "sbatch analysis/slurm/vgp_read_validation/P09_low_coverage.sh  # 1787121, 1787122, corrected 1787124",
                "sbatch --partition=highmem --nodelist=octopus02 analysis/slurm/vgp_read_validation/P07_validate.sh  # failed 1787561; cancelled 1787663; preflight-failed 1787664; corrected 1787665",
                "sacct -X -j 1787121,1787122,1787124,1787557,1787561,1787663,1787664,1787665 --parsable2 --format=JobIDRaw,JobName,Partition,State,Elapsed,ElapsedRaw,Timelimit,AllocCPUS,ReqMem,MaxRSS,CPUTimeRAW,TotalCPU,ExitCode,NodeList,Start,End,Submit",
            ],
        },
    }
    return packet, mask_rows


def write_mask_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = (
        "canonical_vgp_root", "selection_id", "mask_id", "minimum_depth", "maximum_depth",
        "callable_bp_common_mask", "assembly_sites", "read_sites", "shared_sites",
        "assembly_only_sites", "read_only_sites", "assembly_pi_common_mask",
        "read_pi_common_mask", "pi_difference_read_minus_assembly", "pi_ratio_read_over_assembly",
        "concordant_pi_lower_bracket", "union_pi_upper_bracket", "jaccard",
        "strong_homozygous_alt_snp_sites", "mapping_consensus_qv",
        "psmc_assembly_theta_0", "psmc_read_theta_0", "psmc_theta_ratio_read_over_assembly",
        "psmc_lambda_pearson_correlation", "psmc_log_lambda_rmse",
    )
    temporary = path.with_suffix(path.suffix + ".partial")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def write_pair_tsv(path: Path, packet: Mapping[str, Any]) -> None:
    fields = (
        "canonical_vgp_root", "selection_id", "species", "individual", "biosample",
        "assembly_result_status", "raw_validation_status", "assembly_pi", "primary_callable_bp",
        "primary_assembly_pi", "primary_read_pi", "primary_pi_ratio",
        "primary_assembly_pi_to_core_pi_ratio", "primary_read_pi_to_core_pi_ratio", "kmer_qv",
        "kmer_heterozygosity", "kmer_heterozygosity_to_core_pi_ratio",
        "primary_psmc_theta_ratio", "primary_psmc_lambda_correlation", "effect",
    )
    rows = []
    for selection_id, pair in packet["pairs"].items():
        comparison = pair.get("primary_common_mask_comparison", {})
        psmc = pair.get("primary_psmc_comparison", {})
        rows.append({
            "canonical_vgp_root": str(CANONICAL_ROOT), "selection_id": selection_id,
            "species": pair["species"], "individual": pair["individual"], "biosample": pair["biosample"],
            "assembly_result_status": pair["assembly_result_status"],
            "raw_validation_status": pair["raw_validation_status"], "assembly_pi": pair.get("assembly_pi", ""),
            "primary_callable_bp": comparison.get("callable_bp_common_mask", ""),
            "primary_assembly_pi": comparison.get("assembly_pi_common_mask", ""),
            "primary_read_pi": comparison.get("read_pi_common_mask", ""),
            "primary_pi_ratio": comparison.get("pi_ratio_read_over_assembly", ""),
            "primary_assembly_pi_to_core_pi_ratio": pair.get("primary_mask_exclusion_diagnostic", {}).get("primary_assembly_pi_to_core_pi_ratio", ""),
            "primary_read_pi_to_core_pi_ratio": pair.get("primary_mask_exclusion_diagnostic", {}).get("primary_read_pi_to_core_pi_ratio", ""),
            "kmer_qv": pair.get("kmer_qv", {}).get("qv", ""),
            "kmer_heterozygosity": pair.get("kmer_heterozygosity", {}).get("heterozygosity_per_base", ""),
            "kmer_heterozygosity_to_core_pi_ratio": pair.get("kmer_heterozygosity_to_core_pi_ratio", ""),
            "primary_psmc_theta_ratio": psmc.get("theta_ratio_read_over_assembly", ""),
            "primary_psmc_lambda_correlation": psmc.get("lambda_pearson_correlation", ""),
            "effect": pair["effect"],
        })
    temporary = path.with_suffix(path.suffix + ".partial")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def render_report(packet: Mapping[str, Any], mask_rows: list[dict[str, Any]]) -> str:
    p07 = packet["pairs"]["P07"]
    p09 = packet["pairs"]["P09"]
    primary = p07["primary_common_mask_comparison"]
    psmc = p07["primary_psmc_comparison"]
    depth = p07["depth_and_collapse_diagnostics"]
    qv = p07["kmer_qv"]
    het = p07["kmer_heterozygosity"]
    illumina = p07["read_backed_assembly_sites"]["illumina"]
    hifi = p07["read_backed_assembly_sites"]["hifi"]
    mapping_qc = p07["mapping_qc"]
    assessment = p07["validation_assessment"]
    exclusion = p07["primary_mask_exclusion_diagnostic"]
    mask_by_id = {row["mask_id"]: row for row in mask_rows}
    lines = [
        "# VGP pilot raw-read validation v1",
        "",
        f"Date completed: {packet['completed_utc'][:10]} UTC  ",
        "WG task: `validate-vgp-pilot-reads`  ",
        f"Canonical shared VGP root: `{CANONICAL_ROOT}`",
        "",
        "## Decision",
        "",
        f"P07 is validated as an exact-individual paired-method sensitivity result and classified `{assessment['classification']}`. The downstream action is to {assessment['downstream_action']}. This preserves the original artifact and provenance even when validation changes its evidentiary use. P04 remains a valid completed core result with raw validation pending. P09 supplies only a low-coverage compatibility control and no absent estimate is converted to zero.",
        "",
        "| Pair | Raw evidence outcome | Assembly pi | Read/pi comparison | PSMC effect | Decision |",
        "|---|---|---:|---:|---|---|",
        f"| P04 | Exact CLR run still pending | {packet['pairs']['P04']['assembly_pi']:.10g} | Not estimable | Not estimable | Retain core; validation pending |",
        f"| P07 | Exact BioSample Illumina + HiFi complete | {p07['assembly_pi']:.10g} | primary common-mask ratio {primary['pi_ratio_read_over_assembly']:.6g} | theta ratio {psmc['theta_ratio_read_over_assembly']:.6g}; lambda r {psmc['lambda_pearson_correlation']:.6g} | {assessment['classification']} |",
        f"| P09 | One of seven HiFi cells; {p09['mapping']['nominal_diploid_equivalent_coverage']:.3f}x diploid-equivalent | Incomplete at freeze | Not estimable | Not estimable | Compatibility only |",
        "",
        "## P07 common-mask pi concordance",
        "",
        "The inherited final denominator was independently reconstructed as 267,379,237 bp after the exact non-SNP flank subtraction; 574,122 SNPs reproduce pi = 0.0021472198306856562. Every comparison below restricts both callsets to the same assembly coordinates and the same read-depth mask. Thus differences are paired sensitivity estimates, not independent replication.",
        f"The primary depth mask retains {pct(exclusion['primary_callable_bp_fraction_of_core'])} of the inherited callable bases but only {pct(exclusion['primary_assembly_site_fraction_of_core'])} of its assembly differences. The excluded {exclusion['excluded_bp']:,} bp contain {exclusion['excluded_assembly_sites']:,} assembly differences, an assembly-difference density {f(exclusion['excluded_to_primary_assembly_site_density_ratio'], 6)} times that inside the primary common mask. Accordingly, primary common-mask assembly/read pi are {f(exclusion['primary_assembly_pi_to_core_pi_ratio'], 6)}/{f(exclusion['primary_read_pi_to_core_pi_ratio'], 6)} times the inherited core pi. This is a strong mappability/collapse sensitivity signal, not an assertion that every excluded difference is false.",
        "",
        "| Mask (inclusive DP) | Common bp | Assembly pi | Read pi | Read/assembly | Shared | Assembly-only | Read-only | Jaccard |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in mask_rows:
        lines.append(
            f"| {row['mask_id']} ({row['minimum_depth']}-{row['maximum_depth']}) | "
            f"{row['callable_bp_common_mask']:,} | {row['assembly_pi_common_mask']:.8g} | "
            f"{row['read_pi_common_mask']:.8g} | {row['pi_ratio_read_over_assembly']:.6g} | "
            f"{row['shared_sites']:,} | {row['assembly_only_sites']:,} | {row['read_only_sites']:,} | {row['jaccard']:.6g} |"
        )
    lines += [
        "",
        f"Changing only the upper cutoff from DP80 to DP100 changes the read/assembly ratio from {mask_by_id['dp10_80']['pi_ratio_read_over_assembly']:.6g} to {mask_by_id['dp10_100']['pi_ratio_read_over_assembly']:.6g}, whereas raising the lower cutoff from DP5 to DP20 changes it from {mask_by_id['dp5_100']['pi_ratio_read_over_assembly']:.6g} to {mask_by_id['dp20_80']['pi_ratio_read_over_assembly']:.6g}. Low-depth/mappability exclusion, rather than high-depth collapse filtering, is therefore the dominant callable-mask sensitivity in this read set.",
        f"For the predeclared primary DP10-80 mask, the lower concordant bracket is {primary['concordant_pi_lower_bracket']:.8g} and the union upper bracket is {primary['union_pi_upper_bracket']:.8g}. Assembly-only calls give a candidate assembly false-positive upper bound of {pct(primary['candidate_assembly_false_positive_upper_bound_fraction_of_assembly_calls'])}; read-only calls give a candidate assembly false-negative upper bound of {pct(primary['candidate_assembly_false_negative_upper_bound_fraction_of_read_calls'])}. Neither is a proven error rate: the former includes read-caller misses and the latter includes assembly-caller misses. The read heterozygote allele-balance median is {f(primary['read_het_alt_balance']['median'], 6)}; {pct(primary['read_het_alt_balance']['fraction_0_3_to_0_7'])} lie from 0.30 through 0.70.",
        "",
        "## Sequence QV, k-mer heterozygosity, and structural diagnostics",
        "",
        f"Illumina primary mapping is {pct(mapping_qc['illumina']['flagstat']['primary_mapped_fraction'])} with {mapping_qc['illumina']['coverage']['length_weighted_mean_depth']:.3f}x whole-H1 mean depth and {pct(mapping_qc['illumina']['coverage']['breadth_at_least_one'])} breadth. HiFi primary mapping is {pct(mapping_qc['hifi']['flagstat']['primary_mapped_fraction'])} with {mapping_qc['hifi']['coverage']['length_weighted_mean_depth']:.3f}x mean depth and {pct(mapping_qc['hifi']['coverage']['breadth_at_least_one'])} breadth.",
        f"The Illumina 21-mer containment estimate gives H1 QV {f(qv.get('qv'), 7)} (binomial-only 95% interval {f(qv.get('qv_lower_95'), 7)}-{f(qv.get('qv_upper_95'), 7)}), from {qv['assembly_kmer_occurrences']:,} assembly k-mer occurrences and {qv['assembly_kmer_occurrences_below_trusted_read_threshold']:,} below the trusted read threshold. Correlated k-mers, read errors, coverage bias, and assembly/read dependence remain systematic and are not hidden inside that narrow interval.",
        f"The transparent four-component negative-binomial spectrum model estimates k-mer heterozygosity {f(het.get('heterozygosity_per_base'), 8)}, {f(p07['kmer_heterozygosity_to_core_pi_ratio'], 6)} times the inherited core pi, with heterozygous/homozygous peaks at {f(het.get('heterozygous_peak_depth'), 5)}x/{f(het.get('homozygous_peak_depth'), 5)}x and fit R-squared {f(het.get('fit_r_squared'), 6)}. This is explicitly a model-based spectrum estimate, not a substituted published number or a population estimate; its disagreement with stringent mapped-read calling is retained rather than averaged away.",
        f"Across the inherited assembly callable input, the positive-depth mode is {depth['depth_structure']['modal_positive_depth']}x; zero depth is {pct(depth['depth_structure']['zero_depth_fraction'])}, below-half-mode is {pct(depth['depth_structure']['below_half_mode_fraction'])}, above 1.5x is {pct(depth['depth_structure']['above_1_5x_mode_fraction'])}, and above 2x is {pct(depth['depth_structure']['above_2x_mode_fraction'])}. Low-depth excess diagnoses duplication/mappability/dropout and high-depth excess diagnoses collapse/repeats; neither alone proves an assembly error.",
        "",
        "## Direct read support and error bounds",
        "",
        f"At the DP10-80 assembly SNPs, Illumina classifies {illumina['supported_heterozygous']:,} as supported heterozygotes and {illumina['contradicted_homozygous_reference']:,} as depth-qualified homozygous-reference contradictions. The concrete false-positive lower-bound fraction is {pct(illumina['concrete_false_positive_lower_bound_fraction'])} (binomial Wilson 95% interval {pct(illumina['concrete_false_positive_fraction_wilson_95'][0])}-{pct(illumina['concrete_false_positive_fraction_wilson_95'][1])}); counting every ambiguous, out-of-mask, or unobserved site produces the deliberately conservative candidate upper bound {pct(illumina['candidate_false_positive_upper_bound_fraction'])}.",
        f"Independent HiFi pileups classify {hifi['supported_heterozygous']:,} supported and {hifi['contradicted_homozygous_reference']:,} contradicted sites, with concrete lower/upper-candidate fractions {pct(hifi['concrete_false_positive_lower_bound_fraction'])}/{pct(hifi['candidate_false_positive_upper_bound_fraction'])} and a lower-bound Wilson 95% interval {pct(hifi['concrete_false_positive_fraction_wilson_95'][0])}-{pct(hifi['concrete_false_positive_fraction_wilson_95'][1])}. HiFi and Illumina share the individual and reference but differ in chemistry and mapping behavior; their agreement is informative without being fully independent.",
        f"The primary read callset has {primary['strong_homozygous_alt_snp_sites']:,} homozygous-alt SNP discrepancies, corresponding to mapping-consensus QV {f(primary.get('mapping_consensus_qv'), 7)} with scope limited to accessible homozygous-alt SNPs. Structural and inaccessible errors are excluded from that QV.",
        "",
        "## PSMC sensitivity",
        "",
        "| Mask | Theta read/assembly | Lambda Pearson r | log-lambda RMSE | time-grid r |",
        "|---|---:|---:|---:|---:|",
    ]
    for row in mask_rows:
        lines.append(
            f"| {row['mask_id']} | {row['psmc_theta_ratio_read_over_assembly']:.6g} | "
            f"{row['psmc_lambda_pearson_correlation']:.6g} | {row['psmc_log_lambda_rmse']:.6g} | "
            f"{row['psmc_time_2N0_pearson_correlation']:.6g} |"
        )
    lines += [
        "",
        "Only the final completed 64-interval optimization round is compared. Both trajectories use the same H1 coordinates, PSMC parameterization, inherited structural mask, and overlapping biological response; the correlations therefore measure callable/caller sensitivity and must not be treated as an independent demographic replication. Scaling is intentionally not introduced here.",
        "",
        "## P09 low-coverage control",
        "",
        f"The exact P09 cell mapped {p09['mapping']['mapped_queries']:,}/{p09['mapping']['raw_reads_metadata']:,} reads ({pct(p09['mapping']['mapped_query_fraction_of_raw_reads'])}); query-union bases were {pct(p09['mapping']['mapped_query_fraction_of_raw_bases'])} of metadata bases, weighted alignment identity was {pct(p09['mapping']['weighted_alignment_identity'])}, and H1 breadth was {pct(p09['mapping']['reference_breadth_at_least_one_mapping_fraction'])}. Coverage is only {p09['mapping']['nominal_diploid_equivalent_coverage']:.3f}x diploid-equivalent from one of seven cells. QV, k-mer heterozygosity, callable read pi, and PSMC validation are therefore not estimable. These metrics establish compatibility, not callability.",
        "",
        "## Raw evidence, failures, reproducibility, and scope",
        "",
        f"The cumulative acquisition ledger contains {packet['raw_read_acquisition']['summary']['verified']['objects']} verified/reused objects totaling {packet['raw_read_acquisition']['summary']['verified']['bytes']:,} bytes. P07 contributes one HiFi and paired Illumina objects; P09 contributes the retained cell. Every canonical CAS object was rehashed offline. P04's 42,344,746,693-byte exact CLR run remains planned.",
        f"A full-size first P07 R2 transfer failed ENA MD5 and gzip integrity and remains quarantined at `{packet['raw_read_acquisition']['retained_failed_R2_transfer']['path']}` (SHA-256 `{packet['raw_read_acquisition']['retained_failed_R2_transfer']['observed_sha256_from_invocation']}`). A clean independent retry produced the verified canonical SHA-256 `c542f6efd9fc1d8f557c89629743fbe4a39584f24002c84185e479057cb443ac`; the failed payload was never promoted.",
        f"GNU Guix is frozen at channel commit `{packet['environment']['frozen_channel_commit']}` with profile `{packet['environment']['profile_store_path']}` and closure digest `{packet['environment']['closure_sha256']}`. Slurm jobs, failed/cancelled precursors, requested resources, elapsed times, and MaxRSS availability are retained in `analysis/vgp_read_validation_sacct_v1.tsv`; executed workers, tool digests, stdout/stderr snapshots, per-stage GNU time telemetry, input manifests, and output manifests remain below `{CANONICAL_ROOT}/derived/read-validation/runs/`.",
        "",
        f"The result is per individual, not a population mean. Binomial intervals exclude correlated-k-mer and mapping systematics. Read and assembly methods intentionally share the individual, H1 coordinate system, structural mask, and in part the molecular data, so covariance is preserved in every machine comparison. The conservative machine decision rule is: {assessment['decision_rule']}. Validation never deletes the original artifact; only a concrete identity, reconstruction, or severe sequence-error failure changes its downstream quantitative use.",
    ]
    return "\n".join(lines) + "\n"


def write_evidence_manifest(args: argparse.Namespace) -> None:
    """Bind repository summaries to the canonical raw/promoted evidence ledgers."""

    repository_artifacts = [args.results, args.pairs, args.masks, args.report, args.sacct]
    external_ledgers = [
        args.acquisition_manifest,
        args.environment_manifest,
        args.p07_run / "output_manifest.tsv",
        args.p09_run / "output_manifest.tsv",
        args.p09_refined,
    ]
    atomic_json(
        args.evidence_manifest,
        {
            "schema_version": "vgp-read-validation-evidence-manifest-v1",
            "task_id": "validate-vgp-pilot-reads",
            "canonical_vgp_root": str(CANONICAL_ROOT),
            "repository_artifacts": [evidence(path.resolve()) for path in repository_artifacts],
            "canonical_and_external_ledgers": [evidence(path.resolve()) for path in external_ledgers],
            "scope": (
                "SHA-256 binding for final repository summaries and immutable acquisition, "
                "environment, promoted-run, and P09-refinement ledgers; promoted-run ledgers "
                "were independently rehashed before this manifest was emitted"
            ),
        },
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--p07-run", type=Path, required=True)
    parser.add_argument("--p09-run", type=Path, required=True)
    parser.add_argument("--p09-refined", type=Path, required=True)
    parser.add_argument("--acquisition-manifest", type=Path, default=REPOSITORY_ROOT / "analysis/vgp_validation_reads_manifest_v1.json")
    parser.add_argument("--environment-manifest", type=Path, default=REPOSITORY_ROOT / "analysis/vgp_read_validation_environment_v1.json")
    parser.add_argument("--sacct", type=Path, default=REPOSITORY_ROOT / "analysis/vgp_read_validation_sacct_v1.tsv")
    parser.add_argument("--results", type=Path, default=REPOSITORY_ROOT / "analysis/vgp_read_validation_results_v1.json")
    parser.add_argument("--pairs", type=Path, default=REPOSITORY_ROOT / "analysis/vgp_read_validation_per_pair_v1.tsv")
    parser.add_argument("--masks", type=Path, default=REPOSITORY_ROOT / "analysis/vgp_read_validation_mask_sensitivity_v1.tsv")
    parser.add_argument("--report", type=Path, default=REPOSITORY_ROOT / "analysis/vgp_read_validation_report_v1.md")
    parser.add_argument("--evidence-manifest", type=Path, default=REPOSITORY_ROOT / "analysis/vgp_read_validation_evidence_manifest_v1.json")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    packet, mask_rows = build_packet(args)
    atomic_json(args.results, packet)
    write_pair_tsv(args.pairs, packet)
    write_mask_tsv(args.masks, mask_rows)
    atomic_text(args.report, render_report(packet, mask_rows))
    write_evidence_manifest(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
