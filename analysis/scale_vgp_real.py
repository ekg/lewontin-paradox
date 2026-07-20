#!/usr/bin/env python3
"""Reconcile and independently verify the real VGP core scale-out.

The Freeze 1 spreadsheet is a discovery catalog, not an eligibility oracle.
Only the ten pairs in ``vgp_10_pair_manifest.tsv`` have frozen reciprocal,
same-individual provenance and mutually-comparable assembly reports.  This
module therefore accounts for all 716 catalog rows and all 569 catalog links,
while refusing to promote the remaining discovery links to biological pairs.

``VGP_ROOT`` is the single configuration variable for every canonical input,
scratch-promotion target, and derived product.  Repository files are immutable
method/provenance inputs; biological products are always rooted below it.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import os
import re
import subprocess
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Mapping, Sequence


REPO_ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = REPO_ROOT / "analysis"
DEFAULT_VGP_ROOT = Path(os.environ.get("VGP_ROOT", "/moosefs/erikg/vgp"))
CATALOG_NAME = "VGPPhase1-freeze-1.0.commit-dc1b2af5a7741b97d66fb10cb2bce97f41765cdf.tsv"
CATALOG_SHA256 = "9c58420484a8b76a2d6175b7c26bf709e68bdc726a67fc7541b8c2b5a2fc13a4"
CATALOG_COMMIT = "dc1b2af5a7741b97d66fb10cb2bce97f41765cdf"
AUTHORIZATION_ID = "vgp10-auth-20260718-v2"
PAIR_IDS = tuple(f"P{value:02d}" for value in range(1, 11))
ACCESSION_RE = re.compile(r"^GC[AF]_\d{9}\.\d+$")

PAIR_FIELDS = (
    "schema_version", "canonical_vgp_root", "catalog_commit", "catalog_sha256",
    "catalog_line", "catalog_record_ordinal", "catalog_pair_id", "link_class",
    "scientific_name", "taxid", "main_assembly_id", "linked_assembly_id",
    "h1_accession_version", "h2_accession_version", "selection_id", "eligibility",
    "eligibility_reason", "confidence_tier", "annotation_covariate", "qv_covariate",
    "busco_covariate", "raw_read_covariate", "repeat_covariate", "kmer_covariate",
    "execution_disposition", "hard_failure_class", "latest_stage", "latest_job_id",
    "latest_job_state", "attempted_slurm_jobs", "completed_slurm_allocations",
    "callable_bp", "heterozygous_snps", "callable_pi", "primary_theta_0",
    "bootstrap_attempts", "finite_bootstraps", "bootstrap_theta_q025",
    "bootstrap_theta_median", "bootstrap_theta_q975", "psmc_uncertainty_status",
    "scenario_uncertainty_status", "annotation_partition_status",
    "fastga_scratch_contract_status", "fastga_scratch_contract_job_id",
    "fastga_scratch_contract_packet",
    "same_pair_pi_psmc_non_independent", "population_inference_authorized",
    "verification_job_id", "verification_status", "result_packet",
)

CATALOG_FIELDS = (
    "schema_version", "canonical_vgp_root", "catalog_commit", "catalog_sha256",
    "catalog_line", "catalog_record_ordinal", "scientific_name", "taxid",
    "main_accession_version", "catalog_link_count", "eligible_audited_pair_count",
    "catalog_disposition",
)

SACCT_FIELDS = (
    "JobIDRaw", "JobName", "State", "Elapsed", "ElapsedRaw", "Timelimit",
    "AllocCPUS", "ReqMem", "MaxRSS", "CPUTimeRAW", "TotalCPU", "ExitCode",
    "NodeList", "Start", "End", "Submit",
)


class ScaleError(RuntimeError):
    """A closed-world, provenance, result, or scheduler invariant failed."""


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_json(path: Path) -> dict[str, object]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ScaleError(f"expected JSON object: {path}")
    return value


def atomic_text(path: Path, value: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_name(f".{path.name}.partial-{os.getpid()}")
    with partial.open("w", encoding="utf-8", newline="") as handle:
        handle.write(value)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(partial, path)


def atomic_json(path: Path, value: Mapping[str, object]) -> None:
    atomic_text(path, json.dumps(value, indent=2, sort_keys=True) + "\n")


def atomic_tsv(path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, object]]) -> None:
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    atomic_text(path, buffer.getvalue())


def split_values(value: str) -> list[str]:
    return [item.strip() for item in value.split(",") if item.strip()]


def canonical_paths(vgp_root: Path) -> dict[str, Path]:
    root = vgp_root.resolve()
    if str(root) != "/moosefs/erikg/vgp":
        raise ScaleError(f"canonical root must resolve exactly to /moosefs/erikg/vgp: {root}")
    return {
        "root": root,
        "catalog": root / "manifests" / CATALOG_NAME,
        "inputs": root / "pilot" / "inputs",
        "pilot_run": root / "pilot" / "runs" / f"{AUTHORIZATION_ID}-pilot-v1",
        "p07_run": root / "pilot" / "outputs" / AUTHORIZATION_ID / "P07" / "core",
        "scale_root": root / "derived" / "scale-vgp-real-v1",
        "submissions": root / "pilot" / "manifests" / f"{AUTHORIZATION_ID}-pilot-v1.submissions.tsv",
    }


def load_catalog(path: Path) -> list[dict[str, str]]:
    if sha256_file(path) != CATALOG_SHA256:
        raise ScaleError("frozen catalog SHA-256 drift")
    rows = read_tsv(path)
    if len(rows) != 716:
        raise ScaleError(f"frozen catalog row drift: {len(rows)} != 716")
    return rows


def load_roster(path: Path) -> list[dict[str, str]]:
    rows = read_tsv(path)
    if [row["selection_id"] for row in rows] != list(PAIR_IDS):
        raise ScaleError("audited eligible roster must be exactly ordered P01..P10")
    pairs: set[tuple[str, str]] = set()
    for row in rows:
        h1, h2 = row["h1_accession_version"], row["h2_accession_version"]
        if not ACCESSION_RE.fullmatch(h1) or not ACCESSION_RE.fullmatch(h2) or h1 == h2:
            raise ScaleError(f"invalid audited pair accessions: {row['selection_id']}")
        if (h1, h2) in pairs:
            raise ScaleError("duplicate audited H1/H2 tuple")
        pairs.add((h1, h2))
        evidence = row["reciprocal_pair_evidence"].lower()
        if "shared biosample" not in evidence or not any(token in evidence for token in ("reciprocal", "each report")):
            raise ScaleError(f"same-individual reciprocal evidence absent: {row['selection_id']}")
        ratio = int(row["h2_length_bp"]) / int(row["h1_length_bp"])
        if not 0.80 <= ratio <= 1.25:
            raise ScaleError(f"audited assemblies are not mutually comparable: {row['selection_id']}")
    return rows


def links(row: Mapping[str, str]) -> list[tuple[str, int, str, str]]:
    result: list[tuple[str, int, str, str]] = []
    for link_class, id_column, accession_column in (
        ("other_high_quality", "Assembly IDs other high-quality haplotypes", "Accession #s other high-quality haplotypes"),
        ("alternate", "Assembly IDs alternate haplotypes", "Accession #s alternate haplotypes"),
    ):
        identifiers = split_values(row[id_column])
        accessions = split_values(row[accession_column])
        for ordinal, accession in enumerate(accessions, 1):
            linked_id = identifiers[ordinal - 1] if ordinal <= len(identifiers) else "UNRESOLVED_CARDINALITY"
            result.append((link_class, ordinal, linked_id, accession))
    return result


def build_closed_world(
    catalog: Sequence[Mapping[str, str]], roster: Sequence[Mapping[str, str]], vgp_root: Path,
) -> tuple[list[dict[str, object]], list[dict[str, object]], dict[str, str]]:
    roster_by_pair = {(row["h1_accession_version"], row["h2_accession_version"]): row for row in roster}
    matched: dict[str, str] = {}
    pair_rows: list[dict[str, object]] = []
    catalog_rows: list[dict[str, object]] = []
    linked_total = 0
    for ordinal, row in enumerate(catalog, 1):
        catalog_line = ordinal + 1
        h1 = row["Accession # for main haplotype"].strip()
        row_links = links(row)
        eligible_count = 0
        class_counts: Counter[str] = Counter()
        for link_class, _class_ordinal, linked_id, h2 in row_links:
            linked_total += 1
            class_counts[link_class] += 1
            tag = "HQ" if link_class == "other_high_quality" else "ALT"
            catalog_pair_id = f"R{catalog_line:04d}-{tag}{class_counts[link_class]:02d}"
            roster_row = roster_by_pair.get((h1, h2))
            is_self = bool(h1) and h1 == h2
            if roster_row:
                selection_id = roster_row["selection_id"]
                if selection_id in matched:
                    raise ScaleError(f"audited pair matched multiple catalog links: {selection_id}")
                matched[selection_id] = catalog_pair_id
                eligible_count += 1
                eligibility = "eligible_exact_audited"
                reason = "PASS_RECIPROCAL_SAME_BIOSAMPLE_MUTUALLY_COMPARABLE_FROZEN_ACCESSIONS"
                confidence = roster_row["core_confidence_tier"]
            elif is_self:
                selection_id, eligibility, reason, confidence = (
                    "", "excluded", "CATALOG_SELF_LINK_NOT_A_PAIR", "X"
                )
            else:
                selection_id, eligibility, reason, confidence = (
                    "", "not_eligible_discovery_link",
                    "EXACT_SAME_INDIVIDUAL_AND_MUTUAL_COMPARABILITY_NOT_AUDITED", "UNASSIGNED"
                )
            pair_rows.append({
                "schema_version": "vgp-real-scaleout-pair-v1", "canonical_vgp_root": str(vgp_root),
                "catalog_commit": CATALOG_COMMIT, "catalog_sha256": CATALOG_SHA256,
                "catalog_line": catalog_line, "catalog_record_ordinal": ordinal,
                "catalog_pair_id": catalog_pair_id, "link_class": link_class,
                "scientific_name": row["Scientific Name"].strip(), "taxid": row["NCBI taxon ID"].strip(),
                "main_assembly_id": row["Assembly ID main haplotype"].strip(),
                "linked_assembly_id": linked_id, "h1_accession_version": h1,
                "h2_accession_version": h2, "selection_id": selection_id, "eligibility": eligibility,
                "eligibility_reason": reason, "confidence_tier": confidence,
                "annotation_covariate": "NOT_EVALUATED_NONELIGIBLE",
                "qv_covariate": row["QV"].strip() or "UNVALIDATED",
                "busco_covariate": "UNVALIDATED", "raw_read_covariate": "UNVALIDATED",
                "repeat_covariate": "UNVALIDATED", "kmer_covariate": "UNVALIDATED",
                "execution_disposition": "NOT_APPLICABLE_NONELIGIBLE" if not roster_row else "PENDING_RECONCILIATION",
                "hard_failure_class": "", "latest_stage": "", "latest_job_id": "",
                "latest_job_state": "", "attempted_slurm_jobs": 0,
                "completed_slurm_allocations": 0, "callable_bp": "", "heterozygous_snps": "",
                "callable_pi": "", "primary_theta_0": "", "bootstrap_attempts": "",
                "finite_bootstraps": "", "bootstrap_theta_q025": "",
                "bootstrap_theta_median": "", "bootstrap_theta_q975": "",
                "psmc_uncertainty_status": "NOT_APPLICABLE_NONELIGIBLE",
                "scenario_uncertainty_status": "NOT_APPLICABLE_NONELIGIBLE",
                "annotation_partition_status": "NOT_APPLICABLE_NONELIGIBLE",
                "fastga_scratch_contract_status": "NOT_APPLICABLE_NONELIGIBLE",
                "fastga_scratch_contract_job_id": "", "fastga_scratch_contract_packet": "",
                "same_pair_pi_psmc_non_independent": "true" if roster_row else "not_applicable",
                "population_inference_authorized": "false", "verification_job_id": "",
                "verification_status": "", "result_packet": "NOT_APPLICABLE_NONELIGIBLE",
            })
        if not h1:
            disposition = "CATALOG_UNRELEASED_NO_MAIN_ACCESSION"
        elif not row_links:
            disposition = "CATALOG_RELEASED_NO_LINKED_HAPLOTYPE"
        elif eligible_count:
            disposition = "CONTAINS_AUDITED_ELIGIBLE_PAIR"
        else:
            disposition = "LINKS_ACCOUNTED_NONE_AUDITED_ELIGIBLE"
        catalog_rows.append({
            "schema_version": "vgp-real-scaleout-catalog-v1", "canonical_vgp_root": str(vgp_root),
            "catalog_commit": CATALOG_COMMIT, "catalog_sha256": CATALOG_SHA256,
            "catalog_line": catalog_line, "catalog_record_ordinal": ordinal,
            "scientific_name": row["Scientific Name"].strip(), "taxid": row["NCBI taxon ID"].strip(),
            "main_accession_version": h1, "catalog_link_count": len(row_links),
            "eligible_audited_pair_count": eligible_count, "catalog_disposition": disposition,
        })
    if linked_total != 569 or len(pair_rows) != 569 or len(catalog_rows) != 716:
        raise ScaleError("closed-world catalog multiplicity drift")
    if set(matched) != set(PAIR_IDS):
        raise ScaleError(f"audited roster/catalog reconciliation failed: {sorted(matched)}")
    return catalog_rows, pair_rows, matched


def load_submissions(path: Path, require_all_pairs: bool = True) -> list[dict[str, str]]:
    rows = read_tsv(path)
    eligible = [row for row in rows if row.get("selection_id") in PAIR_IDS]
    if not eligible:
        raise ScaleError("zero biological Slurm submissions is forbidden")
    if require_all_pairs and {row["selection_id"] for row in eligible} != set(PAIR_IDS):
        raise ScaleError("not every audited eligible pair has a biological Slurm submission")
    if any(row.get("canonical_vgp_root") != "/moosefs/erikg/vgp" for row in eligible):
        raise ScaleError("pilot submission escaped canonical root")
    return eligible


def require_submission_coverage(rows: Sequence[Mapping[str, str]]) -> None:
    observed = {row["selection_id"] for row in rows if row.get("selection_id") in PAIR_IDS}
    if observed != set(PAIR_IDS):
        raise ScaleError(f"biological Slurm submission coverage is incomplete: {sorted(observed)}")
    if any(row.get("canonical_vgp_root") != "/moosefs/erikg/vgp" for row in rows):
        raise ScaleError("submission manifest escaped canonical root")


def query_sacct(job_ids: Iterable[str]) -> list[dict[str, str]]:
    ids = sorted(set(job_ids), key=lambda value: (int(value.split("_")[0]), value))
    result: list[dict[str, str]] = []
    for start in range(0, len(ids), 80):
        run = subprocess.run(
            ["sacct", "-X", "-n", "-P", "-j", ",".join(ids[start:start + 80]),
             "--format=" + ",".join(SACCT_FIELDS)],
            check=False, text=True, capture_output=True,
        )
        if run.returncode:
            raise ScaleError(f"sacct failed: {run.stderr.strip()}")
        for line in run.stdout.splitlines():
            values = line.split("|")
            if values and values[-1] == "":
                values.pop()
            if len(values) != len(SACCT_FIELDS):
                raise ScaleError("unexpected sacct field count")
            result.append(dict(zip(SACCT_FIELDS, values)))
    return result


def collect_sacct(submissions: Path, scale_submissions: Path | None, output: Path) -> int:
    rows = load_submissions(submissions, require_all_pairs=False)
    scale_rows = read_tsv(scale_submissions) if scale_submissions and scale_submissions.exists() else []
    all_submission_rows = [*rows, *scale_rows]
    require_submission_coverage(all_submission_rows)
    job_ids = {row["job_id"] for row in all_submission_rows if row.get("job_id")}
    if not job_ids:
        raise ScaleError("zero-job completion forbidden")
    sacct_rows = query_sacct(job_ids)
    fields = ("canonical_vgp_root", "authorization_id", *SACCT_FIELDS)
    atomic_tsv(output, fields, ({
        "canonical_vgp_root": "/moosefs/erikg/vgp", "authorization_id": AUTHORIZATION_ID, **row
    } for row in sacct_rows))
    return len(sacct_rows)


def verify_pair(pair: str, vgp_root: Path, output: Path, job_id: str) -> dict[str, object]:
    if pair not in {"P04", "P07"}:
        raise ScaleError("independent scale verification is defined for completed P04/P07 only")
    # Import only here so catalog-only operations remain lightweight.  The
    # independent reviewer itself does not import the production workflow.
    from analysis import review_vgp_real_pilot as reviewer

    paths = canonical_paths(vgp_root)
    reviewer.DEFAULT_VGP_ROOT = paths["root"]
    run_root = paths["pilot_run"] / pair if pair == "P04" else paths["p07_run"]
    audited = reviewer.audit_pair(pair, run_root)
    repair = read_json(ANALYSIS / "vgp_psmc_bootstrap_repair_v1.json")
    repaired = repair["pairs"][pair]
    primary_theta = audited["psmc_recomputation"]["primary"]["theta_0_per_100bp_bin"]
    if primary_theta != repaired["primary"]["theta_0_per_100bp_bin"]:
        raise ScaleError("repaired bootstrap packet disagrees with independently parsed primary PSMC")
    if not repaired["passed"] or not repaired["masked_and_callable_sampling_population_preserved"]:
        raise ScaleError("repaired bootstrap packet did not pass")
    value = {
        "schema_version": "vgp-real-scaleout-independent-verification-v1",
        "canonical_vgp_root": str(paths["root"]), "authorization_id": AUTHORIZATION_ID,
        "selection_id": pair, "verification_job_id": job_id, "verified_at_utc": utc_now(),
        "status": "PASS", "audited_core": audited,
        "environment_capture": ({
            "path": os.environ["VGP_ENVIRONMENT_CAPTURE"],
            "sha256": sha256_file(Path(os.environ["VGP_ENVIRONMENT_CAPTURE"])),
        } if os.environ.get("VGP_ENVIRONMENT_CAPTURE") else {"status": "LOCAL_TEST_ONLY"}),
        "repaired_psmc_uncertainty": repaired,
        "same_pair_pi_psmc_non_independent": True,
        "population_inference_authorized": False,
        "scenario_uncertainty": {
            "status": "PRESERVED_GENERIC_SENSITIVITY_GRID_NOT_SPECIES_CALIBRATION",
            "scenario_count": audited["psmc_recomputation"]["scenario_count"],
            "scenario_ids": audited["psmc_recomputation"]["scenario_ids"],
            "scaling_sources": audited["psmc_recomputation"]["scaling_sources"],
        },
    }
    atomic_json(output, value)
    return value


def latest_attempts(
    submissions: Sequence[Mapping[str, str]], sacct: Sequence[Mapping[str, str]],
) -> dict[str, dict[str, object]]:
    by_job = {row["JobIDRaw"]: row for row in sacct}
    by_pair: dict[str, list[tuple[Mapping[str, str], Mapping[str, str] | None]]] = defaultdict(list)
    for row in submissions:
        by_pair[row["selection_id"]].append((row, by_job.get(row["job_id"])))
    result: dict[str, dict[str, object]] = {}
    for pair in PAIR_IDS:
        attempts = by_pair[pair]
        observed = [(submission, allocation) for submission, allocation in attempts if allocation]
        # Dependency placeholders have larger job IDs than the upstream stage
        # that failed or is still running.  They must not hide the concrete
        # execution state.  Prefer a live allocation, then a terminal error,
        # then the numerically latest observed/submitted stage.
        live = [item for item in observed if item[1]["State"].startswith("RUNNING")]
        executed = [item for item in observed if not item[1]["State"].startswith("PENDING")]
        latest_submission, latest_allocation = max(
            live or executed or observed or attempts,
            key=lambda item: int(item[0]["job_id"].split("_")[0]),
        )
        pair_allocations = [
            row for row in sacct
            if f"-{pair}-" in row.get("JobName", "")
        ]
        result[pair] = {
            "attempted_slurm_jobs": len({row[0]["job_id"] for row in attempts}),
            "completed_slurm_allocations": sum(
                row["State"].startswith("COMPLETED") for row in pair_allocations
            ),
            "latest_stage": latest_submission["stage"], "latest_job_id": latest_submission["job_id"],
            "latest_job_state": latest_allocation["State"] if latest_allocation else "SACCT_PENDING",
        }
    return result


def apply_results(
    pair_rows: list[dict[str, object]], roster: Sequence[Mapping[str, str]],
    review: Mapping[str, object], repair: Mapping[str, object], attempts: Mapping[str, Mapping[str, object]],
    verification_root: Path, scale_sacct: Sequence[Mapping[str, str]],
    fastga_contracts: Mapping[str, Mapping[str, object]],
) -> None:
    roster_by_id = {row["selection_id"]: row for row in roster}
    scale_job_by_pair = {
        row["JobName"].split("-")[2]: row for row in scale_sacct
        if row.get("JobName", "").startswith("vgp-scale-") and len(row["JobName"].split("-")) >= 4
    }
    for row in pair_rows:
        pair = str(row["selection_id"])
        if not pair:
            continue
        roster_row = roster_by_id[pair]
        row.update(attempts[pair])
        row["annotation_covariate"] = roster_row["annotation_status"]
        row["qv_covariate"] = roster_row["qv_evidence"]
        row["busco_covariate"] = roster_row["completeness_evidence"]
        row["raw_read_covariate"] = roster_row["read_hifi_provenance"]
        row["repeat_covariate"] = "UNVALIDATED_NOT_CORE_GATE"
        row["kmer_covariate"] = roster_row["duplication_collapse_evidence"]
        fastga = fastga_contracts.get(pair)
        if fastga:
            row["fastga_scratch_contract_status"] = "PASS_LIVE_PROC_NODE_LOCAL_SCRATCH"
            row["fastga_scratch_contract_job_id"] = fastga["job_id"]
            row["fastga_scratch_contract_packet"] = fastga["contract_packet"]
        else:
            row["fastga_scratch_contract_status"] = "NOT_RUN_NO_VALID_CORE"
        verification_path = verification_root / f"{pair}.json"
        scale_allocation = scale_job_by_pair.get(pair)
        if scale_allocation:
            row["verification_job_id"] = scale_allocation["JobIDRaw"]
        if verification_path.is_file():
            packet = read_json(verification_path)
            if packet.get("canonical_vgp_root") != "/moosefs/erikg/vgp" or packet.get("status") != "PASS":
                raise ScaleError(f"invalid independent verification packet: {pair}")
            core = packet["audited_core"]
            variants = core["variant_recomputation"]
            repaired = packet["repaired_psmc_uncertainty"]
            theta = repaired["bootstrap_theta_0_per_100bp_bin"]
            row.update({
                "execution_disposition": "VERIFIED_CORE_COMPLETE",
                "hard_failure_class": "", "callable_bp": variants["final_callable_bp"],
                "heterozygous_snps": variants["heterozygous_snps"], "callable_pi": variants["pi"],
                "primary_theta_0": repaired["primary"]["theta_0_per_100bp_bin"],
                "bootstrap_attempts": repaired["bootstrap_attempts"],
                "finite_bootstraps": repaired["finite_bootstraps"],
                "bootstrap_theta_q025": theta["q025"], "bootstrap_theta_median": theta["median"],
                "bootstrap_theta_q975": theta["q975"],
                "psmc_uncertainty_status": "PASS_PRIMARY_PSMCFA_BLOCK_BOOTSTRAP_CENTERED",
                "scenario_uncertainty_status": "PRESERVED_9_GENERIC_MU_X_GENERATION_SCENARIOS",
                "annotation_partition_status": (
                    "PASS_EXACT_NATIVE" if pair == "P07" else "NOT_AVAILABLE_NOT_CORE_VETO"
                ),
                "verification_status": "PASS", "result_packet": str(verification_path),
            })
            continue
        disposition = review["pair_dispositions"][pair]
        if disposition.startswith("hard_invalid"):
            row["execution_disposition"] = "HARD_INVALID_PRIMARY"
            row["hard_failure_class"] = disposition
        elif row["latest_job_state"].startswith("RUNNING"):
            row["execution_disposition"] = "RUNNING_RESUMABLE_WAVE"
            row["hard_failure_class"] = ""
        elif row["latest_job_state"].startswith(("FAILED", "CANCELLED", "TIMEOUT")):
            row["execution_disposition"] = "HARD_EXECUTION_ERROR_NO_ESTIMATE"
            row["hard_failure_class"] = "CONCRETE_EXECUTION_ERROR"
        else:
            row["execution_disposition"] = "SUBMITTED_DEPENDENCY_PENDING_OR_NEVER_SATISFIED"
            row["hard_failure_class"] = ""
        row["psmc_uncertainty_status"] = "NOT_COMPUTED_NO_VALID_CORE"
        row["scenario_uncertainty_status"] = "PRESPECIFIED_NOT_COMPUTED_NO_VALID_CORE"
        row["annotation_partition_status"] = "NOT_RUN_NO_VALID_CORE"
        row["verification_status"] = "NOT_APPLICABLE_NO_COMPLETED_CORE"
        row["result_packet"] = "NOT_APPLICABLE_NO_VALID_CORE"


def validate_fastga_contracts(vgp_root: Path, scale_root: Path) -> list[dict[str, object]]:
    """Require live /proc node-local evidence for every completed mapping."""
    results: list[dict[str, object]] = []

    p07_path = scale_root / "fastga/P07/contract.json"
    p07 = read_json(p07_path)
    if (
        p07.get("canonical_vgp_root") != str(vgp_root)
        or p07.get("selection_id") != "P07"
        or p07.get("contract_valid") is not True
        or p07.get("fastga_live_contract_valid") is not True
        or p07.get("multiplicity_contract_valid") is not True
        or not str(p07.get("private_working_directory", "")).startswith("/scratch/")
        or not all(
            str(path).startswith(str(p07.get("resolved_node_local_scratch_root", "")) + os.sep)
            for path in p07.get("observed_fastga_cwds", [])
        )
        or not all(value.get("exact_match") for value in p07["frozen_mapping_exact_reproduction"].values())
    ):
        raise ScaleError("P07 FastGA node-local scratch revalidation failed")
    results.append({
        "schema_version": "vgp-real-scaleout-fastga-contract-ledger-v1",
        "canonical_vgp_root": str(vgp_root), "selection_id": "P07",
        "job_id": p07["job_id"], "status": "PASS_LIVE_PROC_NODE_LOCAL_SCRATCH",
        "requested_scratch_root": p07["private_working_directory"],
        "resolved_node_local_scratch_root": p07["resolved_node_local_scratch_root"],
        "live_fastga_snapshots": p07["fastga_process_snapshot_count"],
        "frozen_mapping_exact_match": "true", "contract_packet": str(p07_path),
        "contract_packet_sha256": sha256_file(p07_path),
    })

    p04_root = vgp_root / "pilot/independent/P04/mapping"
    manifests = sorted(
        p04_root.glob("P04.independent.*.native.1to1.paf.manifest.json"),
        key=lambda path: int(path.name.split(".")[2]),
    )
    if not manifests:
        raise ScaleError("P04 independent FastGA node-local scratch revalidation absent")
    p04_path = manifests[-1]
    p04 = read_json(p04_path)
    manifest_suffix = ".manifest.json"
    guard_path = Path(str(p04_path)[:-len(manifest_suffix)] + ".fastga_scratch_contract.json")
    guard = read_json(guard_path)
    observed = p04["paf"]["sha256"]
    frozen_path = vgp_root / f"pilot/runs/{AUTHORIZATION_ID}-pilot-v1/P04/mapping/h2_to_h1.native.1to1.paf"
    exact = observed == sha256_file(frozen_path)
    requested_scratch = str(guard.get("canonical_requested_scratch_root", ""))
    resolved_scratch = str(guard.get("resolved_node_local_scratch_root", ""))
    if (
        p04.get("canonical_vgp_root") != str(vgp_root)
        or p04.get("selection_id") != "P04"
        or guard.get("contract_valid") is not True
        or not requested_scratch.startswith("/scratch/")
        or not all(
            str(path).startswith(resolved_scratch + os.sep)
            for path in guard.get("observed_cwds", [])
        )
        or not exact
    ):
        raise ScaleError("P04 FastGA node-local scratch revalidation failed")
    results.append({
        "schema_version": "vgp-real-scaleout-fastga-contract-ledger-v1",
        "canonical_vgp_root": str(vgp_root), "selection_id": "P04",
        "job_id": p04["job_id"], "status": "PASS_LIVE_PROC_NODE_LOCAL_SCRATCH",
        "requested_scratch_root": requested_scratch,
        "resolved_node_local_scratch_root": resolved_scratch,
        "live_fastga_snapshots": guard["fastga_snapshot_count"],
        "frozen_mapping_exact_match": "true", "contract_packet": str(p04_path),
        "contract_packet_sha256": sha256_file(p04_path),
    })
    return sorted(results, key=lambda row: row["selection_id"])


def resource_telemetry(
    sacct: Sequence[Mapping[str, str]], scale_sacct: Sequence[Mapping[str, str]], vgp_root: Path,
) -> list[dict[str, object]]:
    result: list[dict[str, object]] = []
    for scope, rows in (("pilot_and_scaleout_attempts", sacct), ("scaleout_jobs", scale_sacct)):
        states = Counter(row["State"].split()[0] for row in rows)
        result.append({
            "schema_version": "vgp-real-scaleout-telemetry-v1", "canonical_vgp_root": str(vgp_root),
            "scope": scope, "allocations": len(rows), "state_counts": json.dumps(dict(sorted(states.items()))),
            "allocated_core_hours": sum(int(row.get("CPUTimeRAW") or 0) for row in rows) / 3600,
            "elapsed_allocation_hours": sum(int(row.get("ElapsedRaw") or 0) for row in rows) / 3600,
            "maxrss_observed_allocations": sum(bool(row.get("MaxRSS")) for row in rows),
            "maxrss_not_imputed": "true", "scratch_measurement": "stage_telemetry_only_not_sacct",
        })
    return result


def write_biological_tables(output: Path, verification_root: Path, vgp_root: Path) -> dict[str, int]:
    scenario_rows: list[dict[str, object]] = []
    trajectory_rows: list[dict[str, object]] = []
    for pair in ("P04", "P07"):
        packet = read_json(verification_root / f"{pair}.json")
        run_root = Path(packet["audited_core"]["canonical_run_root"])
        if not str(run_root.resolve()).startswith(str(vgp_root) + os.sep):
            raise ScaleError("verified biological run root escaped VGP_ROOT")
        scenario_path = run_root / "psmc/finalize/scenario_scaled_trajectories.tsv"
        scenarios = read_tsv(scenario_path)
        if len(scenarios) != 576 or len({row["scenario_id"] for row in scenarios}) != 9:
            raise ScaleError(f"scenario uncertainty multiplicity drift: {pair}")
        for row in scenarios:
            scenario_rows.append({
                "canonical_vgp_root": str(vgp_root), "selection_id": pair,
                "source_path": str(scenario_path), "source_sha256": sha256_file(scenario_path), **row,
            })
        unscaled_path = run_root / "psmc/finalize/unscaled_trajectory.tsv"
        trajectories = read_tsv(unscaled_path)
        repaired = packet["repaired_psmc_uncertainty"]
        if len(trajectories) != 64:
            raise ScaleError(f"unscaled PSMC trajectory multiplicity drift: {pair}")
        for row in trajectories:
            trajectory_rows.append({
                "canonical_vgp_root": str(vgp_root), "selection_id": pair,
                "primary_theta_0_per_100bp_bin": repaired["primary"]["theta_0_per_100bp_bin"],
                "source_path": str(unscaled_path), "source_sha256": sha256_file(unscaled_path), **row,
            })
    annotation_path = vgp_root / "pilot/outputs" / AUTHORIZATION_ID / "P07/annotation/exact_partitions.json"
    annotation = read_json(annotation_path)
    if (
        annotation.get("canonical_vgp_root") != str(vgp_root)
        or annotation.get("annotation_status") != "exact_native"
        or annotation.get("sequence_dictionary_equal") is not True
    ):
        raise ScaleError("P07 exact-native annotation binding failed")
    annotation_rows = [{
        "canonical_vgp_root": str(vgp_root), "selection_id": "P07",
        "assembly_accession_version": annotation["assembly_accession_version"],
        "annotation_accession_version": annotation["annotation_accession_version"],
        "annotation_gff_sha256": annotation["annotation_gff"]["sha256"],
        "sequence_dictionary_equal": "true", "slurm_job_id": annotation["slurm_job_id"],
        **partition,
    } for partition in annotation["partitions"]]
    atomic_tsv(output / "scenario_uncertainty.tsv", tuple(scenario_rows[0]), scenario_rows)
    atomic_tsv(output / "psmc_unscaled_trajectories.tsv", tuple(trajectory_rows[0]), trajectory_rows)
    atomic_tsv(output / "exact_native_annotation_partitions.tsv", tuple(annotation_rows[0]), annotation_rows)
    return {
        "scenario_rows": len(scenario_rows), "psmc_trajectory_rows": len(trajectory_rows),
        "annotation_partition_rows": len(annotation_rows),
    }


def write_per_pair_resource_plan(output: Path, vgp_root: Path) -> dict[str, object]:
    """Bind the executed waves to the telemetry-calibrated per-pair plan."""
    source = ANALYSIS / "vgp_real_pilot_resource_plan_v1.json"
    plan = read_json(source)
    if plan.get("canonical_vgp_root") != str(vgp_root):
        raise ScaleError("resource plan escaped canonical VGP root")
    if set(plan.get("pairs", {})) != set(PAIR_IDS) - {"P07"}:
        raise ScaleError("resource plan must cover the nine non-canary audited pairs")
    if plan.get("basis", {}).get("slurm_job_id") != "1781798":
        raise ScaleError("resource plan lost its measured canary telemetry binding")
    plan["pairs"]["P07"] = {
        "measurement_sequence_bp": plan["basis"]["canary_sequence_bp"],
        "wave": 0,
        "resource_reestimate_basis": (
            "Direct canary allocation job 1781798: 32 CPUs, 128 GiB, 02:27:33 elapsed; "
            "private-scratch mapping revalidated by completed scale job 1787599."
        ),
        "stages": {"mapping": {
            "cpus_per_task": 32, "slurm_mem": "128G", "slurm_partition": "workers",
            "slurm_time": "3-00:00:00", "slurm_exclude": "octopus11",
            "scratch_bytes_high": 96636764160,
        }},
    }
    plan["scaleout_augmentation"] = {
        "selection_id": "P07", "basis": "direct measured canary telemetry",
        "canary_job_id": "1781798", "scratch_revalidation_job_id": "1787599",
    }
    target = output / "per_pair_resource_plan.json"
    atomic_json(target, plan)
    return {
        "schema_version": plan["schema_version"], "canonical_vgp_root": str(vgp_root),
        "pair_count": len(plan["pairs"]), "source_path": str(source),
        "source_sha256": sha256_file(source), "measured_canary_job_id": "1781798",
        "policy": plan["basis"]["policy"],
    }


def render_report(summary: Mapping[str, object]) -> str:
    counts = summary["counts"]
    return f"""# Real VGP core scale-out accounting

Generated: {summary['generated_at_utc']}

This is a closed-world accounting of the immutable 716-row Freeze 1 catalog:
all **569** catalog links are represented exactly once. Exactly **10** links
match the separately frozen reciprocal, same-individual, mutually comparable
roster; the other discovery links were not silently promoted to eligibility.

Real biological Slurm work was nonzero and covered all ten audited pairs.
There were {counts['submitted_biological_jobs']} distinct submitted pilot and
scale-out jobs, including {counts['scale_verification_jobs']} independent
estimate-verification jobs and {counts['fastga_scratch_revalidation_jobs']} live
FastGA `/scratch` revalidations. At this cutoff, {counts['verified_pairs']}
pairs have independently verified callable diversity and assembly-derived PSMC
packets, {counts['hard_invalid_pairs']} are hard-invalid primary executions,
{counts['execution_error_pairs']} end in concrete execution errors, and
{counts['running_pairs']} remain in a resumable wave.

Pi and PSMC are two outcomes of the same H1/H2 pair and are explicitly marked
non-independent. PSMC bootstraps sample blocks of the primary PSMCFA (including
both masked and callable N/K/T bins), and the repaired 200-replicate intervals
pass the predeclared primary-theta centering diagnostic. Absolute histories
remain a nine-scenario generic mutation-rate × generation-time sensitivity
grid; none is promoted to a species calibration. Population inference is not
authorized.

Annotation absence, external Ne, QV, BUSCO, raw reads, repeat reports, and
k-mer audits are retained as confidence covariates rather than universal core
gates. Exact-native annotation partitions are admitted only for P07; P04's
absence is not a veto. All canonical products resolve from `VGP_ROOT` =
`/moosefs/erikg/vgp`.

Both completed mappings were rerun with live `/proc` guards. Their FastGA
working directories, temporary environments, indexes, pair files, and managed
intermediates resolved beneath private node-local `/scratch`; both reruns
exactly reproduced the frozen native PAFs used by downstream pi and PSMC.
The ten-pair wave/resource plan is retained verbatim with its measured canary
job binding and per-pair memory, CPU, walltime, and scratch allocations;
concrete resource failures trigger explicit higher-allocation retries.
"""


def reconcile(args: argparse.Namespace) -> dict[str, object]:
    paths = canonical_paths(args.vgp_root)
    catalog = load_catalog(paths["catalog"])
    roster = load_roster(args.roster)
    catalog_rows, pair_rows, matched = build_closed_world(catalog, roster, paths["root"])
    pilot_submissions = load_submissions(args.pilot_submissions, require_all_pairs=False)
    scale_submissions = read_tsv(args.scale_submissions)
    submissions = [*pilot_submissions, *scale_submissions]
    require_submission_coverage(submissions)
    sacct = read_tsv(args.sacct)
    if any(row.get("canonical_vgp_root") != str(paths["root"]) for row in sacct):
        raise ScaleError("sacct snapshot escaped canonical root")
    attempts = latest_attempts(submissions, sacct)
    review = read_json(args.review)
    repair = read_json(args.repair)
    if review.get("decision") != "CONDITIONAL_GO" or "core-only scale-out" not in str(review.get("decision_scope")):
        raise ScaleError("bounded core scale-out authorization absent")
    if repair.get("passed") is not True or repair.get("canonical_vgp_root") != str(paths["root"]):
        raise ScaleError("PSMC repair authorization absent")
    scale_sacct = read_tsv(args.scale_sacct) if args.scale_sacct.is_file() else []
    if any(row.get("canonical_vgp_root") != str(paths["root"]) for row in scale_sacct):
        raise ScaleError("scale verification sacct escaped canonical root")
    scale_job_sacct = [
        row for row in scale_sacct if row.get("JobName", "").startswith("vgp-scale-")
    ]
    scale_verification_sacct = [
        row for row in scale_job_sacct if row.get("JobName", "").endswith("-verify")
    ]
    scale_fastga_sacct = [
        row for row in scale_job_sacct if row.get("JobName", "").endswith("-fastga")
    ]
    fastga_rows = validate_fastga_contracts(paths["root"], paths["scale_root"])
    fastga_by_pair = {str(row["selection_id"]): row for row in fastga_rows}
    apply_results(
        pair_rows, roster, review, repair, attempts, args.verification_root,
        scale_verification_sacct, fastga_by_pair,
    )
    eligible = [row for row in pair_rows if row["eligibility"] == "eligible_exact_audited"]
    disposition_counts = Counter(str(row["execution_disposition"]) for row in eligible)
    verified = [row for row in eligible if row["execution_disposition"] == "VERIFIED_CORE_COMPLETE"]
    scale_jobs = {row["JobIDRaw"] for row in scale_job_sacct}
    args.output.mkdir(parents=True, exist_ok=True)
    biological_table_counts = write_biological_tables(
        args.output, args.verification_root, paths["root"],
    )
    resource_plan = write_per_pair_resource_plan(args.output, paths["root"])
    summary: dict[str, object] = {
        "schema_version": "vgp-real-scaleout-summary-v1", "generated_at_utc": utc_now(),
        "canonical_vgp_root": str(paths["root"]), "authorization_id": AUTHORIZATION_ID,
        "canonical_manifest_root": str(paths["scale_root"] / "manifests"),
        "authorization": "CONDITIONAL_GO_BOUNDED_CORE_ONLY",
        "catalog": {"rows": 716, "links": 569, "commit": CATALOG_COMMIT, "sha256": CATALOG_SHA256},
        "eligible_roster": {"pairs": 10, "pair_ids": list(PAIR_IDS), "catalog_pair_mapping": matched},
        "counts": {
            "catalog_rows": 716, "catalog_links": 569, "audited_eligible_pairs": 10,
            "not_eligible_or_self_links": 559,
            "submitted_biological_jobs": len({row["job_id"] for row in submissions}),
            "scaleout_jobs": len(scale_jobs),
            "scale_verification_jobs": len({row["JobIDRaw"] for row in scale_verification_sacct}),
            "fastga_scratch_revalidation_jobs": len({row["JobIDRaw"] for row in scale_fastga_sacct}),
            "fastga_scratch_contracts_passed": len(fastga_rows), "verified_pairs": len(verified),
            "hard_invalid_pairs": disposition_counts["HARD_INVALID_PRIMARY"],
            "execution_error_pairs": disposition_counts["HARD_EXECUTION_ERROR_NO_ESTIMATE"],
            "running_pairs": disposition_counts["RUNNING_RESUMABLE_WAVE"],
            "pending_pairs": disposition_counts["SUBMITTED_DEPENDENCY_PENDING_OR_NEVER_SATISFIED"],
            **biological_table_counts,
        },
        "verified_pair_estimates": [{
            "selection_id": row["selection_id"], "callable_bp": row["callable_bp"],
            "heterozygous_snps": row["heterozygous_snps"], "callable_pi": row["callable_pi"],
            "primary_theta_0": row["primary_theta_0"], "bootstrap_attempts": row["bootstrap_attempts"],
            "finite_bootstraps": row["finite_bootstraps"],
            "bootstrap_theta_interval": [row["bootstrap_theta_q025"], row["bootstrap_theta_q975"]],
            "same_pair_pi_psmc_non_independent": True,
        } for row in verified],
        "uncertainty_contract": {
            "psmc": "primary-PSMCFA block bootstrap with masked and callable N/K/T population preserved",
            "absolute_scaling": "nine generic mutation-rate x generation-time scenarios; none preferred",
            "pair_dependence": "pi and PSMC share the same H1/H2 pair and are not independent",
            "population_inference_authorized": False,
        },
        "confidence_covariates_not_core_gates": [
            "annotation", "external_Ne", "final_QV", "BUSCO", "raw_reads",
            "repeat_reports", "kmer_audits",
        ],
        "telemetry_calibrated_per_pair_resource_plan": resource_plan,
    }
    atomic_tsv(args.output / "catalog_accounting.tsv", CATALOG_FIELDS, catalog_rows)
    atomic_tsv(args.output / "pair_accounting.tsv", PAIR_FIELDS, pair_rows)
    atomic_tsv(args.output / "fastga_scratch_contracts.tsv", tuple(fastga_rows[0]), fastga_rows)
    telemetry = resource_telemetry(sacct, scale_job_sacct, paths["root"])
    atomic_tsv(args.output / "scaleout_telemetry.tsv", tuple(telemetry[0]), telemetry)
    atomic_json(args.output / "summary.json", summary)
    atomic_text(args.output / "results.md", render_report(summary))
    return summary


def parser() -> argparse.ArgumentParser:
    value = argparse.ArgumentParser(description=__doc__)
    sub = value.add_subparsers(dest="command", required=True)
    verify = sub.add_parser("verify-pair")
    verify.add_argument("--pair", required=True, choices=("P04", "P07"))
    verify.add_argument("--vgp-root", type=Path, default=DEFAULT_VGP_ROOT)
    verify.add_argument("--output", type=Path, required=True)
    verify.add_argument("--job-id", default=os.environ.get("SLURM_JOB_ID", "local"))
    collect = sub.add_parser("collect-sacct")
    collect.add_argument("--pilot-submissions", type=Path, required=True)
    collect.add_argument("--scale-submissions", type=Path)
    collect.add_argument("--output", type=Path, required=True)
    reconcile_cmd = sub.add_parser("reconcile")
    reconcile_cmd.add_argument("--vgp-root", type=Path, default=DEFAULT_VGP_ROOT)
    reconcile_cmd.add_argument("--roster", type=Path, default=ANALYSIS / "vgp_10_pair_manifest.tsv")
    reconcile_cmd.add_argument("--review", type=Path, default=ANALYSIS / "vgp_real_pilot_independent_review_v1.json")
    reconcile_cmd.add_argument("--repair", type=Path, default=ANALYSIS / "vgp_psmc_bootstrap_repair_v1.json")
    reconcile_cmd.add_argument("--pilot-submissions", type=Path, default=DEFAULT_VGP_ROOT / "pilot/manifests/vgp10-auth-20260718-v2-pilot-v1.submissions.tsv")
    reconcile_cmd.add_argument("--scale-submissions", type=Path, default=DEFAULT_VGP_ROOT / "derived/scale-vgp-real-v1/submissions.tsv")
    reconcile_cmd.add_argument("--sacct", type=Path, required=True)
    reconcile_cmd.add_argument("--scale-sacct", type=Path, required=True)
    reconcile_cmd.add_argument("--verification-root", type=Path, default=DEFAULT_VGP_ROOT / "derived/scale-vgp-real-v1/verification")
    reconcile_cmd.add_argument("--output", type=Path, required=True)
    return value


def main(argv: Sequence[str] | None = None) -> int:
    args = parser().parse_args(argv)
    try:
        if args.command == "verify-pair":
            verify_pair(args.pair, args.vgp_root, args.output, args.job_id)
        elif args.command == "collect-sacct":
            collect_sacct(args.pilot_submissions, args.scale_submissions, args.output)
        else:
            reconcile(args)
        return 0
    except (OSError, ValueError, KeyError, ScaleError) as error:
        print(f"ERROR: {error}")
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
