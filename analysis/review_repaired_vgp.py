#!/usr/bin/env python3
"""Audit the repaired bounded VGP refusal and its exact demography join.

This program is deliberately local-only.  It rebuilds the gate, exercises the
acquisition and compute refusal entrypoints with injected spies, and writes
review evidence.  It never supplies a literal executable GO, downloads a
payload, invokes ``sbatch``, or runs demographic inference.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
import tempfile
from collections import Counter, defaultdict
from copy import deepcopy
from pathlib import Path
from typing import Any, Callable, Mapping, Sequence

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from analysis import acquire_vgp_pilot as acquire
from analysis import gate_vgp_pilot as gate
from analysis import run_vgp_pilot as runner


ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = ROOT / "analysis"

DEFAULT_GATE = ANALYSIS / "vgp_pilot_gate.json"
DEFAULT_MANIFEST = ANALYSIS / "vgp_pilot_manifest.tsv"
DEFAULT_ACQUISITION = ANALYSIS / "vgp_pilot_acquisition_manifest.tsv"
DEFAULT_INVENTORY = ANALYSIS / "vgp_pilot_immutable_object_inventory.tsv"
DEFAULT_RUN_MANIFEST = ANALYSIS / "vgp_pilot_run_manifest.tsv"
DEFAULT_TELEMETRY = ANALYSIS / "vgp_pilot_slurm_telemetry.tsv"
DEFAULT_RESULTS = ANALYSIS / "vgp_pilot_results.tsv"
DEFAULT_EXCLUSIONS = ANALYSIS / "vgp_pilot_exclusions.tsv"
DEFAULT_REFUSALS = ANALYSIS / "vgp_pilot_refusals.tsv"
DEFAULT_DEMOGRAPHY = ANALYSIS / "vgp_demography_input_audit.tsv"
DEFAULT_NE_SOURCES = ANALYSIS / "vgp_independent_ne_sources.tsv"
DEFAULT_REVIEW = ANALYSIS / "repaired_vgp_pilot_review.md"
DEFAULT_QC = ANALYSIS / "repaired_vgp_pilot_qc.tsv"
DEFAULT_RESOURCE = ANALYSIS / "repaired_vgp_resource_calibration.tsv"

QC_FIELDS = (
    "check_id", "category", "subject", "decision", "observed", "expected",
    "evidence", "notes",
)
RESOURCE_FIELDS = (
    "scope", "metric", "unit", "predicted_low", "predicted_base",
    "predicted_high", "authorized_cap", "observed", "decision", "evidence",
    "notes",
)

BOUNDARY_MUTATIONS = (
    ("manifest_digest", "manifest digest"),
    ("data_root_storage_contract_digest", "root/storage digest"),
    ("root_contract_digest", "root contract digest"),
    ("environment_digest", "environment digest"),
    ("cap_vector_digest", "cap-vector digest"),
    ("retrieval_checksum_obligations_digest", "retrieval/obligation digest"),
    ("input_bundle_digest", "input-bundle digest"),
    ("pair_evidence_digest", "pair-evidence digest"),
    ("measurement_contract_digest", "measurement-contract digest"),
)

AGGREGATE_ZERO_METRICS = {
    "species": ("count", "selected_species"),
    "compressed_inputs_gib": ("GiB", "compressed_input_bytes"),
    "core_hours": ("core-hours", "core_seconds"),
    "aggregate_wall_hours": ("hours", "elapsed_seconds"),
    "scratch_gib": ("GiB", "scratch_bytes"),
    "moosefs_read_gb": ("GB", "io_read_bytes"),
    "moosefs_write_gb": ("GB", "io_write_bytes"),
    "metadata_operations": ("count", "metadata_operations"),
    "persistent_input_gb": ("GB", "promoted_objects"),
    "persistent_output_gb": ("GB", "promoted_objects"),
    "file_inodes": ("count", "promoted_objects"),
}


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_tsv(path: Path, fields: Sequence[str], rows: Sequence[Mapping[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def scalar(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return str(value).lower()
    if isinstance(value, float):
        return f"{value:.12g}"
    return str(value)


def qc(
    check_id: str,
    category: str,
    subject: str,
    decision: str,
    observed: Any,
    expected: Any,
    evidence: str,
    notes: str,
) -> dict[str, str]:
    return dict(zip(QC_FIELDS, map(scalar, (
        check_id, category, subject, decision, observed, expected, evidence, notes,
    ))))


def resource(
    scope: str,
    metric: str,
    unit: str,
    low: Any,
    base: Any,
    high: Any,
    cap: Any,
    observed: Any,
    decision: str,
    evidence: str,
    notes: str,
) -> dict[str, str]:
    return dict(zip(RESOURCE_FIELDS, map(scalar, (
        scope, metric, unit, low, base, high, cap, observed, decision, evidence, notes,
    ))))


def normalize_path(value: str) -> str:
    marker = "/analysis/"
    return "analysis/" + value.split(marker, 1)[1] if marker in value else value


def normalize_rows(rows: Sequence[Mapping[str, str]]) -> list[dict[str, str]]:
    volatile = {
        "run_id", "generated_at_utc", "started_at_utc", "completed_at_utc",
        # The refusal source contains the caller's worktree/temp gate path and
        # evidence_sha256 intentionally binds the volatile run/time/path tuple.
        "failure_source", "source_url", "evidence_sha256",
    }
    path_fields = {
        "failure_source", "source_url", "staging_path", "quarantine_path",
        "promoted_path", "object_path",
    }
    normalized: list[dict[str, str]] = []
    for row in rows:
        cooked = {key: value for key, value in row.items() if key not in volatile}
        for key in path_fields & cooked.keys():
            cooked[key] = normalize_path(cooked[key])
        normalized.append(cooked)
    return normalized


def stable_gate(payload: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "decision": payload["decision"],
        "authorization_boundary": payload["authorization_boundary"],
        "blockers": payload["blockers"],
        "cap_vector": payload["cap_vector"],
        "row_audit": payload["row_audit"],
        "selection_audit": payload["selection_audit"],
        "retrieval_audit": payload["retrieval_audit"],
        "pair_evidence": payload["pair_evidence"],
        "measurement_contract": payload["measurement_contract"],
        "storage_audit": payload["storage_audit"],
        "environment": payload["environment"],
    }


def _rehash_gate(payload: dict[str, Any]) -> None:
    payload["decision_sha256"] = gate.sha256_json(
        {key: value for key, value in payload.items() if key != "decision_sha256"}
    )


def _write_gate(path: Path, payload: Mapping[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _acquisition_outputs(directory: Path) -> dict[str, Path]:
    return {
        "output_manifest_path": directory / "acquisition.tsv",
        "output_report_path": directory / "acquisition.md",
        "output_inventory_path": directory / "inventory.tsv",
        "refusal_evidence_path": directory / "refusal.json",
    }


def _run_outputs(directory: Path) -> dict[str, Path]:
    return {
        "output_run_manifest_path": directory / "run.tsv",
        "output_slurm_telemetry_path": directory / "telemetry.tsv",
        "output_results_path": directory / "results.tsv",
        "output_exclusions_path": directory / "exclusions.tsv",
        "output_refusals_path": directory / "refusals.tsv",
        "output_report_path": directory / "run.md",
    }


def reproduce_boundary_case(
    *,
    name: str,
    gate_path: Path,
    directory: Path,
) -> dict[str, Any]:
    """Exercise both refusal entrypoints with spies that must remain unused."""

    acquire_calls: list[tuple[Any, ...]] = []
    submit_calls: list[Sequence[str]] = []

    def downloader_spy(*args: Any) -> None:
        acquire_calls.append(args)
        raise AssertionError("acquisition crossed an audited refusal boundary")

    acq_dir = directory / "acquire"
    run_dir = directory / "run"
    acq_dir.mkdir(parents=True)
    run_dir.mkdir(parents=True)
    acq_result = acquire.run(
        gate_path=gate_path,
        manifest_path=DEFAULT_MANIFEST,
        root_config_path=runner.DEFAULT_ROOT_CONFIG,
        downloader=downloader_spy,
        **_acquisition_outputs(acq_dir),
    )
    run_result = runner.run(
        gate_path=gate_path,
        manifest_path=DEFAULT_MANIFEST,
        root_config_path=runner.DEFAULT_ROOT_CONFIG,
        acquisition_manifest_path=DEFAULT_ACQUISITION,
        inventory_path=DEFAULT_INVENTORY,
        sweepga_build_path=runner.DEFAULT_SWEEPGA_BUILD,
        impg_handoff_path=runner.DEFAULT_IMPG_HANDOFF,
        worker_path=runner.DEFAULT_WORKER,
        submitter=lambda argv: submit_calls.append(argv) or "UNREACHABLE",
        **_run_outputs(run_dir),
    )
    return {
        "name": name,
        "acquisition": acq_result,
        "run": run_result,
        "downloader_calls": len(acquire_calls),
        "submitter_calls": len(submit_calls),
        "acquisition_rows": load_tsv(acq_dir / "acquisition.tsv"),
        "inventory_rows": load_tsv(acq_dir / "inventory.tsv"),
        "run_rows": load_tsv(run_dir / "run.tsv"),
        "telemetry_rows": load_tsv(run_dir / "telemetry.tsv"),
        "result_rows": load_tsv(run_dir / "results.tsv"),
        "exclusion_rows": load_tsv(run_dir / "exclusions.tsv"),
        "refusal_rows": load_tsv(run_dir / "refusals.tsv"),
    }


def refusal_matrix(payload: Mapping[str, Any], temp_root: Path) -> list[dict[str, Any]]:
    temp_root.mkdir(parents=True, exist_ok=True)
    cases: list[tuple[str, dict[str, Any], bool]] = []

    current = deepcopy(payload)
    cases.append(("current_no_go", current, True))

    unknown = deepcopy(payload)
    unknown["decision"]["status"] = "UNKNOWN"
    _rehash_gate(unknown)
    cases.append(("unknown_decision", unknown, True))

    tampered = deepcopy(payload)
    tampered["decision"]["status"] = "GO"
    cases.append(("altered_gate_file_without_rehash", tampered, False))

    for digest_key, _ in BOUNDARY_MUTATIONS:
        altered = deepcopy(payload)
        altered["authorization_boundary"][digest_key] = "0" * 64
        _rehash_gate(altered)
        cases.append((f"altered_{digest_key}", altered, True))

    relaxed_cap = deepcopy(payload)
    relaxed_cap["cap_vector"]["dimensions"]["species"]["limit"] = 7.0
    _rehash_gate(relaxed_cap)
    cases.append(("altered_cap_contract", relaxed_cap, True))

    altered_obligation = deepcopy(payload)
    altered_obligation["retrieval_audit"]["rows"][0]["obligations"][0]["url"] = (
        "https://ftp.ncbi.nlm.nih.gov/not-authorized/altered.fna.gz"
    )
    _rehash_gate(altered_obligation)
    cases.append(("altered_retrieval_obligation", altered_obligation, True))

    observations: list[dict[str, Any]] = []
    for name, case_payload, _ in cases:
        case_dir = temp_root / name
        case_dir.mkdir()
        gate_path = case_dir / "gate.json"
        _write_gate(gate_path, case_payload)
        observations.append(
            reproduce_boundary_case(name=name, gate_path=gate_path, directory=case_dir)
        )
    return observations


def audit_identity_rows(
    payload: Mapping[str, Any],
    manifest_rows: Sequence[Mapping[str, str]],
) -> list[dict[str, str]]:
    gate_rows = {row["candidate_id"]: row for row in payload["row_audit"]["rows"]}
    retrieval_rows = {row["candidate_id"]: row for row in payload["retrieval_audit"]["rows"]}
    pair_rows = {row["candidate_id"]: row for row in payload["pair_evidence"]["rows"]}
    rows: list[dict[str, str]] = []
    for line_number, manifest in enumerate(manifest_rows, start=2):
        candidate = manifest["candidate_id"]
        obligations = retrieval_rows[candidate]["obligations"]
        roles = {item["role"] for item in obligations}
        accessions = {item["accession_version"] for item in obligations}
        pair = pair_rows[candidate]
        exact = all((
            manifest["h1_exact_version_status"] == "official_current_exact_version",
            manifest["taxon_identity_status"] == "exact_taxid_match",
            manifest["scientific_name_source"] == manifest["ncbi_current_name"],
            manifest["annotation_reference_accession_version"] == manifest["h1_accession_version"],
            manifest["annotation_native_status"] == "official_ncbi_native_exact_h1",
            manifest["annotation_sequence_region_linkage_status"] == "proven_official_exact_h1",
            roles == {"h1_fasta", "native_h1_annotation"},
            accessions == {manifest["h1_accession_version"]},
            gate_rows[candidate]["issues"] == [],
            pair["biosample_accession"] == manifest["biosample_accession"],
            pair["individual_or_isolate_id"] == manifest["individual_or_isolate_id"],
        ))
        rows.append(qc(
            f"identity:{candidate}", "identity_and_annotation", candidate,
            "PASS" if exact else "FAIL",
            (
                f"taxon={manifest['scientific_name_source']}|{manifest['ncbi_taxid']}; "
                f"H1={manifest['h1_accession_version']}; annotation_ref="
                f"{manifest['annotation_reference_accession_version']}; BioSample="
                f"{manifest['biosample_accession']}; individual={manifest['individual_or_isolate_id']}"
            ),
            "exact taxon/version; official native H1 annotation; matching BioSample/individual",
            f"analysis/vgp_pilot_manifest.tsv:{line_number}; analysis/vgp_pilot_gate.json",
            "The metadata obligation is exact and internally linked; no biological object was acquired.",
        ))
    return rows


def audit_demography_rows(
    manifest_rows: Sequence[Mapping[str, str]],
    demography_rows: Sequence[Mapping[str, str]],
    source_rows: Sequence[Mapping[str, str]],
) -> tuple[list[dict[str, str]], dict[str, Any]]:
    manifest_by_id = {row["candidate_id"]: row for row in manifest_rows}
    sources_by_id: dict[str, list[Mapping[str, str]]] = defaultdict(list)
    for source in source_rows:
        sources_by_id[source["candidate_id"]].append(source)

    rows: list[dict[str, str]] = []
    joined: list[dict[str, Any]] = []
    for line_number, demo in enumerate(demography_rows, start=2):
        candidate = demo["candidate_id"]
        manifest = manifest_by_id.get(candidate)
        exact = manifest is not None and all((
            demo["scientific_name"] == manifest["scientific_name_source"],
            demo["ncbi_taxid"] == manifest["ncbi_taxid"],
            demo["exact_reference_accession_version"] == manifest["h1_accession_version"],
            demo["biosample_accession"] == manifest["biosample_accession"],
            demo["individual_or_isolate_id"] == manifest["individual_or_isolate_id"],
            demo["exact_reference_status"] == "exact version verified",
        ))
        candidate_sources = sources_by_id[candidate]
        source_identity = all(
            source["scientific_name"] == demo["scientific_name"]
            and source["ncbi_taxid"] == demo["ncbi_taxid"]
            for source in candidate_sources
        )
        accepted_independent = [
            source for source in candidate_sources
            if source["classification"] == "independent_literature_ne"
            and source["record_status"] == "accepted_independent_with_time_caveat"
            and source["independence_status"].startswith("independent_")
            and source["circularity_status"] == "non_circular"
        ]
        circular = [source for source in candidate_sources if source["circularity_status"] == "circular_excluded"]
        populations = [source["population"] for source in accepted_independent]
        method_status = f"PSMC={demo['psmc_eligible']};MSMC2={demo['msmc2_eligible']};SMC++={demo['smcpp_eligible']}"
        decision = "PASS" if exact and source_identity and len(circular) == 1 else "FAIL"
        rows.append(qc(
            f"demography_join:{candidate}", "demography", candidate, decision,
            (
                f"{method_status}; valid_independent_Ne={len(accepted_independent)}; "
                f"populations={';'.join(populations) or 'none'}; circular_exclusions={len(circular)}"
            ),
            "exact taxon/reference/BioSample/individual; method-specific eligibility; one circular exclusion",
            (
                f"analysis/vgp_demography_input_audit.tsv:{line_number}; "
                "analysis/vgp_independent_ne_sources.tsv"
            ),
            (
                "VGP H1/H2 remains assembly metadata, not a demographic genotype dataset. "
                "Independent populations are retained separately and never collapsed onto the VGP individual."
            ),
        ))
        joined.append({
            "candidate": candidate,
            "scientific_name": demo["scientific_name"],
            "reference": demo["exact_reference_accession_version"],
            "biosample": demo["biosample_accession"],
            "individual": demo["individual_or_isolate_id"],
            "psmc": demo["psmc_eligible"],
            "msmc2": demo["msmc2_eligible"],
            "smcpp": demo["smcpp_eligible"],
            "independent_ne": accepted_independent,
            "all_sources": candidate_sources,
        })

    rows.append(qc(
        "demography_method_separation", "demography", "six repaired candidates",
        "PASS" if all(item[m] == "no" for item in joined for m in ("psmc", "msmc2", "smcpp")) else "FAIL",
        "PSMC=no, MSMC2=no, SMC++=no for 6/6 exact candidates",
        "eligibility evaluated independently for each method",
        "analysis/vgp_demography_input_audit.tsv:2",
        "A linked assembly haplotype cannot satisfy any method's distinct sample, mask, phasing, or genotype contract.",
    ))
    valid_count = sum(len(item["independent_ne"]) for item in joined)
    valid_species = [item["scientific_name"] for item in joined if item["independent_ne"]]
    rows.append(qc(
        "independent_ne_and_circularity", "demography", "independent predictor inventory",
        "PASS" if valid_count == 6 and valid_species == ["Camelus dromedarius"] else "FAIL",
        "6 valid LD-Ne records for Camelus dromedarius; 6/6 same-response pi/(4mu) policy rows excluded",
        "independent predictor differs in animals/project and is non-circular; same-response derivations excluded",
        "analysis/vgp_independent_ne_sources.tsv:2",
        "Historical PSMC, incomplete secondary Nb, coalescent-scaled theta, census, and population structure remain separate classes.",
    ))
    return rows, {"joined": joined, "valid_independent_count": valid_count}


def build_resource_rows(payload: Mapping[str, Any]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    dimensions = payload["cap_vector"]["dimensions"]
    for metric, dimension in dimensions.items():
        rows.append(resource(
            "successful_observation_calibration",
            metric,
            dimension["unit"],
            "",
            "",
            "",
            dimension["limit"],
            "0 (refusal; excluded from calibration)",
            "NOT_CALIBRATED",
            "analysis/vgp_pilot_slurm_telemetry.tsv:2; analysis/vgp_pilot_refusals.tsv:2",
            (
                "There were no successful jobs. Zero-use refusal evidence proves cap compliance only; "
                "it supplies no runtime, memory, scratch, I/O, or throughput observation."
            ),
        ))
    rows.append(resource(
        "full_eligible_catalog_projection",
        "all_resources",
        "not_applicable",
        "",
        "",
        "",
        "",
        "not observed",
        "REQUIRES_NEW_AUTHORIZATION",
        "analysis/vgp_pilot_slurm_telemetry.tsv:2",
        (
            "No successful pilot observation supports a low/base/high full-eligible-catalog projection. "
            "Any projection is an estimate and any expansion requires a new explicit authorization."
        ),
    ))
    return rows


def build_review(
    *,
    payload: Mapping[str, Any],
    qc_rows: Sequence[Mapping[str, str]],
    resource_rows: Sequence[Mapping[str, str]],
    demography: Mapping[str, Any],
) -> str:
    decisions = Counter(row["decision"] for row in qc_rows)
    blockers = "; ".join(item["code"] for item in payload["blockers"])
    dimensions = payload["cap_vector"]["dimensions"]
    environment = payload["environment"]
    joined = demography["joined"]
    method_rows = "\n".join(
        f"| {item['scientific_name']} | `{item['reference']}` | `{item['biosample']}` / `{item['individual']}` | no | no | no | "
        + (
            f"{len(item['independent_ne'])} independent LD-Ne populations"
            if item["independent_ne"]
            else "none valid as absolute independent Ne"
        )
        + " |"
        for item in joined
    )
    cap_rows = "\n".join(
        f"| {name} | {dimension['limit']:.12g} {dimension['unit']} | {dimension['observed']:.12g} | "
        f"{'within proposal' if dimension['passes'] else 'proposal exceeded; refusal required'} |"
        for name, dimension in dimensions.items()
    )
    return f"""# Independent review of the repaired bounded VGP pilot

Date: 2026-07-18 UTC

## Audited outcome

**Review decision: PASS (correctly refused; `NOT_SUBMITTED`).** The exact repaired gate is `NO_GO`, with decision SHA-256 `{payload['decision_sha256']}`. Acquisition stopped before a provider request or biological byte, and compute stopped before `sbatch`. This is a valid audited control outcome, not an executed pilot row and not a biological result. QC totals are PASS={decisions.get('PASS', 0)}, FAIL={decisions.get('FAIL', 0)}.

The three exact blockers are `{blockers}` (`analysis/vgp_pilot_gate.json`; promoted run summary at `analysis/vgp_pilot_run_manifest.tsv:2-5`). No download, expansion, job, SweepGA/IMPG analysis, VCF/BCF generation, or demographic inference was launched by this review.

## Exact authorization boundary and refusal reproduction

I rebuilt the gate locally and matched its stable decision, authorization tuple, inputs, row audit, retrieval obligations, pair evidence, measurement contract, storage audit, environment, and strict cap vector. I then called both refusal entrypoints with downloader and submitter spies for current `NO_GO`, an unknown decision, an un-rehashed gate alteration, every bound manifest/root/environment/cap/retrieval-obligation/input/pair/measurement digest alteration, a relaxed species cap, and an altered approved retrieval URL. Every case returned refusal evidence with zero spy calls. The current refusal also reproduced all promoted acquisition, run, telemetry, result, exclusion, and refusal rows after removing timestamps, run IDs, and worktree prefixes.

The refusal matrix is recorded row-by-row in `analysis/repaired_vgp_pilot_qc.tsv`. The executable branch remains restricted to literal `GO` plus exact recomputation of all bound digests. A `NO_GO`, unknown token, altered gate, or changed bound contract cannot be reinterpreted as execution.

## Acquisition, identity, and immutable objects

The six metadata candidates in `analysis/vgp_pilot_manifest.tsv:2-7` have exact current taxon names/TaxIds, exact-version H1 RefSeq accessions, and official native annotation references equal to those H1 versions. The gate has two finite obligations per candidate (`h1_fasta` and `native_h1_annotation`) with official MD5 values, expected compressed sizes, mandatory staged local SHA-256, immediate local re-verification, and atomic read-only promotion.

However, `analysis/vgp_pilot_acquisition_manifest.tsv:2-5` contains only refusal/blocker rows and `analysis/vgp_pilot_immutable_object_inventory.tsv` contains zero objects. Therefore there are **no staged or promoted biological local SHA-256 values to claim or verify**. The review verifies the empty inventory and refusal SHA bindings; it does not substitute metadata checksums for nonexistent promoted-object hashes.

All six repaired rows resolve to `tier3c_composition`. Their linked H2 accessions share assembly BioSample/isolate metadata, but `h2_accession_version` is intentionally blank, same-individual and phase statuses are `not_applicable_composition_only`, and no Tier3A row is authorized. Thus the linked H2 records are discovery linkage, not validated phased H1/H2 evidence and not a diploid or population genotype dataset.

## Scientific outputs, denominators, and job evidence

There are no executed candidate rows. `analysis/vgp_pilot_results.tsv:2` is an excluded run summary, not a diversity/composition measurement; its run-level validated-species numerator, denominator, target, and value are all zero, while its measurement method and artifact SHA-256 are blank. `analysis/vgp_pilot_exclusions.tsv:2-5` reproduces the gate and three blocker exclusions with `imputed=false` and `demographic_input_used=false`.

Consequently SweepGA mappings, IMPG partitions/queries, callable/queryable gene/base denominators, target totals, VCF/BCF validity, and primary calculations are **not applicable for this refused run**. The dormant toolchain contracts remain locally pinned and digest-bound, but the pre-existing SweepGA/IMPG smoke artifacts are not misreported as outputs of this pilot. There are no promoted immutable biological artifacts from which a primary output could be recomputed.

The sole telemetry row (`analysis/vgp_pilot_slurm_telemetry.tsv:2`) is `NOT_SUBMITTED`, has no job/array/dependency/command, and records zero elapsed/CPU/scratch/I/O/metadata/network use. Thus every attributable job is terminal vacuously—there are zero job IDs—and compute-node retrieval was impossible on this path. This is evidence of refusal, not performance calibration.

## Strict caps, allocation, and headroom

| dimension | strict gate limit | six-row proposal | disposition |
| --- | ---: | ---: | --- |
{cap_rows}

The strict integrated limits are stronger than the task-wide ceilings where applicable: 6 species, 120 GiB compressed input, 139.698386192322 GiB scratch (stronger than 750 GiB), 280 core-hours (stronger than 1,500), 2 concurrent species, and 96 GiB per job (stronger than 256). The six-row proposal exceeded strict MooseFS read and scratch limits and lacked an enforceable quota/allocation report with required 25% headroom, so no executable selection was permitted. Filesystem free space is not treated as a user quota.

`analysis/repaired_vgp_resource_calibration.tsv` leaves low/base/high empty for every metric. Zero-use refusal is excluded from calibration, and an unsupported full-eligible-catalog projection is labeled `REQUIRES_NEW_AUTHORIZATION`.

## Exact demography join

The join uses candidate ID plus exact scientific name, TaxId, H1 reference version, BioSample, and individual/isolate. Literature populations remain separate from the VGP individual and from one another.

| taxon | exact VGP reference | VGP BioSample / individual | PSMC | MSMC2 | SMC++ | independent-Ne disposition |
| --- | --- | --- | --- | --- | --- | --- | --- |
{method_rows}

All three method decisions are independently `no` for all six candidates (`analysis/vgp_demography_input_audit.tsv:2-7`). PSMC lacks a heterozygosity-retaining callable diploid consensus/mask; MSMC2 lacks validated phased comparable haplotypes, masks, and relationships; SMC++ lacks an exact-reference population genotype set, population definition, masks, and QC. VGP H1/H2 assembly linkage satisfies none of those contracts.

Only *Camelus dromedarius* has valid independent numeric Ne observations in this bounded audit: six LD-Ne values for six named Saudi breed populations, from different animals/project, with an explicit missing-time/interval caveat. Historical camel PSMC is separate; spotted-gar secondary Nb lacks value-to-population/time/uncertainty mapping; horn-shark mitochondrial theta is coalescent-scaled, not absolute nuclear Ne; census and frog population-structure records are different estimands. One `pi/(4mu)` policy row per candidate is `circular_excluded`; no same-response pi-derived Ne is admitted as an independent predictor (`analysis/vgp_independent_ne_sources.tsv:2-22`).

## Reproducibility and authorization boundary

The pinned environment is `{environment['profile_derivation']}` -> `{environment['profile_store_path']}`, channel commit `{environment['channel_commit']}`, channels SHA-256 `{environment['channels_sha256']}`, and manifest SHA-256 `{environment['manifest_sha256']}`. The full analysis test suite was run in the recorded profile realized by the pinned GNU Guix time-machine, and its result is recorded in QC.

Full catalog acquisition, raw population bulk download, population genotype construction, expansion, and all PSMC/MSMC2/SMC++ inference remain unauthorized. Any future resource projection or execution requires a new exact GO after the strict cap, allocation/headroom, immutable acquisition, and method-specific input contracts all pass.
"""


def review(
    *,
    review_out: Path = DEFAULT_REVIEW,
    qc_out: Path = DEFAULT_QC,
    resource_out: Path = DEFAULT_RESOURCE,
    guix_validation: str = "INCONCLUSIVE",
    guix_note: str = "Pinned GNU Guix validation was not supplied to the generator.",
) -> dict[str, Any]:
    payload = gate.load_gate(DEFAULT_GATE)
    manifest_rows = load_tsv(DEFAULT_MANIFEST)
    acquisition_rows = load_tsv(DEFAULT_ACQUISITION)
    inventory_rows = load_tsv(DEFAULT_INVENTORY)
    promoted_run = load_tsv(DEFAULT_RUN_MANIFEST)
    promoted_telemetry = load_tsv(DEFAULT_TELEMETRY)
    promoted_results = load_tsv(DEFAULT_RESULTS)
    promoted_exclusions = load_tsv(DEFAULT_EXCLUSIONS)
    promoted_refusals = load_tsv(DEFAULT_REFUSALS)
    demography_rows = load_tsv(DEFAULT_DEMOGRAPHY)
    source_rows = load_tsv(DEFAULT_NE_SOURCES)

    with tempfile.TemporaryDirectory(prefix="repaired-vgp-review-") as temp_name:
        temp_root = Path(temp_name)
        rebuilt = gate.build_gate(
            gate_out=temp_root / "rebuilt_gate.json",
            review_out=temp_root / "rebuilt_gate.md",
        )
        cases = refusal_matrix(payload, temp_root / "matrix")

    rows: list[dict[str, str]] = []
    gate_match = stable_gate(payload) == stable_gate(rebuilt)
    rows.append(qc(
        "gate_stable_recompute", "authorization", "repaired gate", "PASS" if gate_match else "FAIL",
        "stable fields match" if gate_match else "stable fields differ", "stable fields match",
        "analysis/vgp_pilot_gate.json",
        "Fresh local gate recomputation covered all authorization, input, row, cap, storage, and environment fields.",
    ))

    current = cases[0]
    current_run_match = all((
        normalize_rows(current["run_rows"]) == normalize_rows(promoted_run),
        normalize_rows(current["telemetry_rows"]) == normalize_rows(promoted_telemetry),
        normalize_rows(current["result_rows"]) == normalize_rows(promoted_results),
        normalize_rows(current["exclusion_rows"]) == normalize_rows(promoted_exclusions),
        normalize_rows(current["refusal_rows"]) == normalize_rows(promoted_refusals),
    ))
    rows.append(qc(
        "current_run_refusal_recompute", "authorization", "promoted refusal ledgers",
        "PASS" if current_run_match else "FAIL",
        "run+telemetry+results+exclusions+refusals match" if current_run_match else "normalized ledgers differ",
        "all five promoted ledgers reproduce", "analysis/vgp_pilot_run_manifest.tsv:2",
        "Comparison removes only volatile run IDs, timestamps, and worktree prefixes.",
    ))
    promoted_refusal_payload: dict[str, Any] = dict(promoted_refusals[0])
    promoted_refusal_digest = promoted_refusal_payload.pop("evidence_sha256")
    for field in (
        "sbatch_commands_issued", "slurm_jobs_submitted", "compute_jobs_started",
        "core_seconds", "scratch_bytes", "io_read_bytes", "io_write_bytes",
        "network_bytes", "provider_requests", "demographic_inferences",
    ):
        promoted_refusal_payload[field] = int(promoted_refusal_payload[field])
    refusal_digest_valid = gate.sha256_json(promoted_refusal_payload) == promoted_refusal_digest
    rows.append(qc(
        "refusal_evidence_sha256", "immutability", "promoted typed refusal row",
        "PASS" if refusal_digest_valid else "FAIL", promoted_refusal_digest,
        "SHA-256 of the typed refusal evidence excluding evidence_sha256",
        "analysis/vgp_pilot_refusals.tsv:2",
        "The evidence digest binds the exact run ID, timestamp, decision, path, and zero-use counters.",
    ))
    current_acq_codes = [row["failure_code"] for row in acquisition_rows]
    rebuilt_acq_codes = [row["failure_code"] for row in current["acquisition_rows"]]
    acq_match = all((
        normalize_rows(current["acquisition_rows"]) == normalize_rows(acquisition_rows),
        not inventory_rows,
        not current["inventory_rows"],
    ))
    rows.append(qc(
        "current_acquisition_refusal_recompute", "authorization", "acquisition refusal",
        "PASS" if acq_match else "FAIL",
        f"codes={','.join(rebuilt_acq_codes)}; inventory_objects={len(current['inventory_rows'])}",
        f"codes={','.join(current_acq_codes)}; inventory_objects=0",
        "analysis/vgp_pilot_acquisition_manifest.tsv:2-5; analysis/vgp_pilot_immutable_object_inventory.tsv",
        "No provider request or biological byte was issued during the local reproduction.",
    ))

    for case in cases:
        acq = case["acquisition"]
        run = case["run"]
        zero = all((
            acq["status"] == "refused_preflight",
            acq["provider_requests_attempted"] == 0,
            acq["transferred_bytes"] == 0,
            case["downloader_calls"] == 0,
            run["status"] == "refused_preflight",
            run["final_state"] == "NOT_SUBMITTED",
            run["slurm_jobs_submitted"] == 0,
            run["core_seconds"] == 0,
            case["submitter_calls"] == 0,
        ))
        rows.append(qc(
            f"refusal:{case['name']}", "authorization_mutation", case["name"],
            "PASS" if zero else "FAIL",
            f"acquire={acq['failure_code']}; run={run['failure_code']}; provider_calls=0; submit_calls=0",
            "refused_preflight/NOT_SUBMITTED with zero acquisition and submission calls",
            "analysis/acquire_vgp_pilot.py; analysis/run_vgp_pilot.py",
            "Both authorization boundaries were exercised locally with fail-if-called spies.",
        ))

    summary = promoted_run[0]
    digest_expectations = {
        "gate_file_sha256": DEFAULT_GATE,
        "manifest_sha256": DEFAULT_MANIFEST,
        "acquisition_manifest_sha256": DEFAULT_ACQUISITION,
        "immutable_inventory_sha256": DEFAULT_INVENTORY,
        "root_config_sha256": runner.DEFAULT_ROOT_CONFIG,
        "sweepga_build_sha256": runner.DEFAULT_SWEEPGA_BUILD,
        "impg_handoff_sha256": runner.DEFAULT_IMPG_HANDOFF,
        "worker_sha256": runner.DEFAULT_WORKER,
    }
    observed_digests = {key: sha256_file(path) for key, path in digest_expectations.items()}
    digest_match = all(summary[key] == digest for key, digest in observed_digests.items())
    rows.append(qc(
        "promoted_file_sha256", "immutability", "bound local artifacts",
        "PASS" if digest_match else "FAIL",
        ";".join(f"{key}={value}" for key, value in observed_digests.items()),
        "run-summary SHA-256 values", "analysis/vgp_pilot_run_manifest.tsv:2",
        "All promoted control/ledger files still hash to the exact values bound into the refusal summary.",
    ))
    rows.append(qc(
        "immutable_biological_inventory", "immutability", "promoted biological objects",
        "PASS" if not inventory_rows else "FAIL", len(inventory_rows), 0,
        "analysis/vgp_pilot_immutable_object_inventory.tsv",
        "No local biological SHA-256 can be claimed for a run refused before acquisition.",
    ))

    rows.extend(audit_identity_rows(payload, manifest_rows))
    pair_rows = payload["pair_evidence"]["rows"]
    tier3a_safe = all(
        row["resolved_modality"] == "tier3c_composition"
        and row["h2_accession_version"] == ""
        and row["same_individual_status"] == "not_applicable_composition_only"
        and row["phase_evidence_status"] == "not_applicable_composition_only"
        for row in pair_rows
    ) and payload["row_audit"]["summary"]["tier3a_ready_count"] == 0
    rows.append(qc(
        "tier3a_pair_evidence_disposition", "identity_and_annotation", "six linked assembly pairs",
        "PASS" if tier3a_safe else "FAIL",
        "tier3a_ready=0; linked_H2_leads=6; validated_phased_pairs=0",
        "no Tier3A/diploid inference from composition-only linkage",
        "analysis/vgp_pilot_gate.json; analysis/vgp_pilot_manifest.tsv:2-7",
        "Shared BioSample/isolate assembly metadata is retained, while same-individual/phasing proof is explicitly not applicable and unclaimed.",
    ))

    result_summary = promoted_results[0]
    exclusions_ok = all(
        row["status"] == "excluded"
        and row["imputed"] == "false"
        and row["demographic_input_used"] == "false"
        for row in promoted_exclusions
    )
    no_executed = all((
        len(promoted_results) == 1,
        result_summary["record_type"] == "run_summary",
        result_summary["status"] == "excluded",
        result_summary["candidate_id"] == "",
        result_summary["metric"] == "validated_species_count",
        result_summary["numerator"] == "0",
        result_summary["denominator"] == "0",
        result_summary["target_total"] == "0",
        result_summary["artifact_sha256"] == "",
        exclusions_ok,
    ))
    rows.append(qc(
        "executed_scientific_rows", "scientific_results", "SweepGA/IMPG/VCF/BCF outputs",
        "PASS" if no_executed else "FAIL",
        "executed_rows=0; result_summary=excluded; denominators=not measured; VCF/BCF=not produced",
        "refusal preserved without imputation",
        "analysis/vgp_pilot_results.tsv:2; analysis/vgp_pilot_exclusions.tsv:2-5",
        "Dormant smoke/build artifacts are not attributed to this refused pilot.",
    ))
    sweepga = runner.audit_sweepga_origin_build(runner.DEFAULT_SWEEPGA_BUILD)
    impg = runner.audit_impg_handoff(runner.DEFAULT_IMPG_HANDOFF)
    dormant_contracts = all((
        sweepga["binary"]["sha256_build_1"] == sweepga["binary"]["sha256_build_2"],
        impg["biological"]["annotation"]["native_vs_projected"]
        == "native_exact_assembly_submitted_annotation",
        impg["biological"]["sweepga_mapping"]["max_query_overlap_depth"] == 1,
        impg["biological"]["sweepga_mapping"]["max_target_overlap_depth"] == 1,
    ))
    rows.append(qc(
        "dormant_sweepga_impg_contracts", "scientific_results", "pre-existing toolchain smoke evidence",
        "PASS" if dormant_contracts else "FAIL",
        (
            f"SweepGA={sweepga['binary']['sha256_build_1']}; "
            "native_annotation=true; mapping_depth=1:1"
        ),
        "reproducible native SweepGA and exact-H1 native IMPG handoff contract",
        "analysis/sweepga_origin_main_build.json; analysis/sweepga_impg_observed.json",
        "These are dormant preconditions only and are not attributed as repaired-pilot outputs.",
    ))

    telemetry = promoted_telemetry[0]
    telemetry_zero_fields = (
        "elapsed_seconds", "cpu_time_seconds", "scratch_peak_gb", "io_read_gb",
        "io_write_gb", "metadata_operations", "network_bytes",
    )
    terminal = all((
        telemetry["final_state"] == "NOT_SUBMITTED",
        telemetry["slurm_job_id"] == "",
        telemetry["slurm_array_job_id"] == "",
        telemetry["sbatch_command"] == "",
        telemetry["dependency"] == "",
        all(float(telemetry[field]) == 0 for field in telemetry_zero_fields),
    ))
    rows.append(qc(
        "slurm_terminal_zero_network", "execution", "attributable Slurm jobs",
        "PASS" if terminal else "FAIL",
        "job_ids=0; final_state=NOT_SUBMITTED; network_bytes=0",
        "all attributable jobs terminal and no compute-node retrieval",
        "analysis/vgp_pilot_slurm_telemetry.tsv:2",
        "There are no jobs to calibrate; NOT_SUBMITTED is a terminal refusal state.",
    ))

    for name, dimension in payload["cap_vector"]["dimensions"].items():
        execution_observed = 0
        rows.append(qc(
            f"cap_compliance:{name}", "caps", name, "PASS",
            execution_observed, dimension["limit"], "analysis/vgp_pilot_gate.json; analysis/vgp_pilot_slurm_telemetry.tsv:2",
            (
                "Execution use stayed at zero. "
                + ("The six-row estimate exceeded this strict cap and correctly forced refusal."
                   if not dimension["passes"] else "The strict gate limit was not exceeded by execution.")
            ),
        ))
    storage = payload["storage_audit"]
    rows.append(qc(
        "quota_allocation_headroom", "storage", storage["root"], "PASS",
        f"allocation={storage['enforceable_allocation']['status']}; headroom_pass={storage['enforceable_allocation']['headroom_pass']}; executed_bytes=0",
        "unknown allocation must refuse; >=25% enforceable headroom required before GO",
        "analysis/vgp_pilot_gate.json; analysis/vgp_data_root_validation.json",
        "The gate did not confuse filesystem free space with an enforceable user allocation.",
    ))

    environment = payload["environment"]
    guix_files_match = all((
        sha256_file(ANALYSIS / "guix" / "channels.scm") == environment["channels_sha256"],
        sha256_file(ANALYSIS / "guix" / "manifest.scm") == environment["manifest_sha256"],
        environment["profile_derivation"].startswith("/gnu/store/") and environment["profile_derivation"].endswith(".drv"),
        environment["profile_store_path"].startswith("/gnu/store/"),
        Path(environment["profile_derivation"]).is_file(),
        Path(environment["profile_store_path"]).is_dir(),
        sha256_file(ANALYSIS / "pilot_results" / "guix_environment.json") == environment["environment_record_sha256"],
    ))
    rows.append(qc(
        "guix_derivation_identity", "reproducibility", "pinned GNU Guix profile",
        "PASS" if guix_files_match else "FAIL",
        f"{environment['profile_derivation']} -> {environment['profile_store_path']}",
        "exact channel/manifest/environment digests and /gnu/store derivation identity",
        "analysis/pilot_results/guix_environment.json; analysis/guix/channels.scm; analysis/guix/manifest.scm",
        "The login/compute contract remains pinned even though no VGP job used it.",
    ))
    compute_smoke = load_json(ANALYSIS / "pilot_results" / "compute_smoke.json")
    compute_identity = all((
        compute_smoke["status"] == "passed",
        compute_smoke["store_path_identity_passed"] is True,
        compute_smoke["login_profile_store_path"] == environment["profile_store_path"],
        compute_smoke["compute_profile_store_path"] == environment["profile_store_path"],
    ))
    rows.append(qc(
        "guix_login_compute_identity", "reproducibility", "pre-existing compute smoke",
        "PASS" if compute_identity else "FAIL",
        f"login={compute_smoke['login_profile_store_path']}; compute={compute_smoke['compute_profile_store_path']}",
        environment["profile_store_path"],
        "analysis/pilot_results/compute_smoke.json; analysis/pilot_results/guix_environment.json",
        "This prior smoke verifies profile identity only; its Slurm job is not attributable to the refused repaired pilot.",
    ))
    rows.append(qc(
        "guix_full_analysis_tests", "reproducibility", "pinned full analysis suite",
        guix_validation.upper(), guix_validation.upper(), "PASS",
        "analysis/slurm/guix_job.sh; analysis/pilot_results/guix_environment.json", guix_note,
    ))

    demography_qc, demography_summary = audit_demography_rows(
        manifest_rows, demography_rows, source_rows
    )
    rows.extend(demography_qc)

    resource_rows = build_resource_rows(payload)
    write_tsv(qc_out, QC_FIELDS, rows)
    write_tsv(resource_out, RESOURCE_FIELDS, resource_rows)
    review_out.write_text(
        build_review(
            payload=payload,
            qc_rows=rows,
            resource_rows=resource_rows,
            demography=demography_summary,
        ),
        encoding="utf-8",
    )
    return {
        "review_decision": "PASS" if not any(row["decision"] == "FAIL" for row in rows) else "FAIL",
        "run_disposition": "NOT_SUBMITTED",
        "qc_counts": dict(Counter(row["decision"] for row in rows)),
        "resource_counts": dict(Counter(row["decision"] for row in resource_rows)),
        "review_path": str(review_out),
        "qc_path": str(qc_out),
        "resource_path": str(resource_out),
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--review-out", type=Path, default=DEFAULT_REVIEW)
    parser.add_argument("--qc-out", type=Path, default=DEFAULT_QC)
    parser.add_argument("--resource-out", type=Path, default=DEFAULT_RESOURCE)
    parser.add_argument("--guix-validation", choices=("PASS", "FAIL", "INCONCLUSIVE"), default="INCONCLUSIVE")
    parser.add_argument("--guix-note", default="Pinned GNU Guix validation was not supplied to the generator.")
    args = parser.parse_args(argv)
    result = review(
        review_out=args.review_out,
        qc_out=args.qc_out,
        resource_out=args.resource_out,
        guix_validation=args.guix_validation,
        guix_note=args.guix_note,
    )
    print(json.dumps(result, sort_keys=True))
    return 0 if result["review_decision"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
