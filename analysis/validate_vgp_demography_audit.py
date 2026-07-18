#!/usr/bin/env python3
"""Fail-closed validation for the six-species demography input audit."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from typing import Any

from analysis.build_vgp_demography_audit import (
    AUDIT, AUDIT_HEADER, CACHE_INDEX, MANIFEST, REPORT, RESOLUTION_INDEX,
    SOURCE_HEADER, SOURCES,
)
from analysis.refresh_vgp_demography_metadata import CHECKPOINT, PRIORITIZED_IDS, canonical, request_key


ROOT = Path(__file__).resolve().parents[1]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_tsv(path: Path, expected: list[str]) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != expected:
            raise RuntimeError(f"{path} header drift: {reader.fieldnames!r}")
        return list(reader)


def validate_cache() -> dict[str, Any]:
    index = json.loads(CACHE_INDEX.read_text(encoding="utf-8"))
    if index["authorization_boundary"] != {
        "biological_payloads_downloaded": 0,
        "demographic_inference_performed": False,
        "jobs_submitted": 0,
        "raw_population_bulk_download_authorized": False,
    }:
        raise RuntimeError("metadata-only authorization boundary changed")
    policy = index["lookup_policy"]
    if not policy["metadata_only"] or policy["ncbi_requests"] != 4:
        raise RuntimeError("cache is not the bounded four-request metadata audit")
    if policy["rate_limit_minimum_interval_seconds"] < 0.34:
        raise RuntimeError("NCBI rate limit interval is too short")
    if not policy["retry_after_aware"] or "batch" not in policy["batching"]:
        raise RuntimeError("cache policy lacks batching or Retry-After awareness")
    if policy["resume_checkpoint"] != "analysis/vgp_demography_cache/checkpoint.json":
        raise RuntimeError("resume checkpoint path drift")
    if len(index["responses"]) != 4:
        raise RuntimeError("expected exactly four batched NCBI responses")
    observed_keys: set[str] = set()
    db_modes: set[tuple[str, str]] = set()
    for entry in index["responses"]:
        request = entry["request"]
        key = request_key(request)
        if key != entry["request_key"] or key in observed_keys:
            raise RuntimeError("request digest mismatch or duplicate")
        observed_keys.add(key)
        response = ROOT / entry["response_path"]
        if not response.is_file() or sha256(response) != entry["response_sha256"]:
            raise RuntimeError("response path/digest mismatch: " + key)
        if response.stat().st_size != entry["response_size_bytes"]:
            raise RuntimeError("response size mismatch: " + key)
        json.loads(response.read_text(encoding="utf-8"))
        params = request["parameters"]
        db_modes.add((params["db"], "summary" if "id" in params else "search"))
        if not entry["retrieved_at_utc"] or not entry["software_environment"]["guix_channels_sha256"]:
            raise RuntimeError("response lacks retrieval/software provenance")
    if db_modes != {("taxonomy", "search"), ("taxonomy", "summary"), ("sra", "search"), ("sra", "summary")}:
        raise RuntimeError("batch endpoint coverage drift")
    checkpoint = json.loads(CHECKPOINT.read_text(encoding="utf-8"))
    if checkpoint["state"] != "complete" or checkpoint["pending_request_keys"]:
        raise RuntimeError("resume state is incomplete")
    if set(checkpoint["completed_request_keys"]) != observed_keys:
        raise RuntimeError("resume checkpoint keys do not match cache")
    if index["manifest"]["sha256"] != sha256(MANIFEST):
        raise RuntimeError("cache manifest digest mismatch")
    return index


def validate() -> dict[str, int]:
    cache = validate_cache()
    audit = read_tsv(AUDIT, AUDIT_HEADER)
    sources = read_tsv(SOURCES, SOURCE_HEADER)
    with MANIFEST.open(newline="", encoding="utf-8") as handle:
        manifest = {row["candidate_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    if set(manifest) != set(PRIORITIZED_IDS):
        raise RuntimeError("manifest denominator is not the repaired six")
    if len(audit) != 6 or {row["candidate_id"] for row in audit} != set(PRIORITIZED_IDS):
        raise RuntimeError("audit does not cover each repaired candidate exactly once")
    if len({row["candidate_id"] for row in audit}) != len(audit):
        raise RuntimeError("duplicate species audit row")

    taxonomy_entry = next(
        item for item in cache["responses"]
        if item["request"]["parameters"].get("db") == "taxonomy"
        and "id" in item["request"]["parameters"]
    )
    taxonomy_payload = json.loads((ROOT / taxonomy_entry["response_path"]).read_text(encoding="utf-8"))["result"]
    expected_taxonomy = {
        row["ncbi_taxid"]: row["scientific_name_source"] for row in manifest.values()
    }
    observed_taxonomy = {
        uid: taxonomy_payload[uid]["scientificname"] for uid in taxonomy_payload["uids"]
    }
    if observed_taxonomy != expected_taxonomy:
        raise RuntimeError("batched NCBI taxonomy does not exactly match manifest names/TaxIds")

    required = set(AUDIT_HEADER)
    for row in audit:
        missing = sorted(key for key in required if row[key] == "")
        if missing:
            raise RuntimeError(f"implicit missingness for {row['candidate_id']}: {missing}")
        parent = manifest[row["candidate_id"]]
        if row["scientific_name"] != parent["scientific_name_source"] or row["ncbi_taxid"] != parent["ncbi_taxid"]:
            raise RuntimeError("taxonomy mismatch: " + row["candidate_id"])
        if row["exact_reference_accession_version"] != parent["h1_accession_version"]:
            raise RuntimeError("exact-reference mismatch: " + row["candidate_id"])
        if parent["h1_exact_version_status"] != "official_current_exact_version":
            raise RuntimeError("manifest exact-reference status is not resolved")
        if row["biosample_accession"] != parent["biosample_accession"]:
            raise RuntimeError("BioSample mismatch")
        if row["psmc_eligible"] not in {"yes", "no"} or row["msmc2_eligible"] not in {"yes", "no"} or row["smcpp_eligible"] not in {"yes", "no"}:
            raise RuntimeError("invalid method-specific eligibility token")
        if row["psmc_eligible"] == "yes" and any("missing" in row[key].lower() for key in ["callable_diploid_genome_status", "compatible_mask_status", "coverage_provenance_status"]):
            raise RuntimeError("PSMC promoted despite missing prerequisite")
        if row["msmc2_eligible"] == "yes" and any("missing" in row[key].lower() for key in ["phasing_accuracy_status", "compatible_mask_status", "individual_population_relationship_status"]):
            raise RuntimeError("MSMC2 promoted despite missing prerequisite")
        if row["smcpp_eligible"] == "yes" and any("missing" in row[key].lower() for key in ["population_genotype_status", "population_definition_status", "population_reference_status", "population_mask_status", "population_qc_status"]):
            raise RuntimeError("SMC++ promoted despite missing prerequisite")
        if not row["psmc_blockers"] or not row["msmc2_blockers"] or not row["smcpp_blockers"]:
            raise RuntimeError("method-specific blockers must be distinct and explicit")
        if "not generated" not in row["coalescent_scaled_output_status"] and row["scientific_name"] != "Heterodontus francisci":
            raise RuntimeError("unexpected coalescent output")
        if "requires explicit mutation-rate and generation-time scenario" not in row["absolute_ne_time_status"] and row["scientific_name"] not in {"Camelus dromedarius", "Heterodontus francisci"}:
            raise RuntimeError("absolute scaling separation missing")
        keys = row["metadata_cache_request_keys"].split(";")
        if not set(item["request_key"] for item in cache["responses"]).issubset(keys):
            raise RuntimeError("audit row omits metadata-cache provenance")

    record_ids = [row["record_id"] for row in sources]
    if len(record_ids) != len(set(record_ids)) or any(not value for value in record_ids):
        raise RuntimeError("duplicate or empty source record_id")
    for row in sources:
        blank = sorted(key for key in SOURCE_HEADER if row[key] == "")
        if blank:
            raise RuntimeError(f"implicit source missingness for {row['record_id']}: {blank}")
        if row["candidate_id"] not in manifest:
            raise RuntimeError("source references non-denominator species")
        parent = manifest[row["candidate_id"]]
        if row["scientific_name"] != parent["scientific_name_source"] or row["ncbi_taxid"] != parent["ncbi_taxid"]:
            raise RuntimeError("source taxonomy mismatch")
        if row["classification"] in {"independent_literature_ne", "independent_literature_nb"}:
            for key in ["estimand", "estimand_definition", "method", "population", "geography", "measurement_time", "uncertainty_status", "source_kind", "source_title", "source_year", "source_locator"]:
                if not row[key]:
                    raise RuntimeError(f"independent estimate lacks {key}: {row['record_id']}")
            if row["response_dataset_overlap"].startswith("derived") or row["circularity_status"] != "non_circular":
                raise RuntimeError("circular estimate promoted as independent")
        if row["classification"] == "coalescent_scaled_not_absolute_ne":
            if row["estimand"] != "theta=4Ne_mu" or "not" not in row["mutation_rate"]:
                raise RuntimeError("coalescent scale was silently made absolute")
        if row["classification"] == "excluded_circular":
            if row["record_status"] != "excluded_if_proposed" or row["circularity_status"] != "circular_excluded":
                raise RuntimeError("circular row not excluded")
            if row["value"] != "not calculated" or row["response_dataset_overlap"] != "derived_from_response":
                raise RuntimeError("circularity guard calculated/promoted a value")
    circular_ids = {row["candidate_id"] for row in sources if row["classification"] == "excluded_circular"}
    if circular_ids != set(PRIORITIZED_IDS):
        raise RuntimeError("each species needs an explicit circularity exclusion")
    represented = {row["candidate_id"] for row in sources}
    if represented != set(PRIORITIZED_IDS):
        raise RuntimeError("source audit lacks an explicit record for a species")

    resolution = json.loads(RESOLUTION_INDEX.read_text(encoding="utf-8"))
    for entry in resolution["responses"]:
        response = ROOT / entry["response_path"]
        if sha256(response) != entry["response_sha256"]:
            raise RuntimeError("inherited exact-reference cache digest failure")

    text = REPORT.read_text(encoding="utf-8")
    for phrase in [
        "None of the six", "Coalescent scale versus absolute scale", "Independent Ne, census, and circularity",
        "zero biological payloads", "159.95 MB", "PSMC", "MSMC2", "SMC++", "QUOTA_UNAVAILABLE",
    ]:
        if phrase not in text:
            raise RuntimeError("report missing required statement: " + phrase)
    return {"species": len(audit), "sources": len(sources), "cache_responses": len(cache["responses"])}


if __name__ == "__main__":
    result = validate()
    print("VGP_DEMOGRAPHY_AUDIT_OK " + " ".join(f"{key}={value}" for key, value in result.items()))
