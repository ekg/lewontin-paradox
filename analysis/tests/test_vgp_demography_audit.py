import csv
import json

import pytest

from analysis import build_vgp_demography_audit as builder
from analysis import validate_vgp_demography_audit as validator
from analysis.refresh_vgp_demography_metadata import PRIORITIZED_IDS, request_key


def _rows(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _write(path, header, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def test_built_audit_covers_repaired_six_and_validates():
    builder.build()
    result = validator.validate()
    assert result == {"species": 6, "sources": 21, "cache_responses": 4}
    rows = _rows(builder.AUDIT)
    assert {row["candidate_id"] for row in rows} == set(PRIORITIZED_IDS)
    assert all(row["psmc_eligible"] == "no" for row in rows)
    assert all(row["msmc2_eligible"] == "no" for row in rows)
    assert all(row["smcpp_eligible"] == "no" for row in rows)


def test_cache_requests_are_unique_batched_digested_and_resumable():
    index = json.loads(builder.CACHE_INDEX.read_text(encoding="utf-8"))
    keys = [item["request_key"] for item in index["responses"]]
    assert len(keys) == len(set(keys)) == 4
    assert all(item["request_key"] == request_key(item["request"]) for item in index["responses"])
    assert index["lookup_policy"]["rate_limit_minimum_interval_seconds"] >= 0.34
    assert index["lookup_policy"]["retry_after_aware"] is True
    checkpoint = json.loads((builder.ROOT / index["lookup_policy"]["resume_checkpoint"]).read_text(encoding="utf-8"))
    assert checkpoint["state"] == "complete"
    assert checkpoint["pending_request_keys"] == []
    assert set(checkpoint["completed_request_keys"]) == set(keys)


def test_validator_rejects_psmc_promotion_from_haploid_h1(tmp_path, monkeypatch):
    rows = _rows(builder.AUDIT)
    rows[0]["psmc_eligible"] = "yes"
    path = tmp_path / "audit.tsv"
    _write(path, builder.AUDIT_HEADER, rows)
    monkeypatch.setattr(validator, "AUDIT", path)
    with pytest.raises(RuntimeError, match="PSMC promoted despite missing prerequisite"):
        validator.validate()


def test_validator_rejects_circular_estimate_promotion(tmp_path, monkeypatch):
    rows = _rows(builder.SOURCES)
    circular = next(row for row in rows if row["classification"] == "excluded_circular")
    circular["record_status"] = "accepted_independent"
    path = tmp_path / "sources.tsv"
    _write(path, builder.SOURCE_HEADER, rows)
    monkeypatch.setattr(validator, "SOURCES", path)
    with pytest.raises(RuntimeError, match="circular row not excluded"):
        validator.validate()


def test_output_separates_coalescent_and_absolute_scaling():
    rows = {row["scientific_name"]: row for row in _rows(builder.AUDIT)}
    horn = rows["Heterodontus francisci"]
    assert "theta=4Ne-mu" in horn["coalescent_scaled_output_status"]
    assert "not converted" in horn["absolute_ne_time_status"]
    camel = rows["Camelus dromedarius"]
    assert camel["mutation_rate_scenario"].startswith("1.1e-8")
    assert camel["generation_time_scenario"].startswith("5 years")
    assert camel["smcpp_eligible"] == "no"


def test_each_species_has_explicit_circularity_guard_and_source_disposition():
    rows = _rows(builder.SOURCES)
    circular = {row["candidate_id"] for row in rows if row["classification"] == "excluded_circular"}
    represented = {row["candidate_id"] for row in rows}
    assert circular == represented == set(PRIORITIZED_IDS)
    assert all(row["value"] == "not calculated" for row in rows if row["classification"] == "excluded_circular")
