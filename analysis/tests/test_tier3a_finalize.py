import csv
import copy
import json
from pathlib import Path

import pytest

from analysis.tier3_common import Tier3ValidationError
from analysis.tier3a_finalize_run import DATA_COLUMNS, finalize_ineligible_run


ROOT = Path(__file__).resolve().parents[2]
PILOTS = ROOT / "analysis/pilot_results/pilot_registry.json"
ENVIRONMENT = ROOT / "analysis/pilot_results/guix_environment.json"
SMOKE = ROOT / "analysis/pilot_results/compute_smoke.json"


def _run(tmp_path, *, pilots=PILOTS, environment=ENVIRONMENT, smoke=SMOKE):
    data = tmp_path / "tier3a_data.tsv"
    failures = tmp_path / "tier3a_failure_ledger.tsv"
    qc = tmp_path / "tier3a_qc_provenance.json"
    result = finalize_ineligible_run(
        pilot_registry_path=pilots,
        environment_path=environment,
        compute_smoke_path=smoke,
        data_path=data,
        failure_ledger_path=failures,
        qc_path=qc,
    )
    return result, data, failures, qc


def test_finalize_vgp_run_emits_structured_missingness_not_zero(tmp_path):
    result, data, failures, qc_path = _run(tmp_path)
    with data.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows == []
    assert data.read_text(encoding="utf-8").rstrip("\n").split("\t") == list(DATA_COLUMNS)

    with failures.open(newline="", encoding="utf-8") as handle:
        failure_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(failure_rows) == 2
    assert {row["job_status"] for row in failure_rows} == {"not_submitted_gate_failure"}
    assert {row["direct_wfmash_status"] for row in failure_rows} == {"not_run_preflight_gate_failure"}
    assert all(row["impg_status"] == "not_run_execution_not_approved" for row in failure_rows)

    committed = json.loads(qc_path.read_text(encoding="utf-8"))
    assert committed == result
    assert result["overall_status"] == "no_eligible_vgp_individual_tuples"
    assert result["included_individual_count"] == 0
    assert result["submitted_jobs"] == []
    assert result["result_table"]["rows"] == 0
    assert "zero heterozygosity" in result["result_table"]["interpretation"]
    assert result["input_audit"]["vgp_phase1_freeze_inventory_sha256"] == (
        "9c58420484a8b76a2d6175b7c26bf709e68bdc726a67fc7541b8c2b5a2fc13a4"
    )
    assert result["input_audit"]["buffalo_core_sha256"] == (
        "df559451dad94b53ba8675e09811708107a57eeb6ffe8f72b944bcbbf3a1f2eb"
    )


def test_candidate_qc_preserves_native_annotation_and_denominator_gates(tmp_path):
    result, *_ = _run(tmp_path)
    by_id = {candidate["candidate_id"]: candidate for candidate in result["candidates"]}
    pilot = by_id["vgp_dual_modality_individual"]
    assert pilot["dual_modality_concordance"]["promotion_passed"] is False
    assert pilot["deposited_modality"]["unique_callable_h1_bases"] is None
    assert pilot["statistics"]["total_denominator"] is None
    assert pilot["statistics"]["fourfold_W_denominator"] is None
    assert pilot["statistics"]["fourfold_S_denominator"] is None
    assert pilot["annotation"]["projected_annotation_used_for_primary"] is False
    assert pilot["annotation"]["sampled_cds_reconstruction"] == "not_run"

    local = by_id["mmyodau2.1_local_inventory"]
    assert local["reference_tuple"]["h1_accession_version"] == "GCF_963259705.1"
    assert local["reference_tuple"]["h2_accession_version"] == "GCA_963242275.1"
    assert local["reference_tuple"]["exact_reference_tuple_frozen"] is False
    assert local["phasing_and_collapse"]["h1_h2_phase_identity_audit_passed"] is False
    assert local["phasing_and_collapse"]["collapse_duplication_qc_passed"] is False
    assert "missing_exact_buffalo_species_covariate" in local["blocking_codes"]


def test_environment_audit_pins_wfmash_and_keeps_impg_disabled(tmp_path):
    result, *_ = _run(tmp_path)
    environment = result["environment"]
    assert environment["wfmash_commit"] == "e040aa10e87cab44ed5a4db005e784be62b0bd21"
    assert environment["impg_commit"] == "101df81eb28a809c8fac97d297acd9fcfbbfa048"
    assert environment["impg_executable_present"] is False
    assert environment["forbidden_tools_used"] == []
    assert environment["compute_smoke"]["wfmash_extended_cigar"]["extended_cigar_passed"] is True
    assert environment["compute_smoke"]["bcf_csi_passed"] is True


@pytest.mark.parametrize("mutation,match", [
    (lambda value: value["pilots"][-1].update(eligibility="eligible"), "must be run"),
    (lambda value: value["pilots"][-1]["optional_impg"].update(execution_approved=True), "IMPG execution"),
])
def test_finalizer_refuses_changed_frozen_pilot(tmp_path, mutation, match):
    value = json.loads(PILOTS.read_text(encoding="utf-8"))
    mutation(value)
    changed = tmp_path / "pilots.json"
    changed.write_text(json.dumps(value), encoding="utf-8")
    with pytest.raises(Tier3ValidationError, match=match):
        _run(tmp_path, pilots=changed)


def test_finalizer_rejects_unpinned_wfmash_environment(tmp_path):
    value = copy.deepcopy(json.loads(ENVIRONMENT.read_text(encoding="utf-8")))
    value["tool_versions"]["wfmash"] = "wfmash 0.12.5"
    changed = tmp_path / "environment.json"
    changed.write_text(json.dumps(value), encoding="utf-8")
    with pytest.raises(Tier3ValidationError, match="pinned WFMASH"):
        _run(tmp_path, environment=changed)


def test_finalization_is_byte_deterministic(tmp_path):
    first = tmp_path / "first"
    second = tmp_path / "second"
    _run(first)
    _run(second)
    for name in ("tier3a_data.tsv", "tier3a_failure_ledger.tsv", "tier3a_qc_provenance.json"):
        assert (first / name).read_bytes() == (second / name).read_bytes()
