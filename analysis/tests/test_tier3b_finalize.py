import csv
import json
from pathlib import Path

import pytest

from analysis.tier3_common import Tier3ValidationError
from analysis.tier3b_finalize_run import finalize_ineligible_run


ROOT = Path(__file__).parents[2]


def _finalize(tmp_path, registry=None):
    registry_path = ROOT / "analysis/pilot_results/pilot_registry.json"
    if registry is not None:
        registry_path = tmp_path / "registry.json"
        registry_path.write_text(json.dumps(registry), encoding="utf-8")
    outputs = {
        "data_path": tmp_path / "tier3b_data.tsv",
        "failure_ledger_path": tmp_path / "tier3b_failure_ledger.tsv",
        "qc_path": tmp_path / "tier3b_qc_provenance.json",
    }
    finalize_ineligible_run(
        pilot_registry_path=registry_path,
        environment_path=ROOT / "analysis/pilot_results/guix_environment.json",
        compute_smoke_path=ROOT / "analysis/pilot_results/compute_smoke.json",
        **outputs,
    )
    return outputs


def test_no_approved_population_tuple_emits_empty_data_and_complete_failure_ledger(tmp_path):
    outputs = _finalize(tmp_path)

    data_lines = outputs["data_path"].read_text(encoding="utf-8").splitlines()
    assert len(data_lines) == 1
    assert "population_pi" in data_lines[0]
    assert "pi_S_over_pi_W" in data_lines[0]
    assert "sfs" not in data_lines[0].lower()

    with outputs["failure_ledger_path"].open(encoding="utf-8", newline="") as handle:
        failures = list(csv.DictReader(handle, dialect="excel-tab"))
    assert {row["candidate_id"] for row in failures} == {
        "dgrp_freeze2",
        "ag1000g_phase3_gambiae",
        "drosophila_simulans_170",
        "drosophila_pseudoobscura_panel",
        "aedes_aegypti_1206",
        "daphnia_pulex_panels",
    }
    assert all(row["job_status"] == "not_submitted_gate_failure" for row in failures)
    assert all(row["raw_read_calling"] == "not_launched" for row in failures)

    qc = json.loads(outputs["qc_path"].read_text(encoding="utf-8"))
    assert qc["overall_status"] == "no_eligible_population_tuples"
    assert qc["included_population_count"] == 0
    assert qc["submitted_jobs"] == []
    assert qc["pilot_gate"]["expansion_allowed"] is False
    assert qc["polarization_gate"]["sfs_b_output"] == "absent"
    assert qc["environment"]["compute_smoke"]["status"] == "passed"
    assert qc["environment"]["profile_store_path"] == qc["environment"]["compute_smoke"][
        "compute_profile_store_path"
    ]
    assert all(
        candidate["annotation"]["primary_4d_status"] == "unavailable_not_audited"
        for candidate in qc["candidates"]
    )


def test_fail_closed_finalization_is_byte_reproducible(tmp_path):
    outputs = _finalize(tmp_path)
    first = {key: path.read_bytes() for key, path in outputs.items()}
    _finalize(tmp_path)
    assert {key: path.read_bytes() for key, path in outputs.items()} == first


def test_empty_finalization_refuses_a_promoted_pilot(tmp_path):
    registry = json.loads(
        (ROOT / "analysis/pilot_results/pilot_registry.json").read_text(encoding="utf-8")
    )
    next(
        pilot for pilot in registry["pilots"] if pilot["pilot_id"] == "population_dgrp_freeze2"
    )["eligibility"] = "eligible"
    with pytest.raises(Tier3ValidationError, match="must be run, not finalized as unavailable"):
        _finalize(tmp_path, registry)
