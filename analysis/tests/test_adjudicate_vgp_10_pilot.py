import csv
import json
from pathlib import Path

import pytest

from analysis import adjudicate_vgp_10_pilot as run


def rows(path):
    with Path(path).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_production_closed_world_has_ten_explicit_preflight_failures_and_no_alternate():
    acquisitions = run.read_tsv(run.DEFAULT_ACQUISITION)
    inventory = run.read_tsv(run.DEFAULT_INVENTORY)
    primaries, alternates = run.validate_closed_world(acquisitions, inventory)
    assert [row["selection_id"] for row in primaries] == list(run.PRIMARY_IDS)
    assert len(alternates) == 6
    reasons = {row["selection_id"]: run.preflight_reason_codes(row) for row in primaries}
    assert all(value[0] == "MISSING_EXACT_FINAL_SEQUENCE_QV" for value in reasons.values())
    assert all("MISSING_H2_BUSCO" in value for value in reasons.values())
    assert all("MISSING_KMER_COPY_NUMBER_AUDIT" in value for value in reasons.values())
    assert all("MISSING_REPEAT_OR_LOW_COMPLEXITY_MASK" in value for value in reasons.values())
    assert {key for key, value in reasons.items() if "MISSING_H1_BUSCO" in value} == {"P07", "P08"}
    assert {key for key, value in reasons.items() if "UNRESOLVED_EXACT_READ_CHEMISTRY" in value} == {"P09"}


def test_generate_and_independently_recompute_terminal_accounting(tmp_path):
    summary = run.generate(output=tmp_path, reverify=False)
    result = run.independently_validate(tmp_path)
    assert summary["primary_slots_accounted"] == 10
    assert summary["primary_completed"] == 0
    assert summary["primary_failed_preflight"] == 10
    assert summary["alternate_activations"] == 0
    assert summary["slurm_jobs_submitted"] == 0
    assert result["primary_slots"] == result["failed_preflight"] == 10
    assert result["completed"] == result["alternate_activations"] == 0
    assert result["slurm_jobs_submitted"] == result["unresolved_identity"] == 0
    assert result["checksum_failures"] == result["bootstrap_attempts"] == 0
    assert result["artifact_packets_verified"] == 10
    assert result["validation_controls_not_run"] == 40
    assert result["reason_counts"]["MISSING_EXACT_FINAL_SEQUENCE_QV"] == 10
    manifest = rows(tmp_path / "vgp_10_pilot_result_manifest.tsv")
    qc = rows(tmp_path / "vgp_10_pilot_qc.tsv")
    telemetry = rows(tmp_path / "vgp_10_pilot_resource_telemetry.tsv")
    assert [row["selection_id"] for row in manifest] == list(run.PRIMARY_IDS)
    assert {row["terminal_state"] for row in manifest} == {"FAILED_PREFLIGHT"}
    assert {row["mapping_status"] for row in manifest} == {"not_run_preflight_failed"}
    assert {row["multiplicity_status"] for row in manifest} == {"not_measured_preflight_failed"}
    assert {row["bootstrap_attempts"] for row in manifest} == {"0"}
    assert {row["core_qc_pass"] for row in qc} == {"false"}
    assert len(telemetry) == 11
    for selection_id in run.PRIMARY_IDS:
        pair = tmp_path / "vgp_10_pilot_pair_artifacts" / selection_id
        assert {path.name for path in pair.iterdir()} == {
            "qc.json", "diversity.tsv", "psmc_trajectory.tsv", "bootstrap.tsv", "scenario.tsv",
            "validation.tsv", "annotation.json", "telemetry.json", "failure.json",
        }
        failure = json.loads((pair / "failure.json").read_text())
        assert failure["primary_retained"] is True
        assert failure["alternate_activated"] is False
        assert failure["slurm_jobs_submitted"] == 0


def test_closed_world_refuses_silent_primary_drop_or_alternate_activation():
    acquisitions = run.read_tsv(run.DEFAULT_ACQUISITION)
    inventory = run.read_tsv(run.DEFAULT_INVENTORY)
    with pytest.raises(run.AdjudicationError, match="P01..P10"):
        run.validate_closed_world(acquisitions[1:], inventory)
    changed = [dict(row) for row in acquisitions]
    alternate = next(row for row in changed if row["roster_type"] == "alternate")
    alternate["activation_status"] = "active_alternate"
    with pytest.raises(run.AdjudicationError, match="standbys"):
        run.validate_closed_world(changed, inventory)


def test_independent_validation_rejects_promoted_failure(tmp_path):
    run.generate(output=tmp_path, reverify=False)
    path = tmp_path / "vgp_10_pilot_qc.tsv"
    value = path.read_text().replace("\tfalse\t", "\ttrue\t", 1)
    path.write_text(value)
    with pytest.raises(run.AdjudicationError, match="promoted"):
        run.independently_validate(tmp_path)
