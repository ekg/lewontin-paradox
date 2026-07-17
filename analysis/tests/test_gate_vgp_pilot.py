import copy
import csv
import json
from pathlib import Path

import pytest

from analysis import gate_vgp_pilot as gate
from analysis.tier3_common import Tier3ValidationError


ROOT = Path(__file__).parents[2]


def _build(tmp_path: Path, **overrides):
    gate_out = tmp_path / "vgp_pilot_gate.json"
    review_out = tmp_path / "vgp_pilot_gate_review.md"
    kwargs = {
        "manifest_path": ROOT / "analysis/vgp_pilot_manifest.tsv",
        "size_budget_path": ROOT / "analysis/vgp_pilot_size_budget.tsv",
        "freeze_provenance_path": ROOT / "analysis/vgp_phase1_freeze_provenance.json",
        "root_config_path": ROOT / "analysis/vgp_data_root_config.json",
        "root_validation_path": ROOT / "analysis/vgp_data_root_validation.json",
        "decisions_path": ROOT / "analysis/vertebrate_scaleout_decisions.tsv",
        "execution_plan_path": ROOT / "analysis/vertebrate_scaleout_execution_plan.md",
        "resource_budget_path": ROOT / "analysis/vertebrate_scaleout_resource_budget.tsv",
        "gate_out": gate_out,
        "review_out": review_out,
    }
    kwargs.update(overrides)
    built = gate.build_gate(**kwargs)
    return built, gate_out, review_out


def _write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def test_build_gate_recomputes_current_artifacts_and_emits_nogo(tmp_path):
    built, gate_out, review_out = _build(tmp_path)
    assert gate_out.is_file()
    assert review_out.is_file()
    assert built["decision"]["status"] == "NO_GO"
    assert built["reproduction"]["manifest_candidate_count"] == 74
    assert built["reproduction"]["selected_row_count"] == 0
    codes = {item["code"] for item in built["blockers"]}
    assert "NO_SELECTED_ROWS" in codes
    assert "ZERO_COMPOSITION_ELIGIBLE_ROWS" in codes
    assert "ZERO_DIVERSITY_ELIGIBLE_ROWS" in codes
    assert "QUOTA_UNAVAILABLE" in codes
    assert built["cap_vector"]["dimensions"]["persistent_input_gb"]["value"] == 0.0
    assert built["cap_vector"]["dimensions"]["aggregate_core_hours"]["value"] == 0.0
    assert built["authorization_boundary"]["cap_vector_digest"] == built["cap_vector"]["sha256"]


def test_authorize_rejects_nogo_for_acquire_and_compute(tmp_path):
    built, gate_out, _review_out = _build(tmp_path)
    assert built["decision"]["status"] == "NO_GO"
    with pytest.raises(Tier3ValidationError, match="gate decision is NO_GO"):
        gate.authorize_gate_action(
            gate_out,
            ROOT / "analysis/vgp_pilot_manifest.tsv",
            ROOT / "analysis/vgp_data_root_config.json",
            "acquire",
        )
    with pytest.raises(Tier3ValidationError, match="gate decision is NO_GO"):
        gate.authorize_gate_action(
            gate_out,
            ROOT / "analysis/vgp_pilot_manifest.tsv",
            ROOT / "analysis/vgp_data_root_config.json",
            "compute",
        )


def test_authorize_rejects_mutated_manifest_and_changed_root(tmp_path):
    _built, gate_out, _review_out = _build(tmp_path)

    mutated_manifest = tmp_path / "mutated_manifest.tsv"
    manifest_text = (ROOT / "analysis/vgp_pilot_manifest.tsv").read_text(encoding="utf-8")
    mutated_manifest.write_text(manifest_text.replace("rejected_unresolved", "rejected_CHANGED", 1), encoding="utf-8")
    with pytest.raises(Tier3ValidationError, match="manifest digest mismatch"):
        gate.authorize_gate_action(
            gate_out,
            mutated_manifest,
            ROOT / "analysis/vgp_data_root_config.json",
            "acquire",
        )

    changed_root = tmp_path / "changed_root.json"
    root_config = json.loads((ROOT / "analysis/vgp_data_root_config.json").read_text(encoding="utf-8"))
    root_config["root"] = root_config["root"] + "-changed"
    changed_root.write_text(json.dumps(root_config, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    with pytest.raises(Tier3ValidationError, match="root contract digest mismatch"):
        gate.authorize_gate_action(
            gate_out,
            ROOT / "analysis/vgp_pilot_manifest.tsv",
            changed_root,
            "compute",
        )


def test_authorize_rejects_tampered_gate_payload(tmp_path):
    _built, gate_out, _review_out = _build(tmp_path)
    tampered = json.loads(gate_out.read_text(encoding="utf-8"))
    tampered["decision"]["status"] = "GO"
    gate_out.write_text(json.dumps(tampered, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    with pytest.raises(Tier3ValidationError, match="gate decision hash does not match"):
        gate.authorize_gate_action(
            gate_out,
            ROOT / "analysis/vgp_pilot_manifest.tsv",
            ROOT / "analysis/vgp_data_root_config.json",
            "acquire",
        )


def test_authorize_rejects_changed_live_cap_vector(tmp_path):
    _built, gate_out, _review_out = _build(tmp_path)
    rows = gate.load_tsv(ROOT / "analysis/vertebrate_scaleout_resource_budget.tsv")
    mutated = copy.deepcopy(rows)
    for row in mutated:
        if row["stage_or_dataset"] == "stratified_pilot" and row["scenario"] == "high":
            row["peak_aggregate_bandwidth_mib_s"] = "119.0"
            break
    mutated_path = tmp_path / "resource_budget_mutated.tsv"
    _write_tsv(mutated_path, mutated)
    with pytest.raises(Tier3ValidationError, match="cap vector digest mismatch"):
        gate.authorize_gate_action(
            gate_out,
            ROOT / "analysis/vgp_pilot_manifest.tsv",
            ROOT / "analysis/vgp_data_root_config.json",
            "acquire",
            resource_budget_path=mutated_path,
        )


def test_build_gate_fails_closed_when_exact_size_is_missing(tmp_path):
    rows = gate.load_tsv(ROOT / "analysis/vgp_pilot_size_budget.tsv")
    mutated = copy.deepcopy(rows)
    mutated[0]["download_bytes_exact"] = ""
    mutated_path = tmp_path / "size_budget_missing.tsv"
    _write_tsv(mutated_path, mutated)

    built, _gate_out, _review_out = _build(tmp_path, size_budget_path=mutated_path)
    codes = {item["code"] for item in built["blockers"]}
    assert "SIZE_BUDGET_MISSING_DOWNLOAD_BYTES" in codes
    assert built["decision"]["status"] == "NO_GO"


def test_build_gate_fails_closed_when_quota_metadata_is_missing(tmp_path):
    validation = json.loads((ROOT / "analysis/vgp_data_root_validation.json").read_text(encoding="utf-8"))
    del validation["system_evidence"]["quota_state"]
    validation_path = tmp_path / "root_validation_missing_quota.json"
    validation_path.write_text(json.dumps(validation, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    built, _gate_out, _review_out = _build(tmp_path, root_validation_path=validation_path)
    codes = {item["code"] for item in built["blockers"]}
    assert "QUOTA_METADATA_MISSING" in codes
    assert "QUOTA_UNAVAILABLE" in codes
    assert built["decision"]["status"] == "NO_GO"
