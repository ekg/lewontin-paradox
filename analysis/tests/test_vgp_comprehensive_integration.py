import json
from pathlib import Path

import pytest

from analysis import validate_vgp_comprehensive_integration as integration


ROOT = Path(__file__).parents[2]


def test_integrated_chain_rosters_tooling_and_history_validate():
    assert integration.validate() == {
        "alternate_pairs": 6,
        "historical_artifacts": 11,
        "historical_runs": 5,
        "primary_pairs": 10,
        "source_commits": 17,
        "tasks": 13,
        "tool_checks": 12,
    }


def test_historical_decision_cannot_be_promoted_as_current_authorization(tmp_path):
    registry = json.loads(integration.HISTORY.read_text(encoding="utf-8"))
    registry["authorization_policy"]["current_execution_authorization"] = (
        "analysis/vgp_10_pilot_review_decision.json"
    )
    changed = tmp_path / "history.json"
    changed.write_text(json.dumps(registry), encoding="utf-8")
    with pytest.raises(integration.IntegrationError, match="must not be current"):
        integration.verify_historical_registry(changed, ROOT)


def test_historical_artifact_digest_drift_fails_closed(tmp_path):
    registry = json.loads(integration.HISTORY.read_text(encoding="utf-8"))
    registry["runs"][0]["artifacts"][0]["sha256"] = "0" * 64
    changed = tmp_path / "history.json"
    changed.write_text(json.dumps(registry), encoding="utf-8")
    with pytest.raises(integration.IntegrationError, match="artifact drift"):
        integration.verify_historical_registry(changed, ROOT)


def test_dependency_order_drift_fails_closed(tmp_path):
    manifest = json.loads(integration.INTEGRATION.read_text(encoding="utf-8"))
    manifest["integration_order"][2]["dependencies"] = ["synthesize-vgp-program"]
    changed = tmp_path / "integration.json"
    changed.write_text(json.dumps(manifest), encoding="utf-8")
    with pytest.raises(integration.IntegrationError, match="dependency order violation"):
        integration.verify_integration_manifest(changed)
