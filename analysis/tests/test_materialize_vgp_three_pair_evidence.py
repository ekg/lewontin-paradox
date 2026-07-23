import csv
import json
from pathlib import Path

import pytest

from analysis.materialize_vgp_three_pair_evidence import (
    materialize_independent_reaudit,
    materialize_report,
    materialize_telemetry,
)
from analysis.vgp_10_pilot import PilotError


def pair(selection_id: str) -> dict:
    return {
        "selection_id": selection_id,
        "species": "Species example",
        "failure_class": "example_failure",
        "actual_core_biological_result": True,
        "diversity": {"heterozygous_snps": 10, "callable_bp": 1000, "pi": 0.01},
        "psmc": {"finite_bootstraps": 200},
        "annotation": {"annotation_status": "exact_native"},
        "independent_stratified_reaudit": [{"stratum": value} for value in ("early", "middle", "late")],
        "coordinate_and_strand_audit": {"invalid_coordinates": 0},
        "ref_alt_reconstruction": {"reconstruction_failures": 0},
        "normalized_concordance": {"h2_reconstruction_failures": 0},
        "graph_identifier_audit": {"unresolved_ids": 0},
        "mask": {"accounting_discrepancy_bp": 0},
        "psmc_population_preserved": True,
        "resource_telemetry": [{
            "selection_id": selection_id, "stage": "consensus", "job_id": "123",
            "disposition": "success", "started_epoch": 100, "ended_epoch": 125,
        }],
    }


@pytest.fixture
def execution(tmp_path: Path) -> tuple[dict, Path]:
    value = {
        "task_id": "run-vgp-three-pair", "run_id": "test", "actual_core_biological_results": 3,
        "completion_gate_passed": True, "pairs": [pair(value) for value in ("P07", "P03", "P02")],
        "controlled_fastga_wfmash_comparison": {
            "fastga_target_bp": 90, "wfmash_target_bp": 95, "overlapping_target_bp": 88,
            "target_coverage_jaccard": 0.907, "exact_variant_jaccard": 0.25,
        },
        "remaining_pipeline_limitations": ["Example limitation."],
    }
    path = tmp_path / "execution.json"
    path.write_text(json.dumps(value))
    return value, path


def test_materializes_independent_and_stage_telemetry(execution, tmp_path: Path):
    value, path = execution
    independent, telemetry = tmp_path / "independent.json", tmp_path / "telemetry.tsv"
    result = materialize_independent_reaudit(value, path, independent)
    assert result["all_three_stratified_reaudits_present"] is True
    assert len(json.loads(independent.read_text())["pairs"]) == 3
    assert materialize_telemetry(value, telemetry) == 3
    with telemetry.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert {row["selection_id"] for row in rows} == {"P07", "P03", "P02"}
    assert {row["elapsed_seconds"] for row in rows} == {"25"}


def test_report_names_three_results_and_backend_limitation(execution, tmp_path: Path):
    value, path = execution
    sacct = tmp_path / "sacct.tsv"
    sacct.write_text(
        "JobIDRaw|JobName|State|ExitCode|Elapsed|AllocCPUS|ReqMem|NodeList\n"
        "123|vgp-P03|COMPLETED|0:0|00:01:00|32|160G|node1\n"
    )
    report = tmp_path / "report.md"
    materialize_report(value, path, sacct, report)
    text = report.read_text()
    assert "three actual core biological results" in text
    assert all(selection_id in text for selection_id in ("P07", "P03", "P02"))
    assert "gap placement" in text and "not biological exclusions" in text


def test_refuses_fewer_than_three_results(execution, tmp_path: Path):
    value, path = execution
    value["pairs"].pop()
    value["actual_core_biological_results"] = 2
    with pytest.raises(PilotError, match="three actual"):
        materialize_independent_reaudit(value, path, tmp_path / "bad.json")
