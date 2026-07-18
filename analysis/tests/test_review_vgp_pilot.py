from pathlib import Path

from analysis import review_vgp_pilot as review


ROOT = Path(__file__).parents[2]


def _read_tsv(path: Path) -> list[dict[str, str]]:
    return review.load_tsv(path)


def test_review_recomputes_current_refusal_artifacts(tmp_path):
    review_out = tmp_path / "review.md"
    qc_out = tmp_path / "qc.tsv"
    resource_out = tmp_path / "resource.tsv"

    result = review.review(
        gate_path=ROOT / "analysis/vgp_pilot_gate.json",
        run_manifest_path=ROOT / "analysis/vgp_pilot_run_manifest.tsv",
        results_path=ROOT / "analysis/vgp_pilot_results.tsv",
        telemetry_path=ROOT / "analysis/vgp_pilot_slurm_telemetry.tsv",
        review_out=review_out,
        qc_out=qc_out,
        resource_out=resource_out,
        guix_validation="INCONCLUSIVE",
        guix_note="test run skipped the external Guix command",
    )

    assert result["overall_decision"] == "FAIL"
    assert review_out.is_file()
    assert qc_out.is_file()
    assert resource_out.is_file()

    qc_rows = _read_tsv(qc_out)
    decisions = {row["check_id"]: row["decision"] for row in qc_rows}
    assert decisions["promoted_gate_recompute"] == "PASS"
    assert decisions["run_manifest_recompute"] == "FAIL"
    assert decisions["telemetry_recompute"] == "PASS"
    assert decisions["results_recompute"] == "FAIL"
    assert decisions["source_catalog_counts"] == "PASS"
    assert decisions["selected_manifest_rows"] == "FAIL"
    assert decisions["quota_interface"] == "FAIL"
    assert decisions["slurm_terminal_state"] == "PASS"
    assert decisions["network_in_arrays"] == "PASS"
    assert decisions["guix_analysis_suite"] == "INCONCLUSIVE"

    resource_rows = _read_tsv(resource_out)
    by_metric = {row["metric"]: row for row in resource_rows}
    assert by_metric["aggregate_core_hours"]["decision"] == "PASS"
    assert by_metric["aggregate_core_hours"]["authorized_cap"] == "280"
    assert by_metric["aggregate_core_hours"]["observed"] == "0"
    assert by_metric["peak_memory_gib_per_job"]["decision"] == "INCONCLUSIVE"

    review_text = review_out.read_text(encoding="utf-8")
    assert "Overall decision: `FAIL`" in review_text
    assert "No Slurm array job ID" in review_text
