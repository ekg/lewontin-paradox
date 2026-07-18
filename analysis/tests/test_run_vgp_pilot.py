import csv
import json
from pathlib import Path

from analysis import run_vgp_pilot as runner


ROOT = Path(__file__).parents[2]


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_run_writes_refusal_outputs_for_current_nogo_gate(tmp_path):
    run_manifest = tmp_path / "run_manifest.tsv"
    slurm_telemetry = tmp_path / "slurm_telemetry.tsv"
    results = tmp_path / "results.tsv"
    result = runner.run(
        gate_path=ROOT / "analysis/vgp_pilot_gate.json",
        manifest_path=ROOT / "analysis/vgp_pilot_manifest.tsv",
        root_config_path=ROOT / "analysis/vgp_data_root_config.json",
        sweepga_build_path=ROOT / "analysis/sweepga_origin_main_build.json",
        impg_handoff_path=ROOT / "analysis/sweepga_impg_observed.json",
        output_run_manifest_path=run_manifest,
        output_slurm_telemetry_path=slurm_telemetry,
        output_results_path=results,
    )
    assert result["status"] == "refused_preflight"
    # Preserve the integrated refusal gate rather than rebinding it in the
    # resolver task; compute rejects the repaired manifest digest until regate.
    assert result["failure_code"] == "MANIFEST_DIGEST_MISMATCH"

    manifest_rows = _read_tsv(run_manifest)
    assert manifest_rows[0]["record_type"] == "run_summary"
    assert manifest_rows[0]["status"] == "refused_preflight"
    assert manifest_rows[0]["gate_decision"] == "NO_GO"
    assert manifest_rows[0]["manifest_digest"]
    assert manifest_rows[0]["cap_vector_digest"]
    assert manifest_rows[0]["sweepga_origin_build_sha256"]
    assert manifest_rows[0]["impg_handoff_sha256"]
    blocker_codes = {row["failure_code"] for row in manifest_rows[1:]}
    assert "NO_SELECTED_ROWS" in blocker_codes
    assert "ZERO_COMPOSITION_ELIGIBLE_ROWS" in blocker_codes
    assert "ZERO_DIVERSITY_ELIGIBLE_ROWS" in blocker_codes

    telemetry_rows = _read_tsv(slurm_telemetry)
    assert telemetry_rows == [
        {
            "run_id": telemetry_rows[0]["run_id"],
            "generated_at_utc": telemetry_rows[0]["generated_at_utc"],
            "record_type": "run_summary",
            "status": "refused_preflight",
            "candidate_id": "",
            "sbatch_command": "",
            "slurm_job_id": "",
            "slurm_array_job_id": "",
            "dependency": "",
            "requested_cpus": "",
            "requested_memory_gib": "",
            "requested_wall_hours": "",
            "requested_scratch_gb": "",
            "requested_read_gb": "",
            "requested_write_gb": "",
            "retry_index": "0",
            "final_state": "NOT_SUBMITTED",
            "max_rss_gib": "",
            "elapsed_seconds": "0",
            "cpu_time_seconds": "0",
            "scratch_peak_gb": "0",
            "io_read_gb": "0",
            "io_write_gb": "0",
            "metadata_operations": "0",
            "failure_code": "MANIFEST_DIGEST_MISMATCH",
            "notes": "gate refusal prevented sbatch submission and sacct telemetry collection",
        }
    ]

    result_rows = _read_tsv(results)
    assert result_rows[0]["metric"] == "validated_species_count"
    assert result_rows[0]["value"] == "0"
    assert result_rows[0]["failure_code"] == "MANIFEST_DIGEST_MISMATCH"
    assert {row["failure_code"] for row in result_rows[1:]} >= {
        "NO_SELECTED_ROWS",
        "ZERO_COMPOSITION_ELIGIBLE_ROWS",
        "ZERO_DIVERSITY_ELIGIBLE_ROWS",
    }


def test_run_refuses_tampered_gate_and_preserves_audited_outputs(tmp_path):
    gate_path = tmp_path / "tampered_gate.json"
    gate_payload = json.loads((ROOT / "analysis/vgp_pilot_gate.json").read_text(encoding="utf-8"))
    gate_payload["decision"]["status"] = "GO"
    gate_path.write_text(json.dumps(gate_payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    run_manifest = tmp_path / "run_manifest.tsv"
    slurm_telemetry = tmp_path / "slurm_telemetry.tsv"
    results = tmp_path / "results.tsv"
    result = runner.run(
        gate_path=gate_path,
        manifest_path=ROOT / "analysis/vgp_pilot_manifest.tsv",
        root_config_path=ROOT / "analysis/vgp_data_root_config.json",
        sweepga_build_path=ROOT / "analysis/sweepga_origin_main_build.json",
        impg_handoff_path=ROOT / "analysis/sweepga_impg_observed.json",
        output_run_manifest_path=run_manifest,
        output_slurm_telemetry_path=slurm_telemetry,
        output_results_path=results,
    )
    assert result["status"] == "refused_preflight"
    assert result["failure_code"] == "GATE_TAMPERED"
    manifest_rows = _read_tsv(run_manifest)
    assert manifest_rows[0]["failure_code"] == "GATE_TAMPERED"
    assert "gate decision hash does not match the gate payload" in manifest_rows[0]["failure_message"]
    telemetry_rows = _read_tsv(slurm_telemetry)
    assert telemetry_rows[0]["failure_code"] == "GATE_TAMPERED"
    result_rows = _read_tsv(results)
    assert result_rows[0]["failure_code"] == "GATE_TAMPERED"
