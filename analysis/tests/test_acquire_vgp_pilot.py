import csv
import json
from pathlib import Path

from analysis import acquire_vgp_pilot as acquire


ROOT = Path(__file__).parents[2]


def _read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_run_writes_refusal_manifest_and_report_for_current_nogo_gate(tmp_path):
    manifest_out = tmp_path / "acquisition_manifest.tsv"
    report_out = tmp_path / "acquisition_report.md"
    result = acquire.run(
        gate_path=ROOT / "analysis/vgp_pilot_gate.json",
        manifest_path=ROOT / "analysis/vgp_pilot_manifest.tsv",
        root_config_path=ROOT / "analysis/vgp_data_root_config.json",
        output_manifest_path=manifest_out,
        output_report_path=report_out,
    )
    assert result["status"] == "refused_preflight"
    assert result["failure_code"] == "GATE_NO_GO"
    rows = _read_manifest(manifest_out)
    assert rows[0]["record_type"] == "run_summary"
    assert rows[0]["status"] == "refused_preflight"
    assert rows[0]["observed_bytes"] == "0"
    blocker_codes = {row["failure_code"] for row in rows[1:]}
    assert "NO_SELECTED_ROWS" in blocker_codes
    assert "ZERO_COMPOSITION_ELIGIBLE_ROWS" in blocker_codes
    report = report_out.read_text(encoding="utf-8")
    assert "Gate decision: `NO_GO`" in report
    assert "Refused before first biological byte: `true`" in report
    assert "Exact-reference/native-annotation linkage validation under pinned GNU Guix was not re-run" in report


def test_run_refuses_tampered_gate_and_records_gate_tampered_failure(tmp_path):
    gate_out = tmp_path / "tampered_gate.json"
    gate_payload = json.loads((ROOT / "analysis/vgp_pilot_gate.json").read_text(encoding="utf-8"))
    gate_payload["decision"]["status"] = "GO"
    gate_out.write_text(json.dumps(gate_payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    manifest_out = tmp_path / "acquisition_manifest.tsv"
    report_out = tmp_path / "acquisition_report.md"
    result = acquire.run(
        gate_path=gate_out,
        manifest_path=ROOT / "analysis/vgp_pilot_manifest.tsv",
        root_config_path=ROOT / "analysis/vgp_data_root_config.json",
        output_manifest_path=manifest_out,
        output_report_path=report_out,
    )
    assert result["status"] == "refused_preflight"
    assert result["failure_code"] == "GATE_TAMPERED"
    rows = _read_manifest(manifest_out)
    assert rows[0]["failure_code"] == "GATE_TAMPERED"
    assert "gate decision hash does not match the gate payload" in rows[0]["failure_message"]
