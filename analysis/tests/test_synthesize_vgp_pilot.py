from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
import sys

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from analysis import synthesize_vgp_pilot as synth


ROOT = Path(__file__).parents[2]


def _read_tsv(path: Path) -> list[dict[str, str]]:
    return synth.load_tsv(path)


class SynthesizeVgpPilotTest(unittest.TestCase):
    def test_synthesis_rebuilds_current_fail_closed_handoff(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            synthesis_out = tmp_path / "vgp_pilot_synthesis.md"
            paper_table_out = tmp_path / "vgp_pilot_paper_table.tsv"
            next_decision_out = tmp_path / "vgp_pilot_next_decision.tsv"

            result = synth.build_synthesis(
                manifest_path=ROOT / "analysis/vgp_pilot_manifest.tsv",
                rejections_path=ROOT / "analysis/vgp_pilot_rejections.tsv",
                results_path=ROOT / "analysis/vgp_pilot_results.tsv",
                telemetry_path=ROOT / "analysis/vgp_pilot_slurm_telemetry.tsv",
                qc_path=ROOT / "analysis/vgp_pilot_qc.tsv",
                resource_path=ROOT / "analysis/vgp_pilot_resource_calibration.tsv",
                ne_sources_path=ROOT / "analysis/vgp_pilot_ne_sources.tsv",
                availability_path=ROOT / "analysis/vgp_pilot_population_data_availability.tsv",
                budget_path=ROOT / "analysis/vertebrate_scaleout_resource_budget.tsv",
                synthesis_out=synthesis_out,
                paper_table_out=paper_table_out,
                next_decision_out=next_decision_out,
            )

            self.assertEqual(result["manifest_rows"], 6)
            self.assertEqual(result["paper_rows"], 6)
            self.assertEqual(result["selected_count"], 0)
            self.assertEqual(result["matched_ne_rows"], 0)
            self.assertEqual(result["matched_availability_rows"], 0)

            paper_rows = _read_tsv(paper_table_out)
            self.assertEqual(len(paper_rows), 6)
            self.assertEqual({row["ne_source_rows"] for row in paper_rows}, {"0"})
            self.assertEqual({row["population_availability_rows"] for row in paper_rows}, {"0"})
            self.assertEqual(
                {row["join_status"] for row in paper_rows},
                {"no_selected_pilot_species_no_inventory_rows"},
            )

            next_rows = _read_tsv(next_decision_out)
            recommendation_rows = [row for row in next_rows if row["row_type"] == "recommendation"]
            self.assertEqual(recommendation_rows[0]["recommendation"], "stop_repair")

            current_core = next(
                row
                for row in next_rows
                if row["scenario_id"] == "current_no_go_executable_cost" and row["metric"] == "aggregate_core_hours"
            )
            self.assertEqual(current_core["base"], "0")

            next_wave = next(
                row
                for row in next_rows
                if row["scenario_id"] == "contingent_repaired_next_wave_leq6"
                and row["metric"] == "aggregate_core_hours"
            )
            self.assertEqual(next_wave["calibration_status"], "observed_historical_proxy")
            self.assertGreater(float(next_wave["base"]), 0.0)

            synthesis_text = synthesis_out.read_text(encoding="utf-8")
            self.assertIn("Recommended next decision: `stop_repair`.", synthesis_text)
            self.assertIn("validated executable species: `0`", synthesis_text)
            self.assertIn(
                "the bounded pilot produced **no** promoted cross-species diversity or composition estimate.",
                synthesis_text,
            )


if __name__ == "__main__":
    unittest.main()
