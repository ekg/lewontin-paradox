from __future__ import annotations

import csv
import tempfile
from pathlib import Path

from analysis import synthesize_repaired_vgp as synth


ROOT = Path(__file__).parents[2]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_repaired_synthesis_preserves_refusal_and_exact_demography_join() -> None:
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        summary = synth.build_synthesis(
            root=ROOT,
            synthesis_out=tmp / "synthesis.md",
            paper_table_out=tmp / "paper.tsv",
            decision_out=tmp / "decision.tsv",
        )

        assert summary == {
            "candidate_count": 6,
            "executed_candidate_count": 0,
            "valid_independent_ne_count": 6,
            "valid_independent_ne_species_count": 1,
            "paper_row_count": 6,
        }

        paper = read_tsv(tmp / "paper.tsv")
        assert len(paper) == 6
        assert {row["exact_join_status"] for row in paper} == {"PASS_EXACT_IDENTITY"}
        assert {row["candidate_execution_status"] for row in paper} == {
            "NOT_EXECUTED_NO_GO"
        }
        assert {row["resolved_modality"] for row in paper} == {"tier3c_composition"}
        assert {row["metadata_composition_eligible"] for row in paper} == {"yes"}
        assert {row["metadata_diversity_eligible"] for row in paper} == {"no"}
        assert {row["composition_measurement_status"] for row in paper} == {
            "NOT_MEASURED_REFUSED"
        }
        assert {row["diversity_measurement_status"] for row in paper} == {
            "NOT_MEASURED_INELIGIBLE_AND_REFUSED"
        }
        assert {row["psmc_eligible"] for row in paper} == {"no"}
        assert {row["msmc2_eligible"] for row in paper} == {"no"}
        assert {row["smcpp_eligible"] for row in paper} == {"no"}

        camel = next(row for row in paper if row["scientific_name"] == "Camelus dromedarius")
        assert camel["valid_independent_ne_record_count"] == "6"
        assert camel["valid_independent_ne_populations"] == (
            "Awarik;Haddana;Majaheem;Sahliah;Shul;Sofor"
        )
        assert camel["valid_independent_ne_values"] == "15;11;37;24;17;23"
        assert camel["valid_independent_ne_interval_types"] == ";".join(
            ["not_applicable"] * 6
        )
        assert camel["valid_independent_ne_uncertainty"] == ";".join(
            ["interval not reported"] * 6
        )
        assert camel["historical_scenario_record_ids"] == "camel_psmc_fitak2020"
        assert camel["historical_scenario_mutation_rates"] == (
            "1.1e-8 changes/site/generation"
        )
        assert camel["historical_scenario_generation_times"] == "5 years"

        shark = next(row for row in paper if row["scientific_name"] == "Heterodontus francisci")
        assert shark["valid_independent_ne_record_count"] == "0"
        assert shark["coalescent_scaled_record_ids"] == "horn_theta_ima3_2022"
        assert shark["coalescent_scaled_units"] == "theta (4Ne-mu)"
        assert "not absolute" in shark["coalescent_scaled_disposition"]

        assert all(row["circular_estimate_record_count"] == "1" for row in paper)
        assert all(row["h1_h2_demography_disposition"].startswith("assembly haplotypes only") for row in paper)

        decisions = read_tsv(tmp / "decision.tsv")
        assert {row["option_id"] for row in decisions} == {
            "repair_remaining_candidate_metadata",
            "request_bounded_expansion_wave",
            "request_population_data_subset",
            "stop",
        }
        assert {row["creates_ready_executable_task"] for row in decisions} == {"no"}
        assert {row["authorization_granted"] for row in decisions} == {"no"}

        projections = [row for row in decisions if row["row_type"] == "resource_projection"]
        assert {row["projection_scope"] for row in projections} == {
            "next_wave",
            "full_eligible_catalog",
        }
        assert all(row["low"] == row["base"] == row["high"] == "" for row in projections)
        assert {row["calibration_status"] for row in projections} == {
            "NOT_ESTIMABLE_NO_SUCCESSFUL_OBSERVATION"
        }

        text = (tmp / "synthesis.md").read_text(encoding="utf-8")
        assert "NO_GO" in text
        assert "zero-byte" in text
        assert "NOT_SUBMITTED" in text
        assert "empty-result" in text
        assert "No biological pilot ran" in text
        assert "Lewontin" in text
        assert "not performance telemetry" in text
        assert "No attributable active jobs or download processes" in text
