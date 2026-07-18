import csv
import json
from pathlib import Path

import pytest

from analysis import review_vgp_10_pilot as review


def rows(path):
    with Path(path).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_review_reconciles_all_slots_and_refuses_scaleout(tmp_path):
    decision = review.review(tmp_path, tmp_path / "independent")
    assert decision["program_decision"] == "CONDITIONAL_GO"
    assert decision["full_biological_scaleout_authorized"] is False
    assert decision["primary_passes"] == 0
    assert decision["primary_slots_reviewed"] == 10
    assert decision["alternates_reviewed"] == 6
    assert decision["alternate_activations"] == 0
    assert decision["slurm_jobs_submitted_by_review"] == 0
    assert decision["classification_counts"] == {
        "technical_failure": 10,
        "low_confidence_usable_core": 0,
        "biological_outlier": 0,
        "biological_not_estimable": 10,
    }
    scaleout = rows(tmp_path / "vgp_10_pilot_scaleout_manifest.tsv")
    assert len(scaleout) == 16
    assert {row["biological_jobs_authorized"] for row in scaleout} == {"false"}
    assert {row["full_scaleout_authorized"] for row in scaleout} == {"false"}
    assert sum(row["roster_type"] == "primary" for row in scaleout) == 10
    assert sum(row["roster_type"] == "alternate" for row in scaleout) == 6


def test_every_hard_and_quantitative_gate_is_explicit(tmp_path):
    review.review(tmp_path, tmp_path / "independent")
    gates = {row["gate_id"]: row for row in rows(tmp_path / "vgp_10_pilot_review_gates.tsv")}
    assert set(("H01", "H02", "H03", "H04", "H05", "H06", "H07", "H08", "H09")) <= set(gates)
    assert gates["P01"]["threshold"] == ">=8/10"
    assert gates["P01"]["status"] == "FAIL"
    assert gates["P02"]["observed"] == "0/6"
    assert gates["P03"]["observed"] == "0/2"
    assert gates["B01"]["threshold"] == ">=100"
    assert gates["B02"]["threshold"] == ">=0.95"
    assert gates["R01"]["threshold"] == "<=0.25"
    assert gates["R02"]["threshold"] == "<=0.50"
    assert gates["R01"]["status"] == gates["R02"]["status"] == "FAIL_NOT_ESTIMABLE"
    assert gates["N01"]["status"] == gates["N02"]["status"] == "PASS"


def test_predeclared_subset_is_live_rehashed_and_nonmaterialization_compared(tmp_path):
    decision = review.review(tmp_path, tmp_path / "independent")
    validation = json.loads((tmp_path / "independent/subset_validation.json").read_text())
    assert validation["subset"] == ["P07", "P08"]
    assert validation["immutable_objects_revalidated"] == 20
    assert validation["checksum_or_size_failures"] == 0
    assert validation["biological_quantities_compared"] == 16
    assert validation["biological_recomputation_status"] == "NOT_ELIGIBLE_PREFLIGHT_FAILED"
    assert validation["slurm_jobs_submitted"] == 0
    assert decision["independent_validation"] == validation
    comparisons = rows(tmp_path / "independent/subset_comparison.tsv")
    assert len(comparisons) == 16
    assert {row["comparison"] for row in comparisons} == {"MATCH_NOT_MATERIALIZED_PRECHECK_FAILED"}


def test_branch_decisions_keep_specialized_state_orthogonal(tmp_path):
    decision = review.review(tmp_path, tmp_path / "independent")
    branches = decision["branches"]
    assert branches["core_diversity_psmc"]["decision"] == "CONDITIONAL_GO"
    assert branches["exact_annotation_partitions"]["decision"] == "CONDITIONAL_GO"
    assert branches["direct_conversion"]["decision"] == "CONDITIONAL_GO"
    assert branches["phylogenetic_substitution_bias"]["decision"] == "CONDITIONAL_GO"
    assert branches["population_gbgc"]["decision"] == "NOT_RUN/DESIGN_ONLY"
    assert branches["non_allelic_conversion"]["decision"] == "NOT_RUN/DESIGN_ONLY"
    assert all(value["scaleout_authorized"] is False for value in branches.values())
    assert decision["interpretation_guards"]["annotation_absence_vetoes_core"] is False
    assert decision["interpretation_guards"]["same_pair_psmc_is_independent_validation"] is False


def test_review_rejects_provenance_or_terminal_promotion():
    design = review.read_tsv(review.ANALYSIS / "vgp_10_pair_manifest.tsv")
    acquisitions = review.read_tsv(review.ANALYSIS / "vgp_10_pilot_acquisition_manifest.tsv")
    results = review.read_tsv(review.ANALYSIS / "vgp_10_pilot_result_manifest.tsv")
    qc = review.read_tsv(review.ANALYSIS / "vgp_10_pilot_qc.tsv")
    telemetry = review.read_tsv(review.ANALYSIS / "vgp_10_pilot_resource_telemetry.tsv")
    summary = json.loads((review.ANALYSIS / "vgp_10_pilot_run_summary.json").read_text())
    changed = [dict(row) for row in results]
    changed[0]["h1_accession_version"] = "GCA_000000000.0"
    with pytest.raises(review.ReviewError, match="provenance drift"):
        review.validate_closed_ledgers(design, acquisitions, changed, qc, telemetry, summary)
    changed = [dict(row) for row in results]
    changed[0]["terminal_state"] = "COMPLETED"
    with pytest.raises(review.ReviewError, match="promoted"):
        review.validate_closed_ledgers(design, acquisitions, changed, qc, telemetry, summary)


def test_resource_estimate_is_observed_where_possible_and_marks_unknown_model():
    telemetry = review.read_tsv(review.ANALYSIS / "vgp_10_pilot_resource_telemetry.tsv")
    value = review.resource_assessment(telemetry)
    assert value["observed_preflight"]["core_objects"] == 100
    assert value["observed_preflight"]["logical_bytes"] == 12_567_760_437
    assert value["observed_preflight"]["logical_write_bytes"] == 20_793
    assert value["resource_model_gate"]["median_ape"] is None
    assert value["resource_model_gate"]["p95_ape"] is None
    assert value["resource_model_gate"]["status"] == "FAIL_NOT_ESTIMABLE_ZERO_CLUSTER_JOBS"
    envelope = value["upper_bound_716_pair_biological_planning_envelope"]
    assert envelope["low"]["core_hours"] < envelope["base"]["core_hours"] < envelope["high"]["core_hours"]
    assert envelope["high"]["scratch_gb_per_job"] == 37.5
    assert envelope["high"]["scratch_gb_aggregate"] == 375
    assert value["headroom"]["storage_and_inode_multiplier"] == 1.25
    assert value["headroom"]["per_job_stop_multiple_of_reviewed_high_estimate"] == 1.5
