import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
ANALYSIS = ROOT / "analysis"


def load(name: str) -> dict:
    return json.loads((ANALYSIS / name).read_text())


def test_fresh_canary_reproduction_records_bounded_equivalence_and_cleanup():
    value = load("vgp_three_pair_review_canary_reproduction_v1.json")

    assert value["range"] == {
        "contig": "CM106587.1",
        "start": 0,
        "end": 5_000_000,
        "native_partitions": 2_500,
    }
    assert value["regional_vcf_count"] == 2_500
    assert value["normalized_variant_keys"] == 3_871
    assert value["normalized_variant_key_sha256"] == (
        "b7f64cd018a693ade59de7daec49b66d"
        "9f57e6ea7baa888c6e9e1c7099da1405"
    )
    assert value["callable_bp"] == 3_116_900
    assert value["matches_clean_impg_subset"] is True
    assert value["local_graph_temporaries_disposed"] is True
    assert value["global_graph_materialized"] is False
    assert value["global_impg_lace_used"] is False


def test_corrected_resource_model_separates_requests_measurements_and_recovery():
    value = load("vgp_three_pair_review_resource_model_v1.json")
    jobs = {row["job_id"]: row for row in value["measured_jobs"]}
    pairs = {row["pair"]: row for row in value["observed_pair_scale"]}

    assert set(jobs) == {
        "1797782",
        "1799306",
        "1799307",
        "1799308",
        "1804024",
        "1804979",
    }
    assert all(row["max_rss"] is None for row in jobs.values())
    assert value["corrected_interpretation"]["memory_is_request_not_usage"] is True
    assert pairs["P02"]["accepted_job_or_lineage"] == "1804024+1804979"
    assert pairs["P02"]["accepted_lineage_elapsed_hours"] == 22.999167
    assert pairs["P02"]["accepted_lineage_allocated_cpu_hours"] == 735.973333
    assert value["architecture"]["aggregate_variant_policy"].startswith(
        "bcftools concat"
    )
    assert value["required_scale_out_telemetry"] == {
        "capture_cgroup_peak_rss_per_stage": True,
        "capture_scratch_high_water_and_file_count_per_stage": True,
        "capture_user_system_and_allocated_cpu_seconds_per_stage": True,
        "record_filesystem_type_and_resolved_scratch_path": True,
        "retain_failed_lineage_cost_separately_from_steady_state_cost": True,
    }


def test_review_issues_four_decisions_and_classifies_every_defect_family():
    text = (ANALYSIS / "vgp_three_pair_independent_review_v1.md").read_text()

    for heading in (
        "## Fresh bounded P07 reproduction",
        "## Recomputed core measurements",
        "## IMPG architecture audit",
        "## Annotation partitions",
        "## PSMC and repaired centering",
        "## Mapping fallback equivalence",
        "## Read sensitivity",
        "## Provenance and scratch",
        "## Decisions",
        "## Remaining defect ledger",
        "## Corrected resource handoff",
    ):
        assert heading in text

    for decision in (
        "| Core assembly inference | **ACCEPT** |",
        "| Annotation partitions | **ACCEPT P07 ONLY** |",
        "| Read sensitivity | **INCONCLUSIVE** |",
        "| Broader scale-out | **NO-GO for biological inference;",
    ):
        assert decision in text

    for defect_prefix in ("| I-1 | implementation |", "| D-1 | data provenance |",
                          "| A-1 | assembly confidence |", "| B-1 | biological limitation |"):
        assert defect_prefix in text

    assert "no symmetric H1/H2/graph comparison exists" in text
    assert "not treated as a generic truth oracle" in text
    assert "16,133 promoted" in text
    assert "files totaling 12,179,399,820 bytes" in text
