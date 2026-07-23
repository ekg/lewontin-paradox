import json
from pathlib import Path

from analysis.materialize_vgp_bounded_evidence import render


def test_bounded_report_requires_and_renders_three_actual_results(tmp_path: Path):
    pairs = []
    for selection_id in ("P07", "P03", "P02"):
        pairs.append({
            "selection_id": selection_id,
            "species": f"species-{selection_id}",
            "failure_class": f"class-{selection_id}",
            "diversity": {"heterozygous_snps": 10, "callable_bp": 1000, "pi": 0.01},
            "bounded_range_plan": {"range_count": 3},
            "psmc": {"finite_bootstraps": 200, "primary_theta_centered": True},
            "annotation": {
                "annotation_status": "exact_native" if selection_id == "P07"
                else "missing_nonblocking"
            },
            "graph_identifier_audit": {"unresolved_ids": 0, "silently_omitted_regions": 0},
            "range_variant_audit": {"exact_duplicate_keys_between_ranges": 0},
            "callable_ownership_audit": {
                "unowned_callable_bp": 0, "multiply_owned_callable_bp": 0
            },
            "independent_stratified_reaudit": [
                {
                    "stratum": label, "range_id": f"r{number}", "contig": "chr1",
                    "start": number * 100, "end": (number + 1) * 100,
                    "variant_records": number + 1, "callable_bp": 90,
                    "consensus_non_N_bp": 90,
                    "psmc_population_bins": {"N": 1, "K": 2, "T": 3},
                }
                for number, label in enumerate(("early", "middle", "late"))
            ],
        })
    execution = tmp_path / "execution.json"
    execution.write_text(json.dumps({
        "actual_core_biological_results": 3,
        "pairs": pairs,
        "canceled_global_jobs": ["1", "2"],
        "controlled_fastga_wfmash_comparison": {
            "overlapping_target_bp": 100, "target_coverage_jaccard": 0.9,
            "exact_variant_jaccard": 0.5,
        },
        "remaining_pipeline_limitations": ["technical limit"],
        "selection_freeze_sha256": "a" * 64,
    }))
    transition = tmp_path / "transition.json"
    transition.write_text(json.dumps({"bounded_canary": {
        "contig": "chr1", "start": 0, "end": 100, "range_bp": 100,
        "native_partition_count": 5, "normalized_variant_keys": 4,
        "exact_variant_key_sha256": "b" * 64, "callable_bp": 90,
        "peak_local_graph_state_bytes": 1234,
    }}))
    sacct = tmp_path / "sacct.tsv"
    sacct.write_text(
        "JobIDRaw\tJobName\tState\tElapsed\tAllocCPUS\tMaxRSS\tNodeList\n"
        "1\tbounded\tCOMPLETED\t00:01:00\t2\t100K\tnode\n"
    )
    report = tmp_path / "report.md"
    render(execution, transition, sacct, report)
    text = report.read_text()
    assert "Three actual same-individual" in text
    assert "no all-genome graph" in text
    assert "P07" in text and "P03" in text and "P02" in text
    assert "200/200 finite; centered=true" in text
    assert "unowned callable bp=0" in text
    assert "technical limit" in text
