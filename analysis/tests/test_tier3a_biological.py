import csv
import math
from pathlib import Path

from analysis.tier3a_biological import (
    DIVERSITY_COLUMNS,
    bootstrap_ratio,
    intersect_by_contig,
    merge_by_contig,
    select_execution_panel,
    validate_diversity_rows,
)


ROOT = Path(__file__).resolve().parents[2]


def test_execution_panel_is_deterministic_and_not_input_order_dependent():
    spans = [
        {"execution_span_id": f"span-{i}", "contig": "chr1", "start_0based": str(i * 100), "end_0based_exclusive": str(i * 100 + 80)}
        for i in range(20)
    ]
    first = select_execution_panel("biological-tuple", spans, 500)
    second = select_execution_panel("biological-tuple", list(reversed(spans)), 500)
    assert [row["execution_span_id"] for row in first] == [row["execution_span_id"] for row in second]
    assert sum(int(row["end_0based_exclusive"]) - int(row["start_0based"]) for row in first) >= 500
    assert len(first) < len(spans)


def test_callable_intersection_merges_overlap_and_preserves_half_open_bounds():
    left = {"chr1": [(0, 10), (8, 20)], "chr2": [(5, 9)]}
    right = {"chr1": [(3, 7), (12, 30)], "chr2": [(0, 6)]}
    assert merge_by_contig(left) == {"chr1": [(0, 20)], "chr2": [(5, 9)]}
    assert intersect_by_contig(left, right) == {
        "chr1": [(3, 7), (12, 20)],
        "chr2": [(5, 6)],
    }


def test_block_bootstrap_has_finite_uncertainty_and_is_seeded():
    blocks = [(100, 10), (200, 12), (150, 11), (250, 17)]
    first = bootstrap_ratio(blocks, replicates=500, seed=17)
    second = bootstrap_ratio(blocks, replicates=500, seed=17)
    assert first == second
    assert first["replicates"] == 500
    assert 0 < first["ci_low"] <= first["estimate"] <= first["ci_high"]
    assert math.isfinite(first["standard_error"])


def test_diversity_validator_requires_positive_biological_rows():
    rows = [{
        "dataset_id": "tuple-1",
        "scientific_name": "Biologica realis",
        "annotation_class": "CDS",
        "eligible_haplotypes": "2",
        "variant_numerator": "3",
        "callable_denominator": "1000",
        "estimate": "0.003",
        "bootstrap_ci_low": "0.001",
        "bootstrap_ci_high": "0.005",
        "bootstrap_standard_error": "0.0007",
        "biological_input": "yes",
    }]
    validate_diversity_rows(rows, {"tuple-1"})
    bad = [dict(rows[0], callable_denominator="0")]
    try:
        validate_diversity_rows(bad, {"tuple-1"})
    except ValueError as error:
        assert "callable" in str(error)
    else:
        raise AssertionError("zero callable denominator was accepted")


def test_diversity_schema_names_reference_conditioned_statistics():
    assert "statistic_label" in DIVERSITY_COLUMNS
    assert DIVERSITY_COLUMNS.index("annotation_class") < DIVERSITY_COLUMNS.index("statistic_label")


def test_committed_biological_results_cover_every_acquired_tuple():
    with (ROOT / "results/tier3a/acquisition_manifest.tsv").open(encoding="utf-8", newline="") as handle:
        acquired = {
            row["dataset_id"] for row in csv.DictReader(handle, delimiter="\t")
            if row["eligibility_status"] == "eligible_biological"
        }
    with (ROOT / "results/tier3a/diploid_diversity.tsv").open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    validate_diversity_rows(rows, acquired)
    assert len(rows) == 4 * len(acquired)
    by_dataset = {dataset: {row["statistic_label"] for row in rows if row["dataset_id"] == dataset} for dataset in acquired}
    assert all(labels == {
        "diploid_haplotype_diversity",
        "pi_S_reference_conditioned",
        "pi_W_reference_conditioned",
    } for labels in by_dataset.values())
    assert all(row["scope_label"].endswith("not_genome_wide") for row in rows)


def test_committed_qc_and_manifest_preserve_separate_tool_provenance():
    qc = (ROOT / "results/tier3a/diploid_qc.md").read_text(encoding="utf-8")
    assert "SweepGA supplied complete whole-H1-versus-H2 mappings" in qc
    assert "IMPG selected" in qc
    assert "both allele checks passed" in qc
    assert "Caps 5/10 are coverage sensitivity only" in qc
    with (ROOT / "results/tier3a/diploid_run_manifest.tsv").open(encoding="utf-8", newline="") as handle:
        runs = list(csv.DictReader(handle, delimiter="\t"))
    assert len(runs) == 3
    assert all(row["status"] == "completed" and row["sweepga_hit_cap"] == "1:1" for row in runs)
    assert all(row["impg_commit"] and row["impg_native_partitions_selected"] for row in runs)
    assert all(row["guix_channel_commit"] and row["guix_profile_store_path"].startswith("/gnu/store/") for row in runs)
