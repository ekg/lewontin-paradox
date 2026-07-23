import json
from pathlib import Path

import pytest

from analysis.vgp_10_pilot import PilotError
from analysis.vgp_bounded_ranges import (
    choose_validation_ranges,
    emit_range_bed,
    freeze_range_plan,
)


def _fixture(tmp_path: Path):
    fai = tmp_path / "h1.fa.fai"
    fai.write_text("chr1\t12000\t0\t0\t0\nchr2\t4500\t0\t0\t0\n")
    partitions = tmp_path / "partitions.bed"
    partitions.write_text(
        "chr1\t0\t2000\tp0\n"
        "h2a\t0\t1900\tp0\n"
        "chr1\t2000\t4000\tp1\n"
        "chr1\t4000\t6000\tp2\n"
        "h2a\t1900\t5800\tp2\n"
        "chr1\t6000\t8000\tp3\n"
        "chr1\t8000\t10000\tp4\n"
        "chr1\t10000\t12000\tp5\n"
        "chr2\t0\t2000\tp6\n"
        "chr2\t2000\t4000\tp7\n"
        "chr2\t4000\t4500\tp8\n"
    )
    return fai, partitions


def test_freeze_plan_is_native_boundary_aligned_disjoint_exhaustive_and_one_owner(tmp_path):
    fai, partitions = _fixture(tmp_path)
    plan_path = tmp_path / "plan.json"
    plan = freeze_range_plan(
        fai, partitions, plan_path, selection_id="P02",
        target_bp=5000, hard_max_bp=8000,
    )
    assert plan["h1_total_bp"] == plan["assigned_bp"] == 16500
    assert plan["native_h1_partition_rows"] == plan["assigned_native_partition_rows"] == 9
    assert plan["range_plan_disjoint"] is True
    assert plan["range_plan_exhaustive"] is True
    assert plan["native_partition_one_owner"] is True
    assert plan["global_impg_lace_created"] is False
    assert [(row["contig"], row["start"], row["end"]) for row in plan["ranges"]] == [
        ("chr1", 0, 4000), ("chr1", 4000, 8000), ("chr1", 8000, 12000),
        ("chr2", 0, 4500),
    ]
    assert json.loads(plan_path.read_text())["global_partition_assignment_ledger_materialized"] is False


def test_freeze_plan_rejects_gap_overlap_and_unbounded_native_partition(tmp_path):
    fai, partitions = _fixture(tmp_path)
    plan = tmp_path / "plan.json"
    partitions.write_text("chr1\t0\t2000\tp0\nchr1\t3000\t12000\tp1\nchr2\t0\t4500\tp2\n")
    with pytest.raises(PilotError, match="gap"):
        freeze_range_plan(fai, partitions, plan, selection_id="P03")
    partitions.write_text("chr1\t0\t12000\tp0\nchr2\t0\t4500\tp1\n")
    with pytest.raises(PilotError, match="exceeds hard"):
        freeze_range_plan(
            fai, partitions, plan, selection_id="P03",
            target_bp=5000, hard_max_bp=10000,
        )


def test_emit_one_range_bed_preserves_exact_frozen_partition_census(tmp_path):
    fai, partitions = _fixture(tmp_path)
    plan_path = tmp_path / "plan.json"
    plan = freeze_range_plan(
        fai, partitions, plan_path, selection_id="P07",
        target_bp=5000, hard_max_bp=8000,
    )
    bed = tmp_path / "r000000.bed"
    result = emit_range_bed(plan_path, partitions, "r000000", bed)
    assert result == {
        "range_id": "r000000", "contig": "chr1", "start": 0, "end": 4000,
        "partition_count": 2, "query_required": True,
        "only_requested_range_materialized": True, "one_owner_census_passed": True,
    }
    assert bed.read_text().splitlines() == [
        "chr1\t0\t2000\tquery000000000\tp0",
        "chr1\t2000\t4000\tquery000000001\tp1",
    ]
    assert choose_validation_ranges(plan) == ["r000000", "r000002", "r000003"]
