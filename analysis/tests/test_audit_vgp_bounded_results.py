import gzip
import json
from pathlib import Path

import pytest

from analysis.audit_vgp_bounded_results import (
    _intersection_bp,
    audit_callable_ownership,
    audit_range_plan,
    audit_range_variants,
)
from analysis.vgp_10_pilot import Interval
from analysis.vgp_10_pilot import PilotError


def _plan():
    return {
        "range_count": 3,
        "ranges": [
            {
                "range_id": "r0", "contig": "chr1", "start": 0, "end": 10,
                "query_required": True,
            },
            {
                "range_id": "r1", "contig": "chr1", "start": 10, "end": 20,
                "query_required": True,
            },
            {
                "range_id": "r2", "contig": "chr2", "start": 0, "end": 5,
                "query_required": False,
            },
        ],
    }


def _write_range(root: Path, range_id: str, rows: list[str]):
    path = root / f"ranges/{range_id}/normalized.vcf.gz"
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as handle:
        handle.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for row in rows:
            handle.write(row + "\n")


def test_range_variant_audit_enforces_half_open_one_owner_and_exact_reduction(tmp_path):
    _write_range(tmp_path, "r0", ["chr1\t10\t.\tA\tC\t.\tPASS\t."])
    _write_range(tmp_path, "r1", ["chr1\t11\t.\tG\tT\t.\tPASS\t."])
    (tmp_path / "range_completion.json").write_text(json.dumps({
        "normalized_variant_records": 2
    }))
    result = audit_range_variants(tmp_path, _plan())
    assert result["unique_normalized_variant_keys"] == 2
    assert result["boundary_ownership_failures"] == 0
    _write_range(tmp_path, "r1", ["chr1\t10\t.\tA\tC\t.\tPASS\t."])
    with pytest.raises(PilotError, match="two range owners"):
        audit_range_variants(tmp_path, _plan())


def test_callable_audit_sums_ranges_once_and_rejects_unowned_bases(tmp_path):
    (tmp_path / "consensus/masks").mkdir(parents=True)
    (tmp_path / "consensus/masks/callable.bed").write_text(
        "chr1\t2\t12\nchr1\t15\t20\n"
    )
    (tmp_path / "execution.json").write_text(json.dumps({
        "diversity": {"callable_bp": 15}
    }))
    result = audit_callable_ownership(tmp_path, _plan())
    assert result["callable_bp"] == 15
    assert result["unowned_callable_bp"] == 0
    (tmp_path / "consensus/masks/callable.bed").write_text("chr2\t0\t6\n")
    (tmp_path / "execution.json").write_text(json.dumps({
        "diversity": {"callable_bp": 6}
    }))
    with pytest.raises(PilotError, match="lack a range owner"):
        audit_callable_ownership(tmp_path, _plan())


def test_independent_plan_audit_recounts_exhaustive_partition_owners(tmp_path):
    plan = _plan()
    plan.update({
        "h1_total_bp": 25,
        "h1_contig_count": 2,
        "native_h1_partition_rows": 2,
        "partitioned_h1_bp": 20,
    })
    for row in plan["ranges"]:
        row["length"] = row["end"] - row["start"]
        row["partition_count"] = int(row["query_required"])
    (tmp_path / "index").mkdir()
    (tmp_path / "index/partitions.bed").write_text(
        "chr1\t0\t10\tp0\nchr1\t10\t20\tp1\nother\t0\t10\tq0\n"
    )
    result = audit_range_plan(tmp_path, plan)
    assert result["range_plan_exhaustive"] is True
    assert result["native_partition_one_owner"] is True
    plan["ranges"][1]["start"] = 9
    plan["ranges"][1]["length"] = 11
    with pytest.raises(PilotError, match="not exhaustive/disjoint"):
        audit_range_plan(tmp_path, plan)


def test_interval_intersection_uses_half_open_merged_sets():
    left = [
        Interval("chr1", 0, 10),
        Interval("chr1", 20, 30),
        Interval("chr2", 0, 5),
    ]
    right = [
        Interval("chr1", 5, 25),
        Interval("chr2", 5, 10),
    ]
    assert _intersection_bp(left, right) == 10
