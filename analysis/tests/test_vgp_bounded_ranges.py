import json
import subprocess
from pathlib import Path

import pytest

from analysis.vgp_10_pilot import PilotError
from analysis.vgp_bounded_ranges import (
    choose_validation_ranges,
    emit_range_bed,
    finalize_callable_masks,
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


def test_slurm_contract_queries_and_laces_only_one_bounded_range_at_a_time():
    root = Path(__file__).parents[2]
    production = (root / "analysis/slurm/run_vgp_bounded_pair.sh").read_text()
    canary = (root / "analysis/slurm/run_vgp_bounded_p07_canary.sh").read_text()
    for script in (production, canary):
        assert "emit-range-bed" in script
        assert "--target-bp 5000000 --hard-max-bp 20000000" in script
        assert "--force-large-region" in script
        assert "--sequence-files" in script
        assert '"$impg" lace' in script
        assert "global IMPG lace" not in script
        assert "emit-range-beds" not in script
        assert "partition_assignments.tsv" not in script
    assert '-l "$work/vcf.list"' in production
    assert "psmc_workers=$cpus" in production
    assert "psmc_workers > 20" not in production
    assert 'rm -rf -- "$work"' in production
    assert '"$bcftools" concat' in production
    assert "convenience genome-wide file" in production
    assert "global_partition_assignment_ledger_materialized" in production
    assert "finalize-callable" in production
    assert "for directory in results plan index" in production
    assert 'cp -a "$scratch/$directory/." "$failure/$directory/"' in production
    assert "failure_preservation.json" in production
    assert '"$bcftools" norm -f "$h1" -c s -m -any' in production
    assert '"$bcftools" norm -f "$h1" -c e -d exact' in production
    assert "ref_alt_swaps_against_exact_h1" in production


def test_pinned_bcftools_reconstructs_graph_ref_against_exact_h1(tmp_path):
    root = Path(__file__).parents[2]
    realization = json.loads(
        (root / "analysis/guix/vgp_10_pilot/realization.json").read_text()
    )
    tools = {row["name"]: row["path"] for row in realization["executables"]}
    fasta = tmp_path / "h1.fa"
    vcf = tmp_path / "graph.vcf"
    reconstructed = tmp_path / "reconstructed.vcf"
    fasta.write_text(">chr1\nACGT\n")
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=4>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t2\t.\tA\tC\t.\tPASS\t.\n"
    )
    subprocess.run([tools["samtools"], "faidx", str(fasta)], check=True)
    subprocess.run(
        [
            tools["bcftools"], "norm", "-f", str(fasta), "-c", "s", "-m", "-any",
            "-Ov", "-o", str(reconstructed), str(vcf),
        ],
        check=True,
    )
    strict = subprocess.run(
        [
            tools["bcftools"], "norm", "-f", str(fasta), "-c", "e", "-d", "exact",
            "-Ov", str(reconstructed),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    record = next(line for line in strict.stdout.splitlines() if not line.startswith("#"))
    assert record.split("\t")[3:5] == ["C", "A"]


def test_finalize_callable_masks_closes_indel_accounting_and_splits_ranges(tmp_path):
    plan = tmp_path / "plan.json"
    plan.write_text(json.dumps({"ranges": [
        {"range_id": "r0", "contig": "chr1", "start": 0, "end": 5},
        {"range_id": "r1", "contig": "chr1", "start": 5, "end": 10},
        {"range_id": "r2", "contig": "chr2", "start": 0, "end": 4},
    ]}))
    consensus = tmp_path / "consensus.fa"
    consensus.write_text(">chr1\nAANNCAAGNN\n>chr2\nNNAT\n")
    root = tmp_path / "consensus"
    (root / "masks").mkdir(parents=True)
    (root / "masks/callable.bed").write_text("chr1\t0\t10\nchr2\t0\t4\n")
    (root / "join_qc.json").write_text(json.dumps({
        "mask": {
            "callable_bp": 14, "callable_fraction": 1.0, "universe_bp": 14,
            "excluded_bp_by_primary_reason": {"not_1to1": 0},
            "reason_order": ["not_1to1"], "accounting_discrepancy_bp": 0,
        },
        "consensus": {"consensus_callable_bp": 8},
    }))
    result = finalize_callable_masks(plan, consensus, root, tmp_path / "ranges")
    assert result["final_callable_bp"] == 8
    assert result["variant_indel_flank_excluded_bp"] == 6
    assert (root / "masks/callable.prevariant.bed").is_file()
    assert (root / "masks/callable.bed").read_text().splitlines() == [
        "chr1\t0\t2", "chr1\t4\t8", "chr2\t2\t4"
    ]
    assert (tmp_path / "ranges/r0/callable.bed").read_text().splitlines() == [
        "chr1\t0\t2", "chr1\t4\t5"
    ]
    assert (tmp_path / "ranges/r1/callable.bed").read_text().strip() == "chr1\t5\t8"
    mask = json.loads((root / "masks/mask_reconciliation.json").read_text())
    assert mask["accounting_discrepancy_bp"] == 0
    assert mask["excluded_bp_by_primary_reason"]["variant_indel_flank"] == 6
