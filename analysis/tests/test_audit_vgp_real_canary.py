import json

from analysis.audit_vgp_real_canary import (
    IntervalIndex,
    independent_mask_reconstruction,
    maximum_axis_depth,
    parse_bed,
    read_regional_vcf_audit,
    summarize_scaled_scenarios,
    summarize_vcf,
)
from analysis.vgp_real_canary_annotation import annotation_partition_counts
from analysis.promote_vgp_real_canary import verify_stage_sentinels


def test_independent_mask_reconstruction_preserves_first_reason_and_subset(tmp_path):
    universe = {"chr1": 20}
    flags = {
        "first": [("chr1", 2, 8)],
        "second": [("chr1", 5, 12)],
    }
    result = independent_mask_reconstruction(universe, ("first", "second"), flags)
    assert result["callable"] == [("chr1", 0, 2), ("chr1", 12, 20)]
    assert result["by_reason"]["first"] == [("chr1", 2, 8)]
    assert result["by_reason"]["second"] == [("chr1", 8, 12)]
    assert result["callable_bp"] == 10
    assert result["accounting_discrepancy_bp"] == 0


def test_independent_mask_sweep_and_interval_index_handle_fragmented_boundaries():
    universe = {"chr1": 30, "chr2": 5}
    flags = {
        "first": [("chr1", -2, 4), ("chr1", 4, 8), ("chr1", 20, 40)],
        "second": [("chr1", 2, 6), ("chr1", 10, 25)],
    }
    result = independent_mask_reconstruction(universe, ("first", "second"), flags)
    assert result["by_reason"]["first"] == [("chr1", 0, 8), ("chr1", 20, 30)]
    assert result["by_reason"]["second"] == [("chr1", 10, 20)]
    assert result["callable"] == [("chr1", 8, 10), ("chr2", 0, 5)]
    assert result["accounting_discrepancy_bp"] == 0

    index = IntervalIndex(result["callable"])
    assert index.contains("chr1", 8, 10)
    assert index.contains("chr2", 0, 5)
    assert not index.contains("chr1", 7, 9)
    assert not index.contains("chr1", 9, 11)


def test_paf_depth_and_vcf_subset_are_reconstructed_without_pipeline_helpers(tmp_path):
    paf = tmp_path / "mapping.paf"
    paf.write_text(
        "q1\t20\t0\t10\t+\tt1\t20\t0\t10\t10\t10\t60\n"
        "q1\t20\t10\t20\t+\tt1\t20\t10\t20\t10\t10\t60\n"
    )
    rows = [line.split("\t") for line in paf.read_text().splitlines()]
    assert maximum_axis_depth(rows, "query") == 1
    assert maximum_axis_depth(rows, "target") == 1

    vcf = tmp_path / "calls.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t4\t.\tA\tG\t.\tPASS\t.\n"
        "chr1\t15\t.\tA\tAT\t.\tPASS\t.\n"
    )
    summary = summarize_vcf(vcf, [("chr1", 0, 10)])
    assert summary["normalized_variant_records"] == 2
    assert summary["normalized_snp_records"] == 1
    assert summary["callable_variant_records"] == 1
    assert summary["callable_snp_records"] == 1
    assert summary["callable_indel_records"] == 0
    assert summary["subset_records"] == [("chr1", 4, "A", "G")]


def test_exact_annotation_partitions_keep_callable_denominators_and_ws_sw_direction():
    rows = annotation_partition_counts(
        [("chr1", 0, 12)],
        [("chr1", 0, 12)],
        {("chr1", 2), ("chr1", 5)},
        {"chr1": "AATCCGAAAAAA"},
        [("chr1", 2, "T", "G"), ("chr1", 5, "G", "A")],
    )
    by_name = {row["partition"]: row for row in rows}
    assert by_name["CDS"]["callable_bp"] == 12
    assert by_name["fourfold"]["callable_bp"] == 2
    assert by_name["WS"]["heterozygous_variants"] == 1
    assert by_name["SW"]["heterozygous_variants"] == 1


def test_atomic_promotion_verifier_rehashes_stage_payloads(tmp_path):
    stage = tmp_path / "mapping"
    stage.mkdir()
    payload = stage / "mapping.paf"
    payload.write_text("biological\n")
    import hashlib
    digest = hashlib.sha256(payload.read_bytes()).hexdigest()
    (stage / ".complete.json").write_text(json.dumps({"files": {"mapping.paf": digest}}))
    audit = verify_stage_sentinels(tmp_path)
    assert audit["stage_sentinels"] == 1
    assert audit["verified_payload_files"] == 1
    payload.write_text("drift\n")
    import pytest
    with pytest.raises(RuntimeError, match="digest mismatch"):
        verify_stage_sentinels(tmp_path)


def test_regional_vcf_audit_replaces_transient_absolute_path_list(tmp_path):
    path = tmp_path / "regional_vcf_audit.json"
    path.write_text(json.dumps({
        "canonical_vgp_root": "/moosefs/erikg/vgp",
        "focus_rows": 203698,
        "unique_query_names": 203698,
        "unique_native_partition_ids": 142198,
        "regional_vcf_count": 203698,
        "regional_vcf_total_bytes": 987654,
        "all_regional_vcfs_nonempty": True,
        "transient_shards_removed_after_lacing": True,
    }))
    value = read_regional_vcf_audit(path, 203698, "/moosefs/erikg/vgp")
    assert value["regional_vcf_count"] == 203698
    value["regional_vcf_count"] = 1
    path.write_text(json.dumps(value))
    import pytest
    with pytest.raises(RuntimeError, match="census"):
        read_regional_vcf_audit(path, 203698, "/moosefs/erikg/vgp")


def test_scaled_scenario_audit_requires_labels_sources_and_frozen_bin_size():
    rows = [{
        "scenario_id": "SENS_MU1E-8_G2Y",
        "mutation_rate_per_generation": "1e-8",
        "generation_time_years": "2",
        "psmc_bin_size_bp": "100",
        "mutation_rate_source": "sensitivity_not_calibration",
        "generation_time_source": "sensitivity_not_calibration",
    }]
    result = summarize_scaled_scenarios(rows)
    assert result["scenario_ids"] == ["SENS_MU1E-8_G2Y"]
    assert result["psmc_bin_size_bp"] == 100
    rows[0]["psmc_bin_size_bp"] = "1"
    import pytest
    with pytest.raises(RuntimeError, match="100-bp"):
        summarize_scaled_scenarios(rows)
