import copy
import csv
import hashlib
import json
import os
import subprocess
from pathlib import Path

import jsonschema
import pytest

from analysis.vgp_10_pilot import (
    ENVIRONMENT_LOCK,
    IUPAC,
    REASON_ORDER,
    Interval,
    PilotError,
    Variant,
    assert_tool_roles,
    atomic_promote,
    audit_sweepga_paf,
    bcftools_commands,
    bootstrap_manifest,
    bootstrap_psmcfa,
    build_diploid_consensus,
    confidence_tier,
    consensus_to_psmcfa,
    construct_reason_mask,
    estimate_resources,
    freeze_bootstrap_units,
    impg_commands,
    interval_bp,
    materialize_mask_consensus_psmc,
    non_acgt_intervals,
    paf_h1_intervals,
    parse_fasta,
    parse_paf,
    parse_psmc_unscaled,
    parse_vcf,
    project_h2_non_acgt_to_h1,
    resource_prediction_errors,
    scale_unscaled_trajectory,
    select_native_partitions,
    sequence_dictionary,
    sha256_file,
    summarize_annotation_partitions,
    subtract_intervals,
    sweepga_command,
    validate_annotation_binding,
    validate_pair_input_manifest,
    validate_ref_and_reconstruct_h2,
    verify_environment_capture,
    verify_environment_lock,
)


ROOT = Path(__file__).parents[2]
SLURM = ROOT / "analysis/slurm/vgp_10_pilot"
GUIX = ROOT / "analysis/guix/vgp_10_pilot"


def write_fasta(path, sequences):
    path.write_text("".join(f">{name}\n{sequence}\n" for name, sequence in sequences.items()))


def apply(sequence, variants):
    out, cursor = [], 0
    for variant in variants:
        out.extend((sequence[cursor:variant.pos0], variant.alt))
        cursor = variant.pos0 + len(variant.ref)
    return "".join(out) + sequence[cursor:]


def truth(tmp_path):
    sequence = "ACGT" * 75
    variants = [
        Variant("h1", 10, "G", "A"),
        Variant("h1", 40, "A", "ATT"),
        Variant("h1", 70, "GTA", "G"),
    ]
    h2_sequence = apply(sequence, variants)
    h1, h2, vcf = tmp_path / "h1.fa", tmp_path / "h2.fa", tmp_path / "calls.vcf"
    write_fasta(h1, {"h1": sequence})
    write_fasta(h2, {"h2": h2_sequence})
    vcf.write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
        "".join(f"{v.contig}\t{v.pos0 + 1}\t.\t{v.ref}\t{v.alt}\t.\tPASS\t.\n" for v in variants)
    )
    return sequence, variants, h1, h2, vcf


def test_sweepga_whole_assembly_native_1to1_and_tool_role_refusals():
    command = sweepga_command("sweepga", "H1.fa", "H2.fa", "whole.paf", 7)
    assert command[:3] == ["sweepga", "H2.fa", "H1.fa"]
    assert command[command.index("--num-mappings") + 1] == "1:1"
    assert command[command.index("--scaffold-jump") + 1] == "0"
    assert_tool_roles("mapping", command)
    with pytest.raises(PilotError, match="variant caller"):
        assert_tool_roles("mapping", command + ["--output-vcf", "calls.vcf"])
    with pytest.raises(PilotError, match="SweepGA"):
        assert_tool_roles("mapping", ["wfmash", "H2.fa", "H1.fa"])
    with pytest.raises(PilotError, match="bcftools"):
        assert_tool_roles("normalize", ["impg", "norm"])


def _paf(query_start, query_end, target_start, target_end):
    length = query_end - query_start
    return (f"h2\t300\t{query_start}\t{query_end}\t+\th1\t300\t{target_start}\t{target_end}"
            f"\t{length}\t{length}\t60\tcg:Z:{length}=\n")


def test_bidirectional_paf_multiplicity_and_orientation_are_fail_closed(tmp_path):
    paf = tmp_path / "x.paf"
    paf.write_text(_paf(0, 100, 0, 100) + _paf(100, 200, 100, 200))
    audit = audit_sweepga_paf(paf, {"h1"}, {"h2"})
    assert audit["maximum_query_overlap_depth"] == audit["maximum_target_overlap_depth"] == 1

    paf.write_text(_paf(0, 100, 0, 100) + _paf(50, 150, 100, 200))
    with pytest.raises(PilotError, match="query=2, target=1"):
        audit_sweepga_paf(paf, {"h1"}, {"h2"})
    paf.write_text(_paf(0, 100, 0, 100) + _paf(100, 200, 50, 150))
    with pytest.raises(PilotError, match="query=1, target=2"):
        audit_sweepga_paf(paf, {"h1"}, {"h2"})
    paf.write_text(_paf(0, 100, 0, 100).replace("h2", "h1", 1))
    with pytest.raises(PilotError, match="orientation"):
        audit_sweepga_paf(paf, {"h1"}, {"h2"})


def test_mapping_derives_h1_1to1_complement_and_projects_h2_N_on_both_strands(tmp_path):
    plus = tmp_path / "plus.paf"
    plus.write_text("h2\t20\t2\t12\t+\th1\t30\t5\t15\t10\t10\t60\tcg:Z:10=\n")
    records = parse_paf(plus)
    assert paf_h1_intervals(records) == [Interval("h1", 5, 15)]
    assert subtract_intervals([Interval("h1", 0, 20)], paf_h1_intervals(records)) == [
        Interval("h1", 0, 5), Interval("h1", 15, 20)]
    h2 = {"h2": "AAANNN" + "A" * 14}
    assert non_acgt_intervals(h2) == [Interval("h2", 3, 6)]
    assert project_h2_non_acgt_to_h1(records, h2) == [Interval("h1", 6, 9)]

    minus = tmp_path / "minus.paf"
    minus.write_text("h2\t20\t2\t12\t-\th1\t30\t5\t15\t10\t10\t60\tcg:Z:10=\n")
    # Forward query [3,6) maps in reverse to target [11,14).
    assert project_h2_non_acgt_to_h1(parse_paf(minus), h2) == [Interval("h1", 11, 14)]

    deletion = tmp_path / "deletion.paf"
    deletion.write_text("h2\t20\t2\t10\t+\th1\t30\t5\t15\t8\t10\t60\tcg:Z:4=2D4=\n")
    assert Interval("h1", 9, 11) in project_h2_non_acgt_to_h1(
        parse_paf(deletion), {"h2": "A" * 20})
    malformed = tmp_path / "malformed.paf"
    malformed.write_text("h2\t20\t2\t12\t+\th1\t30\t5\t15\t10\t10\t60\tcg:Z:9=\n")
    with pytest.raises(PilotError, match="consumption"):
        project_h2_non_acgt_to_h1(parse_paf(malformed), {"h2": "A" * 20})


def test_exact_impg_index_partition_query_lace_order_and_both_sequences():
    commands = impg_commands(
        "impg", "exact.paf", "graph.impg", "partitions", "focus.bed", "calls",
        "vcfs.list", "laced.vcf", "H1.fa", "H2.fa", 6,
    )
    assert [row[1] for row in commands] == ["index", "partition", "query", "lace"]
    assert commands[0][commands[0].index("-a") + 1] == "exact.paf"
    assert "-w" not in commands[1]  # IMPG owns native partition defaults
    assert commands[2][commands[2].index("--sequence-files") + 1:] [:2] == ["H1.fa", "H2.fa"]
    assert commands[3][commands[3].index("--reference") + 1] == "H1.fa"
    for stage, command in zip(("index", "partition", "query", "lace"), commands):
        assert_tool_roles(stage, command)
    with pytest.raises(PilotError, match="both sequence"):
        assert_tool_roles("query", ["impg", "query", "-b", "focus.bed"])


def test_partition_focus_is_exact_subset_of_native_impg_rows(tmp_path):
    native = tmp_path / "partitions.bed"
    focus = tmp_path / "focus.bed"
    native.write_text("h1\t0\t50\tp0\nh1\t50\t100\tp1\nh1\t100\t150\tp2\n")
    selected = select_native_partitions(native, [Interval("h1", 70, 120)], focus)
    assert selected == [(Interval("h1", 50, 100), "p1"), (Interval("h1", 100, 150), "p2")]
    assert focus.read_text().splitlines() == ["h1\t50\t100\tp1", "h1\t100\t150\tp2"]
    malformed = tmp_path / "user-windows.bed"
    malformed.write_text("h1\t0\t10\n")
    with pytest.raises(PilotError, match="partition identifier"):
        select_native_partitions(malformed, [Interval("h1", 0, 10)], focus)


def test_bcftools_owns_decompose_trim_exact_dedup_and_both_indexes():
    commands = bcftools_commands("bcftools", "laced.vcf", "H1.fa", "one.bed", "normalized")
    assert [row[1] for row in commands] == ["norm", "view", "norm", "index", "norm", "index"]
    assert commands[0][commands[0].index("-m") + 1] == "-any"
    assert commands[1][commands[1].index("-R") + 1] == "one.bed"
    assert "exact" in commands[2] and "-Oz" in commands[2]
    assert "-t" in commands[3]
    assert "exact" in commands[4] and "-Ob" in commands[4]


def test_normalized_indels_parse_and_exact_duplicates_are_rejected(tmp_path):
    _, expected, _, _, vcf = truth(tmp_path)
    assert parse_vcf(vcf) == expected
    vcf.write_text(vcf.read_text() + "h1\t41\t.\tA\tATT\t.\tPASS\t.\n")
    with pytest.raises(PilotError, match="exact duplicate"):
        parse_vcf(vcf)


def test_reason_mask_order_disjoint_complement_and_exact_denominator():
    universe = [Interval("h1", 0, 100)]
    flags = {
        "h2_gap_or_N": [Interval("h1", 20, 50)],
        "not_eligible_contig": [Interval("h1", 0, 10), Interval("h1", 30, 40)],
        "repeat_or_low_complexity_primary": [Interval("h1", 45, 60)],
    }
    callable_rows, exclusions, reconciliation = construct_reason_mask(universe, flags)
    assert exclusions["not_eligible_contig"] == [Interval("h1", 0, 10), Interval("h1", 30, 40)]
    assert exclusions["h2_gap_or_N"] == [Interval("h1", 20, 30), Interval("h1", 40, 50)]
    assert exclusions["repeat_or_low_complexity_primary"] == [Interval("h1", 50, 60)]
    assert interval_bp(callable_rows) == 50
    assert reconciliation["universe_bp"] == 100
    assert reconciliation["callable_bp"] + sum(reconciliation["excluded_bp_by_primary_reason"].values()) == 100
    assert reconciliation["reason_order"] == list(REASON_ORDER)
    with pytest.raises(PilotError, match="unregistered"):
        construct_reason_mask(universe, {"made_up": [Interval("h1", 1, 2)]})


def test_h1_ref_and_alternate_sequence_concordance_with_snps_and_indels(tmp_path):
    _, variants, h1_path, h2_path, _ = truth(tmp_path)
    h1, h2 = parse_fasta(h1_path), parse_fasta(h2_path)
    audit = validate_ref_and_reconstruct_h2(h1, h2, variants, {"h1": "h2"})
    assert audit["h1_ref_mismatches"] == audit["h2_reconstruction_failures"] == 0
    invalid = list(variants)
    invalid[0] = Variant("h1", 10, "C", "A")
    with pytest.raises(PilotError, match="H1 REF mismatch"):
        validate_ref_and_reconstruct_h2(h1, h2, invalid, {"h1": "h2"})
    bad_h2 = dict(h2)
    bad_h2["h2"] = bad_h2["h2"][:-1] + "A"
    with pytest.raises(PilotError, match="reconstruction failed"):
        validate_ref_and_reconstruct_h2(h1, bad_h2, variants, {"h1": "h2"})


def test_regional_h2_reconstruction_is_manifest_bound(tmp_path):
    sequence, variants, h1_path, h2_path, _ = truth(tmp_path)
    h1, h2 = parse_fasta(h1_path), parse_fasta(h2_path)
    regions = [{"h1_contig": "h1", "h1_start": 0, "h1_end": len(sequence),
                "h2_contig": "h2", "h2_start": 0, "h2_end": len(h2["h2"]), "strand": "+"}]
    assert validate_ref_and_reconstruct_h2(h1, h2, variants, aligned_regions=regions)[
        "h2_reconstruction_failures"] == 0
    with pytest.raises(PilotError, match="not every normalized variant"):
        validate_ref_and_reconstruct_h2(h1, h2, variants, aligned_regions=[{
            "h1_contig": "h1", "h1_start": 0, "h1_end": 20,
            "h2_contig": "h2", "h2_start": 0, "h2_end": 20, "strand": "+"}])


def test_diploid_consensus_iupac_indel_mask_and_noncallable_N(tmp_path):
    sequence, variants, h1_path, _, _ = truth(tmp_path)
    h1 = parse_fasta(h1_path)
    callable_rows = [Interval("h1", 5, 90), Interval("h1", 100, 290)]
    consensus, qc = build_diploid_consensus(h1, callable_rows, variants, indel_flank=10)
    assert consensus["h1"][10] == IUPAC[frozenset(("G", "A"))]
    assert set(consensus["h1"][:5]) == {"N"}
    assert set(consensus["h1"][90:100]) == {"N"}
    assert set(consensus["h1"][30:51]) == {"N"}  # insertion anchor +/- 10
    assert set(consensus["h1"][60:83]) == {"N"}  # deletion REF span +/- 10
    assert qc["non_callable_bases_encoded_homozygous_reference"] == 0
    with pytest.raises(PilotError, match="callable"):
        build_diploid_consensus(h1, [Interval("h1", 20, 290)], variants)


def test_psmcfa_preserves_mask_and_heterozygous_bins():
    value = consensus_to_psmcfa({"a": "A" * 100 + "R" + "A" * 99 + "N" + "A" * 99})
    assert value == ">a\nTKN\n"


def test_at_least_100_deterministic_boundary_aware_bootstraps():
    callable_rows = [Interval("a", 0, 12), Interval("a", 20, 31), Interval("b", 0, 9)]
    units = freeze_bootstrap_units(callable_rows, 5)
    assert units == [
        Interval("a", 0, 5), Interval("a", 5, 10), Interval("a", 10, 12),
        Interval("a", 20, 25), Interval("a", 25, 30), Interval("a", 30, 31),
        Interval("b", 0, 5), Interval("b", 5, 9),
    ]
    rows = bootstrap_manifest(callable_rows, "P07", "a" * 64, attempts=200, block_bp=5)
    assert len(rows) == 200
    assert rows == bootstrap_manifest(callable_rows, "P07", "a" * 64, attempts=200, block_bp=5)
    assert all(len(row["sampled_unit_indices"]) == len(units) for row in rows)
    with pytest.raises(PilotError, match="at least 100"):
        bootstrap_manifest(callable_rows, "P07", "a" * 64, attempts=99, block_bp=5)
    encoded = bootstrap_psmcfa({"a": "A" * 31, "b": "C" * 9}, units, rows[0]["sampled_unit_indices"], 5)
    assert encoded.count(">block_") == len(units)


def test_unscaled_trajectory_and_scaling_scenarios_remain_separate():
    source = [{"interval": 0, "time_2N0": 0.5, "lambda": 2.0}]
    scenarios = [{"scenario_id": "s1", "mutation_rate_per_generation": 1e-8,
                  "generation_time_years": 5, "theta_0": 0.001,
                  "mutation_rate_source": "doi:mu", "generation_time_source": "doi:g"}]
    unscaled, scaled = scale_unscaled_trajectory(source, scenarios)
    assert unscaled == [{"interval": 0, "time_2N0": 0.5, "lambda": 2.0}]
    assert "scenario_id" not in unscaled[0]
    assert scaled[0]["scenario_id"] == "s1"
    assert scaled[0]["mutation_rate_source"] == "doi:mu"


def test_psmc_parser_selects_final_unscaled_iteration_and_theta(tmp_path):
    output = tmp_path / "truth.psmc"
    output.write_text(
        "TR\t0\t0.001\t0.01\nRS\t0\t0\t0.1\t1.0\n"
        "TR\t25\t0.002\t0.02\nRS\t25\t1\t0.5\t2.0\nRS\t25\t0\t0.2\t1.5\n"
    )
    rows, theta = parse_psmc_unscaled(output)
    assert theta == 0.002
    assert rows == [
        {"interval": 0, "time_2N0": 0.2, "lambda": 1.5},
        {"interval": 1, "time_2N0": 0.5, "lambda": 2.0},
    ]
    output.write_text("RS\t25\t0\t0.2\t1.5\n")
    with pytest.raises(PilotError, match="theta"):
        parse_psmc_unscaled(output)


def test_exact_annotation_accession_dictionary_or_validated_liftover_only(tmp_path):
    dictionary = [{"name": "h1", "length": 300, "md5": "a" * 32}]
    native = {"assembly_accession_version": "GCA_000001.1", "sequence_dictionary": dictionary,
              "annotation_accession_version": "GCA_000001.1-ANN_1", "gff_sha256": "f" * 64}
    assert validate_annotation_binding("GCA_000001.1", dictionary, native)["status"] == "exact_native"
    assert validate_annotation_binding("GCA_000001.1", dictionary, None) == {
        "status": "not_available", "core_eligible": True, "annotation_outputs_allowed": False}
    wrong = copy.deepcopy(native)
    wrong["sequence_dictionary"][0]["length"] = 301
    with pytest.raises(PilotError, match="dictionary"):
        validate_annotation_binding("GCA_000001.1", dictionary, wrong)
    with pytest.raises(PilotError, match="liftover"):
        validate_annotation_binding("GCA_000001.1", dictionary, {
            "assembly_accession_version": "GCA_000002.1", "sequence_dictionary": dictionary,
            "annotation_accession_version": "GCA_000002.1-ANN_1", "gff_sha256": "f" * 64})
    lifted = {"assembly_accession_version": "GCA_000002.1",
              "annotation_accession_version": "GCA_000002.1-ANN_1", "gff_sha256": "f" * 64,
              "validated_liftover": {
        "source_accession_version": "GCA_000002.1", "target_accession_version": "GCA_000001.1",
        "chain_sha256": "1" * 64, "validation_sha256": "2" * 64,
        "manifest_sha256": "3" * 64, "passed": True}}
    assert validate_annotation_binding("GCA_000001.1", dictionary, lifted)["status"] == "validated_liftover"


def test_annotation_cds_fourfold_effect_ws_sw_and_gc3_partitions():
    callable_rows = [Interval("h1", 0, 100)]
    variants = [Variant("h1", 10, "A", "G"), Variant("h1", 20, "G", "A")]
    features = {name: [Interval("h1", 0, 50)] for name in
                ("CDS", "fourfold", "nonsynonymous", "synonymous", "WS", "SW", "GC3")}
    rows = {row["partition"]: row for row in
            summarize_annotation_partitions(callable_rows, variants, features)}
    assert set(rows) == {"CDS", "fourfold", "nonsynonymous", "synonymous", "WS", "SW", "GC3"}
    assert rows["WS"]["variant_records"] == 1
    assert rows["SW"]["variant_records"] == 1
    assert all(row["callable_h1_bp"] == 50 for row in rows.values())


def test_pair_input_provenance_core_qc_and_selective_validation_policy(tmp_path):
    import analysis.vgp_10_pilot as pilot
    row = pilot.load_primary_pair("P01")
    h1 = tmp_path / "h1.fa"
    h2 = tmp_path / "h2.fa"
    write_fasta(h1, {"h1": "A" * 20})
    write_fasta(h2, {"h2": "A" * 20})
    def asset(path):
        return {"path": str(path), "sha256": sha256_file(path), "size_bytes": path.stat().st_size,
                "sequence_dictionary": sequence_dictionary(parse_fasta(path))}
    value = {key: row[key] for key in
             ("biosample", "individual_or_isolate", "h1_accession_version", "h2_accession_version")}
    value.update({"selection_id": "P01", "orientation": "H1_reference_H2_query",
                  "assets": {"h1_fasta": asset(h1), "h2_fasta": asset(h2)},
                  "read_technology_resolved": True, "long_range_phasing_evidence": "Hi-C",
                  "selective_validation": {},
                  "core_qc": {hap: {"qv": 45, "busco_complete_fraction": .95,
                                    "busco_missing_fraction": .02, "busco_duplicated_fraction": .01,
                                    "copy_number_and_kmer_audit_passed": True} for hap in ("h1", "h2")}})
    result = validate_pair_input_manifest(value)
    assert result["selective_validation_is_universal_core_gate"] is False
    bad = copy.deepcopy(value)
    bad["h1_accession_version"] = "GCA_000000001.1"
    with pytest.raises(PilotError, match="provenance"):
        validate_pair_input_manifest(bad)
    bad = copy.deepcopy(value)
    bad["core_qc"]["h2"]["qv"] = 39
    with pytest.raises(PilotError, match="QV"):
        validate_pair_input_manifest(bad)


def test_confidence_tiers_track_core_and_selective_evidence_without_universal_gate():
    hard = {key: True for key in ("exact_pair", "qv_pass", "completeness_pass", "collapse_pass",
            "mapping_1to1_pass", "callability_pass", "consensus_pass", "reproducibility_pass")}
    assert confidence_tier(hard) == "B"
    complete = dict(hard, raw_read_validation=True, kmer_validation=True,
                    published_estimate_validation=True, long_range_switch_validation=True)
    assert confidence_tier(complete) == "A"
    incomplete = dict(hard)
    incomplete["qv_pass"] = None
    assert confidence_tier(incomplete) == "C"
    failed = dict(hard)
    failed["callability_pass"] = False
    assert confidence_tier(failed) == "X"


def test_resource_estimator_uses_measured_size_and_has_no_global_ceiling():
    small = estimate_resources(100, 100, 100_000_000, 100_000_000, 10, 100)
    large = estimate_resources(1000, 1000, 4_000_000_000, 4_000_000_000, 5000, 10000)
    assert large["map_cpu_hours_estimate"] > small["map_cpu_hours_estimate"]
    assert large["scratch_bytes_estimate"] > small["scratch_bytes_estimate"]
    assert large["read_bytes_estimate"] > small["read_bytes_estimate"]
    assert large["map_wall_hours_estimate"] == large["map_cpu_hours_estimate"] / 8
    assert "no global" in large["scheduler_limits"]
    with pytest.raises(PilotError, match="positive sequence"):
        estimate_resources(1, 1, 0, 1, 1)
    errors = resource_prediction_errors(
        {key: 10 for key in ("wall_time", "cpu_hours", "peak_rss", "scratch_high_water", "read_bytes", "write_bytes")},
        {key: 8 for key in ("wall_time", "cpu_hours", "peak_rss", "scratch_high_water", "read_bytes", "write_bytes")},
    )
    assert set(errors.values()) == {0.25}


def test_atomic_promotion_sentinel_is_idempotent(tmp_path):
    partial, final = tmp_path / "stage.partial", tmp_path / "stage"
    partial.mkdir()
    (partial / "result.txt").write_text("complete")
    atomic_promote(partial, final, {"stage": "truth"})
    assert json.loads((final / ".complete.json").read_text())["files"]["result.txt"] == hashlib.sha256(b"complete").hexdigest()
    atomic_promote(tmp_path / "absent", final, {"stage": "truth"})
    broken = tmp_path / "broken"
    broken.mkdir()
    with pytest.raises(PilotError, match="sentinel"):
        atomic_promote(tmp_path / "absent2", broken, {})


def test_guix_lock_pins_channels_manifests_sources_submodules_and_fails_unrealized():
    lock = verify_environment_lock(require_realized=False)
    assert lock["channel_commit"] == "44bbfc24e4bcc48d0e3343cd3d83452721af8c36"
    assert lock["source_identities"]["psmc"]["commit"] == "b37b1cfa05b89c67c2ad1b63c699a27600d5516e"
    assert set(lock["source_identities"]["impg"]["submodules"]) == {"gfaffix", "syng"}
    with pytest.raises(PilotError, match="not been captured"):
        verify_environment_lock(require_realized=True)
    assert "guix time-machine" in (GUIX / "capture_environment.sh").read_text()
    capture = verify_environment_capture(GUIX / "realization.json")
    assert capture["reproducibility"]["psmc_guix_check"] == "passed"
    assert len(capture["executables"]) == 16


def test_materialized_join_proves_concordance_mask_consensus_and_200_bootstraps(tmp_path):
    sequence, _, h1, h2, vcf = truth(tmp_path)
    universe = tmp_path / "universe.bed"
    universe.write_text(f"h1\t0\t{len(sequence)}\n")
    exclusion = tmp_path / "not1.bed"
    exclusion.write_text("h1\t290\t300\n")
    output = tmp_path / "join"
    qc = materialize_mask_consensus_psmc(
        h1, h2, vcf, universe, {"not_1to1": exclusion}, output,
        contig_map={"h1": "h2"}, selection_id="P07", attempts=200,
    )
    assert qc["bootstrap_attempts"] == 200
    assert qc["concordance"]["h2_reconstruction_failures"] == 0
    assert (output / "masks/callable.bed").is_file()
    manifest = list(csv.DictReader((output / "consensus/bootstrap_manifest.tsv").open(), delimiter="\t"))
    assert len(manifest) == 200
    assert (output / "consensus/bootstrap_units.1mb.bed").is_file()
    assert (output / "consensus/bootstrap_units.10mb.bed").is_file()


def test_output_schema_is_valid_and_annotation_absence_is_nonblocking():
    schema = json.loads((ROOT / "analysis/vgp_10_pilot_output_schema.json").read_text())
    jsonschema.Draft202012Validator.check_schema(schema)
    annotation_schema = schema["properties"]["annotation"]
    jsonschema.Draft202012Validator(annotation_schema).validate({"status": "not_available", "core_blocked": False})
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.Draft202012Validator(annotation_schema).validate({"status": "not_available", "core_blocked": True})


def test_slurm_entrypoints_are_resumable_atomic_ordered_and_have_no_global_memory_ceiling():
    pair = (SLURM / "pair_stage.sh").read_text()
    common = (SLURM / "common.sh").read_text()
    submit = (SLURM / "submit.sh").read_text()
    psmc = (SLURM / "psmc_array.sh").read_text()
    for token in ("--num-mappings 1:1", "impg\" index", "impg\" partition", "impg\" query", "impg\" lace",
                  "bcftools\" norm", "-d exact", "materialize_mask_consensus_psmc"):
        assert token in pair
    assert pair.index('"$impg" index') < pair.index('"$impg" partition') < pair.index('"$impg" query') < pair.index('"$impg" lace')
    assert ".complete.json" in common and "atomic_promote" in common and "record_telemetry" in common
    assert "--dry-run" in submit and "--submit" in submit and "resources.json" in submit
    assert "#SBATCH --mem" not in pair and "#SBATCH --mem" not in psmc
    assert "emit-bootstrap" in psmc and "SLURM_ARRAY_TASK_ID" in psmc
    assert "--array=0-200%20" in submit and "psmc_finalize.sh" in submit
    for script in SLURM.glob("*.sh"):
        subprocess.run(["bash", "-n", str(script)], check=True)


def test_cli_plan_is_h1_reference_h2_query_and_has_exact_stage_order(tmp_path):
    result = subprocess.run(
        [os.environ.get("PYTHON", "python3"), "-m", "analysis.vgp_10_pilot", "plan", "P01", str(tmp_path)],
        cwd=ROOT, check=True, capture_output=True, text=True,
    )
    value = json.loads(result.stdout)
    assert value["orientation"] == "H1_reference_H2_query"
    assert [row[1] for row in value["commands"]["impg"]] == ["index", "partition", "query", "lace"]
