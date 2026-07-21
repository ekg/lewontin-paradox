from __future__ import annotations

import io
import json
import math
from pathlib import Path

import pytest

from analysis.vgp_read_validation import (
    DepthMaskSpec,
    assembly_evidence,
    estimate_kmer_heterozygosity,
    kmer_qv,
    parse_psmc,
    parse_pileup_bases,
    stream_depth_masks,
    summarize_paf,
    summarize_site_concordance,
)


def test_depth_masks_are_half_open_and_split_on_coordinate_gaps(tmp_path: Path) -> None:
    stream = io.StringIO(
        "chr1\t1\t8\n"
        "chr1\t2\t12\n"
        "chr1\t3\t20\n"
        "chr1\t5\t20\n"
        "chr2\t1\t12\n"
    )
    summary = stream_depth_masks(
        stream,
        tmp_path,
        [DepthMaskSpec("dp10_15", 10, 15), DepthMaskSpec("dp10_30", 10, 30)],
    )

    assert (tmp_path / "dp10_15.bed").read_text() == "chr1\t1\t2\nchr2\t0\t1\n"
    assert (tmp_path / "dp10_30.bed").read_text() == (
        "chr1\t1\t3\nchr1\t4\t5\nchr2\t0\t1\n"
    )
    assert summary["observed_positions"] == 5
    assert summary["depth_histogram"] == {"8": 1, "12": 2, "20": 2}
    assert summary["masks"]["dp10_30"]["callable_bp"] == 4
    assert summary["depth_structure"]["modal_positive_depth"] == 12
    assert summary["depth_structure"]["below_half_mode_fraction"] == 0
    assert summary["depth_structure"]["above_1_5x_mode_fraction"] == pytest.approx(0.4)


def test_pileup_parser_handles_markers_indels_and_strands() -> None:
    # Four REF observations, four ALT observations, one other base.  The +2tt
    # insertion and start/end markers annotate observations and are not bases.
    counts = parse_pileup_bases(".,Aa^F.,+2ttAa-1a$g", "C", "A")
    assert counts == {"ref": 4, "alt": 4, "other": 1, "deletion": 0}


def test_assembly_evidence_classifies_supported_and_contradicted_sites() -> None:
    assembly = io.StringIO("chr1\t10\tA\tG\nchr1\t20\tC\tT\n")
    pileup = io.StringIO(
        "chr1\t10\tA\t10\t.....GGGGG\tFFFFFFFFFF\n"
        "chr1\t20\tC\t10\t..........\tFFFFFFFFFF\n"
    )
    out = io.StringIO()
    summary = assembly_evidence(assembly, pileup, out, minimum_depth=10, maximum_depth=80)

    assert summary["supported_heterozygous"] == 1
    assert summary["contradicted_homozygous_reference"] == 1
    assert summary["not_observed"] == 0
    assert summary["concrete_false_positive_lower_bound_sites"] == 1
    assert summary["concrete_false_positive_lower_bound_fraction"] == pytest.approx(0.5)
    assert summary["candidate_false_positive_upper_bound_sites"] == 1
    rows = out.getvalue().splitlines()
    assert rows[0].startswith("chrom\tposition_1based")
    assert rows[1].endswith("supported_heterozygous")
    assert rows[2].endswith("contradicted_homozygous_reference")


def test_kmer_qv_uses_per_base_error_probability_and_wilson_bound() -> None:
    result = kmer_qv(total_kmers=1_000_000, error_kmers=2_100, k=21)
    expected_error = 1.0 - (1.0 - 0.0021) ** (1.0 / 21.0)
    assert result["error_rate"] == pytest.approx(expected_error)
    assert result["qv"] == pytest.approx(-10.0 * math.log10(expected_error))
    assert result["qv_lower_95"] < result["qv"] < result["qv_upper_95"]


def test_kmer_heterozygosity_recovers_synthetic_two_copy_spectrum() -> None:
    pytest.importorskip("scipy")
    # Synthetic Poisson mixture: 36 M distinct heterozygous k-mers at 22x,
    # 390 M homozygous k-mers at 44x, and smaller repeat components.
    from scipy.stats import poisson

    hist = {}
    for depth in range(1, 180):
        value = (
            36_000_000 * poisson.pmf(depth, 22)
            + 390_000_000 * poisson.pmf(depth, 44)
            + 12_000_000 * poisson.pmf(depth, 66)
            + 8_000_000 * poisson.pmf(depth, 88)
        )
        hist[depth] = int(round(value))
    result = estimate_kmer_heterozygosity(hist, k=21, minimum_depth=5)
    expected = 1.0 - (2 * 390 / (36 + 2 * 390)) ** (1 / 21)

    assert result["status"] == "estimated"
    assert result["homozygous_peak_depth"] == pytest.approx(44, abs=1.0)
    assert result["heterozygosity_per_base"] == pytest.approx(expected, rel=0.12)
    assert result["fit_r_squared"] > 0.995


def test_site_concordance_preserves_common_denominator_and_covariance() -> None:
    assembly = {
        ("chr1", 10, "A", "G"),
        ("chr1", 20, "C", "T"),
        ("chr1", 30, "G", "A"),
    }
    reads = {
        ("chr1", 10, "A", "G"),
        ("chr1", 20, "C", "T"),
        ("chr1", 40, "T", "C"),
        ("chr1", 50, "A", "C"),
    }
    result = summarize_site_concordance(assembly, reads, callable_bp=1_000)

    assert result["shared_sites"] == 2
    assert result["assembly_only_sites"] == 1
    assert result["read_only_sites"] == 2
    assert result["jaccard"] == pytest.approx(0.4)
    assert result["assembly_pi_common_mask"] == pytest.approx(0.003)
    assert result["read_pi_common_mask"] == pytest.approx(0.004)
    assert result["pi_difference_read_minus_assembly"] == pytest.approx(0.001)
    assert result["concordant_pi_lower_bracket"] == pytest.approx(0.002)
    assert result["union_pi_upper_bracket"] == pytest.approx(0.005)
    assert result["candidate_assembly_false_positive_upper_bound_sites"] == 1
    assert result["candidate_assembly_false_negative_upper_bound_sites"] == 2
    assert result["method_covariance"] == "paired_shared_reference_and_callable_mask"


def test_cli_manifest_requires_canonical_root(tmp_path: Path) -> None:
    from analysis.vgp_read_validation import validate_root_config

    config = tmp_path / "root.json"
    config.write_text(json.dumps({"root": "/wrong", "layout": {}}))
    with pytest.raises(ValueError, match="canonical VGP root"):
        validate_root_config(config)


def test_paf_summary_merges_target_breadth_and_does_not_claim_callability() -> None:
    paf = io.StringIO(
        "r1\t100\t0\t100\t+\tchr1\t1000\t0\t100\t98\t100\t60\n"
        "r1\t100\t80\t100\t+\tchr2\t500\t100\t120\t19\t20\t40\n"
        "r2\t100\t0\t100\t+\tchr1\t1000\t50\t150\t95\t100\t30\n"
        "r3\t100\t0\t80\t+\tchr2\t500\t10\t90\t72\t80\t10\n"
    )
    result = summarize_paf(paf, reference_bp=1_500, raw_bases=300, raw_reads=3)
    assert result["non_secondary_alignment_segments"] == 4
    assert result["mapped_queries"] == 3
    assert result["mapped_query_bases_union"] == 280
    assert result["reference_breadth_at_least_one_mapping_bp"] == 250
    assert result["weighted_alignment_identity"] == pytest.approx((98 + 19 + 95 + 72) / 300)
    assert result["mapped_query_fraction_of_raw_reads"] == 1
    assert result["callability_status"] == "not_established_from_single_coverage_breadth"


def test_psmc_parser_uses_only_final_optimization_round(tmp_path: Path) -> None:
    psmc = tmp_path / "two-round.psmc"
    psmc.write_text(
        "RD\t0\n"
        "TR\t0.100000\t0.010000\n"
        "RS\t0\t0.000000\t99.000000\n"
        "RS\t1\t1.000000\t99.000000\n"
        "//\n"
        "RD\t25\n"
        "TR\t0.020000\t0.003000\n"
        "RS\t0\t0.000000\t1.500000\n"
        "RS\t1\t1.000000\t2.500000\n"
        "//\n"
    )

    theta, rows = parse_psmc(psmc)
    assert theta == pytest.approx(0.02)
    assert rows == [
        {"interval": 0, "time_2N0": 0.0, "lambda": 1.5},
        {"interval": 1, "time_2N0": 1.0, "lambda": 2.5},
    ]


def test_p07_worker_passes_fastq_path_to_fq2psmcfa_and_retains_failed_scratch() -> None:
    repository_root = Path(__file__).resolve().parents[2]
    script = (repository_root / "analysis/slurm/vgp_read_validation/P07_validate.sh").read_text()
    assert 'fq2psmcfa -q 20 "$WORK/$mask_id.consensus.fastq"' in script
    assert 'fasta-to-fastq < "$WORK/$mask_id.consensus.fa" |' not in script
    assert 'FAILED_WORKDIR=$WORK' in script
    assert '"$WORK/$mask_id.vcf.gz" >> "$mask_dir/read.snps.tsv"' in script
    assert '"$WORK/read.norm.bcf" >> "$mask_dir/read.snps.tsv"' not in script
    assert 'row.get("status") not in {"verified", "reused"}' in script
