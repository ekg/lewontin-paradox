from collections import Counter

from analysis.audit_vgp_three_pair_results import (
    fasta_lengths,
    fasta_non_n_counts,
    pipeline_limitations,
    psmc_population_counts,
    select_strata,
    variant_counts,
)


def test_independent_stratified_reaudit_counts_fasta_vcf_and_psmc_populations(tmp_path):
    fasta = tmp_path / "consensus.fa"
    fasta.write_text(
        ">early\n" + "A" * 5_000_000 + "\n>middle\n" + "N" * 2_500_000 +
        "A" * 5_000_000 + "N" * 2_500_000 + "\n>late\n" + "C" * 5_000_000 + "\n"
    )
    regions = select_strata(fasta_lengths(fasta))
    assert [row.contig for row in regions] == ["early", "middle", "late"]
    non_n = fasta_non_n_counts(fasta, regions)
    assert [non_n[row] for row in regions] == [5_000_000, 5_000_000, 5_000_000]

    vcf = tmp_path / "calls.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "early\t1\t.\tA\tG\t.\tPASS\t.\n"
        "middle\t5000001\t.\tA\tAT\t.\tPASS\t.\n"
        "late\t5000000\t.\tC\tT\t.\tPASS\t.\n"
    )
    counts = variant_counts(vcf, regions)
    assert [counts[row]["records"] for row in regions] == [1, 1, 1]
    assert [counts[row]["snps"] for row in regions] == [1, 0, 1]

    psmcfa = tmp_path / "input.psmcfa"
    psmcfa.write_text(
        ">early\n" + "K" * 50_000 + "\n>middle\n" + "N" * 100_000 +
        "\n>late\n" + "T" * 50_000 + "\n"
    )
    populations = psmc_population_counts(psmcfa, regions)
    assert populations[regions[0]] == Counter(K=50_000)
    assert populations[regions[1]] == Counter(N=50_000)
    assert populations[regions[2]] == Counter(T=50_000)


def test_pipeline_limitations_retain_pinned_impg_single_thread_stall():
    limitations = pipeline_limitations()
    assert any("lace -t 1" in limitation and "non-progressing" in limitation
               for limitation in limitations)
