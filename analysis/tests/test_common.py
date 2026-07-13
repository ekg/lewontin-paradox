import hashlib
import json
from pathlib import Path

import pytest

from analysis.tier3_common import (
    Tier3ValidationError,
    allele_pairwise_diversity,
    collect_fourfold_sites,
    intersect_intervals,
    parse_gff,
    parse_gt,
    read_fasta,
    reconstruct_cds,
    resolve_contig_aliases,
    sha256_file,
    traverse_paf,
    verify_file,
)


FIXTURES = Path(__file__).parent / "fixtures"


def test_gff_strand_phase_and_transcript_selection():
    fasta = read_fasta(FIXTURES / "truth.fa")
    annotation = parse_gff(FIXTURES / "truth.gff3")

    # tx_plus has two biological segments; the second phase removes the
    # deliberately prefixed base.  tx_minus tests descending segment order,
    # reverse-complementing, and phase at the biological start.
    assert reconstruct_cds(fasta, annotation.transcripts["tx_plus"]) == "GCTGCC"
    assert reconstruct_cds(fasta, annotation.transcripts["tx_minus"]) == "GGTGGA"
    sites, excluded = collect_fourfold_sites(fasta, annotation)
    assert ("chr1", 2) in sites
    assert ("chr1", 22) in sites
    assert not excluded
    # The provider-designated canonical transcript, not merely the longest,
    # is retained for gene_plus.
    assert annotation.canonical_transcripts["gene_plus"] == "tx_plus"


def test_gff_overlapping_discordant_frames_are_excluded(tmp_path):
    fasta_path = tmp_path / "x.fa"
    fasta_path.write_text(">c\nGCTGCTGCT\n")
    gff_path = tmp_path / "x.gff3"
    gff_path.write_text(
        "##gff-version 3\n##sequence-region c 1 9\n"
        "c\tt\tmRNA\t1\t9\t.\t+\t.\tID=t1;Parent=g1;tag=canonical\n"
        "c\tt\tCDS\t1\t9\t.\t+\t0\tParent=t1\n"
        "c\tt\tmRNA\t2\t7\t.\t+\t.\tID=t2;Parent=g2;tag=canonical\n"
        "c\tt\tCDS\t2\t7\t.\t+\t0\tParent=t2\n"
    )
    sites, excluded = collect_fourfold_sites(read_fasta(fasta_path), parse_gff(gff_path))
    assert ("c", 2) not in sites
    assert ("c", 2) in excluded


def test_contig_aliases_are_explicit_one_to_one_and_length_checked():
    assert resolve_contig_aliases(
        {"chr1": 30, "chr2": 20}, {"1": 30, "2": 20}, {"1": "chr1", "2": "chr2"}
    ) == {"1": "chr1", "2": "chr2"}
    with pytest.raises(Tier3ValidationError, match="undeclared"):
        resolve_contig_aliases({"chr1": 30}, {"1": 30}, {})
    with pytest.raises(Tier3ValidationError, match="one-to-one"):
        resolve_contig_aliases({"chr1": 30}, {"1": 30, "01": 30}, {"1": "chr1", "01": "chr1"})
    with pytest.raises(Tier3ValidationError, match="length"):
        resolve_contig_aliases({"chr1": 30}, {"1": 29}, {"1": "chr1"})


def test_genotype_missingness_ploidy_and_unbiased_diversity():
    assert parse_gt("0|1", expected_ploidy=2) == (0, 1)
    assert parse_gt("./1", expected_ploidy=2) == (None, 1)
    with pytest.raises(Tier3ValidationError, match="ploidy"):
        parse_gt("0/1/1", expected_ploidy=2)
    assert parse_gt("0/1", expected_ploidy=1, haploidize_heterozygous=True) == (None,)
    # Three called copies (0, 0, 1): 1 - 2/(3*2) = 2/3.
    assert allele_pairwise_diversity([(0, 0), (1, None)], minimum_called=3) == pytest.approx(2 / 3)
    assert allele_pairwise_diversity([(0, None), (None, 1)], minimum_called=3) is None


def test_paf_extended_cigar_eq_x_i_d_traversal_ambiguity_and_overlap():
    target_bases = list("A" * 240)
    query_bases = list("A" * 240)
    target_bases[111], query_bases[111] = "C", "G"  # X
    query_bases[117] = "T"  # I relative to H1
    target_bases[120] = "C"  # D relative to H1
    target = {"t": "".join(target_bases)}
    query = {"q": "".join(query_bases)}
    paf = [
        "q\t240\t0\t240\t+\tt\t240\t0\t240\t238\t240\t60\tcg:Z:111=1X5=1I3=1D9=110=",
    ]
    result = traverse_paf(paf, target, query, edge_exclusion_bp=100, indel_flank_bp=1)
    assert ("t", 111) in result.snv_positions
    assert ("t", 110) in result.callable_positions
    # The I/D anchors and one-base flanks are excluded from the denominator.
    assert ("t", 118) not in result.callable_positions
    assert result.operation_counts == {"=": 238, "X": 1, "I": 1, "D": 1}

    ambiguous_query = {"q": query["q"][:112] + "N" + query["q"][113:]}
    ambiguous = traverse_paf(paf, target, ambiguous_query, edge_exclusion_bp=100, indel_flank_bp=1)
    assert ("t", 112) not in ambiguous.callable_positions
    assert ambiguous.exclusion_counts["ambiguous_base"] == 1

    overlap = traverse_paf(paf + paf, target, query, edge_exclusion_bp=100, indel_flank_bp=1)
    assert not overlap.callable_positions
    assert overlap.exclusion_counts["multiple_projection"] > 0


def test_paf_rejects_m_and_supports_reverse_strand():
    with pytest.raises(Tier3ValidationError, match="extended"):
        traverse_paf(["q\t4\t0\t4\t+\tt\t4\t0\t4\t4\t4\t60\tcg:Z:4M"], {"t": "ACGT"}, {"q": "ACGT"}, 0, 0)
    result = traverse_paf(
        ["q\t4\t0\t4\t-\tt\t4\t0\t4\t4\t4\t60\tcg:Z:4="],
        {"t": "ACGT"},
        {"q": "ACGT"},
        edge_exclusion_bp=0,
        indel_flank_bp=0,
    )
    assert len(result.callable_positions) == 4


def test_interval_intersection_is_half_open_and_merges():
    assert intersect_intervals([(0, 5), (4, 8), (10, 12)], [(3, 6), (7, 11)]) == [
        (3, 6),
        (7, 8),
        (10, 11),
    ]


def test_checksums_fail_closed(tmp_path):
    p = tmp_path / "object"
    p.write_bytes(b"tier3\n")
    digest = hashlib.sha256(b"tier3\n").hexdigest()
    assert sha256_file(p) == digest
    verify_file(p, digest, 6)
    with pytest.raises(Tier3ValidationError, match="SHA-256"):
        verify_file(p, "0" * 64, 6)
    with pytest.raises(Tier3ValidationError, match="size"):
        verify_file(p, digest, 7)


def test_fixture_expected_truth_is_independent_of_variant_rows():
    import pysam

    expected = json.loads((FIXTURES / "expected.json").read_text())
    assert sha256_file(FIXTURES / "truth.normalized.bcf") == expected["normalized_bcf_sha256"]
    assert sha256_file(FIXTURES / "truth.normalized.bcf.csi") == expected["normalized_bcf_index_sha256"]
    with pysam.VariantFile(str(FIXTURES / "truth.normalized.bcf")) as variants:
        records = list(variants)
    snvs = [record for record in records if len(record.ref) == 1 and all(len(alt) == 1 for alt in record.alts)]
    indels = [record for record in records if record not in snvs]
    observed_truth = [
        {
            "contig": record.contig,
            "position_1based": record.pos,
            "ref": record.ref,
            "alt": record.alts[0],
            "genotypes": ["|".join(str(allele) for allele in record.samples[sample]["GT"]) for sample in record.samples],
        }
        for record in snvs
    ]
    assert observed_truth == expected["allele_truth"]
    assert len(snvs) == expected["variant_snv_count"]
    assert len(indels) == expected["variant_indel_count"]

    def read_bed(name):
        return [line.rstrip().split("\t") for line in (FIXTURES / name).read_text().splitlines()]

    callable_intervals = read_bed(expected["callable_bed"])
    fourfold_intervals = read_bed(expected["fourfold_intersection_bed"])
    callable_bases = sum(int(end) - int(start) for _contig, start, end in callable_intervals)
    fourfold_bases = sum(int(end) - int(start) for _contig, start, end in fourfold_intervals)
    assert callable_bases == expected["callable_bases"]
    assert fourfold_bases == expected["fourfold_callable_bases"]
    assert expected["callable_bases"] > expected["variant_snv_count"]
    assert expected["fourfold_callable_bases"] == sum(expected["fourfold_by_class"].values())
