import json
import shutil
from pathlib import Path

import pytest

from analysis.tier3_common import Tier3ValidationError, sha256_file
from analysis.tier3a_vgp_collect import (
    audit_phase_orientation,
    collect_deposited_variants,
    collect_direct_alignment,
    impg_core_windows,
    own_normalized_records,
)
from analysis.tier3a_vgp_compute import common_callable_concordance, compute_tier3a


FIXTURES = Path(__file__).parent / "fixtures"


def _store_root(program):
    executable = Path(shutil.which(program) or "").resolve()
    if not executable.is_file() or executable.parent.name != "bin":
        pytest.skip(f"{program} is not available from the pinned Guix profile")
    root = executable.parent.parent
    if not str(root).startswith("/gnu/store/"):
        pytest.skip(f"{program} did not resolve to a Guix store path")
    return str(root)


def _write_fasta(path, name, sequence):
    path.write_text(f">{name}\n{sequence}\n", encoding="utf-8")


def _annotation_provenance(fasta, gff, **updates):
    value = {
        "provider": "truth-provider",
        "release": "1",
        "assembly_accession_version": "GCA_000001.1",
        "status": "native",
        "genetic_code": 1,
        "fasta_sha256": sha256_file(fasta),
        "gff_sha256": sha256_file(gff),
        "contig_mapping": {},
    }
    value.update(updates)
    return value


def test_deposited_truth_uses_explicit_invariant_denominator_and_reference_conditioned_4d():
    result = compute_tier3a(
        dataset_id="truth.deposited",
        reference_fasta=FIXTURES / "truth.fa",
        normalized_bcf=FIXTURES / "truth.normalized.bcf",
        callable_bed=FIXTURES / "expected.callable.bed",
        sample="truth",
        modality="deposited_exact_reference_variants_plus_mask",
        reference_accession="GCA_000001.1",
        annotation_gff=FIXTURES / "truth.gff3",
        annotation_provenance=_annotation_provenance(
            FIXTURES / "truth.fa", FIXTURES / "truth.gff3"
        ),
        minimum_4d_class_sites=1,
    )
    assert result["statistic_label"] == "individual_snv_heterozygosity"
    assert result["population_pi"] is None
    assert result["individual_snv_heterozygosity"] == pytest.approx(2 / 12)
    assert result["heterozygous_snvs"] == 2
    assert result["callable_bases"] == 12
    assert result["fourfold"]["S"]["heterozygous_snvs"] == 0
    assert result["fourfold"]["S"]["callable_bases"] == 1
    assert result["fourfold"]["W"]["heterozygous_snvs"] == 2
    assert result["fourfold"]["W"]["callable_bases"] == 3
    assert result["pi_S_over_pi_W"] == pytest.approx(0.0)
    assert result["annotation"]["status"] == "eligible_native_exact_reference"


def test_variant_only_input_is_missing_never_zero_or_reference_length():
    result = compute_tier3a(
        dataset_id="truth.no-mask",
        reference_fasta=FIXTURES / "truth.fa",
        normalized_bcf=FIXTURES / "truth.normalized.bcf",
        callable_bed=None,
        sample="truth",
        modality="deposited_exact_reference_variants_plus_mask",
        reference_accession="GCA_000001.1",
    )
    assert result["status"] == "unavailable_missing_callable_denominator"
    assert result["individual_snv_heterozygosity"] is None
    assert result["callable_bases"] is None
    assert result["heterozygous_snvs"] is None


def test_deposited_collector_preserves_normalized_exact_reference_tuple(tmp_path):
    result = collect_deposited_variants(
        reference_fasta=FIXTURES / "truth.fa",
        deposited_variants=FIXTURES / "truth.vcf",
        callable_bed=FIXTURES / "expected.callable.bed",
        output_dir=tmp_path / "deposited",
        sample="truth",
        bcftools_store_path=_store_root("bcftools"),
        reference_assembly_accession="GCA_000001.1",
        variant_reference_accession="GCA_000001.1",
    )
    assert result["status"] == "eligible"
    assert result["statistic_label"] == "individual_snv_heterozygosity"
    assert result["population_pi"] is None
    assert result["variant_only_reference_length_assumption_used"] is False
    assert Path(result["outputs"]["original_variants"]).read_bytes() == (
        FIXTURES / "truth.vcf"
    ).read_bytes()
    assert Path(result["outputs"]["normalized_bcf"]).is_file()
    assert Path(result["outputs"]["normalized_bcf_index"]).is_file()

    unavailable = collect_deposited_variants(
        reference_fasta=FIXTURES / "truth.fa",
        deposited_variants=FIXTURES / "truth.vcf",
        callable_bed=None,
        output_dir=tmp_path / "variant-only",
        sample="truth",
        bcftools_store_path=_store_root("bcftools"),
        reference_assembly_accession="GCA_000001.1",
        variant_reference_accession="GCA_000001.1",
    )
    assert unavailable["status"] == "unavailable_missing_callable_denominator"
    assert unavailable["callable_bases"] is None


def test_exact_reference_mismatch_is_rejected(tmp_path):
    wrong = tmp_path / "wrong.fa"
    sequence = (FIXTURES / "truth.fa").read_text().replace("GCT", "GCA", 1)
    wrong.write_text(sequence, encoding="utf-8")
    with pytest.raises(Tier3ValidationError, match="REF mismatch"):
        compute_tier3a(
            dataset_id="truth.mismatch",
            reference_fasta=wrong,
            normalized_bcf=FIXTURES / "truth.normalized.bcf",
            callable_bed=FIXTURES / "expected.callable.bed",
            sample="truth",
            modality="deposited_exact_reference_variants_plus_mask",
            reference_accession="GCA_000001.1",
        )


def test_annotation_degrades_without_native_exact_provenance_but_gc_remains():
    provenance = _annotation_provenance(FIXTURES / "truth.fa", FIXTURES / "truth.gff3")
    provenance["status"] = "projected"
    result = compute_tier3a(
        dataset_id="truth.projected",
        reference_fasta=FIXTURES / "truth.fa",
        normalized_bcf=FIXTURES / "truth.normalized.bcf",
        callable_bed=FIXTURES / "expected.callable.bed",
        sample="truth",
        modality="deposited_exact_reference_variants_plus_mask",
        reference_accession="GCA_000001.1",
        annotation_gff=FIXTURES / "truth.gff3",
        annotation_provenance=provenance,
    )
    assert result["whole_genome_gc"]["value"] is not None
    assert result["gc3"] is None
    assert result["fourfold"] is None
    assert result["annotation"]["status"] == "unavailable_non_native_annotation"


def test_direct_alignment_retains_provenance_qc_bed_bcf_csi_and_store_paths(tmp_path):
    h1 = tmp_path / "h1.fa"
    h2 = tmp_path / "h2.fa"
    target = list("A" * 230)
    query = list(target)
    target[110] = "C"
    query[110] = "G"
    query[116] = "T"
    target[119] = "C"
    _write_fasta(h1, "h1", "".join(target))
    _write_fasta(h2, "h2", "".join(query))
    paf = tmp_path / "raw.input.paf"
    paf.write_text(
        "h2\t230\t0\t230\t+\th1\t230\t0\t230\t228\t230\t60\t"
        "cg:Z:110=1X5=1I3=1D9=101=\n",
        encoding="utf-8",
    )
    out = tmp_path / "out"
    record = collect_direct_alignment(
        h1_fasta=h1,
        h2_fasta=h2,
        paf_path=paf,
        output_dir=out,
        sample="individual",
        wfmash_store_path=_store_root("wfmash"),
        bcftools_store_path=_store_root("bcftools"),
        phase_qc_passed=True,
        collapse_qc_passed=True,
        edge_exclusion_bp=100,
        indel_flank_bp=1,
    )
    assert (out / "raw.wfmash.paf").read_bytes() == paf.read_bytes()
    assert Path(record["outputs"]["accepted_mapping_qc"]).is_file()
    assert Path(record["outputs"]["callable_bed"]).read_text().strip()
    assert Path(record["outputs"]["normalized_snv_bcf"]).is_file()
    assert Path(record["outputs"]["normalized_snv_bcf_index"]).is_file()
    assert record["tools"]["wfmash_store_path"].startswith("/gnu/store/")
    assert record["tools"]["bcftools_store_path"].startswith("/gnu/store/")
    assert record["operation_counts"] == {"=": 228, "X": 1, "I": 1, "D": 1}
    assert record["heterozygous_snvs"] == 1


@pytest.mark.parametrize(
    "cigar,error",
    [(None, "cg:Z"), ("230M", "extended")],
)
def test_direct_alignment_rejects_missing_or_approximate_cigar(tmp_path, cigar, error):
    h1, h2, paf = tmp_path / "h1.fa", tmp_path / "h2.fa", tmp_path / "x.paf"
    _write_fasta(h1, "h1", "A" * 230)
    _write_fasta(h2, "h2", "A" * 230)
    tag = "" if cigar is None else f"\tcg:Z:{cigar}"
    paf.write_text(f"h2\t230\t0\t230\t+\th1\t230\t0\t230\t230\t230\t60{tag}\n")
    with pytest.raises(Tier3ValidationError, match=error):
        collect_direct_alignment(
            h1, h2, paf, tmp_path / "out", "i",
            _store_root("wfmash"), _store_root("bcftools"), True, True,
        )


def test_direct_alignment_reverse_cigar_swap_edges_ambiguity_overlap_and_gaps(tmp_path):
    h1, h2 = tmp_path / "h1.fa", tmp_path / "h2.fa"
    _write_fasta(h1, "target", "ACGT")
    _write_fasta(h2, "query", "ACGT")
    reverse = tmp_path / "reverse.paf"
    reverse.write_text("query\t4\t0\t4\t-\ttarget\t4\t0\t4\t4\t4\t60\tcg:Z:4=\n")
    record = collect_direct_alignment(
        h1, h2, reverse, tmp_path / "reverse-out", "i",
        _store_root("wfmash"), _store_root("bcftools"), True, True, 0, 0,
    )
    assert record["callable_bases"] == 4

    # Swapping H1/H2 changes REF/ALT orientation, but not individual
    # heterozygosity.  This guards against silently treating H2 as reference.
    _write_fasta(h1, "target", "ACGT")
    _write_fasta(h2, "query", "AGGT")
    forward = tmp_path / "forward.paf"
    forward.write_text("query\t4\t0\t4\t+\ttarget\t4\t0\t4\t3\t4\t60\tcg:Z:1=1X2=\n")
    a = collect_direct_alignment(
        h1, h2, forward, tmp_path / "a", "i",
        _store_root("wfmash"), _store_root("bcftools"), True, True, 0, 0,
    )
    swapped = tmp_path / "swapped.paf"
    swapped.write_text("target\t4\t0\t4\t+\tquery\t4\t0\t4\t3\t4\t60\tcg:Z:1=1X2=\n")
    b = collect_direct_alignment(
        h2, h1, swapped, tmp_path / "b", "i",
        _store_root("wfmash"), _store_root("bcftools"), True, True, 0, 0,
    )
    assert a["heterozygous_snvs"] == b["heterozygous_snvs"] == 1
    assert a["variant_alleles"][0]["ref"] == "C"
    assert b["variant_alleles"][0]["ref"] == "G"

    # A duplicated mapping and ambiguous target base eliminate those bases;
    # an insertion gap and flank do not enter the denominator either.
    _write_fasta(h1, "target", "AANAAA")
    _write_fasta(h2, "query", "AATAAAA")
    complex_paf = tmp_path / "complex.paf"
    line = "query\t7\t0\t7\t+\ttarget\t6\t0\t6\t5\t7\t60\tcg:Z:2=1X1I3=\n"
    complex_paf.write_text(line + line)
    c = collect_direct_alignment(
        h1, h2, complex_paf, tmp_path / "c", "i",
        _store_root("wfmash"), _store_root("bcftools"), True, True, 0, 0,
    )
    assert c["callable_bases"] == 0
    assert c["exclusion_counts"]["multiple_projection"] > 0


def test_direct_alignment_rejects_collapse_phase_and_annotation_mismatch(tmp_path):
    h1, h2, paf = tmp_path / "h1.fa", tmp_path / "h2.fa", tmp_path / "x.paf"
    _write_fasta(h1, "h1", "AAAA")
    _write_fasta(h2, "h2", "AAAA")
    paf.write_text("h2\t4\t0\t4\t+\th1\t4\t0\t4\t4\t4\t60\tcg:Z:4=\n")
    common = dict(
        h1_fasta=h1, h2_fasta=h2, paf_path=paf, output_dir=tmp_path / "out",
        sample="i", wfmash_store_path=_store_root("wfmash"),
        bcftools_store_path=_store_root("bcftools"), edge_exclusion_bp=0,
    )
    with pytest.raises(Tier3ValidationError, match="phase"):
        collect_direct_alignment(**common, phase_qc_passed=False, collapse_qc_passed=True)
    with pytest.raises(Tier3ValidationError, match="collapse"):
        collect_direct_alignment(**common, phase_qc_passed=True, collapse_qc_passed=False)
    with pytest.raises(Tier3ValidationError, match="annotation.*H1"):
        collect_direct_alignment(
            **common, phase_qc_passed=True, collapse_qc_passed=True,
            h1_assembly_accession="GCA_1.1", annotation_assembly_accession="GCA_2.1",
        )


def test_deposited_and_haplotype_modalities_agree_on_common_callable_truth(tmp_path):
    h1, h2, paf = tmp_path / "h1.fa", tmp_path / "h2.fa", tmp_path / "truth.paf"
    reference = "GCTAGCCAAAAAAAAAAAATCCACCTAAAAAAAAA"
    alternate = list(reference)
    alternate[2] = "C"
    alternate[19] = "C"
    _write_fasta(h1, "chr1", reference)
    _write_fasta(h2, "chr1", "".join(alternate))
    paf.write_text(
        "chr1\t35\t0\t35\t+\tchr1\t35\t0\t35\t33\t35\t60\t"
        "cg:Z:2=1X16=1X15=\n"
    )
    alignment = collect_direct_alignment(
        h1, h2, paf, tmp_path / "alignment", "truth",
        _store_root("wfmash"), _store_root("bcftools"), True, True, 0, 0,
        h1_accessibility_bed=FIXTURES / "expected.callable.bed",
    )
    alignment_hets = {
        (item["contig"], item["position_1based"] - 1)
        for item in alignment["variant_alleles"]
    }
    deposited_hets = {("chr1", 2), ("chr1", 19)}
    audit = common_callable_concordance(
        reference_fasta=h1,
        left_callable_bed=FIXTURES / "expected.callable.bed",
        right_callable_bed=alignment["outputs"]["callable_bed"],
        left_heterozygous=deposited_hets,
        right_heterozygous=alignment_hets,
        synthetic_fixture=True,
    )
    assert audit["common_callable_bases"] == 12
    assert audit["snv_precision"] == audit["snv_recall"] == 1.0
    assert audit["heterozygous_nonheterozygous_genotype_concordance"] == 1.0
    assert audit["left_individual_snv_heterozygosity"] == pytest.approx(2 / 12)
    assert audit["right_individual_snv_heterozygosity"] == pytest.approx(2 / 12)
    assert audit["passed"] is True


def test_impg_padded_cores_own_boundaries_once_and_phase_is_audited():
    windows = impg_core_windows({"chr1": 2_000_001}, core_size_bp=1_000_000, padding_bp=10_000)
    assert [(w.core_start, w.core_end, w.query_start, w.query_end) for w in windows] == [
        (0, 1_000_000, 0, 1_010_000),
        (1_000_000, 2_000_000, 990_000, 2_000_001),
        (2_000_000, 2_000_001, 1_990_000, 2_000_001),
    ]
    records = [
        {"contig": "chr1", "pos": 1_000_000, "ref": "A", "alt": "C"},
        {"contig": "chr1", "pos": 1_000_001, "ref": "A", "alt": "G"},
        {"contig": "chr1", "pos": 1_000_001, "ref": "A", "alt": "G"},
    ]
    owned = own_normalized_records(records, windows)
    assert [(item["pos"], item["owner_core_start"]) for item in owned] == [
        (1_000_000, 0),
        (1_000_001, 1_000_000),
    ]
    audit = audit_phase_orientation(
        expected={("chr1", 1_000_001, "A", "G"): (0, 1)},
        observed={("chr1", 1_000_001, "A", "G"): (1, 0)},
    )
    assert audit["result"] == "orientation_inverted"
    assert audit["phase_sensitive_eligible"] is False
