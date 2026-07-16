import json
import math
import shutil
import subprocess
from pathlib import Path

import pytest

from analysis.tier3_common import Tier3ValidationError
from analysis.tier3b_popvcf_collect import (
    collect_population_vcf,
    deterministic_sample_selection,
)
from analysis.tier3b_popvcf_compute import compute_population_pi


FIXTURES = Path(__file__).with_name("fixtures")


def _samples(count=20):
    return ["s{:02d}".format(index) for index in range(count)]


def _write_fasta(tmp_path, sequences):
    path = tmp_path / "reference.fa"
    path.write_text(
        "".join(">{}\n{}\n".format(contig, sequence) for contig, sequence in sequences.items()),
        encoding="utf-8",
    )
    return path


def _write_vcf(tmp_path, sequences, samples, rows, name="calls.vcf"):
    path = tmp_path / name
    header = [
        "##fileformat=VCFv4.2",
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    ]
    header.extend(
        "##contig=<ID={},length={}>".format(contig, len(sequence))
        for contig, sequence in sequences.items()
    )
    header.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples))
    lines = []
    for row in rows:
        contig, pos, ref, alt, genotypes = row[:5]
        filter_value = row[5] if len(row) > 5 else "PASS"
        info = row[6] if len(row) > 6 else "."
        lines.append(
            "{}\t{}\t.\t{}\t{}\t60\t{}\t{}\tGT\t{}".format(
                contig, pos, ref, alt, filter_value, info, "\t".join(genotypes)
            )
        )
    path.write_text("\n".join(header + lines) + "\n", encoding="utf-8")
    return path


def _write_bed(tmp_path, intervals, name="callable.bed"):
    path = tmp_path / name
    path.write_text(
        "".join("{}\t{}\t{}\n".format(contig, start, end) for contig, start, end in intervals),
        encoding="utf-8",
    )
    return path


def _compute(tmp_path, vcf, fasta, samples, bed=None, kind="cohort_callable_mask", design="wild_diploid", **kwargs):
    return compute_population_pi(
        dataset_id=kwargs.pop("dataset_id", "synthetic"),
        vcf_path=vcf,
        fasta_path=fasta,
        selected_samples=samples,
        design=design,
        denominator_kind=kind,
        callable_bed_path=bed,
        bootstrap_replicates=kwargs.pop("bootstrap_replicates", 32),
        **kwargs
    )


def test_variant_only_input_never_manufactures_invariant_denominator(tmp_path):
    samples = _samples()
    sequences = {"chr1": "A" * 10}
    fasta = _write_fasta(tmp_path, sequences)
    genotypes = ["0/1"] + ["0/0"] * 19
    vcf = _write_vcf(tmp_path, sequences, samples, [("chr1", 2, "A", "G", genotypes)])
    with pytest.raises(Tier3ValidationError, match="variant-only"):
        _compute(tmp_path, vcf, fasta, samples, kind="all_sites_vcf")
    with pytest.raises(Tier3ValidationError, match="explicit denominator"):
        _compute(tmp_path, vcf, fasta, samples, kind="unavailable")


def test_mask_and_all_sites_denominators_include_invariant_sites_and_match_hand_calculation(tmp_path):
    samples = _samples()
    sequences = {"chr1": "A" * 10}
    fasta = _write_fasta(tmp_path, sequences)
    variant_genotypes = ["0/1"] + ["0/0"] * 19
    variant = ("chr1", 2, "A", "G", variant_genotypes)
    sparse = _write_vcf(tmp_path, sequences, samples, [variant], "sparse.vcf")
    bed = _write_bed(tmp_path, [("chr1", 0, 10)])
    mask_result = _compute(tmp_path, sparse, fasta, samples, bed)

    rows = [variant]
    invariant = ["0/0"] * 20
    rows.extend(("chr1", pos, "A", ".", invariant) for pos in range(1, 11) if pos != 2)
    rows.sort(key=lambda row: row[1])
    all_sites = _write_vcf(tmp_path, sequences, samples, rows, "all-sites.vcf")
    all_sites_result = _compute(tmp_path, all_sites, fasta, samples, kind="all_sites_vcf")

    # One alternate among 40 chromosomes: site pi=2*1*39/(40*39)=0.05;
    # nine explicit invariant sites make the genome-wide estimate 0.005.
    for result in (mask_result, all_sites_result):
        assert result["population_pi"]["diversity_sum"] == pytest.approx(0.05)
        assert result["population_pi"]["callable_sites"] == 10
        assert result["population_pi"]["point_estimate"] == pytest.approx(0.005)


def test_invariant_only_all_sites_vcf_has_zero_pi(tmp_path):
    samples = _samples()
    sequences = {"chr1": "ACGT"}
    fasta = _write_fasta(tmp_path, sequences)
    invariant = ["0/0"] * 20
    rows = [("chr1", pos + 1, base, ".", invariant) for pos, base in enumerate(sequences["chr1"])]
    vcf = _write_vcf(tmp_path, sequences, samples, rows)
    result = _compute(tmp_path, vcf, fasta, samples, kind="all_sites_vcf")
    assert result["population_pi"]["callable_sites"] == 4
    assert result["population_pi"]["diversity_sum"] == 0
    assert result["population_pi"]["point_estimate"] == 0


def test_gvcf_reference_blocks_are_explicit_and_failed_blocks_are_wholly_removed(tmp_path):
    samples = _samples()
    sequences = {"chr1": "A" * 10}
    fasta = _write_fasta(tmp_path, sequences)
    called = ["0/0"] * 20
    under_called = ["./."] * 20
    rows = [
        ("chr1", 1, "A", "<NON_REF>", called, "PASS", "END=5"),
        ("chr1", 6, "A", "<NON_REF>", under_called, "PASS", "END=10"),
    ]
    vcf = _write_vcf(tmp_path, sequences, samples, rows)
    result = _compute(tmp_path, vcf, fasta, samples, kind="gvcf")
    assert result["population_pi"]["callable_sites"] == 5
    assert result["population_pi"]["point_estimate"] == 0
    assert result["exclusion_counts"]["insufficient_called_chromosomes"] == 1


def test_per_site_missing_chromosome_count_uses_actual_n_and_threshold(tmp_path):
    samples = _samples()
    sequences = {"chr1": "AAA"}
    fasta = _write_fasta(tmp_path, sequences)
    bed = _write_bed(tmp_path, [("chr1", 0, 3)])
    # Exactly 36 called chromosomes is retained and the formula uses n=36.
    retained = ["0/1"] + ["0/0"] * 17 + ["./."] * 2
    # 35 called chromosomes is excluded from both numerator and denominator.
    excluded = ["0/1"] + ["0/0"] * 16 + ["0/."] + ["./."] * 2
    vcf = _write_vcf(
        tmp_path,
        sequences,
        samples,
        [("chr1", 1, "A", "G", retained), ("chr1", 2, "A", "G", excluded)],
    )
    result = _compute(tmp_path, vcf, fasta, samples, bed)
    expected = 1.0 - (35 * 34) / (36 * 35)
    assert result["population_pi"]["diversity_sum"] == pytest.approx(expected)
    assert result["population_pi"]["callable_sites"] == 2
    assert result["population_pi"]["called_chromosome_histogram_at_records"] == {"35": 1, "36": 1}
    assert result["exclusion_counts"]["insufficient_called_chromosomes"] == 1


def test_multiallelic_snv_uses_all_allele_counts_once(tmp_path):
    samples = _samples()
    sequences = {"chr1": "A"}
    fasta = _write_fasta(tmp_path, sequences)
    bed = _write_bed(tmp_path, [("chr1", 0, 1)])
    genotypes = ["0/1", "2/2"] + ["0/0"] * 18
    vcf = _write_vcf(tmp_path, sequences, samples, [("chr1", 1, "A", "C,G", genotypes)])
    result = _compute(tmp_path, vcf, fasta, samples, bed)
    assert result["population_pi"]["diversity_sum"] == pytest.approx(113 / 780)
    assert result["population_pi"]["callable_sites"] == 1


def test_haploid_diploid_and_inbred_genotype_semantics(tmp_path):
    sequences = {"chr1": "A"}
    fasta = _write_fasta(tmp_path, sequences)
    bed = _write_bed(tmp_path, [("chr1", 0, 1)])

    diploid_samples = _samples()
    diploid_vcf = _write_vcf(
        tmp_path, sequences, diploid_samples,
        [("chr1", 1, "A", "G", ["0/1"] + ["0/0"] * 19)], "diploid.vcf",
    )
    assert _compute(tmp_path, diploid_vcf, fasta, diploid_samples, bed)["population_pi"]["diversity_sum"] == pytest.approx(0.05)

    haploid_samples = _samples()
    haploid_vcf = _write_vcf(
        tmp_path, sequences, haploid_samples,
        [("chr1", 1, "A", "G", ["1"] + ["0"] * 19)], "haploid.vcf",
    )
    haploid = _compute(tmp_path, haploid_vcf, fasta, haploid_samples, bed, design="haploid")
    assert haploid["population_pi"]["diversity_sum"] == pytest.approx(0.1)

    # Inbred diploid-format homozygotes yield one consensus allele.  Two
    # heterozygous/uncertain lines are missing; the remaining 18 meet the gate.
    inbred_vcf = _write_vcf(
        tmp_path, sequences, diploid_samples,
        [("chr1", 1, "A", "G", ["1/1"] + ["0/0"] * 17 + ["0/1", "./."])], "inbred.vcf",
    )
    inbred = _compute(
        tmp_path, inbred_vcf, fasta, diploid_samples, bed, design="inbred_lines_haploidized"
    )
    expected = 1.0 - (17 * 16) / (18 * 17)
    assert inbred["population_pi"]["diversity_sum"] == pytest.approx(expected)
    assert inbred["population_pi"]["called_chromosome_histogram_at_records"] == {"18": 1}


def test_filtered_indels_and_ambiguous_reference_do_not_enter_denominator(tmp_path):
    samples = _samples()
    sequences = {"chr1": "AANA"}
    fasta = _write_fasta(tmp_path, sequences)
    bed = _write_bed(tmp_path, [("chr1", 0, 4)])
    calls = ["0/1"] + ["0/0"] * 19
    rows = [
        ("chr1", 1, "A", "G", calls, "LowQual"),
        ("chr1", 2, "A", "AT", calls),
        ("chr1", 3, "N", ".", ["0/0"] * 20),
    ]
    vcf = _write_vcf(tmp_path, sequences, samples, rows)
    result = _compute(tmp_path, vcf, fasta, samples, bed)
    assert result["denominator"]["initial_callable_acgt_sites"] == 3
    assert result["population_pi"]["callable_sites"] == 1
    assert result["exclusion_counts"] == {"ambiguous_reference": 1, "filtered": 1, "non_snv": 1}


def test_exact_reference_dictionary_and_ref_mismatch_are_fatal(tmp_path):
    samples = _samples()
    sequences = {"chr1": "AAAA"}
    fasta = _write_fasta(tmp_path, sequences)
    bed = _write_bed(tmp_path, [("chr1", 0, 4)])
    genotypes = ["0/0"] * 20
    wrong_ref = _write_vcf(tmp_path, sequences, samples, [("chr1", 1, "C", ".", genotypes)])
    with pytest.raises(Tier3ValidationError, match="REF mismatch"):
        _compute(tmp_path, wrong_ref, fasta, samples, bed)

    wrong_dictionary = _write_vcf(tmp_path, {"chr1": "AAAAA"}, samples, [], "wrong-dict.vcf")
    with pytest.raises(Tier3ValidationError, match="contig length"):
        _compute(tmp_path, wrong_dictionary, fasta, samples, bed)


def test_population_structure_and_downsampling_are_deterministic():
    rows = []
    for population, count in (("popB", 22), ("popA", 22), ("smaller", 21)):
        for index in range(count):
            rows.append({
                "sample_id": "{}_{}".format(population, index),
                "population_id": population,
                "qc_pass": "true",
                "unrelated": "true",
            })
    rows.append({"sample_id": "bad", "population_id": "popA", "qc_pass": "false"})
    first = deterministic_sample_selection("dataset", rows)
    second = deterministic_sample_selection("dataset", list(reversed(rows)))
    assert first == second
    assert first["population_id"] == "popA"  # equal-largest tie is bytewise
    assert len(first["selected_samples"]) == 20
    assert all(sample.startswith("popA_") for sample in first["selected_samples"])
    assert first["exclusion_list_sha256"] != "0" * 64


def _annotation_metadata(status="native", assembly="GCA_000001.1", fasta_assembly="GCA_000001.1"):
    return {
        "provider": "synthetic-provider",
        "release": "fixture-1",
        "assembly_accession": assembly,
        "fasta_assembly_accession": fasta_assembly,
        "status": status,
        "genetic_code": 1,
    }


def test_reference_conditioned_W_S_classes_and_annotation_provenance(tmp_path):
    samples = _samples()
    fasta = tmp_path / "truth.fa"
    fasta.write_bytes((FIXTURES / "truth.fa").read_bytes())
    bed = tmp_path / "callable.bed"
    bed.write_bytes((FIXTURES / "expected.callable.bed").read_bytes())
    sequences = {"chr1": "GCTAGCCAAAAAAAAAAAATCCACCTAAAAAAAAA"}
    calls = ["0/1"] * 20
    vcf = _write_vcf(
        tmp_path, sequences, samples,
        [("chr1", 3, "T", "C", calls), ("chr1", 7, "C", "T", calls)],
    )
    result = _compute(
        tmp_path,
        vcf,
        fasta,
        samples,
        bed,
        gff_path=FIXTURES / "truth.gff3",
        annotation_metadata=_annotation_metadata(),
        minimum_4d_class_sites=1,
    )
    ratio = result["pi_S_over_pi_W"]
    assert ratio["S"]["callable_sites"] == 1
    assert ratio["W"]["callable_sites"] == 3
    assert ratio["S"]["diversity_sum"] == pytest.approx(20 / 39)
    assert ratio["W"]["diversity_sum"] == pytest.approx(20 / 39)
    assert ratio["point_estimate"] == pytest.approx(3.0)
    assert ratio["bootstrap"]["unavailable_reason"] == "fewer_than_20_eligible_blocks"
    assert result["annotation"]["status"] == "native"
    assert result["annotation"]["provider"] == "synthetic-provider"
    assert len(result["annotation"]["fasta_sha256"]) == 64
    assert len(result["annotation"]["gff_sha256"]) == 64
    assert result["annotation"]["contig_dictionary_passed"] is True


def test_explicit_native_transcript_exclusion_is_audited_without_editing_gff(tmp_path):
    samples = _samples()
    fasta = tmp_path / "truth.fa"
    fasta.write_bytes((FIXTURES / "truth.fa").read_bytes())
    bed = tmp_path / "callable.bed"
    bed.write_bytes((FIXTURES / "expected.callable.bed").read_bytes())
    sequences = {"chr1": "GCTAGCCAAAAAAAAAAAATCCACCTAAAAAAAAA"}
    vcf = _write_vcf(tmp_path, sequences, samples, [])
    metadata = _annotation_metadata()
    metadata["excluded_transcripts"] = {"tx_minus": "provider CDS ambiguity documented upstream"}
    result = _compute(
        tmp_path,
        vcf,
        fasta,
        samples,
        bed,
        gff_path=FIXTURES / "truth.gff3",
        annotation_metadata=metadata,
    )
    assert result["annotation"]["excluded_transcripts"] == metadata["excluded_transcripts"]
    assert result["annotation"]["declared_excluded_transcripts"] == metadata["excluded_transcripts"]
    assert result["annotation"]["retained_transcripts"] == 1
    assert result["annotation"]["gff_sha256"]

    metadata["excluded_transcripts"] = {"absent": "not allowed"}
    with pytest.raises(Tier3ValidationError, match="absent transcripts"):
        _compute(
            tmp_path,
            vcf,
            fasta,
            samples,
            bed,
            gff_path=FIXTURES / "truth.gff3",
            annotation_metadata=metadata,
        )


def test_invalid_canonical_transcripts_can_be_excluded_with_exact_audit(tmp_path):
    samples = _samples()
    sequences = {"chr1": "GCN"}
    fasta = _write_fasta(tmp_path, sequences)
    bed = _write_bed(tmp_path, [("chr1", 0, 3)])
    vcf = _write_vcf(tmp_path, sequences, samples, [])
    gff = tmp_path / "native.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "##sequence-region chr1 1 3\n"
        "chr1\tx\tgene\t1\t3\t.\t+\t.\tID=g\n"
        "chr1\tx\tmRNA\t1\t3\t.\t+\t.\tID=t;Parent=g;tag=canonical\n"
        "chr1\tx\tCDS\t1\t3\t.\t+\t0\tParent=t\n",
        encoding="utf-8",
    )
    metadata = _annotation_metadata()
    metadata["invalid_transcript_policy"] = "exclude_with_audit"
    with pytest.raises(Tier3ValidationError, match="remove every retained transcript"):
        _compute(
            tmp_path,
            vcf,
            fasta,
            samples,
            bed,
            gff_path=gff,
            annotation_metadata=metadata,
        )


def test_delete_one_individual_jackknife_reports_complete_uncertainty(tmp_path):
    samples = _samples()
    sequences = {"chr1": "A" * 10}
    fasta = _write_fasta(tmp_path, sequences)
    bed = _write_bed(tmp_path, [("chr1", 0, 10)])
    genotypes = ["0/1"] + ["0/0"] * 19
    vcf = _write_vcf(tmp_path, sequences, samples, [("chr1", 2, "A", "G", genotypes)])
    result = _compute(
        tmp_path,
        vcf,
        fasta,
        samples,
        bed,
        sampling_unit_jackknife=True,
    )
    uncertainty = result["population_pi"]["uncertainty"]
    assert uncertainty["method"] == "delete_one_sampling_unit_jackknife"
    assert uncertainty["replicates"] == 20
    assert uncertainty["standard_error"] > 0
    assert len(uncertainty["leave_one_out_estimates"]) == 20
    assert 0 <= uncertainty["interval"][0] <= result["population_pi"]["point_estimate"]
    assert uncertainty["interval"][1] >= result["population_pi"]["point_estimate"]


@pytest.mark.parametrize(
    "metadata, message",
    [
        (_annotation_metadata(status="projected"), "requires native"),
        (_annotation_metadata(fasta_assembly="GCA_000002.1"), "assembly accession"),
    ],
)
def test_primary_4d_rejects_projected_or_wrong_assembly_annotation(tmp_path, metadata, message):
    samples = _samples()
    fasta = tmp_path / "truth.fa"
    fasta.write_bytes((FIXTURES / "truth.fa").read_bytes())
    bed = tmp_path / "callable.bed"
    bed.write_bytes((FIXTURES / "expected.callable.bed").read_bytes())
    sequences = {"chr1": "GCTAGCCAAAAAAAAAAAATCCACCTAAAAAAAAA"}
    vcf = _write_vcf(tmp_path, sequences, samples, [])
    with pytest.raises(Tier3ValidationError, match=message):
        _compute(
            tmp_path, vcf, fasta, samples, bed,
            gff_path=FIXTURES / "truth.gff3", annotation_metadata=metadata,
        )


def test_frozen_stratified_block_bootstrap_is_deterministic_and_reports_ratio_uncertainty(tmp_path):
    samples = _samples()
    sequences = {}
    gff_lines = ["##gff-version 3"]
    rows = []
    intervals = []
    calls = ["0/1"] * 20
    for index in range(40):
        contig = "c{:02d}".format(index)
        sequence = "GCC" if index % 2 == 0 else "GCT"
        sequences[contig] = sequence
        gff_lines.extend([
            "##sequence-region {} 1 3".format(contig),
            "{}\tx\tgene\t1\t3\t.\t+\t.\tID=g{}".format(contig, index),
            "{}\tx\tmRNA\t1\t3\t.\t+\t.\tID=t{};Parent=g{};tag=canonical".format(contig, index, index),
            "{}\tx\tCDS\t1\t3\t.\t+\t0\tParent=t{}".format(contig, index),
        ])
        rows.append((contig, 3, sequence[2], "A" if sequence[2] != "A" else "G", calls))
        intervals.append((contig, 0, 3))
    fasta = _write_fasta(tmp_path, sequences)
    gff = tmp_path / "native.gff3"
    gff.write_text("\n".join(gff_lines) + "\n", encoding="utf-8")
    vcf = _write_vcf(tmp_path, sequences, samples, rows)
    bed = _write_bed(tmp_path, intervals)
    kwargs = dict(
        gff_path=gff,
        annotation_metadata=_annotation_metadata(),
        minimum_4d_class_sites=20,
        bootstrap_replicates=64,
        dataset_id="bootstrap-fixture",
    )
    first = _compute(tmp_path, vcf, fasta, samples, bed, **kwargs)
    second = _compute(tmp_path, vcf, fasta, samples, bed, **kwargs)
    ratio = first["pi_S_over_pi_W"]
    assert ratio == second["pi_S_over_pi_W"]
    assert ratio["S"]["callable_sites"] == ratio["W"]["callable_sites"] == 20
    assert ratio["point_estimate"] == pytest.approx(1.0)
    assert ratio["bootstrap"]["eligible_blocks"] == 40
    assert ratio["bootstrap"]["interval"] == pytest.approx([1.0, 1.0])
    assert ratio["bootstrap"]["rng"] == "sha256-counter-v1"
    assert len(ratio["bootstrap"]["seed_digest"]) == 64


def test_sfs_B_is_absent_even_when_optional_gate_is_requested(tmp_path):
    samples = _samples()
    sequences = {"chr1": "AA"}
    fasta = _write_fasta(tmp_path, sequences)
    bed = _write_bed(tmp_path, [("chr1", 0, 2)])
    vcf = _write_vcf(tmp_path, sequences, samples, [])
    result = compute_population_pi(
        dataset_id="no-polarization",
        vcf_path=vcf,
        fasta_path=fasta,
        selected_samples=samples,
        design="wild_diploid",
        denominator_kind="cohort_callable_mask",
        callable_bed_path=bed,
        bootstrap_replicates=1,
        polarization_gate={
            "frozen_outgroup": True,
            "sample_size": 20,
            "polarization_error_model": "proposed",
            "demographic_model": "proposed",
        },
    )
    assert "polarized_sfs_B" not in result
    assert result["polarization_gate"] == {
        "requested": True, "passed": False, "reason": "deferred_by_tier3-decisions-v1"
    }
    assert "SFS" not in " ".join(result.keys())


@pytest.mark.parametrize("pilot,design", [("dgrp", "inbred_lines_haploidized"), ("ag1000g", "wild_diploid")])
def test_pilot_collection_is_idempotent_and_has_provenance(tmp_path, pilot, design):
    if not shutil.which("bcftools"):
        pytest.skip("bcftools is supplied by the pinned pure Guix environment")
    samples = _samples(22)
    sequences = {"chr1": "A"}
    fasta = _write_fasta(tmp_path, sequences)
    vcf = _write_vcf(tmp_path, sequences, samples, [])
    metadata = tmp_path / "samples.tsv"
    metadata.write_text(
        "sample_id\tpopulation_id\tqc_pass\tunrelated\n"
        + "".join("{}\tlocality\ttrue\ttrue\n".format(sample) for sample in samples),
        encoding="utf-8",
    )
    bed = _write_bed(tmp_path, [("chr1", 0, 1)])
    arguments = dict(
        dataset_id=pilot + "-fixture",
        input_vcf=vcf,
        fasta_path=fasta,
        sample_metadata_path=metadata,
        output_dir=tmp_path / ("out-" + pilot),
        design=design,
        denominator_kind="cohort_callable_mask",
        callable_bed_path=bed,
        population_id=None,
        assembly_accession="GCA_000001.1",
        pilot=pilot,
    )
    first = collect_population_vcf(**arguments)
    bcf_mtime = Path(first["outputs"]["normalized_bcf"]["path"]).stat().st_mtime_ns
    second = collect_population_vcf(**arguments)
    assert first["idempotent_reuse"] is False
    assert second["idempotent_reuse"] is True
    assert Path(second["outputs"]["normalized_bcf"]["path"]).stat().st_mtime_ns == bcf_mtime
    assert first["selection"]["selected_sample_list_sha256"] == first["denominator"]["sample_list_sha256"]
    assert first["raw_read_mapping_or_joint_calling"] == "not_implemented_by_policy"
    assert first["outputs"]["normalized_bcf"]["sha256"] == second["outputs"]["normalized_bcf"]["sha256"]


def test_hand_calculation_has_independent_bcftools_vcftools_crosscheck(tmp_path):
    if not shutil.which("bcftools") or not shutil.which("vcftools"):
        pytest.skip("cross-check tools are supplied by the pinned pure Guix environment")
    samples = _samples()
    sequences = {"chr1": "A" * 10}
    fasta = _write_fasta(tmp_path, sequences)
    genotypes = ["0/1"] + ["0/0"] * 19
    vcf = _write_vcf(tmp_path, sequences, samples, [("chr1", 2, "A", "G", genotypes)])
    bed = _write_bed(tmp_path, [("chr1", 0, 10)])
    result = _compute(tmp_path, vcf, fasta, samples, bed)

    stats = subprocess.run(
        ["bcftools", "stats", str(vcf)], check=True, text=True, stdout=subprocess.PIPE
    ).stdout
    assert any(line.startswith("SN\t0\tnumber of SNPs:\t1") for line in stats.splitlines())
    prefix = tmp_path / "vcftools"
    subprocess.run(
        ["vcftools", "--vcf", str(vcf), "--site-pi", "--out", str(prefix)],
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    values = (tmp_path / "vcftools.sites.pi").read_text(encoding="utf-8").splitlines()
    site_pi = float(values[1].split("\t")[-1])
    assert site_pi == pytest.approx(0.05)
    assert site_pi / 10 == pytest.approx(result["population_pi"]["point_estimate"])
