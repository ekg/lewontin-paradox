from pathlib import Path

import pytest

from analysis import review_vgp_real_pilot as review


def test_variant_recount_reapplies_callable_and_indel_flank_policy(tmp_path: Path):
    vcf = tmp_path / "variants.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t2\t.\tA\tG\t.\tPASS\t.\n"
        "chr1\t10\t.\tA\tAT\t.\tPASS\t.\n"
        "chr1\t20\t.\tC\tT\t.\tPASS\t.\n"
        "chr1\t30\t.\tG\tA\t.\tPASS\t.\n"
        "chr1\t32\t.\tA\tT\t.\tPASS\t.\n",
        encoding="utf-8",
    )
    result = review.audit_variants(vcf, [("chr1", 0, 30)], indel_flank=2)
    assert result == {
        "method": "two-pass VCF recount plus BED containment and independent merged +/-2-bp callable indel mask",
        "normalized_variant_records": 5,
        "normalized_snp_records": 4,
        "normalized_indel_records": 1,
        "pre_indel_callable_bp": 30,
        "callable_variant_records": 4,
        "callable_snp_records": 3,
        "callable_indel_records": 1,
        "callable_snps_masked_by_indel_flank": 0,
        "indel_masked_h1_bp": 5,
        "final_callable_bp": 25,
        "heterozygous_snps": 3,
        "pi": 3 / 25,
    }


def test_variant_recount_masks_snp_overlapping_indel_window(tmp_path: Path):
    vcf = tmp_path / "variants.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t9\t.\tA\tG\t.\tPASS\t.\n"
        "chr1\t10\t.\tA\tAT\t.\tPASS\t.\n",
        encoding="utf-8",
    )
    result = review.audit_variants(vcf, [("chr1", 0, 20)], indel_flank=2)
    assert result["callable_snps_masked_by_indel_flank"] == 1
    assert result["heterozygous_snps"] == 0
    assert result["indel_masked_h1_bp"] == 5


@pytest.mark.parametrize(
    ("value", "expected"),
    [("16G", 16), ("134217728K", 128), ("4096M", 4), ("1T", 1024)],
)
def test_request_memory_gib_parses_sacct_units(value: str, expected: float):
    assert review.request_memory_gib(value) == expected


def test_allocation_summary_never_imputes_missing_maxrss():
    rows = [
        {
            "State": "COMPLETED",
            "ElapsedRaw": "3600",
            "AllocCPUS": "4",
            "ReqMem": "16G",
            "MaxRSS": "",
            "CPUTimeRAW": "14400",
        }
    ]
    result = review.allocation_summary(rows)
    assert result["allocated_core_hours"] == 4
    assert result["requested_memory_gib_hours"] == 16
    assert result["maxrss_observed_rows"] == 0


def test_parse_psmc_uses_final_native_iteration(tmp_path: Path):
    psmc = tmp_path / "test.psmc"
    psmc.write_text(
        "RD\t1\nTR\t1.0\t2.0\nRS\t0\t0.0\t2.0\n"
        "RD\t2\nTR\t3.5\t4.0\nRS\t0\t0.0\t1.0\nRS\t1\t0.1\t0.5\n",
        encoding="utf-8",
    )
    rows, theta = review.parse_psmc(psmc)
    assert theta == 3.5
    assert [row["interval"] for row in rows] == [0, 1]


def test_parse_psmc_rejects_nonfinite_output(tmp_path: Path):
    psmc = tmp_path / "test.psmc"
    psmc.write_text("RD\t1\nTR\t1.0\t2.0\nRS\t0\t0.0\tnan\n", encoding="utf-8")
    with pytest.raises(review.ReviewError, match="non-finite"):
        review.parse_psmc(psmc)
