"""Unit tests for analysis/vgp_symmetric_read_test.py (synthetic fixtures only)."""

import csv
import json
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from analysis import vgp_symmetric_read_test as vsrt  # noqa: E402


def write_paf(tmp_path, rows):
    """rows: list of (q, qlen, qs, qe, strand, t, tlen, ts, te, cigar)."""
    path = tmp_path / "chain.paf"
    with path.open("w") as handle:
        for q, qlen, qs, qe, strand, t, tlen, ts, te, cigar in rows:
            handle.write(
                f"{q}\t{qlen}\t{qs}\t{qe}\t{strand}\t{t}\t{tlen}\t{ts}\t{te}\t"
                f"{qe - qs}\t{te - ts}\t255\tcg:Z:{cigar}\n"
            )
    return path


@pytest.fixture
def simple_chain(tmp_path):
    # Forward row A[100,300) -> B[500,700) with a 3bp insertion after 50 M
    # and a 2bp deletion after another 40 M: query 100 maps to target 500.
    # '+' row:  50M 3I 40M 2D 105M  (query span 100+198=198? see test)
    # query span: 50+3+40+105 = 198 -> [100,298)
    # target span: 50+40+2+105 = 197 -> [500,697)
    return write_paf(tmp_path, [
        ("A", 1000, 100, 298, "+", "B", 1200, 500, 697, "50M3I40M2D105M"),
    ])


class TestCigarWalk:
    def test_forward_plus_strand_plain_match(self, simple_chain):
        index, _ = vsrt.load_paf(simple_chain)
        # A:100 (0-based) -> B:500
        assert index.lift_site("A", 100, forward=True) == ("B", 500)
        # A:149 -> B:549 (50M block)
        assert index.lift_site("A", 149, forward=True) == ("B", 549)
        # A:150 is inside the 3I -> unmappable
        assert index.lift_site("A", 150, forward=True) is None
        assert index.lift_site("A", 151, forward=True) is None
        assert index.lift_site("A", 152, forward=True) is None
        # A:153 -> after 3I: target continues 550
        assert index.lift_site("A", 153, forward=True) == ("B", 550)
        # A:192 (40 more M) -> B:589; then 2D: A:193 -> B:592
        assert index.lift_site("A", 192, forward=True) == ("B", 589)
        assert index.lift_site("A", 193, forward=True) == ("B", 592)

    def test_reverse_plus_strand_round_trip(self, simple_chain):
        index, _ = vsrt.load_paf(simple_chain)
        # target-side gaps (2D region) are unmappable from B
        assert index.lift_site("B", 589, forward=False) == ("A", 192)
        assert index.lift_site("B", 590, forward=False) is None  # inside 2D
        assert index.lift_site("B", 592, forward=False) == ("A", 193)
        assert index.round_trip_site("A", 193) == ("B", 592)

    def test_forward_minus_strand(self, tmp_path):
        # '-' row: C[50,100) -> D[10,60), CIGAR 50M on revcomp(C) vs forward D
        paf = write_paf(tmp_path, [
            ("C", 200, 50, 100, "-", "D", 300, 10, 60, "50M"),
        ])
        index, _ = vsrt.load_paf(paf)
        # C:99 (last base) maps to revcomp offset 0 -> D:10
        assert index.lift_site("C", 99, forward=True) == ("D", 10)
        # C:50 maps to revcomp offset 49 -> D:59
        assert index.lift_site("C", 50, forward=True) == ("D", 59)
        # reverse: D:10 -> C:99
        assert index.lift_site("D", 10, forward=False) == ("C", 99)
        assert index.round_trip_site("C", 50) == ("D", 59)

    def test_allele_complement_on_minus(self):
        assert vsrt.complement_alleles("ACGT", "TGCA") == ("TGCA", "ACGT")
        assert vsrt.complement_alleles("A", "G") == ("T", "C")


class TestReciprocalRejection:
    def test_inconsistent_inverse_rejected(self, tmp_path):
        # Forward row E[0,100)->F[0,100) 50M... and an inverse-ish second row
        # that maps F[0,100) back to E[200,300): round trips must fail.
        paf = write_paf(tmp_path, [
            ("E", 500, 0, 100, "+", "F", 500, 0, 100, "100M"),
            ("F", 500, 0, 100, "+", "E", 500, 200, 300, "100M"),
        ])
        index, _ = vsrt.load_paf(paf)
        # forward lift works, but round trip is inconsistent
        assert index.lift_site("E", 10, forward=True) == ("F", 10)
        assert index.round_trip_site("E", 10) is None
        assert index.round_trip_site("F", 10) is None

    def test_overlapping_positions_unliftable(self, tmp_path):
        paf = write_paf(tmp_path, [
            ("E", 500, 0, 100, "+", "F", 500, 0, 100, "100M"),
            ("E", 500, 50, 150, "+", "G", 500, 0, 100, "100M"),
        ])
        index, _ = vsrt.load_paf(paf)
        # unique zone of E lifts normally
        assert index.lift_site("E", 10, forward=True) == ("F", 10)
        assert index.round_trip_site("E", 10) == ("F", 10)
        # overlap zone [50,100) is ambiguous in both directions
        assert index.lift_site("E", 60, forward=True) is None
        assert index.lift_site("F", 60, forward=False) is None
        assert index.round_trip_site("E", 60) is None
        # interval straddling the overlap is rejected wholesale
        assert index.lift_interval("E", 10, 80, forward=True) is None
        # clean interval still round-trips
        assert index.round_trip_interval("E", 10, 40, forward=True) == ("F", 10, 40)

    def test_interval_partial_rejection(self, simple_chain):
        index, _ = vsrt.load_paf(simple_chain)
        # interval crossing the 3I: both endpoints lift, but the pair straddles
        # a query-only gap; endpoints A[149,153): 149->B:549, 152->None
        assert index.lift_interval("A", 149, 153, forward=True) is None
        # fully contained clean interval round-trips exactly
        assert index.round_trip_interval("A", 120, 140, forward=True) == ("B", 520, 540)
        # interval extending past the row end is rejected
        assert index.round_trip_interval("A", 290, 320, forward=True) is None


class TestIntervalLiftInteriorRows:
    """Regression tests for the backward-only scan blocker: interval queries
    must see rows that START inside the interval, not only rows that start at
    or before it."""

    def test_interior_overlap_zone_rejects_interval(self, tmp_path):
        # E[0,100)->F[0,100) plus interior E[50,60)->G[0,10): lifting E[10,80)
        # must be REJECTED (two candidate rows intersect the interval), never
        # silently lifted whole through row F.
        paf = write_paf(tmp_path, [
            ("E", 500, 0, 100, "+", "F", 500, 0, 100, "100M"),
            ("E", 500, 50, 60, "+", "G", 500, 0, 10, "10M"),
        ])
        index, _ = vsrt.load_paf(paf)
        assert index.lift_interval("E", 10, 80, forward=True) is None
        assert index.round_trip_interval("E", 10, 80, forward=True) is None
        # a sub-interval clear of the overlap zone still lifts
        assert index.lift_interval("E", 10, 40, forward=True) == ("F", 10, 40)

    def test_seam_stitch_rejects_interval(self, tmp_path):
        # E[0,100)->F[0,100) and E[150,250)->F[120,220): E[80,200) straddles
        # the seam between two rows and must be REJECTED, never stitched into
        # a single F interval (the backward-only scan missed the second row).
        paf = write_paf(tmp_path, [
            ("E", 500, 0, 100, "+", "F", 500, 0, 100, "100M"),
            ("E", 500, 150, 250, "+", "F", 500, 120, 220, "100M"),
        ])
        index, _ = vsrt.load_paf(paf)
        assert index.lift_interval("E", 80, 200, forward=True) is None
        assert index.round_trip_interval("E", 80, 200, forward=True) is None
        # each side of the seam alone still lifts correctly
        assert index.lift_interval("E", 60, 90, forward=True) == ("F", 60, 90)
        assert index.lift_interval("E", 160, 200, forward=True) == ("F", 130, 170)


class TestTransferCommand:
    def test_transfer_end_to_end(self, tmp_path, simple_chain):
        bed = tmp_path / "pi.callable.bed"
        bed.write_text("A\t120\t140\nA\t150\t151\n")  # second sits in the 3I gap
        sites = tmp_path / "sites.tsv"
        sites.write_text(
            "chrom\tposition_1based\tref\talt\n"
            "A\t121\tC\tT\n"     # lifts cleanly
            "A\t151\tG\tA\n"     # inside 3I -> rejected
        )
        out = tmp_path / "transfer"
        result = subprocess.run(
            [sys.executable, "-m", "analysis.vgp_symmetric_read_test", "transfer",
             "--paf", str(simple_chain), "--bed", str(bed), "--sites", str(sites),
             "--target-label", "h2", "--output-dir", str(out)],
            check=True, capture_output=True, text=True,
        )
        lifted_bed = (out / "h2.pi.callable.bed").read_text().splitlines()
        assert lifted_bed == ["B\t520\t540"]
        lifted_sites = (out / "h2.assembly.snps.tsv").read_text().splitlines()
        assert lifted_sites == ["chrom\tposition_1based\tref\talt", "B\t521\tC\tT"]
        stats = json.loads((out / "transfer_stats.json").read_text())
        assert stats["bed"]["input_intervals"] == 2
        assert stats["bed"]["lifted_intervals"] == 1
        assert stats["sites"]["total"] == 2
        assert stats["sites"]["lifted"] == 1
        assert stats["roundtrip_fidelity"]["fidelity"] == 0.5
        # symmetric H1-frame subset mirrors exactly what lifted
        assert (out / "h1.symmetric.pi.callable.bed").read_text().splitlines() == ["A\t120\t140"]
        assert (out / "h1.symmetric.assembly.snps.tsv").read_text().splitlines() == [
            "chrom\tposition_1based\tref\talt", "A\t121\tC\tT"]

    def test_transfer_minus_strand_complements(self, tmp_path):
        paf = write_paf(tmp_path, [("C", 200, 50, 100, "-", "D", 300, 10, 60, "50M")])
        bed = tmp_path / "pi.callable.bed"
        bed.write_text("C\t60\t70\n")
        sites = tmp_path / "sites.tsv"
        sites.write_text("chrom\tposition_1based\tref\talt\nC\t61\tAC\tGT\n")
        out = tmp_path / "transfer"
        subprocess.run(
            [sys.executable, "-m", "analysis.vgp_symmetric_read_test", "transfer",
             "--paf", str(paf), "--bed", str(bed), "--sites", str(sites),
             "--target-label", "h2", "--output-dir", str(out)],
            check=True, capture_output=True, text=True,
        )
        # C:60 -> D:49; C:69 -> D:40  => interval D[40,50)
        assert (out / "h2.pi.callable.bed").read_text().splitlines() == ["D\t40\t50"]
        # site C:60 (0-based) -> D:49; alleles complement AC->TG, GT->CA
        assert (out / "h2.assembly.snps.tsv").read_text().splitlines() == [
            "chrom\tposition_1based\tref\talt", "D\t50\tTG\tCA",
        ]

    def test_transfer_rejection_reasons(self, tmp_path):
        paf = write_paf(tmp_path, [
            ("E", 500, 0, 100, "+", "F", 500, 0, 100, "100M"),
            ("E", 500, 50, 150, "+", "G", 500, 0, 100, "100M"),
        ])
        bed = tmp_path / "pi.callable.bed"
        bed.write_text("E\t10\t20\n")
        sites = tmp_path / "sites.tsv"
        sites.write_text(
            "chrom\tposition_1based\tref\talt\n"
            "E\t20\tC\tT\n"    # 0-based 19: unique zone -> lifted
            "E\t80\tA\tG\n"    # 0-based 79: overlap zone [50,100) -> overlap
            "E\t300\tA\tG\n"   # 0-based 299: no covering row -> unaligned
        )
        out = tmp_path / "transfer"
        subprocess.run(
            [sys.executable, "-m", "analysis.vgp_symmetric_read_test", "transfer",
             "--paf", str(paf), "--bed", str(bed), "--sites", str(sites),
             "--target-label", "h2", "--output-dir", str(out)],
            check=True, capture_output=True, text=True,
        )
        stats = json.loads((out / "transfer_stats.json").read_text())
        assert stats["sites"]["total"] == 3
        assert stats["sites"]["lifted"] == 1
        assert stats["sites"]["rejected"] == 2
        assert stats["sites"]["rejection_reasons"] == {"overlap": 1, "unaligned": 1}

    def test_transfer_bp_delta_internal_indel(self, tmp_path, simple_chain):
        # A[149,193) spans the 3I query-only gap: endpoints lift (149->B:549,
        # 192->B:589) but the destination interval is 41bp vs 44bp source,
        # so bp_delta_source_minus_dest must be exactly +3.
        bed = tmp_path / "pi.callable.bed"
        bed.write_text("A\t149\t193\n")
        sites = tmp_path / "sites.tsv"
        sites.write_text("chrom\tposition_1based\tref\talt\nA\t170\tC\tT\n")
        out = tmp_path / "transfer"
        subprocess.run(
            [sys.executable, "-m", "analysis.vgp_symmetric_read_test", "transfer",
             "--paf", str(simple_chain), "--bed", str(bed), "--sites", str(sites),
             "--target-label", "h2", "--output-dir", str(out)],
            check=True, capture_output=True, text=True,
        )
        stats = json.loads((out / "transfer_stats.json").read_text())
        assert stats["bed"]["lifted_intervals"] == 1
        assert stats["bed"]["bp_delta_source_minus_dest"] == 3
        assert (out / "h2.pi.callable.bed").read_text().splitlines() == ["B\t549\t590"]


class TestMetricsCommand:
    def make_frame(self, tmp_path, name, sites, reads, evidence):
        sites_path = tmp_path / f"{name}.sites.tsv"
        sites_path.write_text("chrom\tposition_1based\tref\talt\n" + "".join(sites))
        reads_path = tmp_path / f"{name}.reads.tsv"
        reads_path.write_text(
            "chrom\tposition_1based\tref\talt\tquality\tgenotype\n" + "".join(reads))
        evidence_path = tmp_path / f"{name}.evidence.tsv"
        evidence_path.write_text(
            "chrom\tposition_1based\tclassification\n" + "".join(evidence))
        bed_path = tmp_path / f"{name}.callable.bed"
        bed_path.write_text("chr1\t0\t3000000\n")
        return sites_path, reads_path, evidence_path, bed_path

    def test_metrics_math_and_bins(self, tmp_path):
        # 3 assembly sites at 1Mb-bin 0; classifications: supported, contradicted, ambiguous
        sites = ["chr1\t100\tA\tG\n", "chr1\t101\tA\tG\n", "chr1\t102\tA\tG\n"]
        reads = [
            "chr1\t100\tA\tG\t50\t0/1\n",    # het, at assembly site
            "chr1\t101\tA\tA\t50\t0/0\n",    # hom — filtered from read_het
            "chr1\t2000000\tC\tT\t40\t0/1\n",  # het NOT at assembly site (direction B), bin 2
            "chr1\t50\tT\tC\t10\t0/1\n",     # het but low QUAL (still counted in read_het)
        ]
        evidence = [
            "chr1\t100\tsupported_heterozygous\n",
            "chr1\t101\tcontradicted_homozygous_reference\n",
            "chr1\t102\tambiguous\n",
        ]
        sites_path, reads_path, evidence_path, bed_path = self.make_frame(
            tmp_path, "f1", sites, reads, evidence)
        out = tmp_path / "metrics.json"
        subprocess.run(
            [sys.executable, "-m", "analysis.vgp_symmetric_read_test", "metrics",
             "--frame", "h1", "--platform", "illumina",
             "--assembly-sites", str(sites_path), "--read-variants", str(reads_path),
             "--assembly-evidence", str(evidence_path),
             "--frame-callable-bed", str(bed_path),
             "--callable-bp", "1000", "--bin-size", "1000000",
             "--min-sites-per-bin", "1", "--output", str(out)],
            check=True, capture_output=True, text=True,
        )
        payload = json.loads(out.read_text())
        # read_het = sites 100, 2000000, 50 -> 3 het (101 is hom)
        assert payload["read_het_snps_on_mask"] == 3
        assert payload["pi_read"] == pytest.approx(3 / 1000)
        assert payload["concordance"]["supported_heterozygous"] == 1
        assert payload["concordance"]["contradicted_homozygous_reference"] == 1
        # B: het at QUAL>=30 not at assembly sites: only chr1:2000000
        # (chr1:50 is het but QUAL 10 < 30; chr1:100 is an assembly site)
        assert payload["both_direction_contradictions"]["B_readHet_assemblyHomRef"] == 1
        # bin 0 (sites at 100-102 + het at 50); bin 1 holds the 0-based pos
        # 1999999 (1-based 2000000): [1000000, 2000000)
        bins = (tmp_path / "metrics.bins.tsv").read_text().splitlines()
        assert any(row.startswith("chr1\t0\t1000000") for row in bins)
        assert any(row.startswith("chr1\t1000000\t2000000") for row in bins)
        # bin 0: 3 assembly sites, 1 contradicted -> rate 1/3; genome-wide 1/3
        # with flag factor 2.0 nothing is flagged; verify flagged empty:
        assert payload["flagged_bins"] == []

    def test_metrics_bin_flagging(self, tmp_path):
        # 4 sites in bin 0: 2 contradicted (rate .5); 4 in bin 1: 0 contradicted
        sites = [f"chr1\t{p}\tA\tG\n" for p in (101, 102, 103, 104, 2000101, 2000102, 2000103, 2000104)]
        reads = [f"chr1\t{p}\tA\tG\t50\t0/1\n" for p in (101, 102, 103, 104, 2000101, 2000102, 2000103, 2000104)]
        evidence = (
            ["chr1\t101\tcontradicted_homozygous_reference\n",
             "chr1\t102\tcontradicted_homozygous_reference\n",
             "chr1\t103\tsupported_heterozygous\n",
             "chr1\t104\tsupported_heterozygous\n"] +
            [f"chr1\t{p}\tsupported_heterozygous\n" for p in (2000101, 2000102, 2000103, 2000104)]
        )
        sites_path, reads_path, evidence_path, bed_path = self.make_frame(tmp_path, "f2", sites, reads, evidence)
        out = tmp_path / "metrics2.json"
        subprocess.run(
            [sys.executable, "-m", "analysis.vgp_symmetric_read_test", "metrics",
             "--frame", "h2", "--platform", "hifi",
             "--assembly-sites", str(sites_path), "--read-variants", str(reads_path),
             "--assembly-evidence", str(evidence_path),
             "--frame-callable-bed", str(bed_path),
             "--callable-bp", "100", "--bin-size", "1000000",
             "--min-sites-per-bin", "1", "--flag-factor", "1.5",
             "--output", str(out)],
            check=True, capture_output=True, text=True,
        )
        payload = json.loads(out.read_text())
        # genome rate 2/8 = 0.25; bin 0 rate 0.5 >= 1.5*0.25 -> flagged
        assert payload["concordance"]["contradiction_rate_resolved"] == pytest.approx(0.25)
        flagged = payload["flagged_bins"]
        assert len(flagged) == 1
        assert flagged[0]["chrom"] == "chr1"
        assert flagged[0]["bin_start_0based"] == 0


class TestReportCommand:
    def test_report_merges_frames(self, tmp_path):
        def frame(name, frame_id, platform, pi, a, b):
            payload = {
                "schema_version": vsrt.SCHEMA_METRICS,
                "frame": frame_id, "platform": platform,
                "callable_bp": 1000, "read_het_snps_on_mask": pi,
                "pi_read": pi / 1000,
                "assembly_sites_observed": 10,
                "evidence_classifications": {},
                "concordance": {"supported_heterozygous": 5,
                                "contradicted_homozygous_reference": a,
                                "contradiction_rate_resolved": a / 10},
                "both_direction_contradictions": {
                    "A_assemblyHet_readsHomRef": a, "B_readHet_assemblyHomRef": b},
                "flagged_bins": [{"chrom": "chrX", "bin_start_0based": 0,
                                  "bin_end_0based": 1000000, "contradiction_rate": 0.9,
                                  "assembly_sites": 60}],
                "flag_policy": {"bin_size": 1000000, "min_sites_per_bin": 50,
                                "flag_factor": 2.0, "genome_wide_rate": 0.25},
            }
            path = tmp_path / name
            path.write_text(json.dumps(payload))
            return path

        m1 = frame("m1.json", "h1", "illumina", 3, 2, 4)
        m2 = frame("m2.json", "h2", "illumina", 2, 1, 6)
        out = tmp_path / "report.md"
        subprocess.run(
            [sys.executable, "-m", "analysis.vgp_symmetric_read_test", "report",
             "--metrics-json", str(m1), "--metrics-json", str(m2),
             "--output", str(out)],
            check=True, capture_output=True, text=True,
        )
        text = out.read_text()
        assert "| h1 | illumina |" in text
        assert "| h2 | illumina |" in text
        assert "chrX" in text
        assert "Symmetric read-vs-assembly validation" in text


class TestAlleleAwareEvidence:
    """Regression tests for the deep-dive finding A parser defect: '.'/','
    pileup symbols must resolve against the pileup REF column (the frame
    fasta base), never unconditionally to the site's H1 allele."""

    def test_lifted_frame_ref_symbol_resolves_to_h2(self):
        # Site H1 allele A, H2 allele G. Lifted h2 frame: pileup REF (the
        # H2 fasta base) equals G — the ALT of the other frame. Three
        # ref-symbol reads '.'/',' are H2-matching reads.
        counts = vsrt.parse_pileup_bases_allele_aware(".,.", "G", "A", "G")
        assert counts["h1"] == 0
        assert counts["h2"] == 3
        assert counts["frame_ref_symbol"] == 3
        assert vsrt.classify_allele_evidence(3, counts,
                                             minimum_depth=1, maximum_depth=120) == "h2_only"
        # The legacy parser would have called all three reads H1 support.
        from analysis import vgp_read_validation as legacy
        legacy_counts = legacy.parse_pileup_bases(".,.", "A", "G")
        assert legacy_counts["ref"] == 3  # the defect, documented

    def test_h1_frame_balanced(self):
        counts = vsrt.parse_pileup_bases_allele_aware(".G.GG.", "A", "A", "G")
        assert counts["h1"] == 3
        assert counts["h2"] == 3
        assert (vsrt.classify_allele_evidence(6, counts,
                                              minimum_depth=1, maximum_depth=80)
                == "balanced_heterozygous")

    def test_minus_strand_complemented_alleles(self):
        # minus-strand lift: sites TSV already carries complemented alleles
        # (ref=T, alt=C); pileup REF is the h2-frame fasta base T.
        counts = vsrt.parse_pileup_bases_allele_aware(",,CCC,", "T", "T", "C")
        assert counts["h1"] == 3
        assert counts["h2"] == 3
        assert (vsrt.classify_allele_evidence(5, counts,
                                              minimum_depth=1, maximum_depth=80)
                == "balanced_heterozygous")

    def test_frame_ref_carries_neither_allele(self):
        # allele-origin disagreement: frame fasta base T matches neither the
        # H1 allele (A) nor the H2 allele (G).
        counts = vsrt.parse_pileup_bases_allele_aware(".T.", "T", "A", "G")
        assert counts["h1"] == 0
        assert counts["h2"] == 0
        assert counts["frame_ref_mismatch"] == 2
        assert counts["other"] == 1
        assert (vsrt.classify_allele_evidence(3, counts,
                                              minimum_depth=1, maximum_depth=80)
                == "not_observed")

    def test_legacy_classification_mapping(self):
        assert vsrt.legacy_classification_of("balanced_heterozygous") == "supported_heterozygous"
        assert (vsrt.legacy_classification_of("h1_only")
                == "contradicted_homozygous_reference")
        assert vsrt.legacy_classification_of("h2_only") == "ambiguous"
        assert vsrt.legacy_classification_of("skewed") == "ambiguous"
        assert vsrt.legacy_classification_of("outside_depth_mask") == "outside_depth_mask"

    def test_evidence_command_records_allele_origin_metadata(self, tmp_path):
        sites = tmp_path / "sites.tsv"
        sites.write_text("chrom\tposition_1based\tref\talt\n"
                         "chr1\t100\tA\tG\n"
                         "chr1\t101\tA\tG\n"
                         "chr1\t102\tA\tG\n")
        pileup = tmp_path / "pileup.txt"
        pileup.write_text("chr1\t100\tG\t3\t.,.\tIII\n"
                          "chr1\t101\tA\t6\t.G.GG.\tIIIIII\n"
                          "chr1\t102\tT\t3\t.T.\tIII\n")
        out = tmp_path / "evidence.tsv"
        summary = tmp_path / "evidence.json"
        subprocess.run(
            [sys.executable, "-m", "analysis.vgp_symmetric_read_test", "evidence",
             "--sites", str(sites), "--pileup", str(pileup),
             "--frame", "h2", "--platform", "hifi",
             "--minimum-depth", "1", "--maximum-depth", "80",
             "--output", str(out), "--summary", str(summary)],
            check=True, capture_output=True, text=True)
        rows = list(csv.DictReader(out.open(), delimiter="\t"))
        by_pos = {int(row["position_1based"]): row for row in rows}
        # lifted-frame hom-H2 site: '.' reads are H2 support, never H1
        assert by_pos[100]["classification"] == "h2_only"
        assert by_pos[100]["h1_reads"] == "0"
        assert by_pos[100]["h2_reads"] == "3"
        assert by_pos[100]["frame_ref_base"] == "G"
        assert by_pos[100]["site_h1_allele"] == "A"
        assert by_pos[100]["site_h2_allele"] == "G"
        assert by_pos[100]["allele_origin"] == "h2_allele"
        assert by_pos[100]["legacy_classification"] == "ambiguous"
        assert by_pos[101]["classification"] == "balanced_heterozygous"
        assert by_pos[101]["legacy_classification"] == "supported_heterozygous"
        assert by_pos[102]["allele_origin"] == "neither"
        payload = json.loads(summary.read_text())
        assert payload["classifications"]["h2_only"] == 1
        assert payload["classifications"]["balanced_heterozygous"] == 1
        assert payload["allele_origin"]["h2_allele"] == 1
        assert payload["allele_origin"]["neither"] == 1

    def test_metrics_accepts_v2_labels_and_v1_labels(self, tmp_path):
        # the same underlying evidence expressed with v2 vs v1 labels must
        # map to identical metrics numbers (old promoted files keep parsing)
        def run(name, labels):
            sites = [f"chr1\t{p}\tA\tG\n" for p in (101, 102, 103, 104)]
            reads = [f"chr1\t{p}\tA\tG\t50\t0/1\n" for p in (101, 102, 103, 104)]
            evidence = [f"chr1\t{p}\t{label}\n"
                        for p, label in zip((101, 102, 103, 104), labels)]
            (tmp_path / f"{name}.sites.tsv").write_text(
                "chrom\tposition_1based\tref\talt\n" + "".join(sites))
            (tmp_path / f"{name}.reads.tsv").write_text(
                "chrom\tposition_1based\tref\talt\tquality\tgenotype\n" + "".join(reads))
            (tmp_path / f"{name}.evidence.tsv").write_text(
                "chrom\tposition_1based\tclassification\n" + "".join(evidence))
            (tmp_path / f"{name}.bed").write_text("chr1\t0\t3000000\n")
            out = tmp_path / f"{name}.metrics.json"
            subprocess.run(
                [sys.executable, "-m", "analysis.vgp_symmetric_read_test", "metrics",
                 "--frame", "h1", "--platform", "hifi",
                 "--assembly-sites", str(tmp_path / f"{name}.sites.tsv"),
                 "--read-variants", str(tmp_path / f"{name}.reads.tsv"),
                 "--assembly-evidence", str(tmp_path / f"{name}.evidence.tsv"),
                 "--frame-callable-bed", str(tmp_path / f"{name}.bed"),
                 "--callable-bp", "100", "--output", str(out)],
                check=True, capture_output=True, text=True)
            return json.loads(out.read_text())

        v2 = run("v2", ["balanced_heterozygous", "h1_only", "h2_only", "skewed"])
        v1 = run("v1", ["supported_heterozygous",
                        "contradicted_homozygous_reference",
                        "ambiguous", "ambiguous"])
        for payload in (v2, v1):
            assert payload["concordance"]["supported_heterozygous"] == 1
            assert payload["concordance"]["contradicted_homozygous_reference"] == 1
        assert v2["h2_only_reads_carry_h2_allele_exclusively"] == 1
        assert v1["h2_only_reads_carry_h2_allele_exclusively"] == 0  # invisible to v1


class TestTriage:
    """Assembly-only paralog triage on a synthetic three-bin genome."""

    BIN = 1000

    def build(self, tmp_path, bin_sites=(2, 30, 3), gc_fractions=(0.6, 0.3, 0.6),
             disagree_in_hotspot=True):
        """Contig h1c: three 1kb bins; site positions chosen inside callable
        windows [0,950), [1000,1950), [2000,2950). H2 agrees with ALT except
        in the hotspot bin, where it carries the H1 allele (disagreement)."""
        import random
        rng = random.Random(7)
        sequences = []
        for gc in gc_fractions:
            gc_count = int(gc * self.BIN)
            at_count = self.BIN - gc_count
            bases = ["G"] * gc_count + ["A"] * at_count
            rng.shuffle(bases)
            sequences.append(bases)
        h1 = [base for sequence in sequences for base in sequence]

        sites = []
        h2 = list(h1)
        hotspot_bins = {i for i, (n, gc) in enumerate(zip(bin_sites, gc_fractions))
                        if n >= 10}
        for bin_index, count in enumerate(bin_sites):
            for k in range(count):
                pos = bin_index * self.BIN + 100 + k * 8
                ref, alt = "A", "G"
                if h1[pos] != ref:
                    # pick the actual base pair deterministically
                    ref = h1[pos]
                    alt = {"A": "G", "G": "A", "C": "T", "T": "C"}[ref]
                sites.append((pos, ref, alt))
                h2[pos] = ref if (disagree_in_hotspot and bin_index in hotspot_bins) else alt
        h1_fa = tmp_path / "h1.fa"
        h2_fa = tmp_path / "h2.fa"
        for path, seq in ((h1_fa, h1), (h2_fa, h2)):
            with path.open("w") as handle:
                handle.write(">h1c\n" if path is h1_fa else ">h2c\n")
                for i in range(0, len(seq), 60):
                    handle.write("".join(seq[i:i + 60]) + "\n")
        paf = write_paf(tmp_path, [
            ("h1c", 3000, 0, 3000, "+", "h2c", 3000, 0, 3000, "3000M"),
        ])
        sites_tsv = tmp_path / "sites.tsv"
        sites_tsv.write_text(
            "chrom\tposition_1based\tref\talt\n"
            + "".join(f"h1c\t{pos + 1}\t{ref}\t{alt}\n" for pos, ref, alt in sites))
        callable_bed = tmp_path / "callable.bed"
        callable_bed.write_text("h1c\t0\t950\nh1c\t1000\t1950\nh1c\t2000\t2950\n")
        return h1_fa, h2_fa, paf, sites_tsv, callable_bed, sites

    def run_triage(self, tmp_path, fixture, name="t", **overrides):
        h1_fa, h2_fa, paf, sites_tsv, callable_bed, sites = fixture
        out_dir = tmp_path / f"triage-{name}"
        argv = [sys.executable, "-m", "analysis.vgp_symmetric_read_test", "triage",
                "--pair", "TEST", "--paf", str(paf),
                "--h1-fasta", str(h1_fa), "--h2-fasta", str(h2_fa),
                "--sites-tsv", str(sites_tsv), "--callable-bed", str(callable_bed),
                "--bin-size", str(self.BIN), "--min-bin-bases", "900",
                "--output-dir", str(out_dir)]
        for key, value in overrides.items():
            argv += [f"--{key.replace('_', '-')}", str(value)]
        subprocess.run(argv, check=True, capture_output=True, text=True)
        summary = json.loads((out_dir / "triage_summary.json").read_text())
        return summary, (out_dir / "triage_bins.tsv"), (out_dir / "paralog_bins.bed")

    def test_hotspot_flagged_and_corrected_pi(self, tmp_path):
        fixture = self.build(tmp_path)
        summary, bins_tsv, paralog_bed = self.run_triage(tmp_path, fixture)
        # bin divergence: 2.0, 30.0, 3.0 sites/kb -> median 3.0, threshold 6.0
        # GC: 0.6, 0.3, 0.6 -> median 0.6, threshold 0.58 -> only bin 1 flags
        assert summary["thresholds"]["divergence_threshold_sites_per_kb"] == pytest.approx(6.0)
        assert summary["bins"]["flagged_paralog"] == 1
        flagged = summary["flagged_bins"][0]
        assert (flagged["chrom"], flagged["start"], flagged["end"]) == ("h1c", 1000, 2000)
        assert paralog_bed.read_text() == "h1c\t1000\t2000\n"
        # 35 callable sites, 30 in the paralog bin; 2850 callable bp, 950 excluded
        assert summary["sites"]["callable_sites"] == 35
        assert summary["sites"]["callable_sites_in_paralog_bins"] == 30
        assert summary["pi"]["callable_bp"] == 2850
        assert summary["pi"]["callable_bp_in_paralog_bins"] == 950
        assert summary["pi"]["raw_pi"] == pytest.approx(35 / 2850)
        assert summary["pi"]["corrected_pi"] == pytest.approx(5 / 1900)
        assert summary["sites"]["artifact_site_fraction"] == pytest.approx(30 / 35)
        assert summary["pi"]["artifact_bp_fraction"] == pytest.approx(950 / 2850)

    def test_allele_origin_localizes_disagreement(self, tmp_path):
        fixture = self.build(tmp_path)
        summary, bins_tsv, _ = self.run_triage(tmp_path, fixture)
        origin = summary["allele_origin"]
        assert origin["paf_orientation"] == "h1_query"
        assert origin["h2_base_matches_alt"] == 5
        assert origin["h2_base_disagrees"] == 30
        assert origin["disagreement_fraction"] == pytest.approx(30 / 35)
        rows = list(csv.DictReader(bins_tsv.open(), delimiter="\t"))
        by_bin = {int(row["bin_start_0based"]): row for row in rows}
        assert float(by_bin[1000]["allele_origin_disagreement_fraction"]) == pytest.approx(1.0)
        assert float(by_bin[0]["allele_origin_disagreement_fraction"]) == pytest.approx(0.0)

    def test_strict_threshold_not_flagged_at_equality(self, tmp_path):
        # divergences 1, 2, 2, 4, 8 sites/kb -> median 2.0, threshold 4.0;
        # the bin at exactly 4.0 must NOT flag (strict >); 8.0 does.
        fixture = self.build(tmp_path, bin_sites=(1, 2, 2, 4, 8),
                             gc_fractions=(0.6, 0.6, 0.6, 0.3, 0.3))
        summary, _, paralog_bed = self.run_triage(tmp_path, fixture, name="edge")
        assert summary["bins"]["flagged_paralog"] == 1
        assert summary["flagged_bins"][0]["start"] == 4000
        assert paralog_bed.read_text() == "h1c\t4000\t5000\n"

    def test_fasta_scan_gc_and_random_access(self, tmp_path):
        h1_fa, h2_fa, _, _, _, sites = self.build(tmp_path)
        scan = vsrt.FastaScan(h1_fa, 1000)
        bins = scan.bins
        assert set(key[0] for key in bins) == {"h1c"}
        gc_fractions = [scan.bins[("h1c", i)][1] / scan.bins[("h1c", i)][0] for i in range(3)]
        assert gc_fractions[0] == pytest.approx(0.6, abs=0.01)
        assert gc_fractions[1] == pytest.approx(0.3, abs=0.01)
        for pos, ref, alt in sites[:20]:
            assert scan.base_at("h1c", pos) == ref
        scan.close()
