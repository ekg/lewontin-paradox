"""Unit tests for analysis/vgp_symmetric_read_test.py (synthetic fixtures only)."""

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
