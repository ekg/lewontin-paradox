import importlib.util
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "results/tier3a/acquisition_build_queries.py"
SPEC = importlib.util.spec_from_file_location("acquisition_build_queries", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)

SELECTOR_SCRIPT = ROOT / "results/tier3a/acquisition_select_impg_partitions.py"
SELECTOR_SPEC = importlib.util.spec_from_file_location("acquisition_select_impg_partitions", SELECTOR_SCRIPT)
SELECTOR = importlib.util.module_from_spec(SELECTOR_SPEC)
assert SELECTOR_SPEC.loader is not None
sys.modules[SELECTOR_SPEC.name] = SELECTOR
SELECTOR_SPEC.loader.exec_module(SELECTOR)


def test_parse_attributes_decodes_gff3_values():
    observed = MODULE.parse_attributes("ID=gene-a;description=alpha%2Cbeta;Parent=rna-1,rna-2")
    assert observed == {"ID": "gene-a", "description": "alpha,beta", "Parent": "rna-1,rna-2"}


def test_merge_only_overlapping_or_touching_intervals():
    assert MODULE.merge([(10, 20), (20, 30), (31, 40), (5, 8)]) == [
        (5, 8), (10, 30), (31, 40)
    ]


def test_transcript_phase_chain_accepts_all_three_phases():
    transcript = MODULE.Transcript(1, "rna-x", ("gene-x",), "chr1", 0, 100, "+", False)
    rows = [
        MODULE.CDS(2, "chr1", 0, 5, "+", "0", "cds-x", ("rna-x",), "p", "l"),
        MODULE.CDS(3, "chr1", 10, 12, "+", "1", "cds-x", ("rna-x",), "p", "l"),
        MODULE.CDS(4, "chr1", 20, 24, "+", "2", "cds-x", ("rna-x",), "p", "l"),
        MODULE.CDS(5, "chr1", 30, 33, "+", "1", "cds-x", ("rna-x",), "p", "l"),
        MODULE.CDS(6, "chr1", 40, 44, "+", "1", "cds-x", ("rna-x",), "p", "l"),
        MODULE.CDS(7, "chr1", 50, 55, "+", "0", "cds-x", ("rna-x",), "p", "l"),
    ]
    assert MODULE.transcript_phase_qc(transcript, rows, {"chr1": 100}) == (True, "passed")


def test_transcript_phase_chain_rejects_mismatch():
    transcript = MODULE.Transcript(1, "rna-x", ("gene-x",), "chr1", 0, 100, "+", False)
    rows = [
        MODULE.CDS(2, "chr1", 0, 5, "+", "0", "cds-x", ("rna-x",), "p", "l"),
        MODULE.CDS(3, "chr1", 10, 14, "+", "2", "cds-x", ("rna-x",), "p", "l"),
    ]
    assert MODULE.transcript_phase_qc(transcript, rows, {"chr1": 100}) == (
        False, "phase_chain_mismatch_expected_1_observed_2"
    )


def test_fully_covered_requires_one_gapless_mapping_union_interval():
    assert MODULE.fully_covered([(0, 20), (30, 50)], 5, 20)
    assert not MODULE.fully_covered([(0, 20), (30, 50)], 5, 35)


def test_read_h1_coverage_accepts_query_axis(tmp_path):
    paf = tmp_path / "mapping.paf"
    paf.write_text("h1\t100\t10\t90\t+\th2\t110\t20\t100\t70\t80\t60\tcg:Z:80=\n")
    assert MODULE.read_h1_coverage(paf, {"h1": 100}) == ({"h1": [(10, 90)]}, "query")


def test_haplotype_contig_map_preserves_pairs_and_unmapped_contigs(tmp_path):
    paf = tmp_path / "mapping.paf"
    paf.write_text(
        "h1a\t100\t10\t90\t+\th2a\t110\t20\t100\t70\t80\t60\n"
        "h1a\t100\t20\t80\t-\th2a\t110\t30\t90\t50\t60\t60\n"
    )
    output = tmp_path / "map.tsv"
    MODULE.write_haplotype_contig_map(
        paf, {"h1a": 100, "h1_unmapped": 25}, {"h2a": 110, "h2_unmapped": 30}, output
    )
    rows = output.read_text().splitlines()
    assert rows[1].split("\t")[1:10] == [
        "h1a", "100", "h2a", "110", "+,-", "2", "80", "80", "selected_bounded_mapping"
    ]
    assert rows[2].endswith("H1_no_selected_mapping")
    assert rows[3].endswith("H2_no_selected_mapping")


def test_haplotype_contig_map_rejects_axis_specific_length_mismatch(tmp_path):
    paf = tmp_path / "mapping.paf"
    paf.write_text("h1\t999\t10\t90\t+\th2\t110\t20\t100\t70\t80\t60\n")
    output = tmp_path / "map.tsv"
    try:
        MODULE.write_haplotype_contig_map(paf, {"h1": 100}, {"h2": 110}, output)
    except SystemExit as error:
        assert "length mismatch" in str(error)
    else:
        raise AssertionError("mismatched H1 length was accepted")


def test_impg_native_partition_selector_uses_half_open_overlap_and_stable_order():
    partitions = [
        ("chr2", 0, 10, "p4"),
        ("chr1", 20, 30, "p3"),
        ("chr1", 0, 10, "p1"),
        ("chr1", 10, 20, "p2"),
    ]
    targets = [
        ("chr1", 9, 11, "span_b"),
        ("chr1", 2, 5, "span_a"),
        ("chr1", 30, 40, "touches_but_does_not_overlap"),
    ]
    assert SELECTOR.select_partitions(partitions, targets) == [
        ("chr1", 0, 10, "p1", ["span_a", "span_b"]),
        ("chr1", 10, 20, "p2", ["span_b"]),
    ]
