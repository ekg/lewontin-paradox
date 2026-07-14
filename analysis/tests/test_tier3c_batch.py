import gzip
from pathlib import Path

import pytest

import analysis.tier3c_batch as tier3c_batch
from analysis.tier3_common import Tier3ValidationError
from analysis.tier3c_batch import (
    BATCH_SCHEMA,
    MINIMUM_NUCLEAR_REFERENCE_BASES,
    FAILURE_COLUMNS,
    MANIFEST_COLUMNS,
    RESULT_COLUMNS,
    _atomic_tsv,
    _bgzip_replace,
    _excluded_non_nuclear,
    _filter_fasta,
    _filter_gff,
    _validate_nuclear_reference_size,
    _slug,
    atomic_write_json,
    validate_collected,
)


def test_nuclear_filter_uses_assembly_unit_and_preserves_unplaced_contigs(tmp_path):
    report = tmp_path / "assembly_report.txt"
    report.write_text(
        "# header\n"
        "1\tassembled-molecule\t1\tChromosome\tGB1.1\t=\tNC_1.1\tPrimary Assembly\t6\tchr1\n"
        "MT\tassembled-molecule\tMT\tMitochondrion\tGBM.1\t=\tNC_M.1\tnon-nuclear\t6\tchrM\n"
        "scaffold1\tunplaced-scaffold\tna\tna\tGBS.1\t=\tNW_1.1\tPrimary Assembly\t4\tna\n"
        "mito\tassembled-molecule\tmitochondrial_genome\tChromosome\tGBM2.1\t=\tNC_M2.1\tPrimary Assembly\t6\tna\n",
        encoding="utf-8",
    )
    excluded, counts = _excluded_non_nuclear(report)
    assert {"MT", "GBM.1", "NC_M.1", "chrM"} <= excluded
    assert "NW_1.1" not in excluded
    assert counts == {"assembly_unit_non_nuclear": 1, "assigned_molecule_non_nuclear": 1}

    source = tmp_path / "genomic.fna.gz"
    with gzip.open(source, "wt", encoding="utf-8") as handle:
        handle.write(">NC_1.1 chromosome\nACGTAC\n>NC_M.1 mitochondrion\nAAAAAA\n>NW_1.1 scaffold\nCCCC\n")
    destination = tmp_path / "nuclear.fna"
    assert _filter_fasta(source, destination, excluded) == 2
    assert destination.read_text(encoding="utf-8") == (
        ">NC_1.1 chromosome\nACGTAC\n>NW_1.1 scaffold\nCCCC\n"
    )


def test_nuclear_reference_minimum_rejects_symbiont_sized_host_assembly(tmp_path):
    fai = tmp_path / "nuclear.fna.gz.fai"
    fai.write_text("symbiont_contig\t1819390\t0\t0\t0\n", encoding="utf-8")
    with pytest.raises(Tier3ValidationError, match=str(MINIMUM_NUCLEAR_REFERENCE_BASES)):
        _validate_nuclear_reference_size(fai)


def test_gff_filter_removes_non_nuclear_dictionary_and_features(tmp_path):
    source = tmp_path / "genomic.gff.gz"
    with gzip.open(source, "wt", encoding="utf-8") as handle:
        handle.write(
            "##gff-version 3\n"
            "##sequence-region NC_1.1 1 6\n"
            "##sequence-region NC_M.1 1 6\n"
            "NC_1.1\tRefSeq\tgene\t1\t3\t.\t+\t.\tID=g1\n"
            "NC_M.1\tRefSeq\tgene\t1\t3\t.\t+\t.\tID=gm\n"
            "##FASTA\n>ignored\nAAA\n"
        )
    destination = tmp_path / "nuclear.gff3"
    regions, features = _filter_gff(source, destination, {"NC_M.1"})
    assert regions == [("NC_1.1", 6)]
    assert features == 1
    assert "NC_M.1" not in destination.read_text(encoding="utf-8")
    assert "##FASTA" not in destination.read_text(encoding="utf-8")


def test_dataset_slugs_are_stable_and_shell_safe():
    assert _slug("Drosophila melanogaster") == "drosophila.melanogaster"
    assert _slug("Aedes albopictus") == "aedes.albopictus"


def test_bgzip_staging_derivative_is_byte_identical_on_rerun(tmp_path):
    source = tmp_path / "nuclear.fna"
    source.write_text(">chr1\n" + "ACGT" * 1000 + "\n", encoding="utf-8")
    compressed = _bgzip_replace(source)
    first = compressed.read_bytes()
    source.write_text(">chr1\n" + "ACGT" * 1000 + "\n", encoding="utf-8")
    assert _bgzip_replace(source).read_bytes() == first


def test_collected_empty_schema_and_summary_validate(tmp_path):
    _atomic_tsv(tmp_path / "tier3c_manifest.tsv", MANIFEST_COLUMNS, [])
    _atomic_tsv(tmp_path / "tier3c_data.tsv", RESULT_COLUMNS, [])
    _atomic_tsv(tmp_path / "tier3c_failure_ledger.tsv", FAILURE_COLUMNS, [])
    atomic_write_json({"completed": 0, "failures": 0}, tmp_path / "tier3c_qc_summary.json")
    assert validate_collected(tmp_path) == {
        "manifest_rows": 0,
        "result_rows": 0,
        "failure_rows": 0,
        "schema_and_checksum_fields_valid": True,
    }


def test_failed_rerun_removes_stale_scientific_result(tmp_path, monkeypatch):
    batch = tmp_path / "batch.json"
    atomic_write_json(
        {
            "schema_version": BATCH_SCHEMA,
            "datasets": [{"dataset_id": "species.tier3c", "species": {"scientific_name": "Species"}}],
            "environment": {},
        },
        batch,
    )
    output = tmp_path / "results"
    output.mkdir()
    stale = output / "species.tier3c.json"
    stale.write_text('{"stale":true}\n', encoding="utf-8")

    def fail(*_args, **_kwargs):
        raise RuntimeError("current run failed")

    monkeypatch.setattr(tier3c_batch, "analyze_dataset", fail)
    with pytest.raises(RuntimeError, match="current run failed"):
        tier3c_batch.run_one(batch, 0, output)
    assert not stale.exists()
    assert (output / "species.tier3c.failure.json").is_file()
