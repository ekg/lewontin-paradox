import hashlib
from pathlib import Path

import pytest

from analysis.tier3_common import Tier3ValidationError
from analysis.tier3c_ncbi_gc import (
    AssemblyCandidate,
    PILOT_RANGES,
    acquire_verified,
    analyze_dataset,
    atomic_write_json,
    ncbi_artifact_urls,
    parse_ncbi_esummary,
    rank_assemblies,
    select_assembly,
    validate_pilot,
)


def _artifact(path: Path) -> dict:
    data = path.read_bytes()
    return {
        "logical_name": path.name,
        "uri": path.resolve().as_uri(),
        "sha256": hashlib.sha256(data).hexdigest(),
        "size_bytes": len(data),
    }


def _dataset(tmp_path: Path, *, annotation=True, genetic_code=1, annotation_accession=None):
    fasta = tmp_path / "exact.fa"
    # tx_plus -> GCT GCC; tx_minus -> GGT GGA.  The rest is valid nuclear
    # sequence and contributes independently to whole-genome GC.
    fasta.write_text(
        ">chr1 exact nuclear chromosome\nGCTAGCCAAAAAAAAAAAATCCACCTAAAAAAAAA\n",
        encoding="utf-8",
    )
    dataset = {
        "dataset_id": "testus.exactus.tier3c",
        "species": {"scientific_name": "Testus exactus", "taxon_id": 1},
        "reference": {
            "assembly_accession": "GCF_000001.2",
            "provider": "NCBI RefSeq",
            "release": "2026-01-02",
            "fasta": _artifact(fasta),
            "nuclear_contigs_only": True,
        },
        "annotation": None,
        "notes": "synthetic exact-reference test",
    }
    if annotation:
        gff = tmp_path / "exact.gff3"
        gff.write_text(
            "##gff-version 3\n"
            "##sequence-region chr1 1 35\n"
            "chr1\ttest\tmRNA\t1\t7\t.\t+\t.\tID=tx_plus;Parent=gene_plus;tag=canonical\n"
            "chr1\ttest\tCDS\t1\t3\t.\t+\t0\tParent=tx_plus\n"
            "chr1\ttest\tCDS\t4\t7\t.\t+\t1\tParent=tx_plus\n"
            "chr1\ttest\tmRNA\t20\t26\t.\t-\t.\tID=tx_minus;Parent=gene_minus\n"
            "chr1\ttest\tCDS\t23\t26\t.\t-\t1\tParent=tx_minus\n"
            "chr1\ttest\tCDS\t20\t22\t.\t-\t0\tParent=tx_minus\n",
            encoding="utf-8",
        )
        mapping = tmp_path / "contigs.tsv"
        mapping.write_text("annotation_contig\tfasta_contig\nchr1\tchr1\n", encoding="utf-8")
        dataset["annotation"] = {
            "provider": "NCBI RefSeq",
            "release": "2026-01-02",
            "assembly_accession": annotation_accession or "GCF_000001.2",
            "status": "native",
            "file": _artifact(gff),
            "contig_mapping": _artifact(mapping),
            "genetic_code": genetic_code,
            "exact_reference_assertion": True,
        }
    return dataset


def _environment() -> dict:
    return {
        "manager": "gnu-guix",
        "guix_environment": "/gnu/store/00000000000000000000000000000000-tier3-test",
        "channel_commit": "44bbfc24e4bcc48d0e3343cd3d83452721af8c36",
        "tools": {
            "python3": {
                "version": "3.test",
                "executable": "/gnu/store/00000000000000000000000000000000-python/bin/python3",
                "store_path": "/gnu/store/00000000000000000000000000000000-python",
            }
        },
    }


def test_assembly_ranking_is_frozen_and_deterministic():
    candidates = [
        AssemblyCandidate("GCA_000003.1", 7, 7, "representative genome", "Complete Genome", "2025-01-01", "ftp://gca"),
        AssemblyCandidate("GCF_000002.1", 7, 7, "reference genome", "Chromosome", "2024-01-01", "https://gcf-old"),
        AssemblyCandidate("GCF_000002.2", 7, 7, "reference genome", "Chromosome", "2026-01-01", "https://gcf-new"),
        AssemblyCandidate("GCF_000001.1", 8, 7, "reference genome", "Complete Genome", "2026-01-01", "https://wrong-taxon"),
    ]
    assert [item.accession for item in rank_assemblies(candidates, taxon_id=7)] == [
        "GCF_000002.2",
        "GCF_000002.1",
        "GCA_000003.1",
        "GCF_000001.1",
    ]


def test_exact_accession_match_includes_prefix_and_version():
    candidates = [
        AssemblyCandidate("GCF_000001.1", 1, 1, "reference genome", "Chromosome", "2025-01-01", "https://one"),
        AssemblyCandidate("GCA_000001.2", 1, 1, "reference genome", "Chromosome", "2026-01-01", "https://two"),
        AssemblyCandidate("GCF_000001.2", 1, 1, "reference genome", "Chromosome", "2024-01-01", "https://exact"),
    ]
    assert select_assembly(candidates, taxon_id=1, exact_accession="GCF_000001.2").ftp_path == "https://exact"
    with pytest.raises(Tier3ValidationError, match="exact assembly accession"):
        select_assembly(candidates, taxon_id=1, exact_accession="GCF_000001.3")


def test_ncbi_esummary_preserves_gca_gcf_identity_and_builds_https_urls():
    payload = {
        "result": {
            "uids": ["1"],
            "1": {
                "assemblyaccession": "GCA_000001.2",
                "taxid": 1,
                "species_taxid": 1,
                "refseq_category": "na",
                "assemblystatus": "Scaffold",
                "asmreleasedate_genbank": "2026-02-03",
                "ftppath_genbank": "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000001.2_Test",
                "ftppath_refseq": "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000001.2_Test",
            },
        }
    }
    candidate = parse_ncbi_esummary(payload)[0]
    assert candidate.accession == "GCA_000001.2"
    urls = ncbi_artifact_urls(candidate)
    assert urls["fasta"].startswith("https://ftp.ncbi.nlm.nih.gov/")
    assert urls["fasta"].endswith("GCA_000001.2_Test_genomic.fna.gz")


def test_gc3_uses_canonical_cds_phase_and_reverse_strand(tmp_path):
    output = tmp_path / "result.json"
    result = analyze_dataset(_dataset(tmp_path), output, environment=_environment())
    assert result["gc3"] == {
        "status": "available",
        "value": pytest.approx(1 / 4),
        "gc_bases": 1,
        "callable_third_positions": 4,
        "genes": 2,
        "transcripts": 2,
        "terminal_stop_codons_excluded": 0,
    }
    assert result["whole_genome_gc"]["callable_bases"] == 35
    assert result["annotation_provenance"]["contig_mapping"] == {"chr1": "chr1"}
    assert result["annotation_provenance"]["all_retained_cds_validated"] is True
    assert result["annotation_provenance"]["fasta_sha256"] == result["reference"]["fasta_sha256"]
    assert result["annotation_provenance"]["sequence_regions"] == {"chr1": 35}
    assert result["reference"]["contig_dictionary"] == {"chr1": 35}


def test_gff3_edge_cases_and_contig_dictionary_fail_closed(tmp_path):
    dataset = _dataset(tmp_path)
    gff = Path(dataset["annotation"]["file"]["uri"].removeprefix("file://"))
    text = gff.read_text(encoding="utf-8").replace("##sequence-region chr1 1 35", "##sequence-region chr1 1 34")
    gff.write_text(text, encoding="utf-8")
    dataset["annotation"]["file"] = _artifact(gff)
    output = tmp_path / "should-not-exist.json"
    with pytest.raises(Tier3ValidationError, match="contig length mismatch"):
        analyze_dataset(dataset, output, environment=_environment())
    assert not output.exists()


def test_duplicate_gff3_cds_segment_fails_closed(tmp_path):
    dataset = _dataset(tmp_path)
    gff = Path(dataset["annotation"]["file"]["uri"].removeprefix("file://"))
    with gff.open("a", encoding="utf-8") as handle:
        handle.write("chr1\ttest\tCDS\t1\t3\t.\t+\t0\tParent=tx_plus\n")
    dataset["annotation"]["file"] = _artifact(gff)
    output = tmp_path / "result.json"
    with pytest.raises(Tier3ValidationError, match="duplicate CDS segment"):
        analyze_dataset(dataset, output, environment=_environment())
    assert not output.exists()


def test_unsupported_genetic_code_is_structured_missingness(tmp_path):
    result = analyze_dataset(_dataset(tmp_path, genetic_code=2), tmp_path / "result.json", environment=_environment())
    assert result["gc3"] == {
        "status": "unavailable",
        "reason": "unsupported_nuclear_genetic_code",
        "genetic_code": 2,
    }
    assert result["whole_genome_gc"]["status"] == "available"


def test_missing_annotation_preserves_exact_fasta_gc_only(tmp_path):
    result = analyze_dataset(_dataset(tmp_path, annotation=False), tmp_path / "result.json", environment=_environment())
    assert result["annotation_provenance"] is None
    assert result["gc3"] == {"status": "unavailable", "reason": "native_annotation_absent"}
    assert result["whole_genome_gc"]["status"] == "available"


def test_projected_annotation_never_fills_primary_gc3(tmp_path):
    dataset = _dataset(tmp_path)
    dataset["annotation"]["status"] = "projected"
    result = analyze_dataset(dataset, tmp_path / "result.json", environment=_environment())
    assert result["gc3"]["status"] == "unavailable"
    assert result["gc3"]["reason"] == "annotation_not_native"
    assert result["annotation_provenance"]["status"] == "projected"


def test_sampled_provider_cds_mismatch_blocks_annotation_result(tmp_path):
    dataset = _dataset(tmp_path)
    provider = tmp_path / "provider-cds.fa"
    provider.write_text(">tx_plus\nGCTGCA\n>tx_minus\nGGTGGA\n", encoding="utf-8")
    dataset["annotation"]["provider_cds"] = _artifact(provider)
    output = tmp_path / "result.json"
    with pytest.raises(Tier3ValidationError, match="sampled CDS mismatch"):
        analyze_dataset(dataset, output, environment=_environment())
    assert not output.exists()


def test_mismatched_fasta_gff_accessions_fail_before_output(tmp_path):
    output = tmp_path / "result.json"
    with pytest.raises(Tier3ValidationError, match="FASTA/GFF assembly accession mismatch"):
        analyze_dataset(
            _dataset(tmp_path, annotation_accession="GCF_000001.1"),
            output,
            environment=_environment(),
        )
    assert not output.exists()


def test_interrupted_download_never_publishes_partial_file(tmp_path):
    destination = tmp_path / "downloaded.fa.gz"

    class InterruptedResponse:
        def __enter__(self):
            return self

        def __exit__(self, *_args):
            return False

        def read(self, _size):
            raise OSError("connection reset")

    with pytest.raises(Tier3ValidationError, match="download failed"):
        acquire_verified(
            "https://example.invalid/genome.fa.gz",
            destination,
            expected_sha256="0" * 64,
            expected_size=1,
            opener=lambda *_args, **_kwargs: InterruptedResponse(),
        )
    assert not destination.exists()
    assert not list(tmp_path.glob(".downloaded.fa.gz.*"))


def test_completed_download_with_wrong_checksum_is_not_published(tmp_path):
    destination = tmp_path / "downloaded.fa.gz"

    class Response:
        def __init__(self):
            self.done = False

        def __enter__(self):
            return self

        def __exit__(self, *_args):
            return False

        def read(self, _size):
            if self.done:
                return b""
            self.done = True
            return b"wrong"

    with pytest.raises(Tier3ValidationError, match="SHA-256"):
        acquire_verified(
            "https://example.invalid/genome.fa.gz",
            destination,
            expected_sha256=hashlib.sha256(b"right").hexdigest(),
            expected_size=5,
            opener=lambda *_args, **_kwargs: Response(),
        )
    assert not destination.exists()


def test_atomic_idempotent_rerun_and_provenance_denominators(tmp_path):
    output = tmp_path / "result.json"
    first = analyze_dataset(_dataset(tmp_path), output, environment=_environment())
    stat = output.stat()
    second = analyze_dataset(_dataset(tmp_path), output, environment=_environment())
    assert second == first
    assert output.stat().st_mtime_ns == stat.st_mtime_ns
    assert first["reference"]["accession"] == "GCF_000001.2"
    assert len(first["reference"]["fasta_sha256"]) == 64
    assert len(first["annotation_provenance"]["gff_sha256"]) == 64
    assert first["environment"]["manager"] == "gnu-guix"
    assert first["environment"]["tools"]["python3"]["store_path"].startswith("/gnu/store/")
    assert first["notes"]

    # The low-level writer also leaves an already canonical result untouched.
    assert atomic_write_json(first, output) is False
    assert not list(tmp_path.glob(".result.json.*"))


@pytest.mark.parametrize(
    "species,gc3,genome_gc",
    [("Drosophila melanogaster", 0.55, 0.42), ("Homo sapiens", 0.52, 0.41)],
)
def test_predeclared_independent_pilot_ranges(species, gc3, genome_gc):
    assert species in PILOT_RANGES
    accession = PILOT_RANGES[species]["assembly_accession"]
    result = {
        "reference": {"accession": accession},
        "annotation_provenance": {"status": "native", "assembly_accession": accession},
        "gc3": {"status": "available", "value": gc3},
        "whole_genome_gc": {"status": "available", "value": genome_gc},
    }
    assert PILOT_RANGES[species]["source"]
    assert validate_pilot(species, result) == []
    result["gc3"]["value"] = 0.01
    assert validate_pilot(species, result)
