import hashlib
import json
from pathlib import Path

import pytest

from analysis.tier3c_control_audit import AuditError, apply_to_collected, audit_control, combine


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_control_fixture(tmp_path: Path) -> tuple[Path, Path, Path]:
    fasta = tmp_path / "control.fna"
    # gene_gc is GCC (GC3=1).  The reverse complement of bases 7..15 is
    # GGTGGAGGA (GC3=0), making pooled-CDS-base GC3=1/4 but the unweighted
    # mean of the two per-gene GC3 values=1/2.
    fasta.write_text(">chr1\nGCCAAATCCTCCACC\n", encoding="utf-8")
    gff = tmp_path / "control.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "##sequence-region chr1 1 15\n"
        "chr1\tRefSeq\tmRNA\t1\t3\t.\t+\t.\tID=tx_gc;Parent=gene_gc;tag=MANE Select\n"
        "chr1\tRefSeq\tCDS\t1\t3\t.\t+\t0\tParent=tx_gc\n"
        "chr1\tRefSeq\tmRNA\t7\t15\t.\t-\t.\tID=tx_at;Parent=gene_at\n"
        "chr1\tRefSeq\tCDS\t7\t15\t.\t-\t0\tParent=tx_at\n"
        "chr1\tRefSeq\tmRNA\t1\t3\t.\t+\t.\tID=tx_pseudo;Parent=gene_pseudo;pseudo=true\n"
        "chr1\tRefSeq\tCDS\t1\t3\t.\t+\t0\tParent=tx_pseudo\n",
        encoding="utf-8",
    )
    mapping = tmp_path / "annotation-contig-map.tsv"
    mapping.write_text(
        "annotation_contig\tfasta_contig\nchr1\tchr1\n", encoding="utf-8"
    )
    return fasta, gff, mapping


def test_independent_audit_reproduces_pooled_and_states_gene_weighting(tmp_path):
    fasta, gff, mapping = _write_control_fixture(tmp_path)
    result = audit_control(
        dataset_id="test.control.tier3c",
        scientific_name="Test control",
        assembly_accession="GCF_000001.2",
        provider="NCBI RefSeq",
        release="2026-01-01",
        fasta_path=fasta,
        gff_path=gff,
        contig_mapping_path=mapping,
        expected_fasta_sha256=_sha256(fasta),
        expected_gff_sha256=_sha256(gff),
        expected_mapping_sha256=_sha256(mapping),
        genetic_code=1,
        production_gc3={
            "value": 0.25,
            "gc_bases": 1,
            "callable_third_positions": 4,
            "genes": 2,
        },
    )

    assert result["metrics"]["pooled_cds_base_gc3"] == {
        "definition": "sum(G_or_C_third_bases)/sum(callable_third_bases)",
        "value": 0.25,
        "gc_bases": 1,
        "callable_third_positions": 4,
    }
    assert result["metrics"]["gene_weighted_gc3"]["mean"] == 0.5
    assert result["metrics"]["gene_weighted_gc3"]["median"] == 0.5
    assert result["selection"]["retained_genes"] == 2
    assert result["selection"]["exclusions"] == {
        "annotated_translation_exception": 1,
        "gene_without_valid_cds": 1,
    }
    assert result["production_comparison"]["all_exact"] is True
    assert result["provenance"]["dictionary_validation"]["passed"] is True
    assert result["method"]["imports_production_estimator"] is False


def test_independent_audit_fails_closed_on_checksum_and_dictionary(tmp_path):
    fasta, gff, mapping = _write_control_fixture(tmp_path)
    common = dict(
        dataset_id="test.control.tier3c",
        scientific_name="Test control",
        assembly_accession="GCF_000001.2",
        provider="NCBI RefSeq",
        release="2026-01-01",
        fasta_path=fasta,
        gff_path=gff,
        contig_mapping_path=mapping,
        expected_fasta_sha256=_sha256(fasta),
        expected_gff_sha256=_sha256(gff),
        expected_mapping_sha256=_sha256(mapping),
        genetic_code=1,
        production_gc3=None,
    )
    with pytest.raises(AuditError, match="FASTA SHA-256"):
        audit_control(**{**common, "expected_fasta_sha256": "0" * 64})

    gff.write_text(
        gff.read_text(encoding="utf-8").replace(
            "##sequence-region chr1 1 15", "##sequence-region chr1 1 14"
        ),
        encoding="utf-8",
    )
    with pytest.raises(AuditError, match="sequence-region length"):
        audit_control(**{**common, "expected_gff_sha256": _sha256(gff)})


def test_independent_audit_rejects_unsupported_genetic_code(tmp_path):
    fasta, gff, mapping = _write_control_fixture(tmp_path)
    with pytest.raises(AuditError, match="genetic code 1"):
        audit_control(
            dataset_id="test.control.tier3c",
            scientific_name="Test control",
            assembly_accession="GCF_000001.2",
            provider="NCBI RefSeq",
            release="2026-01-01",
            fasta_path=fasta,
            gff_path=gff,
            contig_mapping_path=mapping,
            expected_fasta_sha256=_sha256(fasta),
            expected_gff_sha256=_sha256(gff),
            expected_mapping_sha256=_sha256(mapping),
            genetic_code=2,
            production_gc3=None,
        )


def test_posthoc_decision_preserves_failed_legacy_bands_and_applies_audit(tmp_path):
    inputs = []
    for slug, name, value in (
        ("drosophila.melanogaster", "Drosophila melanogaster", 0.631),
        ("homo.sapiens", "Homo sapiens", 0.585),
    ):
        path = tmp_path / f"{slug}.json"
        path.write_text(
            json.dumps(
                {
                    "schema_version": "tier3c-control-audit-v1",
                    "dataset_id": f"{slug}.tier3c",
                    "scientific_name": name,
                    "provenance": {
                        "annotation_status": "native",
                        "exact_reference_assertion": True,
                        "dictionary_validation": {"passed": True},
                    },
                    "selection": {"all_retained_cds_reconstructed": True},
                    "metrics": {
                        "pooled_cds_base_gc3": {"value": value},
                        "gene_weighted_gc3": {"mean": value, "median": value},
                    },
                    "production_comparison": {"all_exact": True},
                }
            ),
            encoding="utf-8",
        )
        inputs.append(path)
    combined = tmp_path / "combined.json"
    combine(inputs, combined)
    decision = json.loads(combined.read_text(encoding="utf-8"))
    assert decision["original_control_gate"]["passed"] is False
    assert decision["original_control_gate"]["bands_rewritten"] is False
    assert decision["audited_control_gate"]["passed"] is True
    assert decision["promotion_decision"]["decision"] == "promote_tier3c_composition"

    qc_dir = tmp_path / "qc"
    qc_dir.mkdir()
    for control in decision["controls"].values():
        (qc_dir / f"{control['dataset_id']}.json").write_text(
            json.dumps({"dataset_id": control["dataset_id"], "pilot_failures": ["legacy"]}),
            encoding="utf-8",
        )
    summary = tmp_path / "summary.json"
    summary.write_text(
        json.dumps({"completed": 135, "failures": 38, "control_gate_passed": False}),
        encoding="utf-8",
    )
    apply_to_collected(combined, qc_dir, summary)
    updated = json.loads(summary.read_text(encoding="utf-8"))
    assert updated["original_control_gate_passed"] is False
    assert updated["audited_control_gate_passed"] is True
    assert updated["control_gate_passed"] is True
    assert updated["promotion_decision"]["decision"] == "promote_tier3c_composition"
