import hashlib
import json
from pathlib import Path

import pytest

from analysis.vgp_10_pilot import PilotError, sha256_file
from analysis.vgp_three_pair import (
    audit_graph_identifiers,
    compare_controlled_pafs,
    materialize_exact_staged_fastas,
)


ROOT = Path(__file__).parents[2]


def _asset(path: Path, records: dict[str, str]) -> dict[str, object]:
    return {
        "path": str(path),
        "size_bytes": path.stat().st_size,
        "sha256": sha256_file(path),
        "sequence_dictionary": [
            {
                "name": name,
                "length": len(sequence),
                "md5": hashlib.md5(sequence.encode()).hexdigest(),
            }
            for name, sequence in records.items()
        ],
    }


def _fixture(tmp_path: Path):
    h1_records = {"H1.chr": "ACGTN" * 7, "H1.extra": "TGCA" * 5}
    h2_records = {"H2.chr": "ACGTA" * 7, "H2.extra": "TGCA" * 5}
    h1, h2 = tmp_path / "source.h1.fa", tmp_path / "source.h2.fa"
    h1_sequence = h1_records["H1.chr"]
    h1.write_text(
        ">H1.chr description\n" + h1_sequence[:7].lower() + "\n" + h1_sequence[7:].lower() +
        "\n>H1.extra\n" + h1_records["H1.extra"] + "\n"
    )
    h2.write_text("".join(f">{name}\n{sequence}\n" for name, sequence in h2_records.items()))
    manifest = tmp_path / "manifest.json"
    manifest.write_text(json.dumps({
        "selection_id": "P02",
        "assets": {"h1_fasta": _asset(h1, h1_records), "h2_fasta": _asset(h2, h2_records)},
    }))
    return manifest, h1_records, h2_records


def test_exact_staging_rewraps_but_preserves_every_digest_and_universe(tmp_path):
    manifest, h1_records, h2_records = _fixture(tmp_path)
    staged, audit = tmp_path / "staged", tmp_path / "dictionary.json"
    value = materialize_exact_staged_fastas(manifest, staged, audit)
    assert value["logical_sequences_equal_to_frozen_source"] is True
    assert value["roles"]["h1_fasta"]["record_count"] == len(h1_records)
    assert value["roles"]["h2_fasta"]["record_count"] == len(h2_records)
    assert (staged / "h1.fa").read_text().startswith(">H1.chr\nACGTN")
    assert (staged / "h1_universe.bed").read_text().splitlines() == [
        "H1.chr\t0\t35", "H1.extra\t0\t20"
    ]
    assert json.loads(audit.read_text())["h1_universe_sha256"] == sha256_file(
        staged / "h1_universe.bed"
    )


def test_graph_id_census_validates_orientation_lengths_coordinates_and_no_omission(tmp_path):
    manifest, _, _ = _fixture(tmp_path)
    staged, dictionary = tmp_path / "staged", tmp_path / "dictionary.json"
    materialize_exact_staged_fastas(manifest, staged, dictionary)
    paf, partitions, output = tmp_path / "map.paf", tmp_path / "partitions.bed", tmp_path / "audit.json"
    paf.write_text("H2.chr\t35\t0\t20\t+\tH1.chr\t35\t0\t20\t20\t20\t60\tcg:Z:20=\n")
    partitions.write_text("H1.chr\t0\t20\tpart0\nH1.extra\t0\t20\tpart1\n")
    result = audit_graph_identifiers(dictionary, paf, partitions, output)
    assert result["unresolved_ids"] == 0
    assert result["silently_omitted_regions"] == 0
    assert result["partition_rows"] == 2

    partitions.write_text("missing\t0\t20\tpart0\n")
    with pytest.raises(PilotError, match="unresolved partition ID"):
        audit_graph_identifiers(dictionary, paf, partitions, output)


def test_alias_requires_matching_sequence_digest_and_source_fasta_digest(tmp_path):
    manifest, _, _ = _fixture(tmp_path)
    staged, dictionary = tmp_path / "staged", tmp_path / "dictionary.json"
    value = materialize_exact_staged_fastas(manifest, staged, dictionary)
    canonical = value["roles"]["h1_fasta"]["records"][0]
    paf, partitions, output, aliases = (
        tmp_path / "map.paf", tmp_path / "partitions.bed", tmp_path / "audit.json",
        tmp_path / "aliases.json",
    )
    paf.write_text("H2.chr\t35\t0\t20\t+\tH1.chr\t35\t0\t20\t20\t20\t60\tcg:Z:20=\n")
    partitions.write_text("graph-H1\t0\t20\tpart0\n")
    row = {
        "observed_id": "graph-H1", "canonical_id": "H1.chr", "length": 35,
        "sequence_sha256": canonical["sequence_sha256"], "source_fasta_sha256": "a" * 64,
    }
    aliases.write_text(json.dumps({"aliases": [row]}))
    result = audit_graph_identifiers(dictionary, paf, partitions, output, aliases)
    assert result["digest_validated_aliases_used"] == ["graph-H1"]
    row["sequence_sha256"] = "0" * 64
    aliases.write_text(json.dumps({"aliases": [row]}))
    with pytest.raises(PilotError, match="alias sequence digest mismatch"):
        audit_graph_identifiers(dictionary, paf, partitions, output, aliases)


def test_slurm_contract_has_dictionary_rebuild_and_pinned_wfmash_fallback():
    pair = (ROOT / "analysis/slurm/vgp_10_pilot/pair_stage.sh").read_text()
    mapping = (ROOT / "analysis/slurm/vgp_10_pilot/mapping_stage.sh").read_text()
    assert "stage-fastas" in pair and "audit-graph-ids" in pair
    assert 'h1_universe="$SLURM_TMPDIR/inputs/h1_universe.bed"' in pair
    assert "fastga|wfmash" in mapping
    assert "--aligner wfmash" in mapping
    assert "same_staged_fasta_bytes_as_corrected_retry" in mapping
    assert "enforce-paf" in mapping and "audit-paf" in mapping


def test_controlled_backend_comparison_uses_common_coverage_and_exact_alleles(tmp_path):
    h1, h2 = tmp_path / "h1.fa", tmp_path / "h2.fa"
    h1.write_text(">H1\n" + "A" * 100 + "\n")
    h2.write_text(">H2\n" + "A" * 20 + "C" + "A" * 79 + "\n")
    fastga, wfmash = tmp_path / "fastga.paf", tmp_path / "wfmash.paf"
    fastga.write_text("H2\t100\t0\t100\t+\tH1\t100\t0\t100\t99\t100\t60\tcg:Z:20=1X79=\n")
    wfmash.write_text("H2\t100\t10\t100\t+\tH1\t100\t10\t100\t89\t90\t60\tcg:Z:10=1X79=\n")
    result = compare_controlled_pafs(fastga, wfmash, h1, h2, tmp_path / "comparison")
    assert result["overlapping_target_bp"] == 90
    assert result["target_coverage_jaccard"] == 0.9
    assert result["shared_exact_variants"] == 1
    assert result["exact_variant_jaccard"] == 1.0
    assert result["fastga_ref_alt_reconstruction"]["reconstruction_failures"] == 0
