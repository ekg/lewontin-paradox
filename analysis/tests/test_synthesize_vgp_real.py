from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import pytest

from analysis import synthesize_vgp_real as synthesis


ROOT = Path(__file__).parents[2]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


@pytest.fixture(scope="module")
def generated(tmp_path_factory):
    output = tmp_path_factory.mktemp("vgp-real-synthesis")
    manifest = synthesis.generate(output)
    return output, manifest


def test_real_estimates_replace_the_old_zero_estimate_synthesis(generated):
    output, manifest = generated
    observed = manifest["biological_evidence"]["assembly_derived_diversity"]
    assert observed == {
        "completed_pair_count": 2,
        "selection_ids": ["P04", "P07"],
        "minimum_pi": pytest.approx(0.0021472198306856562),
        "maximum_pi": pytest.approx(0.004604184795871289),
        "max_to_min_ratio": pytest.approx(2.1442540396066794),
        "validated_quantitative_pair_count": 0,
        "retained_raw_pending_pair_count": 1,
        "invalidated_by_raw_reads_pair_count": 1,
    }
    assert manifest["closed_world"]["audited_pairs"] == 10
    assert manifest["closed_world"]["verified_core_results"] == 2
    assert "produced two assembly-derived individual diversity estimates" in (
        output / synthesis.REPORT.name
    ).read_text().lower()


def test_pair_table_preserves_confidence_and_validation_dispositions(generated):
    output, _ = generated
    pairs = read_tsv(output / synthesis.PAPER_PAIRS.name)
    assert len(pairs) == 10
    by_id = {row["selection_id"]: row for row in pairs}

    assert by_id["P04"]["assembly_pi"] == "0.004604184795871289"
    assert by_id["P04"]["synthesis_confidence_tier"] == "T1_RETAINED_RAW_VALIDATION_PENDING"
    assert by_id["P04"]["quantitative_disposition"] == "RETAINED_ASSEMBLY_DERIVED_PENDING_RAW_VALIDATION"
    assert by_id["P07"]["assembly_pi"] == "0.0021472198306856562"
    assert by_id["P07"]["synthesis_confidence_tier"] == "T2_INVALIDATED_BY_EXACT_READS"
    assert by_id["P07"]["quantitative_disposition"] == "PRESERVED_NOT_ADMITTED_QUANTITATIVELY"
    assert by_id["P07"]["primary_read_pi"] == "0.00038064456641950404"
    assert by_id["P07"]["primary_read_over_assembly_pi"] == "0.4363316186618393"

    assert {by_id[p]["synthesis_confidence_tier"] for p in ("P01", "P02", "P03", "P05")} == {
        "X_HARD_INVALID_PRIMARY"
    }
    assert {by_id[p]["synthesis_confidence_tier"] for p in ("P06", "P09", "P10")} == {
        "X_EXECUTION_ERROR_NO_ESTIMATE"
    }
    assert by_id["P08"]["synthesis_confidence_tier"] == "P_RESUMABLE_RUNNING_AT_FREEZE"
    assert all(row["canonical_vgp_root"] == "/moosefs/erikg/vgp" for row in pairs)


def test_psmc_is_descriptive_unscaled_and_nonindependent(generated):
    output, manifest = generated
    histories = read_tsv(output / synthesis.PSMC_HISTORIES.name)
    assert {row["selection_id"] for row in histories} == {"P04", "P07"}
    assert {row["pi_psmc_independence"] for row in histories} == {
        "SAME_PAIR_NONINDEPENDENT"
    }
    by_id = {row["selection_id"]: row for row in histories}
    assert by_id["P04"]["primary_theta_0_per_100bp_bin"] == "0.52445"
    assert by_id["P04"]["trajectory_lambda_min"] == "0.492194"
    assert by_id["P04"]["trajectory_lambda_max"] == "2.874536"
    assert by_id["P04"]["quantitative_disposition"] == "RETAINED_DESCRIPTIVE_UNSCALED_RAW_PENDING"
    assert by_id["P07"]["trajectory_lambda_max"] == "159.317994"
    assert by_id["P07"]["quantitative_disposition"] == "INVALIDATED_BY_EXACT_READ_VALIDATION"
    assert {row["absolute_history_status"] for row in histories} == {
        "BOUNDED_9_GENERIC_SCENARIOS_NOT_SPECIES_CALIBRATION"
    }
    assert manifest["analysis_contract"]["same_pair_pi_psmc_independent"] is False
    assert manifest["analysis_contract"]["species_calibrated_absolute_psmc"] is False
    assert manifest["analysis_contract"]["population_inference_authorized"] is False


def test_annotation_is_exact_native_but_not_gene_conversion_evidence(generated):
    output, _ = generated
    annotations = read_tsv(output / synthesis.ANNOTATION_TABLE.name)
    assert len(annotations) == 6
    assert {row["selection_id"] for row in annotations} == {"P07"}
    assert {row["sequence_dictionary_equal"] for row in annotations} == {"true"}
    assert {row["quantitative_disposition"] for row in annotations} == {
        "DESCRIPTIVE_ONLY_PARENT_PAIR_INVALIDATED"
    }
    assert {row["gene_conversion_disposition"] for row in annotations} == {
        "NOT_A_CONFORMING_GENE_CONVERSION_ESTIMATE"
    }

    branches = read_tsv(output / synthesis.GENE_CONVERSION.name)
    assert [row["branch"] for row in branches] == [
        "direct_pedigree_or_gamete",
        "population_allele_frequency_spectrum",
        "historical_phylogenetic_substitution",
        "non_allelic_paralog",
    ]
    assert {row["vgp_integration_status"] for row in branches} == {
        "SEPARATE_NO_ACTUAL_CONFORMING_VGP_ESTIMATE"
    }
    assert {row["estimate"] for row in branches} == {"NOT_ESTIMABLE"}


def test_supported_bounded_suggestive_and_unidentifiable_claims_are_explicit(generated):
    output, _ = generated
    claims = read_tsv(output / synthesis.CLAIMS.name)
    assert {row["classification"] for row in claims} == {
        "supported",
        "bounded",
        "suggestive",
        "unidentifiable",
    }
    assert len(claims) >= 14
    for row in claims:
        assert row["canonical_vgp_root"] == "/moosefs/erikg/vgp"
        assert row["evidence_artifacts"]
        assert row["sampling_unit"]
        assert row["uncertainty_covariance"]
        assert row["forbidden_inference"]
    lewontin = next(row for row in claims if row["claim_id"] == "LR-IMPLICATION")
    assert lewontin["classification"] == "unidentifiable"
    assert "does not test" in lewontin["conclusion"]


def test_every_scheduler_record_and_failure_is_reconciled(generated):
    output, manifest = generated
    jobs = read_tsv(output / synthesis.JOB_LEDGER.name)
    assert len(jobs) == 658
    assert len({row["job_id"] for row in jobs}) == 658
    assert sum(row["packet"] == "scale_vgp_real" for row in jobs) == 650
    assert sum(row["packet"] == "validate_vgp_pilot_reads" for row in jobs) == 8
    assert manifest["jobs"]["scale_vgp_real"] == {
        "allocations": 650,
        "states": {
            "CANCELLED by 1001": 346,
            "COMPLETED": 240,
            "FAILED": 42,
            "PENDING": 21,
            "RUNNING": 1,
        },
    }
    assert manifest["jobs"]["validate_vgp_pilot_reads"] == {
        "allocations": 8,
        "states": {
            "CANCELLED by 1001": 2,
            "COMPLETED": 2,
            "FAILED": 4,
        },
    }
    assert sum(row["failure_or_nonterminal"] == "true" for row in jobs) == 416
    assert all(row["canonical_vgp_root"] == "/moosefs/erikg/vgp" for row in jobs)


def test_inputs_environments_and_digests_are_closed_world(generated):
    output, manifest = generated
    digests = read_tsv(output / synthesis.DIGEST_LEDGER.name)
    repository = [row for row in digests if row["binding_class"] == "REPOSITORY_INPUT"]
    embedded = [row for row in digests if row["binding_class"] == "UPSTREAM_EMBEDDED_DIGEST"]
    assert len(repository) == len(synthesis.INPUT_PATHS)
    assert embedded
    assert {row["algorithm"] for row in digests} == {"sha256", "md5"}
    assert all(
        len(row["digest"]) == (64 if row["algorithm"] == "sha256" else 32)
        for row in digests
    )
    assert {row["verification"] for row in repository} == {"PASS_REHASHED"}
    assert manifest["environment"]["guix_channel_commit"] == "44bbfc24e4bcc48d0e3343cd3d83452721af8c36"
    assert manifest["environment"]["core_guix_closure_sha256"] == "8fcdb32021f1cd8eac839509cff47ab6bdd63b656b30e243fdf78d3c4ba24f9d"
    assert manifest["environment"]["read_validation_guix_closure_sha256"] == "ac0cb3601e56ef62b9ef99419de3659b2a2ba59b2aead29bc5f1928b50c83da2"
    assert manifest["raw_reads"]["verified_or_reused_objects"] == 4
    assert manifest["raw_reads"]["verified_or_reused_bytes"] == 31_058_137_613
    assert manifest["raw_reads"]["planned_pending_bytes"] == 42_344_746_693
    assert manifest["raw_reads"]["quarantined_failed_transfer_count"] == 1
    assert manifest["canonical_vgp_root"] == "/moosefs/erikg/vgp"


def test_static_figures_and_committed_packet_reproduce(generated):
    output, manifest = generated
    for name in (synthesis.DIVERSITY_FIGURE.name, synthesis.PSMC_FIGURE.name):
        text = (output / name).read_text(encoding="utf-8")
        assert text.startswith("<svg")
        assert "<title>" in text and "<desc>" in text
        assert "/moosefs/erikg/vgp" in text
        assert "P07" in text and "P04" in text

    assert synthesis.validate_outputs(output) == []
    for relative, digest in manifest["output_digests"].items():
        assert hashlib.sha256((output / Path(relative).name).read_bytes()).hexdigest() == digest

    committed = json.loads(synthesis.MANIFEST.read_text(encoding="utf-8"))
    assert committed == manifest
    for relative, digest in committed["output_digests"].items():
        assert hashlib.sha256((ROOT / relative).read_bytes()).hexdigest() == digest


def test_validator_rejects_promotion_of_invalidated_p07(generated, tmp_path):
    output, _ = generated
    for name in synthesis.OUTPUT_FILENAMES:
        (tmp_path / name).write_bytes((output / name).read_bytes())
    rows = read_tsv(tmp_path / synthesis.PAPER_PAIRS.name)
    next(row for row in rows if row["selection_id"] == "P07")["quantitative_disposition"] = "ADMITTED"
    with (tmp_path / synthesis.PAPER_PAIRS.name).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    errors = synthesis.validate_outputs(tmp_path, verify_digests=False)
    assert any("P07 invalidated quantitative disposition" in error for error in errors)
