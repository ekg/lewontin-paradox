"""Tests for the corrected VGP evidence publication (correct-vgp-evidence).

These tests pin the correction contract: reviewed bounded values are carried
through unchanged, P07 is reclassified to unresolved read/assembly discordance,
historical failures are pipeline defects rather than biological exclusions, the
annotation catalog is published with explicit non-scale-out semantics, and
historical artifacts are preserved with explicit supersession links.
"""

import csv
import json
import re
import subprocess
import sys
from pathlib import Path

import pytest

from analysis import correct_vgp_evidence as cve

REPO_ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = REPO_ROOT / "analysis" / "vgp_evidence_correction_v1"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


@pytest.fixture(scope="module")
def execution() -> dict:
    return json.loads(
        (REPO_ROOT / "analysis" / "vgp_three_pair_execution_v2.json").read_text(encoding="utf-8")
    )


@pytest.fixture(scope="module")
def result_pairs() -> list[dict[str, str]]:
    return read_tsv(OUT_DIR / "result_pairs.tsv")


@pytest.fixture(scope="module")
def claims() -> list[dict[str, str]]:
    return read_tsv(OUT_DIR / "claims_ledger.tsv")


@pytest.fixture(scope="module")
def reliability() -> list[dict[str, str]]:
    return read_tsv(OUT_DIR / "pipeline_reliability.tsv")


@pytest.fixture(scope="module")
def annotation_paths() -> list[dict[str, str]]:
    return read_tsv(OUT_DIR / "annotation_paths.tsv")


@pytest.fixture(scope="module")
def supersession() -> list[dict[str, str]]:
    return read_tsv(OUT_DIR / "supersession_ledger.tsv")


@pytest.fixture(scope="module")
def report() -> str:
    return (OUT_DIR / "report.md").read_text(encoding="utf-8")


# ---------------------------------------------------------------------------
# Corrected core results
# ---------------------------------------------------------------------------


def test_accepted_rows_match_reviewed_execution(result_pairs, execution):
    by_id = {pair["selection_id"]: pair for pair in execution["pairs"]}
    accepted = [row for row in result_pairs if row["row_kind"] == "accepted_bounded_core"]
    assert {row["selection_id"] for row in accepted} == {"P07", "P03", "P02"}
    for row in accepted:
        pair = by_id[row["selection_id"]]
        diversity = pair["diversity"]
        psmc = pair["psmc"]
        assert int(row["heterozygous_snps"]) == diversity["heterozygous_snps"]
        assert int(row["callable_bp"]) == diversity["callable_bp"]
        assert float(row["pi"]) == pytest.approx(diversity["pi"], rel=0, abs=0)
        assert float(row["psmc_primary_theta"]) == psmc["primary_theta"]
        lo, hi = psmc["nearest_index_central_95pct"]
        assert row["psmc_nearest_index_central_95pct"] == f"[{lo}, {hi}]"
        assert int(row["psmc_finite_bootstraps"]) == 200


def test_clean_canary_control_row_preserved(result_pairs):
    controls = [row for row in result_pairs if row["row_kind"] == "clean_canary_reproduction_control"]
    assert len(controls) == 1
    row = controls[0]
    assert row["selection_id"] == "P07"
    assert int(row["heterozygous_snps"]) == 574122
    assert int(row["callable_bp"]) == 267379237
    assert float(row["pi"]) == pytest.approx(0.0021472198306856562)
    assert row["corrected_disposition"] == "RETAINED_REPRODUCTION_CONTROL"
    assert "NOT_ADMISSIBLE_CORE" in row["corrected_classification"]


def test_p04_retained_not_strengthened(result_pairs):
    p04 = [row for row in result_pairs if row["selection_id"] == "P04"]
    assert len(p04) == 1
    assert p04[0]["corrected_disposition"] == "RETAINED_PRIOR_LINEAGE_RESULT"
    assert float(p04[0]["pi"]) == pytest.approx(0.004604184795871289)


def test_nonrun_pairs_have_no_imputed_estimates(result_pairs):
    for selection_id in ("P01", "P05", "P06", "P08", "P09", "P10"):
        rows = [row for row in result_pairs if row["selection_id"] == selection_id]
        assert len(rows) == 1
        row = rows[0]
        assert row["heterozygous_snps"] == ""
        assert row["callable_bp"] == ""
        assert row["pi"] == ""
        # Pipeline defects are not biological exclusions and not zeros.
        assert row["corrected_disposition"] == "PIPELINE_DEFECT_NOT_BIOLOGICAL_EXCLUSION"
        assert "biological exclusion" in row["disposition_reason"].lower() or "not biologically excluded" in row["disposition_reason"].lower()


def test_trio_pi_values_and_report_table(result_pairs, report):
    accepted = {row["selection_id"]: float(row["pi"]) for row in result_pairs if row["row_kind"] == "accepted_bounded_core"}
    assert set(accepted) == {"P07", "P03", "P02"}
    assert accepted["P02"] > accepted["P03"] > accepted["P07"] > 0
    for value in accepted.values():
        assert repr(value) in report


# ---------------------------------------------------------------------------
# P07 reclassification
# ---------------------------------------------------------------------------


def test_p07_reclassified_to_unresolved_discordance(claims, result_pairs, report):
    claim = {row["claim_id"]: row for row in claims}["P07-DISPOSITION"]
    assert claim["corrected_classification"] == "reclassified_unresolved_discordance"
    statement = claim["corrected_statement"]
    assert "UNRESOLVED READ/ASSEMBLY DISCORDANCE" in statement
    # The old label may appear only as the reclassified prior, never as the
    # corrected disposition of the bounded result.
    p07 = [row for row in result_pairs if row["selection_id"] == "P07" and row["row_kind"] == "accepted_bounded_core"]
    assert len(p07) == 1
    assert "invalid" not in p07[0]["corrected_disposition"].lower()
    assert "excluded" not in p07[0]["corrected_disposition"].lower()
    assert "INCONCLUSIVE" in report
    assert "symmetric" in report


def test_p07_diagnostics_retained_not_transferred(claims):
    claim = {row["claim_id"]: row for row in claims}["P07-DISCORDANCE-EVIDENCE"]
    assert claim["corrected_classification"] == "retained_as_diagnostic_uncertainty"
    assert "0.436332" in claim["corrected_statement"]
    assert "50.122" in claim["corrected_statement"] and "53.166" in claim["corrected_statement"]


# ---------------------------------------------------------------------------
# Pipeline defect reclassification and reliability
# ---------------------------------------------------------------------------


REQUIRED_DEFECT_FAMILIES = {
    "fastga_execution",
    "impg_global_architecture",
    "compression_integrity",
    "sequence_dictionary_binding",
    "annotation_discovery",
    "scratch_containment",
}


def test_required_defect_families_reclassified(reliability):
    defects = [row for row in reliability if row["record_kind"] == "defect_reclassification"]
    families = {row["defect_family"] for row in defects}
    assert REQUIRED_DEFECT_FAMILIES <= families
    for row in defects:
        assert row["corrected_classification"] == "PIPELINE_DEFECT_NOT_BIOLOGICAL_EXCLUSION"
        assert row["historical_evidence"].strip()
    # P02/P03 must be cited as resolved by the bounded lineage.
    resolved = [row for row in defects if row["resolution_status"] == "resolved_in_bounded_lineage"]
    affected = {sid for row in resolved for sid in row["affected_selection_ids"].split(";")}
    assert {"P02", "P03"} <= affected


def test_review_defect_ledger_carried_through(reliability):
    review_rows = [row for row in reliability if row["record_kind"] == "review_defect_ledger"]
    ids = {row["record_id"] for row in review_rows}
    assert len(review_rows) >= 12
    for expected in ("REVIEW-I-1", "REVIEW-D-3", "REVIEW-A-1", "REVIEW-B-1"):
        assert expected in ids


def test_scheduler_metrics_recomputed(reliability):
    metrics = {row["metric"]: row["value"] for row in reliability if row["record_kind"] == "scheduler_reliability_metric"}
    prior = read_tsv(REPO_ROOT / "analysis" / "vgp_real_synthesis_v1" / "job_ledger.tsv")
    sacct = read_tsv(REPO_ROOT / "analysis" / "vgp_three_pair_bounded_sacct_v1.tsv")

    def count(rows, key, value):
        return sum(1 for row in rows if row[key] == value)

    assert int(metrics["prior_scale_packet_allocations_total"]) == len(prior)
    assert int(metrics["prior_scale_packet_completed"]) == count(prior, "state", "COMPLETED")
    assert int(metrics["prior_scale_packet_failed"]) == count(prior, "state", "FAILED")
    assert int(metrics["prior_scale_packet_cancelled"]) == count(prior, "state", "CANCELLED by 1001")
    assert int(metrics["bounded_wave_allocations_total"]) == len(sacct)
    assert int(metrics["bounded_wave_completed"]) == count(sacct, "State", "COMPLETED")
    assert int(metrics["bounded_wave_accepted_results"]) == 3
    nonterminal = (
        len(prior)
        - count(prior, "state", "COMPLETED")
        - count(prior, "state", "FAILED")
        - count(prior, "state", "CANCELLED by 1001")
    )
    assert int(metrics["prior_scale_packet_nonterminal_at_freeze"]) == nonterminal


def test_pipeline_defect_claim(claims):
    claim = {row["claim_id"]: row for row in claims}["PIPELINE-DEFECTS-NOT-EXCLUSIONS"]
    assert claim["corrected_classification"] == "supported"
    statement = claim["corrected_statement"]
    for family in ("FastGA", "IMPG", "compression", "dictionary", "annotation-discovery", "scratch"):
        assert family.lower() in statement.lower()
    assert "HARD_INVALID_PRIMARY" in statement


# ---------------------------------------------------------------------------
# Annotation catalog publication and scale semantics
# ---------------------------------------------------------------------------


def test_annotation_paths_cover_ten_pilots(annotation_paths):
    assert {row["selection_id"] for row in annotation_paths} == {
        f"P{number:02d}" for number in range(1, 11)
    }
    for row in annotation_paths:
        assert row["annotation_path"].startswith("/")
        assert row["catalog_binding_status"] == "EXACT_DICTIONARY"
        assert "NOT a biological scale-out" in row["catalog_scale_status"]


def test_only_p07_has_computed_partitions(annotation_paths):
    computed = {
        row["selection_id"]: row["annotation_partitions_computed"]
        for row in annotation_paths
    }
    assert computed["P07"] == "computed_from_bounded_callset"
    for selection_id in ("P02", "P03"):
        assert computed[selection_id] == "missing_nonblocking"
    for selection_id in ("P01", "P04", "P05", "P06", "P08", "P09", "P10"):
        assert computed[selection_id] == "not_computed_no_accepted_result"


def test_reference_bindings_are_not_claimed_exact_native_to_run_h1(annotation_paths):
    by_id = {row["selection_id"]: row for row in annotation_paths}
    # P07 is exact-native to the accepted run H1.
    assert by_id["P07"]["exact_native_to_run_h1"] == "true"
    assert by_id["P07"]["roster_native_annotation_accession"].startswith(
        by_id["P07"]["pair_h1_accession_version"] + "-"
    )
    # P02/P03 roster annotations bind to RefSeq reference assemblies, not to
    # the run H1 accessions, and no dictionary audit against the run H1 exists.
    for selection_id in ("P02", "P03"):
        row = by_id[selection_id]
        assert row["exact_native_to_run_h1"] == "false_not_audited_vs_run_h1"
        assert row["annotation_reference_accession_version"].startswith("GCF_")
        assert row["annotation_reference_accession_version"] != row["pair_h1_accession_version"]
        assert "not been performed" in row["binding_semantics"]


def test_catalog_metrics_match_catalog_summary():
    manifest = json.loads((OUT_DIR / "manifest.json").read_text(encoding="utf-8"))
    summary = json.loads(
        (REPO_ROOT / "analysis" / "vgp_annotation_catalog_summary.json").read_text(encoding="utf-8")
    )
    metrics = manifest["annotation_catalog_metrics"]
    assert metrics["freeze1_assemblies_with_preferred_exact_dictionary_annotation"] == summary["assembly_accounting"]["total"] == 581
    assert metrics["accepted_annotations"] == summary["accepted_annotations"] == 1833
    assert metrics["pilot_reference_bindings"] == summary["pilot_accounting"]["total"] == 10
    assert metrics["is_biological_scale_out"] is False
    assert metrics["biological_annotation_partitions_computed"] == 1
    assert metrics["pilot_exact_native_to_run_h1"] == 1


def test_report_states_preservation_and_links(report):
    assert "byte-identically" in report
    assert "supersession_ledger.tsv" in report
    assert "No historical file content was rewritten" in report


def test_scale_statement_explicit(report, claims):
    claim = {row["claim_id"]: row for row in claims}["SCALE-STATUS"]
    statement = claim["corrected_statement"]
    assert "What HAS been scaled" in statement
    assert "What has NOT been scaled" in statement
    assert "No full scale-out was launched" in statement or "No scale-out was launched" in statement
    assert "NOT a biological scale-out" in report or "not a biological scale-out" in report
    assert "NOT scaled" in report


def test_annotation_catalog_claim(claims):
    claim = {row["claim_id"]: row for row in claims}["ANNOT-CATALOG-SCALE"]
    assert claim["corrected_classification"] == "supported_catalog_level_only"
    assert "NOT a biological scale-out" in claim["corrected_statement"]
    assert "581" in claim["corrected_statement"]


# ---------------------------------------------------------------------------
# Supersession and historical preservation
# ---------------------------------------------------------------------------


def test_supersession_ledger_covers_prior_synthesis(supersession):
    for row in supersession:
        if row["supersession_kind"] in (
            "superseded_as_primary_evidence",
            "superseded_dispositions",
            "superseded_claims",
            "superseded_values",
            "superseded_as_core_diagnostic_retained",
            "role_clarified_reproduction_control",
        ):
            assert "byte-identical" in row["historical_preservation"], row["superseded_artifact"]
    superseded = {row["superseded_artifact"] for row in supersession}
    for path in (
        "analysis/vgp_real_synthesis_v1/report.md",
        "analysis/vgp_real_synthesis_v1/paper_pairs.tsv",
        "analysis/vgp_real_synthesis_v1/claim_ledger.tsv",
        "analysis/vgp_real_synthesis_v1/annotation_partitions.tsv",
        "analysis/vgp_real_canary_report_v1.md",
        "analysis/vgp_clean_canary_report_v1.md",
    ):
        assert path in superseded
    # The reviewed evidence itself is not superseded.
    roles = {row["superseded_artifact"]: row["supersession_kind"] for row in supersession}
    assert roles["analysis/vgp_three_pair_bounded_report_v1.md"] == "current_accepted_evidence"
    assert roles["analysis/vgp_three_pair_independent_review_v1.md"] == "current_governing_review"


def test_superseded_reports_are_preserved_not_annotated(supersession):
    """Historical artifacts stay byte-identical; supersession lives in the ledger.

    Several superseded artifacts are digest-bound by frozen evidence ledgers;
    editing them (e.g. adding banners) would silently invalidate those
    bindings, so the correction must leave them untouched.
    """
    for relative in (
        "analysis/vgp_real_synthesis_v1/report.md",
        "analysis/vgp_real_canary_report_v1.md",
        "analysis/vgp_clean_canary_report_v1.md",
    ):
        path = REPO_ROOT / relative
        assert path.is_file(), relative
        text = path.read_text(encoding="utf-8")
        assert "correct-vgp-evidence" not in text, f"historical file was annotated: {relative}"
        assert "SUPERSEDED AS" not in text, f"historical file was annotated: {relative}"


def test_preservation_digest_bindings_hold():
    # The prior synthesis manifest binds every packet output by SHA-256 and
    # the correction verified each binding at generation time.
    manifest = json.loads((REPO_ROOT / "analysis" / "vgp_real_synthesis_v1" / "manifest.json").read_text(encoding="utf-8"))
    correction = json.loads((OUT_DIR / "manifest.json").read_text(encoding="utf-8"))
    statuses = correction["preservation_statuses"]
    for relative, digest in manifest["output_digests"].items():
        assert cve.sha256_file(REPO_ROOT / relative) == digest, relative
        assert statuses.get(relative) == "digest_preserved_synthesis_manifest", relative
    for relative in (
        "analysis/vgp_real_canary_report_v1.md",
        "analysis/vgp_clean_canary_report_v1.md",
    ):
        assert statuses.get(relative) == "preserved_unmodified", relative


def test_prior_synthesis_files_still_exist_untouched():
    # The superseded TSVs are never modified by the correction: their digests
    # must match the manifest input bindings.
    manifest = json.loads((OUT_DIR / "manifest.json").read_text(encoding="utf-8"))
    for relative in (
        "analysis/vgp_real_synthesis_v1/paper_pairs.tsv",
        "analysis/vgp_real_synthesis_v1/claim_ledger.tsv",
    ):
        assert manifest["inputs"][relative] == cve.sha256_file(REPO_ROOT / relative)


def test_no_scaleout_launched(manifest=None):
    manifest = json.loads((OUT_DIR / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["no_scaleout_launched"] is True


# ---------------------------------------------------------------------------
# Determinism and manifest binding
# ---------------------------------------------------------------------------


def test_manifest_binds_all_outputs():
    manifest = json.loads((OUT_DIR / "manifest.json").read_text(encoding="utf-8"))
    expected = {
        "result_pairs.tsv",
        "annotation_paths.tsv",
        "pipeline_reliability.tsv",
        "claims_ledger.tsv",
        "supersession_ledger.tsv",
        "report.md",
    }
    # The generator binds exactly its own outputs; hand-audited companions
    # (integration_audit.md) may sit alongside without being regenerated.
    assert expected <= set(manifest["outputs"])
    for name, digest in manifest["outputs"].items():
        assert digest == cve.sha256_file(OUT_DIR / name), name


def test_regeneration_is_deterministic(tmp_path):
    exit_code = cve.main(
        [
            "--repo-root",
            str(REPO_ROOT),
            "--output-dir",
            str(tmp_path),
        ]
    )
    assert exit_code == 0
    for name in (
        "result_pairs.tsv",
        "annotation_paths.tsv",
        "pipeline_reliability.tsv",
        "claims_ledger.tsv",
        "supersession_ledger.tsv",
        "report.md",
    ):
        assert (tmp_path / name).read_bytes() == (OUT_DIR / name).read_bytes(), name


def test_broken_evidence_is_rejected(tmp_path, monkeypatch):
    # A non-passing execution record must not produce a correction.
    bad = json.loads((REPO_ROOT / "analysis" / "vgp_three_pair_execution_v2.json").read_text(encoding="utf-8"))
    bad["actual_core_biological_results"] = 2
    fake_root = tmp_path / "repo"
    (fake_root / "analysis").mkdir(parents=True)
    (fake_root / "analysis" / "vgp_three_pair_execution_v2.json").write_text(json.dumps(bad), encoding="utf-8")
    with pytest.raises(cve.CorrectionError):
        cve.load_inputs(fake_root)
