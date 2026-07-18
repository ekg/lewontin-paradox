import csv
import json
from collections import Counter
from pathlib import Path

import pytest

from analysis import vgp_core_scaleout as scale


ROOT = Path(__file__).parents[2]


def rows(path):
    with Path(path).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


@pytest.fixture(scope="module")
def generated(tmp_path_factory):
    output = tmp_path_factory.mktemp("vgp-core-scaleout")
    summary = scale.generate(output)
    return output, summary


def test_frozen_catalog_discovers_every_linked_haplotype_once():
    catalog = scale.verify_catalog(scale.CATALOG)
    links = [link for row in catalog for link in scale.linked_records(row)]
    assert len(catalog) == 716
    assert len(links) == 569
    assert Counter(link[0] for link in links) == {"other_high_quality": 264, "alternate": 305}
    # Five catalog rows have unequal assembly-ID/accession list cardinality;
    # together they contain eleven accession entries whose assembly ID is not
    # safely assignable by position.
    assert sum(not link[4] for link in links) == 11
    assert sum(
        row["Accession # for main haplotype"].strip() == link[3]
        for row in catalog for link in scale.linked_records(row)
    ) == 3


def test_closed_world_manifest_has_716_rows_plus_569_links(generated):
    output, summary = generated
    manifest = rows(output / "vgp_core_scaleout_manifest.tsv")
    catalog_rows = [row for row in manifest if row["record_type"] == "catalog_row"]
    pairs = [row for row in manifest if row["record_type"] == "linked_pair"]
    assert len(manifest) == 1_285
    assert len(catalog_rows) == 716
    assert len({row["catalog_row"] for row in catalog_rows}) == 716
    assert len(pairs) == 569
    assert len({row["pair_id"] for row in pairs}) == 569
    assert summary["catalog_row_dispositions"] == {"excluded": 187, "failed": 529}
    assert summary["pair_dispositions"] == {"excluded": 3, "failed": 566}
    assert summary["distinct_nonself_pair_candidates"] == 566


def test_pair_qc_is_one_to_one_and_no_pair_is_silently_eligible(generated):
    output, _ = generated
    manifest = rows(output / "vgp_core_scaleout_manifest.tsv")
    pairs = [row for row in manifest if row["record_type"] == "linked_pair"]
    qc = rows(output / "vgp_core_scaleout_qc.tsv")
    assert {row["pair_id"] for row in pairs} == {row["pair_id"] for row in qc}
    assert len(qc) == 569
    assert {row["core_eligible"] for row in qc} == {"false"}
    assert {row["psmc_bootstrap_attempts"] for row in qc} == {"0"}
    assert {row["annotation_absence_core_veto"] for row in qc} == {"false"}
    assert {row["same_pair_psmc_independent_evidence"] for row in qc} == {"false"}
    assert all(";" not in row["primary_reason_code"] for row in qc)


def test_zero_job_boundary_and_non_authorizing_wave_templates(generated):
    output, summary = generated
    telemetry = rows(output / "vgp_core_scaleout_telemetry.tsv")
    waves = rows(output / "vgp_core_scaleout_wave_manifest.tsv")
    assert summary["eligible_pairs"] == summary["completed_pairs"] == 0
    assert summary["slurm_jobs_submitted"] == summary["atomic_promotions"] == 0
    assert sum(int(row["slurm_jobs_submitted"]) for row in telemetry) == 0
    gate = next(row for row in waves if row["wave_id"] == "GATE-0")
    assert gate["authorization"] == "NOT_AUTHORIZED"
    assert gate["max_pair_jobs"] == gate["max_array_concurrency"] == "0"
    templates = [row for row in waves if row["record_type"] == "planning_template"]
    assert {row["scenario"] for row in templates} == {"low", "base", "high"}
    assert {row["authorization"] for row in templates} == {"NOT_AUTHORIZED_PLANNING_ONLY"}
    assert all(int(row["max_pair_jobs"]) == 25 for row in templates)
    assert all(int(row["max_transient_retries"]) == 2 for row in templates)
    assert all(int(row["max_scientific_relaxation_retries"]) == 0 for row in templates)


def test_paper_tables_report_nonmaterialization_and_nonindependence(generated):
    output, _ = generated
    pairs = rows(output / "vgp_core_scaleout_paper_pairs.tsv")
    sensitivity = rows(output / "vgp_core_scaleout_sensitivity.tsv")
    scenarios = rows(output / "vgp_core_scaleout_scaling_scenarios.tsv")
    independent = rows(output / "vgp_core_scaleout_independent_validation.tsv")
    text = (output / "vgp_core_scaleout_results.md").read_text()
    assert len(pairs) == 569
    assert {row["core_callable_diversity"] for row in pairs} == {"NOT_ESTIMABLE"}
    assert {row["same_pair_psmc_independent_evidence"] for row in pairs} == {"false"}
    assert {row["status"] for row in sensitivity} == {"NOT_ESTIMABLE_ZERO_APPROVED_COMPLETED_PAIRS"}
    assert scenarios[0]["scenario_id"] == "UNSCALED_PRIMARY"
    assert len(independent) == 5
    assert {row["selected_eligible_pairs"] for row in independent} == {"0"}
    assert {row["completed_checks"] for row in independent} == {"0"}
    assert {row["slurm_jobs_submitted"] for row in independent} == {"0"}
    assert "is not statistically independent evidence" in text
    assert "Technical non-execution is not a low-diversity result" in " ".join(text.split())
    assert "No valid non-annotated core result was deleted" in text


def test_summary_binds_atomic_output_digests(generated):
    output, summary = generated
    committed = json.loads((output / "vgp_core_scaleout_summary.json").read_text())
    assert committed == summary
    for name, digest in summary["output_digests"].items():
        filename = {
            "manifest": "vgp_core_scaleout_manifest.tsv",
            "qc": "vgp_core_scaleout_qc.tsv",
            "telemetry": "vgp_core_scaleout_telemetry.tsv",
            "waves": "vgp_core_scaleout_wave_manifest.tsv",
            "paper_pairs": "vgp_core_scaleout_paper_pairs.tsv",
            "paper_summary": "vgp_core_scaleout_paper_summary.tsv",
            "sensitivity": "vgp_core_scaleout_sensitivity.tsv",
            "scaling": "vgp_core_scaleout_scaling_scenarios.tsv",
            "results": "vgp_core_scaleout_results.md",
            "independent_validation": "vgp_core_scaleout_independent_validation.tsv",
        }[name]
        assert scale.sha256_file(output / filename) == digest
    assert not list(output.glob(".*.partial-*"))


def test_authorization_drift_refuses_instead_of_launching(tmp_path):
    decision = json.loads(scale.DECISION.read_text())
    decision["full_biological_scaleout_authorized"] = True
    changed = tmp_path / "decision.json"
    changed.write_text(json.dumps(decision))
    with pytest.raises(scale.ScaleoutError, match="authorization boundary drift"):
        scale.verify_authorization(changed, scale.REVIEW_SCALEOUT)


def test_mirror_launch_or_state_drift_requires_new_review(tmp_path):
    summary = json.loads(scale.MIRROR_SUMMARY.read_text())
    summary["bulk_launch"]["launched"] = True
    changed = tmp_path / "mirror.json"
    changed.write_text(json.dumps(summary))
    with pytest.raises(scale.ScaleoutError, match="mirror launch state drift"):
        scale.verify_mirror(changed, scale.MIRROR_MANIFEST)


def test_committed_artifacts_reconcile_and_use_exact_frozen_release():
    summary = json.loads((ROOT / "analysis/vgp_core_scaleout_summary.json").read_text())
    assert summary["catalog"] == {
        "commit": scale.CATALOG_COMMIT, "sha256": scale.CATALOG_SHA256, "rows": 716,
    }
    assert summary["linked_haplotype_entries"] == 569
    assert summary["full_biological_scaleout_authorized"] is False
    assert summary["mirror_verified_or_reused_objects"] == 0
    assert summary["biological_estimates"] == 0
