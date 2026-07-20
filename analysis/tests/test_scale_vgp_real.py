import csv
import json
from pathlib import Path

import pytest

from analysis import scale_vgp_real as scale


def read_rows(path: Path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_closed_world_accounts_all_catalog_links_and_only_audited_roster():
    paths = scale.canonical_paths(Path("/moosefs/erikg/vgp"))
    catalog = scale.load_catalog(paths["catalog"])
    roster = scale.load_roster(scale.ANALYSIS / "vgp_10_pair_manifest.tsv")
    catalog_rows, pair_rows, matched = scale.build_closed_world(catalog, roster, paths["root"])
    assert len(catalog_rows) == 716
    assert len(pair_rows) == 569
    assert set(matched) == set(scale.PAIR_IDS)
    eligible = [row for row in pair_rows if row["eligibility"] == "eligible_exact_audited"]
    assert len(eligible) == 10
    assert {row["selection_id"] for row in eligible} == set(scale.PAIR_IDS)
    assert sum(row["eligibility"] != "eligible_exact_audited" for row in pair_rows) == 559
    assert all(row["population_inference_authorized"] == "false" for row in pair_rows)


def test_roster_rejects_missing_same_individual_evidence(tmp_path: Path):
    rows = scale.load_roster(scale.ANALYSIS / "vgp_10_pair_manifest.tsv")
    rows[0]["reciprocal_pair_evidence"] = "labels only"
    path = tmp_path / "roster.tsv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    with pytest.raises(scale.ScaleError, match="same-individual"):
        scale.load_roster(path)


def test_zero_job_completion_is_forbidden(tmp_path: Path):
    path = tmp_path / "empty.tsv"
    path.write_text("selection_id\tjob_id\tcanonical_vgp_root\n", encoding="utf-8")
    with pytest.raises(scale.ScaleError, match="zero biological"):
        scale.load_submissions(path)


def test_submission_coverage_can_be_union_of_pilot_and_scale_wave():
    rows = [
        {"selection_id": pair, "canonical_vgp_root": "/moosefs/erikg/vgp"}
        for pair in scale.PAIR_IDS
    ]
    scale.require_submission_coverage(rows[:9] + rows[9:])
    with pytest.raises(scale.ScaleError, match="coverage is incomplete"):
        scale.require_submission_coverage(rows[:9])


def test_latest_attempts_preserves_real_stage_state():
    submissions = [
        {"selection_id": pair, "job_id": str(100 + index), "stage": "mapping"}
        for index, pair in enumerate(scale.PAIR_IDS)
    ]
    sacct = [{"JobIDRaw": row["job_id"], "State": "COMPLETED"} for row in submissions]
    result = scale.latest_attempts(submissions, sacct)
    assert set(result) == set(scale.PAIR_IDS)
    assert all(value["attempted_slurm_jobs"] == 1 for value in result.values())
    assert all(value["latest_job_state"] == "COMPLETED" for value in result.values())


def test_latest_attempts_does_not_hide_failure_behind_dependency_placeholder():
    submissions = []
    sacct = []
    for index, pair in enumerate(scale.PAIR_IDS):
        mapping = str(200 + index * 2)
        downstream = str(201 + index * 2)
        submissions.extend((
            {"selection_id": pair, "job_id": mapping, "stage": "mapping"},
            {"selection_id": pair, "job_id": downstream, "stage": "impg"},
        ))
        sacct.extend((
            {"JobIDRaw": mapping, "JobName": f"vgp10-{pair}-mapping", "State": "FAILED"},
            {"JobIDRaw": downstream, "JobName": f"vgp10-{pair}-impg", "State": "PENDING"},
        ))
    result = scale.latest_attempts(submissions, sacct)
    assert all(value["latest_stage"] == "mapping" for value in result.values())
    assert all(value["latest_job_state"] == "FAILED" for value in result.values())


def test_latest_attempts_prefers_later_successful_retry():
    submissions = []
    sacct = []
    for index, pair in enumerate(scale.PAIR_IDS):
        failed = str(400 + index * 2)
        passed = str(401 + index * 2)
        submissions.extend((
            {"selection_id": pair, "job_id": failed, "stage": "verification"},
            {"selection_id": pair, "job_id": passed, "stage": "verification"},
        ))
        sacct.extend((
            {"JobIDRaw": failed, "JobName": f"vgp-scale-{pair}-verify", "State": "FAILED"},
            {"JobIDRaw": passed, "JobName": f"vgp-scale-{pair}-verify", "State": "COMPLETED"},
        ))
    result = scale.latest_attempts(submissions, sacct)
    assert all(value["latest_job_state"] == "COMPLETED" for value in result.values())


def test_worker_enforces_private_node_local_scratch_contract():
    text = (scale.ANALYSIS / "slurm/scale_vgp_real/verify_pair.sh").read_text()
    for token in ("mktemp -d \"/scratch/", "TMPDIR=$scratch", "TMP=$scratch", "TEMP=$scratch", "cd \"$scratch\"", "/proc/$$/cwd", "scratch_resolved"):
        assert token in text
    assert "canonical_vgp_root" in text
    assert "verify-capture" in text
    assert 'profile.startswith("/gnu/store/")' in text


def test_fastga_revalidation_worker_has_live_fail_closed_scratch_contract():
    worker = (scale.ANALYSIS / "slurm/scale_vgp_real/revalidate_p07_fastga.sh").read_text()
    for token in (
        'mktemp -d "/scratch/vgp-scale-fastga-P07-',
        "TMPDIR=$scratch TMP=$scratch TEMP=$scratch",
        'cd -- "$scratch"',
        "/proc/$$/cwd",
        "fastga_scratch_guard.py\" check",
        "pkill -TERM",
        "exit 70",
        "fastga_scratch_guard.py\" finalize",
        "--num-mappings 1:1",
        "Promotion starts only after alignment success",
        "node_local_base_resolved=$(readlink -f -- /scratch)",
        'scratch_resolved != "$node_local_base_resolved"/vgp-scale-fastga-P07-',
        "refusing cleanup outside validated requested/resolved scratch roots",
    ):
        assert token in worker
    assert "/tmp/" not in worker
    assert 'cp -a -- "$partial" "$promote"' in worker
    assert '--scratch "$scratch"' in worker
    assert '--scratch "$scratch_resolved"' not in worker


def test_fastga_input_resolver_is_canonical_and_integrity_checked():
    resolver = (scale.ANALYSIS / "slurm/scale_vgp_real/resolve_p07_inputs.py").read_text()
    assert '"/moosefs/erikg/vgp"' in resolver
    assert "compressed_sha256" in resolver
    assert "FastGA amendment mismatch" in resolver
    assert "lewontin-paradox-data" not in resolver


def test_scale_submitter_materializes_canonical_binding_and_submits_fastga():
    submitter = (scale.ANALYSIS / "run_scale_vgp_real_slurm.sh").read_text()
    assert "vgp-real-scaleout-input-binding-v1" in submitter
    assert '"canonical_vgp_root":str(root)' in submitter
    assert '"migration_action":"verified hard-link reuse; no redownload"' in submitter
    assert "fastga_scratch_revalidation" in submitter
    assert "--cpus-per-task=32 --mem=128G" in submitter


def test_summary_preserves_same_pair_dependence_and_scenario_uncertainty():
    text = (scale.ANALYSIS / "scale_vgp_real.py").read_text()
    assert "same_pair_pi_psmc_non_independent" in text
    assert "nine generic mutation-rate x generation-time scenarios; none preferred" in text
    assert '"population_inference_authorized": False' in text


def test_repaired_psmc_packet_is_centered_and_rooted_canonically():
    packet = json.loads((scale.ANALYSIS / "vgp_psmc_bootstrap_repair_v1.json").read_text())
    assert packet["canonical_vgp_root"] == "/moosefs/erikg/vgp"
    assert packet["passed"] is True
    for pair in ("P04", "P07"):
        row = packet["pairs"][pair]
        assert row["finite_bootstraps"] == row["bootstrap_attempts"] == 200
        assert row["masked_and_callable_sampling_population_preserved"] is True
        assert row["centering_diagnostic"]["passed"] is True


def test_committed_scaleout_is_closed_world_nonzero_and_biological():
    output = scale.ANALYSIS / "vgp_real_scaleout_v1"
    summary = json.loads((output / "summary.json").read_text())
    assert summary["canonical_vgp_root"] == "/moosefs/erikg/vgp"
    assert summary["canonical_manifest_root"] == "/moosefs/erikg/vgp/derived/scale-vgp-real-v1/manifests"
    assert summary["counts"]["catalog_rows"] == 716
    assert summary["counts"]["catalog_links"] == 569
    assert summary["counts"]["audited_eligible_pairs"] == 10
    assert summary["counts"]["submitted_biological_jobs"] > 0
    assert summary["counts"]["scale_verification_jobs"] == 6
    assert summary["counts"]["fastga_scratch_revalidation_jobs"] == 4
    assert summary["counts"]["scaleout_jobs"] == 10
    assert summary["counts"]["fastga_scratch_contracts_passed"] == 2
    assert summary["counts"]["verified_pairs"] == 2
    assert summary["counts"]["scenario_rows"] == 1152
    assert summary["counts"]["psmc_trajectory_rows"] == 128
    assert summary["counts"]["annotation_partition_rows"] == 6
    assert summary["telemetry_calibrated_per_pair_resource_plan"]["pair_count"] == 10
    assert summary["telemetry_calibrated_per_pair_resource_plan"]["measured_canary_job_id"] == "1781798"
    assert {row["selection_id"] for row in summary["verified_pair_estimates"]} == {"P04", "P07"}
    assert all(row["same_pair_pi_psmc_non_independent"] for row in summary["verified_pair_estimates"])


def test_committed_pair_ledger_accounts_every_link_and_hard_failure_class():
    rows = read_rows(scale.ANALYSIS / "vgp_real_scaleout_v1/pair_accounting.tsv")
    assert len(rows) == 569
    assert all(row["canonical_vgp_root"] == "/moosefs/erikg/vgp" for row in rows)
    eligible = {row["selection_id"]: row for row in rows if row["selection_id"]}
    assert set(eligible) == set(scale.PAIR_IDS)
    assert {eligible[pair]["execution_disposition"] for pair in ("P04", "P07")} == {"VERIFIED_CORE_COMPLETE"}
    assert all(
        eligible[pair]["fastga_scratch_contract_status"] == "PASS_LIVE_PROC_NODE_LOCAL_SCRATCH"
        for pair in ("P04", "P07")
    )
    assert all(eligible[pair]["hard_failure_class"] for pair in ("P01", "P02", "P03", "P05"))
    assert all(
        eligible[pair]["execution_disposition"] == "HARD_EXECUTION_ERROR_NO_ESTIMATE"
        for pair in ("P06", "P09", "P10")
    )
    assert eligible["P08"]["execution_disposition"] == "RUNNING_RESUMABLE_WAVE"
    assert all(row["population_inference_authorized"] == "false" for row in eligible.values())


def test_all_trajectory_scenario_annotation_rows_name_canonical_root():
    output = scale.ANALYSIS / "vgp_real_scaleout_v1"
    scenarios = read_rows(output / "scenario_uncertainty.tsv")
    trajectories = read_rows(output / "psmc_unscaled_trajectories.tsv")
    annotation = read_rows(output / "exact_native_annotation_partitions.tsv")
    assert len(scenarios) == 1152 and {row["selection_id"] for row in scenarios} == {"P04", "P07"}
    assert len(trajectories) == 128 and {row["selection_id"] for row in trajectories} == {"P04", "P07"}
    assert len(annotation) == 6 and {row["selection_id"] for row in annotation} == {"P07"}
    assert {row["scenario_id"] for row in scenarios if row["selection_id"] == "P04"} == {
        "SENS_MU1.0E-08_G1Y", "SENS_MU1.0E-08_G2Y", "SENS_MU1.0E-08_G4Y",
        "SENS_MU2.0E-08_G1Y", "SENS_MU2.0E-08_G2Y", "SENS_MU2.0E-08_G4Y",
        "SENS_MU5.0E-09_G1Y", "SENS_MU5.0E-09_G2Y", "SENS_MU5.0E-09_G4Y",
    }
    assert all(
        row["canonical_vgp_root"] == "/moosefs/erikg/vgp"
        for row in [*scenarios, *trajectories, *annotation]
    )


def test_fastga_contract_ledger_has_both_completed_pairs_and_private_scratch():
    rows = read_rows(scale.ANALYSIS / "vgp_real_scaleout_v1/fastga_scratch_contracts.tsv")
    assert {row["selection_id"] for row in rows} == {"P04", "P07"}
    assert all(row["canonical_vgp_root"] == "/moosefs/erikg/vgp" for row in rows)
    assert all(row["status"] == "PASS_LIVE_PROC_NODE_LOCAL_SCRATCH" for row in rows)
    assert all(row["requested_scratch_root"].startswith("/scratch/") for row in rows)
    assert all(row["resolved_node_local_scratch_root"] for row in rows)
    assert all(row["frozen_mapping_exact_match"] == "true" for row in rows)


def test_emitted_resource_plan_is_per_pair_and_telemetry_calibrated():
    plan = json.loads((scale.ANALYSIS / "vgp_real_scaleout_v1/per_pair_resource_plan.json").read_text())
    assert plan["canonical_vgp_root"] == "/moosefs/erikg/vgp"
    assert set(plan["pairs"]) == set(scale.PAIR_IDS)
    assert plan["basis"]["slurm_job_id"] == "1781798"
    assert "Per-genome allocations scaled" in plan["basis"]["policy"]
