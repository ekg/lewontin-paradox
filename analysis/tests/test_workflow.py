import copy
import json
import os
import signal
import shutil
import subprocess
import sys
import time
from pathlib import Path

import pytest

from analysis.run_tier3 import (
    WorkflowError,
    audit_environment_record,
    load_workflow,
    run_array_task,
    sanitized_stage_environment,
    validate_dataset_preflight,
)
from analysis.tier3_common import sha256_file
from analysis.tier3a_vgp_compute import common_callable_concordance


FIXTURES = Path(__file__).parent / "fixtures"
ROOT = Path(__file__).parents[2]


def _annotation_provenance():
    return {
        "provider": "truth-provider",
        "release": "truth-release-1",
        "assembly_accession_version": "GCA_000001.1",
        "fasta_assembly_accession": "GCA_000001.1",
        "status": "native",
        "genetic_code": 1,
        "fasta_sha256": sha256_file(FIXTURES / "truth.fa"),
        "gff_sha256": sha256_file(FIXTURES / "truth.gff3"),
        "contig_mapping": {},
        "native_vs_projected": "native",
    }


def _stage(name, marker):
    code = (
        "from pathlib import Path; import os; "
        "p=Path(__import__('sys').argv[1]); p.write_text(os.environ.get('LEAK_ME','absent')); "
        "c=Path(__import__('sys').argv[2]); "
        "c.write_text(str(int(c.read_text())+1) if c.exists() else '1')"
    )
    return {
        "name": name,
        "argv": [sys.executable, "-c", code, "{stage_dir}/" + marker, "{dataset_dir}/counter-" + name],
        "outputs": [marker],
    }


def _direct_workflow(tmp_path):
    stages = [
        _stage("alignment", "aligned.txt"),
        _stage("mapping", "unique.txt"),
        _stage("normalized_bcf", "calls.txt"),
        _stage("annotation_4d", "fourfold.txt"),
        _stage("qc", "qc.txt"),
    ]
    return {
        "schema_version": "1.0",
        "decision_version": "tier3-decisions-v1",
        "scratch_root": str(tmp_path / "scratch"),
        "output_root": str(tmp_path / "out"),
        "datasets": [
            {
                "dataset_id": "truth.vgp.dual",
                "tier": "3a",
                "pilot": "vgp_dual_modality",
                "mode": "direct_wfmash",
                "reference_accession": "GCA_000001.1",
                "preflight": {
                    "reference_fasta": str(FIXTURES / "truth.fa"),
                    "query_fasta": str(FIXTURES / "truth.fa"),
                    "paf": str(tmp_path / "truth.paf"),
                    "callable_mask": str(FIXTURES / "expected.callable.bed"),
                    "annotation_gff": str(FIXTURES / "truth.gff3"),
                    "annotation_provenance": _annotation_provenance(),
                    "phase_identity_audit_passed": True,
                    "collapse_qc_passed": True,
                },
                "stages": stages,
            }
        ],
    }


def _write_valid_paf(path):
    sequence = "".join(
        line.strip()
        for line in (FIXTURES / "truth.fa").read_text(encoding="utf-8").splitlines()
        if not line.startswith(">")
    )
    path.write_text(
        f"chr1\t{len(sequence)}\t0\t{len(sequence)}\t+\tchr1\t{len(sequence)}\t0\t{len(sequence)}\t"
        f"{len(sequence)}\t{len(sequence)}\t60\tcg:Z:{len(sequence)}=\n",
        encoding="utf-8",
    )


def test_atomic_stage_resume_is_idempotent_and_discards_interrupted_partial(tmp_path, monkeypatch):
    workflow_value = _direct_workflow(tmp_path)
    _write_valid_paf(Path(workflow_value["datasets"][0]["preflight"]["paf"]))
    workflow_path = tmp_path / "workflow.json"
    workflow_path.write_text(json.dumps(workflow_value), encoding="utf-8")
    monkeypatch.setenv("LEAK_ME", "credential-that-must-not-propagate")

    first = run_array_task(workflow_path, tier="3a", array_index=0)
    assert first["status"] == "complete"
    dataset_dir = Path(first["dataset_directory"])
    for stage in workflow_value["datasets"][0]["stages"]:
        assert (dataset_dir / ("counter-" + stage["name"])).read_text() == "1"
    assert (dataset_dir / "stages/00-alignment/aligned.txt").read_text() == "absent"

    # A killed command can leave arbitrary partial bytes, but never a completed
    # checkpoint.  Resubmission removes that directory and leaves completed
    # stages byte-for-byte untouched.
    stale = dataset_dir / "stages/02-normalized_bcf.partial-dead-job"
    stale.mkdir(parents=True)
    (stale / "calls.txt").write_text("truncated", encoding="utf-8")
    second = run_array_task(workflow_path, tier="3a", array_index=0)
    assert second["status"] == "complete"
    assert not stale.exists()
    for stage in workflow_value["datasets"][0]["stages"]:
        assert (dataset_dir / ("counter-" + stage["name"])).read_text() == "1"
    assert all(item["resumed"] for item in second["stages"])


def test_killed_array_process_resubmits_from_last_atomic_checkpoint(tmp_path):
    value = _direct_workflow(tmp_path)
    _write_valid_paf(Path(value["datasets"][0]["preflight"]["paf"]))
    ready, release = tmp_path / "slow-stage.ready", tmp_path / "slow-stage.release"
    slow_code = (
        "from pathlib import Path; import sys,time; "
        "out,ready,release=map(Path,sys.argv[1:]); ready.write_text('ready'); "
        "time.sleep(30) if not release.exists() else None; out.write_text('complete')"
    )
    value["datasets"][0]["stages"][0] = {
        "name": "alignment",
        "argv": [sys.executable, "-c", slow_code, "{stage_dir}/aligned.txt", str(ready), str(release)],
        "outputs": ["aligned.txt"],
    }
    workflow = tmp_path / "workflow.kill.json"
    workflow.write_text(json.dumps(value), encoding="utf-8")
    process = subprocess.Popen(
        [sys.executable, "-m", "analysis.run_tier3", "run-array", str(workflow), "--tier", "3a", "--index", "0"],
        cwd=ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        start_new_session=True,
    )
    for _ in range(200):
        if ready.is_file():
            break
        time.sleep(0.01)
    assert ready.is_file(), "slow stage never started"
    os.killpg(process.pid, signal.SIGKILL)
    process.wait(timeout=5)
    assert process.returncode < 0
    dataset_dir = Path(value["scratch_root"]) / "3a/truth.vgp.dual"
    assert list((dataset_dir / "stages").glob("00-alignment.partial-*"))
    assert not (dataset_dir / "stages/00-alignment/.complete.json").exists()

    release.write_text("resume", encoding="utf-8")
    resumed = run_array_task(workflow, tier="3a", array_index=0)
    assert resumed["status"] == "complete"
    assert not list((dataset_dir / "stages").glob("*.partial-*"))
    assert (dataset_dir / "stages/00-alignment/aligned.txt").read_text() == "complete"
    rerun = run_array_task(workflow, tier="3a", array_index=0)
    assert all(stage["resumed"] for stage in rerun["stages"])


def test_preflight_rejects_missing_mask_reference_mismatch_and_bad_or_multiple_cigar(tmp_path):
    workflow_value = _direct_workflow(tmp_path)
    paf = Path(workflow_value["datasets"][0]["preflight"]["paf"])
    _write_valid_paf(paf)
    dataset = workflow_value["datasets"][0]
    validate_dataset_preflight(dataset)

    no_mask = copy.deepcopy(dataset)
    no_mask["mode"] = "deposited_vcf"
    del no_mask["preflight"]["callable_mask"]
    with pytest.raises(WorkflowError, match="callable mask"):
        validate_dataset_preflight(no_mask)

    wrong_accession = copy.deepcopy(dataset)
    wrong_accession["preflight"]["annotation_provenance"]["assembly_accession_version"] = "GCA_999999.1"
    with pytest.raises(WorkflowError, match="assembly"):
        validate_dataset_preflight(wrong_accession)

    paf.write_text(paf.read_text().replace("cg:Z:35=", "cg:Z:35M"), encoding="utf-8")
    with pytest.raises(WorkflowError, match="extended"):
        validate_dataset_preflight(dataset)

    _write_valid_paf(paf)
    paf.write_text(paf.read_text() * 2, encoding="utf-8")
    with pytest.raises(WorkflowError, match="multiple|unique"):
        validate_dataset_preflight(dataset)


def test_native_annotation_provenance_is_a_hard_primary_gate(tmp_path):
    workflow_value = _direct_workflow(tmp_path)
    _write_valid_paf(Path(workflow_value["datasets"][0]["preflight"]["paf"]))
    dataset = workflow_value["datasets"][0]
    for mutation, error in [
        (("status", "projected"), "native"),
        (("gff_sha256", "0" * 64), "GFF checksum"),
        (("fasta_sha256", "0" * 64), "FASTA checksum"),
    ]:
        invalid = copy.deepcopy(dataset)
        invalid["preflight"]["annotation_provenance"][mutation[0]] = mutation[1]
        with pytest.raises(WorkflowError, match=error):
            validate_dataset_preflight(invalid)

    unavailable = copy.deepcopy(dataset)
    unavailable["annotation_derived"] = "unavailable_missing_native"
    del unavailable["preflight"]["annotation_gff"]
    del unavailable["preflight"]["annotation_provenance"]
    audit = validate_dataset_preflight(unavailable)
    assert audit["annotation"] == {
        "status": "unavailable_missing_native_exact_reference_annotation",
        "gc3": "unavailable",
        "fourfold": "unavailable",
        "whole_genome_gc": "eligible",
    }


def test_dual_modality_vgp_pilot_passes_on_common_callable_truth():
    deposited = {("chr1", 2), ("chr1", 19)}
    direct = {("chr1", 2), ("chr1", 19)}
    result = common_callable_concordance(
        reference_fasta=FIXTURES / "truth.fa",
        left_callable_bed=FIXTURES / "expected.callable.bed",
        right_callable_bed=FIXTURES / "expected.callable.bed",
        left_heterozygous=deposited,
        right_heterozygous=direct,
        synthetic_fixture=True,
    )
    assert result["passed"] is True
    assert result["common_callable_bases"] == 12
    assert result["snv_precision"] == result["snv_recall"] == 1.0


def test_only_declared_scheduler_tier3_and_scratch_variables_survive(monkeypatch, tmp_path):
    monkeypatch.setenv("SLURM_JOB_ID", "123")
    monkeypatch.setenv("TIER3_DATA_ROOT", "/data")
    monkeypatch.setenv("SCRATCH", "/scratch")
    monkeypatch.setenv("TMPDIR", "/caller-tmp")
    monkeypatch.setenv("AWS_SECRET_ACCESS_KEY", "secret")
    monkeypatch.setenv("CONDA_PREFIX", "/forbidden")
    monkeypatch.setenv("LEAK_ME", "secret")
    monkeypatch.setenv("SSL_CERT_DIR", "/gnu/store/nss-certs/etc/ssl/certs")
    monkeypatch.setenv("SSL_CERT_FILE", "/gnu/store/nss-certs/etc/ssl/certs/ca-certificates.crt")
    result = sanitized_stage_environment(tmp_path)
    assert result["SLURM_JOB_ID"] == "123"
    assert result["TIER3_DATA_ROOT"] == "/data"
    assert result["SCRATCH"] == "/scratch"
    assert result["TMPDIR"] == str(tmp_path)
    assert result["SSL_CERT_DIR"].startswith("/gnu/store/")
    assert result["SSL_CERT_FILE"].startswith("/gnu/store/")
    assert "AWS_SECRET_ACCESS_KEY" not in result
    assert "CONDA_PREFIX" not in result
    assert "LEAK_ME" not in result


def test_workflow_rejects_shell_commands_forbidden_managers_and_wrong_stage_order(tmp_path):
    value = _direct_workflow(tmp_path)
    _write_valid_paf(Path(value["datasets"][0]["preflight"]["paf"]))
    path = tmp_path / "workflow.json"
    for argv, error in [
        (["bash", "-c", "echo unsafe"], "shell"),
        (["conda", "run", "python3"], "Conda|forbidden"),
    ]:
        invalid = copy.deepcopy(value)
        invalid["datasets"][0]["stages"][0]["argv"] = argv
        path.write_text(json.dumps(invalid), encoding="utf-8")
        with pytest.raises(WorkflowError, match=error):
            load_workflow(path)
    invalid = copy.deepcopy(value)
    invalid["datasets"][0]["stages"][0:2] = reversed(invalid["datasets"][0]["stages"][0:2])
    path.write_text(json.dumps(invalid), encoding="utf-8")
    with pytest.raises(WorkflowError, match="stage order"):
        load_workflow(path)


def test_environment_record_requires_resolved_store_paths_derivations_and_versions(tmp_path):
    profile = tmp_path / "profile"
    profile.symlink_to("/gnu/store/00000000000000000000000000000000-tier3-profile")
    record = {
        "channel_commit": "44bbfc24e4bcc48d0e3343cd3d83452721af8c36",
        "profile_store_path": str(profile.resolve(strict=False)),
        "profile_gc_root": str(profile.absolute()),
        "derivations": ["/gnu/store/" + "1" * 32 + "-tier3.drv"],
        "store_paths": ["/gnu/store/" + "2" * 32 + "-bcftools-1.14"],
        "tool_versions": {
            "python3": "Python 3",
            "pytest": "pytest",
            "samtools": "samtools",
            "bcftools": "bcftools 1.14",
            "bgzip": "bgzip",
            "tabix": "tabix",
            "bedtools": "bedtools",
            "vcftools": "vcftools",
            "wfmash": "wfmash",
        },
        "resolved_channels_scm": "channel 44bbfc24e4bcc48d0e3343cd3d83452721af8c36\n",
        "pack_fallback": {"required": False},
    }
    record["resolved_channels_sha256"] = __import__("hashlib").sha256(
        record["resolved_channels_scm"].encode()
    ).hexdigest()
    assert audit_environment_record(record, require_existing_profile=False) == record
    for key in ("derivations", "store_paths", "tool_versions"):
        invalid = copy.deepcopy(record)
        invalid[key] = [] if key != "tool_versions" else {}
        with pytest.raises(WorkflowError, match=key.replace("_", " ")):
            audit_environment_record(invalid, require_existing_profile=False)


def test_all_slurm_jobs_use_pinned_pure_guix_and_pass_sbatch_test_only():
    slurm = ROOT / "analysis" / "slurm"
    jobs = sorted(slurm.glob("*_array.sh")) + [slurm / "compute_smoke.sh"]
    assert jobs
    sbatch = shutil.which("sbatch") or "/usr/local/bin/sbatch"
    if not Path(sbatch).is_file():
        pytest.skip("Slurm client is external to the Guix scientific closure")
    forbidden = ("conda", "micromamba", "source activate")
    for job in jobs:
        text = job.read_text(encoding="utf-8").lower()
        assert "#sbatch --export=none" in text
        assert "guix_job.sh" in text
        assert not any(word in text for word in forbidden)
        completed = subprocess.run(
                [sbatch, "--test-only", str(job), "/tmp/workflow.json", str(ROOT)],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
        )
        assert completed.returncode == 0, completed.stdout


def test_pilot_registry_names_all_frozen_pilots_without_promoting_pending_diversity():
    registry = json.loads((ROOT / "analysis/pilot_results/pilot_registry.json").read_text())
    assert {item["pilot_id"] for item in registry["pilots"]} == {
        "composition_drosophila_melanogaster",
        "composition_homo_sapiens",
        "population_dgrp_freeze2",
        "population_ag1000g_phase3",
        "vgp_dual_modality_individual",
    }
    candidates = {item["pilot_id"]: item for item in registry["pilots"]}
    assert candidates["population_dgrp_freeze2"]["eligibility"].startswith("ineligible_pending")
    assert candidates["population_ag1000g_phase3"]["eligibility"].startswith("ineligible_pending")
    assert candidates["vgp_dual_modality_individual"]["required_modalities"] == [
        "deposited_exact_reference_variants_plus_mask",
        "direct_wfmash_extended_cigar",
    ]
