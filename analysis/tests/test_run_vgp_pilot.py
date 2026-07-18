import csv
import hashlib
import json
import stat
from copy import deepcopy
from pathlib import Path

import pytest

from analysis import acquire_vgp_pilot as acquisition
from analysis import run_vgp_pilot as runner
from analysis.tier3_common import Tier3ValidationError


ROOT = Path(__file__).parents[2]


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _outputs(tmp_path: Path) -> dict:
    return {
        "output_run_manifest_path": tmp_path / "run_manifest.tsv",
        "output_slurm_telemetry_path": tmp_path / "slurm.tsv",
        "output_results_path": tmp_path / "results.tsv",
        "output_exclusions_path": tmp_path / "exclusions.tsv",
        "output_refusals_path": tmp_path / "refusals.tsv",
        "output_report_path": tmp_path / "report.md",
    }


def _run(tmp_path: Path, **overrides):
    kwargs = {
        "gate_path": ROOT / "analysis/vgp_pilot_gate.json",
        "manifest_path": ROOT / "analysis/vgp_pilot_manifest.tsv",
        "root_config_path": ROOT / "analysis/vgp_data_root_config.json",
        "acquisition_manifest_path": ROOT / "analysis/vgp_pilot_acquisition_manifest.tsv",
        "inventory_path": ROOT / "analysis/vgp_pilot_immutable_object_inventory.tsv",
        "sweepga_build_path": ROOT / "analysis/sweepga_origin_main_build.json",
        "impg_handoff_path": ROOT / "analysis/sweepga_impg_observed.json",
        "worker_path": ROOT / "analysis/slurm/run_repaired_vgp.sh",
        **_outputs(tmp_path),
    }
    kwargs.update(overrides)
    return runner.run(**kwargs)


def _gate() -> dict:
    return json.loads((ROOT / "analysis/vgp_pilot_gate.json").read_text(encoding="utf-8"))


def _write_gate(tmp_path: Path, payload: dict, name: str = "gate.json") -> Path:
    path = tmp_path / name
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path


def test_current_repaired_nogo_is_not_submitted_with_exact_zero_use(tmp_path):
    calls = []
    result = _run(tmp_path, submitter=lambda argv: calls.append(argv) or "999")
    assert result["status"] == "refused_preflight"
    assert result["final_state"] == "NOT_SUBMITTED"
    assert result["failure_code"] == "GATE_NO_GO"
    assert result["slurm_jobs_submitted"] == result["compute_jobs_started"] == 0
    assert result["core_seconds"] == result["network_bytes"] == 0
    assert calls == []

    manifest = _read_tsv(tmp_path / "run_manifest.tsv")
    summary = manifest[0]
    assert summary["gate_decision"] == "NO_GO"
    assert summary["final_state"] == "NOT_SUBMITTED"
    for field in (
        "gate_file_sha256", "decision_sha256", "authorization_tuple_digest", "manifest_sha256",
        "acquisition_manifest_sha256", "immutable_inventory_sha256", "root_config_sha256",
        "input_bundle_digest", "root_contract_digest", "cap_vector_digest", "retrieval_digest",
        "pair_evidence_digest", "measurement_contract_digest", "environment_digest",
        "sweepga_build_sha256", "impg_handoff_sha256", "worker_sha256",
    ):
        assert summary[field], field
    assert {row["failure_code"] for row in manifest[1:]} == {
        "CAP_MOOSEFS_READ_GB_EXCEEDED", "CAP_SCRATCH_GIB_EXCEEDED", "QUOTA_UNAVAILABLE"
    }

    telemetry = _read_tsv(tmp_path / "slurm.tsv")
    assert telemetry[0]["final_state"] == "NOT_SUBMITTED"
    for field in ("elapsed_seconds", "cpu_time_seconds", "scratch_peak_gb", "io_read_gb", "io_write_gb", "metadata_operations", "network_bytes"):
        assert telemetry[0][field] == "0"
    assert telemetry[0]["sbatch_command"] == telemetry[0]["slurm_job_id"] == ""

    refusal = _read_tsv(tmp_path / "refusals.tsv")[0]
    for field in ("sbatch_commands_issued", "slurm_jobs_submitted", "compute_jobs_started", "core_seconds", "scratch_bytes", "io_read_bytes", "io_write_bytes", "network_bytes", "provider_requests", "demographic_inferences"):
        assert refusal[field] == "0"
    assert refusal["full_catalog_downloaded"] == refusal["population_bulk_downloaded"] == "false"
    assert len(refusal["evidence_sha256"]) == 64
    assert "No callable/queryable denominator or biological result was imputed" in (tmp_path / "report.md").read_text()


def test_unknown_and_missing_gate_refuse_before_submitter(tmp_path):
    unknown = _gate()
    unknown["decision"]["status"] = "MAYBE"
    calls = []
    result = _run(tmp_path / "unknown", gate_path=_write_gate(tmp_path, unknown), submitter=lambda argv: calls.append(argv) or "1")
    assert result["failure_code"] == "GATE_DECISION_UNKNOWN"
    assert result["final_state"] == "NOT_SUBMITTED"
    assert calls == []

    result = _run(tmp_path / "missing", gate_path=tmp_path / "absent.json", submitter=lambda argv: calls.append(argv) or "2")
    assert result["failure_code"] == "BOUND_INPUT_MISSING"
    assert result["final_state"] == "NOT_SUBMITTED"
    assert calls == []


@pytest.mark.parametrize("digest_key", [
    "authorization_tuple_digest", "cap_vector_digest", "catalog_provenance_digest",
    "data_root_storage_contract_digest", "environment_digest", "input_bundle_digest",
    "manifest_digest", "measurement_contract_digest", "pair_evidence_digest",
    "retrieval_checksum_obligations_digest", "root_contract_digest", "row_dispositions_digest",
])
def test_every_altered_bound_digest_refuses_zero_use(tmp_path, digest_key):
    payload = _gate()
    payload["authorization_boundary"][digest_key] = "0" * 64
    calls = []
    gate_path = _write_gate(tmp_path, payload, f"{digest_key}.json")
    result = _run(tmp_path / "out", gate_path=gate_path, submitter=lambda argv: calls.append(argv) or "1")
    assert result["status"] == "refused_preflight"
    assert result["failure_code"] == "BOUND_DIGEST_MISMATCH"
    assert result["slurm_jobs_submitted"] == result["core_seconds"] == 0
    assert calls == []


def _acquisition_fixture(tmp_path: Path):
    tmp_path.mkdir(parents=True, exist_ok=True)
    content = b"exact immutable payload"
    obj = tmp_path / "object"
    obj.write_bytes(content)
    obj.chmod(stat.S_IRUSR | stat.S_IRGRP)
    digest = hashlib.sha256(content).hexdigest()
    candidate = "candidate_a"
    accession = "GCF_000000001.1"
    url = f"https://example.invalid/{accession}/x.gz"
    expected = [{
        "candidate_id": candidate, "role": "h1_fasta", "accession_version": accession,
        "url": url, "expected_size_bytes": len(content),
    }]
    payload = {
        "decision_sha256": "1" * 64,
        "authorization_boundary": {"authorization_tuple_digest": "2" * 64},
    }
    ledger_fields = acquisition.MANIFEST_FIELDS
    ledger_row = {field: "" for field in ledger_fields}
    ledger_row.update({
        "record_type": "asset", "status": "promoted", "candidate_id": candidate,
        "asset_role": "h1_fasta", "accession_version": accession, "source_url": url,
        "observed_bytes": len(content), "observed_sha256": digest,
        "validation_outcomes_json": json.dumps(["local_sha256_reverified", "atomic_read_only_promotion"]),
    })
    inventory_row = {
        "run_id": "r", "promoted_at_utc": "now", "candidate_id": candidate,
        "asset_role": "h1_fasta", "accession_version": accession, "source_url": url,
        "bytes": len(content), "local_sha256": digest, "source_checksum_algorithm": "",
        "source_checksum_value": "", "source_checksum_verified": "not_available",
        "object_path": str(obj), "mode_octal": "0o440", "authorization_tuple_digest": "2" * 64,
        "decision_sha256": "1" * 64,
    }
    ledger = tmp_path / "acquisition.tsv"
    inventory = tmp_path / "inventory.tsv"
    runner.atomic_write_tsv(ledger, ledger_fields, [ledger_row])
    runner.atomic_write_tsv(inventory, acquisition.INVENTORY_FIELDS, [inventory_row])
    return payload, expected, ledger, inventory, obj


def test_altered_local_payload_digest_refuses_acquisition_audit(tmp_path):
    payload, expected, ledger, inventory, obj = _acquisition_fixture(tmp_path)
    obj.chmod(stat.S_IRUSR | stat.S_IWUSR)
    altered = bytearray(obj.read_bytes())
    altered[0] ^= 1
    obj.write_bytes(altered)
    obj.chmod(stat.S_IRUSR | stat.S_IRGRP)
    with pytest.raises(Tier3ValidationError, match="local payload SHA-256"):
        runner.audit_acquisition(payload, expected, ledger, inventory)


def test_incomplete_acquisition_and_accession_substitution_refuse(tmp_path):
    payload, expected, ledger, inventory, _ = _acquisition_fixture(tmp_path)
    runner.atomic_write_tsv(inventory, acquisition.INVENTORY_FIELDS, [])
    with pytest.raises(Tier3ValidationError, match="incomplete"):
        runner.audit_acquisition(payload, expected, ledger, inventory)

    payload, expected, ledger, inventory, _ = _acquisition_fixture(tmp_path / "sub")
    rows = _read_tsv(inventory)
    rows[0]["accession_version"] = "GCF_000000002.1"
    runner.atomic_write_tsv(inventory, acquisition.INVENTORY_FIELDS, rows)
    with pytest.raises(Tier3ValidationError, match="substitution"):
        runner.audit_acquisition(payload, expected, ledger, inventory)


def test_native_annotation_linkage_is_required(tmp_path):
    payload, expected, ledger, inventory, _ = _acquisition_fixture(tmp_path)
    expected[0]["role"] = "native_h1_annotation"
    ledger_rows = _read_tsv(ledger)
    ledger_rows[0]["asset_role"] = "native_h1_annotation"
    inventory_rows = _read_tsv(inventory)
    inventory_rows[0]["asset_role"] = "native_h1_annotation"
    runner.atomic_write_tsv(ledger, acquisition.MANIFEST_FIELDS, ledger_rows)
    runner.atomic_write_tsv(inventory, acquisition.INVENTORY_FIELDS, inventory_rows)
    with pytest.raises(Tier3ValidationError, match="native annotation"):
        runner.audit_acquisition(payload, expected, ledger, inventory)


def test_tier3a_pair_must_be_same_individual_and_correctly_phased():
    payload = {
        "row_audit": {"rows": [{"candidate_id": "a", "tier3a_required": True}]},
        "pair_evidence": {"rows": [{
            "candidate_id": "a", "h1_accession_version": "GCA_1.1", "h2_accession_version": "GCA_2.1",
            "same_individual_status": "unknown", "phase_evidence_status": "unknown",
        }]},
    }
    with pytest.raises(Tier3ValidationError, match="Tier3A pair mismatch"):
        runner.audit_pair_contract(payload, ["a"])
    payload["pair_evidence"]["rows"][0].update({"same_individual_status": "affirmed", "phase_evidence_status": "correctly_phased"})
    runner.audit_pair_contract(payload, ["a"])


def test_weakened_cap_and_storage_contract_refuse():
    payload = _gate()
    weakened = deepcopy(payload)
    weakened["cap_vector"]["dimensions"]["species"]["limit"] = 7
    with pytest.raises(Tier3ValidationError, match="cap vector relaxed"):
        acquisition.verify_cap_contract(weakened)
    weakened = deepcopy(payload)
    weakened["cap_vector"]["operational_thresholds"]["quota_headroom_fraction_minimum"] = 0.1
    with pytest.raises(Tier3ValidationError, match="storage headroom"):
        acquisition.verify_cap_contract(weakened)


def test_every_integrated_io_storage_and_network_dimension_is_mechanically_enforced():
    payload = _gate()
    for dimension in payload["cap_vector"]["dimensions"].values():
        dimension["observed"] = min(float(dimension["observed"]), float(dimension["limit"]))
        dimension["passes"] = True
    payload["cap_vector"]["proposed_metadata_ready"].update({
        "core_hours": 11.0, "scratch_gib": 100.0, "memory_per_job_gib": 32.0,
        "cpus_per_element": 8.0, "aggregate_wall_hours": 8.0,
    })
    payload["storage_audit"]["adequate"] = True
    payload["storage_audit"]["enforceable_allocation"].update({"status": "known", "headroom_pass": True})
    plan = runner.derive_resource_plan(payload, [f"candidate_{i}" for i in range(6)])
    assert plan["concurrency"] == 2
    assert plan["total_core_hours"] == 11.0

    for dimension_name in ("moosefs_read_gb", "moosefs_write_gb", "metadata_operations", "peak_bandwidth_mib_s", "persistent_output_gb", "file_inodes"):
        weakened = deepcopy(payload)
        dimension = weakened["cap_vector"]["dimensions"][dimension_name]
        dimension["observed"] = float(dimension["limit"]) + 1
        dimension["passes"] = False
        with pytest.raises(Tier3ValidationError, match=dimension_name):
            runner.derive_resource_plan(weakened, [f"candidate_{i}" for i in range(6)])

    weakened = deepcopy(payload)
    weakened["storage_audit"]["enforceable_allocation"]["headroom_pass"] = False
    with pytest.raises(Tier3ValidationError, match="storage contract"):
        runner.derive_resource_plan(weakened, [f"candidate_{i}" for i in range(6)])


def test_denominators_and_target_totals_are_measured_thresholded_and_never_imputed():
    contract = _gate()["measurement_contract"]
    packet = {
        "callable_bases": 0, "callable_fraction": None, "queryable_gene_bases": 999999,
        "queryable_gene_count": 999, "target_gene_total": 0, "target_base_total": None,
        "measurement_method": "", "artifact_sha256": "",
    }
    exclusions = runner.validate_denominator_packet(packet, contract)
    assert {row["metric"] for row in exclusions} == {
        "callable_bases", "callable_fraction", "queryable_gene_bases", "queryable_gene_count",
        "target_gene_total", "target_base_total", "denominator_provenance",
    }
    accepted = {
        "callable_bases": 10_000_000, "callable_fraction": 0.5, "queryable_gene_bases": 1_000_000,
        "queryable_gene_count": 1000, "target_gene_total": 1200, "target_base_total": 2_000_000,
        "measurement_method": "post-alignment PAF/callability and native-CDS union",
        "artifact_sha256": "a" * 64,
    }
    assert runner.validate_denominator_packet(accepted, contract) == []


def test_approved_sweepga_impg_and_worker_contracts_are_exact():
    runner.audit_sweepga_origin_build(ROOT / "analysis/sweepga_origin_main_build.json")
    runner.audit_impg_handoff(ROOT / "analysis/sweepga_impg_observed.json")
    worker = (ROOT / "analysis/slurm/run_repaired_vgp.sh").read_text(encoding="utf-8")
    compute = (ROOT / "analysis/vgp_pilot_compute.py").read_text(encoding="utf-8")
    assert "SLURM_TMPDIR" in worker
    assert "unset http_proxy https_proxy" in worker
    assert "sha256sum" in worker and "atomic promotion" in worker
    assert "--num-mappings\", \"1:1" in compute
    # The success-sentinel requirements live in the worker, while literal IMPG
    # commands live in the compute boundary.
    assert all(token in worker for token in ("normalized_vcf_tbi", "normalized_bcf_csi", "denominators_measured", "no_demographic_inference"))
    assert all(token in compute for token in ("impg", "partition", "query", "vcf:poa"))
