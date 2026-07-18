import copy
import csv
import hashlib
import json
from pathlib import Path

import pytest

from analysis import gate_vgp_pilot as gate
from analysis.tier3_common import Tier3ValidationError


ROOT = Path(__file__).parents[2]


def _paths(tmp_path: Path) -> dict:
    return {
        "manifest_path": ROOT / "analysis/vgp_pilot_manifest.tsv",
        "rejections_path": ROOT / "analysis/vgp_pilot_rejections.tsv",
        "size_budget_path": ROOT / "analysis/vgp_pilot_size_budget.tsv",
        "resolution_index_path": ROOT / "analysis/vgp_resolution_cache/index.json",
        "freeze_provenance_path": ROOT / "analysis/vgp_phase1_freeze_provenance.json",
        "root_config_path": ROOT / "analysis/vgp_data_root_config.json",
        "root_validation_path": ROOT / "analysis/vgp_data_root_validation.json",
        "decisions_path": ROOT / "analysis/vertebrate_scaleout_decisions.tsv",
        "execution_plan_path": ROOT / "analysis/vertebrate_scaleout_execution_plan.md",
        "resource_budget_path": ROOT / "analysis/vertebrate_scaleout_resource_budget.tsv",
        "guix_channels_path": ROOT / "analysis/guix/channels.scm",
        "guix_manifest_path": ROOT / "analysis/guix/manifest.scm",
        "guix_environment_path": ROOT / "analysis/pilot_results/guix_environment.json",
        "gate_out": tmp_path / "vgp_pilot_gate.json",
        "review_out": tmp_path / "vgp_pilot_gate_review.md",
    }


def _build(tmp_path: Path, **overrides):
    kwargs = _paths(tmp_path)
    kwargs.update(overrides)
    return gate.build_gate(**kwargs), kwargs


def _write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _mutated_tsv(tmp_path: Path, source: Path, mutate) -> Path:
    rows = gate.load_tsv(source)
    mutate(rows)
    target = tmp_path / ("mutated-" + source.name)
    _write_tsv(target, rows)
    return target


def _mutated_json(tmp_path: Path, source: Path, mutate, name: str | None = None) -> Path:
    payload = json.loads(source.read_text(encoding="utf-8"))
    mutate(payload)
    target = tmp_path / (name or ("mutated-" + source.name))
    target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return target


def _authorize(kwargs: dict, **overrides):
    call = {
        key: value
        for key, value in kwargs.items()
        if key not in {"gate_out", "review_out"}
    }
    call.update(overrides)
    return gate.authorize_gate_action(
        kwargs["gate_out"],
        call.pop("manifest_path"),
        call.pop("root_config_path"),
        "acquire",
        **call,
    )


def test_regenerated_gate_reproduces_catalog_rows_duplicates_and_dispositions(tmp_path):
    built, kwargs = _build(tmp_path)
    assert kwargs["gate_out"].is_file()
    assert kwargs["review_out"].is_file()
    assert built["decision"]["status"] == "NO_GO"
    assert built["catalog_audit"]["statistics"] == {
        "physical_lines": 717,
        "header_lines": 1,
        "data_rows": 716,
        "unique_species": 714,
        "data_row_excess_over_unique_species": 2,
        "duplicated_species": [
            {"scientific_name": "Lophostoma evotis", "multiplicity": 2},
            {"scientific_name": "Micronycteris microtis", "multiplicity": 2},
        ],
    }
    assert built["row_audit"]["summary"]["metadata_ready_count"] == 6
    assert built["row_audit"]["summary"]["tier3a_ready_count"] == 0
    assert built["selection_audit"]["observed_selected_candidate_ids"] == []
    assert built["selection_audit"]["independently_expected_selected_candidate_ids"] == []
    assert built["disposition_audit"]["seed_row_count"] == 74
    assert built["disposition_audit"]["rejection_row_count"] == 74
    assert built["disposition_audit"]["all_rows_match"] is True
    assert len(built["disposition_audit"]["rows"]) == 74
    review = kwargs["review_out"].read_text(encoding="utf-8")
    assert "717 physical lines" in review
    assert "Lophostoma evotis" in review
    assert "NO_GO" in review


def test_gate_recomputes_strictest_caps_and_does_not_promote_free_space_to_quota(tmp_path):
    built, _kwargs = _build(tmp_path)
    caps = built["cap_vector"]["dimensions"]
    assert caps["species"]["limit"] == 6.0
    assert caps["compressed_inputs_gib"]["limit"] == 120.0
    assert caps["scratch_gib"]["limit"] == pytest.approx(150 / gate.BYTES_PER_GIB * 1_000_000_000)
    assert caps["core_hours"]["limit"] == 280.0
    assert caps["concurrent_species"]["limit"] == 2.0
    assert caps["memory_per_job_gib"]["limit"] == 96.0
    assert caps["scratch_gib"]["observed"] == pytest.approx(205.5066)
    assert caps["moosefs_read_gb"]["observed"] == pytest.approx(406.7274)
    codes = {item["code"] for item in built["blockers"]}
    assert {"QUOTA_UNAVAILABLE", "CAP_SCRATCH_GIB_EXCEEDED", "CAP_MOOSEFS_READ_GB_EXCEEDED"} <= codes
    storage = built["storage_audit"]
    assert storage["filesystem"]["free_bytes"] > 0
    assert storage["filesystem"]["free_inodes"] > 0
    assert storage["enforceable_allocation"]["status"] == "unknown"
    assert storage["adequate"] is False
    assert storage["headroom_fraction_required"] == 0.25


@pytest.mark.parametrize(
    ("mutation", "expected_code"),
    [
        (lambda rows: rows[0].__setitem__("h1_accession_version", "GCF_036321535"), "H1_ACCESSION_INVALID"),
        (
            lambda rows: rows[0].__setitem__("annotation_reference_accession_version", "GCF_000000001.1"),
            "ANNOTATION_REFERENCE_MISMATCH",
        ),
        (
            lambda rows: rows[0].__setitem__("annotation_sequence_region_linkage_status", "not_proven"),
            "ANNOTATION_LINKAGE_NOT_PROVEN",
        ),
    ],
)
def test_missing_exact_h1_or_native_annotation_linkage_is_nogo(tmp_path, mutation, expected_code):
    manifest = _mutated_tsv(tmp_path, ROOT / "analysis/vgp_pilot_manifest.tsv", mutation)
    built, _kwargs = _build(tmp_path, manifest_path=manifest)
    assert built["decision"]["status"] == "NO_GO"
    codes = {issue["code"] for row in built["row_audit"]["rows"] for issue in row["issues"]}
    assert expected_code in codes


def test_tier3a_requires_exact_h2_same_individual_and_phasing_evidence(tmp_path):
    def mutation(rows):
        rows[0]["resolved_modality"] = "tier3a_diversity"
        rows[0]["h2_accession_version"] = ""
        rows[0]["same_individual_status"] = "not_proven"
        rows[0]["phase_evidence_status"] = ""
        rows[0]["pair_evidence_url"] = ""

    manifest = _mutated_tsv(tmp_path, ROOT / "analysis/vgp_pilot_manifest.tsv", mutation)
    built, _kwargs = _build(tmp_path, manifest_path=manifest)
    issues = {
        issue["code"]
        for row in built["row_audit"]["rows"]
        if row["candidate_id"] == "camelus_dromedarius_gca_036321535_1"
        for issue in row["issues"]
    }
    assert {
        "H2_ACCESSION_INVALID",
        "PAIR_NOT_SAME_INDIVIDUAL",
        "PAIR_PHASE_EVIDENCE_MISSING",
        "PAIR_EVIDENCE_URL_MISSING",
    } <= issues
    assert built["decision"]["status"] == "NO_GO"


def test_finite_sizes_and_staged_local_sha256_obligations_are_enforced(tmp_path):
    built, _kwargs = _build(tmp_path)
    assert built["retrieval_audit"]["ready_count"] == 6
    for row in built["retrieval_audit"]["rows"]:
        assert row["pre_download_ready"] is True
        for obligation in row["obligations"]:
            assert obligation["local_sha256_after_staging_required"] is True
            assert obligation["local_sha256_reverification_required"] is True
            assert obligation["remote_checksum_required_for_pre_download"] is False
            assert obligation["expected_size_bytes"] > 0
            if obligation["source_checksum"]:
                assert obligation["source_checksum_verified_against_official_catalog"] is True

    def mutation(rows):
        obligations = json.loads(rows[0]["acquisition_obligations"])
        obligations[0]["expected_size_bytes"] = None
        rows[0]["acquisition_obligations"] = json.dumps(obligations)

    manifest = _mutated_tsv(tmp_path, ROOT / "analysis/vgp_pilot_manifest.tsv", mutation)
    failed, _kwargs = _build(tmp_path, manifest_path=manifest)
    codes = {issue["code"] for row in failed["row_audit"]["rows"] for issue in row["issues"]}
    assert "RETRIEVAL_SIZE_INVALID" in codes
    assert failed["decision"]["status"] == "NO_GO"

    missing_source_checksum = _mutated_tsv(
        tmp_path,
        ROOT / "analysis/vgp_pilot_manifest.tsv",
        lambda rows: rows[0].__setitem__("h1_provider_md5", ""),
    )
    checksum_failed, _kwargs = _build(tmp_path, manifest_path=missing_source_checksum)
    checksum_codes = {issue["code"] for row in checksum_failed["row_audit"]["rows"] for issue in row["issues"]}
    assert "RETRIEVAL_SOURCE_CHECKSUM_NOT_VERIFIED" in checksum_codes


def test_denominators_are_post_alignment_acceptance_not_predownload_prerequisites(tmp_path):
    built, _kwargs = _build(tmp_path)
    contract = built["measurement_contract"]
    assert contract["pre_download_prerequisite"] is False
    assert contract["phase"] == "post_alignment_pre_result_acceptance"
    assert built["row_audit"]["summary"]["metadata_ready_count"] == 6
    assert not any(
        issue["code"].startswith(("CALLABLE", "QUERYABLE"))
        for row in built["row_audit"]["rows"]
        for issue in row["issues"]
    )
    missing = gate.evaluate_post_alignment_measurements(contract, {})
    assert missing["accepted"] is False
    assert missing["result_disposition"] == "exclude_downstream_result"
    inadequate = gate.evaluate_post_alignment_measurements(
        contract,
        {"callable_bases": 20_000_000, "callable_fraction": 0.49, "queryable_gene_count": 2_000, "queryable_gene_bases": 2_000_000},
    )
    assert inadequate["accepted"] is False
    assert "callable_fraction_below_minimum" in inadequate["failed_thresholds"]
    adequate = gate.evaluate_post_alignment_measurements(
        contract,
        {"callable_bases": 20_000_000, "callable_fraction": 0.75, "queryable_gene_count": 2_000, "queryable_gene_bases": 2_000_000},
    )
    assert adequate["accepted"] is True


def test_missing_storage_safety_and_cap_compliance_are_nogo(tmp_path):
    root_validation = _mutated_json(
        tmp_path,
        ROOT / "analysis/vgp_data_root_validation.json",
        lambda payload: payload["smoke_tests"]["atomic_promotion"].__setitem__("status", "fail"),
    )
    built, _kwargs = _build(tmp_path, root_validation_path=root_validation)
    assert "ROOT_ATOMIC_PROMOTION_UNSAFE" in {item["code"] for item in built["blockers"]}
    assert built["decision"]["status"] == "NO_GO"


@pytest.mark.parametrize(
    ("kind", "expected"),
    [
        ("manifest", "manifest digest mismatch"),
        ("catalog", "catalog provenance digest mismatch"),
        ("root", "root/storage contract digest mismatch"),
        ("environment", "environment digest mismatch"),
        ("cap", "cap vector digest mismatch"),
        ("retrieval", "retrieval/checksum obligations digest mismatch"),
        ("pair", "pair evidence digest mismatch"),
        ("measurement", "measurement contract digest mismatch"),
    ],
)
def test_authorizer_refuses_every_altered_bound_component(tmp_path, kind, expected):
    _built, kwargs = _build(tmp_path)
    overrides = {}
    if kind == "manifest":
        overrides["manifest_path"] = _mutated_tsv(
            tmp_path,
            kwargs["manifest_path"],
            lambda rows: rows[0].__setitem__("uncertainty_status", "altered"),
        )
    elif kind == "catalog":
        provenance = json.loads(kwargs["freeze_provenance_path"].read_text(encoding="utf-8"))
        source = Path(provenance["source_catalog"]["path"])
        catalog = tmp_path / "catalog.tsv"
        catalog.write_text(source.read_text(encoding="utf-8").replace("Camelus dromedarius", "Camelus altered", 1), encoding="utf-8")
        provenance["source_catalog"]["path"] = str(catalog)
        provenance["source_catalog"]["sha256"] = hashlib.sha256(catalog.read_bytes()).hexdigest()
        changed = tmp_path / "freeze.json"
        changed.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        overrides["freeze_provenance_path"] = changed
    elif kind == "root":
        overrides["root_validation_path"] = _mutated_json(
            tmp_path,
            kwargs["root_validation_path"],
            lambda payload: payload["root_owner"].__setitem__("world_writable", True),
        )
    elif kind == "environment":
        overrides["guix_environment_path"] = _mutated_json(
            tmp_path,
            kwargs["guix_environment_path"],
            lambda payload: payload.__setitem__("channel_commit", "0" * 40),
        )
    elif kind == "cap":
        overrides["decisions_path"] = _mutated_tsv(
            tmp_path,
            kwargs["decisions_path"],
            lambda rows: next(row for row in rows if row["decision_id"] == "D012").__setitem__(
                "resolution_hard_gate",
                next(row for row in rows if row["decision_id"] == "D012")["resolution_hard_gate"].replace("280 core-h", "279 core-h"),
            ),
        )
    elif kind == "retrieval":
        def mutate_retrieval(rows):
            obligations = json.loads(rows[0]["acquisition_obligations"])
            obligations[0]["steps"].append("altered_step")
            rows[0]["acquisition_obligations"] = json.dumps(obligations)
        overrides["manifest_path"] = _mutated_tsv(tmp_path, kwargs["manifest_path"], mutate_retrieval)
    elif kind == "pair":
        overrides["manifest_path"] = _mutated_tsv(
            tmp_path,
            kwargs["manifest_path"],
            lambda rows: rows[0].__setitem__("same_individual_evidence", "altered"),
        )
    elif kind == "measurement":
        overrides["resolution_index_path"] = _mutated_json(
            tmp_path,
            kwargs["resolution_index_path"],
            lambda payload: payload["post_alignment_measurement_contract"]["minimum_thresholds"].__setitem__("callable_bases", 9_999_999),
        )
    with pytest.raises(Tier3ValidationError, match=expected):
        _authorize(kwargs, **overrides)


def test_tampered_gate_and_nogo_are_both_refused(tmp_path):
    built, kwargs = _build(tmp_path)
    assert built["decision"]["status"] == "NO_GO"
    with pytest.raises(Tier3ValidationError, match="gate decision is NO_GO"):
        _authorize(kwargs)

    payload = json.loads(kwargs["gate_out"].read_text(encoding="utf-8"))
    payload["decision"]["status"] = "GO"
    kwargs["gate_out"].write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    with pytest.raises(Tier3ValidationError, match="gate decision hash does not match"):
        _authorize(kwargs)
