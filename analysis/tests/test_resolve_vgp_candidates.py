import json
import csv
from pathlib import Path

import pytest

from analysis import resolve_vgp_candidates as resolver


def _catalog_text():
    names = ["Species {0:04d}".format(index) for index in range(712)]
    names.extend(["Lophostoma evotis", "Micronycteris microtis"])
    rows = list(names)
    rows.extend(["Lophostoma evotis", "Micronycteris microtis"])
    return "Scientific Name\tStatus\n" + "\n".join(
        "{0}\t4".format(name) for name in rows
    ) + "\n"


def test_catalog_units_distinguish_lines_rows_and_unique_species():
    stats = resolver.catalog_statistics(_catalog_text())
    assert stats == {
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


def test_normalized_request_cache_is_immutable_and_provenance_stamped(tmp_path):
    cache = resolver.ImmutableResponseCache(tmp_path)
    request = resolver.NormalizedRequest(
        source="NCBI Datasets",
        source_version="18.33.1",
        method="GET",
        endpoint="https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession",
        parameters={"accessions": ["GCF_2.1", "GCF_1.1"], "page_size": 20},
    )
    equivalent = resolver.NormalizedRequest(
        source="NCBI Datasets",
        source_version="18.33.1",
        method="get",
        endpoint="https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/",
        parameters={"page_size": 20, "accessions": ["GCF_1.1", "GCF_2.1"]},
    )
    assert request.cache_key == equivalent.cache_key

    entry = cache.store(
        request,
        b'{"reports":[]}',
        response_headers={"x-datasets-version": "18.33.1"},
        retrieved_at_utc="2026-07-18T00:00:00Z",
        software_environment={"software": "test-resolver/1", "python": "pinned"},
    )
    assert entry["request"] == request.as_dict()
    assert entry["endpoint"] == request.endpoint
    assert entry["retrieved_at_utc"] == "2026-07-18T00:00:00Z"
    assert len(entry["response_sha256"]) == 64
    assert entry["software_environment"]["software"] == "test-resolver/1"
    assert cache.load(equivalent)[0] == b'{"reports":[]}'

    with pytest.raises(resolver.CacheConflictError, match="immutable cache conflict"):
        cache.store(
            request,
            b'{"reports":[{"changed":true}]}',
            response_headers={},
            retrieved_at_utc="2026-07-18T00:00:01Z",
            software_environment={"software": "test-resolver/1"},
        )


class _Response:
    def __init__(self, status, body=b"", headers=None):
        self.status = status
        self.body = body
        self.headers = headers or {}


def test_batched_retry_after_checkpoint_resume_and_cache(tmp_path):
    calls = []
    sleeps = []
    responses = [
        _Response(429, headers={"Retry-After": "2"}),
        _Response(503),
        _Response(200, b'{"reports":[1,2]}', {"x-datasets-version": "18.33.1"}),
        _Response(200, b'{"reports":[3]}', {"x-datasets-version": "18.33.1"}),
    ]

    def transport(request):
        calls.append(request)
        return responses.pop(0)

    client = resolver.BatchedCachedClient(
        cache=resolver.ImmutableResponseCache(tmp_path / "cache"),
        checkpoint_path=tmp_path / "checkpoint.json",
        transport=transport,
        sleep=sleeps.append,
        jitter=lambda _attempt: 0.0,
        minimum_interval_seconds=0.0,
        maximum_attempts=4,
        maximum_backoff_seconds=8.0,
        software_environment={"software": "test-resolver/1"},
    )
    results = client.fetch_accession_batches(
        accessions=["GCF_3.1", "GCF_1.1", "GCF_2.1"],
        batch_size=2,
        request_factory=lambda batch: resolver.NormalizedRequest(
            source="NCBI Datasets",
            source_version="18.33.1",
            method="GET",
            endpoint="https://example.test/dataset_report",
            parameters={"accessions": list(batch)},
        ),
    )
    assert len(results) == 2
    assert [call.parameters["accessions"] for call in calls] == [
        ["GCF_1.1", "GCF_2.1"],
        ["GCF_1.1", "GCF_2.1"],
        ["GCF_1.1", "GCF_2.1"],
        ["GCF_3.1"],
    ]
    assert sleeps == [2.0, 2.0]
    checkpoint = json.loads((tmp_path / "checkpoint.json").read_text())
    assert len(checkpoint["completed_request_keys"]) == 2

    def must_not_call(_request):
        raise AssertionError("resume must use immutable cached responses")

    resumed = resolver.BatchedCachedClient(
        cache=resolver.ImmutableResponseCache(tmp_path / "cache"),
        checkpoint_path=tmp_path / "checkpoint.json",
        transport=must_not_call,
        sleep=sleeps.append,
        jitter=lambda _attempt: 0.0,
        minimum_interval_seconds=0.0,
        software_environment={"software": "test-resolver/1"},
    ).fetch_accession_batches(
        accessions=["GCF_1.1", "GCF_2.1", "GCF_3.1"],
        batch_size=2,
        request_factory=lambda batch: resolver.NormalizedRequest(
            source="NCBI Datasets",
            source_version="18.33.1",
            method="GET",
            endpoint="https://example.test/dataset_report",
            parameters={"accessions": list(batch)},
        ),
    )
    assert [body for body, _entry in resumed] == [b'{"reports":[1,2]}', b'{"reports":[3]}']


def test_missing_remote_checksum_becomes_staged_local_sha256_obligation():
    obligations = resolver.payload_acquisition_obligations(
        accession_version="GCF_000001.2",
        url="https://ftp.ncbi.nlm.nih.gov/exact/GCF_000001.2_genomic.fna.gz",
        expected_size_bytes=12345,
        official_checksum_algorithm=None,
        official_checksum=None,
    )
    assert obligations["pre_download_eligible"] is True
    assert obligations["remote_checksum_required_for_pre_download"] is False
    assert obligations["official_checksum_verification_required"] is False
    assert obligations["steps"] == [
        "stage_full_payload",
        "verify_expected_size_bytes:12345",
        "compute_local_sha256",
        "reverify_local_sha256_before_promotion",
        "atomic_promote_read_only",
    ]

    with_checksum = resolver.payload_acquisition_obligations(
        accession_version="GCF_000001.2",
        url="https://ftp.ncbi.nlm.nih.gov/exact/GCF_000001.2_genomic.gff.gz",
        expected_size_bytes=987,
        official_checksum_algorithm="md5",
        official_checksum="0" * 32,
    )
    assert with_checksum["official_checksum_verification_required"] is True
    assert "verify_official_md5:" + "0" * 32 in with_checksum["steps"]


def test_denominators_are_post_alignment_acceptance_obligations():
    contract = resolver.POST_ALIGNMENT_MEASUREMENT_CONTRACT
    assert contract["phase"] == "post_alignment_pre_result_acceptance"
    assert contract["executable_command"]

    predownload = resolver.pre_download_eligibility(
        {
            "h1_accession_version": "GCF_000001.2",
            "ncbi_taxid": "123",
            "annotation_accession_version": "GCF_000001.2-RS_2026_01",
            "annotation_reference_accession_version": "GCF_000001.2",
            "annotation_sequence_region_linkage_status": "proven_official_exact_h1",
        },
        require_diversity=False,
    )
    assert predownload["eligible"] is True
    assert not any("CALLABLE" in code or "QUERYABLE" in code for code in predownload["blockers"])

    absent = resolver.evaluate_post_alignment_measurements({})
    assert absent["accepted"] is False
    assert absent["result_disposition"] == "exclude"
    assert set(absent["failed_thresholds"]) == {
        "callable_bases_missing",
        "callable_fraction_missing",
        "queryable_gene_count_missing",
        "queryable_gene_bases_missing",
    }

    adequate = resolver.evaluate_post_alignment_measurements(
        {
            "callable_bases": 20_000_000,
            "callable_fraction": 0.75,
            "queryable_gene_count": 12_000,
            "queryable_gene_bases": 15_000_000,
        }
    )
    assert adequate == {
        "accepted": True,
        "failed_thresholds": [],
        "result_disposition": "include",
    }


def test_tier3a_pre_download_gate_requires_exact_phased_same_individual_pair():
    base = {
        "h1_accession_version": "GCA_000001.2",
        "ncbi_taxid": "123",
        "annotation_accession_version": "GCA_000001.2-ANN_1",
        "annotation_reference_accession_version": "GCA_000001.2",
        "annotation_sequence_region_linkage_status": "proven_official_exact_h1",
        "h2_accession_version": "GCA_000002.1",
        "same_individual_status": "yes",
        "phase_evidence_status": "affirmative_correctly_phased",
        "pair_evidence_url": "https://api.ncbi.nlm.nih.gov/pair",
    }
    assert resolver.pre_download_eligibility(base, require_diversity=True)["eligible"] is True
    broken = dict(base)
    broken["phase_evidence_status"] = "inferred_from_neighboring_accession"
    result = resolver.pre_download_eligibility(broken, require_diversity=True)
    assert result["eligible"] is False
    assert "TIER3A_PHASE_EVIDENCE_NOT_AFFIRMATIVE" in result["blockers"]


def test_integrated_cache_index_and_repaired_outputs_are_reproducible():
    index_path = Path("analysis/vgp_resolution_cache/index.json")
    assert index_path.is_file()
    index = resolver.validate_resolution_cache_index(index_path)
    assert index["catalog_statistics"]["physical_lines"] == 717
    assert index["catalog_statistics"]["header_lines"] == 1
    assert index["catalog_statistics"]["data_rows"] == 716
    assert index["catalog_statistics"]["unique_species"] == 714
    assert index["catalog_statistics"]["data_row_excess_over_unique_species"] == 2
    assert index["selection_comparison"]["old_selected_count"] == 0
    repaired_selected = index["selection_comparison"]["repaired_selected_count"]
    assert 0 <= repaired_selected <= 6
    if repaired_selected == 0:
        assert index["selection_comparison"]["zero_selection_irreducible_blockers"] == [
            "QUOTA_UNAVAILABLE"
        ]
    assert index["authorization_boundary"]["biological_payloads_downloaded"] == 0
    assert index["authorization_boundary"]["jobs_submitted"] == 0
    assert index["authorization_boundary"]["demographic_inference_performed"] is False
    with Path("analysis/vgp_pilot_manifest.tsv").open(newline="") as handle:
        manifest = list(csv.DictReader(handle, delimiter="\t"))
    assert len(manifest) == 6
    assert all(row["h1_accession_version"].startswith("GCF_") for row in manifest)
    assert all(
        row["annotation_reference_accession_version"] == row["h1_accession_version"]
        for row in manifest
    )
    assert all(
        row["annotation_sequence_region_linkage_status"] == "proven_official_exact_h1"
        for row in manifest
    )
    assert all(row["resolved_modality"] == "tier3c_composition" for row in manifest)
    assert all(row["assembly_diversity_eligible"] == "no" for row in manifest)
