from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import pytest

from analysis import vgp_10_pilot_acquisition as acq
from analysis import validate_vgp_10_pilot_acquisition as validator


def _object(tmp_path: Path, payload: bytes = b"immutable pilot object") -> dict[str, object]:
    return {
        "object_id": "P01:h1:genome_fasta",
        "selection_id": "P01",
        "object_role": "genome_fasta",
        "source_url": "https://example.invalid/GCA_000000001.1_genomic.fna.gz",
        "expected_bytes": len(payload),
        "upstream_checksum_algorithm": "md5",
        "upstream_checksum": hashlib.md5(payload).hexdigest(),  # noqa: S324
        "expected_local_sha256": hashlib.sha256(payload).hexdigest(),
    }


def test_promote_reuses_mirror_only_after_local_sha256_revalidation(tmp_path: Path) -> None:
    payload = b"verified mirror bytes"
    row = _object(tmp_path, payload)
    cas = tmp_path / "objects" / "sha256"
    target = acq.cas_path(cas, row["expected_local_sha256"])
    target.parent.mkdir(parents=True)
    target.write_bytes(payload)

    outcome = acq.acquire_one(row, tmp_path, downloader=lambda *_: pytest.fail("downloaded"))

    assert outcome["status"] == "reused"
    assert outcome["local_sha256"] == row["expected_local_sha256"]
    assert outcome["revalidation_count"] == 1
    assert outcome["transferred_bytes"] == 0


def test_resume_partial_then_verify_and_atomically_promote(tmp_path: Path) -> None:
    payload = b"0123456789abcdef"
    row = _object(tmp_path, payload)
    partial = acq.partial_path(tmp_path, row["object_id"])
    partial.parent.mkdir(parents=True)
    partial.write_bytes(payload[:5])
    calls: list[int] = []

    def resume(_row: dict[str, object], path: Path, offset: int) -> int:
        calls.append(offset)
        with path.open("ab") as handle:
            handle.write(payload[offset:])
        return len(payload) - offset

    outcome = acq.acquire_one(row, tmp_path, downloader=resume)

    assert calls == [5]
    assert outcome["status"] == "verified"
    assert outcome["resume_from_bytes"] == 5
    assert outcome["transferred_bytes"] == len(payload) - 5
    promoted = Path(outcome["local_path"])
    assert promoted.read_bytes() == payload
    assert not partial.exists()


def test_checksum_mismatch_is_quarantined_and_never_promoted(tmp_path: Path) -> None:
    payload = b"expected"
    row = _object(tmp_path, payload)

    def corrupt(_row: dict[str, object], path: Path, offset: int) -> int:
        assert offset == 0
        path.write_bytes(b"corrupt!")
        return 8

    outcome = acq.acquire_one(row, tmp_path, downloader=corrupt)

    assert outcome["status"] == "quarantined"
    assert Path(outcome["quarantine_path"]).is_file()
    assert not outcome["local_path"]
    assert not acq.cas_path(tmp_path / "objects" / "sha256", row["expected_local_sha256"]).exists()


def test_pair_identity_fails_closed_but_annotation_failure_is_branch_local() -> None:
    roster = {
        "species": "Testudo exacta",
        "resolved_taxid": "123",
        "biosample": "SAMN1",
        "individual_or_isolate": "animal1",
        "h1_accession_version": "GCA_000000001.1",
        "h2_accession_version": "GCA_000000002.1",
        "h1_role": "haplotype_1",
        "h2_role": "haplotype_2",
    }
    h1 = acq.synthetic_report(roster, "h1")
    h2 = acq.synthetic_report(roster, "h2")
    status = acq.validate_pair_identity(roster, h1, h2)
    assert status["core_identity_status"] == "pass"

    mismatched_annotation = acq.validate_annotation_binding(
        roster["h1_accession_version"], "GCF_000000099.1", "GCF_000000099.1-RS_1"
    )
    assert mismatched_annotation["annotation_branch_status"] == "failed_mismatch"
    assert status["core_identity_status"] == "pass"

    h2["assembly_info"]["biosample"]["accession"] = "SAMN_OTHER"
    with pytest.raises(acq.AcquisitionRefusal, match="BioSample"):
        acq.validate_pair_identity(roster, h1, h2)


def test_summary_uses_exclusive_terminal_dispositions(tmp_path: Path) -> None:
    rows = []
    for index, (status, size) in enumerate(
        [("verified", 10), ("reused", 20), ("missing", 30), ("superseded", 40), ("quarantined", 50)]
    ):
        rows.append(
            {
                "object_id": str(index),
                "status": status,
                "expected_bytes": size,
                "observed_bytes": size if status != "missing" else 0,
                "transferred_bytes": 5 if status in {"verified", "quarantined"} else 0,
            }
        )
    summary = acq.summarize_inventory(rows)
    assert summary["planned"]["objects"] == 5
    assert summary["planned"]["bytes"] == 150
    assert summary["terminal_disposition_reconciliation"]["objects"] == 5
    assert summary["terminal_disposition_reconciliation"]["bytes"] == 150
    assert summary["terminal_disposition_reconciliation"]["matches_planned"] is True
    assert summary["transferred"]["bytes"] == 10


def test_repository_deliverables_are_closed_world_and_no_bulk_reads() -> None:
    root = Path(__file__).resolve().parents[2]
    required = [
        root / "analysis/vgp_10_pilot_acquisition_manifest.tsv",
        root / "analysis/vgp_10_pilot_object_inventory.tsv",
        root / "analysis/vgp_direct_control_acquisition_manifest.tsv",
        root / "analysis/vgp_10_pilot_acquisition_summary.json",
        root / "analysis/vgp_10_pilot_acquisition_handoff.md",
    ]
    if not all(path.exists() for path in required):
        pytest.skip("generated acquisition deliverables are created by the pinned-Guix execution")

    with required[0].open(newline="", encoding="utf-8") as handle:
        pairs = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["selection_id"] for row in pairs if row["roster_type"] == "primary"] == [
        f"P{i:02d}" for i in range(1, 11)
    ]
    assert len([row for row in pairs if row["roster_type"] == "alternate"]) >= 5

    with required[1].open(newline="", encoding="utf-8") as handle:
        objects = list(csv.DictReader(handle, delimiter="\t"))
    assert len({row["object_id"] for row in objects}) == len(objects)
    assert not any(row["object_role"] in acq.BULK_RAW_ROLES for row in objects)
    assert all(row["source_url"].startswith("https://") for row in objects)
    assert all(row["expected_bytes"].isdigit() for row in objects)
    assert all(row["license_or_access"] for row in objects)

    summary = json.loads(required[3].read_text(encoding="utf-8"))
    assert summary["scope_proof"]["unmanifested_bulk_raw_read_objects_acquired"] == 0
    assert summary["scope_proof"]["slurm_jobs_submitted"] == 0
    assert summary["accounting"]["terminal_disposition_reconciliation"]["matches_planned"] is True


def test_repository_deliverables_pass_closed_world_validator() -> None:
    if not acq.SUMMARY_OUTPUT.is_file():
        pytest.skip("generated acquisition deliverables are created by the pinned-Guix execution")
    assert validator.validate(reverify_local=False) == []
