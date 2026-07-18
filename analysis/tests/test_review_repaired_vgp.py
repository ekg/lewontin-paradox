import csv
from pathlib import Path

from analysis import review_repaired_vgp as review


ROOT = Path(__file__).parents[2]


def read_tsv(path: Path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_repaired_review_reproduces_refusal_and_exact_demography_join(tmp_path):
    result = review.review(
        review_out=tmp_path / "review.md",
        qc_out=tmp_path / "qc.tsv",
        resource_out=tmp_path / "resource.tsv",
        guix_validation="PASS",
        guix_note="test fixture records the externally run pinned suite",
    )
    assert result["review_decision"] == "PASS"
    assert result["run_disposition"] == "NOT_SUBMITTED"

    qc_rows = read_tsv(tmp_path / "qc.tsv")
    by_id = {row["check_id"]: row for row in qc_rows}
    assert all(row["decision"] == "PASS" for row in qc_rows)
    assert by_id["gate_stable_recompute"]["decision"] == "PASS"
    assert by_id["current_run_refusal_recompute"]["decision"] == "PASS"
    assert by_id["refusal_evidence_sha256"]["decision"] == "PASS"
    assert by_id["promoted_file_sha256"]["decision"] == "PASS"
    refusal_rows = [row for row in qc_rows if row["check_id"].startswith("refusal:")]
    assert len(refusal_rows) == 14
    assert all("provider_calls=0; submit_calls=0" in row["observed"] for row in refusal_rows)
    assert by_id["refusal:current_no_go"]["observed"].startswith("acquire=GATE_NO_GO; run=GATE_NO_GO")
    assert by_id["refusal:unknown_decision"]["observed"].startswith("acquire=GATE_DECISION_UNKNOWN; run=GATE_DECISION_UNKNOWN")
    assert by_id["refusal:altered_manifest_digest"]["observed"].startswith("acquire=BOUND_DIGEST_MISMATCH; run=BOUND_DIGEST_MISMATCH")
    assert by_id["refusal:altered_retrieval_obligation"]["observed"].startswith(
        "acquire=RETRIEVAL_CONTRACT_MISMATCH; run=BOUND_DIGEST_MISMATCH"
    )
    assert by_id["immutable_biological_inventory"]["observed"] == "0"
    assert by_id["tier3a_pair_evidence_disposition"]["observed"] == "tier3a_ready=0; linked_H2_leads=6; validated_phased_pairs=0"
    assert by_id["demography_method_separation"]["decision"] == "PASS"
    assert "6 valid LD-Ne records for Camelus dromedarius" in by_id["independent_ne_and_circularity"]["observed"]

    text = (tmp_path / "review.md").read_text(encoding="utf-8")
    assert "PASS (correctly refused; `NOT_SUBMITTED`)" in text
    assert "no staged or promoted biological local SHA-256 values" in text
    assert "not performance calibration" in text
    assert "PSMC" in text and "MSMC2" in text and "SMC++" in text


def test_refusal_zero_use_is_not_resource_calibration(tmp_path):
    review.review(
        review_out=tmp_path / "review.md",
        qc_out=tmp_path / "qc.tsv",
        resource_out=tmp_path / "resource.tsv",
        guix_validation="PASS",
        guix_note="passed",
    )
    rows = read_tsv(tmp_path / "resource.tsv")
    calibration = [row for row in rows if row["scope"] == "successful_observation_calibration"]
    assert calibration
    assert all(row["predicted_low"] == row["predicted_base"] == row["predicted_high"] == "" for row in calibration)
    assert all(row["decision"] == "NOT_CALIBRATED" for row in calibration)
    assert all("excluded from calibration" in row["observed"] for row in calibration)
    projection = [row for row in rows if row["scope"] == "full_eligible_catalog_projection"]
    assert len(projection) == 1
    assert projection[0]["decision"] == "REQUIRES_NEW_AUTHORIZATION"


def test_repository_outputs_are_current_repaired_review():
    assert review.DEFAULT_REVIEW.name == "repaired_vgp_pilot_review.md"
    assert review.DEFAULT_QC.name == "repaired_vgp_pilot_qc.tsv"
    assert review.DEFAULT_RESOURCE.name == "repaired_vgp_resource_calibration.tsv"
