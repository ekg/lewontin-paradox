from __future__ import annotations

import hashlib
from pathlib import Path

import pytest

from analysis.compile_vgp_read_validation_results import p07_assessment, qv_effect, verify_output_manifest


def test_verify_output_manifest_rehashes_every_promoted_object(tmp_path: Path) -> None:
    payload = tmp_path / "result.json"
    payload.write_text('{"canonical_vgp_root":"/moosefs/erikg/vgp"}\n')
    digest = hashlib.sha256(payload.read_bytes()).hexdigest()
    (tmp_path / "output_manifest.tsv").write_text(
        f"sha256\trelative_path\n{digest}\tresult.json\n"
    )

    result = verify_output_manifest(tmp_path)
    assert result["verified_objects"] == 1
    assert result["all_digests_match"] is True

    payload.write_text("tampered\n")
    with pytest.raises(ValueError, match="output digest mismatch"):
        verify_output_manifest(tmp_path)


def test_effect_retains_core_and_reports_paired_ratio() -> None:
    assessment = p07_assessment(
        {"qv": 42.5},
        {"pi_ratio_read_over_assembly": 1.02},
        {"lambda_pearson_correlation": 0.9},
        {"concrete_false_positive_lower_bound_fraction": 0.02},
        {"concrete_false_positive_lower_bound_fraction": 0.01},
    )
    result = qv_effect(assessment)
    assert "retain the core result" in result
    assert "1.02" in result
    assert assessment["classification"] == "paired_validation_supports_core_with_measured_sensitivity"


def test_cross_technology_majority_contradiction_is_concrete_failure() -> None:
    assessment = p07_assessment(
        {"qv": 35.9},
        {"pi_ratio_read_over_assembly": 0.44},
        {"lambda_pearson_correlation": -0.14},
        {"concrete_false_positive_lower_bound_fraction": 0.55},
        {"concrete_false_positive_lower_bound_fraction": 0.61},
    )
    assert assessment["classification"] == "concrete_haplotype_reconstruction_failure"
    assert "preserve the core artifact" in assessment["downstream_action"]
    assert "do not use" in assessment["downstream_action"]


def test_caller_disagreement_without_direct_majority_is_not_deletion_trigger() -> None:
    assessment = p07_assessment(
        {"qv": 35.9},
        {"pi_ratio_read_over_assembly": 0.44},
        {"lambda_pearson_correlation": -0.14},
        {"concrete_false_positive_lower_bound_fraction": 0.20},
        {"concrete_false_positive_lower_bound_fraction": 0.10},
    )
    assert assessment["classification"] == "material_method_discordance_without_concrete_failure"
    assert "retain the core result with low validation confidence" in assessment["downstream_action"]
