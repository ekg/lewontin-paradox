import csv
from pathlib import Path

import pytest

from analysis.tier3_common import Tier3ValidationError
from analysis.tier3b_population_recovery import (
    DIVERSITY_FIELDS,
    _extract_subset_vcf,
    _uncertainty_fields,
    _write_tsv,
    scheduler_telemetry,
    validate,
)


def _manifest(path: Path) -> Path:
    fields = ["tuple_id", "status", "biological"]
    rows = [
        {"tuple_id": "biological_a", "status": "approved", "biological": "true"},
        {"tuple_id": "biological_b", "status": "approved", "biological": "true"},
    ]
    _write_tsv(path, fields, rows)
    return path


def _diversity(path: Path, numerator: str = "2.5", uncertainty_method: str = "delete_one_sampling_unit_jackknife") -> Path:
    rows = []
    for tuple_id in ("biological_a", "biological_b"):
        for statistic in ("population_pi", "pi_S", "pi_W", "pi_S_over_pi_W"):
            bootstrap = statistic in ("population_pi", "pi_S_over_pi_W")
            rows.append(
                {
                    "tuple_id": tuple_id,
                    "biological": "true",
                    "statistic": statistic,
                    "variant_count": "5",
                    "numerator": numerator,
                    "callable_site_denominator": "100",
                    "estimate": "0.025",
                    "uncertainty_method": (
                        "chromosome_stratified_block_bootstrap" if bootstrap else uncertainty_method
                    ),
                    "uncertainty_unit": (
                        "1-Mb genomic block resampled within chromosome"
                        if bootstrap
                        else "one selected biological individual"
                    ),
                    "uncertainty_replicates": "10000" if bootstrap else "20",
                    "uncertainty_standard_error": "" if bootstrap else "0.001",
                    "interval_low": "0.023",
                    "interval_high": "0.027",
                    "interval_type": (
                        "percentile_95_percent" if bootstrap else "normal_approximation_95_percent"
                    ),
                    "frozen_genomic_bootstrap_status": (
                        "available:21_eligible_1Mb_blocks"
                        if bootstrap
                        else "not_applicable:class_component_uses_sampling_unit_jackknife"
                    ),
                }
            )
    _write_tsv(path, DIVERSITY_FIELDS, rows)
    return path


def _independent(path: Path) -> Path:
    _write_tsv(
        path,
        ("tuple_id", "status"),
        [
            {"tuple_id": "biological_a", "status": "PASS"},
            {"tuple_id": "biological_b", "status": "PASS"},
        ],
    )
    return path


def test_final_assertion_requires_two_complete_biological_tuples(tmp_path, capsys):
    validate(
        _diversity(tmp_path / "diversity.tsv"),
        _manifest(tmp_path / "manifest.tsv"),
        _independent(tmp_path / "independent.tsv"),
    )
    assert "PASS: 2 biological tuples, 8 estimates" in capsys.readouterr().out


@pytest.mark.parametrize(
    "numerator,uncertainty,error",
    [
        ("0", "delete_one_sampling_unit_jackknife", "non-positive numerator"),
        ("2.5", "", "incomplete uncertainty"),
    ],
)
def test_final_assertion_fails_nonpositive_or_missing_uncertainty(
    tmp_path, numerator, uncertainty, error
):
    with pytest.raises(Tier3ValidationError, match=error):
        validate(
            _diversity(tmp_path / "diversity.tsv", numerator, uncertainty),
            _manifest(tmp_path / "manifest.tsv"),
            _independent(tmp_path / "independent.tsv"),
        )


def test_uncertainty_fields_prefer_powered_genomic_bootstrap():
    fields = _uncertainty_fields(
        {
            "bootstrap": {
                "interval": [0.8, 1.2],
                "replicates": 10000,
                "eligible_blocks": 21,
            },
            "uncertainty": {
                "method": "delete_one_sampling_unit_jackknife",
                "unit": "one selected biological individual",
                "replicates": 20,
                "standard_error": 0.1,
                "interval": [0.7, 1.3],
                "interval_type": "normal_approximation_95_percent",
            },
        }
    )
    assert fields == {
        "uncertainty_method": "chromosome_stratified_block_bootstrap",
        "uncertainty_unit": "1-Mb genomic block resampled within chromosome",
        "uncertainty_replicates": 10000,
        "uncertainty_standard_error": "",
        "interval_low": 0.8,
        "interval_high": 1.2,
        "interval_type": "percentile_95_percent",
        "frozen_genomic_bootstrap_status": "available:21_eligible_1Mb_blocks",
    }


def test_independent_subset_extraction_uses_exact_indexed_coordinates(tmp_path, monkeypatch):
    source = tmp_path / "source.vcf.gz"
    source.write_bytes(b"source")
    destination = tmp_path / "subset.vcf.gz"
    observed = {}

    def fake_run(command, check):
        observed["command"] = command
        observed["check"] = check
        destination.write_bytes(b"subset")

    monkeypatch.setattr("analysis.tier3b_population_recovery.subprocess.run", fake_run)
    _extract_subset_vcf(source, destination, "3R", 9_999_999, 10_099_999)
    assert observed == {
        "command": [
            "bcftools",
            "view",
            "--no-version",
            "--regions",
            "3R:10000000-10099999",
            "--output-type",
            "z",
            "--output",
            str(destination),
            str(source),
        ],
        "check": True,
    }


def test_scheduler_telemetry_maps_array_display_ids_and_batch_rss(tmp_path):
    manifest = _manifest(tmp_path / "manifest.tsv")
    sacct = tmp_path / "sacct.psv"
    sacct.write_text(
        "JobID|JobIDRaw|State|Elapsed|TotalCPU|MaxRSS|ReqMem|AllocCPUS|ExitCode\n"
        "42_0|43|COMPLETED|01:02:03|01:01:00||12G|2|0:0\n"
        "42_0.batch|43.batch|COMPLETED|01:02:03|01:01:00|8123M||2|0:0\n"
        "42_1|42|COMPLETED|01:03:04|01:02:00||12G|2|0:0\n"
        "42_1.batch|42.batch|COMPLETED|01:03:04|01:02:00||||0:0\n",
        encoding="utf-8",
    )
    output = tmp_path / "telemetry.tsv"
    scheduler_telemetry(manifest, sacct, "42", output)
    with output.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["tuple_id"] for row in rows] == ["biological_a", "biological_b"]
    assert [row["job_id"] for row in rows] == ["42_0", "42_1"]
    assert rows[0]["max_rss"] == "8123M"
    assert rows[1]["max_rss"] == "not_reported_by_sacct"
