import json
from pathlib import Path

from analysis.run_vgp_clean_canary import _nearest_quantile, _strata


ROOT = Path(__file__).resolve().parents[2]


def test_selection_is_predeclared_small_exact_p07():
    value = json.loads((ROOT / "analysis/vgp_clean_canary_selection_v1.json").read_text())
    assert value["selection_id"] == "P07"
    assert value["annotation"]["binding_status"] == "EXACT_DICTIONARY"
    assert value["comparison_target"]["callable_bp_absolute_tolerance"] == 0
    assert value["comparison_target"]["pi_absolute_tolerance"] == 1e-12
    assert all(row["derived_bgzf_path"].endswith(".fa.gz") for row in value["immutable_bgzf_inputs"].values())


def test_strata_reconstruct_early_middle_late_without_crossing_contigs():
    dictionary = [
        {"name": "a", "length": 7_000_000},
        {"name": "b", "length": 8_000_000},
        {"name": "c", "length": 9_000_000},
    ]
    rows = _strata(dictionary)
    assert [row["stratum"] for row in rows] == ["early", "middle", "late"]
    assert rows[0] == {"stratum": "early", "contig": "a", "start": 0, "end": 5_000_000}
    assert rows[1]["contig"] == "b"
    assert rows[2]["contig"] == "c"
    assert all(0 <= row["start"] < row["end"] for row in rows)


def test_bootstrap_centering_uses_predeclared_nearest_index_quantiles():
    values = [float(value) for value in range(200)]
    assert _nearest_quantile(values, 0.025) == 5.0
    assert _nearest_quantile(values, 0.5) == 100.0
    assert _nearest_quantile(values, 0.975) == 194.0
