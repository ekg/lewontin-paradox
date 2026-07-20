import json
from pathlib import Path

from analysis import repair_vgp_psmc_bootstrap as repair


def write_psmc(path: Path, theta: float) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        f"RD\t25\nTR\t{theta}\t1.0\nRS\t0\t0.0\t1.0\n",
        encoding="utf-8",
    )


def fixture_root(tmp_path: Path) -> Path:
    root = tmp_path / "vgp"
    for pair, relative in repair.PAIR_RELATIVE_ROOTS.items():
        run = root / relative
        source = run / "consensus/consensus/input.psmcfa"
        source.parent.mkdir(parents=True, exist_ok=True)
        source.write_text(">chr1\nNNKTTNKT\n>chr2\nTNNK\n", encoding="utf-8")
        write_psmc(run / "psmc/replicate-000/unscaled.psmc", 1.0)
    return root


def test_prepare_pair_freezes_primary_psmcfa_population_and_contig_boundaries(tmp_path: Path):
    root = fixture_root(tmp_path)
    output = tmp_path / "diagnostics"
    result = repair.prepare_pair(root, output, "P07")
    assert result["blocks_cross_contig_boundaries"] is False
    assert result["masked_and_callable_sampling_population_preserved"] is True
    assert result["primary_psmcfa_symbols"] == result["frozen_unit_symbols"]
    assert result["primary_psmcfa_bins"] == result["frozen_unit_bins"] == 12
    assert result["unit_count"] == 2
    manifest = (output / "P07/design/bootstrap_manifest.tsv").read_text()
    assert len(manifest.splitlines()) == 201
    assert "block_bins" in manifest.splitlines()[0]


def test_audit_requires_200_finite_outputs_and_predeclared_theta_centering(tmp_path: Path):
    root = fixture_root(tmp_path)
    output = tmp_path / "diagnostics"
    diagnostic = {
        "schema_version": "vgp-psmc-bootstrap-centering-diagnostic-v1",
        "repair_id": repair.REPAIR_ID,
        "canonical_vgp_root": str(root),
        "diagnostic": dict(repair.CENTERING_DIAGNOSTIC),
    }
    repair.atomic_text(
        output / "centering_diagnostic.predeclared.json", repair.canonical_json(diagnostic)
    )
    for pair in repair.PAIR_RELATIVE_ROOTS:
        repair.prepare_pair(root, output, pair)
        for replicate in range(1, 201):
            theta = 0.9 + 0.2 * (replicate - 1) / 199
            write_psmc(
                output / pair / f"replicate-{replicate:03d}/bootstrap.unscaled.psmc",
                theta,
            )
    result = repair.audit(root, output)
    assert result["passed"] is True
    assert json.loads(
        (output / "centering_diagnostic.predeclared.json").read_text()
    )["diagnostic"]["predeclared_before_execution"] is True
    for pair in result["pairs"].values():
        assert pair["finite_bootstraps"] == pair["bootstrap_attempts"] == 200
        assert pair["centering_diagnostic"]["primary_inside_bounds"] is True
        assert pair["passed"] is True
