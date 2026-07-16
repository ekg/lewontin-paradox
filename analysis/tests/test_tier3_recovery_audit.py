from pathlib import Path

from analysis import tier3_recovery_audit as audit


def test_recovery_provenance_rejects_superseded_lineage_and_pins_origin_main():
    provenance = audit.validate_provenance()
    assert len(provenance["sweepga_commit"]) == 40
    assert len(provenance["sweepga_binary_sha256"]) == 64
    assert provenance["tier3a_guix_profile"].startswith("/gnu/store/")
    assert provenance["tier3b_guix_profile"].startswith("/gnu/store/")


def test_recovery_ledger_and_independent_numerical_audit():
    ledger = audit.evidence_ledger()
    assert len(ledger) == 23
    assert {row["tier"] for row in ledger} == {"3A", "3B"}
    checks = audit.headline_audit(ledger)
    assert len(checks) == 23
    assert {row["status"] for row in checks} == {"PASS"}


def test_committed_synthesis_and_figure_cover_every_recovered_observation():
    ledger = audit.evidence_ledger()
    audit.verify_synthesis_table(ledger)
    assert (Path(__file__).resolve().parents[2] / "analysis/fig_tier3.png").stat().st_size > 1000


def test_manuscript_reports_recovery_definitions_uncertainty_and_claim_boundary():
    root = Path(__file__).resolve().parents[2]
    text = (root / "manuscript.typ").read_text(encoding="utf-8")
    normalized = " ".join(text.split())
    for required in (
        "three biological H1/H2 comparisons",
        "H1/H2 alternative alleles per",
        "28,233/2,038,234",
        "20 wild diploid individuals",
        "184,914.842",
        "10,000-replicate chromosome-stratified 1-Mb block-bootstrap",
        "not polarized SFS-$B$",
        "not pooled with the assembly comparisons",
    ):
        assert required in normalized
    assert "have zero eligible tuples" not in text
    assert "tables are header-only" not in text
