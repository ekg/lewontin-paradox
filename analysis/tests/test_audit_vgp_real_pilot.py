import pytest

from analysis.audit_vgp_real_pilot import (
    PAIR_IDS,
    PilotAuditError,
    apply_terminal_failures,
    validate_closed_world,
)


def row(pair, complete=False):
    value = {
        "selection_id": pair,
        "status": "complete" if complete else "submitted_incomplete",
        "submitted_stages": ["mapping", "impg", "variants", "consensus", "psmc", "psmc_finalize"],
    }
    if complete:
        value.update({
            "diversity": {"pi": 0.001, "heterozygous_snps": 10},
            "callability": {"final_callable_bp": 10000},
            "psmc": {"finite_bootstraps": 200, "unscaled_intervals": 64,
                     "scenario_rows": 576, "unscaled_primary_preserved": True},
        })
    return value


def test_closed_world_requires_a_new_positive_biological_completion():
    rows = [row(pair, pair == "P07") for pair in PAIR_IDS]
    with pytest.raises(PilotAuditError, match="no newly completed"):
        validate_closed_world(rows)
    rows[7] = row("P08", True)
    validate_closed_world(rows)


def test_closed_world_rejects_zero_estimate_and_missing_submission_stage():
    rows = [row(pair, pair in {"P07", "P08"}) for pair in PAIR_IDS]
    rows[7]["diversity"]["pi"] = 0
    with pytest.raises(PilotAuditError, match="zero diversity"):
        validate_closed_world(rows)
    rows[7] = row("P08", True)
    rows[0]["submitted_stages"].remove("impg")
    with pytest.raises(PilotAuditError, match="lacks real biological submissions"):
        validate_closed_world(rows)


def test_closed_world_attaches_only_retained_terminal_primary_failures():
    rows = [row(pair, pair == "P07") for pair in PAIR_IDS]
    ledger = {"failures": [{
        "selection_id": "P03",
        "classification": "reproducible_hard_primary_execution_failure",
        "primary_preserved": True,
        "alternate_activated": False,
        "reason": "reproduced on independent healthy nodes",
    }, {
        "selection_id": "P04",
        "classification": "intermediate_scratch_scalability_failure",
        "primary_preserved": True,
        "alternate_activated": False,
    }]}
    apply_terminal_failures(rows, ledger)
    assert rows[2]["status"] == "failed_primary"
    assert rows[2]["terminal_failure"]["reason"].startswith("reproduced")
    assert rows[3]["status"] == "submitted_incomplete"


def test_closed_world_does_not_terminalize_reclassified_fastga_infrastructure_failure():
    rows = [row(pair, pair == "P07") for pair in PAIR_IDS]
    ledger = {"failures": [{
        "selection_id": "P03",
        "classification": "reproducible_hard_primary_execution_failure",
        "prior_reason_retained": "FastGA child exited during pair-file writes",
        "terminal_primary_failure": False,
        "reclassified_as": "fastga_scratch_workdir_infrastructure_failure",
        "primary_preserved": True,
        "alternate_activated": False,
    }]}
    apply_terminal_failures(rows, ledger)
    assert rows[2]["status"] == "submitted_incomplete"
