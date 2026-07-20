import csv

from analysis import collect_vgp_real_pilot_sacct as telemetry
from analysis.collect_vgp_real_pilot_sacct import SACCT_FIELDS, submitted_job_ids


def test_submission_job_census_is_exact_and_nonempty(tmp_path):
    path=tmp_path/"submissions.tsv"
    with path.open("w",newline="") as handle:
        writer=csv.DictWriter(handle,fieldnames=("job_id","selection_id"),delimiter="\t")
        writer.writeheader()
        writer.writerows([
            {"job_id":"101","selection_id":"P01"},
            {"job_id":"101","selection_id":"P01"},
            {"job_id":"202","selection_id":"P02"},
        ])
    assert submitted_job_ids(path) == {"101","202"}
    assert SACCT_FIELDS[0] == "JobIDRaw" and SACCT_FIELDS[-1] == "Submit"


def test_collect_retains_array_task_allocations(monkeypatch):
    def fake_query(job_ids):
        assert job_ids == ["101", "202"]
        return [
            {"JobIDRaw": "101"},
            {"JobIDRaw": "202"},
            {"JobIDRaw": "202_0"},
            {"JobIDRaw": "202_1"},
        ]

    monkeypatch.setattr(telemetry, "query_sacct", fake_query)
    assert [row["JobIDRaw"] for row in telemetry.collect({"202", "101"})] == [
        "101", "202", "202_0", "202_1",
    ]
