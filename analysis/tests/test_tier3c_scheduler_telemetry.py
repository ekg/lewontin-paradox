import json
from pathlib import Path

from analysis.tier3c_scheduler_telemetry import merge_telemetry, parse_sacct


ROOT = Path(__file__).resolve().parents[2]


def test_tier3c_standard_and_outlier_slurm_profiles_are_frozen():
    standard = (ROOT / "analysis/slurm/tier3c_run_batch_array.sh").read_text(encoding="utf-8")
    audit = (ROOT / "analysis/slurm/tier3c_control_audit_array.sh").read_text(encoding="utf-8")
    retry = (ROOT / "analysis/slurm/tier3c_run_outlier_retry.sh").read_text(encoding="utf-8")
    for text in (standard, audit):
        assert "#SBATCH --cpus-per-task=2" in text
        assert "#SBATCH --mem=32G" in text
        assert "#SBATCH --time=02:00:00" in text
    assert "#SBATCH --array=0-0%8" in standard
    assert "#SBATCH --cpus-per-task=2" in retry
    assert "#SBATCH --mem=64G" in retry
    assert "#SBATCH --time=04:00:00" in retry
    assert "#SBATCH --array=0-0%1" in retry


def test_sacct_telemetry_merges_scheduler_fields_with_measured_maxrss(tmp_path: Path):
    sidecar = tmp_path / "control.job.json"
    sidecar.write_text(
        json.dumps(
            {
                "dataset_id": "control.tier3c",
                "purpose": "independent_control_audit",
                "slurm_job_id": "101",
                "slurm_array_job_id": "100",
                "slurm_array_task_id": "1",
                "max_rss_kib": 12345,
                "wall_seconds": 2.5,
            }
        ),
        encoding="utf-8",
    )
    sacct = parse_sacct(
        "101|32G||00:00:03|3|COMPLETED|0:0|2|node1|02:00:00\n"
        "101.batch||12345K|00:00:03|3|COMPLETED|0:0|2|node1|\n"
    )
    rows = merge_telemetry([], [sidecar], sacct)
    assert rows == [
        {
            "dataset_id": "control.tier3c",
            "purpose": "independent_control_audit",
            "slurm_job_id": "101",
            "slurm_array_job_id": "100",
            "slurm_array_task_id": "1",
            "requested_cpus": "2",
            "requested_memory": "32G",
            "time_limit": "02:00:00",
            "max_rss_kib": "12345",
            "max_rss_source": "process_getrusage_and_sacct_agree",
            "sacct_max_rss": "12345K",
            "elapsed": "00:00:03",
            "elapsed_seconds": "3",
            "state": "COMPLETED",
            "exit_code": "0:0",
            "node": "node1",
        }
    ]


def test_sacct_parser_preserves_blank_maxrss_for_cluster_without_gather_plugin():
    records = parse_sacct(
        "202|12G||00:10:31|631|COMPLETED|0:0|2|node2|01:00:00\n"
        "202.batch|||00:10:31|631|COMPLETED|0:0|2|node2|\n"
    )
    assert records["202"]["max_rss"] == ""


def test_production_sidecar_is_labeled_as_selective_checksum_rerun(tmp_path: Path):
    sidecar = tmp_path / "production.job.json"
    sidecar.write_text(
        json.dumps(
            {
                "schema_version": "tier3c-job-v1",
                "dataset_id": "control.tier3c",
                "slurm_job_id": "303",
                "max_rss_kib": 42,
            }
        ),
        encoding="utf-8",
    )
    sacct = parse_sacct(
        "303|32G||00:00:01|1|COMPLETED|0:0|2|node3|02:00:00\n"
        "303.batch|||00:00:01|1|COMPLETED|0:0|2|node3|\n"
    )
    rows = merge_telemetry([], [sidecar], sacct)
    assert rows[0]["purpose"] == "selective_control_checksum_rerun"
