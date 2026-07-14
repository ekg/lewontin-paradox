#!/usr/bin/env python3
"""Collect Tier 3c scheduler and process telemetry into one auditable TSV."""

from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import tempfile
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence


FIELDS = (
    "dataset_id",
    "purpose",
    "slurm_job_id",
    "slurm_array_job_id",
    "slurm_array_task_id",
    "requested_cpus",
    "requested_memory",
    "time_limit",
    "max_rss_kib",
    "max_rss_source",
    "sacct_max_rss",
    "elapsed",
    "elapsed_seconds",
    "state",
    "exit_code",
    "node",
)
SACCT_FIELDS = (
    "JobIDRaw",
    "ReqMem",
    "MaxRSS",
    "Elapsed",
    "ElapsedRaw",
    "State",
    "ExitCode",
    "AllocCPUS",
    "NodeList",
    "Timelimit",
)


class TelemetryError(ValueError):
    """Missing, malformed, or failed scheduler telemetry."""


def _rss_kib(text: str) -> Optional[int]:
    if not text:
        return None
    suffixes = {"K": 1, "M": 1024, "G": 1024 * 1024, "T": 1024 * 1024 * 1024}
    unit = text[-1].upper()
    if unit in suffixes:
        return int(float(text[:-1]) * suffixes[unit])
    return int(text)


def parse_sacct(text: str) -> Dict[str, Dict[str, str]]:
    """Parse pipe-delimited sacct rows and fold step MaxRSS into each job."""

    jobs: Dict[str, Dict[str, str]] = {}
    step_rss: Dict[str, List[str]] = {}
    for raw in text.splitlines():
        if not raw.strip():
            continue
        values = raw.rstrip("\n").split("|")
        if len(values) != len(SACCT_FIELDS):
            raise TelemetryError(f"sacct row has {len(values)} fields, expected {len(SACCT_FIELDS)}")
        row = dict(zip(SACCT_FIELDS, values))
        job_raw = row["JobIDRaw"]
        base = job_raw.split(".", 1)[0]
        if "." not in job_raw:
            jobs[base] = {
                "requested_memory": row["ReqMem"],
                "max_rss": row["MaxRSS"],
                "elapsed": row["Elapsed"],
                "elapsed_seconds": row["ElapsedRaw"],
                "state": row["State"],
                "exit_code": row["ExitCode"],
                "requested_cpus": row["AllocCPUS"],
                "node": row["NodeList"],
                "time_limit": row["Timelimit"],
            }
        if row["MaxRSS"]:
            step_rss.setdefault(base, []).append(row["MaxRSS"])
    for job_id, rows in step_rss.items():
        if job_id in jobs:
            jobs[job_id]["max_rss"] = max(rows, key=lambda value: _rss_kib(value) or -1)
    return jobs


def _baseline_jobs(qc_dir: Path) -> List[Dict[str, Any]]:
    jobs: List[Dict[str, Any]] = []
    for path in sorted(qc_dir.glob("*.json")):
        qc = json.loads(path.read_text(encoding="utf-8"))
        job = qc.get("job")
        if not job:
            continue
        jobs.append(
            {
                **job,
                "purpose": "frozen_primary_analysis",
                "slurm_array_job_id": job.get("slurm_array_job_id", ""),
            }
        )
    return jobs


def _read_sidecars(paths: Iterable[Path]) -> List[Dict[str, Any]]:
    records: List[Dict[str, Any]] = []
    for path in paths:
        record = json.loads(path.read_text(encoding="utf-8"))
        if "purpose" not in record:
            if record.get("schema_version") == "tier3c-job-v1":
                record["purpose"] = "selective_control_checksum_rerun"
            else:
                raise TelemetryError(f"sidecar has no purpose: {path}")
        records.append(record)
    return records


def merge_telemetry(
    baseline_jobs: Sequence[Mapping[str, Any]],
    sidecars: Sequence[Path],
    sacct: Mapping[str, Mapping[str, str]],
) -> List[Dict[str, str]]:
    records = list(baseline_jobs) + _read_sidecars(sidecars)
    rows: List[Dict[str, str]] = []
    for record in records:
        job_id = str(record["slurm_job_id"])
        if job_id == "login":
            continue
        if job_id not in sacct:
            raise TelemetryError(f"sacct did not return job {job_id}")
        scheduler = sacct[job_id]
        state = scheduler["state"]
        exit_code = scheduler["exit_code"]
        if state != "COMPLETED" or exit_code != "0:0":
            raise TelemetryError(f"job {job_id} did not complete cleanly: {state} {exit_code}")
        process_rss = int(record.get("max_rss_kib", 0))
        sacct_rss_text = scheduler["max_rss"]
        sacct_rss = _rss_kib(sacct_rss_text)
        if not process_rss and sacct_rss is None:
            raise TelemetryError(f"job {job_id} has no MaxRSS measurement")
        if process_rss and sacct_rss is not None:
            source = (
                "process_getrusage_and_sacct_agree"
                if process_rss == sacct_rss
                else "process_getrusage_and_sacct_both_reported"
            )
        elif process_rss:
            source = "process_getrusage;_sacct_blank_jobacct_gather_disabled"
        else:
            source = "sacct"
        rows.append(
            {
                "dataset_id": str(record["dataset_id"]),
                "purpose": str(record.get("purpose", "unspecified")),
                "slurm_job_id": job_id,
                "slurm_array_job_id": str(record.get("slurm_array_job_id", "")),
                "slurm_array_task_id": str(record.get("slurm_array_task_id", "")),
                "requested_cpus": scheduler["requested_cpus"],
                "requested_memory": scheduler["requested_memory"],
                "time_limit": scheduler["time_limit"],
                "max_rss_kib": str(process_rss or sacct_rss),
                "max_rss_source": source,
                "sacct_max_rss": sacct_rss_text,
                "elapsed": scheduler["elapsed"],
                "elapsed_seconds": scheduler["elapsed_seconds"],
                "state": state,
                "exit_code": exit_code,
                "node": scheduler["node"],
            }
        )
    return sorted(
        rows,
        key=lambda row: (
            row["purpose"].encode("utf-8"),
            row["dataset_id"].encode("utf-8"),
            int(row["slurm_job_id"]),
        ),
    )


def _query_sacct(job_ids: Sequence[str], starttime: str) -> Dict[str, Dict[str, str]]:
    records: Dict[str, Dict[str, str]] = {}
    for offset in range(0, len(job_ids), 100):
        chunk = job_ids[offset : offset + 100]
        completed = subprocess.run(
            [
                "sacct",
                "-n",
                "-P",
                "-j",
                ",".join(chunk),
                "--starttime",
                starttime,
                "--format=" + ",".join(SACCT_FIELDS),
            ],
            check=True,
            stdout=subprocess.PIPE,
            text=True,
        )
        records.update(parse_sacct(completed.stdout))
    return records


def _atomic_tsv(destination: Path, rows: Sequence[Mapping[str, str]]) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    handle, temporary = tempfile.mkstemp(prefix=f".{destination.name}.", dir=destination.parent)
    try:
        with os.fdopen(handle, "w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, destination)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def collect(qc_dir: Path, sidecars: Sequence[Path], destination: Path, starttime: str) -> None:
    baseline = _baseline_jobs(qc_dir)
    extra = _read_sidecars(sidecars)
    job_ids = sorted(
        {
            str(record["slurm_job_id"])
            for record in baseline + extra
            if str(record["slurm_job_id"]) != "login"
        },
        key=int,
    )
    scheduler = _query_sacct(job_ids, starttime)
    _atomic_tsv(destination, merge_telemetry(baseline, sidecars, scheduler))


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("qc_dir", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--job-sidecar", action="append", default=[], type=Path)
    parser.add_argument("--starttime", default="2026-07-13")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _parser().parse_args(argv)
    collect(args.qc_dir, args.job_sidecar, args.output, args.starttime)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (TelemetryError, OSError, subprocess.CalledProcessError, json.JSONDecodeError) as error:
        print(f"tier3c-scheduler-telemetry: {error}", file=os.sys.stderr)
        raise SystemExit(2)
