#!/usr/bin/env python3
"""Collect authoritative login-node sacct telemetry for the real VGP pilot."""

from __future__ import annotations

import argparse
import csv
import subprocess
import time
from pathlib import Path
from typing import Sequence


SACCT_FIELDS = (
    "JobIDRaw", "JobName", "State", "Elapsed", "ElapsedRaw", "Timelimit",
    "AllocCPUS", "ReqMem", "MaxRSS", "CPUTimeRAW", "TotalCPU", "ExitCode",
    "NodeList", "Start", "End", "Submit",
)


class TelemetryError(RuntimeError):
    """Authoritative scheduler accounting could not be closed."""


def submitted_job_ids(path: Path) -> set[str]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows=list(csv.DictReader(handle, delimiter="\t"))
    result={row["job_id"] for row in rows if row.get("job_id")}
    if not result:
        raise TelemetryError("submission manifest has no Slurm jobs")
    return result


def query_sacct(job_ids: Sequence[str], attempts: int = 4) -> list[dict[str, str]]:
    command = [
        "sacct", "-X", "-n", "-P", "-j", ",".join(job_ids),
        "--format=" + ",".join(SACCT_FIELDS),
    ]
    error = ""
    for attempt in range(attempts):
        run = subprocess.run(command, text=True, capture_output=True)
        if run.returncode == 0:
            rows=[]
            for raw in run.stdout.splitlines():
                fields=raw.split("|")
                if fields and fields[-1] == "": fields.pop()
                if len(fields) != len(SACCT_FIELDS):
                    raise TelemetryError(f"unexpected sacct field count: {len(fields)}")
                rows.append(dict(zip(SACCT_FIELDS, fields)))
            return rows
        error=run.stderr.strip()
        if attempt + 1 < attempts: time.sleep(2 ** attempt)
    raise TelemetryError(f"sacct failed after {attempts} attempts: {error}")


def collect(job_ids: set[str]) -> list[dict[str, str]]:
    rows=[]
    ordered=sorted(job_ids, key=int)
    for start in range(0, len(ordered), 80):
        rows.extend(query_sacct(ordered[start:start + 80]))
    # A requested array parent can expand to its parent plus every array task.
    # Retain those task allocations: reducing back to only the requested IDs
    # would discard the scheduler evidence for the 200 biological bootstraps.
    by_id={row["JobIDRaw"]:row for row in rows}
    missing=job_ids-set(by_id)
    if missing:
        raise TelemetryError(f"sacct lacks submitted jobs: {sorted(missing, key=int)}")
    return list(by_id.values())


def write(args: argparse.Namespace) -> int:
    ids=submitted_job_ids(args.submissions)
    ids.update(args.extra_job_id)
    rows=collect(ids)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    partial=args.output.with_suffix(args.output.suffix+".partial")
    fields=("authorization_id","canonical_vgp_root",*SACCT_FIELDS)
    with partial.open("w",newline="",encoding="utf-8") as handle:
        writer=csv.DictWriter(handle,fieldnames=fields,delimiter="\t",lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({
                "authorization_id":"vgp10-auth-20260718-v2",
                "canonical_vgp_root":args.canonical_root,
                **row,
            })
    partial.replace(args.output)
    complete=sum(row["State"].startswith("COMPLETED") for row in rows)
    print(f"rows={len(rows)} completed={complete} output={args.output}")
    return 0


def parser() -> argparse.ArgumentParser:
    value=argparse.ArgumentParser(description=__doc__)
    value.add_argument("--submissions",type=Path,required=True)
    value.add_argument("--output",type=Path,required=True)
    value.add_argument("--canonical-root",default="/moosefs/erikg/vgp")
    value.add_argument("--extra-job-id",action="append",default=[])
    return value


def main(argv: Sequence[str] | None = None) -> int:
    args=parser().parse_args(argv)
    try:
        return write(args)
    except (OSError,ValueError,TelemetryError) as error:
        print(f"ERROR: {error}")
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
