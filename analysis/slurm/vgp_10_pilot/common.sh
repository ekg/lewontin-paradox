#!/usr/bin/env bash
set -euo pipefail

VGP_REPO_ROOT=$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../../.." && pwd)
readonly VGP_REPO_ROOT

fail() {
    echo "ERROR: $*" >&2
    exit 2
}

require_selection_id() {
    [[ ${1:-} =~ ^P(0[1-9]|10)$ ]] || fail "selection id must be P01..P10"
}

require_runtime() {
    : "${VGP_DATA_ROOT:?VGP_DATA_ROOT is required}"
    : "${VGP_RUN_ID:?VGP_RUN_ID is required}"
    : "${VGP_ENVIRONMENT_CAPTURE:?VGP_ENVIRONMENT_CAPTURE is required}"
    [[ $VGP_DATA_ROOT = /* ]] || fail "VGP_DATA_ROOT must be absolute"
    [[ $VGP_RUN_ID =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || fail "unsafe run id"
    [[ -f $VGP_ENVIRONMENT_CAPTURE ]] || fail "environment capture is absent"
    VGP_SACCT=$(command -v sacct || true)
    export VGP_SACCT
    export PYTHONPATH="$VGP_REPO_ROOT"
    export VGP_PROFILE
    VGP_PROFILE=$(python3 - "$VGP_ENVIRONMENT_CAPTURE" <<'PY'
import json,sys
value=json.load(open(sys.argv[1]))
required=("profile","derivation","closure_sha256","executables")
if any(not value.get(key) for key in required):
    raise SystemExit("uncaptured Guix environment")
print(value["profile"])
PY
)
    [[ $VGP_PROFILE == /gnu/store/* ]] || fail "profile is not a Guix store path"
    VGP_CAPTURED_TOOL_PATH=$(python3 - "$VGP_ENVIRONMENT_CAPTURE" <<'PY'
import json,sys
from pathlib import Path
rows=json.load(open(sys.argv[1]))["executables"]
print(":".join(sorted({str(Path(row["path"]).parent) for row in rows})))
PY
)
    export PATH="$VGP_PROFILE/bin:$VGP_CAPTURED_TOOL_PATH"
    python3 -m analysis.vgp_10_pilot verify-capture "$VGP_ENVIRONMENT_CAPTURE" >/dev/null
    unset PYTHONHOME CONDA_PREFIX VIRTUAL_ENV LD_LIBRARY_PATH
}

tool_path() {
    python3 - "$VGP_ENVIRONMENT_CAPTURE" "$1" <<'PY'
import json,sys
rows=json.load(open(sys.argv[1]))["executables"]
matches=[row for row in rows if row["name"] == sys.argv[2]]
if len(matches) != 1:
    raise SystemExit(f"tool identity does not resolve once: {sys.argv[2]}")
print(matches[0]["path"])
PY
}

verify_tool() {
    python3 - "$VGP_ENVIRONMENT_CAPTURE" "$1" <<'PY'
import hashlib,json,sys
from pathlib import Path
rows=json.load(open(sys.argv[1]))["executables"]
matches=[row for row in rows if row["name"] == sys.argv[2]]
if len(matches) != 1:
    raise SystemExit("tool missing from capture")
row=matches[0]
path=Path(row["path"])
if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != row["sha256"]:
    raise SystemExit(f"tool digest mismatch: {sys.argv[2]}")
PY
}

begin_stage() {
    local stage=$1
    VGP_STAGE_TAG=${stage//\//-}
    export VGP_STAGE_TAG
    VGP_PAIR_RUN="$VGP_DATA_ROOT/pilot/runs/$VGP_RUN_ID/$VGP_SELECTION_ID"
    VGP_STAGE_FINAL="$VGP_PAIR_RUN/$stage"
    if [[ -f $VGP_STAGE_FINAL/.complete.json ]]; then
        echo "RESUME: $VGP_SELECTION_ID $stage already complete"
        exit 0
    fi
    [[ ! -e $VGP_STAGE_FINAL ]] || fail "final stage exists without sentinel: $VGP_STAGE_FINAL"
    local scratch_root=${SLURM_TMPDIR:-${VGP_SCRATCH_ROOT:?VGP_SCRATCH_ROOT is required outside Slurm}}
    VGP_STAGE_PARTIAL="$scratch_root/vgp-$VGP_RUN_ID-$VGP_SELECTION_ID-$VGP_STAGE_TAG-${SLURM_JOB_ID:-local}.partial"
    [[ $VGP_STAGE_PARTIAL = "$scratch_root"/* ]] || fail "unsafe scratch stage path"
    if [[ -e $VGP_STAGE_PARTIAL ]]; then
        mv "$VGP_STAGE_PARTIAL" "$VGP_STAGE_PARTIAL.stale.$(date -u +%Y%m%dT%H%M%SZ)"
    fi
    mkdir -p "$VGP_STAGE_PARTIAL" "$VGP_PAIR_RUN/telemetry"
    export VGP_PAIR_RUN VGP_STAGE_FINAL VGP_STAGE_PARTIAL
    VGP_STAGE_STARTED_EPOCH=$(date +%s)
    export VGP_STAGE_STARTED_EPOCH
    trap 'record_telemetry failure' ERR INT TERM
}

record_telemetry() {
    local disposition=${1:-success}
    local output="$VGP_PAIR_RUN/telemetry/${VGP_STAGE_TAG}.${SLURM_JOB_ID:-local}.json"
    python3 - "$output" "$disposition" <<'PY'
import json,os,resource,sys,time
from pathlib import Path
usage=resource.getrusage(resource.RUSAGE_CHILDREN)
scratch=Path(os.environ["VGP_STAGE_PARTIAL"])
files=[path for path in scratch.rglob("*") if path.is_file()] if scratch.exists() else []
value={
 "selection_id":os.environ["VGP_SELECTION_ID"], "stage":os.environ["VGP_STAGE_NAME"],
 "job_id":os.environ.get("SLURM_JOB_ID","local"), "array_task_id":os.environ.get("SLURM_ARRAY_TASK_ID"),
 "started_epoch":int(os.environ["VGP_STAGE_STARTED_EPOCH"]), "ended_epoch":int(time.time()),
 "disposition":sys.argv[2], "maximum_rss_kib":usage.ru_maxrss,
 "child_cpu_seconds":usage.ru_utime+usage.ru_stime,
 "filesystem_read_bytes":usage.ru_inblock*512,"filesystem_write_bytes":usage.ru_oublock*512,
 "scratch_high_water_bytes":sum(path.stat().st_size for path in files),
 "scratch_file_count":len(files),
 "scratch_path":os.environ["VGP_STAGE_PARTIAL"], "retry":int(os.environ.get("VGP_RETRY","0")),
}
path=sys.argv[1]
partial=path+".partial"
open(partial,"w").write(json.dumps(value,sort_keys=True)+"\n")
os.replace(partial,path)
PY
    if [[ -n ${SLURM_JOB_ID:-} && -n ${VGP_SACCT:-} && -x $VGP_SACCT ]]; then
        local sacct_output="$VGP_PAIR_RUN/telemetry/${VGP_STAGE_TAG}.${SLURM_JOB_ID}.sacct.tsv"
        "$VGP_SACCT" -j "$SLURM_JOB_ID" --parsable2 --noheader \
            -o JobID,ElapsedRaw,TotalCPU,MaxRSS,MaxDiskRead,MaxDiskWrite,State,ExitCode \
            >"$sacct_output.partial" || true
        mv "$sacct_output.partial" "$sacct_output"
    fi
}

promote_stage() {
    trap - ERR INT TERM
    record_telemetry success
    python3 - "$VGP_STAGE_PARTIAL" "$VGP_STAGE_FINAL" "$VGP_SELECTION_ID" "$VGP_STAGE_NAME" <<'PY'
import sys
from pathlib import Path
from analysis.vgp_10_pilot import atomic_promote
atomic_promote(Path(sys.argv[1]),Path(sys.argv[2]),{
 "selection_id":sys.argv[3],"stage":sys.argv[4],"atomic_promotion":True,
})
PY
}
