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

prepare_node_local_scratch() {
    local stage_key=${VGP_STAGE_NAME%%/*}
    [[ $stage_key == psmc ]] && stage_key=psmc
    local resource="$VGP_DATA_ROOT/pilot/inputs/$VGP_SELECTION_ID/resources.json"
    local required_bytes
    required_bytes=$(python3 - "$resource" "$stage_key" "${VGP_RESOURCE_PLAN:-}" \
        "$VGP_SELECTION_ID" <<'PY'
import json,sys
value=json.load(open(sys.argv[1]))
stage=value.get("stages",{}).get(sys.argv[2],{})
if sys.argv[3]:
    plan=json.load(open(sys.argv[3]))
    override=plan.get("pairs",{}).get(sys.argv[4],{}).get("stages",{}).get(sys.argv[2],{})
    stage={**stage,**override}
required=int(stage.get("scratch_bytes_high",0))
if required <= 0: raise SystemExit(f"missing positive scratch estimate for {sys.argv[2]}")
print(required)
PY
)
    VGP_NODE_LOCAL_BASE=${VGP_NODE_LOCAL_BASE:-/scratch}
    [[ -d $VGP_NODE_LOCAL_BASE && -w $VGP_NODE_LOCAL_BASE ]] || \
        fail "node-local scratch base is absent or unwritable: $VGP_NODE_LOCAL_BASE"
    local filesystem_type available_bytes
    filesystem_type=$(stat -f -c %T -- "$VGP_NODE_LOCAL_BASE")
    case "$filesystem_type" in
        nfs|nfs4|fuse*|lustre|gpfs|ceph) fail "scratch is not node-local: $filesystem_type" ;;
    esac
    available_bytes=$(df -PB1 -- "$VGP_NODE_LOCAL_BASE" | awk 'NR==2 {print $4}')
    [[ $available_bytes =~ ^[0-9]+$ ]] || fail "could not measure node-local scratch capacity"
    (( available_bytes >= required_bytes )) || \
        fail "insufficient node-local scratch: $available_bytes < measured $required_bytes bytes"
    if [[ -z ${SLURM_TMPDIR:-} ]]; then
        SLURM_TMPDIR=$(mktemp -d -- \
            "${VGP_NODE_LOCAL_BASE%/}/vgp-${VGP_RUN_ID}-${VGP_SELECTION_ID}-${SLURM_JOB_ID:-local}-XXXXXX")
        VGP_SLURM_TMPDIR_CREATED=1
    else
        VGP_SLURM_TMPDIR_CREATED=0
    fi
    [[ -d $SLURM_TMPDIR && $SLURM_TMPDIR = "${VGP_NODE_LOCAL_BASE%/}"/* ]] || \
        fail "SLURM_TMPDIR is outside the verified node-local base: $SLURM_TMPDIR"
    export SLURM_TMPDIR TMPDIR="$SLURM_TMPDIR" VGP_NODE_LOCAL_BASE VGP_SLURM_TMPDIR_CREATED
    export VGP_SCRATCH_REQUIRED_BYTES=$required_bytes VGP_SCRATCH_AVAILABLE_BYTES=$available_bytes
}

cleanup_stage_scratch() {
    if [[ ${VGP_SLURM_TMPDIR_CREATED:-0} == 1 && -n ${SLURM_TMPDIR:-} ]]; then
        [[ $SLURM_TMPDIR == "${VGP_NODE_LOCAL_BASE%/}/vgp-${VGP_RUN_ID}-${VGP_SELECTION_ID}-${SLURM_JOB_ID:-local}-"* ]] || \
            fail "refusing unsafe scratch cleanup target: $SLURM_TMPDIR"
        rm -rf -- "$SLURM_TMPDIR"
    fi
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
    prepare_node_local_scratch
    local scratch_root=$SLURM_TMPDIR
    VGP_STAGE_PARTIAL="$scratch_root/vgp-$VGP_RUN_ID-$VGP_SELECTION_ID-$VGP_STAGE_TAG-${SLURM_JOB_ID:-local}.partial"
    [[ $VGP_STAGE_PARTIAL = "$scratch_root"/* ]] || fail "unsafe scratch stage path"
    if [[ -e $VGP_STAGE_PARTIAL ]]; then
        mv "$VGP_STAGE_PARTIAL" "$VGP_STAGE_PARTIAL.stale.$(date -u +%Y%m%dT%H%M%SZ)"
    fi
    mkdir -p "$VGP_STAGE_PARTIAL" "$VGP_PAIR_RUN/telemetry"
    export VGP_PAIR_RUN VGP_STAGE_FINAL VGP_STAGE_PARTIAL
    VGP_STAGE_STARTED_EPOCH=$(date +%s)
    export VGP_STAGE_STARTED_EPOCH
    trap 'record_telemetry failure || true; cleanup_stage_scratch || true' ERR INT TERM
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
diagnostics={}
for diagnostic in files:
    if diagnostic.suffix in {".stderr",".stdout",".log"}:
        try:
            diagnostics[str(diagnostic.relative_to(scratch))]=diagnostic.read_text(
                encoding="utf-8",errors="replace")[-32768:]
        except OSError:
            pass
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
 "canonical_vgp_root":os.environ["VGP_DATA_ROOT"],
 "scratch_required_bytes":int(os.environ["VGP_SCRATCH_REQUIRED_BYTES"]),
 "scratch_available_bytes_at_start":int(os.environ["VGP_SCRATCH_AVAILABLE_BYTES"]),
 "node_local_scratch_base":os.environ["VGP_NODE_LOCAL_BASE"],
 "diagnostic_tails":diagnostics,
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
    python3 - "$VGP_STAGE_PARTIAL" "$VGP_STAGE_FINAL" "$VGP_SELECTION_ID" "$VGP_STAGE_NAME" \
        "${SLURM_JOB_ID:-local}" "$VGP_DATA_ROOT" <<'PY'
import sys
from pathlib import Path
from analysis.vgp_10_pilot import promote_stage_cross_filesystem
promote_stage_cross_filesystem(Path(sys.argv[1]),Path(sys.argv[2]),{
 "selection_id":sys.argv[3],"stage":sys.argv[4],"atomic_promotion":True,
 "canonical_vgp_root":sys.argv[6],"slurm_job_id":sys.argv[5],
},sys.argv[5])
PY
    cleanup_stage_scratch
}
