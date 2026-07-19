#!/usr/bin/env bash
set -euo pipefail

PACKET_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
REPO_ROOT=$(CDPATH= cd -- "$PACKET_DIR/../../../../.." && pwd)
AUTHORIZATION="$REPO_ROOT/analysis/vgp_pilot_authorization_v2.json"
SELECTION_ID=P07
AUTHORIZATION_ID=vgp10-auth-20260718-v2
: "${SLURM_JOB_ID:?the authorized canary must execute inside Slurm}"
: "${VGP_DURABLE_ROOT:=/moosefs/erikg/vgp}"
: "${VGP_ENVIRONMENT_CAPTURE:=$REPO_ROOT/analysis/guix/vgp_10_pilot/realization.json}"
: "${VGP_RETRY:=0}"

# This cluster declares Slurm's TmpFS as /tmp but does not synthesize a
# per-job SLURM_TMPDIR.  Create a private job directory there when needed and
# reject network filesystems so the high-I/O working set remains node-local.
: "${VGP_NODE_LOCAL_BASE:=/tmp}"
[[ -d $VGP_NODE_LOCAL_BASE && $VGP_NODE_LOCAL_BASE = /* ]] || {
    echo "invalid node-local scratch base: $VGP_NODE_LOCAL_BASE" >&2; exit 2;
}
node_local_filesystem_type=$(stat -f -c %T -- "$VGP_NODE_LOCAL_BASE")
case "$node_local_filesystem_type" in
    nfs|nfs4|fuse*|lustre|gpfs|ceph)
        echo "refusing non-local scratch filesystem: $node_local_filesystem_type" >&2
        exit 2
        ;;
esac
SLURM_TMPDIR_CREATED=0
if [[ -z ${SLURM_TMPDIR:-} ]]; then
    SLURM_TMPDIR=$(mktemp -d -- "${VGP_NODE_LOCAL_BASE%/}/vgp-slurm-${SLURM_JOB_ID}-XXXXXX")
    chmod 700 "$SLURM_TMPDIR"
    SLURM_TMPDIR_CREATED=1
fi
[[ -d $SLURM_TMPDIR && $SLURM_TMPDIR = "${VGP_NODE_LOCAL_BASE%/}"/* ]] || {
    echo "SLURM_TMPDIR is not inside the verified node-local base: $SLURM_TMPDIR" >&2
    exit 2
}
export SLURM_TMPDIR VGP_NODE_LOCAL_BASE VGP_DURABLE_ROOT node_local_filesystem_type SLURM_TMPDIR_CREATED
export TMPDIR="$SLURM_TMPDIR"

WORK_ROOT="$SLURM_TMPDIR/$AUTHORIZATION_ID-$SELECTION_ID"
LOCAL_DATA="$WORK_ROOT/data"
CHECKPOINT="$VGP_DURABLE_ROOT/pilot/authorized-checkpoints/$AUTHORIZATION_ID/$SELECTION_ID"
RUN_ID="$AUTHORIZATION_ID-canary"
mkdir -p "$LOCAL_DATA/pilot/runs/$RUN_ID" "$CHECKPOINT"

checkpoint() {
    if [[ -d $LOCAL_DATA/pilot/runs/$RUN_ID/$SELECTION_ID ]]; then
        mkdir -p "$CHECKPOINT/run"
        cp -a "$LOCAL_DATA/pilot/runs/$RUN_ID/$SELECTION_ID/." "$CHECKPOINT/run/"
    fi
    python3 - "$CHECKPOINT/checkpoint.json" <<'PY'
import json,os,sys,time
from pathlib import Path
path=Path(sys.argv[1])
value={
  "authorization_id":"vgp10-auth-20260718-v2", "selection_id":"P07",
  "slurm_job_id":os.environ.get("SLURM_JOB_ID"), "retry":int(os.environ.get("VGP_RETRY","0")),
  "checkpoint_epoch":int(time.time()), "node_local_scratch_staged":True,
  "node_local_filesystem_type":os.environ["node_local_filesystem_type"],
  "canonical_vgp_root":os.environ["VGP_DURABLE_ROOT"],
}
partial=path.with_suffix(".json.partial")
partial.write_text(json.dumps(value,sort_keys=True)+"\n")
partial.replace(path)
PY
}
finish() {
    local scratch=${SLURM_TMPDIR:-}
    checkpoint
    if [[ $SLURM_TMPDIR_CREATED == 1 ]]; then
        [[ $scratch == "${VGP_NODE_LOCAL_BASE%/}/vgp-slurm-${SLURM_JOB_ID}-"* ]] || {
            echo "refusing unsafe scratch cleanup target: $scratch" >&2
            return 2
        }
        rm -rf -- "$scratch"
    fi
}
trap finish EXIT
trap 'checkpoint; exit 99' USR1 TERM

if [[ -d $CHECKPOINT/run ]]; then
    mkdir -p "$LOCAL_DATA/pilot/runs/$RUN_ID/$SELECTION_ID"
    cp -a "$CHECKPOINT/run/." "$LOCAL_DATA/pilot/runs/$RUN_ID/$SELECTION_ID/"
fi

export PYTHONPATH="$REPO_ROOT"
python3 -m analysis.vgp_pilot_authorization validate --authorization "$AUTHORIZATION"
python3 -m analysis.vgp_pilot_authorization materialize-input \
    --authorization "$AUTHORIZATION" --selection-id "$SELECTION_ID" --data-root "$LOCAL_DATA"

export VGP_DATA_ROOT="$LOCAL_DATA"
export VGP_RUN_ID="$RUN_ID"
export VGP_ENVIRONMENT_CAPTURE VGP_RETRY

for stage in preflight mapping impg variants consensus; do
    bash "$REPO_ROOT/analysis/slurm/vgp_10_pilot/pair_stage.sh" "$stage" "$SELECTION_ID"
    checkpoint
done

: "${P07_CANARY_PSMC_PARALLELISM:=20}"
[[ $P07_CANARY_PSMC_PARALLELISM =~ ^[1-9][0-9]*$ ]] || {
    echo "invalid PSMC parallelism: $P07_CANARY_PSMC_PARALLELISM" >&2; exit 2;
}
(( P07_CANARY_PSMC_PARALLELISM <= ${SLURM_CPUS_PER_TASK:-1} )) || {
    echo "PSMC parallelism exceeds the Slurm CPU allocation" >&2; exit 2;
}
for batch_start in $(seq 0 "$P07_CANARY_PSMC_PARALLELISM" 200); do
    batch_end=$((batch_start + P07_CANARY_PSMC_PARALLELISM - 1))
    (( batch_end > 200 )) && batch_end=200
    declare -a psmc_pids=()
    for replicate in $(seq "$batch_start" "$batch_end"); do
        (
            export SLURM_ARRAY_TASK_ID=$replicate
            bash "$REPO_ROOT/analysis/slurm/vgp_10_pilot/psmc_array.sh" "$SELECTION_ID"
        ) &
        psmc_pids+=("$!")
    done
    psmc_failed=0
    for pid in "${psmc_pids[@]}"; do
        if ! wait "$pid"; then psmc_failed=1; fi
    done
    (( psmc_failed == 0 )) || { echo "one or more PSMC replicates failed" >&2; exit 2; }
    checkpoint
done
bash "$REPO_ROOT/analysis/slurm/vgp_10_pilot/psmc_finalize.sh" "$SELECTION_ID"
checkpoint
