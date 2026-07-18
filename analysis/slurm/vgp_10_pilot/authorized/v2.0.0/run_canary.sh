#!/usr/bin/env bash
set -euo pipefail

PACKET_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
REPO_ROOT=$(CDPATH= cd -- "$PACKET_DIR/../../../../.." && pwd)
AUTHORIZATION="$REPO_ROOT/analysis/vgp_pilot_authorization_v2.json"
SELECTION_ID=P07
AUTHORIZATION_ID=vgp10-auth-20260718-v2
: "${SLURM_TMPDIR:?the authorized canary requires node-local SLURM_TMPDIR}"
: "${VGP_DURABLE_ROOT:=/moosefs/erikg/vgp}"
: "${VGP_ENVIRONMENT_CAPTURE:=$REPO_ROOT/analysis/guix/vgp_10_pilot/realization.json}"
: "${VGP_RETRY:=0}"

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
}
partial=path.with_suffix(".json.partial")
partial.write_text(json.dumps(value,sort_keys=True)+"\n")
partial.replace(path)
PY
}
trap checkpoint EXIT
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

for replicate in $(seq 0 200); do
    export SLURM_ARRAY_TASK_ID=$replicate
    bash "$REPO_ROOT/analysis/slurm/vgp_10_pilot/psmc_array.sh" "$SELECTION_ID"
    if (( replicate % 10 == 0 )); then checkpoint; fi
done
unset SLURM_ARRAY_TASK_ID
bash "$REPO_ROOT/analysis/slurm/vgp_10_pilot/psmc_finalize.sh" "$SELECTION_ID"
checkpoint
