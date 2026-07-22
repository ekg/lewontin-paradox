#!/usr/bin/env bash
#SBATCH --job-name=vgp-clean-P07
#SBATCH --partition=highmem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=3-00:00:00
set -euo pipefail

: "${SLURM_JOB_ID:?this clean canary must run under Slurm}"
REPO_ROOT=${SLURM_SUBMIT_DIR:?submit from repository root}
SELECTION="$REPO_ROOT/analysis/vgp_clean_canary_selection_v1.json"
CAPTURE="$REPO_ROOT/analysis/guix/vgp_10_pilot/realization.json"
RUN_ID=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["run_id"])' "$SELECTION")
SELECTION_ID=P07
DURABLE_BASE=/moosefs/erikg/vgp/pilot/clean-canary
DURABLE_TARGET="$DURABLE_BASE/$RUN_ID/$SELECTION_ID"
[[ ! -e $DURABLE_TARGET ]] || { echo "clean target already exists: $DURABLE_TARGET" >&2; exit 2; }
scratch=$(mktemp -d -- "/scratch/vgp-clean-$SELECTION_ID-${SLURM_JOB_ID}-XXXXXX")
scratch_resolved=$(readlink -f -- "$scratch")
case $(stat -f -c %T -- "$scratch") in nfs|nfs4|fuse*|lustre|gpfs|ceph) echo "scratch is not node-local" >&2; exit 2;; esac
export SLURM_TMPDIR="$scratch" TMPDIR="$scratch" TMP="$scratch" TEMP="$scratch"
export VGP_NODE_LOCAL_BASE=/scratch VGP_ENVIRONMENT_CAPTURE="$CAPTURE"
export VGP_DATA_ROOT="$scratch/data" VGP_RUN_ID="$RUN_ID" VGP_SELECTION_ID="$SELECTION_ID"
export VGP_FASTGA_GUARD_INTERVAL_SECONDS=2 PYTHONPATH="$REPO_ROOT"
export VGP_TASK_ID=run-vgp-clean-canary
export VGP_MAPPING_SCHEMA_VERSION=vgp-clean-canary-mapping-execution-v1
export VGP_ANNOTATION_SCHEMA_VERSION=vgp-clean-canary-exact-annotation-v1
mkdir -p "$VGP_DATA_ROOT/pilot/runs/$RUN_ID" "$DURABLE_BASE/logs"
started=$(date +%s)

cleanup() {
    status=$?
    if (( status != 0 )); then
        mkdir -p "$DURABLE_BASE/failures/$RUN_ID-$SLURM_JOB_ID"
        cp -a "$VGP_DATA_ROOT/pilot/runs/$RUN_ID/$SELECTION_ID/telemetry" \
            "$DURABLE_BASE/failures/$RUN_ID-$SLURM_JOB_ID/" 2>/dev/null || true
    fi
    [[ $scratch == "/scratch/vgp-clean-$SELECTION_ID-${SLURM_JOB_ID}-"* ]] && rm -rf -- "$scratch"
    exit "$status"
}
trap cleanup EXIT

cd -- "$scratch"
[[ $(readlink -f -- /proc/$$/cwd) == "$scratch_resolved" ]] || exit 70
python3 "$REPO_ROOT/analysis/run_vgp_clean_canary.py" materialize \
    --selection "$SELECTION" --data-root "$VGP_DATA_ROOT" >"$scratch/materialize.json"

export VGP_STAGE_NAME=preflight
bash "$REPO_ROOT/analysis/slurm/vgp_10_pilot/pair_stage.sh" preflight "$SELECTION_ID"
bash "$REPO_ROOT/analysis/slurm/vgp_10_pilot/mapping_stage.sh" "$SELECTION_ID"
for stage in impg variants consensus; do
    bash "$REPO_ROOT/analysis/slurm/vgp_10_pilot/pair_stage.sh" "$stage" "$SELECTION_ID"
done
bash "$REPO_ROOT/analysis/slurm/vgp_10_pilot/annotation_stage.sh" "$SELECTION_ID"

parallelism=20
for batch_start in $(seq 0 "$parallelism" 200); do
    batch_end=$((batch_start + parallelism - 1)); (( batch_end > 200 )) && batch_end=200
    pids=()
    for replicate in $(seq "$batch_start" "$batch_end"); do
        (export SLURM_ARRAY_TASK_ID=$replicate; bash "$REPO_ROOT/analysis/slurm/vgp_10_pilot/psmc_array.sh" "$SELECTION_ID") &
        pids+=("$!")
    done
    failed=0
    for pid in "${pids[@]}"; do wait "$pid" || failed=1; done
    (( failed == 0 )) || { echo "PSMC batch failed" >&2; exit 2; }
done
bash "$REPO_ROOT/analysis/slurm/vgp_10_pilot/psmc_finalize.sh" "$SELECTION_ID"

run_root="$VGP_DATA_ROOT/pilot/runs/$RUN_ID/$SELECTION_ID"
bcftools=$(python3 -c 'import json,sys; print(next(x["path"] for x in json.load(open(sys.argv[1]))["executables"] if x["name"]=="bcftools"))' "$CAPTURE")
python3 "$REPO_ROOT/analysis/run_vgp_clean_canary.py" audit \
    --selection "$SELECTION" --run-root "$run_root" --bcftools "$bcftools" \
    --output "$run_root/execution.json" >"$run_root/audit.stdout.json"
cp "$scratch/materialize.json" "$run_root/materialize.json"
python3 - "$run_root/job.json" "$started" <<'PY'
import json,os,sys,time
from pathlib import Path
Path(sys.argv[1]).write_text(json.dumps({
 "schema_version":"vgp-clean-canary-job-v1","task_id":"run-vgp-clean-canary",
 "job_id":os.environ["SLURM_JOB_ID"],"node":os.environ.get("SLURMD_NODENAME"),
 "started_epoch":int(sys.argv[2]),"ended_epoch":int(time.time()),
 "private_scratch":os.environ["SLURM_TMPDIR"],"TMPDIR":os.environ["TMPDIR"],
 "TMP":os.environ["TMP"],"TEMP":os.environ["TEMP"],"prior_intermediates_reused":False,
},sort_keys=True)+"\n")
PY

mkdir -p "$DURABLE_BASE/$RUN_ID"
staging="$DURABLE_BASE/$RUN_ID/.P07.${SLURM_JOB_ID}.partial"
[[ ! -e $staging ]] || exit 2
cp -a "$run_root" "$staging"
mv -- "$staging" "$DURABLE_TARGET"
printf '%s\n' "$DURABLE_TARGET"
