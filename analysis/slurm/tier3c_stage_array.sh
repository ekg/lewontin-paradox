#!/usr/bin/env bash
#SBATCH --job-name=tier3c-stage
#SBATCH --partition=workers
#SBATCH --array=0-0
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --export=NONE
set -euo pipefail

DISCOVERY=${1:?discovery manifest required}
STAGE_ROOT=${2:?stage root required}
TIER3_ROOT=${3:?repository root required}
ENVIRONMENT_RECORD=${4:?environment record required}
export TIER3_ROOT TIER3_SCRATCH_ROOT="$STAGE_ROOT"

exec "$TIER3_ROOT/analysis/slurm/guix_job.sh" "$ENVIRONMENT_RECORD" \
    python3 "$TIER3_ROOT/analysis/tier3c_batch.py" stage-one \
    "$DISCOVERY" "${SLURM_ARRAY_TASK_ID:?array index missing}" "$STAGE_ROOT" "$ENVIRONMENT_RECORD"
