#!/usr/bin/env bash
#SBATCH --job-name=tier3c-run
#SBATCH --partition=workers
#SBATCH --array=0-0%8
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --export=NONE
set -euo pipefail

BATCH=${1:?frozen batch manifest required}
OUTPUT_ROOT=${2:?output root required}
TIER3_ROOT=${3:?repository root required}
ENVIRONMENT_RECORD=${4:?environment record required}
export TIER3_ROOT TIER3_SCRATCH_ROOT="$OUTPUT_ROOT"

exec "$TIER3_ROOT/analysis/slurm/guix_job.sh" "$ENVIRONMENT_RECORD" \
    python3 "$TIER3_ROOT/analysis/tier3c_batch.py" run-one \
    "$BATCH" "${SLURM_ARRAY_TASK_ID:?array index missing}" "$OUTPUT_ROOT"
