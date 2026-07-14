#!/usr/bin/env bash
#SBATCH --job-name=tier3c-outlier-retry
#SBATCH --partition=workers
#SBATCH --array=0-0%1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --export=NONE
set -euo pipefail

BATCH=${1:?frozen batch manifest required}
OUTPUT_ROOT=${2:?output root required}
TIER3_ROOT=${3:?repository root required}
ENVIRONMENT_RECORD=${4:?environment record required}
INDEX_FILE=${5:?newline-delimited frozen batch indices required}
export TIER3_ROOT TIER3_SCRATCH_ROOT="$OUTPUT_ROOT"

INDEX=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$INDEX_FILE")
case "$INDEX" in
    ''|*[!0-9]*) echo "outlier index file yielded an invalid batch index: $INDEX" >&2; exit 64 ;;
esac

exec "$TIER3_ROOT/analysis/slurm/guix_job.sh" "$ENVIRONMENT_RECORD" \
    python3 "$TIER3_ROOT/analysis/tier3c_batch.py" run-one \
    "$BATCH" "$INDEX" "$OUTPUT_ROOT"
