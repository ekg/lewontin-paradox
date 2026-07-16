#!/usr/bin/env bash
#SBATCH --job-name=tier3b-recovery
#SBATCH --partition=workers
#SBATCH --array=0-1%1
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=04:00:00
#SBATCH --export=NONE
#SBATCH --output=results/tier3b/run_logs/population-%A_%a.out
#SBATCH --error=results/tier3b/run_logs/population-%A_%a.err
set -euo pipefail

TIER3_ROOT=${1:?absolute repository root required}
MANIFEST=${2:?absolute acquisition manifest required}
OUTPUT_DIR=${3:?absolute raw output directory required}
ENVIRONMENT_RECORD=${4:?absolute Guix environment record required}
TIER3_SCRATCH_ROOT=${5:?absolute scratch root required}
export TIER3_ROOT TIER3_ENVIRONMENT_RECORD TIER3_SCRATCH_ROOT

exec "$TIER3_ROOT/analysis/slurm/guix_job.sh" "$ENVIRONMENT_RECORD" \
    python3 "$TIER3_ROOT/analysis/tier3b_population_recovery.py" run-one \
    --manifest "$MANIFEST" \
    --index "${SLURM_ARRAY_TASK_ID:?array index missing}" \
    --output-dir "$OUTPUT_DIR"
