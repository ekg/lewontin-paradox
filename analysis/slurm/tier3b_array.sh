#!/usr/bin/env bash
#SBATCH --job-name=tier3b-popvcf
#SBATCH --partition=workers
#SBATCH --array=0-0
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --export=NONE
#SBATCH --output=analysis/pilot_results/logs/tier3b-%A_%a.out
#SBATCH --error=analysis/pilot_results/logs/tier3b-%A_%a.err
set -euo pipefail

WORKFLOW=${1:?workflow path required}
TIER3_ROOT=${2:?repository root required}
ENVIRONMENT_RECORD=${3:-${TIER3_ENVIRONMENT_RECORD:?environment record required}}
TIER3_SCRATCH_ROOT=$(/usr/bin/python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["scratch_root"])' "$WORKFLOW")
export TIER3_ROOT TIER3_WORKFLOW="$WORKFLOW" TIER3_ENVIRONMENT_RECORD="$ENVIRONMENT_RECORD" TIER3_SCRATCH_ROOT

exec "$TIER3_ROOT/analysis/slurm/guix_job.sh" "$ENVIRONMENT_RECORD" \
    python3 "$TIER3_ROOT/analysis/run_tier3.py" run-array "$WORKFLOW" \
    --tier 3b --index "${SLURM_ARRAY_TASK_ID:?array index missing}"
