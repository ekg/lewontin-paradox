#!/usr/bin/env bash
#SBATCH --job-name=tier3-guix-smoke
#SBATCH --partition=workers
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --export=NONE
#SBATCH --output=analysis/pilot_results/logs/compute-smoke-%j.out
#SBATCH --error=analysis/pilot_results/logs/compute-smoke-%j.err
set -euo pipefail

ENVIRONMENT_RECORD=${1:?environment record required}
TIER3_ROOT=${2:?repository root required}
OUTPUT=${3:-$TIER3_ROOT/analysis/pilot_results/compute_smoke.json}
TIER3_SCRATCH_ROOT=${4:-$(dirname "$ENVIRONMENT_RECORD")/compute-smoke-scratch}
mkdir -p "$TIER3_SCRATCH_ROOT"
export TIER3_ROOT TIER3_ENVIRONMENT_RECORD="$ENVIRONMENT_RECORD" TIER3_SCRATCH_ROOT

exec "$TIER3_ROOT/analysis/slurm/guix_job.sh" "$ENVIRONMENT_RECORD" \
    python3 "$TIER3_ROOT/analysis/run_tier3.py" compute-smoke \
    --environment-record "$ENVIRONMENT_RECORD" --output "$OUTPUT"
