#!/usr/bin/env bash
#SBATCH --job-name=tier3b-independent
#SBATCH --partition=workers
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=03:00:00
#SBATCH --export=NONE
#SBATCH --output=results/tier3b/run_logs/population-independent-%j.out
#SBATCH --error=results/tier3b/run_logs/population-independent-%j.err
set -euo pipefail

TIER3_ROOT=${1:?absolute repository root required}
MANIFEST=${2:?absolute acquisition manifest required}
OUTPUT_DIR=${3:?absolute subset output directory required}
ENVIRONMENT_RECORD=${4:?absolute Guix environment record required}
TIER3_SCRATCH_ROOT=${5:?absolute scratch root required}
export TIER3_ROOT TIER3_ENVIRONMENT_RECORD TIER3_SCRATCH_ROOT

exec "$TIER3_ROOT/analysis/slurm/guix_job.sh" "$ENVIRONMENT_RECORD" \
    python3 "$TIER3_ROOT/analysis/tier3b_population_recovery.py" independent-check \
    --manifest "$MANIFEST" \
    --output-dir "$OUTPUT_DIR"
