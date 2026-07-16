#!/usr/bin/env bash
#SBATCH --job-name=tier3b-acq-power
#SBATCH --partition=workers
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --export=NONE
#SBATCH --output=/moosefs/erikg/tier3scratch/tier3b-acquisition/repair-logs/acquisition-power-%j.out
#SBATCH --error=/moosefs/erikg/tier3scratch/tier3b-acquisition/repair-logs/acquisition-power-%j.err
set -euo pipefail

REPO=${1:?absolute repository root required}
ENVIRONMENT_RECORD=${2:?absolute Guix environment record required}
ROOT=${3:?absolute acquisition root required}
MODE=${4:?current or repair required}
CANDIDATE=${5:?candidate path required}
OUTPUT=${6:?output JSON path required}
export TIER3_SCRATCH_ROOT=$ROOT/environment

exec "$REPO/analysis/slurm/guix_job.sh" "$ENVIRONMENT_RECORD" \
  python3 "$REPO/results/tier3b/acquisition_repair_power.py" \
  --root "$ROOT" --mode "$MODE" --candidate "$CANDIDATE" --output "$OUTPUT"
