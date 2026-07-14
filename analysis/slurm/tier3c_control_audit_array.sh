#!/usr/bin/env bash
#SBATCH --job-name=tier3c-control-audit
#SBATCH --partition=workers
#SBATCH --array=0-1%2
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --export=NONE
set -euo pipefail

RESULT_ROOT=${1:?result root required}
AUDIT_ROOT=${2:?audit output root required}
TIER3_ROOT=${3:?repository root required}
ENVIRONMENT_RECORD=${4:?environment record required}
export TIER3_ROOT TIER3_SCRATCH_ROOT="$AUDIT_ROOT"

case "${SLURM_ARRAY_TASK_ID:?array index missing}" in
    0) DATASET=drosophila.melanogaster.tier3c ;;
    1) DATASET=homo.sapiens.tier3c ;;
    *) echo "control audit index must be 0 or 1" >&2; exit 64 ;;
esac

exec "$TIER3_ROOT/analysis/slurm/guix_job.sh" "$ENVIRONMENT_RECORD" \
    python3 "$TIER3_ROOT/analysis/tier3c_control_audit.py" audit \
    "$RESULT_ROOT/$DATASET.json" "$AUDIT_ROOT/$DATASET.audit.json"
