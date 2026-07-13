#!/usr/bin/env bash
set -euo pipefail

ROOT=$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)
WORKFLOW=${1:?usage: submit_tier3.sh ABSOLUTE_WORKFLOW ABSOLUTE_STATE_DIRECTORY}
STATE_DIR=${2:?usage: submit_tier3.sh ABSOLUTE_WORKFLOW ABSOLUTE_STATE_DIRECTORY}
case "$WORKFLOW:$STATE_DIR" in
    /*:/*) ;;
    *) echo "workflow and state directory must be absolute" >&2; exit 64 ;;
esac

umask 077
mkdir -p "$ROOT/analysis/pilot_results/logs"
ENVIRONMENT_RECORD=$("$ROOT/analysis/slurm/prepare_guix.sh" "$STATE_DIR")
export TIER3_ROOT="$ROOT" TIER3_ENVIRONMENT_RECORD="$ENVIRONMENT_RECORD"

"$ROOT/analysis/slurm/guix_job.sh" "$ENVIRONMENT_RECORD" \
    python3 "$ROOT/analysis/run_tier3.py" preflight "$WORKFLOW" >/dev/null

SCRATCH_ROOT=$(
    "$ROOT/analysis/slurm/guix_job.sh" "$ENVIRONMENT_RECORD" \
        python3 - "$WORKFLOW" <<'PY'
import json
import sys
from pathlib import Path
print(Path(json.load(open(sys.argv[1], encoding="utf-8"))["scratch_root"]).resolve())
PY
)
mkdir -p "$SCRATCH_ROOT"
FREE_BYTES=$(df -PB1 "$SCRATCH_ROOT" | awk 'NR==2 {print $4}')
MINIMUM_FREE=536870912000
(( FREE_BYTES >= MINIMUM_FREE )) || {
    echo "Tier 3 submission requires at least $MINIMUM_FREE free bytes; observed $FREE_BYTES" >&2
    exit 75
}
export TIER3_SCRATCH_ROOT="$SCRATCH_ROOT"

SMOKE_JOB=$(sbatch --parsable \
    "$ROOT/analysis/slurm/compute_smoke.sh" "$ENVIRONMENT_RECORD" "$ROOT")
printf 'compute_smoke_job=%s\n' "$SMOKE_JOB"

for TIER in 3a 3b 3c; do
    COUNT=$(
        "$ROOT/analysis/slurm/guix_job.sh" "$ENVIRONMENT_RECORD" \
            python3 "$ROOT/analysis/run_tier3.py" count "$WORKFLOW" --tier "$TIER"
    )
    if (( COUNT == 0 )); then
        continue
    fi
    SCRIPT=$ROOT/analysis/slurm/tier${TIER}_array.sh
    JOB=$(sbatch --parsable --dependency="afterok:$SMOKE_JOB" --array="0-$((COUNT - 1))" \
        "$SCRIPT" "$WORKFLOW" "$ROOT" "$ENVIRONMENT_RECORD")
    printf 'tier%s_array_job=%s count=%s\n' "$TIER" "$JOB" "$COUNT"
done
