#!/usr/bin/env bash
set -euo pipefail

ROOT=$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
CHANNELS="$ROOT/analysis/guix/channels.scm"
MANIFEST="$ROOT/analysis/guix/manifest.scm"
LOAD_PATH="$ROOT/analysis/guix"

if [[ -n ${SLURM_JOB_ID:-}${SLURM_ARRAY_JOB_ID:-}${SLURM_ARRAY_TASK_ID:-} ]]; then
    echo "refusing: final synthesis is reconciliation, not biological compute" >&2
    exit 65
fi
if ! command -v guix >/dev/null 2>&1; then
    echo "refusing: pinned GNU Guix is required" >&2
    exit 69
fi

if [[ ${1:-} == --all-tests ]]; then
    shift
    command=(python3 -m pytest -q "$ROOT/analysis/tests" "$@")
elif [[ ${1:-} == --tests ]]; then
    shift
    command=(python3 -m pytest -q "$ROOT/analysis/tests/test_vgp_comprehensive_synthesis.py" "$@")
else
    command=(python3 "$ROOT/analysis/synthesize_vgp_comprehensive.py" "$@")
fi

exec guix time-machine -C "$CHANNELS" -- \
    shell -L "$LOAD_PATH" -m "$MANIFEST" --pure -- \
    bash -c 'export PYTHONPATH="$1"; shift; exec "$@"' \
    bash "$ROOT" "${command[@]}"
