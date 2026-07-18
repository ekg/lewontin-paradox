#!/usr/bin/env bash
set -euo pipefail

ROOT=$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
CHANNELS="$ROOT/analysis/guix/vgp_10_pilot/channels.scm"
MANIFEST="$ROOT/analysis/guix/vgp_10_pilot/manifest.scm"
LOAD_PATH="$ROOT/analysis/guix"

if [[ -n ${SLURM_JOB_ID:-}${SLURM_ARRAY_JOB_ID:-}${SLURM_ARRAY_TASK_ID:-} ]]; then
    echo "refusing: preflight adjudication runs before the Slurm submission boundary" >&2
    exit 65
fi
if ! command -v guix >/dev/null 2>&1; then
    echo "refusing: GNU Guix is required" >&2
    exit 69
fi

# Resolve the authenticated channel commit and production manifest explicitly.
# The adjudicator then verifies every captured store/executable identity again.
exec guix time-machine -C "$CHANNELS" -- \
    shell -L "$LOAD_PATH" -m "$MANIFEST" --pure -- \
    bash -c 'export PYTHONPATH="$1"; shift; exec python3 "$@"' \
    bash "$ROOT" "$ROOT/analysis/adjudicate_vgp_10_pilot.py" "$@"
