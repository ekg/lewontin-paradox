#!/usr/bin/env bash
set -euo pipefail

ROOT=$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
CHANNELS="$ROOT/analysis/guix/channels.scm"
MANIFEST="$ROOT/analysis/guix/manifest.scm"
LOAD_PATH="$ROOT/analysis/guix"

if [[ -n ${SLURM_JOB_ID:-}${SLURM_ARRAY_JOB_ID:-}${SLURM_ARRAY_TASK_ID:-} ]]; then
    echo "refusing: VGP acquisition must not run in Slurm" >&2
    exit 65
fi
if ! command -v guix >/dev/null 2>&1; then
    echo "refusing: GNU Guix is required" >&2
    exit 69
fi

exec guix time-machine -C "$CHANNELS" -- \
    shell -L "$LOAD_PATH" -m "$MANIFEST" --pure -- \
    bash -c 'export VGPPILOT_PINNED_GUIX=1 PYTHONPATH=$1 SSL_CERT_FILE=/etc/ssl/certs/ca-certificates.crt; shift; exec python3 "$@"' \
    bash "$ROOT" "$ROOT/analysis/vgp_10_pilot_acquisition.py" "$@"
