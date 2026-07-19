#!/usr/bin/env bash
set -euo pipefail

ROOT=$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
MANIFEST="$ROOT/analysis/guix/manifest.scm"
MIRROR="$ROOT/analysis/mirror_vgp_freeze1.py"
WORKER_UNIT=vgp-freeze1-mirror.service
ADJUDICATOR_UNIT=vgp-freeze1-conflict-adjudicator.service
ADJUDICATOR_TIMER=vgp-freeze1-conflict-adjudicator.timer
GUIX=$(command -v guix || true)

if [[ -z $GUIX ]]; then
    echo "refusing: GNU Guix is required" >&2
    exit 69
fi

command=${1:-status}
case "$command" in
    start)
        concurrency=${VGP_MIRROR_CONCURRENCY:-4}
        if [[ ! $concurrency =~ ^[1-8]$ ]]; then
            echo "refusing: VGP_MIRROR_CONCURRENCY must be between 1 and 8" >&2
            exit 64
        fi
        if systemctl --user is-active --quiet "$WORKER_UNIT"; then
            echo "$WORKER_UNIT is already active" >&2
            exit 0
        fi
        systemd-run --user --unit="${WORKER_UNIT%.service}" --collect \
            --description="Restart-safe official VGP Freeze 1 mirror" \
            --working-directory="$ROOT" \
            --property=Restart=on-failure \
            --property=RestartSec=60s \
            --property=TimeoutStopSec=180s \
            "$GUIX" shell -m "$MANIFEST" -- \
            python3 "$MIRROR" worker --concurrency "$concurrency"
        if ! systemctl --user is-active --quiet "$ADJUDICATOR_TIMER"; then
            systemd-run --user --unit="${ADJUDICATOR_UNIT%.service}" \
                --description="Adjudicate isolated VGP Freeze 1 metadata conflicts" \
                --working-directory="$ROOT" \
                --on-active=2m --on-unit-active=5m \
                "$GUIX" shell -m "$MANIFEST" -- \
                python3 "$MIRROR" adjudicate-conflicts
        fi
        ;;
    status)
        systemctl --user --no-pager --full status "$WORKER_UNIT" "$ADJUDICATOR_TIMER" || true
        "$GUIX" shell -m "$MANIFEST" -- python3 "$MIRROR" status
        ;;
    stop)
        systemctl --user stop "$ADJUDICATOR_TIMER" "$ADJUDICATOR_UNIT" "$WORKER_UNIT"
        ;;
    *)
        echo "usage: $0 {start|status|stop}" >&2
        exit 64
        ;;
esac
