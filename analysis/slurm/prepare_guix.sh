#!/usr/bin/env bash
set -euo pipefail

ROOT=$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)
STATE_DIR=${1:?usage: prepare_guix.sh ABSOLUTE_STATE_DIRECTORY}
case "$STATE_DIR" in
    /*) ;;
    *) echo "state directory must be absolute" >&2; exit 64 ;;
esac

case "$(hostname -s)" in
    octopus01) ;;
    *) echo "the frozen topology permits Guix realization only on octopus01" >&2; exit 69 ;;
esac

umask 077
mkdir -p "$STATE_DIR/gcroots"
CHANNELS=$ROOT/analysis/guix/channels.scm
MANIFEST=$ROOT/analysis/guix/manifest.scm
LOAD_PATH=$ROOT/analysis/guix
PROFILE_LINK=$STATE_DIR/profile
GC_ROOT=$STATE_DIR/gcroots/tier3-profile
DERIVATIONS=$STATE_DIR/derivations.txt
STORE_PATHS=$STATE_DIR/store-paths.txt
RESOLVED_CHANNELS=$STATE_DIR/resolved-channels.scm
RECORD=$STATE_DIR/environment.json

guix time-machine -C "$CHANNELS" -- \
    package -L "$LOAD_PATH" -p "$PROFILE_LINK" -m "$MANIFEST" --no-grafts
PROFILE_STORE=$(readlink -f "$PROFILE_LINK")
case "$PROFILE_STORE" in
    /gnu/store/????????????????????????????????-*) ;;
    *) echo "profile did not resolve into the Guix store" >&2; exit 70 ;;
esac

guix gc --derivers "$PROFILE_STORE" >"$STATE_DIR/profile-deriver.txt"
[[ $(wc -l <"$STATE_DIR/profile-deriver.txt") == 1 ]] || {
    echo "realized profile does not have exactly one recorded derivation" >&2
    exit 70
}
PROFILE_DERIVER=$(<"$STATE_DIR/profile-deriver.txt")
rm -f "$GC_ROOT"
guix build --root="$GC_ROOT" "$PROFILE_DERIVER" >/dev/null
[[ $(readlink -f "$GC_ROOT") == "$PROFILE_STORE" ]] || {
    echo "persistent GC root does not resolve to the realized profile" >&2
    exit 70
}

guix time-machine -C "$CHANNELS" -- \
    build -L "$LOAD_PATH" -m "$MANIFEST" --no-grafts --derivations \
    >"$DERIVATIONS.partial"
sort -u "$DERIVATIONS.partial" "$STATE_DIR/profile-deriver.txt" >"$DERIVATIONS"
rm -f "$DERIVATIONS.partial"
guix gc --requisites "$PROFILE_STORE" | sort -u >"$STORE_PATHS.partial"
mv "$STORE_PATHS.partial" "$STORE_PATHS"
guix time-machine -C "$CHANNELS" -- describe --format=channels >"$RESOLVED_CHANNELS.partial"
mv "$RESOLVED_CHANNELS.partial" "$RESOLVED_CHANNELS"

PYTHON_SITE=$(find "$PROFILE_STORE/lib" -maxdepth 2 -type d -path '*/python*/site-packages' -print -quit)
env \
    PATH="$PROFILE_STORE/bin:/usr/bin:/bin" \
    GUIX_PROFILE="$PROFILE_STORE" \
    GUIX_PYTHONPATH="$PYTHON_SITE" \
    PYTHONPATH="$ROOT" \
    "$PROFILE_STORE/bin/python3" "$ROOT/analysis/run_tier3.py" record-environment \
        --profile "$GC_ROOT" \
        --derivations "$DERIVATIONS" \
        --store-paths "$STORE_PATHS" \
        --resolved-channels "$RESOLVED_CHANNELS" \
        --output "$RECORD" \
        >/dev/null

echo "$RECORD"
