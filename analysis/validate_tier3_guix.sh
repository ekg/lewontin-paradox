#!/usr/bin/env bash
set -euo pipefail

ROOT=$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
CHANNELS="$ROOT/analysis/guix/channels.scm"
MANIFEST="$ROOT/analysis/guix/manifest.scm"
LOAD_PATH="$ROOT/analysis/guix"
TMP_ROOT=$(mktemp -d "${TMPDIR:-/tmp}/tier3-guix.XXXXXXXX")
trap 'rm -rf "$TMP_ROOT"' EXIT

if ! command -v guix >/dev/null 2>&1; then
    echo "GNU Guix is required" >&2
    exit 1
fi

build_once() {
    local profile=$1
    guix time-machine -C "$CHANNELS" -- \
        package -L "$LOAD_PATH" -p "$profile" -m "$MANIFEST" --no-grafts \
        >/dev/null
    readlink -f "$profile"
}

build_once "$TMP_ROOT/profile-a" >"$TMP_ROOT/paths-a"
build_once "$TMP_ROOT/profile-b" >"$TMP_ROOT/paths-b"
diff -u "$TMP_ROOT/paths-a" "$TMP_ROOT/paths-b"

PROFILE_PATH=$(cat "$TMP_ROOT/paths-a")
case "$PROFILE_PATH" in
    /gnu/store/*) ;;
    *) echo "manifest did not realize to a Guix store path: $PROFILE_PATH" >&2; exit 1 ;;
esac

guix time-machine -C "$CHANNELS" -- \
    shell -L "$LOAD_PATH" -m "$MANIFEST" --pure -- \
    bash -c '
        set -euo pipefail
        for tool in python3 pytest samtools bcftools bgzip tabix bedtools vcftools wfmash; do
            command -v "$tool" >/dev/null
        done
        python3 - <<"PY"
from importlib.metadata import version
import Bio, numpy, pandas, pyfaidx, pysam, scipy
print("python scientific versions:", Bio.__version__, version("jsonschema"), numpy.__version__, pandas.__version__, pyfaidx.__version__, pysam.__version__, scipy.__version__)
PY
        python3 -m pytest -q analysis/tests/test_common.py analysis/tests/test_manifest.py
        wfmash --version
        bcftools --version
        samtools --version
        bgzip --version
        bedtools --version
        vcftools --version
        if command -v impg >/dev/null; then
            echo "unapproved IMPG executable leaked into the pure profile" >&2
            exit 1
        fi
    '

echo "reproducible Tier 3 profile: $PROFILE_PATH"
