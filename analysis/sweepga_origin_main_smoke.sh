#!/usr/bin/env bash
# Hard gate: prove the exact rebuilt origin/main binary accepts literal -n.
set -euo pipefail

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
CHANNELS="$ROOT/analysis/guix/channels.scm"
MANIFEST="$ROOT/analysis/guix/sweepga_impg_smoke_manifest.scm"
SCRATCH=${SWEEPGA_BUILD_SCRATCH:-/moosefs/erikg/tier3scratch/sweepga-origin-main-018e4ce}
BIN=${SWEEPGA_ORIGIN_MAIN_BIN:-$SCRATCH/bin-1/sweepga}
OUT=${SWEEPGA_ORIGIN_MAIN_SMOKE_OUT:-$SCRATCH/smoke-literal-n}
BIO=${TIER3A_BIOLOGICAL_ROOT:-/moosefs/erikg/tier3data/tier3a-acquisition-20260716/spinachia_spinachia_SK-2024b}
EXPECTED_SHA=fa7f0edb9b7e275c288db254046020e136d4267dd5ee043379227ef80da0573b

if [[ ${1:-} != --inside-guix ]]; then
    exec guix time-machine -C "$CHANNELS" -- \
        shell -m "$MANIFEST" --pure -- \
        env ROOT="$ROOT" CHANNELS="$CHANNELS" MANIFEST="$MANIFEST" \
            SCRATCH="$SCRATCH" BIN="$BIN" OUT="$OUT" BIO="$BIO" \
            EXPECTED_SHA="$EXPECTED_SHA" bash "$0" --inside-guix
fi

rm -rf "$OUT"
mkdir -p "$OUT"

sha256sum "$BIN" > "$OUT/binary.sha256"
[[ $(awk '{print $1}' "$OUT/binary.sha256") == "$EXPECTED_SHA" ]] || {
    echo "wrong SweepGA binary" >&2
    exit 1
}
realpath "$BIN" > "$OUT/binary.realpath"
PATH="$(dirname -- "$BIN"):$PATH"
export PATH
type -a sweepga > "$OUT/type-a.txt"
command -v sweepga > "$OUT/command-v.txt"
[[ $(command -v sweepga) == "$BIN" ]] || {
    echo "PATH resolves a different SweepGA" >&2
    exit 1
}

"$BIN" --version > "$OUT/version.stdout" 2> "$OUT/version.stderr"
"$BIN" --help > "$OUT/help.stdout" 2> "$OUT/help.stderr"

# Real, staged Spinachia H1/H2 excerpts.  The literal option is intentionally
# passed to the native binary without any wrapper, alias, or CLI patch.
samtools faidx "$BIO/h1.fna" 'CM106590.1:50001-70000' > "$OUT/h1.fa"
samtools faidx "$BIO/h2.fna" 'CM106672.1:45001-65000' > "$OUT/h2.fa"
set +e
"$BIN" "$OUT/h1.fa" "$OUT/h2.fa" \
    --output-file "$OUT/literal-n-1to1.paf" \
    -n 1:1 --scaffold-jump 0 --threads 2 \
    > "$OUT/literal-n.stdout" 2> "$OUT/literal-n.stderr"
status=$?
set -e
printf '%s\n' "$status" > "$OUT/literal-n.exit-status"

if [[ $status -ne 0 ]]; then
    echo "HARD GATE FAILED: origin/main SweepGA rejected literal -n 1:1" >&2
    exit "$status"
fi
[[ -s $OUT/literal-n-1to1.paf ]] || {
    echo "literal -n command succeeded but emitted no mapping" >&2
    exit 1
}
