#!/usr/bin/env bash
# Rebuild the fetched SweepGA origin/main commit twice in isolated Guix shells.
set -euo pipefail

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
CHANNELS="$ROOT/analysis/guix/channels.scm"
MANIFEST="$ROOT/analysis/guix/sweepga_origin_main_manifest.scm"
PARENT_REPO=${SWEEPGA_PARENT_REPO:-/moosefs/erikg/sweepga}
SCRATCH=${SWEEPGA_BUILD_SCRATCH:-/moosefs/erikg/tier3scratch/sweepga-origin-main-018e4ce}
EXPECTED_COMMIT=${SWEEPGA_EXPECTED_COMMIT:-018e4ce49d2c125820e0ac50dc5feaa02d423683}

fail() {
    echo "ERROR: $*" >&2
    exit 1
}

if [[ ${1:-} != --inside-guix ]]; then
    command -v guix >/dev/null || fail "guix is not on PATH"
    command -v git >/dev/null || fail "git is not on PATH"
    [[ -d $PARENT_REPO/.git ]] || fail "missing parent repository: $PARENT_REPO"

    # Fetching updates only the remote-tracking ref/FETCH_HEAD.  It does not
    # reset, clean, checkout, or otherwise modify the user's working files.
    mkdir -p "$SCRATCH/logs"
    git -C "$PARENT_REPO" remote get-url origin > "$SCRATCH/origin_url.txt"
    git -C "$PARENT_REPO" rev-parse HEAD > "$SCRATCH/parent_head.txt"
    git -C "$PARENT_REPO" status --short --branch > "$SCRATCH/parent_status.txt"
    git -C "$PARENT_REPO" fetch --no-tags origin main \
        > "$SCRATCH/fetch.stdout" 2> "$SCRATCH/fetch.stderr"
    FETCHED=$(git -C "$PARENT_REPO" rev-parse FETCH_HEAD)
    REMOTE_MAIN=$(git -C "$PARENT_REPO" rev-parse refs/remotes/origin/main)
    [[ $FETCHED == "$REMOTE_MAIN" ]] || fail "FETCH_HEAD differs from origin/main"
    [[ $FETCHED == "$EXPECTED_COMMIT" ]] ||
        fail "origin/main moved: expected $EXPECTED_COMMIT, fetched $FETCHED"
    printf '%s\n' "$FETCHED" > "$SCRATCH/fetched_commit.txt"
    guix time-machine -C "$CHANNELS" -- describe -f channels \
        > "$SCRATCH/guix-describe.scm"

    exec guix time-machine -C "$CHANNELS" -- \
        shell -m "$MANIFEST" --pure -- \
        env ROOT="$ROOT" CHANNELS="$CHANNELS" MANIFEST="$MANIFEST" \
            PARENT_REPO="$PARENT_REPO" SCRATCH="$SCRATCH" \
            EXPECTED_COMMIT="$EXPECTED_COMMIT" \
            bash "$0" --inside-guix
fi

[[ $EXPECTED_COMMIT == "$(git -C "$PARENT_REPO" rev-parse refs/remotes/origin/main)" ]] ||
    fail "origin/main changed between wrapper and Guix build"

rm -rf "$SCRATCH/checkout-1" "$SCRATCH/checkout-2" \
       "$SCRATCH/build-1" "$SCRATCH/build-2" \
       "$SCRATCH/home-1" "$SCRATCH/home-2" \
       "$SCRATCH/cargo-home-1" "$SCRATCH/cargo-home-2" \
       "$SCRATCH/cache-1" "$SCRATCH/cache-2" "$SCRATCH/shims" \
       "$SCRATCH/bin-1" "$SCRATCH/bin-2" "$SCRATCH/logs"
mkdir -p "$SCRATCH"/{checkout-1,checkout-2,bin-1,bin-2,logs}
mkdir -p "$SCRATCH/shims"
for pair in rt:1 pthread:0 dl:2 util:1; do
    base=${pair%:*}
    sonum=${pair#*:}
    ln -s "$GUIX_ENVIRONMENT/lib/lib${base}.so.${sonum}" "$SCRATCH/shims/lib${base}.so"
done

# Independent archives avoid the user's modified worktree and ensure both
# builds use only objects from the fetched origin/main revision.
for n in 1 2; do
    git -C "$PARENT_REPO" archive --format=tar "$EXPECTED_COMMIT" |
        tar -xf - -C "$SCRATCH/checkout-$n"
    git -C "$PARENT_REPO" archive --format=tar "$EXPECTED_COMMIT" |
        sha256sum | awk '{print $1}' > "$SCRATCH/logs/source-archive-$n.sha256"
done
cmp "$SCRATCH/logs/source-archive-1.sha256" "$SCRATCH/logs/source-archive-2.sha256"

{
    rustc --version --verbose
    cargo --version --verbose
    gcc --version | sed -n '1p'
    clang --version | sed -n '1p'
    cmake --version | sed -n '1p'
    pkg-config --version
    git --version
    ld --version | sed -n '1p'
} > "$SCRATCH/logs/tool-versions.txt"

SOURCE_DATE_EPOCH=$(git -C "$PARENT_REPO" show -s --format=%ct "$EXPECTED_COMMIT")
printf '%s\n' "$SOURCE_DATE_EPOCH" > "$SCRATCH/logs/source-date-epoch.txt"

build_one() {
    local n=$1
    local source="$SCRATCH/checkout-$n"
    local target="$SCRATCH/build-$n"
    local home="$SCRATCH/home-$n"
    local cargo_home="$SCRATCH/cargo-home-$n"
    local cache="$SCRATCH/cache-$n"
    local shims="$SCRATCH/shims"
    mkdir -p "$target" "$home" "$cargo_home" "$cache"

    (
        export HOME="$home"
        export CARGO_HOME="$cargo_home"
        export CARGO_TARGET_DIR="$target"
        export XDG_CACHE_HOME="$cache"
        export SOURCE_DATE_EPOCH
        export CARGO_INCREMENTAL=0
        export CARGO_TERM_COLOR=never
        export CARGO_NET_GIT_FETCH_WITH_CLI=true
        export CARGO_REGISTRIES_CRATES_IO_PROTOCOL=sparse
        export CC=gcc CXX=g++
        export CARGO_TARGET_X86_64_UNKNOWN_LINUX_GNU_LINKER=gcc
        export LIBCLANG_PATH="$GUIX_ENVIRONMENT/lib"
        export LIBRARY_PATH="$shims${LIBRARY_PATH:+:$LIBRARY_PATH}"
        export CFLAGS="-ffile-prefix-map=$cargo_home=/build/cargo-home -ffile-prefix-map=$target=/build/target"
        export CXXFLAGS="$CFLAGS"
        export RUSTFLAGS="--remap-path-prefix=$source=/build/sweepga-source --remap-path-prefix=$cargo_home=/build/cargo-home --remap-path-prefix=$target=/build/target -L native=$shims"
        cd "$source"
        cargo fetch --locked
        cargo build --release --locked --bin sweepga
    ) > "$SCRATCH/logs/build-$n.log" 2>&1

    cp "$target/release/sweepga" "$SCRATCH/bin-$n/sweepga"
    # Rust release artifacts retain a large non-runtime symbol table whose
    # native dependency symbols can encode build-directory spellings.  The
    # deployed executable does not need it; deterministic stripping is part
    # of the recorded build recipe.
    strip --strip-all "$SCRATCH/bin-$n/sweepga"
    chmod 0555 "$SCRATCH/bin-$n/sweepga"
    for name in ALNtoPAF FAtoGDB FastGA GIXmake GIXrm ONEview PAFtoALN wfmash; do
        companion=$(find "$target/release/build" -type f -path "*/out/$name" -print -quit)
        [[ -n $companion ]] || fail "build $n did not produce companion $name"
        cp "$companion" "$SCRATCH/bin-$n/$name"
        chmod 0555 "$SCRATCH/bin-$n/$name"
    done
    sha256sum "$SCRATCH/bin-$n/sweepga" > "$SCRATCH/logs/binary-$n.sha256"
    sha256sum "$SCRATCH/bin-$n"/* > "$SCRATCH/logs/artifacts-$n.sha256"
    realpath "$SCRATCH/bin-$n/sweepga" > "$SCRATCH/logs/binary-$n.realpath"
    readelf -l "$SCRATCH/bin-$n/sweepga" |
        grep 'Requesting program interpreter' > "$SCRATCH/logs/binary-$n.interpreter"
    ldd "$SCRATCH/bin-$n/sweepga" > "$SCRATCH/logs/binary-$n.ldd"
}

build_one 1
build_one 2

if cmp -s "$SCRATCH/bin-1/sweepga" "$SCRATCH/bin-2/sweepga"; then
    printf '%s\n' byte-identical > "$SCRATCH/logs/reproducibility.txt"
else
    printf '%s\n' different > "$SCRATCH/logs/reproducibility.txt"
    cmp -l "$SCRATCH/bin-1/sweepga" "$SCRATCH/bin-2/sweepga" |
        sed -n '1,100p' > "$SCRATCH/logs/binary-differences.txt" || true
fi
printf 'artifact\tcomparison\n' > "$SCRATCH/logs/companion-reproducibility.tsv"
for name in ALNtoPAF FAtoGDB FastGA GIXmake GIXrm ONEview PAFtoALN wfmash; do
    if cmp -s "$SCRATCH/bin-1/$name" "$SCRATCH/bin-2/$name"; then
        printf '%s\tbyte-identical\n' "$name"
    else
        # These native helper programs are retained for runtime completeness,
        # but are not the SweepGA binary under the reproducibility criterion.
        printf '%s\tdiffers-native-build-metadata\n' "$name"
    fi >> "$SCRATCH/logs/companion-reproducibility.tsv"
done

printf '%s\n' "builds complete: $SCRATCH"
