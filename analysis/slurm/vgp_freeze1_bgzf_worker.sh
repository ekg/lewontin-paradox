#!/usr/bin/env bash
set -euo pipefail

: "${SLURM_ARRAY_TASK_ID:?worker requires a Slurm array task}"
: "${VGP_BGZF_PROFILE:?pinned Guix profile is required}"
: "${VGP_BGZF_INVENTORY:?inventory is required}"
: "${VGP_BGZF_DERIVED_ROOT:?derived root is required}"
: "${VGP_BGZF_PIPELINE:?pipeline path is required}"

export PATH="$VGP_BGZF_PROFILE/bin"
export LC_ALL=C
export LANG=C
umask 0027

scratch_created=0
if [[ -z ${SLURM_TMPDIR:-} ]]; then
    : "${VGP_NODE_LOCAL_BASE:=/scratch}"
    [[ -d $VGP_NODE_LOCAL_BASE && -w $VGP_NODE_LOCAL_BASE ]] || {
        echo "node-local scratch base is absent or unwritable: $VGP_NODE_LOCAL_BASE" >&2
        exit 2
    }
    scratch_fs=$(stat -f -c %T -- "$VGP_NODE_LOCAL_BASE")
    case "$scratch_fs" in
        nfs|nfs4|fuse*|moosefs|lustre|gpfs)
            echo "refusing shared scratch filesystem: $scratch_fs" >&2
            exit 2 ;;
    esac
    SLURM_TMPDIR=$(mktemp -d -- "${VGP_NODE_LOCAL_BASE%/}/vgp-bgzf-${SLURM_JOB_ID}-${SLURM_ARRAY_TASK_ID}-XXXXXX")
    scratch_created=1
    export SLURM_TMPDIR
fi
[[ -d $SLURM_TMPDIR ]] || { echo "invalid SLURM_TMPDIR: $SLURM_TMPDIR" >&2; exit 2; }

cleanup_scratch() {
    if [[ $scratch_created == 1 && -n $SLURM_TMPDIR && $SLURM_TMPDIR = "${VGP_NODE_LOCAL_BASE%/}"/vgp-bgzf-* ]]; then
        rm -rf -- "$SLURM_TMPDIR"
    fi
}
trap cleanup_scratch EXIT

# Moderate array concurrency is the primary shared-I/O control.  Low-priority
# best-effort I/O and CPU scheduling keep these streaming conversions behind
# interactive and pilot workloads on each node.
ionice -c 2 -n 6 nice -n 5 \
    python3 "$VGP_BGZF_PIPELINE" worker \
        --inventory "$VGP_BGZF_INVENTORY" \
        --derived-root "$VGP_BGZF_DERIVED_ROOT" \
        --task-index "$SLURM_ARRAY_TASK_ID" \
        --scratch-root "$SLURM_TMPDIR"
