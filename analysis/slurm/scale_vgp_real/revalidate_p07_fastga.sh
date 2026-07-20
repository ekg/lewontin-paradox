#!/usr/bin/env bash
# Recompute the frozen P07 SweepGA mapping while continuously proving that
# FastGA and every companion process remain inside private node-local /scratch.
set -euo pipefail

: "${VGP_ROOT:?VGP_ROOT is required}"
: "${VGP_REPO_ROOT:?VGP_REPO_ROOT is required}"
: "${VGP_SCALE_ROOT:?VGP_SCALE_ROOT is required}"
: "${VGP_ENVIRONMENT_CAPTURE:?VGP_ENVIRONMENT_CAPTURE is required}"
: "${VGP_FASTGA_AMENDMENT:?VGP_FASTGA_AMENDMENT is required}"
[[ $VGP_ROOT == /moosefs/erikg/vgp ]] || { echo "ERROR: noncanonical VGP_ROOT" >&2; exit 2; }
[[ $VGP_SCALE_ROOT == "$VGP_ROOT"/derived/scale-vgp-real-v1 ]] || {
    echo "ERROR: scale root does not derive from VGP_ROOT" >&2; exit 2;
}
[[ -d /scratch && -w /scratch ]] || { echo "ERROR: /scratch unavailable" >&2; exit 2; }
case $(stat -f -c %T -- /scratch) in
    nfs|nfs4|fuse*|lustre|gpfs|ceph) echo "ERROR: /scratch is not node-local" >&2; exit 2 ;;
esac

input_manifest=$VGP_SCALE_ROOT/inputs/P07.json
prior_mapping=$VGP_ROOT/pilot/outputs/vgp10-auth-20260718-v2/P07/core/mapping
[[ -f $input_manifest && -f $prior_mapping/h2_to_h1.1to1.paf ]] || {
    echo "ERROR: canonical inputs or frozen P07 mapping absent" >&2; exit 2;
}

readarray -t resolved < <("${VGP_REPO_ROOT}/analysis/slurm/scale_vgp_real/resolve_p07_inputs.py" \
    "$input_manifest" "$VGP_ENVIRONMENT_CAPTURE" "$VGP_FASTGA_AMENDMENT")
h1_gz=${resolved[0]}
h2_gz=${resolved[1]}
h1_fa_sha256=${resolved[2]}
h2_fa_sha256=${resolved[3]}
sweepga=${resolved[4]}
capture_sha256=${resolved[5]}
amendment_sha256=${resolved[6]}

scratch=$(mktemp -d "/scratch/vgp-scale-fastga-P07-${SLURM_JOB_ID:?}-XXXXXX")
scratch_resolved=$(readlink -f -- "$scratch")
node_local_base_resolved=$(readlink -f -- /scratch)
cleanup() {
    status=$?
    if [[ $scratch != /scratch/vgp-scale-fastga-P07-"$SLURM_JOB_ID"-* \
          || $scratch_resolved != "$node_local_base_resolved"/vgp-scale-fastga-P07-"$SLURM_JOB_ID"-* ]]; then
        echo "ERROR: refusing cleanup outside validated requested/resolved scratch roots" >&2
        exit 70
    fi
    rm -rf -- "$scratch_resolved"
    exit "$status"
}
trap cleanup EXIT INT TERM
mkdir -p "$scratch/inputs" "$scratch/result.partial"
export TMPDIR=$scratch TMP=$scratch TEMP=$scratch PYTHONPATH=$VGP_REPO_ROOT
cd -- "$scratch"
cwd_resolved=$(readlink -f -- /proc/$$/cwd)
[[ $cwd_resolved == "$scratch_resolved" ]] || {
    echo "ERROR: batch cwd outside private node-local /scratch: $cwd_resolved" >&2; exit 70;
}
for variable in TMPDIR TMP TEMP; do
    [[ $(readlink -f -- "${!variable}") == "$scratch_resolved" ]] || {
        echo "ERROR: $variable outside private node-local /scratch" >&2; exit 70;
    }
done

# The immutable compressed CAS objects remain on MooseFS.  FASTAs are created
# only inside this job's private scratch directory.
gzip -dc -- "$h1_gz" > "$scratch/inputs/h1.fa"
gzip -dc -- "$h2_gz" > "$scratch/inputs/h2.fa"
[[ $(sha256sum "$scratch/inputs/h1.fa" | awk '{print $1}') == "$h1_fa_sha256" ]] || {
    echo "ERROR: staged H1 digest mismatch" >&2; exit 2;
}
[[ $(sha256sum "$scratch/inputs/h2.fa" | awk '{print $1}') == "$h2_fa_sha256" ]] || {
    echo "ERROR: staged H2 digest mismatch" >&2; exit 2;
}

partial=$scratch/result.partial
started=$(date -u +%Y-%m-%dT%H:%M:%SZ)
"$sweepga" "$scratch/inputs/h2.fa" "$scratch/inputs/h1.fa" \
    --output-file "$partial/h2_to_h1.native.1to1.paf" \
    --num-mappings 1:1 --scaffold-jump 0 --overlap 0 \
    --scoring log-length-ani --threads "${SLURM_CPUS_PER_TASK:-1}" \
    >"$partial/sweepga.stdout" 2>"$partial/sweepga.stderr" &
sweepga_pid=$!
guard_log=$partial/fastga_scratch_snapshots.jsonl
guard_failed=0
while [[ $(ps -o stat= -p "$sweepga_pid" 2>/dev/null) != *Z* ]] && kill -0 "$sweepga_pid" 2>/dev/null; do
    if ! python3 "$VGP_REPO_ROOT/analysis/fastga_scratch_guard.py" check \
        --parent-pid "$sweepga_pid" --scratch "$scratch" --audit-jsonl "$guard_log"; then
        guard_failed=1
        pkill -TERM -P "$sweepga_pid" 2>/dev/null || true
        kill -TERM "$sweepga_pid" 2>/dev/null || true
        break
    fi
    sleep "${VGP_FASTGA_GUARD_INTERVAL_SECONDS:-2}"
done
if wait "$sweepga_pid"; then sweepga_status=0; else sweepga_status=$?; fi
(( guard_failed == 0 )) || {
    echo "ERROR: hard infrastructure error: FastGA escaped private node-local /scratch" >&2
    exit 70
}
(( sweepga_status == 0 )) || exit "$sweepga_status"
python3 "$VGP_REPO_ROOT/analysis/fastga_scratch_guard.py" finalize \
    --scratch "$scratch" --audit-jsonl "$guard_log" \
    --output "$partial/fastga_scratch_contract.json"
python3 -m analysis.vgp_10_pilot enforce-paf \
    "$partial/h2_to_h1.native.1to1.paf" "$partial/h2_to_h1.1to1.paf" \
    > "$partial/exact_multiplicity_filter.json"
python3 -m analysis.vgp_10_pilot audit-paf \
    "$partial/h2_to_h1.1to1.paf" "$scratch/inputs/h1.fa" "$scratch/inputs/h2.fa" \
    > "$partial/multiplicity.json"

python3 "$VGP_REPO_ROOT/analysis/slurm/scale_vgp_real/finalize_p07_fastga.py" \
    --partial "$partial" --input-manifest "$input_manifest" \
    --prior-mapping "$prior_mapping" --scratch "$scratch" \
    --job-id "$SLURM_JOB_ID" --started "$started" \
    --capture-sha256 "$capture_sha256" --amendment-sha256 "$amendment_sha256"

# Promotion starts only after alignment success, multiplicity audit, live guard
# finalization, and exact equality to the frozen mapping used downstream.
final=$VGP_SCALE_ROOT/fastga/P07
[[ ! -e $final ]] || { echo "ERROR: final P07 FastGA validation already exists" >&2; exit 2; }
mkdir -p "$(dirname "$final")"
promote=$VGP_SCALE_ROOT/fastga/.P07."$SLURM_JOB_ID".partial
[[ ! -e $promote ]] || { echo "ERROR: stale promotion path" >&2; exit 2; }
cp -a -- "$partial" "$promote"
sync "$promote"
mv -- "$promote" "$final"
cmp -- "$partial/contract.json" "$final/contract.json"
