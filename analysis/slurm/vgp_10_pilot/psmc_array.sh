#!/usr/bin/env bash
#SBATCH --job-name=vgp10-psmc
set -euo pipefail

if [[ -n ${SLURM_JOB_ID:-} ]]; then
    VGP_STAGE_REPO_ROOT=${SLURM_SUBMIT_DIR:?submit from the repository root}
else
    VGP_STAGE_REPO_ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/../../.." && pwd)
fi
source "$VGP_STAGE_REPO_ROOT/analysis/slurm/vgp_10_pilot/common.sh"
VGP_SELECTION_ID=${1:?usage: psmc_array.sh SELECTION_ID}
export VGP_SELECTION_ID
require_selection_id "$VGP_SELECTION_ID"
require_runtime
replicate=${SLURM_ARRAY_TASK_ID:-0}
[[ $replicate =~ ^([0-9]|[1-9][0-9]|1[0-9][0-9]|200)$ ]] || fail "replicate must be 0..200"
VGP_STAGE_NAME=$(printf 'psmc/replicate-%03d' "$replicate")
export VGP_STAGE_NAME
begin_stage "$VGP_STAGE_NAME"

[[ -f $VGP_PAIR_RUN/consensus/.complete.json ]] || fail "consensus sentinel absent"
verify_tool psmc
psmc=$(tool_path psmc)
threads=${SLURM_CPUS_PER_TASK:-1}
input="$VGP_PAIR_RUN/consensus/consensus/input.psmcfa"
if [[ $replicate -eq 0 ]]; then
    "$psmc" -N25 -t15 -r5 -p '4+25*2+4+6' -o "$VGP_STAGE_PARTIAL/unscaled.psmc" "$input"
else
    unit="$VGP_STAGE_PARTIAL/replicate-$(printf '%03d' "$replicate").psmcfa"
    python3 -m analysis.vgp_10_pilot emit-bootstrap \
        "$VGP_PAIR_RUN/consensus/consensus/consensus.fa" \
        "$VGP_PAIR_RUN/consensus/consensus/bootstrap_units.5mb.bed" \
        "$VGP_PAIR_RUN/consensus/consensus/bootstrap_manifest.tsv" "$replicate" "$unit" \
        >"$VGP_STAGE_PARTIAL/bootstrap_input.json"
    "$psmc" -b -N25 -t15 -r5 -p '4+25*2+4+6' \
        -o "$VGP_STAGE_PARTIAL/bootstrap.unscaled.psmc" "$unit"
    rm -- "$unit"
fi
promote_stage
