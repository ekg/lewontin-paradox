#!/usr/bin/env bash
set -euo pipefail

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/../../.." && pwd)
mode=${1:---dry-run}
[[ $mode == --dry-run || $mode == --submit ]] || {
    echo "usage: submit.sh [--dry-run|--submit]" >&2
    exit 2
}
: "${VGP_DATA_ROOT:?VGP_DATA_ROOT is required}"
: "${VGP_RUN_ID:?VGP_RUN_ID is required}"
: "${VGP_ENVIRONMENT_CAPTURE:?VGP_ENVIRONMENT_CAPTURE is required}"

submit() {
    if [[ $mode == --dry-run ]]; then
        printf '%q ' "$@"
        printf '\n'
    else
        "$@"
    fi
}

read_limits() {
    local resource=$1
    local stage=$2
    python3 - "$resource" "$stage" <<'PY'
import json,sys
v=json.load(open(sys.argv[1]))
stage=v.get("stages",{}).get(sys.argv[2])
if not isinstance(stage,dict): raise SystemExit(f"missing stage-specific estimate: {sys.argv[2]}")
for key in ("cpus_per_task","slurm_time","slurm_mem","scratch_bytes_high"):
    if not stage.get(key): raise SystemExit(f"missing scheduler-specific limit: {sys.argv[2]}.{key}")
    print(stage[key])
PY
}

for selection_id in P01 P02 P03 P04 P05 P06 P07 P08 P09 P10; do
    resource="$VGP_DATA_ROOT/pilot/inputs/$selection_id/resources.json"
    [[ -f $resource ]] || {
        echo "ERROR: measured resource estimate absent: $resource" >&2
        exit 2
    }
    dependency=
    for stage in preflight mapping impg variants consensus; do
        readarray -t limits < <(read_limits "$resource" "$stage")
        command=(sbatch --parsable --cpus-per-task="${limits[0]}" --time="${limits[1]}" \
            --mem="${limits[2]}" --export="ALL,VGP_DATA_ROOT,VGP_RUN_ID,VGP_ENVIRONMENT_CAPTURE" \
            "$ROOT/analysis/slurm/vgp_10_pilot/pair_stage.sh" "$stage" "$selection_id")
        if [[ -n $dependency ]]; then
            command=(sbatch --parsable --dependency="afterok:$dependency" \
                --cpus-per-task="${limits[0]}" --time="${limits[1]}" --mem="${limits[2]}" \
                --export="ALL,VGP_DATA_ROOT,VGP_RUN_ID,VGP_ENVIRONMENT_CAPTURE" \
                "$ROOT/analysis/slurm/vgp_10_pilot/pair_stage.sh" "$stage" "$selection_id")
        fi
        if [[ $mode == --submit ]]; then
            dependency=$("${command[@]}")
        else
            submit "${command[@]}"
            dependency="<${selection_id}-${stage}-jobid>"
        fi
    done
    consensus_dependency=$dependency
    readarray -t limits < <(read_limits "$resource" psmc)
    command=(sbatch --parsable --dependency="afterok:$consensus_dependency" --array=0-200%20 \
        --cpus-per-task="${limits[0]}" --time="${limits[1]}" --mem="${limits[2]}" \
        --export="ALL,VGP_DATA_ROOT,VGP_RUN_ID,VGP_ENVIRONMENT_CAPTURE" \
        "$ROOT/analysis/slurm/vgp_10_pilot/psmc_array.sh" "$selection_id")
    if [[ $mode == --submit ]]; then
        psmc_dependency=$("${command[@]}")
    else
        submit "${command[@]}"
        psmc_dependency="<${selection_id}-psmc-array-jobid>"
    fi
    readarray -t limits < <(read_limits "$resource" psmc_finalize)
    command=(sbatch --parsable --dependency="afterok:$psmc_dependency" \
        --cpus-per-task="${limits[0]}" --time="${limits[1]}" --mem="${limits[2]}" \
        --export="ALL,VGP_DATA_ROOT,VGP_RUN_ID,VGP_ENVIRONMENT_CAPTURE" \
        "$ROOT/analysis/slurm/vgp_10_pilot/psmc_finalize.sh" "$selection_id")
    submit "${command[@]}"

    if python3 - "$VGP_DATA_ROOT/pilot/inputs/$selection_id/input-manifest.json" <<'PY'
import json,sys
raise SystemExit(0 if json.load(open(sys.argv[1])).get("annotation") is not None else 1)
PY
    then
        readarray -t limits < <(read_limits "$resource" annotation)
        command=(sbatch --parsable --dependency="afterok:$consensus_dependency" \
            --cpus-per-task="${limits[0]}" --time="${limits[1]}" --mem="${limits[2]}" \
            --export="ALL,VGP_DATA_ROOT,VGP_RUN_ID,VGP_ENVIRONMENT_CAPTURE" \
            "$ROOT/analysis/slurm/vgp_10_pilot/annotation_stage.sh" "$selection_id")
        submit "${command[@]}"
    fi
done
