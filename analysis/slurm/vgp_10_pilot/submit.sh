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
: "${VGP_SELECTION_IDS:=P01 P02 P03 P04 P05 P06 P08 P09 P10}"
: "${VGP_NODE_LOCAL_BASE:=/scratch}"
: "${VGP_SUBMISSION_MANIFEST:=$VGP_DATA_ROOT/pilot/manifests/$VGP_RUN_ID.submissions.tsv}"
: "${VGP_RESOURCE_PLAN:=$ROOT/analysis/vgp_real_pilot_resource_plan_v1.json}"
export VGP_DATA_ROOT VGP_RUN_ID VGP_ENVIRONMENT_CAPTURE VGP_NODE_LOCAL_BASE \
    VGP_SUBMISSION_MANIFEST VGP_RESOURCE_PLAN
mkdir -p "$VGP_DATA_ROOT/pilot/logs" "$(dirname -- "$VGP_SUBMISSION_MANIFEST")"

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
    local selection=$3
    python3 - "$resource" "$stage" "$selection" "$VGP_RESOURCE_PLAN" <<'PY'
import json,sys
v=json.load(open(sys.argv[1]))
stage=v.get("stages",{}).get(sys.argv[2])
if not isinstance(stage,dict): raise SystemExit(f"missing stage-specific estimate: {sys.argv[2]}")
try:
    plan=json.load(open(sys.argv[4]))
except FileNotFoundError:
    plan={}
override=plan.get("pairs",{}).get(sys.argv[3],{}).get("stages",{}).get(sys.argv[2],{})
stage={**stage,**override}
for key in ("cpus_per_task","slurm_time","slurm_mem","scratch_bytes_high"):
    if not stage.get(key): raise SystemExit(f"missing scheduler-specific limit: {sys.argv[2]}.{key}")
    print(stage[key])
print(override.get("slurm_partition",plan.get("stage_partitions",{}).get(sys.argv[2],"workers")))
print(override.get("slurm_exclude",plan.get("stage_excludes",{}).get(sys.argv[2],"")))
PY
}

record_submission() {
    local job_id=$1 selection=$2 stage=$3 dependency=$4
    python3 - "$VGP_SUBMISSION_MANIFEST" "$job_id" "$selection" "$stage" "$dependency" <<'PY'
import csv,datetime,os,sys
from pathlib import Path
path=Path(sys.argv[1]); exists=path.exists()
with path.open("a",newline="") as handle:
    fields=("authorization_id","run_id","canonical_vgp_root","selection_id","stage",
            "job_id","dependency","submitted_at_utc","resource_plan")
    writer=csv.DictWriter(handle,fieldnames=fields,delimiter="\t",lineterminator="\n")
    if not exists: writer.writeheader()
    writer.writerow({
        "authorization_id":"vgp10-auth-20260718-v2","run_id":os.environ["VGP_RUN_ID"],
        "canonical_vgp_root":os.environ["VGP_DATA_ROOT"],"selection_id":sys.argv[3],
        "stage":sys.argv[4],"job_id":sys.argv[2],"dependency":sys.argv[5],
        "submitted_at_utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "resource_plan":os.environ["VGP_RESOURCE_PLAN"],
    })
PY
}

submit_job() {
    local selection=$1 stage=$2 dependency=$3
    shift 3
    local job_id
    job_id=$("$@")
    record_submission "$job_id" "$selection" "$stage" "$dependency"
    printf '%s\n' "$job_id"
}

declare -A seen=()
for selection_id in $VGP_SELECTION_IDS; do
    require_pattern='^P(0[1-9]|10)$'
    [[ $selection_id =~ $require_pattern && $selection_id != P07 ]] || {
        echo "ERROR: remaining-primary selection must be P01..P10 except completed P07: $selection_id" >&2
        exit 2
    }
    [[ -z ${seen[$selection_id]:-} ]] || { echo "ERROR: duplicate selection: $selection_id" >&2; exit 2; }
    seen[$selection_id]=1
    resource="$VGP_DATA_ROOT/pilot/inputs/$selection_id/resources.json"
    [[ -f $resource ]] || {
        echo "ERROR: measured resource estimate absent: $resource" >&2
        exit 2
    }
    dependency=
    for stage in preflight mapping impg variants consensus; do
        sentinel="$VGP_DATA_ROOT/pilot/runs/$VGP_RUN_ID/$selection_id/$stage/.complete.json"
        if [[ -f $sentinel ]]; then
            echo "RESUME-SUBMIT: $selection_id $stage already complete" >&2
            continue
        fi
        readarray -t limits < <(read_limits "$resource" "$stage" "$selection_id")
        log_prefix="$VGP_DATA_ROOT/pilot/logs/$VGP_RUN_ID-$selection_id-$stage"
        if [[ $stage == mapping ]]; then
            stage_script="$ROOT/analysis/slurm/vgp_10_pilot/mapping_stage.sh"
            stage_args=("$selection_id")
        else
            stage_script="$ROOT/analysis/slurm/vgp_10_pilot/pair_stage.sh"
            stage_args=("$stage" "$selection_id")
        fi
        node_args=()
        [[ -n ${limits[5]:-} ]] && node_args=(--exclude="${limits[5]}")
        command=(sbatch --parsable --cpus-per-task="${limits[0]}" --time="${limits[1]}" \
            --mem="${limits[2]}" --partition="${limits[4]}" \
            "${node_args[@]}" \
            --job-name="vgp10-$selection_id-$stage" --output="$log_prefix-%j.out" --error="$log_prefix-%j.err" \
            --export="ALL,VGP_DATA_ROOT,VGP_RUN_ID,VGP_ENVIRONMENT_CAPTURE,VGP_NODE_LOCAL_BASE,VGP_RESOURCE_PLAN" \
            "$stage_script" "${stage_args[@]}")
        if [[ -n $dependency ]]; then
            command=(sbatch --parsable --dependency="afterok:$dependency" \
                --cpus-per-task="${limits[0]}" --time="${limits[1]}" --mem="${limits[2]}" \
                --partition="${limits[4]}" \
                "${node_args[@]}" \
                --job-name="vgp10-$selection_id-$stage" --output="$log_prefix-%j.out" --error="$log_prefix-%j.err" \
                --export="ALL,VGP_DATA_ROOT,VGP_RUN_ID,VGP_ENVIRONMENT_CAPTURE,VGP_NODE_LOCAL_BASE,VGP_RESOURCE_PLAN" \
                "$stage_script" "${stage_args[@]}")
        fi
        if [[ $mode == --submit ]]; then
            dependency=$(submit_job "$selection_id" "$stage" "$dependency" "${command[@]}")
        else
            submit "${command[@]}"
            dependency="<${selection_id}-${stage}-jobid>"
        fi
    done
    consensus_dependency=$dependency
    psmc_finalize_sentinel="$VGP_DATA_ROOT/pilot/runs/$VGP_RUN_ID/$selection_id/psmc/finalize/.complete.json"
    if [[ -f $psmc_finalize_sentinel ]]; then
        echo "RESUME-SUBMIT: $selection_id PSMC finalize already complete" >&2
    else
      readarray -t limits < <(read_limits "$resource" psmc "$selection_id")
      command=(sbatch --parsable --array=0-200%20 \
        --cpus-per-task="${limits[0]}" --time="${limits[1]}" --mem="${limits[2]}" --partition="${limits[4]}" \
        --job-name="vgp10-$selection_id-psmc" \
        --output="$VGP_DATA_ROOT/pilot/logs/$VGP_RUN_ID-$selection_id-psmc-%A_%a.out" \
        --error="$VGP_DATA_ROOT/pilot/logs/$VGP_RUN_ID-$selection_id-psmc-%A_%a.err" \
        --export="ALL,VGP_DATA_ROOT,VGP_RUN_ID,VGP_ENVIRONMENT_CAPTURE,VGP_NODE_LOCAL_BASE,VGP_RESOURCE_PLAN" \
        "$ROOT/analysis/slurm/vgp_10_pilot/psmc_array.sh" "$selection_id")
      if [[ -n $consensus_dependency ]]; then
        command=(sbatch --parsable --dependency="afterok:$consensus_dependency" "${command[@]:2}")
      fi
    if [[ $mode == --submit ]]; then
        psmc_dependency=$(submit_job "$selection_id" psmc "$consensus_dependency" "${command[@]}")
    else
        submit "${command[@]}"
        psmc_dependency="<${selection_id}-psmc-array-jobid>"
    fi
    readarray -t limits < <(read_limits "$resource" psmc_finalize "$selection_id")
    command=(sbatch --parsable --dependency="afterok:$psmc_dependency" \
        --cpus-per-task="${limits[0]}" --time="${limits[1]}" --mem="${limits[2]}" --partition="${limits[4]}" \
        --job-name="vgp10-$selection_id-finalize" \
        --output="$VGP_DATA_ROOT/pilot/logs/$VGP_RUN_ID-$selection_id-finalize-%j.out" \
        --error="$VGP_DATA_ROOT/pilot/logs/$VGP_RUN_ID-$selection_id-finalize-%j.err" \
        --export="ALL,VGP_DATA_ROOT,VGP_RUN_ID,VGP_ENVIRONMENT_CAPTURE,VGP_NODE_LOCAL_BASE,VGP_RESOURCE_PLAN" \
        "$ROOT/analysis/slurm/vgp_10_pilot/psmc_finalize.sh" "$selection_id")
    if [[ $mode == --submit ]]; then
        submit_job "$selection_id" psmc_finalize "$psmc_dependency" "${command[@]}" >/dev/null
    else
        submit "${command[@]}"
    fi
    fi

    if python3 - "$VGP_DATA_ROOT/pilot/inputs/$selection_id/input-manifest.json" <<'PY'
import json,sys
raise SystemExit(0 if json.load(open(sys.argv[1])).get("annotation") is not None else 1)
PY
    then
        annotation_sentinel="$VGP_DATA_ROOT/pilot/runs/$VGP_RUN_ID/$selection_id/annotation/.complete.json"
        if [[ -f $annotation_sentinel ]]; then
            echo "RESUME-SUBMIT: $selection_id annotation already complete" >&2
            continue
        fi
        readarray -t limits < <(read_limits "$resource" annotation "$selection_id")
        command=(sbatch --parsable \
            --cpus-per-task="${limits[0]}" --time="${limits[1]}" --mem="${limits[2]}" --partition="${limits[4]}" \
            --job-name="vgp10-$selection_id-annotation" \
            --output="$VGP_DATA_ROOT/pilot/logs/$VGP_RUN_ID-$selection_id-annotation-%j.out" \
            --error="$VGP_DATA_ROOT/pilot/logs/$VGP_RUN_ID-$selection_id-annotation-%j.err" \
            --export="ALL,VGP_DATA_ROOT,VGP_RUN_ID,VGP_ENVIRONMENT_CAPTURE,VGP_NODE_LOCAL_BASE,VGP_RESOURCE_PLAN" \
            "$ROOT/analysis/slurm/vgp_10_pilot/annotation_stage.sh" "$selection_id")
        if [[ -n $consensus_dependency ]]; then
            command=(sbatch --parsable --dependency="afterok:$consensus_dependency" "${command[@]:2}")
        fi
        if [[ $mode == --submit ]]; then
            submit_job "$selection_id" annotation "$consensus_dependency" "${command[@]}" >/dev/null
        else
            submit "${command[@]}"
        fi
    fi
done
