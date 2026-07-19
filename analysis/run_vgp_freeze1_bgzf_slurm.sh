#!/usr/bin/env bash
set -euo pipefail

repo_root=$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
readonly repo_root
: "${VGP_BGZF_DERIVED_ROOT:=/moosefs/erikg/vgp/derived/freeze1-bgzf}"
: "${VGP_BGZF_MIRROR_MANIFEST:=$repo_root/analysis/vgp_freeze1_mirror_manifest.tsv}"
: "${VGP_BGZF_INVENTORY:=$repo_root/analysis/vgp_freeze1_bgzf_inventory.tsv}"
: "${VGP_BGZF_PROFILE:=$VGP_BGZF_DERIVED_ROOT/environment/profile}"
: "${VGP_BGZF_PARTITION:=workers}"
: "${VGP_BGZF_RETRY_LIMIT:=3}"
: "${VGP_BGZF_NODE_LOCAL_BASE:=/scratch}"

channels="$repo_root/analysis/guix/channels.scm"
manifest="$repo_root/analysis/guix/vgp-freeze1-bgzf-manifest.scm"
pipeline="$repo_root/analysis/build_vgp_freeze1_bgzf.py"

usage() {
    printf 'usage: %s inventory|realize|submit|retry|finalize|status\n' "$0" >&2
    exit 64
}

inventory() {
    python3 "$pipeline" inventory \
        --mirror-manifest "$VGP_BGZF_MIRROR_MANIFEST" \
        --derived-root "$VGP_BGZF_DERIVED_ROOT" \
        --output "$VGP_BGZF_INVENTORY"
}

realize() {
    local environment_id marker
    mkdir -p "$(dirname -- "$VGP_BGZF_PROFILE")"
    environment_id=$(sha256sum "$channels" "$manifest" | sha256sum | cut -d ' ' -f 1)
    marker="${VGP_BGZF_PROFILE}.environment-id"
    if [[ ! -x "$VGP_BGZF_PROFILE/bin/bgzip" || ! -x "$VGP_BGZF_PROFILE/bin/samtools" || ! -f $marker || $(<"$marker") != "$environment_id" ]]; then
        guix time-machine -C "$channels" -- package \
            --profile="$VGP_BGZF_PROFILE" --manifest="$manifest"
        printf '%s\n' "$environment_id" >"$marker"
    fi
    "$VGP_BGZF_PROFILE/bin/bgzip" --version
    "$VGP_BGZF_PROFILE/bin/samtools" --version | head -n 1
}

indices_for_class() {
    local class=$1 mode=$2
    python3 - "$VGP_BGZF_INVENTORY" "$VGP_BGZF_DERIVED_ROOT" "$class" "$mode" "$VGP_BGZF_RETRY_LIMIT" <<'PY'
import csv,json,sys
from pathlib import Path
inventory,root,wanted,mode,retry_limit=sys.argv[1:]
chosen=[]
with open(inventory,newline='') as handle:
    for row in csv.DictReader(handle,delimiter='\t'):
        if row['resource_class'] != wanted:
            continue
        status=Path(root)/'status'/f"{int(row['task_index']):06d}.{row['accession_version']}.json"
        data={}
        if status.is_file():
            try: data=json.loads(status.read_text())
            except json.JSONDecodeError: pass
        if data.get('state') in {'converted','reused'}:
            continue
        attempts=int(data.get('attempts',0))
        if mode == 'retry' and (not data or attempts >= int(retry_limit)):
            continue
        chosen.append(row['task_index'])
print(','.join(chosen))
PY
}

resource_value() {
    local class=$1 field=$2
    python3 - "$VGP_BGZF_INVENTORY" "$class" "$field" <<'PY'
import csv,sys
with open(sys.argv[1],newline='') as h:
    for row in csv.DictReader(h,delimiter='\t'):
        if row['resource_class']==sys.argv[2]:
            print(row[sys.argv[3]]); break
PY
}

submit_arrays() {
    local mode=$1 class indices cpus memory walltime parallel job_id release_id control
    realize >/dev/null
    mkdir -p "$VGP_BGZF_DERIVED_ROOT/logs"
    release_id=$(
        sha256sum "$pipeline" "$repo_root/analysis/slurm/vgp_freeze1_bgzf_worker.sh" "$VGP_BGZF_INVENTORY" \
            | sha256sum | cut -d ' ' -f 1
    )
    control="$VGP_BGZF_DERIVED_ROOT/control/$release_id"
    mkdir -p "$control"
    install -m 0555 "$pipeline" "$control/build_vgp_freeze1_bgzf.py"
    install -m 0555 "$repo_root/analysis/slurm/vgp_freeze1_bgzf_worker.sh" "$control/worker.sh"
    install -m 0444 "$VGP_BGZF_INVENTORY" "$control/inventory.tsv"
    for class in tiny small medium large; do
        indices=$(indices_for_class "$class" "$mode")
        [[ -n $indices ]] || continue
        cpus=$(resource_value "$class" cpus)
        memory=$(resource_value "$class" memory_mb)
        walltime=$(resource_value "$class" walltime)
        case "$class" in tiny) parallel=6;; small) parallel=8;; medium) parallel=4;; large) parallel=2;; esac
        job_id=$(sbatch --parsable \
            --job-name="vgp-bgzf-$class" --partition="$VGP_BGZF_PARTITION" \
            --array="$indices%$parallel" --cpus-per-task="$cpus" --mem="${memory}M" \
            --time="$walltime" --nice=500 \
            --output="$VGP_BGZF_DERIVED_ROOT/logs/%A_%a.out" \
            --error="$VGP_BGZF_DERIVED_ROOT/logs/%A_%a.err" \
            --export="ALL,VGP_NODE_LOCAL_BASE=$VGP_BGZF_NODE_LOCAL_BASE,VGP_BGZF_PROFILE=$VGP_BGZF_PROFILE,VGP_BGZF_INVENTORY=$control/inventory.tsv,VGP_BGZF_DERIVED_ROOT=$VGP_BGZF_DERIVED_ROOT,VGP_BGZF_PIPELINE=$control/build_vgp_freeze1_bgzf.py" \
            "$control/worker.sh")
        printf '%s\t%s\t%s\n' "$class" "$job_id" "$indices"
    done
}

finalize() {
    python3 "$pipeline" finalize --inventory "$VGP_BGZF_INVENTORY" \
        --derived-root "$VGP_BGZF_DERIVED_ROOT" \
        --shared-manifest "$repo_root/analysis/vgp_freeze1_bgzf_manifest.tsv" \
        --summary "$repo_root/analysis/vgp_freeze1_bgzf_summary.json"
}

status_report() {
    python3 - "$VGP_BGZF_INVENTORY" "$VGP_BGZF_DERIVED_ROOT" <<'PY'
import csv,json,collections,sys
from pathlib import Path
rows=list(csv.DictReader(open(sys.argv[1]),delimiter='\t')); root=Path(sys.argv[2]); states=collections.Counter(); reasons=collections.Counter()
for row in rows:
 p=root/'status'/f"{int(row['task_index']):06d}.{row['accession_version']}.json"
 if not p.is_file(): states['pending']+=1; continue
 try: d=json.loads(p.read_text())
 except json.JSONDecodeError: states['invalid_status']+=1; continue
 states[d.get('state','unknown')]+=1; reasons[d.get('reason_code','')]+=1
print(json.dumps({'total':len(rows),'states':states,'reasons':reasons},indent=2,sort_keys=True))
PY
}

[[ $# -eq 1 ]] || usage
case "$1" in
    inventory) inventory ;;
    realize) realize ;;
    submit) inventory; submit_arrays submit ;;
    retry) submit_arrays retry ;;
    finalize) finalize ;;
    status) status_report ;;
    *) usage ;;
esac
