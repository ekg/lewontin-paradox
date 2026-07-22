#!/usr/bin/env bash
# Minimal canary-style whole-assembly mapping entry point for scale-out.
set -euo pipefail

VGP_SELECTION_ID=${1:?usage: mapping_stage.sh SELECTION_ID}
[[ $VGP_SELECTION_ID =~ ^P(0[1-9]|10)$ ]] || { echo "invalid selection" >&2; exit 2; }
VGP_REPO_ROOT=${SLURM_SUBMIT_DIR:?submit from the repository root}
: "${VGP_DATA_ROOT:?VGP_DATA_ROOT is required}"
: "${VGP_RUN_ID:?VGP_RUN_ID is required}"
: "${VGP_ENVIRONMENT_CAPTURE:?VGP_ENVIRONMENT_CAPTURE is required}"
: "${VGP_FASTGA_AMENDMENT:=$VGP_REPO_ROOT/analysis/vgp_real_pilot_fastga_amendment_v1.json}"
: "${VGP_NODE_LOCAL_BASE:=/scratch}"
: "${VGP_RESOURCE_PLAN:=}"
export PYTHONPATH=$VGP_REPO_ROOT

input_dir="$VGP_DATA_ROOT/pilot/inputs/$VGP_SELECTION_ID"
pair_run="$VGP_DATA_ROOT/pilot/runs/$VGP_RUN_ID/$VGP_SELECTION_ID"
final="$pair_run/mapping"
telemetry="$pair_run/telemetry"
if [[ -f $final/.complete.json ]]; then
    echo "RESUME: $VGP_SELECTION_ID mapping already complete"
    exit 0
fi
[[ -f $pair_run/preflight/.complete.json ]] || { echo "preflight sentinel absent" >&2; exit 2; }
[[ ! -e $final ]] || { echo "mapping final exists without sentinel" >&2; exit 2; }
mkdir -p "$telemetry"

readarray -t verified < <(python3 - "$VGP_ENVIRONMENT_CAPTURE" "$VGP_FASTGA_AMENDMENT" <<'PY'
import hashlib,json,sys
from pathlib import Path
capture=json.load(open(sys.argv[1])); amendment=json.load(open(sys.argv[2]))
def verify(row,name):
    path=Path(row.get("path",""))
    if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest()!=row.get("sha256"):
        raise SystemExit(f"digest mismatch: {name}")
    return str(path)
rows={row["name"]:row for row in capture["executables"]}
print(verify(rows["sweepga"],"sweepga"))
for name,row in amendment["companions"].items(): verify(row,name)
print(hashlib.sha256(Path(sys.argv[1]).read_bytes()).hexdigest())
print(hashlib.sha256(Path(sys.argv[2]).read_bytes()).hexdigest())
PY
)
sweepga=${verified[0]}

required_bytes=$(python3 - "$input_dir/resources.json" "$VGP_RESOURCE_PLAN" \
    "$VGP_SELECTION_ID" <<'PY'
import json,sys
stage=json.load(open(sys.argv[1]))["stages"]["mapping"]
if sys.argv[2]:
    plan=json.load(open(sys.argv[2]))
    stage={**stage,**plan.get("pairs",{}).get(sys.argv[3],{}).get("stages",{}).get("mapping",{})}
print(int(stage["scratch_bytes_high"]))
PY
)
[[ -d $VGP_NODE_LOCAL_BASE && -w $VGP_NODE_LOCAL_BASE ]] || { echo "scratch unavailable" >&2; exit 2; }
case $(stat -f -c %T -- "$VGP_NODE_LOCAL_BASE") in
    nfs|nfs4|fuse*|lustre|gpfs|ceph) echo "scratch is not node-local" >&2; exit 2 ;;
esac
available_bytes=$(df -PB1 -- "$VGP_NODE_LOCAL_BASE" | awk 'NR==2 {print $4}')
(( available_bytes >= required_bytes )) || { echo "insufficient measured scratch" >&2; exit 2; }
scratch=$(mktemp -d -- "$VGP_NODE_LOCAL_BASE/vgp-map-$VGP_SELECTION_ID-${SLURM_JOB_ID:?}-XXXXXX")
partial="$scratch/mapping.partial"
mkdir -p "$partial" "$scratch/inputs"
export TMPDIR="$scratch"
export TMP="$scratch"
export TEMP="$scratch"
started=$(date +%s)

finish() {
    status=$?
    if (( status != 0 )); then
        python3 - "$telemetry/mapping.${SLURM_JOB_ID}.json" "$partial" "$status" \
            "$started" "$required_bytes" "$available_bytes" "$VGP_DATA_ROOT" "$VGP_SELECTION_ID" <<'PY' || true
import json,sys,time
from pathlib import Path
out,partial=Path(sys.argv[1]),Path(sys.argv[2]); diagnostics={}
if partial.exists():
    for path in partial.rglob("*"):
        if path.is_file() and path.suffix in {".stderr",".stdout",".log"}:
            diagnostics[str(path.relative_to(partial))]=path.read_text(errors="replace")[-32768:]
out.write_text(json.dumps({
 "selection_id":sys.argv[8],"stage":"mapping","job_id":out.stem.split(".")[-1],
 "disposition":"failure","exit_status":int(sys.argv[3]),"started_epoch":int(sys.argv[4]),
 "ended_epoch":int(time.time()),"scratch_required_bytes":int(sys.argv[5]),
 "scratch_available_bytes_at_start":int(sys.argv[6]),"canonical_vgp_root":sys.argv[7],
 "diagnostic_tails":diagnostics,
},sort_keys=True)+"\n")
PY
    fi
    [[ $scratch == "$VGP_NODE_LOCAL_BASE/vgp-map-$VGP_SELECTION_ID-${SLURM_JOB_ID}-"* ]] && rm -rf -- "$scratch"
    exit "$status"
}
trap finish EXIT

cp "$input_dir/h1.fa" "$scratch/inputs/h1.fa"
cp "$input_dir/h2.fa" "$scratch/inputs/h2.fa"
cd -- "$scratch"
scratch_resolved=$(readlink -f -- "$scratch")
cwd_resolved=$(readlink -f -- /proc/$$/cwd)
[[ $cwd_resolved == "$scratch_resolved" ]] || {
    echo "batch cwd outside private node-local scratch: $cwd_resolved" >&2
    exit 70
}
threads=${SLURM_CPUS_PER_TASK:-1}
"$sweepga" "$scratch/inputs/h2.fa" "$scratch/inputs/h1.fa" \
    --output-file "$partial/h2_to_h1.native.1to1.paf" \
    --num-mappings 1:1 --scaffold-jump 0 --overlap 0 \
    --scoring log-length-ani --threads "$threads" \
    >"$partial/sweepga.stdout" 2>"$partial/sweepga.stderr" &
sweepga_pid=$!
guard_log="$partial/fastga_scratch_snapshots.jsonl"
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
if wait "$sweepga_pid"; then
    sweepga_status=0
else
    sweepga_status=$?
fi
(( guard_failed == 0 )) || {
    echo "hard infrastructure error: FastGA escaped private node-local scratch" >&2
    exit 70
}
(( sweepga_status == 0 )) || exit "$sweepga_status"
python3 "$VGP_REPO_ROOT/analysis/fastga_scratch_guard.py" finalize \
    --scratch "$scratch" --audit-jsonl "$guard_log" \
    --output "$partial/fastga_scratch_contract.json"
python3 -m analysis.vgp_10_pilot enforce-paf \
    "$partial/h2_to_h1.native.1to1.paf" "$partial/h2_to_h1.1to1.paf" \
    >"$partial/exact_multiplicity_filter.json"
python3 -m analysis.vgp_10_pilot audit-paf \
    "$partial/h2_to_h1.1to1.paf" "$scratch/inputs/h1.fa" "$scratch/inputs/h2.fa" \
    >"$partial/multiplicity.json"
python3 - "$partial/h2_to_h1.1to1.paf" "$scratch/inputs/h1.fa" \
    "$scratch/inputs/h2.fa" "$input_dir/h1_universe.bed" "$partial" <<'PY'
import sys
from pathlib import Path
from analysis.vgp_10_pilot import (
 low_complexity_intervals,non_acgt_intervals,paf_h1_intervals,parse_fasta,parse_paf,
 project_h2_non_acgt_to_h1,read_bed,subtract_intervals,write_bed,
)
paf,h1,h2,universe,out=sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],Path(sys.argv[5])
records=parse_paf(paf); one_to_one=paf_h1_intervals(records); h1_sequences=parse_fasta(h1)
write_bed(out/"h1.1to1.bed",one_to_one)
write_bed(out/"not_1to1.bed",subtract_intervals(read_bed(universe),one_to_one))
write_bed(out/"h1_gap_or_N.bed",non_acgt_intervals(h1_sequences))
write_bed(out/"h2_gap_or_N.bed",project_h2_non_acgt_to_h1(records,parse_fasta(h2)))
write_bed(out/"repeat_or_low_complexity_primary.bed",low_complexity_intervals(h1_sequences))
PY
python3 - "$partial/mapping_execution.json" "$VGP_DATA_ROOT" "$VGP_SELECTION_ID" \
    "$VGP_ENVIRONMENT_CAPTURE" "${verified[1]}" "$VGP_FASTGA_AMENDMENT" "${verified[2]}" <<'PY'
import json,os,sys
from pathlib import Path
Path(sys.argv[1]).write_text(json.dumps({
 "schema_version":os.environ.get("VGP_MAPPING_SCHEMA_VERSION","vgp-real-pilot-mapping-execution-v1"),
 "task_id":os.environ.get("VGP_TASK_ID","run-vgp-real-pilot"),
 "authorization_id":"vgp10-auth-20260718-v2","canonical_vgp_root":sys.argv[2],
 "selection_id":sys.argv[3],"slurm_job_id":os.environ["SLURM_JOB_ID"],
 "required_option":"--num-mappings 1:1","environment_capture":{"path":sys.argv[4],"sha256":sys.argv[5]},
 "fastga_amendment":{"path":sys.argv[6],"sha256":sys.argv[7]},
},sort_keys=True)+"\n")
PY
python3 - "$partial" "$final" "$VGP_SELECTION_ID" "$SLURM_JOB_ID" "$VGP_DATA_ROOT" <<'PY'
import sys
from pathlib import Path
from analysis.vgp_10_pilot import promote_stage_cross_filesystem
promote_stage_cross_filesystem(Path(sys.argv[1]),Path(sys.argv[2]),{
 "selection_id":sys.argv[3],"stage":"mapping","atomic_promotion":True,
 "slurm_job_id":sys.argv[4],"canonical_vgp_root":sys.argv[5],
},sys.argv[4])
PY
python3 - "$telemetry/mapping.${SLURM_JOB_ID}.json" "$started" "$required_bytes" \
    "$available_bytes" "$VGP_DATA_ROOT" "$VGP_SELECTION_ID" <<'PY'
import json,os,sys,time
from pathlib import Path
Path(sys.argv[1]).write_text(json.dumps({
 "selection_id":sys.argv[6],"stage":"mapping","job_id":os.environ["SLURM_JOB_ID"],
 "disposition":"success","started_epoch":int(sys.argv[2]),"ended_epoch":int(time.time()),
 "scratch_required_bytes":int(sys.argv[3]),"scratch_available_bytes_at_start":int(sys.argv[4]),
 "canonical_vgp_root":sys.argv[5],"node_local_scratch_base":os.environ["VGP_NODE_LOCAL_BASE"],
},sort_keys=True)+"\n")
PY
