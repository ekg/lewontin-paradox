#!/usr/bin/env bash
# Independently recompute a frozen primary's whole-assembly SweepGA mapping.
set -euo pipefail

selection_id=${1:?usage: independent_mapping.sh SELECTION_ID}
[[ $selection_id =~ ^P(0[1-9]|10)$ ]] || { echo "invalid selection id" >&2; exit 2; }
VGP_REPO_ROOT=${SLURM_SUBMIT_DIR:?submit from the repository root}
: "${VGP_DATA_ROOT:?VGP_DATA_ROOT is required}"
: "${VGP_ENVIRONMENT_CAPTURE:?VGP_ENVIRONMENT_CAPTURE is required}"
: "${VGP_FASTGA_AMENDMENT:?VGP_FASTGA_AMENDMENT is required}"
: "${VGP_NODE_LOCAL_BASE:=/scratch}"
[[ $VGP_DATA_ROOT == /moosefs/erikg/vgp ]] || {
    echo "independent recomputation must use the canonical VGP root" >&2
    exit 2
}

input_dir="$VGP_DATA_ROOT/pilot/inputs/$selection_id"
output_dir="$VGP_DATA_ROOT/pilot/independent/$selection_id/mapping"
mkdir -p "$output_dir" "$VGP_DATA_ROOT/pilot/logs"

readarray -t verified < <(python3 - "$VGP_ENVIRONMENT_CAPTURE" \
    "$VGP_FASTGA_AMENDMENT" <<'PY'
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
capture_sha256=${verified[1]}
amendment_sha256=${verified[2]}

[[ -d $VGP_NODE_LOCAL_BASE && -w $VGP_NODE_LOCAL_BASE ]] || {
    echo "node-local scratch is unavailable" >&2; exit 2;
}
case $(stat -f -c %T -- "$VGP_NODE_LOCAL_BASE") in
    nfs|nfs4|fuse*|lustre|gpfs|ceph) echo "scratch is not node-local" >&2; exit 2 ;;
esac
scratch=$(mktemp -d -- "$VGP_NODE_LOCAL_BASE/vgp-independent-$selection_id-${SLURM_JOB_ID:?}-XXXXXX")
cleanup() {
    status=$?
    [[ $scratch == "$VGP_NODE_LOCAL_BASE/vgp-independent-$selection_id-${SLURM_JOB_ID}-"* ]] && rm -rf -- "$scratch"
    exit "$status"
}
trap cleanup EXIT
mkdir -p "$scratch/inputs"
cp "$input_dir/h1.fa" "$scratch/inputs/h1.fa"
cp "$input_dir/h2.fa" "$scratch/inputs/h2.fa"
export TMPDIR="$scratch"
export TMP="$scratch"
export TEMP="$scratch"
cd -- "$scratch"
scratch_resolved=$(readlink -f -- "$scratch")
cwd_resolved=$(readlink -f -- /proc/$$/cwd)
[[ $cwd_resolved == "$scratch_resolved" ]] || {
    echo "batch cwd outside private node-local scratch: $cwd_resolved" >&2
    exit 70
}

started=$(date -u +%Y-%m-%dT%H:%M:%SZ)
"$sweepga" "$scratch/inputs/h2.fa" "$scratch/inputs/h1.fa" \
    --output-file "$scratch/h2_to_h1.native.1to1.paf" \
    --num-mappings 1:1 --scaffold-jump 0 --overlap 0 \
    --scoring log-length-ani --threads "${SLURM_CPUS_PER_TASK:-1}" \
    >"$scratch/sweepga.stdout" 2>"$scratch/sweepga.stderr" &
sweepga_pid=$!
guard_log="$scratch/fastga_scratch_snapshots.jsonl"
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
(( guard_failed == 0 )) || exit 70
(( sweepga_status == 0 )) || exit "$sweepga_status"
python3 "$VGP_REPO_ROOT/analysis/fastga_scratch_guard.py" finalize \
    --scratch "$scratch" --audit-jsonl "$guard_log" \
    --output "$scratch/fastga_scratch_contract.json"

final="$output_dir/${selection_id}.independent.${SLURM_JOB_ID}.native.1to1.paf"
partial="$final.partial"
cp "$scratch/h2_to_h1.native.1to1.paf" "$partial"
mv "$partial" "$final"
cp "$scratch/fastga_scratch_contract.json" "$final.fastga_scratch_contract.json.partial"
mv "$final.fastga_scratch_contract.json.partial" "$final.fastga_scratch_contract.json"
python3 - "$final" "$final.manifest.json" "$selection_id" "$started" \
    "$VGP_DATA_ROOT" "$capture_sha256" "$amendment_sha256" <<'PY'
import datetime,hashlib,json,os,sys
from pathlib import Path
p=Path(sys.argv[1])
manifest={
 "schema_version":"vgp-independent-sweepga-mapping-v1",
 "selection_id":sys.argv[3],"canonical_vgp_root":sys.argv[5],
 "job_id":os.environ["SLURM_JOB_ID"],"started_at_utc":sys.argv[4],
 "completed_at_utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),
 "command_contract":{"num_mappings":"1:1","scaffold_jump":0,"overlap":0,
                     "scoring":"log-length-ani","orientation":"H2_query_to_H1_reference"},
 "paf":{"path":str(p),"size_bytes":p.stat().st_size,
        "sha256":hashlib.sha256(p.read_bytes()).hexdigest(),
        "records":sum(1 for _ in p.open())},
 "environment_capture_sha256":sys.argv[6],
 "fastga_amendment_sha256":sys.argv[7],
 "independent_of_primary_output":True,
}
Path(sys.argv[2]).write_text(json.dumps(manifest,sort_keys=True,indent=2)+"\n")
PY
