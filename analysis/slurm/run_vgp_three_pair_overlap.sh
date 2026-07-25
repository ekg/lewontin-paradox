#!/usr/bin/env bash
#SBATCH --job-name=vgp3-P03-overlap
#SBATCH --partition=highmem
#SBATCH --exclude=octopus11
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=12:00:00
set -euo pipefail

: "${SLURM_JOB_ID:?run under Slurm}"
REPO_ROOT=${SLURM_SUBMIT_DIR:?submit from repository root}
SELECTION="$REPO_ROOT/analysis/vgp_three_pair_selection_v1.json"
CAPTURE="$REPO_ROOT/analysis/guix/vgp_10_pilot/realization.json"
AMENDMENT="$REPO_ROOT/analysis/vgp_real_pilot_fastga_amendment_v1.json"
RUN_ID=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["run_id"])' "$SELECTION")
MANIFEST=/moosefs/erikg/vgp/pilot/inputs/P03/input-manifest.json
DURABLE=/moosefs/erikg/vgp/pilot/three-pair/$RUN_ID/P03/backend-comparison
[[ ! -e $DURABLE ]] || { echo "comparison target already exists" >&2; exit 2; }
scratch=$(mktemp -d -- "/scratch/vgp3-overlap-${SLURM_JOB_ID}-XXXXXX")
trap 'status=$?; [[ $scratch == /scratch/vgp3-overlap-${SLURM_JOB_ID}-* ]] && rm -rf -- "$scratch"; exit $status' EXIT
export TMPDIR="$scratch" TMP="$scratch" TEMP="$scratch" PYTHONPATH="$REPO_ROOT"
cd -- "$scratch"
started=$(date +%s)

readarray -t tools < <(python3 - "$CAPTURE" "$AMENDMENT" <<'PY'
import hashlib,json,sys
from pathlib import Path
capture=json.load(open(sys.argv[1])); amendment=json.load(open(sys.argv[2]))
def verified(row):
    path=Path(row["path"])
    if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest()!=row["sha256"]:
        raise SystemExit(f"digest mismatch: {path}")
    return str(path)
rows={row["name"]:row for row in capture["executables"]}
print(verified(rows["sweepga"]))
for name in ("FAtoGDB","GIXmake","wfmash"):
    print(verified(amendment["companions"][name]))
PY
)
sweepga=${tools[0]}; fatogdb=${tools[1]}; gixmake=${tools[2]}

python3 - "$MANIFEST" "$scratch" <<'PY'
import hashlib,json,sys
from pathlib import Path
from analysis.vgp_10_pilot import sha256_file
manifest=json.load(open(sys.argv[1])); out=Path(sys.argv[2]); limit=20_000_000
audit={"schema_version":"vgp-three-pair-controlled-subset-v1","subset_bp":limit,"roles":{}}
for role,filename in (("h1_fasta","h1.subset.fa"),("h2_fasta","h2.subset.fa")):
    asset=manifest["assets"][role]; source=Path(asset["path"])
    if source.stat().st_size != asset["size_bytes"] or sha256_file(source) != asset["sha256"]:
        raise SystemExit(f"immutable source mismatch: {role}")
    expected=asset["sequence_dictionary"][0]
    parts=[]; name=None; collected=0
    with source.open() as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None: break
                name=line[1:].split()[0]
            elif name is not None and collected < limit:
                part=line.strip().upper(); parts.append(part); collected+=len(part)
            elif name is not None:
                break
    sequence="".join(parts)[:limit]
    if name != expected["name"] or len(sequence) != limit:
        raise SystemExit(f"controlled first-sequence extraction failed: {role}")
    target=out/filename
    with target.open("w") as handle:
        handle.write(f">{name}\n")
        for offset in range(0,len(sequence),80): handle.write(sequence[offset:offset+80]+"\n")
    audit["roles"][role]={"source_sha256":asset["sha256"],"name":name,"start":0,"end":limit,
        "subset_sequence_sha256":hashlib.sha256(sequence.encode()).hexdigest(),"staged_fasta_sha256":sha256_file(target)}
(out/"subset.json").write_text(json.dumps(audit,sort_keys=True)+"\n")
PY

for fasta in "$scratch/h2.subset.fa" "$scratch/h1.subset.fa"; do
    "$fatogdb" "$fasta"
    "$gixmake" -T"${SLURM_CPUS_PER_TASK}" -P"$scratch" "${fasta%.fa}"
done
common=("$scratch/h2.subset.fa" "$scratch/h1.subset.fa" --num-mappings 1:1
    --scaffold-jump 0 --overlap 0 --scoring log-length-ani --threads "$SLURM_CPUS_PER_TASK")
"$sweepga" "${common[@]}" --output-file "$scratch/fastga.native.paf" \
    >"$scratch/fastga.stdout" 2>"$scratch/fastga.stderr"
"$sweepga" "${common[@]}" --aligner wfmash --map-pct-identity 90 --min-aln-length 25000 \
    --output-file "$scratch/wfmash.native.paf" \
    >"$scratch/wfmash.stdout" 2>"$scratch/wfmash.stderr"
for backend in fastga wfmash; do
    python3 -m analysis.vgp_10_pilot enforce-paf \
        "$scratch/$backend.native.paf" "$scratch/$backend.exact.paf" \
        >"$scratch/$backend.filter.json"
    python3 -m analysis.vgp_10_pilot audit-paf "$scratch/$backend.exact.paf" \
        "$scratch/h1.subset.fa" "$scratch/h2.subset.fa" >"$scratch/$backend.multiplicity.json"
done
python3 -m analysis.vgp_three_pair compare-pafs \
    "$scratch/fastga.exact.paf" "$scratch/wfmash.exact.paf" \
    "$scratch/h1.subset.fa" "$scratch/h2.subset.fa" "$scratch/comparison" \
    >"$scratch/comparison.stdout.json"
python3 - "$scratch/job.json" "$started" <<'PY'
import json,os,resource,sys,time
from pathlib import Path
u=resource.getrusage(resource.RUSAGE_CHILDREN)
Path(sys.argv[1]).write_text(json.dumps({
 "schema_version":"vgp-three-pair-overlap-job-v1","job_id":os.environ["SLURM_JOB_ID"],
 "node":os.environ.get("SLURMD_NODENAME"),"started_epoch":int(sys.argv[2]),"ended_epoch":int(time.time()),
 "maximum_child_rss_kib":u.ru_maxrss,"child_cpu_seconds":u.ru_utime+u.ru_stime,
 "private_node_local_scratch":True,
},sort_keys=True)+"\n")
PY
mkdir -p "$(dirname -- "$DURABLE")"
partial="$DURABLE.${SLURM_JOB_ID}.partial"
mkdir "$partial"
cp subset.json fastga.exact.paf wfmash.exact.paf fastga.filter.json wfmash.filter.json \
    fastga.multiplicity.json wfmash.multiplicity.json fastga.stderr wfmash.stderr job.json \
    comparison.stdout.json "$partial/"
cp -a comparison "$partial/"
mv -- "$partial" "$DURABLE"
printf '%s\n' "$DURABLE"
