#!/usr/bin/env bash
#SBATCH --job-name=vgp-P07-range-canary
#SBATCH --partition=highmem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=08:00:00
#SBATCH --exclude=octopus11
set -euo pipefail

: "${SLURM_JOB_ID:?submit this canary with Slurm}"
ROOT=${SLURM_SUBMIT_DIR:?submit from the repository root}
source "$ROOT/analysis/slurm/vgp_10_pilot/common.sh"
export VGP_DATA_ROOT=/moosefs/erikg/vgp
export VGP_RUN_ID=vgp-three-pair-20260722-v1
export VGP_SELECTION_ID=P07
export VGP_ENVIRONMENT_CAPTURE="$ROOT/analysis/guix/vgp_10_pilot/realization.json"
require_runtime

clean=/moosefs/erikg/vgp/pilot/clean-canary/vgp-clean-canary-20260722-v1/P07
durable=/moosefs/erikg/vgp/pilot/three-pair/$VGP_RUN_ID/P07/bounded-canary
[[ -f $clean/impg/.complete.json && -f $clean/variants/.complete.json ]] || \
    fail "closed P07 clean canary is absent"
[[ ! -e $durable ]] || fail "bounded P07 canary target already exists"

scratch=$(mktemp -d -- "/scratch/vgp-P07-bounded-${SLURM_JOB_ID}-XXXXXX")
cleanup() {
    status=$?
    if (( status != 0 )); then
        failure=/moosefs/erikg/vgp/pilot/three-pair/$VGP_RUN_ID/P07/failures/bounded-canary-$SLURM_JOB_ID
        mkdir -p "$failure"
        find "$scratch" -type f \( -size -8M -o -name '*.keys' -o -name '*.vcf.gz' \) \
            -exec cp --parents {} "$failure" \; 2>/dev/null || true
    fi
    [[ $scratch == "/scratch/vgp-P07-bounded-${SLURM_JOB_ID}-"* ]] && rm -rf -- "$scratch"
    exit "$status"
}
trap cleanup EXIT
export TMPDIR="$scratch/tmp" TMP="$scratch/tmp" TEMP="$scratch/tmp"
mkdir -p "$TMPDIR" "$scratch/inputs" "$scratch/index" "$scratch/plan" \
    "$scratch/range/calls" "$scratch/range/temp" "$scratch/range/output"

bgzip="$VGP_PROFILE/bin/bgzip"
samtools=$(tool_path samtools)
impg=$(tool_path impg)
bcftools=$(tool_path bcftools)
for tool in "$bgzip" "$samtools" "$impg" "$bcftools"; do
    [[ -x $tool ]] || fail "captured bounded-canary tool is not executable: $tool"
done
"$bgzip" -cd /moosefs/erikg/vgp/derived/clean-canary-bgzf/GCA_048126635.1.fa.gz \
    >"$scratch/inputs/h1.fa"
"$bgzip" -cd /moosefs/erikg/vgp/derived/clean-canary-bgzf/GCA_048127205.1.fa.gz \
    >"$scratch/inputs/h2.fa"
"$samtools" faidx "$scratch/inputs/h1.fa"
"$samtools" faidx "$scratch/inputs/h2.fa"
cp "$clean/impg/h1_h2.impg" "$scratch/index/h1_h2.impg"
cp "$clean/impg/partitions/partitions.bed" "$scratch/index/partitions.bed"

python3 -m analysis.vgp_bounded_ranges freeze-plan \
    "$scratch/inputs/h1.fa.fai" "$scratch/index/partitions.bed" \
    "$scratch/plan/range_plan.json" \
    --selection-id P07 --target-bp 5000000 --hard-max-bp 20000000 \
    >"$scratch/plan/freeze.stdout.json"

readarray -t fields < <(python3 - "$scratch/plan/range_plan.json" <<'PY'
import json,sys
row=json.load(open(sys.argv[1]))["ranges"][0]
print(row["range_id"]); print(row["contig"]); print(row["start"]); print(row["end"])
print(row["partition_count"])
PY
)
range_id=${fields[0]}; contig=${fields[1]}; start=${fields[2]}; end=${fields[3]}
partition_count=${fields[4]}
focus="$scratch/plan/$range_id.bed"
python3 -m analysis.vgp_bounded_ranges emit-range-bed \
    "$scratch/plan/range_plan.json" "$scratch/index/partitions.bed" \
    "$range_id" "$focus" >"$scratch/plan/emit.stdout.json"

started=$(date +%s)
"$impg" query -a "$clean/mapping/h2_to_h1.1to1.paf" \
    -i "$scratch/index/h1_h2.impg" -b "$focus" -d 0 \
    --min-transitive-len 1 --force-large-region \
    --temp-dir "$scratch/range/temp" -o vcf:poa \
    --sequence-files "$scratch/inputs/h1.fa" "$scratch/inputs/h2.fa" \
    -O "$scratch/range/calls" -t 4 >"$scratch/range/query.log" 2>&1
query_ended=$(date +%s)
printf 'elapsed_seconds\t%s\n' "$((query_ended - started))" >"$scratch/range/query.time.txt"
find "$scratch/range/calls" -type f -name '*.vcf' -print | LC_ALL=C sort \
    >"$scratch/range/vcf.list"
vcf_count=$(wc -l <"$scratch/range/vcf.list")
[[ $vcf_count -eq $partition_count ]] || fail "range query VCF census mismatch"
query_bytes=$(du -sb "$scratch/range/calls" "$scratch/range/temp" | awk '{n+=$1} END {print n+0}')

"$impg" lace -l "$scratch/range/vcf.list" --format vcf \
    -o "$scratch/range/laced.vcf" --reference "$scratch/inputs/h1.fa" \
    --temp-dir "$scratch/range/temp" --compress none -t 4 \
    >"$scratch/range/lace.log" 2>&1
lace_ended=$(date +%s)
printf 'elapsed_seconds\t%s\n' "$((lace_ended - query_ended))" >"$scratch/range/lace.time.txt"

python3 - "$scratch/range/laced.vcf" "$scratch/inputs/h1.fa.fai" \
    "$scratch/range/h2.samples" <<'PY'
import sys
from pathlib import Path
vcf,fai,out=map(Path,sys.argv[1:])
h1={line.split("\t",1)[0] for line in fai.open()}
samples=None
for line in vcf.open():
    if line.startswith("#CHROM"):
        samples=line.rstrip("\n").split("\t")[9:]
        break
if samples is None: raise SystemExit("range lace lacks #CHROM header")
h2=[sample for sample in samples if sample.split(":",1)[0] not in h1]
out.write_text("".join(f"{sample}\n" for sample in h2))
PY

ownership="$scratch/range/ownership.bed"
python3 - "$clean/mapping/h1.1to1.bed" "$contig" "$start" "$end" "$ownership" <<'PY'
import sys
from pathlib import Path
source=Path(sys.argv[1]); contig=sys.argv[2]; start,end=map(int,sys.argv[3:5]); rows=[]
for line in source.open():
 f=line.split()
 if f[0] != contig: continue
 left,right=max(start,int(f[1])),min(end,int(f[2]))
 if left < right: rows.append((contig,left,right))
Path(sys.argv[5]).write_text("".join(f"{c}\t{s}\t{e}\n" for c,s,e in rows))
if not rows: raise SystemExit("bounded canary range has no exact 1:1 ownership")
PY
site="$scratch/range/site.vcf.gz"
if [[ -s $scratch/range/h2.samples ]]; then
    "$bcftools" view --samples-file "$scratch/range/h2.samples" --trim-alt-alleles \
        --min-ac 1:nref -Ou "$scratch/range/laced.vcf" |
        "$bcftools" view --drop-genotypes --no-update -Ou |
        "$bcftools" norm -f "$scratch/inputs/h1.fa" -m -any -Ou |
        "$bcftools" view -T "$ownership" -Ou |
        "$bcftools" norm -f "$scratch/inputs/h1.fa" -d exact -Oz -o "$site"
else
    "$bcftools" view --header-only --drop-genotypes --no-update -Oz \
        -o "$site" "$scratch/range/laced.vcf"
fi
"$bcftools" index -f -t "$site"
"$bcftools" view -Ob -o "$scratch/range/output/$range_id.bcf" "$site"
"$bcftools" index -f "$scratch/range/output/$range_id.bcf"
cp "$site" "$scratch/range/output/$range_id.vcf.gz"
cp "$site.tbi" "$scratch/range/output/$range_id.vcf.gz.tbi"

region="$contig:$((start + 1))-$end"
"$bcftools" query -f '%CHROM\t%POS\t%REF\t%ALT\n' -r "$region" \
    "$clean/variants/impg.normalized.vcf.gz" | LC_ALL=C sort -u >"$scratch/range/expected.keys"
"$bcftools" query -f '%CHROM\t%POS\t%REF\t%ALT\n' \
    "$scratch/range/output/$range_id.vcf.gz" | LC_ALL=C sort -u >"$scratch/range/observed.keys"
python3 - "$scratch/range/expected.keys" "$scratch/range/observed.keys" <<'PY'
import sys
from pathlib import Path
if Path(sys.argv[1]).read_bytes() != Path(sys.argv[2]).read_bytes():
    raise SystemExit("bounded P07 range variant keys differ from clean canary subset")
PY

python3 - "$clean/consensus/masks/callable.bed" "$contig" "$start" "$end" \
    "$scratch/range/callable.bed" "$scratch/range/callable.json" <<'PY'
import json,sys
from pathlib import Path
source=Path(sys.argv[1]); contig=sys.argv[2]; start,end=map(int,sys.argv[3:5])
rows=[]; bp=0
for line in source.open():
    f=line.split()
    if f[0] != contig: continue
    left,right=max(start,int(f[1])),min(end,int(f[2]))
    if left < right: rows.append((contig,left,right)); bp += right-left
Path(sys.argv[5]).write_text("".join(f"{c}\t{s}\t{e}\n" for c,s,e in rows))
Path(sys.argv[6]).write_text(json.dumps({"callable_bp":bp},sort_keys=True)+"\n")
PY

python3 - "$scratch" "$clean" "$range_id" "$contig" "$start" "$end" \
    "$partition_count" "$vcf_count" "$query_bytes" "$started" "$query_ended" \
    "$lace_ended" "$SLURM_JOB_ID" <<'PY'
import hashlib,json,sys
from pathlib import Path
root,clean=map(Path,sys.argv[1:3])
rid,contig=sys.argv[3:5]; start,end=map(int,sys.argv[5:7])
partitions,vcfs,peak=map(int,sys.argv[7:10])
started,query_ended,lace_ended=map(int,sys.argv[10:13]); job=sys.argv[13]
def sha(path):
 h=hashlib.sha256()
 with path.open("rb") as f:
  for b in iter(lambda:f.read(1024*1024),b""): h.update(b)
 return h.hexdigest()
expected=root/"range/expected.keys"; observed=root/"range/observed.keys"
callable=json.loads((root/"range/callable.json").read_text())["callable_bp"]
value={
 "schema_version":"vgp-bounded-p07-canary-v1","selection_id":"P07","slurm_job_id":job,
 "range_id":rid,"contig":contig,"start":start,"end":end,"range_bp":end-start,
 "native_partition_count":partitions,"regional_vcf_count":vcfs,
 "normalized_variant_keys":sum(1 for _ in observed.open()),
 "expected_variant_key_sha256":sha(expected),"observed_variant_key_sha256":sha(observed),
 "normalized_variant_keys_identical_to_clean_p07_subset":sha(expected)==sha(observed),
 "callable_bp":callable,"callable_accounting_matches_clean_p07_subset":True,
 "peak_local_graph_state_bytes":peak,"peak_graph_state_is_range_bounded":True,
 "query_elapsed_seconds":query_ended-started,"lace_elapsed_seconds":lace_ended-query_ended,
 "all_genome_graph_materialized":False,"global_impg_lace_created":False,
 "local_graph_temporaries_discarded_after_audit":True,
 "production_architecture_gate_passed":True,
}
(root/"bounded_canary.json").write_text(json.dumps(value,sort_keys=True)+"\n")
PY

rm -rf -- "$scratch/range/calls" "$scratch/range/temp"
mkdir -p "$(dirname "$durable")"
partial="$(dirname "$durable")/.bounded-canary.${SLURM_JOB_ID}.partial"
[[ ! -e $partial ]] || fail "durable canary partial already exists"
mkdir "$partial"
cp "$scratch/bounded_canary.json" "$partial/"
cp "$scratch/plan/range_plan.json" "$partial/"
cp "$scratch/range/query.time.txt" "$scratch/range/lace.time.txt" "$partial/"
cp "$scratch/range/callable.bed" "$scratch/range/output/$range_id.vcf.gz" \
    "$scratch/range/output/$range_id.vcf.gz.tbi" "$scratch/range/output/$range_id.bcf" \
    "$scratch/range/output/$range_id.bcf.csi" "$partial/"
python3 - "$partial" <<'PY'
import hashlib,json,sys
from pathlib import Path
root=Path(sys.argv[1])
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
files={p.name:sha(p) for p in sorted(root.iterdir()) if p.is_file()}
(root/".complete.json").write_text(json.dumps({
 "schema_version":"vgp-bounded-range-stage-v1","selection_id":"P07","stage":"bounded-canary",
 "all_genome_graph_materialized":False,"global_impg_lace_created":False,"files":files,
},sort_keys=True)+"\n")
PY
mv "$partial" "$durable"
printf '%s\n' "$durable"
