#!/usr/bin/env bash
#SBATCH --job-name=vgp10-P07-impg-rescue
#SBATCH --partition=highmem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=1-00:00:00
#SBATCH --output=/moosefs/erikg/vgp/pilot/logs/vgp10-auth-v2-P07-impg-rescue-%j.out
#SBATCH --error=/moosefs/erikg/vgp/pilot/logs/vgp10-auth-v2-P07-impg-rescue-%j.err
set -euo pipefail

SUBMIT_ROOT=${SLURM_SUBMIT_DIR:?submit from the repository root}
CANONICAL_ROOT=/moosefs/erikg/vgp
AUTHORIZATION_ID=vgp10-auth-20260718-v2
SELECTION_ID=P07
RESCUE="$CANONICAL_ROOT/pilot/authorized-checkpoints/$AUTHORIZATION_ID/$SELECTION_ID/rescue/impg-partial-1781666.tar"
RESCUE_SHA256=1db1573e24991d3cbebfe9e28e3ca99421632e87ef73c327e9a2909ab65563d8
FINAL="$CANONICAL_ROOT/pilot/authorized-checkpoints/$AUTHORIZATION_ID/$SELECTION_ID/run/impg"
CAPTURE="$SUBMIT_ROOT/analysis/guix/vgp_10_pilot/realization.json"
H1_CAS="$CANONICAL_ROOT/objects/sha256/2b/75/2b758e606304f7cb5e795d7939979b08c21bf4f3eac7ea3cf1c6ab0a463733c7"

[[ -f $RESCUE && -f $CAPTURE && -f $H1_CAS ]] || {
    echo "rescue, environment capture, or H1 CAS object absent" >&2; exit 2;
}
[[ ! -e $FINAL ]] || { echo "refusing existing IMPG final stage: $FINAL" >&2; exit 2; }

SCRATCH_BASE=/scratch
[[ -d $SCRATCH_BASE && -w $SCRATCH_BASE ]] || {
    echo "required node-local scratch base is absent or unwritable: $SCRATCH_BASE" >&2; exit 2;
}
available_scratch_kib=$(df -Pk -- "$SCRATCH_BASE" | awk 'NR==2 {print $4}')
minimum_scratch_kib=$((20 * 1024 * 1024))
[[ $available_scratch_kib =~ ^[0-9]+$ ]] || {
    echo "could not measure node-local scratch capacity" >&2; exit 2;
}
(( available_scratch_kib >= minimum_scratch_kib )) || {
    echo "insufficient node-local scratch: ${available_scratch_kib} KiB available" >&2; exit 2;
}
scratch=$(mktemp -d -- "$SCRATCH_BASE/vgp-P07-impg-rescue-${SLURM_JOB_ID:?}-XXXXXX")
filesystem_type=$(stat -f -c %T -- "$scratch")
case "$filesystem_type" in
    nfs|nfs4|fuse*|lustre|gpfs|ceph)
        echo "rescue scratch is not node-local: $filesystem_type" >&2; exit 2 ;;
esac
cleanup() {
    [[ $scratch == "$SCRATCH_BASE/vgp-P07-impg-rescue-${SLURM_JOB_ID}-"* ]] || return 2
    rm -rf -- "$scratch"
}
trap cleanup EXIT
export SLURM_TMPDIR=$scratch TMPDIR=$scratch
export VGP_DATA_ROOT="$scratch/data" VGP_RUN_ID=vgp10-auth-20260718-v2-canary
export VGP_ENVIRONMENT_CAPTURE=$CAPTURE PYTHONPATH=$SUBMIT_ROOT
mkdir -p "$VGP_DATA_ROOT"
source "$SUBMIT_ROOT/analysis/slurm/vgp_10_pilot/common.sh"
require_runtime
verify_tool impg
verify_tool bcftools
impg=$(tool_path impg)
bcftools=$(tool_path bcftools)

printf '%s  %s\n' "$RESCUE_SHA256" "$RESCUE" | sha256sum -c -
python3 - "$RESCUE" "$scratch" <<'PY'
import sys,tarfile
from pathlib import Path
archive,destination=Path(sys.argv[1]),Path(sys.argv[2]).resolve()
with tarfile.open(archive,"r:") as handle:
    members=handle.getmembers()
    for member in members:
        resolved=(destination/member.name).resolve()
        if not resolved.is_relative_to(destination):
            raise SystemExit(f"rescue tar member escapes scratch: {member.name}")
        if member.issym() or member.islnk():
            raise SystemExit(f"rescue tar contains a link: {member.name}")
    handle.extractall(destination,members=members)
PY
partial="$scratch/vgp-vgp10-auth-20260718-v2-canary-P07-impg-1781666.partial"
[[ -d $partial/calls && -s $partial/vcf.list && -s $partial/regional_vcf_audit.json ]] || {
    echo "rescued IMPG shard set is incomplete" >&2; exit 2;
}
python3 - "$partial/vcf.list" "$partial/calls" "$partial/rescue_path_rebinding.json" <<'PY'
import json,sys
from pathlib import Path

vcf_list,calls,output=map(Path,sys.argv[1:])
original=[Path(line.strip()) for line in vcf_list.read_text().splitlines() if line.strip()]
if not original:
    raise SystemExit("rescued IMPG VCF list is empty")
names=[path.name for path in original]
if len(names) != len(set(names)):
    raise SystemExit("cannot safely rebind rescued IMPG VCF list: duplicate basenames")
rebound=[calls/path.name for path in original]
missing=[str(path) for path in rebound if not path.is_file() or path.stat().st_size == 0]
if missing:
    raise SystemExit(f"cannot safely rebind rescued IMPG VCF list: {len(missing)} shards absent or empty")
temporary=vcf_list.with_suffix(".list.rebound")
temporary.write_text("".join(f"{path}\n" for path in rebound))
temporary.replace(vcf_list)
output.write_text(json.dumps({
  "schema_version":"vgp-impg-rescue-path-rebinding-v1",
  "task_id":"run-vgp-real-canary",
  "canonical_vgp_root":"/moosefs/erikg/vgp",
  "source_slurm_job_id":"1781666",
  "input_count":len(original),
  "unique_basename_count":len(set(names)),
  "original_parent_count":len({str(path.parent) for path in original}),
  "rebound_calls_directory":str(calls),
  "rescued_vcf_list_rebound":True,
  "all_rebound_shards_nonempty":True,
},sort_keys=True)+"\n")
PY
rm -f -- "$partial/laced.vcf"

cp "$H1_CAS" "$scratch/h1.fa.gz"
printf '%s  %s\n' \
    2b758e606304f7cb5e795d7939979b08c21bf4f3eac7ea3cf1c6ab0a463733c7 \
    "$scratch/h1.fa.gz" | sha256sum -c -
python3 - "$scratch/h1.fa.gz" "$scratch/h1.fa" <<'PY'
import gzip,shutil,sys
with gzip.open(sys.argv[1],"rb") as source, open(sys.argv[2],"wb") as target:
    shutil.copyfileobj(source,target)
PY

worker_threads=2
available_cpus=${SLURM_CPUS_PER_TASK:-1}
(( available_cpus >= worker_threads )) || {
    echo "IMPG hierarchical lace requires at least $worker_threads CPUs" >&2; exit 2;
}
chunk_count=$(( available_cpus / worker_threads ))
(( chunk_count > 16 )) && chunk_count=16
(( chunk_count > 0 )) || { echo "invalid lace chunk count: $chunk_count" >&2; exit 2; }
mkdir -p "$scratch/chunk-lists" "$scratch/lace-chunks" "$scratch/chunk-temp" "$scratch/chunk-logs"
python3 - "$partial/vcf.list" "$scratch/chunk-lists" "$chunk_count" <<'PY'
import json,math,sys
from pathlib import Path

source,output,requested=Path(sys.argv[1]),Path(sys.argv[2]),int(sys.argv[3])
rows=[line.strip() for line in source.read_text().splitlines() if line.strip()]
if not rows:
    raise SystemExit("rescued IMPG VCF list is empty")
count=min(requested,len(rows))
width=math.ceil(len(rows)/count)
chunks=[]
for number,start in enumerate(range(0,len(rows),width)):
    values=rows[start:start+width]
    path=output/f"{number:03d}.list"
    path.write_text("".join(f"{value}\n" for value in values))
    chunks.append({"chunk":number,"first_input_index":start,
                   "last_input_index":start+len(values)-1,"input_count":len(values)})
if sum(row["input_count"] for row in chunks) != len(rows):
    raise SystemExit("hierarchical lace split did not preserve every input exactly once")
(output/"split_manifest.json").write_text(json.dumps({
    "schema_version":"vgp-impg-hierarchical-lace-split-v1",
    "canonical_vgp_root":"/moosefs/erikg/vgp",
    "input_count":len(rows),"chunk_count":len(chunks),"chunks":chunks,
},sort_keys=True)+"\n")
PY

declare -a lace_pids=()
declare -a lace_ids=()
for chunk_list in "$scratch/chunk-lists"/[0-9][0-9][0-9].list; do
    chunk_id=${chunk_list##*/}
    chunk_id=${chunk_id%.list}
    mkdir -p "$scratch/chunk-temp/$chunk_id"
    "$impg" lace -l "$chunk_list" --format vcf \
        -o "$scratch/lace-chunks/$chunk_id.vcf" --reference "$scratch/h1.fa" \
        --temp-dir "$scratch/chunk-temp/$chunk_id" --compress none -t "$worker_threads" \
        >"$scratch/chunk-logs/$chunk_id.log" 2>&1 &
    lace_pids+=("$!")
    lace_ids+=("$chunk_id")
done
lace_failed=0
for index in "${!lace_pids[@]}"; do
    pid=${lace_pids[$index]}
    if wait "$pid"; then
        [[ -s "$scratch/lace-chunks/${lace_ids[$index]}.vcf" ]] || lace_failed=1
    else
        lace_failed=1
    fi
done
if (( lace_failed != 0 )); then
    for log in "$scratch/chunk-logs"/*.log; do
        echo "===== $log =====" >&2
        tail -n 40 -- "$log" >&2 || true
    done
    echo "one or more hierarchical IMPG lace chunks failed" >&2
    exit 2
fi
python3 - "$scratch/lace-chunks" "$scratch/lace-chunks.list" <<'PY'
import sys
from pathlib import Path
root,output=Path(sys.argv[1]),Path(sys.argv[2])
chunks=sorted(root.glob("[0-9][0-9][0-9].vcf"))
if not chunks or any(path.stat().st_size == 0 for path in chunks):
    raise SystemExit("hierarchical lace chunk output is absent or empty")
output.write_text("".join(f"{path}\n" for path in chunks))
PY
mkdir -p "$scratch/h2-samples"
python3 - "$scratch/h1.fa" "$scratch/lace-chunks" "$scratch/h2-samples" \
    "$partial/h2_sample_projection.json" <<'PY'
import json,sys
from pathlib import Path

h1,chunk_root,output_root,manifest=map(Path,sys.argv[1:])
h1_contigs=set()
with h1.open() as handle:
    for line in handle:
        if line.startswith(">"):
            h1_contigs.add(line[1:].split()[0])
if not h1_contigs:
    raise SystemExit("H1 reference contains no contigs")
rows=[]
all_h2_samples=[]
for chunk in sorted(chunk_root.glob("[0-9][0-9][0-9].vcf")):
    samples=None
    with chunk.open() as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                samples=line.rstrip("\n").split("\t")[9:]
                break
    if samples is None:
        raise SystemExit(f"IMPG chunk lacks a VCF sample header: {chunk}")
    h2_samples=[sample for sample in samples if sample.split(":",1)[0] not in h1_contigs]
    if not h2_samples:
        raise SystemExit(f"IMPG chunk has no H2 region samples: {chunk}")
    sample_path=output_root/f"{chunk.stem}.txt"
    sample_path.write_text("".join(f"{sample}\n" for sample in h2_samples))
    all_h2_samples.extend(h2_samples)
    rows.append({"chunk":chunk.stem,"total_sample_count":len(samples),
                 "h2_sample_count":len(h2_samples),"h1_sample_count":len(samples)-len(h2_samples)})
manifest.write_text(json.dumps({
  "schema_version":"vgp-impg-h2-sample-projection-v1",
  "task_id":"run-vgp-real-canary",
  "canonical_vgp_root":"/moosefs/erikg/vgp",
  "h1_contig_count":len(h1_contigs),
  "chunk_count":len(rows),
  "h2_sample_count":len(all_h2_samples),
  "unique_h2_sample_count":len(set(all_h2_samples)),
  "chunks":rows,
  "projection_rule":"sample contig prefix absent from the H1 reference sequence dictionary",
},sort_keys=True)+"\n")
PY
mkdir -p "$scratch/lace-sites"
for chunk_vcf in "$scratch/lace-chunks"/[0-9][0-9][0-9].vcf; do
    chunk_id=${chunk_vcf##*/}
    chunk_stem=${chunk_id%.vcf}
    site_vcf="$scratch/lace-sites/$chunk_id.gz"
    "$bcftools" view --samples-file "$scratch/h2-samples/$chunk_stem.txt" \
        --trim-alt-alleles --min-ac 1:nref -Ou "$chunk_vcf" | \
        "$bcftools" view --drop-genotypes --no-update -Oz -o "$site_vcf"
    "$bcftools" index -f -t "$site_vcf"
    [[ -s $site_vcf && -s $site_vcf.tbi ]] || {
        echo "site-only projection is absent for $chunk_id" >&2; exit 2;
    }
done
python3 - "$scratch/lace-sites" "$scratch/lace-sites.list" <<'PY'
import sys
from pathlib import Path
root,output=Path(sys.argv[1]),Path(sys.argv[2])
sites=sorted(root.glob("[0-9][0-9][0-9].vcf.gz"))
if not sites or any(path.stat().st_size == 0 for path in sites):
    raise SystemExit("site-only IMPG lace projection is absent or empty")
output.write_text("".join(f"{path}\n" for path in sites))
PY
mkdir -p "$scratch/final-temp"
"$bcftools" concat -a -d exact -f "$scratch/lace-sites.list" -Ov \
    -o "$scratch/laced.concatenated.vcf"
"$bcftools" sort -m 8G --temp-dir "$scratch/final-temp" -Ov \
    -o "$partial/laced.vcf" "$scratch/laced.concatenated.vcf"
rm -- "$scratch/laced.concatenated.vcf"
[[ -s $partial/laced.vcf ]] || { echo "rescued IMPG lace emitted no VCF" >&2; exit 2; }
python3 - "$scratch/chunk-lists/split_manifest.json" "$scratch/lace-chunks" \
    "$scratch/lace-sites" "$partial/laced.vcf" "$partial/hierarchical_lace_audit.json" <<'PY'
import hashlib,json,sys
from pathlib import Path

split_path,chunk_root,site_root,final_path,output=map(Path,sys.argv[1:])
split=json.loads(split_path.read_text())
def digest(path):
    value=hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda:handle.read(1024*1024),b""):
            value.update(block)
    return value.hexdigest()
chunks=[]
for path in sorted(chunk_root.glob("[0-9][0-9][0-9].vcf")):
    chunks.append({"name":path.name,"bytes":path.stat().st_size,"sha256":digest(path)})
if len(chunks) != split["chunk_count"]:
    raise SystemExit("hierarchical lace output count differs from split count")
sites=[]
for path in sorted(site_root.glob("[0-9][0-9][0-9].vcf.gz")):
    index=Path(f"{path}.tbi")
    if not index.is_file() or index.stat().st_size == 0:
        raise SystemExit(f"site-only projection index absent: {index}")
    sites.append({"name":path.name,"bytes":path.stat().st_size,"sha256":digest(path),
                  "index_bytes":index.stat().st_size,"index_sha256":digest(index)})
if len(sites) != split["chunk_count"]:
    raise SystemExit("site-only projection count differs from split count")
value={
  "schema_version":"vgp-impg-hierarchical-lace-audit-v1",
  "task_id":"run-vgp-real-canary",
  "canonical_vgp_root":"/moosefs/erikg/vgp",
  "source_regional_vcf_count":split["input_count"],
  "chunk_count":len(chunks),
  "split_is_disjoint_and_exhaustive":True,
  "site_only_projection":True,
  "h2_nonreference_projection":True,
  "site_projection_rationale":"retain ALT alleles carried by H2 region samples, then remove region-window bookkeeping genotypes",
  "boundary_reconciliation_engine":"pinned bcftools concat -a -d exact followed by coordinate sort",
  "final_boundary_reconciliation":True,
  "chunk_outputs":chunks,
  "site_projection_outputs":sites,
  "final_laced_vcf":{"bytes":final_path.stat().st_size,"sha256":digest(final_path)},
}
output.write_text(json.dumps(value,sort_keys=True)+"\n")
PY
rm -rf -- "$partial/calls"
rm -- "$partial/vcf.list"
python3 - "$partial/rescue_provenance.json" "$filesystem_type" "$chunk_count" "$worker_threads" \
    "$SCRATCH_BASE" "$available_scratch_kib" <<'PY'
import json,os,sys
from pathlib import Path
value={
  "schema_version":"vgp-impg-rescue-v1",
  "task_id":"run-vgp-real-canary",
  "authorization_id":"vgp10-auth-20260718-v2",
  "selection_id":"P07",
  "canonical_vgp_root":"/moosefs/erikg/vgp",
  "source_slurm_job_id":"1781666",
  "rescue_slurm_job_id":os.environ["SLURM_JOB_ID"],
  "rescue_tar_sha256":"1db1573e24991d3cbebfe9e28e3ca99421632e87ef73c327e9a2909ab65563d8",
  "node_local_filesystem_type":sys.argv[2],
  "node_local_scratch_base":sys.argv[5],
  "node_local_scratch_available_kib_at_start":int(sys.argv[6]),
  "rescued_vcf_list_rebound":True,
  "lace_temp_explicitly_node_local":True,
  "hierarchical_lace_chunk_count":int(sys.argv[3]),
  "hierarchical_lace_worker_threads":int(sys.argv[4]),
  "site_only_projection":True,
  "h2_nonreference_projection":True,
  "boundary_reconciliation_engine":"pinned bcftools concat -a -d exact followed by coordinate sort",
  "final_boundary_reconciliation":True,
  "regional_shards_removed_after_verified_lacing":True,
}
Path(sys.argv[1]).write_text(json.dumps(value,sort_keys=True)+"\n")
PY

durable_partial="${FINAL}.rescue-${SLURM_JOB_ID}.partial"
[[ ! -e $durable_partial ]] || { echo "durable rescue partial already exists" >&2; exit 2; }
cp -a "$partial" "$durable_partial"
python3 - "$durable_partial" "$FINAL" "$SLURM_JOB_ID" <<'PY'
import sys
from pathlib import Path
from analysis.vgp_10_pilot import atomic_promote
atomic_promote(Path(sys.argv[1]),Path(sys.argv[2]),{
  "selection_id":"P07",
  "stage":"impg",
  "atomic_promotion":True,
  "rescued_from_completed_regional_shards":True,
  "source_slurm_job_id":"1781666",
  "rescue_slurm_job_id":sys.argv[3],
  "canonical_vgp_root":"/moosefs/erikg/vgp",
})
PY
