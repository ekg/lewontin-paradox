#!/usr/bin/env bash
# Boundary-safe, resumable scale-out lacing for an already queried IMPG stage.
set -euo pipefail

: "${VGP_STAGE_PARTIAL:?VGP_STAGE_PARTIAL is required}"
: "${VGP_SELECTION_ID:?VGP_SELECTION_ID is required}"
: "${VGP_DATA_ROOT:?VGP_DATA_ROOT is required}"
: "${SLURM_TMPDIR:?SLURM_TMPDIR is required}"
: "${h1:?h1 is required}"
: "${impg:?impg is required}"
: "${bcftools:?bcftools is required}"

lace_threads=${VGP_IMPG_LACE_THREADS:-2}
max_inputs_per_chunk=${VGP_IMPG_LACE_MAX_INPUTS:-4096}
available_cpus=${SLURM_CPUS_PER_TASK:-1}
(( lace_threads >= 2 && available_cpus >= lace_threads && max_inputs_per_chunk >= 1 )) || \
    fail "hierarchical IMPG lace has an invalid CPU/thread allocation"
worker_count=$((available_cpus / lace_threads))
work="$SLURM_TMPDIR/impg-hierarchical-${VGP_SELECTION_ID}-${SLURM_JOB_ID:-local}"
mkdir -p "$work/chunk-lists" "$work/lace-chunks" "$work/chunk-temp" \
    "$work/chunk-logs" "$work/h2-samples" "$work/lace-sites" "$work/final-temp"

python3 - "$VGP_STAGE_PARTIAL/vcf.list" "$work/chunk-lists" "$worker_count" \
    "$max_inputs_per_chunk" "$VGP_DATA_ROOT" "$VGP_SELECTION_ID" <<'PY'
import json,math,sys
from pathlib import Path
source,output=Path(sys.argv[1]),Path(sys.argv[2])
workers,maximum=int(sys.argv[3]),int(sys.argv[4])
canonical_root,selection=sys.argv[5],sys.argv[6]
rows=[line.strip() for line in source.read_text().splitlines() if line.strip()]
if not rows: raise SystemExit("IMPG VCF list is empty")
count=max(workers,math.ceil(len(rows)/maximum)); count=min(count,len(rows))
if count > 1000: raise SystemExit("hierarchical lace requires more than 1000 bounded chunks")
width=math.ceil(len(rows)/count); chunks=[]
for number,start in enumerate(range(0,len(rows),width)):
    values=rows[start:start+width]
    path=output/f"{number:03d}.list"
    path.write_text("".join(f"{value}\n" for value in values))
    chunks.append({"chunk":number,"first_input_index":start,
                   "last_input_index":start+len(values)-1,"input_count":len(values)})
if sum(row["input_count"] for row in chunks) != len(rows):
    raise SystemExit("hierarchical split is not exhaustive")
(output/"split_manifest.json").write_text(json.dumps({
    "schema_version":"vgp-impg-hierarchical-lace-split-v1",
    "canonical_vgp_root":canonical_root,"selection_id":selection,
    "input_count":len(rows),"chunk_count":len(chunks),"chunks":chunks,
    "parallel_worker_count":workers,"maximum_inputs_per_chunk":maximum,
    "split_is_disjoint_and_exhaustive":True,
},sort_keys=True)+"\n")
PY

mapfile -t chunk_lists < <(printf '%s\n' "$work/chunk-lists"/[0-9][0-9][0-9].list)
lace_failed=0
for ((batch_start=0; batch_start<${#chunk_lists[@]}; batch_start+=worker_count)); do
    declare -a lace_pids=()
    declare -a lace_ids=()
    for ((index=batch_start; index<${#chunk_lists[@]} && index<batch_start+worker_count; index++)); do
        chunk_list=${chunk_lists[$index]}
        chunk_id=${chunk_list##*/}; chunk_id=${chunk_id%.list}
        mkdir -p "$work/chunk-temp/$chunk_id"
        "$impg" lace -l "$chunk_list" --format vcf \
            -o "$work/lace-chunks/$chunk_id.vcf.bz2" --reference "$h1" \
            --temp-dir "$work/chunk-temp/$chunk_id" --compress bgzip -t "$lace_threads" \
            >"$work/chunk-logs/$chunk_id.log" 2>&1 &
        lace_pids+=("$!"); lace_ids+=("$chunk_id")
    done
    for index in "${!lace_pids[@]}"; do
        if wait "${lace_pids[$index]}"; then
            [[ -s "$work/lace-chunks/${lace_ids[$index]}.vcf.bz2" ]] || lace_failed=1
        else
            lace_failed=1
        fi
    done
    (( lace_failed == 0 )) || break
done
if (( lace_failed != 0 )); then
    for log in "$work/chunk-logs"/*.log; do tail -n 40 -- "$log" >&2 || true; done
    fail "one or more hierarchical IMPG lace chunks failed"
fi

python3 - "$h1" "$work/lace-chunks" "$work/h2-samples" \
    "$VGP_STAGE_PARTIAL/h2_sample_projection.json" "$VGP_DATA_ROOT" "$VGP_SELECTION_ID" <<'PY'
import bz2,json,sys
from pathlib import Path
h1,chunk_root,output_root,manifest=map(Path,sys.argv[1:5])
canonical_root,selection=sys.argv[5:7]
h1_contigs=set()
with h1.open() as handle:
    for line in handle:
        if line.startswith(">"): h1_contigs.add(line[1:].split()[0])
if not h1_contigs: raise SystemExit("H1 reference has no contigs")
rows=[]; all_h2=[]
for chunk in sorted(chunk_root.glob("[0-9][0-9][0-9].vcf.bz2")):
    samples=None
    with bz2.open(chunk,"rt") as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                samples=line.rstrip("\n").split("\t")[9:]; break
    if samples is None: raise SystemExit(f"IMPG chunk lacks sample header: {chunk}")
    h2_samples=[sample for sample in samples if sample.split(":",1)[0] not in h1_contigs]
    if not h2_samples: raise SystemExit(f"IMPG chunk has no H2 region sample: {chunk}")
    chunk_id=chunk.name.removesuffix(".vcf.bz2")
    (output_root/f"{chunk_id}.txt").write_text("".join(f"{sample}\n" for sample in h2_samples))
    all_h2.extend(h2_samples)
    rows.append({"chunk":chunk_id,"total_sample_count":len(samples),
                 "h2_sample_count":len(h2_samples),"h1_sample_count":len(samples)-len(h2_samples)})
manifest.write_text(json.dumps({
    "schema_version":"vgp-impg-h2-sample-projection-v1",
    "canonical_vgp_root":canonical_root,"selection_id":selection,
    "h1_contig_count":len(h1_contigs),"chunk_count":len(rows),
    "h2_sample_count":len(all_h2),"unique_h2_sample_count":len(set(all_h2)),
    "chunks":rows,"projection_rule":"sample contig prefix absent from H1 dictionary",
    "chunk_compression":"bzip2 (observed BZ stream from pinned IMPG --compress bgzip)",
},sort_keys=True)+"\n")
PY

for chunk_vcf in "$work/lace-chunks"/[0-9][0-9][0-9].vcf.bz2; do
    chunk_name=${chunk_vcf##*/}; chunk_stem=${chunk_name%.vcf.bz2}
    site_vcf="$work/lace-sites/$chunk_stem.vcf.gz"
    python3 -c \
        'import bz2,shutil,sys; source=bz2.open(sys.argv[1],"rb"); shutil.copyfileobj(source,sys.stdout.buffer,1024*1024); source.close()' \
        "$chunk_vcf" |
        "$bcftools" view --samples-file "$work/h2-samples/$chunk_stem.txt" \
        --trim-alt-alleles --min-ac 1:nref -Ou - | \
        "$bcftools" view --drop-genotypes --no-update -Oz -o "$site_vcf"
    "$bcftools" index -f -t "$site_vcf"
done
python3 - "$work/lace-sites" "$work/lace-sites.list" <<'PY'
import sys
from pathlib import Path
root,output=map(Path,sys.argv[1:])
rows=sorted(root.glob("[0-9][0-9][0-9].vcf.gz"))
if not rows or any(not path.stat().st_size for path in rows):
    raise SystemExit("site-only lace projection is absent or empty")
output.write_text("".join(f"{path}\n" for path in rows))
PY

"$bcftools" concat -a -d exact -f "$work/lace-sites.list" -Ov \
    -o "$work/laced.concatenated.vcf"
"$bcftools" sort -m 8G --temp-dir "$work/final-temp" -Ov \
    -o "$VGP_STAGE_PARTIAL/laced.vcf" "$work/laced.concatenated.vcf"
[[ -s $VGP_STAGE_PARTIAL/laced.vcf ]] || fail "hierarchical IMPG lace emitted no VCF"

python3 - "$work/chunk-lists/split_manifest.json" "$work/lace-chunks" \
    "$work/lace-sites" "$VGP_STAGE_PARTIAL/laced.vcf" \
    "$VGP_STAGE_PARTIAL/hierarchical_lace_audit.json" "$VGP_DATA_ROOT" \
    "$VGP_SELECTION_ID" "$lace_threads" "$worker_count" "$max_inputs_per_chunk" <<'PY'
import hashlib,json,sys
from pathlib import Path
split_path,chunk_root,site_root,final_path,output=map(Path,sys.argv[1:6])
canonical_root,selection=sys.argv[6],sys.argv[7]
threads,workers,maximum=map(int,sys.argv[8:11])
split=json.loads(split_path.read_text())
def record(path):
    digest=hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda:handle.read(1024*1024),b""): digest.update(block)
    return {"name":path.name,"bytes":path.stat().st_size,"sha256":digest.hexdigest()}
chunks=[record(path) for path in sorted(chunk_root.glob("[0-9][0-9][0-9].vcf.bz2"))]
sites=[record(path) for path in sorted(site_root.glob("[0-9][0-9][0-9].vcf.gz"))]
if len(chunks) != split["chunk_count"] or len(sites) != split["chunk_count"]:
    raise SystemExit("hierarchical output census differs from split manifest")
output.write_text(json.dumps({
    "schema_version":"vgp-impg-hierarchical-lace-audit-v1",
    "canonical_vgp_root":canonical_root,"selection_id":selection,
    "source_regional_vcf_count":split["input_count"],"chunk_count":len(chunks),
    "worker_threads":threads,"parallel_worker_count":workers,
    "maximum_inputs_per_chunk":maximum,"split_is_disjoint_and_exhaustive":True,
    "chunk_compression":"bzip2 (observed BZ stream from pinned IMPG --compress bgzip)",
    "site_only_projection":True,"h2_nonreference_projection":True,
    "boundary_reconciliation_engine":"pinned bcftools concat -a -d exact then coordinate sort",
    "final_boundary_reconciliation":True,
    "regional_shards_removed_after_verified_lacing":True,
    "chunk_outputs":chunks,"site_projection_outputs":sites,"final_laced_vcf":record(final_path),
},sort_keys=True)+"\n")
PY

rm -rf -- "$VGP_STAGE_PARTIAL/calls" "$work"
rm -- "$VGP_STAGE_PARTIAL/vcf.list"
