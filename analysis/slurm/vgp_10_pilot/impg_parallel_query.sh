#!/usr/bin/env bash
# Disjoint/exhaustive concurrent IMPG query extraction within one allocation.
set -euo pipefail

: "${VGP_STAGE_PARTIAL:?VGP_STAGE_PARTIAL is required}"
: "${VGP_SELECTION_ID:?VGP_SELECTION_ID is required}"
: "${VGP_DATA_ROOT:?VGP_DATA_ROOT is required}"
: "${SLURM_TMPDIR:?SLURM_TMPDIR is required}"
: "${h1:?h1 is required}"
: "${h2:?h2 is required}"
: "${paf:?paf is required}"
: "${impg:?impg is required}"

available_cpus=${SLURM_CPUS_PER_TASK:-1}
query_threads=2
if (( available_cpus < query_threads )); then query_threads=$available_cpus; fi
query_workers=$((available_cpus / query_threads))
(( query_workers > 16 )) && query_workers=16
(( query_workers > 0 )) || fail "IMPG query requires at least one CPU"
query_root="$SLURM_TMPDIR/impg-query-$VGP_SELECTION_ID-${SLURM_JOB_ID:-local}"
mkdir -p "$query_root/beds" "$query_root/temp" "$query_root/logs" \
    "$VGP_STAGE_PARTIAL/calls"

python3 - "$VGP_STAGE_PARTIAL/focus.native.bed" "$query_root/beds" \
    "$query_workers" "$VGP_DATA_ROOT" "$VGP_SELECTION_ID" <<'PY'
import json,math,sys
from pathlib import Path
source,out=Path(sys.argv[1]),Path(sys.argv[2]); requested=int(sys.argv[3])
rows=[line for line in source.read_text().splitlines(keepends=True) if line.strip()]
if not rows: raise SystemExit("focused IMPG BED is empty")
count=min(requested,len(rows)); width=math.ceil(len(rows)/count); chunks=[]
for number,start in enumerate(range(0,len(rows),width)):
    values=rows[start:start+width]; path=out/f"{number:03d}.bed"
    path.write_text("".join(values))
    chunks.append({"chunk":number,"first_query_index":start,
                   "last_query_index":start+len(values)-1,"query_count":len(values)})
if sum(row["query_count"] for row in chunks) != len(rows):
    raise SystemExit("parallel IMPG query split is not exhaustive")
(out/"split_manifest.json").write_text(json.dumps({
 "schema_version":"vgp-impg-parallel-query-split-v1",
 "canonical_vgp_root":sys.argv[4],"selection_id":sys.argv[5],
 "query_count":len(rows),"chunk_count":len(chunks),"chunks":chunks,
 "split_is_disjoint_and_exhaustive":True,
},sort_keys=True)+"\n")
PY

declare -a query_pids=()
declare -a query_ids=()
for query_bed in "$query_root/beds"/[0-9][0-9][0-9].bed; do
    query_id=${query_bed##*/}; query_id=${query_id%.bed}
    mkdir -p "$query_root/temp/$query_id" "$VGP_STAGE_PARTIAL/calls/$query_id"
    "$impg" query -a "$paf" -i "$VGP_STAGE_PARTIAL/h1_h2.impg" \
        -b "$query_bed" -d 0 --min-transitive-len 1 \
        --force-large-region \
        --temp-dir "$query_root/temp/$query_id" -o vcf:poa \
        --sequence-files "$h1" "$h2" -O "$VGP_STAGE_PARTIAL/calls/$query_id" \
        -t "$query_threads" >"$query_root/logs/$query_id.log" 2>&1 &
    query_pids+=("$!"); query_ids+=("$query_id")
done
query_failed=0
for index in "${!query_pids[@]}"; do
    if ! wait "${query_pids[$index]}"; then query_failed=1; fi
done
if (( query_failed != 0 )); then
    for log in "$query_root/logs"/*.log; do tail -n 60 -- "$log" >&2 || true; done
    fail "one or more disjoint IMPG query chunks failed"
fi
cp "$query_root/beds/split_manifest.json" "$VGP_STAGE_PARTIAL/parallel_query_audit.json"
[[ $query_root == "$SLURM_TMPDIR/impg-query-$VGP_SELECTION_ID-"* ]] || \
    fail "refusing unsafe IMPG query scratch cleanup: $query_root"
rm -rf -- "$query_root"
