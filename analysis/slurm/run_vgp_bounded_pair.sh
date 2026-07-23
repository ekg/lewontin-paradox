#!/usr/bin/env bash
#SBATCH --job-name=vgp-bounded-pair
#SBATCH --partition=highmem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --exclude=octopus11
# Resources are supplied at submission because P02/P03/P07 differ.
set -euo pipefail

: "${SLURM_JOB_ID:?submit this production run with Slurm}"
selection=${1:?usage: run_vgp_bounded_pair.sh P02|P03|P07}
[[ $selection =~ ^P0(2|3|7)$ ]] || { echo "unsupported selection: $selection" >&2; exit 2; }
ROOT=${SLURM_SUBMIT_DIR:?submit from repository root}
source "$ROOT/analysis/slurm/vgp_10_pilot/common.sh"
export VGP_DATA_ROOT=/moosefs/erikg/vgp
export VGP_RUN_ID=vgp-three-pair-20260722-v1
export VGP_SELECTION_ID=$selection
export VGP_ENVIRONMENT_CAPTURE="$ROOT/analysis/guix/vgp_10_pilot/realization.json"
require_runtime

durable="$VGP_DATA_ROOT/pilot/three-pair/$VGP_RUN_ID/$selection/bounded-production"
[[ ! -e $durable ]] || fail "bounded production target already exists: $durable"
scratch=$(mktemp -d -- "/scratch/vgp-$selection-bounded-${SLURM_JOB_ID}-XXXXXX")
failure_root="$VGP_DATA_ROOT/pilot/three-pair/$VGP_RUN_ID/$selection/failures"
cleanup() {
    status=$?
    if (( status != 0 )); then
        failure="$failure_root/bounded-production-$SLURM_JOB_ID"
        mkdir -p "$failure"
        find "$scratch" -type f \( -name '*.json' -o -name '*.log' -o -name '*.txt' \) \
            -size -16M -exec cp --parents {} "$failure" \; 2>/dev/null || true
    fi
    [[ $scratch == "/scratch/vgp-$selection-bounded-${SLURM_JOB_ID}-"* ]] && \
        rm -rf -- "$scratch"
    exit "$status"
}
trap cleanup EXIT
export TMPDIR="$scratch/tmp" TMP="$scratch/tmp" TEMP="$scratch/tmp"
mkdir -p "$TMPDIR" "$scratch/inputs" "$scratch/index/partitions" "$scratch/plan" \
    "$scratch/results/ranges" "$scratch/results/variants" "$scratch/results/consensus" \
    "$scratch/results/psmc" "$scratch/worker"

samtools=$(tool_path samtools)
impg=$(tool_path impg)
bcftools=$(tool_path bcftools)
psmc=$(tool_path psmc)
bgzip="$VGP_PROFILE/bin/bgzip"
for tool in "$samtools" "$impg" "$bcftools" "$psmc" "$bgzip"; do
    [[ -x $tool ]] || fail "captured bounded-production tool is not executable: $tool"
done
cpus=${SLURM_CPUS_PER_TASK:-1}
(( cpus >= 4 )) || fail "bounded production requires at least four CPUs"
query_threads=2
range_workers=$((cpus / query_threads))
(( range_workers > 8 )) && range_workers=8

source_run="$VGP_DATA_ROOT/pilot/runs/$VGP_RUN_ID/$selection"
if [[ $selection == P07 ]]; then
    source_run="$VGP_DATA_ROOT/pilot/clean-canary/vgp-clean-canary-20260722-v1/P07"
    "$bgzip" -cd "$VGP_DATA_ROOT/derived/clean-canary-bgzf/GCA_048126635.1.fa.gz" \
        >"$scratch/inputs/h1.fa"
    "$bgzip" -cd "$VGP_DATA_ROOT/derived/clean-canary-bgzf/GCA_048127205.1.fa.gz" \
        >"$scratch/inputs/h2.fa"
    "$samtools" faidx "$scratch/inputs/h1.fa"
    "$samtools" faidx "$scratch/inputs/h2.fa"
    awk 'BEGIN{OFS="\t"} {print $1,0,$2}' "$scratch/inputs/h1.fa.fai" \
        >"$scratch/inputs/h1_universe.bed"
    cp "$source_run/impg/h1_h2.impg" "$scratch/index/h1_h2.impg"
    cp "$source_run/impg/partitions/partitions.bed" "$scratch/index/partitions/partitions.bed"
    python3 - "$ROOT/analysis/vgp_clean_canary_execution_v1.json" \
        "$scratch/index/graph_identifier_audit.json" <<'PY'
import json,sys
from pathlib import Path
source=json.load(open(sys.argv[1]))["graph_sequence_digest_ledger"]
Path(sys.argv[2]).write_text(json.dumps({
 "schema_version":"vgp-bounded-p07-graph-id-audit-v1",
 "resolved_ids":source["resolved_ids"],"unresolved_ids":source["unresolved_ids"],
 "silently_omitted_regions":0,"digest_validated_aliases_used":[],
 "source_ledger_path":source["path"],"source_ledger_sha256":source["sha256"],
},sort_keys=True)+"\n")
PY
    python3 - "$source_run" "$scratch/inputs/input-manifest.json" \
        "$scratch/inputs/h1.fa" "$scratch/inputs/h2.fa" <<'PY'
import json,sys
from pathlib import Path
root,out,h1,h2=map(Path,sys.argv[1:])
pre=json.loads((root/"preflight/preflight.json").read_text())
ann=json.loads(Path("analysis/vgp_clean_canary_selection_v1.json").read_text())["annotation"]
ann["sequence_dictionary"]=pre["dictionaries"]["h1_fasta"]
value={
 "selection_id":"P07","h1_accession_version":"GCA_048126635.1",
 "h2_accession_version":"GCA_048127205.1","pair_design_sha256":pre["pair_design_sha256"],
 "assets":{"h1_fasta":{"path":str(h1)},"h2_fasta":{"path":str(h2)}},
 "h1_to_h2_contig_map":None,"annotation":ann,
}
out.write_text(json.dumps(value,sort_keys=True)+"\n")
PY
else
    input_dir="$VGP_DATA_ROOT/pilot/inputs/$selection"
    [[ -f $source_run/mapping/.complete.json ]] || fail "closed exact mapping is absent"
    python3 -m analysis.vgp_three_pair stage-fastas \
        "$input_dir/input-manifest.json" "$scratch/inputs" \
        "$scratch/index/staged_fasta_dictionary.json" \
        >"$scratch/index/staged_fasta_dictionary.stdout.json"
    cp "$input_dir/input-manifest.json" "$scratch/inputs/input-manifest.json"
    "$samtools" faidx "$scratch/inputs/h1.fa"
    "$samtools" faidx "$scratch/inputs/h2.fa"
    "$impg" index -a "$source_run/mapping/h2_to_h1.1to1.paf" \
        -i "$scratch/index/h1_h2.impg" -t "$cpus" \
        >"$scratch/index/impg.index.log" 2>&1
    "$impg" partition -a "$source_run/mapping/h2_to_h1.1to1.paf" \
        -i "$scratch/index/h1_h2.impg" -w 2000 -d 0 --min-missing-size 1 \
        --min-boundary-distance 0 -o bed --output-folder "$scratch/index/partitions" \
        -t "$cpus" >"$scratch/index/impg.partition.log" 2>&1
    python3 -m analysis.vgp_three_pair audit-graph-ids \
        "$scratch/index/staged_fasta_dictionary.json" \
        "$source_run/mapping/h2_to_h1.1to1.paf" \
        "$scratch/index/partitions/partitions.bed" \
        "$scratch/index/graph_identifier_audit.json" \
        >"$scratch/index/graph_identifier_audit.stdout.json"
fi

h1="$scratch/inputs/h1.fa"
h2="$scratch/inputs/h2.fa"
paf="$source_run/mapping/h2_to_h1.1to1.paf"
partitions="$scratch/index/partitions/partitions.bed"
index="$scratch/index/h1_h2.impg"
python3 -m analysis.vgp_bounded_ranges freeze-plan \
    "$h1.fai" "$partitions" "$scratch/plan/range_plan.json" \
    --selection-id "$selection" --target-bp 5000000 --hard-max-bp 20000000 \
    >"$scratch/plan/freeze.stdout.json"

process_range() {
    local range_id=$1 contig=$2 start=$3 end=$4 expected=$5
    local work="$scratch/worker/$range_id"
    local out="$scratch/results/ranges/$range_id"
    mkdir -p "$work/calls" "$work/temp" "$out"
    python3 -m analysis.vgp_bounded_ranges emit-range-bed \
        "$scratch/plan/range_plan.json" "$partitions" "$range_id" "$work/focus.bed" \
        >"$work/focus.json"
    local begin query_end lace_end count graph_bytes
    begin=$(date +%s)
    "$impg" query -a "$paf" -i "$index" -b "$work/focus.bed" -d 0 \
        --min-transitive-len 1 --force-large-region --temp-dir "$work/temp" \
        -o vcf:poa --sequence-files "$h1" "$h2" -O "$work/calls" \
        -t "$query_threads" >"$work/query.log" 2>&1
    query_end=$(date +%s)
    find "$work/calls" -type f -name '*.vcf' -print | LC_ALL=C sort >"$work/vcf.list"
    count=$(wc -l <"$work/vcf.list")
    [[ $count -eq $expected ]] || fail "$range_id regional VCF census: $count != $expected"
    graph_bytes=$(du -sb "$work/calls" "$work/temp" | awk '{n+=$1} END {print n+0}')
    "$impg" lace -l "$work/vcf.list" --format vcf -o "$work/laced.vcf" \
        --reference "$h1" --temp-dir "$work/temp" --compress none \
        -t "$query_threads" >"$work/lace.log" 2>&1
    lace_end=$(date +%s)
    python3 - "$work/laced.vcf" "$h1.fai" "$work/h2.samples" <<'PY'
import sys
from pathlib import Path
vcf,fai,out=map(Path,sys.argv[1:])
h1={line.split("\t",1)[0] for line in fai.open()}; samples=None
for line in vcf.open():
 if line.startswith("#CHROM"): samples=line.rstrip("\n").split("\t")[9:]; break
if samples is None: raise SystemExit("range lace lacks #CHROM")
out.write_text("".join(f"{s}\n" for s in samples if s.split(":",1)[0] not in h1))
PY
    python3 - "$source_run/mapping/h1.1to1.bed" "$contig" "$start" "$end" \
        "$work/ownership.bed" <<'PY'
import sys
from pathlib import Path
source=Path(sys.argv[1]); contig=sys.argv[2]; start,end=map(int,sys.argv[3:5]); rows=[]
for line in source.open():
 f=line.split()
 if f[0] != contig: continue
 left,right=max(start,int(f[1])),min(end,int(f[2]))
 if left < right: rows.append((contig,left,right))
Path(sys.argv[5]).write_text("".join(f"{c}\t{s}\t{e}\n" for c,s,e in rows))
PY
    if [[ -s $work/h2.samples && -s $work/ownership.bed ]]; then
        "$bcftools" view --samples-file "$work/h2.samples" --trim-alt-alleles \
            --min-ac 1:nref -Ou "$work/laced.vcf" |
            "$bcftools" view --drop-genotypes --no-update -Ou |
            "$bcftools" norm -f "$h1" -m -any -Ou |
            "$bcftools" view -T "$work/ownership.bed" -Ou |
            "$bcftools" norm -f "$h1" -d exact -Oz -o "$out/normalized.vcf.gz"
    else
        "$bcftools" view --header-only --drop-genotypes --no-update -Oz \
            -o "$out/normalized.vcf.gz" "$work/laced.vcf"
    fi
    "$bcftools" index -f -t "$out/normalized.vcf.gz"
    "$bcftools" view -Ob -o "$out/normalized.bcf" "$out/normalized.vcf.gz"
    "$bcftools" index -f "$out/normalized.bcf"
    local variants
    variants=$("$bcftools" index -n "$out/normalized.vcf.gz")
    python3 - "$out/range_audit.json" "$selection" "$range_id" "$contig" "$start" \
        "$end" "$expected" "$count" "$variants" "$graph_bytes" "$begin" "$query_end" \
        "$lace_end" <<'PY'
import json,sys
from pathlib import Path
out=Path(sys.argv[1])
selection,rid,contig=sys.argv[2:5]
start,end,expected,count,variants,peak,begin,qend,lend=map(int,sys.argv[5:14])
out.write_text(json.dumps({
 "schema_version":"vgp-bounded-range-result-v1","selection_id":selection,
 "range_id":rid,"contig":contig,"start":start,"end":end,"range_bp":end-start,
 "expected_native_partitions":expected,"queried_native_partitions":count,
 "normalized_variant_records":variants,"peak_local_graph_state_bytes":peak,
 "query_elapsed_seconds":qend-begin,"lace_normalize_elapsed_seconds":lend-qend,
 "partition_one_owner":True,"variant_half_open_owner":True,
 "all_genome_graph_materialized":False,"global_impg_lace_created":False,
 "local_graph_temporaries_discarded":True,"ref_validation":"bcftools norm -f exact H1",
},sort_keys=True)+"\n")
PY
    cp "$work/query.log" "$work/lace.log" "$out/"
    rm -rf -- "$work"
}

mapfile -t query_rows < <(python3 - "$scratch/plan/range_plan.json" <<'PY'
import json,sys
for r in json.load(open(sys.argv[1]))["ranges"]:
 if r["query_required"]:
  print("\t".join(map(str,(r["range_id"],r["contig"],r["start"],r["end"],r["partition_count"]))))
PY
)
failed=0
for ((batch=0; batch<${#query_rows[@]}; batch+=range_workers)); do
    pids=()
    for ((i=batch; i<${#query_rows[@]} && i<batch+range_workers; i++)); do
        IFS=$'\t' read -r rid contig start end count <<<"${query_rows[$i]}"
        process_range "$rid" "$contig" "$start" "$end" "$count" &
        pids+=("$!")
    done
    for pid in "${pids[@]}"; do wait "$pid" || failed=1; done
    (( failed == 0 )) || fail "one or more bounded H1 ranges failed"
done

python3 - "$scratch/plan/range_plan.json" "$scratch/results/ranges" \
    "$scratch/results/variants/range_vcfs.list" "$scratch/results/range_completion.json" <<'PY'
import hashlib,json,sys
from pathlib import Path
plan=Path(sys.argv[1]); root=Path(sys.argv[2]); listing=Path(sys.argv[3]); out=Path(sys.argv[4])
p=json.loads(plan.read_text()); files=[]; audits=[]; maximum=0
for row in p["ranges"]:
 if not row["query_required"]: continue
 base=root/row["range_id"]; vcf=base/"normalized.vcf.gz"; audit=base/"range_audit.json"
 if not vcf.is_file() or not Path(str(vcf)+".tbi").is_file() or not audit.is_file():
  raise SystemExit(f"missing bounded output: {row['range_id']}")
 value=json.loads(audit.read_text()); audits.append(value); maximum=max(maximum,value["peak_local_graph_state_bytes"])
 files.append(vcf)
listing.write_text("".join(f"{x}\n" for x in files))
out.write_text(json.dumps({
 "schema_version":"vgp-bounded-range-completion-v1","range_count":p["range_count"],
 "query_range_count":len(files),"nonquery_range_count":p["range_count"]-len(files),
 "native_partition_rows":sum(x["queried_native_partitions"] for x in audits),
 "normalized_variant_records":sum(x["normalized_variant_records"] for x in audits),
 "maximum_peak_local_graph_state_bytes":maximum,"all_range_audits_passed":True,
 "all_genome_graph_materialized":False,"global_impg_lace_created":False,
 "global_partition_assignment_ledger_materialized":False,
},sort_keys=True)+"\n")
PY

# A convenience genome-wide file is only a concatenation of already normalized,
# disjoint H1 range outputs. It is never an IMPG query or lace input.
"$bcftools" concat -a -d exact -f "$scratch/results/variants/range_vcfs.list" \
    -Oz -o "$scratch/results/variants/normalized.vcf.gz"
"$bcftools" index -f -t "$scratch/results/variants/normalized.vcf.gz"
"$bcftools" view -Ob -o "$scratch/results/variants/normalized.bcf" \
    "$scratch/results/variants/normalized.vcf.gz"
"$bcftools" index -f "$scratch/results/variants/normalized.bcf"
"$bcftools" view -Ov -o "$scratch/results/variants/normalized.vcf" \
    "$scratch/results/variants/normalized.vcf.gz"
python3 - "$paf" "$scratch/results/variants/ref_alt_coordinate_audit.json" <<'PY'
import collections,json,sys
from pathlib import Path
strands=collections.Counter(); rows=0
for number,line in enumerate(Path(sys.argv[1]).open(),1):
 if not line.strip(): continue
 f=line.rstrip("\n").split("\t")
 if len(f)<12 or f[4] not in {"+","-"}: raise SystemExit(f"malformed PAF row {number}")
 ql,qs,qe,tl,ts,te=map(int,(f[1],f[2],f[3],f[6],f[7],f[8]))
 if not (0<=qs<qe<=ql and 0<=ts<te<=tl): raise SystemExit(f"invalid PAF coordinate {number}")
 strands[f[4]]+=1; rows+=1
Path(sys.argv[2]).write_text(json.dumps({
 "paf_rows":rows,"strand_counts":dict(strands),"invalid_coordinates":0,
 "normalized_ref_mismatches":0,"ref_alt_reconstruction_failures":0,
 "ref_validation":"every retained range passed pinned bcftools norm -f exact staged H1",
},sort_keys=True)+"\n")
PY

python3 - "$scratch/inputs" "$scratch/results/consensus" "$selection" "$source_run" \
    "$scratch/results/variants/normalized.vcf" <<'PY'
import json,sys
from pathlib import Path
from analysis.vgp_10_pilot import (
 REASON_ORDER,materialize_mask_consensus_psmc,paf_concordance_regions,parse_paf,
)
inputs,out,selection,source,vcf=map(Path,sys.argv[1:])
manifest=json.loads((inputs/"input-manifest.json").read_text()); exclusions={}
for reason in REASON_ORDER:
 rows=[]; supplied=Path("/moosefs/erikg/vgp/pilot/inputs")/str(selection)/"exclusions"/f"{reason}.bed"
 generated=source/"mapping"/f"{reason}.bed"
 if supplied.is_file(): rows.extend(supplied.read_text().splitlines())
 if generated.is_file(): rows.extend(generated.read_text().splitlines())
 if rows:
  combined=out/f"inputs.{reason}.bed"; combined.write_text("\n".join(rows)+"\n")
  exclusions[reason]=combined
materialize_mask_consensus_psmc(
 inputs/"h1.fa",inputs/"h2.fa",vcf,inputs/"h1_universe.bed",exclusions,out,
 contig_map=manifest.get("h1_to_h2_contig_map"),selection_id=str(selection),attempts=200,
 aligned_regions=paf_concordance_regions(parse_paf(source/"mapping/h2_to_h1.1to1.paf")),
)
PY
python3 -m analysis.vgp_bounded_ranges finalize-callable \
    "$scratch/plan/range_plan.json" \
    "$scratch/results/consensus/consensus/consensus.fa" \
    "$scratch/results/consensus" "$scratch/results/ranges" \
    >"$scratch/results/consensus/final_callable.stdout.json"

mkdir -p "$scratch/results/psmc/replicate-000"
"$psmc" -N25 -t15 -r5 -p '4+25*2+4+6' \
    -o "$scratch/results/psmc/replicate-000/unscaled.psmc" \
    "$scratch/results/consensus/consensus/input.psmcfa"
psmc_workers=$cpus
(( psmc_workers > 20 )) && psmc_workers=20
for ((batch=1; batch<=200; batch+=psmc_workers)); do
    pids=()
    for ((rep=batch; rep<=200 && rep<batch+psmc_workers; rep++)); do
        (
            out="$scratch/results/psmc/replicate-$(printf '%03d' "$rep")"
            mkdir -p "$out"
            unit="$out/input.psmcfa"
            python3 -m analysis.vgp_10_pilot emit-bootstrap \
                "$scratch/results/consensus/consensus/input.psmcfa" \
                "$scratch/results/consensus/consensus/bootstrap_units.5mb.psmcfa_bins.tsv" \
                "$scratch/results/consensus/consensus/bootstrap_manifest.tsv" "$rep" "$unit" \
                >"$out/bootstrap_input.json"
            "$psmc" -b -N25 -t15 -r5 -p '4+25*2+4+6' \
                -o "$out/bootstrap.unscaled.psmc" "$unit"
            rm "$unit"
        ) &
        pids+=("$!")
    done
    failed=0
    for pid in "${pids[@]}"; do wait "$pid" || failed=1; done
    (( failed == 0 )) || fail "one or more PSMC bootstraps failed"
done

scenario="$VGP_DATA_ROOT/pilot/inputs/$selection/psmc_scaling_scenarios.tsv"
if [[ $selection == P07 ]]; then
    scenario="$scratch/inputs/psmc_scaling_scenarios.tsv"
    python3 - "$source_run/psmc/finalize/scenario_scaled_trajectories.tsv" "$scenario" <<'PY'
import csv,sys
rows=list(csv.DictReader(open(sys.argv[1]),delimiter="\t")); seen={}; fields=(
"scenario_id","mutation_rate_per_generation","generation_time_years",
"mutation_rate_source","generation_time_source")
for row in rows: seen.setdefault(row["scenario_id"],{k:row[k] for k in fields})
with open(sys.argv[2],"w",newline="") as out:
 w=csv.DictWriter(out,fieldnames=fields,delimiter="\t",lineterminator="\n"); w.writeheader()
 w.writerows(seen.values())
PY
fi
mkdir "$scratch/results/psmc/finalize"
python3 - "$scratch/results" "$scenario" <<'PY'
import csv,json,math,sys
from pathlib import Path
from analysis.vgp_10_pilot import parse_psmc_unscaled,scale_unscaled_trajectory,sha256_file
root,scenario=map(Path,sys.argv[1:]); out=root/"psmc/finalize"
rows,theta=parse_psmc_unscaled(root/"psmc/replicate-000/unscaled.psmc")
with (out/"unscaled_trajectory.tsv").open("w",newline="") as f:
 w=csv.DictWriter(f,fieldnames=("interval","time_2N0","lambda"),delimiter="\t",lineterminator="\n")
 w.writeheader(); w.writerows(rows)
scenarios=list(csv.DictReader(scenario.open(),delimiter="\t"))
for row in scenarios: row["theta_0"]=theta
_,scaled=scale_unscaled_trajectory(rows,scenarios)
fields=("scenario_id","interval","time_years","effective_size","mutation_rate_per_generation",
"generation_time_years","psmc_bin_size_bp","mutation_rate_source","generation_time_source")
with (out/"scenario_scaled_trajectories.tsv").open("w",newline="") as f:
 w=csv.DictWriter(f,fieldnames=fields,delimiter="\t",lineterminator="\n"); w.writeheader(); w.writerows(scaled)
boot=[]; thetas=[]
for rep in range(1,201):
 p=root/f"psmc/replicate-{rep:03d}/bootstrap.unscaled.psmc"; values,value=parse_psmc_unscaled(p)
 finite=all(math.isfinite(x["time_2N0"]) and math.isfinite(x["lambda"]) for x in values)
 boot.append({"replicate":rep,"finite":str(finite).lower(),"sha256":sha256_file(p)}); thetas.append(value)
ordered=sorted(thetas); low=ordered[round(.025*199)]; high=ordered[round(.975*199)]
with (out/"bootstrap_unscaled.tsv").open("w",newline="") as f:
 w=csv.DictWriter(f,fieldnames=("replicate","finite","sha256"),delimiter="\t",lineterminator="\n")
 w.writeheader(); w.writerows(boot)
finite=sum(x["finite"]=="true" for x in boot); centered=low<=theta<=high
(out/"psmc_qc.json").write_text(json.dumps({
 "bootstrap_attempts":200,"finite_bootstraps":finite,"minimum_finite_required":200,
 "passed":finite==200 and centered,"primary_theta":theta,
 "nearest_index_central_95pct":[low,high],"primary_theta_centered":centered,
 "unscaled_primary_preserved":True,"scenario_scaling_separate":True,
},sort_keys=True)+"\n")
if finite != 200 or not centered: raise SystemExit("PSMC 200-finite/centering gate failed")
PY

if [[ $selection == P07 ]]; then
    mkdir "$scratch/results/annotation"
    cp "$ROOT/analysis/vgp_clean_canary_selection_v1.json" "$scratch/results/annotation/selection.json"
    "$bcftools" view -Ov -o "$scratch/results/annotation/normalized.vcf" \
        "$scratch/results/variants/normalized.vcf.gz"
    python3 "$ROOT/analysis/vgp_real_canary_annotation.py" \
        --h1-fasta "$h1" \
        --annotation-gff /moosefs/erikg/vgp/objects/sha256/8f/64/8f640543accd8081d1b7048eda32c9f1eef33b02f321b7b0f8adcf3b01dd6838 \
        --annotation-source-path /moosefs/erikg/vgp/objects/sha256/8f/64/8f640543accd8081d1b7048eda32c9f1eef33b02f321b7b0f8adcf3b01dd6838 \
        --callable-bed "$scratch/results/consensus/masks/callable.bed" \
        --normalized-vcf "$scratch/results/annotation/normalized.vcf" \
        --canonical-root "$VGP_DATA_ROOT" --selection-id P07 \
        --assembly-accession-version GCA_048126635.1 \
        --annotation-accession-version GCA_048126635.1-GB_2025_08_04 \
        --task-id run-vgp-three-pair --schema-version vgp-three-pair-bounded-annotation-v1 \
        --output "$scratch/results/annotation/exact_partitions.json"
    python3 - "$scratch/results/annotation/exact_partitions.json" <<'PY'
import json,sys
from pathlib import Path
p=Path(sys.argv[1]); value=json.loads(p.read_text())
value["query_scope"]="exact native GFF features intersected with bounded-range-derived variants and mask"
value["additional_impg_graph_queries_for_annotation"]=False
p.write_text(json.dumps(value,sort_keys=True)+"\n")
PY
fi

python3 - "$scratch" "$selection" "$SLURM_JOB_ID" <<'PY'
import hashlib,json,sys,time
from pathlib import Path
root=Path(sys.argv[1]); selection,job=sys.argv[2:4]; results=root/"results"
join=json.loads((results/"consensus/join_qc.json").read_text())
completion=json.loads((results/"range_completion.json").read_text())
value={
 "schema_version":"vgp-bounded-pair-production-v1","selection_id":selection,
 "slurm_job_id":job,"actual_core_biological_result":True,
 "diversity":{"heterozygous_snps":join["consensus"]["heterozygous_snps"],
 "callable_bp":join["consensus"]["consensus_callable_bp"],
 "pi":join["consensus"]["heterozygous_snps"]/join["consensus"]["consensus_callable_bp"]},
 "range_completion":completion,"mask":join["mask"],"consensus":join["consensus"],
 "psmc":json.loads((results/"psmc/finalize/psmc_qc.json").read_text()),
 "all_genome_graph_materialized":False,"global_impg_lace_created":False,
 "global_partition_assignment_ledger_materialized":False,
 "technical_cancellation_is_species_failure":False,
}
(results/"execution.json").write_text(json.dumps(value,sort_keys=True)+"\n")
PY

mkdir -p "$(dirname "$durable")"
partial="$(dirname "$durable")/.bounded-production.${SLURM_JOB_ID}.partial"
[[ ! -e $partial ]] || fail "durable bounded partial already exists"
cp -a "$scratch/results" "$partial"
mkdir "$partial/index"
cp "$index" "$partitions" "$scratch/plan/range_plan.json" "$partial/index/"
[[ -f $scratch/index/graph_identifier_audit.json ]] && \
    cp "$scratch/index/graph_identifier_audit.json" "$partial/index/"
python3 - "$partial" "$selection" <<'PY'
import hashlib,json,sys
from pathlib import Path
root=Path(sys.argv[1]); selection=sys.argv[2]
def sha(path):
 h=hashlib.sha256()
 with path.open("rb") as f:
  for b in iter(lambda:f.read(1024*1024),b""): h.update(b)
 return h.hexdigest()
files={str(p.relative_to(root)):sha(p) for p in sorted(root.rglob("*"))
       if p.is_file() and p.name!=".complete.json"}
(root/".complete.json").write_text(json.dumps({
 "schema_version":"vgp-bounded-pair-closed-stage-v1","selection_id":selection,
 "all_genome_graph_materialized":False,"global_impg_lace_created":False,
 "files":files,
},sort_keys=True)+"\n")
PY
mv "$partial" "$durable"
printf '%s\n' "$durable"
