#!/usr/bin/env bash
#SBATCH --job-name=vgp10-pair
# Scheduler time/CPU/memory/scratch limits are supplied per job by submit.sh
# from the measured estimator record.  There is intentionally no global cap.
set -euo pipefail

source "$(dirname -- "$0")/common.sh"
VGP_STAGE_NAME=${1:?usage: pair_stage.sh STAGE SELECTION_ID}
VGP_SELECTION_ID=${2:?usage: pair_stage.sh STAGE SELECTION_ID}
export VGP_STAGE_NAME VGP_SELECTION_ID
require_selection_id "$VGP_SELECTION_ID"
require_runtime
begin_stage "$VGP_STAGE_NAME"

input_dir="$VGP_DATA_ROOT/pilot/inputs/$VGP_SELECTION_ID"
h1="$input_dir/h1.fa"
h2="$input_dir/h2.fa"
threads=${SLURM_CPUS_PER_TASK:-1}

case "$VGP_STAGE_NAME" in
preflight)
    [[ -f $input_dir/input-manifest.json ]] || fail "pair input manifest absent"
    python3 - "$input_dir/input-manifest.json" "$VGP_STAGE_PARTIAL/preflight.json" <<'PY'
import json,sys
from pathlib import Path
from analysis.vgp_10_pilot import canonical_json,validate_pair_input_manifest
value=json.load(open(sys.argv[1]))
result=validate_pair_input_manifest(value)
Path(sys.argv[2]).write_text(canonical_json(result))
PY
    ;;
mapping)
    [[ -f $VGP_PAIR_RUN/preflight/.complete.json ]] || fail "preflight sentinel absent"
    verify_tool sweepga
    sweepga=$(tool_path sweepga)
    "$sweepga" "$h2" "$h1" --output-file "$VGP_STAGE_PARTIAL/h2_to_h1.1to1.paf" \
        --num-mappings 1:1 --scaffold-jump 0 --overlap 0.95 \
        --scoring log-length-ani --threads "$threads" \
        >"$VGP_STAGE_PARTIAL/sweepga.stdout" 2>"$VGP_STAGE_PARTIAL/sweepga.stderr"
    python3 -m analysis.vgp_10_pilot audit-paf \
        "$VGP_STAGE_PARTIAL/h2_to_h1.1to1.paf" "$h1" "$h2" \
        >"$VGP_STAGE_PARTIAL/multiplicity.json"
    python3 - "$VGP_STAGE_PARTIAL/h2_to_h1.1to1.paf" "$h1" "$h2" \
        "$input_dir/h1_universe.bed" "$VGP_STAGE_PARTIAL" <<'PY'
import sys
from pathlib import Path
from analysis.vgp_10_pilot import (
    low_complexity_intervals,non_acgt_intervals,paf_h1_intervals,parse_fasta,parse_paf,
    project_h2_non_acgt_to_h1,read_bed,subtract_intervals,write_bed,
)
paf,h1,h2,universe,out=sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],Path(sys.argv[5])
records=parse_paf(paf)
one_to_one=paf_h1_intervals(records)
h1_sequences=parse_fasta(h1)
write_bed(out/"h1.1to1.bed",one_to_one)
write_bed(out/"not_1to1.bed",subtract_intervals(read_bed(universe),one_to_one))
write_bed(out/"h1_gap_or_N.bed",non_acgt_intervals(h1_sequences))
write_bed(out/"h2_gap_or_N.bed",project_h2_non_acgt_to_h1(records,parse_fasta(h2)))
write_bed(out/"repeat_or_low_complexity_primary.bed",low_complexity_intervals(h1_sequences))
PY
    ;;
impg)
    [[ -f $VGP_PAIR_RUN/mapping/.complete.json ]] || fail "mapping sentinel absent"
    verify_tool impg
    impg=$(tool_path impg)
    paf="$VGP_PAIR_RUN/mapping/h2_to_h1.1to1.paf"
    "$impg" index -a "$paf" -i "$VGP_STAGE_PARTIAL/h1_h2.impg" -t "$threads"
    mkdir "$VGP_STAGE_PARTIAL/partitions"
    "$impg" partition -a "$paf" -i "$VGP_STAGE_PARTIAL/h1_h2.impg" -d 0 \
        -o bed --output-folder "$VGP_STAGE_PARTIAL/partitions" -t "$threads"
    python3 - "$VGP_STAGE_PARTIAL/partitions/partitions.bed" \
        "$input_dir/eligible_query_regions.bed" "$VGP_STAGE_PARTIAL/focus.native.bed" <<'PY'
import sys
from analysis.vgp_10_pilot import read_bed,select_native_partitions
select_native_partitions(sys.argv[1],read_bed(sys.argv[2]),sys.argv[3])
PY
    mkdir "$VGP_STAGE_PARTIAL/calls"
    "$impg" query -a "$paf" -i "$VGP_STAGE_PARTIAL/h1_h2.impg" \
        -b "$VGP_STAGE_PARTIAL/focus.native.bed" -d 0 -o vcf:poa \
        --sequence-files "$h1" "$h2" -O "$VGP_STAGE_PARTIAL/calls" -t "$threads"
    find "$VGP_STAGE_PARTIAL/calls" -type f -name '*.vcf' -print | LC_ALL=C sort \
        >"$VGP_STAGE_PARTIAL/vcf.list"
    [[ -s $VGP_STAGE_PARTIAL/vcf.list ]] || fail "IMPG query emitted no regional VCF"
    "$impg" lace -l "$VGP_STAGE_PARTIAL/vcf.list" --format vcf \
        -o "$VGP_STAGE_PARTIAL/laced.vcf" --reference "$h1" --compress none -t "$threads"
    ;;
variants)
    [[ -f $VGP_PAIR_RUN/impg/.complete.json ]] || fail "IMPG sentinel absent"
    verify_tool bcftools
    bcftools=$(tool_path bcftools)
    "$bcftools" norm -f "$h1" -m -any -Oz -o "$VGP_STAGE_PARTIAL/split.vcf.gz" \
        "$VGP_PAIR_RUN/impg/laced.vcf"
    "$bcftools" view -R "$VGP_PAIR_RUN/mapping/h1.1to1.bed" -Oz \
        -o "$VGP_STAGE_PARTIAL/trimmed.vcf.gz" "$VGP_STAGE_PARTIAL/split.vcf.gz"
    "$bcftools" norm -f "$h1" -d exact -Oz -o "$VGP_STAGE_PARTIAL/normalized.vcf.gz" \
        "$VGP_STAGE_PARTIAL/trimmed.vcf.gz"
    "$bcftools" index -f -t "$VGP_STAGE_PARTIAL/normalized.vcf.gz"
    "$bcftools" norm -f "$h1" -d exact -Ob -o "$VGP_STAGE_PARTIAL/normalized.bcf" \
        "$VGP_STAGE_PARTIAL/trimmed.vcf.gz"
    "$bcftools" index -f "$VGP_STAGE_PARTIAL/normalized.bcf"
    ;;
consensus)
    [[ -f $VGP_PAIR_RUN/variants/.complete.json ]] || fail "variant sentinel absent"
    verify_tool bcftools
    bcftools=$(tool_path bcftools)
    "$bcftools" view -Ov -o "$VGP_STAGE_PARTIAL/normalized.vcf" \
        "$VGP_PAIR_RUN/variants/normalized.vcf.gz"
    python3 - "$input_dir" "$VGP_STAGE_PARTIAL" "$VGP_SELECTION_ID" "$VGP_PAIR_RUN" <<'PY'
import json,sys
from pathlib import Path
from analysis.vgp_10_pilot import REASON_ORDER,materialize_mask_consensus_psmc
inputs,output,selection,pair_run=Path(sys.argv[1]),Path(sys.argv[2]),sys.argv[3],Path(sys.argv[4])
manifest=json.loads((inputs/"input-manifest.json").read_text())
exclusions={}
derived=pair_run/"mapping"
for reason in REASON_ORDER:
    rows=[]
    supplied=inputs/"exclusions"/f"{reason}.bed"
    generated=derived/f"{reason}.bed"
    if supplied.is_file(): rows.extend(supplied.read_text().splitlines())
    if generated.is_file(): rows.extend(generated.read_text().splitlines())
    if rows:
        combined=output/f"inputs.{reason}.bed"
        combined.write_text("\n".join(rows)+"\n")
        exclusions[reason]=combined
materialize_mask_consensus_psmc(
    inputs/"h1.fa",inputs/"h2.fa",output/"normalized.vcf",inputs/"h1_universe.bed",
    exclusions,output,contig_map=manifest.get("h1_to_h2_contig_map"),
    selection_id=selection,attempts=200,aligned_regions=manifest.get("concordance_regions"),
)
reconciliation=json.loads((output/"masks/mask_reconciliation.json").read_text())
result_gates=manifest.get("result_gates",{})
minimum_bp=int(result_gates.get("minimum_callable_bp",100_000_000))
minimum_fraction=float(result_gates.get("minimum_callable_fraction",0.60))
if reconciliation["callable_bp"] < minimum_bp:
    raise SystemExit(f"callable sequence hard gate failed: {reconciliation['callable_bp']} < {minimum_bp}")
if reconciliation["callable_fraction"] < minimum_fraction:
    raise SystemExit(
        f"callable fraction hard gate failed: {reconciliation['callable_fraction']} < {minimum_fraction}")
PY
    ;;
*)
    fail "unsupported pair stage: $VGP_STAGE_NAME"
    ;;
esac

promote_stage
