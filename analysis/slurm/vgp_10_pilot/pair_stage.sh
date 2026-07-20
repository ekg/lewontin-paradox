#!/usr/bin/env bash
#SBATCH --job-name=vgp10-pair
# Scheduler time/CPU/memory/scratch limits are supplied per job by submit.sh
# from the measured estimator record.  There is intentionally no global cap.
set -euo pipefail

if [[ -n ${SLURM_JOB_ID:-} ]]; then
    VGP_STAGE_REPO_ROOT=${SLURM_SUBMIT_DIR:?submit from the repository root}
else
    VGP_STAGE_REPO_ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/../../.." && pwd)
fi
source "$VGP_STAGE_REPO_ROOT/analysis/slurm/vgp_10_pilot/common.sh"
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

# FastGA/IMPG create indexes and companion files relative to their sequence
# inputs.  Canonical inputs are immutable, so every biological stage receives
# private node-local working copies sized by resources.json.
if [[ $VGP_STAGE_NAME != preflight ]]; then
    mkdir -p "$SLURM_TMPDIR/inputs"
    cp "$h1" "$SLURM_TMPDIR/inputs/h1.fa"
    cp "$h2" "$SLURM_TMPDIR/inputs/h2.fa"
    h1="$SLURM_TMPDIR/inputs/h1.fa"
    h2="$SLURM_TMPDIR/inputs/h2.fa"
fi

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
    : "${VGP_FASTGA_AMENDMENT:=$VGP_REPO_ROOT/analysis/vgp_real_pilot_fastga_amendment_v1.json}"
    readarray -t fastga_companions < <(python3 - "$VGP_FASTGA_AMENDMENT" <<'PY'
import hashlib,json,sys
from pathlib import Path
value=json.load(open(sys.argv[1])); companions=value.get("companions",{})
required=("ALNtoPAF","FAtoGDB","FastGA","GIXmake","GIXrm","ONEview","PAFtoALN","wfmash")
for name in required:
    row=companions.get(name,{}); path=Path(row.get("path",""))
    if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest()!=row.get("sha256"):
        raise SystemExit(f"FastGA amendment digest mismatch: {name}")
print(companions["FAtoGDB"]["path"]); print(companions["GIXmake"]["path"])
PY
)
    fatogdb=${fastga_companions[0]}
    gixmake=${fastga_companions[1]}
    for fasta in "$h2" "$h1"; do
        "$fatogdb" "$fasta"
        "$gixmake" -T"$threads" -P"$SLURM_TMPDIR" "${fasta%.fa}"
    done
    python3 - "$VGP_STAGE_PARTIAL/fastga_sidecar_prebuild.json" "$fatogdb" "$gixmake" \
        "$VGP_SELECTION_ID" "$VGP_DATA_ROOT" "$VGP_FASTGA_AMENDMENT" <<'PY'
import hashlib,json,sys
from pathlib import Path
def digest(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
Path(sys.argv[1]).write_text(json.dumps({
 "schema_version":"vgp-fastga-sidecar-prebuild-v1","selection_id":sys.argv[4],
 "canonical_vgp_root":sys.argv[5],"reason":"avoid nested child signal observed in scale-out",
 "sweepga_remains_mapping_owner":True,"num_mappings":"1:1",
 "FAtoGDB":{"path":sys.argv[2],"sha256":digest(sys.argv[2])},
 "GIXmake":{"path":sys.argv[3],"sha256":digest(sys.argv[3])},
 "runtime_amendment":{"path":sys.argv[6],"sha256":digest(sys.argv[6])},
},sort_keys=True)+"\n")
PY
    "$sweepga" "$h2" "$h1" --output-file "$VGP_STAGE_PARTIAL/h2_to_h1.native.1to1.paf" \
        --num-mappings 1:1 --scaffold-jump 0 --overlap 0 \
        --scoring log-length-ani --threads "$threads" \
        >"$VGP_STAGE_PARTIAL/sweepga.stdout" 2>"$VGP_STAGE_PARTIAL/sweepga.stderr"
    python3 -m analysis.vgp_10_pilot enforce-paf \
        "$VGP_STAGE_PARTIAL/h2_to_h1.native.1to1.paf" \
        "$VGP_STAGE_PARTIAL/h2_to_h1.1to1.paf" \
        >"$VGP_STAGE_PARTIAL/exact_multiplicity_filter.json"
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
    # IMPG emits one INFO row per 2-kb partition (hundreds of thousands to
    # millions per genome).  Preserve warnings/errors and structured audits
    # without turning scheduler stderr into a second large biological object.
    export RUST_LOG=${VGP_IMPG_LOG_LEVEL:-warn}
    verify_tool impg
    verify_tool bcftools
    verify_tool samtools
    impg=$(tool_path impg)
    bcftools=$(tool_path bcftools)
    samtools=$(tool_path samtools)
    paf="$VGP_PAIR_RUN/mapping/h2_to_h1.1to1.paf"
    "$impg" index -a "$paf" -i "$VGP_STAGE_PARTIAL/h1_h2.impg" -t "$threads"
    mkdir "$VGP_STAGE_PARTIAL/partitions"
    "$impg" partition -a "$paf" -i "$VGP_STAGE_PARTIAL/h1_h2.impg" \
        -w 2000 -d 0 --min-missing-size 1 --min-boundary-distance 0 \
        -o bed --output-folder "$VGP_STAGE_PARTIAL/partitions" -t "$threads" \
        2> >(awk '!/ INFO  /' >"$VGP_STAGE_PARTIAL/impg.partition.stderr")
    python3 - "$VGP_STAGE_PARTIAL/partitions/partitions.bed" \
        "$input_dir/h1_universe.bed" "$VGP_STAGE_PARTIAL/focus.native.bed" <<'PY'
import sys
from analysis.vgp_10_pilot import read_bed,select_native_partitions
select_native_partitions(sys.argv[1],read_bed(sys.argv[2]),sys.argv[3])
PY
    # Every IMPG query process otherwise races to create the same .fai files.
    # Build the private staged indexes once, before starting concurrent readers.
    "$samtools" faidx "$h1"
    "$samtools" faidx "$h2"
    [[ -s $h1.fai && -s $h2.fai ]] || fail "staged FASTA indexes are absent"
    source "$VGP_REPO_ROOT/analysis/slurm/vgp_10_pilot/impg_parallel_query.sh"
    find "$VGP_STAGE_PARTIAL/calls" -type f -name '*.vcf' -print | LC_ALL=C sort \
        >"$VGP_STAGE_PARTIAL/vcf.list"
    [[ -s $VGP_STAGE_PARTIAL/vcf.list ]] || fail "IMPG query emitted no regional VCF"
    python3 - "$VGP_STAGE_PARTIAL/focus.native.bed" "$VGP_STAGE_PARTIAL/vcf.list" \
        "$VGP_STAGE_PARTIAL/regional_vcf_audit.json" "$VGP_SELECTION_ID" \
        "$VGP_DATA_ROOT" <<'PY'
import json,sys
from pathlib import Path
focus,listing,output,selection,canonical_root=(
    Path(sys.argv[1]),Path(sys.argv[2]),Path(sys.argv[3]),sys.argv[4],sys.argv[5])
rows=[line.rstrip("\n").split("\t") for line in focus.open() if line.strip() and not line.startswith("#")]
if not rows or any(len(row) < 5 for row in rows):
    raise SystemExit("focused IMPG BED lacks query/native identifiers")
query_names=[row[3] for row in rows]
if len(query_names) != len(set(query_names)):
    raise SystemExit("focused IMPG query names are not unique")
files=[Path(line.strip()) for line in listing.open() if line.strip()]
if len(files) != len(rows):
    raise SystemExit(f"IMPG regional VCF census mismatch: {len(files)} != {len(rows)}")
missing=[str(path) for path in files if not path.is_file() or path.stat().st_size == 0]
if missing:
    raise SystemExit(f"IMPG regional VCFs missing or empty: {missing[:3]}")
file_names={path.stem for path in files}
if file_names != set(query_names):
    raise SystemExit("IMPG regional VCF names do not match focused query rows")
value={
    "canonical_vgp_root":canonical_root,
    "selection_id":selection,
    "focus_rows":len(rows),
    "unique_query_names":len(set(query_names)),
    "unique_native_partition_ids":len({row[4] for row in rows}),
    "regional_vcf_count":len(files),
    "regional_vcf_total_bytes":sum(path.stat().st_size for path in files),
    "all_regional_vcfs_nonempty":True,
    "transient_shards_removed_after_lacing":True,
}
output.write_text(json.dumps(value,sort_keys=True)+"\n")
PY
    # The P07 canary established that a single monolithic lace can deadlock on
    # thousands of regional VCFs.  Preserve every query shard through a
    # parallel hierarchical lace and exact boundary reconciliation instead.
    source "$VGP_REPO_ROOT/analysis/slurm/vgp_10_pilot/impg_hierarchical_lace.sh"
    ;;
variants)
    [[ -f $VGP_PAIR_RUN/impg/.complete.json ]] || fail "IMPG sentinel absent"
    verify_tool bcftools
    bcftools=$(tool_path bcftools)
    "$bcftools" norm -f "$h1" -m -any -Oz -o "$VGP_STAGE_PARTIAL/impg.split.vcf.gz" \
        "$VGP_PAIR_RUN/impg/laced.vcf"
    "$bcftools" index -f -t "$VGP_STAGE_PARTIAL/impg.split.vcf.gz"
    "$bcftools" view -R "$VGP_PAIR_RUN/mapping/h1.1to1.bed" -Oz \
        -o "$VGP_STAGE_PARTIAL/impg.trimmed.vcf.gz" "$VGP_STAGE_PARTIAL/impg.split.vcf.gz"
    "$bcftools" norm -f "$h1" -d exact -Oz -o "$VGP_STAGE_PARTIAL/impg.normalized.vcf.gz" \
        "$VGP_STAGE_PARTIAL/impg.trimmed.vcf.gz"
    "$bcftools" index -f -t "$VGP_STAGE_PARTIAL/impg.normalized.vcf.gz"
    python3 -m analysis.vgp_10_pilot paf-vcf \
        "$VGP_PAIR_RUN/mapping/h2_to_h1.1to1.paf" "$h1" "$h2" \
        "$VGP_STAGE_PARTIAL/paf.raw.vcf" "$VGP_STAGE_PARTIAL/paf_variant_audit.json"
    "$bcftools" norm -f "$h1" -m -any -Oz -o "$VGP_STAGE_PARTIAL/split.vcf.gz" \
        "$VGP_STAGE_PARTIAL/paf.raw.vcf"
    "$bcftools" index -f -t "$VGP_STAGE_PARTIAL/split.vcf.gz"
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
from analysis.vgp_10_pilot import (
    REASON_ORDER,materialize_mask_consensus_psmc,paf_concordance_regions,parse_paf,
)
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
    selection_id=selection,attempts=200,
    aligned_regions=paf_concordance_regions(parse_paf(pair_run/"mapping/h2_to_h1.1to1.paf")),
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
