#!/usr/bin/env bash
#SBATCH --job-name=tier3a-origin-rerun
#SBATCH --array=0-2
#SBATCH --partition=workers
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --export=NONE
#SBATCH --output=results/tier3a/logs/origin-rerun-%A_%a.out
#SBATCH --error=results/tier3a/logs/origin-rerun-%A_%a.err

set -euo pipefail

if [[ -n ${SLURM_SUBMIT_DIR:-} ]]; then
    ROOT=$(CDPATH= cd -- "$SLURM_SUBMIT_DIR" && pwd)
else
    ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
fi
CHANNELS="$ROOT/analysis/guix/channels.scm"
GUIX_MANIFEST="$ROOT/analysis/guix/manifest.scm"
ENVIRONMENT_RECORD="$ROOT/analysis/pilot_results/guix_environment.json"
SUPPLEMENTAL_PROFILE=/gnu/store/8x4hx7d9hnv187yprjrzqyg0kxj2z32k-profile
ACQUISITION="$ROOT/results/tier3a/acquisition_corrected_manifest.tsv"
CORRECTED_COMMANDS="$ROOT/results/tier3a/acquisition_corrected_commands.tsv"
SUPERSESSION="$ROOT/results/tier3a/acquisition_sweepga_supersession.tsv"
BUILD_PROVENANCE="$ROOT/analysis/sweepga_origin_main_build.json"
WORK_ROOT=${TIER3A_WORK_ROOT:-/moosefs/erikg/tier3data/tier3a-origin-biological-rerun-20260716}
PANEL_BASES=${TIER3A_PANEL_BASES:-2000000}

if [[ ${1:-} != --inside-guix ]]; then
    export TIER3_ROOT="$ROOT"
    exec "$ROOT/analysis/slurm/guix_job.sh" "$ENVIRONMENT_RECORD" \
        bash "$ROOT/analysis/slurm/tier3a_biological_array.sh" --inside-guix
fi

[[ -d $SUPPLEMENTAL_PROFILE ]] || { echo "Pinned supplemental Guix profile is unavailable" >&2; exit 69; }
export PATH="$SUPPLEMENTAL_PROFILE/bin:$PATH"

PREFLIGHT="$WORK_ROOT/preflight.${SLURM_ARRAY_TASK_ID:-0}.json"
mkdir -p "$WORK_ROOT"
python3 "$ROOT/analysis/tier3a_biological.py" preflight \
    --manifest "$ACQUISITION" --commands "$CORRECTED_COMMANDS" \
    --supersession-ledger "$SUPERSESSION" --build-provenance "$BUILD_PROVENANCE" \
    --output "$PREFLIGHT"

INDEX=${SLURM_ARRAY_TASK_ID:-0}
mapfile -t META < <(python3 - "$ACQUISITION" "$INDEX" <<'PY'
import csv
import sys
csv.field_size_limit(sys.maxsize)
with open(sys.argv[1], encoding="utf-8", newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
index = int(sys.argv[2])
if not 0 <= index < len(rows):
    raise SystemExit("array index outside acquisition manifest")
row = rows[index]
for key in ("dataset_id", "sweepga_bounded_paf_path", "h1_fasta_path", "h2_fasta_path", "impg_binary_path", "annotation_gff_path", "annotation_release_date"):
    print(row[key])
PY
)
[[ ${#META[@]} -eq 7 ]] || { echo "Incomplete acquisition metadata at array index $INDEX" >&2; exit 2; }
DATASET=${META[0]}
PAF=${META[1]}
H1=${META[2]}
H2=${META[3]}
IMPG=${META[4]}
GFF=${META[5]}
ANNOTATION_RELEASE=${META[6]}
TOOL_DIR=${IMPG%/*}
export PATH="$TOOL_DIR:$PATH"
OUT="$WORK_ROOT/$DATASET"
[[ ! -e $OUT ]] || { echo "Refusing to reuse existing rerun tuple directory: $OUT" >&2; exit 70; }
mkdir -p "$OUT/native_annotation" "$OUT/partitions" "$OUT/calls"
cp "$PREFLIGHT" "$OUT/preflight_audit.json"

python3 "$ROOT/results/tier3a/acquisition_build_queries.py" \
    --dataset-id "$DATASET" --gff "$GFF" --fai "$H1.fai" --h2-fai "$H2.fai" \
    --bounded-paf "$PAF" --output-dir "$OUT/native_annotation" \
    --annotation-release "$ANNOTATION_RELEASE" \
    --mapping-status origin_main_native_1to1_corrected_mapping \
    > "$OUT/native_annotation/build.json"

python3 "$ROOT/analysis/tier3a_biological.py" prepare \
    --manifest "$ACQUISITION" --dataset-id "$DATASET" --output-dir "$OUT" \
    --annotation-dir "$OUT/native_annotation" --preflight-audit "$OUT/preflight_audit.json" \
    --panel-bases "$PANEL_BASES"

START_EPOCH=$(date +%s)
"$IMPG" index -a "$PAF" -i "$OUT/graph.impg" -t "$SLURM_CPUS_PER_TASK" \
    > "$OUT/index.stdout" 2> "$OUT/index.stderr"
"$IMPG" partition -a "$PAF" -i "$OUT/graph.impg" -w 2000 -d 0 \
    --min-missing-size 1 --min-boundary-distance 0 -o bed \
    --output-folder "$OUT/partitions" -t "$SLURM_CPUS_PER_TASK" -v 0 \
    > "$OUT/partition.stdout" 2> "$OUT/partition.stderr"

python3 "$ROOT/results/tier3a/acquisition_select_impg_partitions.py" \
    --partitions "$OUT/partitions/partitions.bed" \
    --targets "$OUT/mapping_callable.bed" \
    --focus-bed "$OUT/focus.bed" --mapping-tsv "$OUT/partition_annotation_map.tsv"

"$IMPG" query -a "$PAF" -i "$OUT/graph.impg" -b "$OUT/focus.bed" -d 0 \
    -o vcf:poa --force-large-region --min-transitive-len 1 \
    --sequence-files "$H1" "$H2" -O "$OUT/calls" \
    -t "$SLURM_CPUS_PER_TASK" > "$OUT/query.stdout" 2> "$OUT/query.stderr"
find "$OUT/calls" -name '*.vcf' -type f | sort > "$OUT/vcf.list"
[[ -s $OUT/vcf.list ]] || { echo "IMPG emitted no regional VCFs for $DATASET" >&2; exit 3; }

"$IMPG" lace -l "$OUT/vcf.list" --format vcf -o "$OUT/laced.vcf" \
    --reference "$H1" --compress none -t "$SLURM_CPUS_PER_TASK" \
    > "$OUT/lace.stdout" 2> "$OUT/lace.stderr"
bcftools norm -f "$H1" -m -any -Ob -o "$OUT/normalized.untrimmed.bcf" "$OUT/laced.vcf" \
    > "$OUT/norm.stdout" 2> "$OUT/norm.stderr"
bcftools index -f "$OUT/normalized.untrimmed.bcf"
bcftools view -R "$OUT/mapping_callable.bed" -Ou "$OUT/normalized.untrimmed.bcf" | \
    bcftools norm -d exact -Ob -o "$OUT/normalized.bcf"
bcftools index -f "$OUT/normalized.bcf"

python3 "$ROOT/analysis/tier3a_biological.py" summarize \
    --output-dir "$OUT" --focus-bed "$OUT/focus.bed" \
    --partitions "$OUT/partitions/partitions.bed" --index "$OUT/graph.impg" \
    --bcf "$OUT/normalized.bcf" --bcftools "$(command -v bcftools)" \
    --block-size 50000 --bootstrap-replicates 1000

END_EPOCH=$(date +%s)
printf '%s\t%s\t%s\t%s\n' "$DATASET" "$START_EPOCH" "$END_EPOCH" "$GUIX_ENVIRONMENT" > "$OUT/completed.tsv"
