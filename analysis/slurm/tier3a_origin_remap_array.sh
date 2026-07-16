#!/usr/bin/env bash
#SBATCH --job-name=t3a-origin-remap
#SBATCH --array=0-2
#SBATCH --partition=workers
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --export=NONE
#SBATCH --output=results/tier3a/logs/origin-remap-%A_%a.out
#SBATCH --error=results/tier3a/logs/origin-remap-%A_%a.err

set -euo pipefail

if [[ -n ${SLURM_SUBMIT_DIR:-} ]]; then
    ROOT=$(CDPATH= cd -- "$SLURM_SUBMIT_DIR" && pwd)
else
    ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
fi
BASE_MANIFEST="$ROOT/results/tier3a/acquisition_manifest.tsv"
ENVIRONMENT_RECORD="$ROOT/analysis/pilot_results/guix_environment.json"
SUPPLEMENTAL_PROFILE=/gnu/store/8x4hx7d9hnv187yprjrzqyg0kxj2z32k-profile
WORK_ROOT=${TIER3A_ORIGIN_WORK_ROOT:-/moosefs/erikg/tier3data/tier3a-origin-remap-20260716}
PINNED_DIR=/moosefs/erikg/tier3scratch/sweepga-origin-main-018e4ce/bin-1
PINNED_SWEEPGA="$PINNED_DIR/sweepga"
EXPECTED_SWEEPGA_SHA256=fa7f0edb9b7e275c288db254046020e136d4267dd5ee043379227ef80da0573b
PANEL_BASES=${TIER3A_PANEL_BASES:-2000000}

if [[ ${1:-} != --inside-guix ]]; then
    export TIER3_ROOT="$ROOT"
    exec "$ROOT/analysis/slurm/guix_job.sh" "$ENVIRONMENT_RECORD" \
        bash "$ROOT/analysis/slurm/tier3a_origin_remap_array.sh" --inside-guix
fi

[[ -d $SUPPLEMENTAL_PROFILE ]] || { echo "Pinned supplemental Guix profile is unavailable" >&2; exit 69; }
export PATH="$SUPPLEMENTAL_PROFILE/bin:$PATH"
export SWEEPGA_PINNED_REALPATH
SWEEPGA_PINNED_REALPATH=$(realpath "$PINNED_SWEEPGA")
[[ $SWEEPGA_PINNED_REALPATH == "$PINNED_SWEEPGA" ]] || { echo "SweepGA realpath changed" >&2; exit 70; }
[[ $(sha256sum "$SWEEPGA_PINNED_REALPATH" | awk '{print $1}') == "$EXPECTED_SWEEPGA_SHA256" ]] || {
    echo "Pinned SweepGA checksum mismatch" >&2
    exit 70
}

INDEX=${SLURM_ARRAY_TASK_ID:-0}
mapfile -t META < <(python3 - "$BASE_MANIFEST" "$INDEX" <<'PY'
import csv, sys
csv.field_size_limit(sys.maxsize)
with open(sys.argv[1], encoding="utf-8", newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
row = rows[int(sys.argv[2])]
for key in ("dataset_id", "h1_fasta_path", "h2_fasta_path", "annotation_gff_path", "annotation_accession_version", "impg_binary_path"):
    print(row[key])
PY
)
[[ ${#META[@]} -eq 6 ]] || { echo "Incomplete acquisition metadata for array index $INDEX" >&2; exit 2; }
DATASET=${META[0]}
H1=${META[1]}
H2=${META[2]}
GFF=${META[3]}
ANNOTATION_RELEASE=${META[4]}
IMPG=${META[5]}
IMPG_DIR=${IMPG%/*}
OUT="$WORK_ROOT/$DATASET"
MAPPING="$OUT/mapping"
QUERY="$OUT/annotation_query"
SCRATCH="$WORK_ROOT/scratch/$DATASET"
START_EPOCH=$(date +%s)

# This correction has an independent root and cannot accidentally resume from
# the superseded acquisition or biological trees.
if [[ -e $OUT/success.json ]]; then
    echo "Refusing to overwrite an already published success marker: $OUT" >&2
    exit 73
fi
install -d -m 700 "$MAPPING" "$QUERY" "$OUT/partitions" "$OUT/calls" "$SCRATCH"

PAF_PART="$MAPPING/production.1to1.paf.part"
COMMAND=(
    "$SWEEPGA_PINNED_REALPATH" "$H1" "$H2"
    --output-file "$PAF_PART"
    --aligner wfmash
    --map-pct-identity 90
    --min-aln-length 25k
    --num-mappings 1:1
    --scaffold-jump 0
    --overlap 0.95
    --scoring log-length-ani
    --threads "${SLURM_CPUS_PER_TASK:-8}"
    --temp-dir "$SCRATCH"
)
export PATH="$PINNED_DIR:$IMPG_DIR:$PATH" WFMASH_BIN_DIR="$PINNED_DIR"
if [[ ! -s $MAPPING/production.1to1.paf ]]; then
    printf '%q ' "${COMMAND[@]}" > "$MAPPING/command.txt"
    printf '\n' >> "$MAPPING/command.txt"
    printf '%s  %s\n' "$EXPECTED_SWEEPGA_SHA256" "$SWEEPGA_PINNED_REALPATH" > "$MAPPING/sweepga.sha256"
    "${COMMAND[@]}" > "$MAPPING/sweepga.stdout" 2> "$MAPPING/sweepga.stderr"
    [[ -s $PAF_PART ]] || { echo "SweepGA emitted an empty PAF for $DATASET" >&2; exit 65; }
    mv "$PAF_PART" "$MAPPING/production.1to1.paf"
else
    grep -F -- '--num-mappings 1:1' "$MAPPING/command.txt" >/dev/null
    grep -F "$EXPECTED_SWEEPGA_SHA256" "$MAPPING/sweepga.sha256" >/dev/null
fi

# Native read-only multiplicity audit.  This output is evidence only: it is
# compared with, and never promoted over, the production mapping.
RECHECK_PART="$MAPPING/native_recheck.1to1.paf.part"
"$SWEEPGA_PINNED_REALPATH" "$MAPPING/production.1to1.paf" \
    --output-file "$RECHECK_PART" --num-mappings 1:1 --scaffold-jump 0 \
    --overlap 0.95 --scoring log-length-ani --threads 2 \
    > "$MAPPING/native_recheck.stdout" 2> "$MAPPING/native_recheck.stderr"
[[ -s $RECHECK_PART ]] || { echo "Native 1:1 audit emitted an empty PAF for $DATASET" >&2; exit 65; }
mv "$RECHECK_PART" "$MAPPING/native_recheck.1to1.paf"

python3 "$ROOT/analysis/tier3a_origin_remap.py" audit-mapping \
    --paf "$MAPPING/production.1to1.paf" \
    --native-recheck-paf "$MAPPING/native_recheck.1to1.paf" \
    --h1-fai "$H1.fai" --h2-fai "$H2.fai" \
    --output "$MAPPING/native_multiplicity_audit.json" \
    --contig-map "$MAPPING/contig_coverage.tsv"

python3 "$ROOT/results/tier3a/acquisition_build_queries.py" \
    --dataset-id "$DATASET" --gff "$GFF" --fai "$H1.fai" --h2-fai "$H2.fai" \
    --bounded-paf "$MAPPING/production.1to1.paf" --output-dir "$QUERY" \
    --annotation-release "$ANNOTATION_RELEASE" > "$QUERY/build.json"

python3 "$ROOT/analysis/tier3a_origin_remap.py" stage-manifest \
    --base-manifest "$BASE_MANIFEST" --dataset-id "$DATASET" \
    --work-dir "$OUT" --output "$OUT/staging_manifest.tsv"

python3 "$ROOT/analysis/tier3a_biological.py" prepare \
    --manifest "$OUT/staging_manifest.tsv" --dataset-id "$DATASET" \
    --output-dir "$OUT" --panel-bases "$PANEL_BASES"

"$IMPG" index -a "$MAPPING/production.1to1.paf" -i "$OUT/graph.impg" \
    -t "${SLURM_CPUS_PER_TASK:-8}" > "$OUT/index.stdout" 2> "$OUT/index.stderr"
"$IMPG" partition -a "$MAPPING/production.1to1.paf" -i "$OUT/graph.impg" \
    -w 2000 -d 0 --min-missing-size 1 --min-boundary-distance 0 -o bed \
    --output-folder "$OUT/partitions" -t "${SLURM_CPUS_PER_TASK:-8}" -v 0 \
    > "$OUT/partition.stdout" 2> "$OUT/partition.stderr"

python3 "$ROOT/results/tier3a/acquisition_select_impg_partitions.py" \
    --partitions "$OUT/partitions/partitions.bed" --targets "$OUT/mapping_callable.bed" \
    --focus-bed "$OUT/focus.bed" --mapping-tsv "$OUT/partition_annotation_map.tsv"

"$IMPG" query -a "$MAPPING/production.1to1.paf" -i "$OUT/graph.impg" \
    -b "$OUT/focus.bed" -d 0 -o vcf:poa --force-large-region \
    --min-transitive-len 1 --sequence-files "$H1" "$H2" -O "$OUT/calls" \
    -t "${SLURM_CPUS_PER_TASK:-8}" > "$OUT/query.stdout" 2> "$OUT/query.stderr"
find "$OUT/calls" -name '*.vcf' -type f | sort > "$OUT/vcf.list"
[[ -s $OUT/vcf.list ]] || { echo "IMPG emitted no regional VCFs for $DATASET" >&2; exit 3; }

"$IMPG" lace -l "$OUT/vcf.list" --format vcf -o "$OUT/laced.vcf" \
    --reference "$H1" --compress none -t "${SLURM_CPUS_PER_TASK:-8}" \
    > "$OUT/lace.stdout" 2> "$OUT/lace.stderr"
bcftools norm -f "$H1" -m -any -Ob -o "$OUT/normalized.untrimmed.bcf" "$OUT/laced.vcf" \
    > "$OUT/norm.stdout" 2> "$OUT/norm.stderr"
bcftools index -f "$OUT/normalized.untrimmed.bcf"
bcftools view -R "$OUT/mapping_callable.bed" -Ou "$OUT/normalized.untrimmed.bcf" | \
    bcftools norm -d exact -Ob -o "$OUT/normalized.bcf"
bcftools index -f "$OUT/normalized.bcf"
bcftools view -Oz -o "$OUT/normalized.vcf.gz" "$OUT/normalized.bcf"
bcftools index -f -t "$OUT/normalized.vcf.gz"

python3 "$ROOT/analysis/tier3a_biological.py" summarize \
    --output-dir "$OUT" --focus-bed "$OUT/focus.bed" \
    --partitions "$OUT/partitions/partitions.bed" --index "$OUT/graph.impg" \
    --bcf "$OUT/normalized.bcf" --bcftools "$(command -v bcftools)" \
    --block-size 50000 --bootstrap-replicates 1000

END_EPOCH=$(date +%s)
export DATASET START_EPOCH END_EPOCH
python3 - "$OUT/success.json" <<'PY'
import json, os, sys
value = {
    "dataset_id": os.environ["DATASET"],
    "status": "completed",
    "slurm_job_id": os.environ.get("SLURM_ARRAY_JOB_ID", os.environ.get("SLURM_JOB_ID", "local")),
    "slurm_array_task_id": os.environ.get("SLURM_ARRAY_TASK_ID", "0"),
    "node": os.environ.get("SLURMD_NODENAME", os.environ.get("HOSTNAME", "unknown")),
    "requested_cpus": os.environ.get("SLURM_CPUS_PER_TASK", "8"),
    "requested_memory": os.environ.get("SLURM_MEM_PER_NODE", "65536M"),
    "start_epoch": int(os.environ["START_EPOCH"]),
    "end_epoch": int(os.environ["END_EPOCH"]),
    "elapsed_seconds": int(os.environ["END_EPOCH"]) - int(os.environ["START_EPOCH"]),
}
with open(sys.argv[1] + ".tmp", "w", encoding="utf-8") as handle:
    json.dump(value, handle, indent=2, sort_keys=True)
    handle.write("\n")
os.replace(sys.argv[1] + ".tmp", sys.argv[1])
PY
