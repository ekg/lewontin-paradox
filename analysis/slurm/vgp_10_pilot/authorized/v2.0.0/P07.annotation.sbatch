#!/usr/bin/env bash
#SBATCH --job-name=vgp10-P07-annotation
#SBATCH --partition=highmem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=1-00:00:00
#SBATCH --output=/moosefs/erikg/vgp/pilot/logs/vgp10-auth-v2-P07-annotation-%j.out
#SBATCH --error=/moosefs/erikg/vgp/pilot/logs/vgp10-auth-v2-P07-annotation-%j.err
set -euo pipefail

SUBMIT_ROOT="${SLURM_SUBMIT_DIR:?submit from the repository root}"
CANONICAL_ROOT=/moosefs/erikg/vgp
RUN_ROOT="$CANONICAL_ROOT/pilot/authorized-checkpoints/vgp10-auth-20260718-v2/P07/run"
CAPTURE="$SUBMIT_ROOT/analysis/guix/vgp_10_pilot/realization.json"
H1_CAS="$CANONICAL_ROOT/objects/sha256/2b/75/2b758e606304f7cb5e795d7939979b08c21bf4f3eac7ea3cf1c6ab0a463733c7"
GFF_CAS="$CANONICAL_ROOT/objects/sha256/8f/64/8f640543accd8081d1b7048eda32c9f1eef33b02f321b7b0f8adcf3b01dd6838"
OUTPUT="$CANONICAL_ROOT/pilot/outputs/vgp10-auth-20260718-v2/P07/annotation/exact_partitions.json"
for required in "$CAPTURE" "$H1_CAS" "$GFF_CAS" \
    "$RUN_ROOT/consensus/masks/callable.bed" "$RUN_ROOT/consensus/normalized.vcf"; do
    [[ -f $required ]] || { echo "required annotation input absent: $required" >&2; exit 2; }
done

scratch=$(mktemp -d -- "/tmp/vgp-P07-annotation-${SLURM_JOB_ID:?}-XXXXXX")
cleanup() {
    [[ $scratch == "/tmp/vgp-P07-annotation-${SLURM_JOB_ID}-"* ]] || return 2
    rm -rf -- "$scratch"
}
trap cleanup EXIT
cp "$H1_CAS" "$scratch/h1.fa.gz"
cp "$GFF_CAS" "$scratch/annotation.gff.gz"
cp "$RUN_ROOT/consensus/masks/callable.bed" "$scratch/callable.bed"
cp "$RUN_ROOT/consensus/normalized.vcf" "$scratch/normalized.vcf"
printf '%s  %s\n' \
    2b758e606304f7cb5e795d7939979b08c21bf4f3eac7ea3cf1c6ab0a463733c7 "$scratch/h1.fa.gz" \
    8f640543accd8081d1b7048eda32c9f1eef33b02f321b7b0f8adcf3b01dd6838 "$scratch/annotation.gff.gz" \
    | sha256sum -c -

export VGP_DATA_ROOT="$CANONICAL_ROOT" VGP_RUN_ID=vgp10-auth-20260718-v2-canary
export VGP_ENVIRONMENT_CAPTURE="$CAPTURE" PYTHONPATH="$SUBMIT_ROOT"
source "$SUBMIT_ROOT/analysis/slurm/vgp_10_pilot/common.sh"
require_runtime
python3 "$SUBMIT_ROOT/analysis/vgp_real_canary_annotation.py" \
    --h1-fasta "$scratch/h1.fa.gz" \
    --annotation-gff "$scratch/annotation.gff.gz" \
    --annotation-source-path "$GFF_CAS" \
    --callable-bed "$scratch/callable.bed" \
    --normalized-vcf "$scratch/normalized.vcf" \
    --canonical-root "$CANONICAL_ROOT" \
    --output "$OUTPUT"
