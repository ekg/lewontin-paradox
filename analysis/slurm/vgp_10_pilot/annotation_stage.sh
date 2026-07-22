#!/usr/bin/env bash
#SBATCH --job-name=vgp10-annotation
set -euo pipefail

if [[ -n ${SLURM_JOB_ID:-} ]]; then
    VGP_STAGE_REPO_ROOT=${SLURM_SUBMIT_DIR:?submit from the repository root}
else
    VGP_STAGE_REPO_ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/../../.." && pwd)
fi
source "$VGP_STAGE_REPO_ROOT/analysis/slurm/vgp_10_pilot/common.sh"
VGP_SELECTION_ID=${1:?usage: annotation_stage.sh SELECTION_ID}
VGP_STAGE_NAME=annotation
export VGP_SELECTION_ID VGP_STAGE_NAME
require_selection_id "$VGP_SELECTION_ID"
require_runtime
begin_stage "$VGP_STAGE_NAME"

input_dir="$VGP_DATA_ROOT/pilot/inputs/$VGP_SELECTION_ID"
[[ -f $VGP_PAIR_RUN/consensus/.complete.json ]] || fail "core consensus sentinel absent"
verify_tool bcftools
bcftools=$(tool_path bcftools)
"$bcftools" view -Ov -o "$VGP_STAGE_PARTIAL/normalized.vcf" \
    "$VGP_PAIR_RUN/variants/normalized.vcf.gz"
python3 - "$input_dir/input-manifest.json" "$VGP_STAGE_PARTIAL" <<'PY'
import json,sys
from pathlib import Path
from analysis.vgp_10_pilot import (
    canonical_json,parse_fasta,sequence_dictionary,sha256_file,validate_annotation_binding,
)
manifest=json.load(open(sys.argv[1]))
out=Path(sys.argv[2])
dictionary=sequence_dictionary(parse_fasta(manifest["assets"]["h1_fasta"]["path"]))
annotation=manifest.get("annotation")
audit=validate_annotation_binding(manifest["h1_accession_version"],dictionary,annotation)
(out/"annotation_dictionary_audit.json").write_text(canonical_json(audit))
if not audit["annotation_outputs_allowed"]:
    raise SystemExit("annotation stage was submitted without an eligible exact binding")
gff=Path(annotation["gff_path"])
if not gff.is_file() or sha256_file(gff) != annotation["gff_sha256"]:
    raise SystemExit("exact annotation GFF digest mismatch")
PY
readarray -t annotation_fields < <(python3 - "$input_dir/input-manifest.json" <<'PY'
import json,sys
value=json.load(open(sys.argv[1])); annotation=value["annotation"]
print(value["assets"]["h1_fasta"]["path"])
print(annotation["gff_path"])
print(value["h1_accession_version"])
print(annotation["annotation_accession_version"])
PY
)
cp "${annotation_fields[1]}" "$VGP_STAGE_PARTIAL/annotation.gff.gz"
python3 "$VGP_REPO_ROOT/analysis/vgp_real_canary_annotation.py" \
    --h1-fasta "${annotation_fields[0]}" \
    --annotation-gff "$VGP_STAGE_PARTIAL/annotation.gff.gz" \
    --annotation-source-path "${annotation_fields[1]}" \
    --callable-bed "$VGP_PAIR_RUN/consensus/masks/callable.bed" \
    --normalized-vcf "$VGP_STAGE_PARTIAL/normalized.vcf" \
    --canonical-root "$VGP_DATA_ROOT" \
    --selection-id "$VGP_SELECTION_ID" \
    --assembly-accession-version "${annotation_fields[2]}" \
    --annotation-accession-version "${annotation_fields[3]}" \
    --task-id "${VGP_TASK_ID:-run-vgp-real-pilot}" \
    --schema-version "${VGP_ANNOTATION_SCHEMA_VERSION:-vgp-real-pilot-exact-annotation-v1}" \
    --output "$VGP_STAGE_PARTIAL/exact_partitions.json"
promote_stage
