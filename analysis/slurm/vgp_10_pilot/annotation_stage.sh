#!/usr/bin/env bash
#SBATCH --job-name=vgp10-annotation
set -euo pipefail

source "$(dirname -- "$0")/common.sh"
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
python3 - "$input_dir/input-manifest.json" "$VGP_PAIR_RUN/consensus/masks/callable.bed" \
    "$VGP_STAGE_PARTIAL/normalized.vcf" "$VGP_STAGE_PARTIAL" <<'PY'
import csv,json,sys
from pathlib import Path
from analysis.vgp_10_pilot import (
    canonical_json,parse_fasta,parse_vcf,read_bed,sequence_dictionary,
    sha256_file,summarize_annotation_partitions,validate_annotation_binding,
)
manifest=json.load(open(sys.argv[1]))
out=Path(sys.argv[4])
dictionary=sequence_dictionary(parse_fasta(manifest["assets"]["h1_fasta"]["path"]))
annotation=manifest.get("annotation")
audit=validate_annotation_binding(manifest["h1_accession_version"],dictionary,annotation)
(out/"annotation_dictionary_audit.json").write_text(canonical_json(audit))
if audit["annotation_outputs_allowed"]:
    gff=Path(annotation["gff_path"])
    if not gff.is_file() or sha256_file(gff) != annotation["gff_sha256"]:
        raise SystemExit("exact annotation GFF digest mismatch")
    feature_paths=annotation.get("feature_beds",{})
    features={}
    for name,record in feature_paths.items():
        path=Path(record["path"])
        if not path.is_file() or sha256_file(path) != record["sha256"]:
            raise SystemExit(f"annotation feature BED digest mismatch: {name}")
        features[name]=read_bed(path)
    rows=summarize_annotation_partitions(read_bed(sys.argv[2]),parse_vcf(sys.argv[3]),features)
    with (out/"feature_diversity.tsv").open("w",newline="") as handle:
        writer=csv.DictWriter(handle,fieldnames=("partition","callable_h1_bp","variant_records",
                              "diversity_per_callable_h1_bp"),delimiter="\t",lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
else:
    (out/"annotation_not_available.txt").write_text(
        "Annotation absence does not block genome-wide diversity or PSMC.\n")
PY
promote_stage
