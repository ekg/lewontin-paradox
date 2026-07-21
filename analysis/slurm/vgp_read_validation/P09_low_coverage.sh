#!/usr/bin/env bash
#SBATCH --job-name=vgp-readval-P09
#SBATCH --partition=workers
#SBATCH --cpus-per-task=32
#SBATCH --mem=160G
#SBATCH --time=1-00:00:00
#SBATCH --output=/moosefs/erikg/vgp/logs/vgp-read-validation-P09-%j.out
#SBATCH --error=/moosefs/erikg/vgp/logs/vgp-read-validation-P09-%j.err
set -euo pipefail

readonly TASK_ID=validate-vgp-pilot-reads
readonly SELECTION_ID=P09
readonly PROFILE=/moosefs/erikg/vgp/derived/read-validation/environment/profile
readonly NODE_LOCAL_BASE=/scratch
readonly REQUIRED_SCRATCH_BYTES=160000000000
if [[ -n ${SLURM_JOB_ID:-} ]]; then
    readonly REPOSITORY_ROOT=${SLURM_SUBMIT_DIR:?submit from the repository root}
else
    REPOSITORY_ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/../../.." && pwd)
    readonly REPOSITORY_ROOT
fi
cd "$REPOSITORY_ROOT"

readonly VGP_ROOT=$(python3 - <<'PY'
import json
value=json.load(open("analysis/vgp_data_root_config.json"))["root"]
if value != "/moosefs/erikg/vgp": raise SystemExit("noncanonical VGP root")
print(value)
PY
)
[[ -x $PROFILE/bin/python3 ]] || { echo "missing pinned profile" >&2; exit 2; }
export GUIX_PROFILE="$PROFILE"
# Guix's generated profile uses `${GUIX_PYTHONPATH:+...}`; initialize it
# before sourcing under `set -u` on clean Slurm nodes.
export GUIX_PYTHONPATH=${GUIX_PYTHONPATH:-}
# shellcheck disable=SC1091
source "$GUIX_PROFILE/etc/profile"
export PYTHONPATH="$REPOSITORY_ROOT"
export LC_ALL=C LANG=C
[[ -n ${SLURM_JOB_ID:-} && -d $NODE_LOCAL_BASE && -w $NODE_LOCAL_BASE ]] || {
    echo "biological computation requires Slurm node-local scratch" >&2; exit 2;
}
case $(stat -f -c %T -- "$NODE_LOCAL_BASE") in
    nfs|nfs4|fuse*|lustre|gpfs|ceph) echo "scratch is not node-local" >&2; exit 2 ;;
esac
available_scratch=$(df -PB1 -- "$NODE_LOCAL_BASE" | awk 'NR==2 {print $4}')
(( available_scratch >= REQUIRED_SCRATCH_BYTES )) || {
    echo "insufficient node-local scratch: $available_scratch < $REQUIRED_SCRATCH_BYTES" >&2
    exit 2
}
for tool in python3 minimap2 sha256sum time; do command -v "$tool" >/dev/null; done

readonly STARTED_UTC=$(date -u +%Y-%m-%dT%H:%M:%SZ)
readonly WORK=$(mktemp -d "$NODE_LOCAL_BASE/vgp-readval-P09-${SLURM_JOB_ID}-XXXXXX")
readonly OUT="$WORK/output"
mkdir -p "$OUT/telemetry"
trap 'rm -rf -- "$WORK"' EXIT
exec > >(tee "$OUT/telemetry/worker.stdout") 2> >(tee "$OUT/telemetry/worker.stderr" >&2)
cp -- "$0" "$OUT/executed_worker.sh"
cp -- analysis/vgp_read_validation.py "$OUT/executed_validation_module.py"
cp -- analysis/vgp_data_root_config.json "$OUT/root_config.json"
run_timed() {
    local label=$1; shift
    echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] START $label: $*"
    command time -v -o "$OUT/telemetry/${label}.time.txt" -- "$@"
    echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] DONE $label"
}

readonly READ_MANIFEST="$VGP_ROOT/pilot/manifests/vgp_validation_reads_manifest_v1.json"
readonly READ_VIEW="$VGP_ROOT/views/accession/P09/SRR29944135/SRR29944135_subreads.fastq.gz"
readonly INPUT_MANIFEST="$VGP_ROOT/pilot/inputs/P09/input-manifest.json"
readonly H1="$VGP_ROOT/pilot/inputs/P09/h1.fa"
for path in "$READ_MANIFEST" "$READ_VIEW" "$INPUT_MANIFEST" "$H1"; do
    [[ -r $path ]] || { echo "missing P09 input: $path" >&2; exit 2; }
done
run_timed stage_reads cp --reflink=auto -- "$READ_VIEW" "$WORK/reads.fastq.gz"
run_timed stage_h1 cp --reflink=auto -- "$H1" "$WORK/h1.fa"

python3 - "$READ_MANIFEST" "$INPUT_MANIFEST" "$WORK" "$OUT/input_manifest.json" <<'PY'
import hashlib,json,sys
from pathlib import Path
def digest(path):
    value=hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda:handle.read(8*1024*1024),b""): value.update(block)
    return value.hexdigest()
read_manifest,input_manifest,work,out=map(Path,sys.argv[1:])
reads=json.loads(read_manifest.read_text())
row=next(row for row in reads["objects"] if row["selection_id"]=="P09")
if reads["canonical_root"]!="/moosefs/erikg/vgp" or row["status"]!="verified":
    raise SystemExit("P09 read manifest is not canonical and verified")
staged=work/"reads.fastq.gz"
read_digest=digest(staged)
if read_digest != row["local_sha256"] or staged.stat().st_size != row["expected_bytes"]:
    raise SystemExit("P09 staged read mismatch")
pair=json.loads(input_manifest.read_text())
h1=work/"h1.fa"
h1_digest=digest(h1)
expected=pair["assets"]["h1_fasta"]
if pair["canonical_vgp_root"]!="/moosefs/erikg/vgp" or h1_digest!=expected["sha256"]:
    raise SystemExit("P09 H1 identity mismatch")
result={"schema_version":"vgp-read-validation-input-manifest-v1","task_id":"validate-vgp-pilot-reads",
        "selection_id":"P09","canonical_vgp_root":"/moosefs/erikg/vgp","exact_individual":"sHetFra1",
        "biosample":"SAMN39432692","inputs":[
            {"role":"one_complete_HiFi_SMRT_cell_of_seven","source":row["accession_view_path"],
             "bytes":staged.stat().st_size,"sha256":read_digest,"raw_bases_metadata":row["base_count"]},
            {"role":"H1_reference","source":expected["path"],"bytes":h1.stat().st_size,
             "sha256":h1_digest,"accession_version":"GCA_036365525.1"}]}
out.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
PY

run_timed map_hifi minimap2 -t 32 -I 8G -x map-hifi --secondary=no -c \
    -o "$OUT/hifi_to_H1.primary.paf" "$WORK/h1.fa" "$WORK/reads.fastq.gz"
python3 -m analysis.vgp_read_validation paf-summary \
    --paf "$OUT/hifi_to_H1.primary.paf" --reference-bp 6013065386 --raw-bases 9742646618 \
    --raw-reads 495557 \
    --selection-id P09 --output "$OUT/mapping_summary.json"
python3 - "$OUT/mapping_summary.json" "$OUT/validation_limits.json" "$STARTED_UTC" "$SLURM_JOB_ID" <<'PY'
import json,sys,time
mapping_path,out,started,job=sys.argv[1:]
mapping=json.load(open(mapping_path))
coverage=mapping["nominal_physical_coverage"]
result={"schema_version":"vgp-read-validation-limits-v1","task_id":"validate-vgp-pilot-reads",
        "selection_id":"P09","canonical_vgp_root":"/moosefs/erikg/vgp","slurm_job_id":job,
        "started_utc":started,"completed_utc":time.strftime("%Y-%m-%dT%H:%M:%SZ",time.gmtime()),
        "raw_scope":"one complete HiFi SMRT cell from a seven-run assembly library",
        "nominal_coverage":coverage,"assembly_qv_status":"not_estimable",
        "kmer_heterozygosity_status":"not_estimable","read_pi_status":"not_estimable",
        "psmc_validation_status":"not_estimable",
        "reasons":["nominal coverage is below 2x", "one SMRT cell is not the complete assembly read library",
                   "diploid reads mapped to one phased haplotype confound heterozygosity with sequence error",
                   "the inherited P09 assembly-derived pi/PSMC product was incomplete at validation freeze"],
        "effect":"assigns low raw-read validation confidence without deleting or converting any future core result to zero"}
out=open(out,"w"); json.dump(result,out,indent=2,sort_keys=True); out.write("\n")
PY

cp -- "$OUT/telemetry/worker.stdout" "$OUT/telemetry/worker.stdout.snapshot"
cp -- "$OUT/telemetry/worker.stderr" "$OUT/telemetry/worker.stderr.snapshot"
find "$OUT" -type f ! -name output_manifest.tsv ! -name worker.stdout ! -name worker.stderr -print0 | sort -z | \
    xargs -0 sha256sum | sed "s#  $OUT/#\t#" | awk 'BEGIN{OFS="\t";print "sha256","relative_path"}{print $1,$2}' \
    > "$OUT/output_manifest.tsv"
readonly PARTIAL="$VGP_ROOT/staging/outputs/vgp-read-validation-P09-${SLURM_JOB_ID}.partial"
readonly FINAL="$VGP_ROOT/derived/read-validation/runs/P09/slurm-${SLURM_JOB_ID}"
[[ ! -e $PARTIAL && ! -e $FINAL ]] || { echo "refusing overwrite" >&2; exit 2; }
mkdir -p "${PARTIAL%/*}" "${FINAL%/*}"
cp -a -- "$OUT" "$PARTIAL"
(cd "$PARTIAL" && tail -n +2 output_manifest.tsv | while IFS=$'\t' read -r digest relative; do
    printf '%s  %s\n' "$digest" "$relative"
done | sha256sum -c -)
mv -- "$PARTIAL" "$FINAL"
echo "PROMOTED=$FINAL"
