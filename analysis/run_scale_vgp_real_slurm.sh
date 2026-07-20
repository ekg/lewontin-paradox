#!/usr/bin/env bash
set -euo pipefail

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
mode=${1:---dry-run}
[[ $mode == --dry-run || $mode == --submit ]] || { echo "usage: $0 [--dry-run|--submit]" >&2; exit 2; }
VGP_ROOT=${VGP_ROOT:-/moosefs/erikg/vgp}
[[ $VGP_ROOT == /moosefs/erikg/vgp ]] || { echo "ERROR: noncanonical VGP_ROOT" >&2; exit 2; }
VGP_SCALE_ROOT=$VGP_ROOT/derived/scale-vgp-real-v1
VGP_ENVIRONMENT_CAPTURE=$ROOT/analysis/guix/vgp_10_pilot/realization.json
VGP_FASTGA_AMENDMENT=$ROOT/analysis/vgp_real_pilot_fastga_amendment_v1.json
[[ -f $VGP_ENVIRONMENT_CAPTURE ]] || { echo "ERROR: Guix realization capture absent" >&2; exit 2; }
[[ -f $VGP_FASTGA_AMENDMENT ]] || { echo "ERROR: FastGA amendment absent" >&2; exit 2; }
manifest=$VGP_SCALE_ROOT/submissions.tsv
mkdir -p "$VGP_SCALE_ROOT/logs"

python3 - "$VGP_SCALE_ROOT/inputs/P07.json" "$VGP_ROOT" <<'PY'
import hashlib,json,os,sys
from pathlib import Path
out=Path(sys.argv[1]); root=Path(sys.argv[2])
rows={
 "h1":("2b758e606304f7cb5e795d7939979b08c21bf4f3eac7ea3cf1c6ab0a463733c7",
       123263591,"438faaebe34180e563b700d911eb80973589f8ad4d5a70861747067621aaf6ba"),
 "h2":("6e7e9d88b88a3d030d80009191a06c525fbcc50044db05da4d81a8e9ad97ed40",
       122248201,"0bbd50ea5954e53e47e4bd80f6e01a22e0b80bfef50dd4781b4acaf5f5fb9418"),
}
inputs={}
for side,(digest,size,decompressed) in rows.items():
    path=root/"derived/scale-vgp-real-v1/cas/sha256"/digest[:2]/digest[2:4]/digest
    if not path.is_file() or path.stat().st_size != size:
        raise SystemExit(f"canonical {side} CAS object absent or wrong size: {path}")
    if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
        raise SystemExit(f"canonical {side} CAS digest mismatch")
    inputs[side]={"canonical_cas_path":str(path),"compressed_sha256":digest,
                  "compressed_bytes":size,"decompressed_sha256":decompressed}
value={"schema_version":"vgp-real-scaleout-input-binding-v1",
       "canonical_vgp_root":str(root),"selection_id":"P07",
       "migration_action":"verified hard-link reuse; no redownload","inputs":inputs}
out.parent.mkdir(parents=True,exist_ok=True)
partial=out.with_name(f".{out.name}.{os.getpid()}.partial")
partial.write_text(json.dumps(value,indent=2,sort_keys=True)+"\n")
os.replace(partial,out)
PY

if [[ $mode == --dry-run ]]; then
    for pair in P04 P07; do
        printf 'sbatch --parsable --job-name=vgp-scale-%s-verify --partition=workers --cpus-per-task=1 --mem=4G --time=04:00:00 --exclude=octopus11 ... verify_pair.sh %s\n' "$pair" "$pair"
    done
    printf 'sbatch --parsable --job-name=vgp-scale-P04-fastga --partition=workers --cpus-per-task=32 --mem=160G --time=3-00:00:00 --exclude=octopus11 ... independent_mapping.sh P04\n'
    printf 'sbatch --parsable --job-name=vgp-scale-P07-fastga --partition=workers --cpus-per-task=32 --mem=128G --time=3-00:00:00 --exclude=octopus11 ... revalidate_p07_fastga.sh\n'
    exit 0
fi

for pair in P04 P07; do
    output=$VGP_SCALE_ROOT/verification/$pair.json
    [[ ! -f $output ]] || { echo "RESUME: $pair already verified" >&2; continue; }
    job_id=$(sbatch --parsable --job-name="vgp-scale-$pair-verify" --partition=workers \
        --cpus-per-task=1 --mem=4G --time=04:00:00 --exclude=octopus11 \
        --output="$VGP_SCALE_ROOT/logs/$pair-%j.out" --error="$VGP_SCALE_ROOT/logs/$pair-%j.err" \
        --export="ALL,VGP_ROOT=$VGP_ROOT,VGP_REPO_ROOT=$ROOT,VGP_SCALE_ROOT=$VGP_SCALE_ROOT,VGP_ENVIRONMENT_CAPTURE=$VGP_ENVIRONMENT_CAPTURE" \
        "$ROOT/analysis/slurm/scale_vgp_real/verify_pair.sh" "$pair")
    python3 - "$manifest" "$job_id" "$pair" "$VGP_ROOT" "$ROOT" <<'PY'
import csv,datetime,hashlib,subprocess,sys
from pathlib import Path
path=Path(sys.argv[1]); exists=path.exists()
fields=("schema_version","canonical_vgp_root","authorization_id","wave_id","selection_id",
        "stage","job_id","submitted_at_utc","repository_commit","worker_sha256")
worker=Path(sys.argv[5])/"analysis/slurm/scale_vgp_real/verify_pair.sh"
with path.open("a",newline="",encoding="utf-8") as handle:
    writer=csv.DictWriter(handle,fieldnames=fields,delimiter="\t",lineterminator="\n")
    if not exists: writer.writeheader()
    writer.writerow({"schema_version":"vgp-real-scaleout-submission-v1",
      "canonical_vgp_root":sys.argv[4],"authorization_id":"vgp10-auth-20260718-v2",
      "wave_id":"SCALE-VERIFY-1","selection_id":sys.argv[3],"stage":"independent_verification",
      "job_id":sys.argv[2],"submitted_at_utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),
      "repository_commit":subprocess.check_output(
          ["git", "-C", sys.argv[5], "rev-parse", "HEAD"], text=True).strip(),
      "worker_sha256":hashlib.sha256(worker.read_bytes()).hexdigest()})
PY
    echo "$pair $job_id"
done

p04_manifests=("$VGP_ROOT"/pilot/independent/P04/mapping/P04.independent.*.native.1to1.paf.manifest.json)
if [[ -e ${p04_manifests[0]} ]]; then
    echo "RESUME: P04 FastGA /scratch contract already verified" >&2
elif awk -F '\t' '$5=="P04" && $6=="fastga_scratch_revalidation" {found=1} END {exit !found}' "$manifest"; then
    echo "RESUME: P04 FastGA /scratch revalidation already submitted" >&2
else
    worker=$ROOT/analysis/slurm/vgp_10_pilot/independent_mapping.sh
    job_id=$(sbatch --parsable --job-name=vgp-scale-P04-fastga --partition=workers \
        --cpus-per-task=32 --mem=160G --time=3-00:00:00 --exclude=octopus11 \
        --output="$VGP_SCALE_ROOT/logs/P04-fastga-%j.out" \
        --error="$VGP_SCALE_ROOT/logs/P04-fastga-%j.err" \
        --export="ALL,VGP_DATA_ROOT=$VGP_ROOT,VGP_ENVIRONMENT_CAPTURE=$VGP_ENVIRONMENT_CAPTURE,VGP_FASTGA_AMENDMENT=$VGP_FASTGA_AMENDMENT" \
        "$worker" P04)
    python3 - "$manifest" "$job_id" "$VGP_ROOT" "$ROOT" "$worker" <<'PY'
import csv,datetime,hashlib,subprocess,sys
from pathlib import Path
path=Path(sys.argv[1]); worker=Path(sys.argv[5])
fields=("schema_version","canonical_vgp_root","authorization_id","wave_id","selection_id",
        "stage","job_id","submitted_at_utc","repository_commit","worker_sha256")
with path.open("a",newline="",encoding="utf-8") as handle:
    writer=csv.DictWriter(handle,fieldnames=fields,delimiter="\t",lineterminator="\n")
    writer.writerow({"schema_version":"vgp-real-scaleout-submission-v1",
      "canonical_vgp_root":sys.argv[3],"authorization_id":"vgp10-auth-20260718-v2",
      "wave_id":"SCALE-FASTGA-REVALIDATE-1","selection_id":"P04",
      "stage":"fastga_scratch_revalidation","job_id":sys.argv[2],
      "submitted_at_utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),
      "repository_commit":subprocess.check_output(
          ["git", "-C", sys.argv[4], "rev-parse", "HEAD"], text=True).strip(),
      "worker_sha256":hashlib.sha256(worker.read_bytes()).hexdigest()})
PY
    echo "P04-fastga $job_id"
fi

fastga_output=$VGP_SCALE_ROOT/fastga/P07/contract.json
latest_p07_fastga_job=$(awk -F '\t' '$5=="P07" && $6=="fastga_scratch_revalidation" {job=$7} END {print job}' "$manifest")
active_p07_fastga_state=$(
    if [[ -n $latest_p07_fastga_job ]]; then squeue -h -j "$latest_p07_fastga_job" -o '%T'; fi
)
if [[ -f $fastga_output ]]; then
    echo "RESUME: P07 FastGA /scratch contract already verified" >&2
elif [[ $active_p07_fastga_state =~ ^(PENDING|RUNNING|CONFIGURING)$ ]]; then
    echo "RESUME: P07 FastGA /scratch revalidation job $latest_p07_fastga_job is $active_p07_fastga_state" >&2
else
    worker=$ROOT/analysis/slurm/scale_vgp_real/revalidate_p07_fastga.sh
    job_id=$(sbatch --parsable --job-name=vgp-scale-P07-fastga --partition=workers \
        --cpus-per-task=32 --mem=128G --time=3-00:00:00 --exclude=octopus11 \
        --output="$VGP_SCALE_ROOT/logs/P07-fastga-%j.out" \
        --error="$VGP_SCALE_ROOT/logs/P07-fastga-%j.err" \
        --export="ALL,VGP_ROOT=$VGP_ROOT,VGP_REPO_ROOT=$ROOT,VGP_SCALE_ROOT=$VGP_SCALE_ROOT,VGP_ENVIRONMENT_CAPTURE=$VGP_ENVIRONMENT_CAPTURE,VGP_FASTGA_AMENDMENT=$VGP_FASTGA_AMENDMENT" \
        "$worker")
    python3 - "$manifest" "$job_id" "$VGP_ROOT" "$ROOT" "$worker" <<'PY'
import csv,datetime,hashlib,subprocess,sys
from pathlib import Path
path=Path(sys.argv[1]); exists=path.exists()
fields=("schema_version","canonical_vgp_root","authorization_id","wave_id","selection_id",
        "stage","job_id","submitted_at_utc","repository_commit","worker_sha256")
worker=Path(sys.argv[5])
with path.open("a",newline="",encoding="utf-8") as handle:
    writer=csv.DictWriter(handle,fieldnames=fields,delimiter="\t",lineterminator="\n")
    if not exists: writer.writeheader()
    writer.writerow({"schema_version":"vgp-real-scaleout-submission-v1",
      "canonical_vgp_root":sys.argv[3],"authorization_id":"vgp10-auth-20260718-v2",
      "wave_id":"SCALE-FASTGA-REVALIDATE-1","selection_id":"P07",
      "stage":"fastga_scratch_revalidation","job_id":sys.argv[2],
      "submitted_at_utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),
      "repository_commit":subprocess.check_output(
          ["git", "-C", sys.argv[4], "rev-parse", "HEAD"], text=True).strip(),
      "worker_sha256":hashlib.sha256(worker.read_bytes()).hexdigest()})
PY
    echo "P07-fastga $job_id"
fi
