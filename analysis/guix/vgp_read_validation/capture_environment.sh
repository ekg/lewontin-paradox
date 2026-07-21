#!/usr/bin/env bash
set -euo pipefail

REPOSITORY_ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/../../.." && pwd)
readonly REPOSITORY_ROOT
readonly DIR="$REPOSITORY_ROOT/analysis/guix/vgp_read_validation"
readonly CHANNELS="$REPOSITORY_ROOT/analysis/guix/vgp_10_pilot/channels.scm"
readonly MANIFEST="$DIR/manifest.scm"
readonly ROOT_CONFIG="$REPOSITORY_ROOT/analysis/vgp_data_root_config.json"
readonly OUTPUT=${1:?usage: capture_environment.sh OUTPUT-JSON}
readonly PROFILE_LINK=${VGP_READ_VALIDATION_PROFILE:-/moosefs/erikg/vgp/derived/read-validation/environment/profile}

VGP_ROOT=$(python3 - "$ROOT_CONFIG" <<'PY'
import json,sys
value=json.load(open(sys.argv[1]))["root"]
if value != "/moosefs/erikg/vgp": raise SystemExit("noncanonical VGP root")
print(value)
PY
)
readonly VGP_ROOT
readonly PROFILE=$(readlink -f "$PROFILE_LINK")
[[ $PROFILE == /gnu/store/*-profile && -d $PROFILE ]] || {
    echo "invalid realized profile: $PROFILE" >&2
    exit 2
}

tmp=$(mktemp -d)
trap 'rm -rf -- "$tmp"' EXIT
guix time-machine -C "$CHANNELS" -- describe -f channels > "$tmp/channels.resolved.scm"
guix package -p "$PROFILE_LINK" --list-installed > "$tmp/packages.tsv"
guix gc --requisites "$PROFILE" | LC_ALL=C sort -u > "$tmp/closure.txt"
derivation=$(guix gc --derivers "$PROFILE")
[[ $derivation == /gnu/store/*.drv ]] || { echo "profile derivation absent" >&2; exit 2; }

env VGP_ROOT="$VGP_ROOT" PROFILE="$PROFILE" PROFILE_LINK="$PROFILE_LINK" \
DERIVATION="$derivation" CHANNELS="$CHANNELS" MANIFEST="$MANIFEST" \
TMP_CAPTURE="$tmp" OUTPUT="$OUTPUT" python3 - <<'PY'
import hashlib,json,os,shutil
from pathlib import Path

def sha(path):
    value=hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda:handle.read(1024*1024),b""): value.update(block)
    return value.hexdigest()

profile=Path(os.environ["PROFILE"])
tmp=Path(os.environ["TMP_CAPTURE"])
executables=[]
for name in ("python3","pytest","minimap2","jellyfish","samtools","bcftools",
             "bgzip","bedtools","psmc","fq2psmcfa"):
    path=profile/"bin"/name
    if not path.is_file(): raise SystemExit(f"required executable absent: {path}")
    executables.append({"name":name,"path":str(path.resolve()),"sha256":sha(path)})
packages=[]
for line in (tmp/"packages.tsv").read_text().splitlines():
    name,version,output,path=line.split("\t")
    packages.append({"name":name.strip(),"version":version.strip(),
                     "output":output.strip(),"store_path":path.strip()})
result={
    "schema_version":"vgp-read-validation-environment-v1",
    "task_id":"validate-vgp-pilot-reads",
    "canonical_vgp_root":os.environ["VGP_ROOT"],
    "frozen_channel_commit":"44bbfc24e4bcc48d0e3343cd3d83452721af8c36",
    "channels_file":{"path":os.environ["CHANNELS"],"sha256":sha(os.environ["CHANNELS"])},
    "resolved_channels_sha256":sha(tmp/"channels.resolved.scm"),
    "manifest":{"path":os.environ["MANIFEST"],"sha256":sha(os.environ["MANIFEST"])},
    "profile_link":os.environ["PROFILE_LINK"],
    "profile":os.environ["PROFILE"],
    "profile_derivation":os.environ["DERIVATION"],
    "closure_sha256":sha(tmp/"closure.txt"),
    "closure_store_items":(tmp/"closure.txt").read_text().splitlines(),
    "packages":packages,
    "executables":executables,
    "realization_command":"guix time-machine -C analysis/guix/vgp_10_pilot/channels.scm -- package -L analysis/guix -p /moosefs/erikg/vgp/derived/read-validation/environment/profile -m analysis/guix/vgp_read_validation/manifest.scm --no-grafts",
    "capture_commands":[
        "guix time-machine -C analysis/guix/vgp_10_pilot/channels.scm -- describe -f channels",
        "guix package -p /moosefs/erikg/vgp/derived/read-validation/environment/profile --list-installed",
        "guix gc --requisites PROFILE | LC_ALL=C sort -u",
        "guix gc --derivers PROFILE",
    ],
}
output=Path(os.environ["OUTPUT"])
output.parent.mkdir(parents=True,exist_ok=True)
partial=output.with_name(output.name+".partial")
partial.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
partial.replace(output)
canonical=Path(os.environ["VGP_ROOT"])/"derived/read-validation/environment/environment_capture_v1.json"
canonical.parent.mkdir(parents=True,exist_ok=True)
canonical_partial=canonical.with_name(canonical.name+".partial")
shutil.copyfile(output,canonical_partial)
canonical_partial.replace(canonical)
PY
