#!/usr/bin/env bash
set -euo pipefail

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
SELECTION=${1:-$ROOT/analysis/vgp_clean_canary_selection_v1.json}
CAPTURE="$ROOT/analysis/guix/vgp_10_pilot/realization.json"
PROFILE=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["profile"])' "$CAPTURE")
BGZIP="$PROFILE/bin/bgzip"
SAMTOOLS=$(python3 -c 'import json,sys; print(next(x["path"] for x in json.load(open(sys.argv[1]))["executables"] if x["name"]=="samtools"))' "$CAPTURE")
[[ -x $BGZIP && -x $SAMTOOLS ]] || { echo "pinned BGZF tools unavailable" >&2; exit 2; }
output_dir=/moosefs/erikg/vgp/derived/clean-canary-bgzf
mkdir -p "$output_dir"
build_dir=$(mktemp -d -- "$output_dir/.build.XXXXXX")
cleanup() { [[ $build_dir == "$output_dir/.build."* ]] && rm -rf -- "$build_dir"; }
trap cleanup EXIT

for side in h1 h2; do
    readarray -t fields < <(python3 - "$SELECTION" "$side" <<'PY'
import json,sys
row=json.load(open(sys.argv[1]))["immutable_bgzf_inputs"][sys.argv[2]]
print(row["accession_version"]); print(row["source_cas_path"]); print(row["source_cas_sha256"]); print(row["derived_bgzf_path"])
PY
)
    accession=${fields[0]}; source=${fields[1]}; expected=${fields[2]}; target=${fields[3]}
    printf '%s  %s\n' "$expected" "$source" | sha256sum -c -
    [[ $target == "$output_dir/$accession.fa.gz" ]] || { echo "unsafe BGZF target" >&2; exit 2; }
    if [[ ! -f $target ]]; then
        gzip -cd -- "$source" | "$BGZIP" -@ 8 -c >"$build_dir/$accession.fa.gz"
        "$SAMTOOLS" faidx "$build_dir/$accession.fa.gz"
        "$BGZIP" -r "$build_dir/$accession.fa.gz"
        mv "$build_dir/$accession.fa.gz" "$target"
        mv "$build_dir/$accession.fa.gz.fai" "$target.fai"
        mv "$build_dir/$accession.fa.gz.gzi" "$target.gzi"
    fi
done

python3 - "$SELECTION" "$output_dir/manifest.json" "$BGZIP" "$SAMTOOLS" <<'PY'
import datetime,hashlib,json,sys
from pathlib import Path
selection=json.load(open(sys.argv[1])); assets={}
def digest(path):
 h=hashlib.sha256()
 with path.open("rb") as handle:
  for block in iter(lambda:handle.read(16*1024*1024),b""): h.update(block)
 return h.hexdigest()
for side,row in selection["immutable_bgzf_inputs"].items():
 path=Path(row["derived_bgzf_path"])
 assets[side]={"accession_version":row["accession_version"],"path":str(path),
  "sha256":digest(path),"bytes":path.stat().st_size,"fai_sha256":digest(Path(str(path)+".fai")),
  "gzi_sha256":digest(Path(str(path)+".gzi")),"source_cas_sha256":row["source_cas_sha256"]}
value={"schema_version":"vgp-clean-canary-bgzf-source-v1","task_id":"run-vgp-clean-canary",
 "created_at_utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),"selection_id":"P07",
 "bgzip":sys.argv[3],"samtools":sys.argv[4],"assets":assets}
Path(sys.argv[2]).write_text(json.dumps(value,sort_keys=True)+"\n")
PY
sha256sum "$output_dir/manifest.json"
