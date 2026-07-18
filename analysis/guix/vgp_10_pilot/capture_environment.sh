#!/usr/bin/env bash
set -euo pipefail

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/../../.." && pwd)
DIR="$ROOT/analysis/guix/vgp_10_pilot"
CHANNELS="$DIR/channels.scm"
MANIFEST="$DIR/manifest.scm"
LOCK="$DIR/environment-lock.json"
OUTPUT=${1:?usage: capture_environment.sh OUTPUT-JSON}

[[ $OUTPUT = /* ]] || OUTPUT="$PWD/$OUTPUT"
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

# No implicit update, module, container, or network installer is permitted.
# Guix may use its configured substitutes/source cache while realizing the
# exact authenticated channel and source hashes.
guix time-machine -C "$CHANNELS" -- describe -f channels >"$tmp/channels.resolved.scm"
guix time-machine -C "$CHANNELS" -- package -L "$ROOT/analysis/guix" \
    -p "$tmp/profile" -m "$MANIFEST" --no-grafts
profile=$(readlink -f "$tmp/profile")
derivation=$(guix gc --derivers "$profile")
[[ $derivation == /gnu/store/*.drv ]] || {
    echo "failed to resolve the profile derivation" >&2
    exit 2
}
guix gc --requisites "$profile" | LC_ALL=C sort -u >"$tmp/closure.txt"
psmc_output=$(guix time-machine -C "$CHANNELS" -- build -L "$ROOT/analysis/guix" \
    -e '(@ (vgp_10_pilot packages psmc) psmc-vgp-pinned)' --no-grafts --check)
psmc_derivation=$(guix gc --derivers "$psmc_output")

PROFILE="$profile" DERIVATION="$derivation" PSMC_OUTPUT="$psmc_output" \
PSMC_DERIVATION="$psmc_derivation" OUTPUT="$OUTPUT" TMP_CAPTURE="$tmp" \
LOCK="$LOCK" python3 - <<'PY'
import hashlib
import json
import os
from pathlib import Path

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

profile = Path(os.environ["PROFILE"])
tmp = Path(os.environ["TMP_CAPTURE"])
executables = []
external = {
    "sweepga": os.environ.get("VGP_SWEEPGA_BIN", ""),
    "impg": os.environ.get("VGP_IMPG_BIN", ""),
    "gfaffix": os.environ.get("VGP_GFAFFIX_BIN", ""),
}
impg_tools = Path(external["impg"]).parent if external["impg"] else Path("/")
for name in ("ALNtoPAF", "FAtoGDB", "FastGA", "GIXmake", "GIXrm", "ONEview", "PAFtoALN", "wfmash"):
    external[name] = str(impg_tools / name)
for name in ("sweepga", "impg", "gfaffix", "ALNtoPAF", "FAtoGDB", "FastGA", "GIXmake",
             "GIXrm", "ONEview", "PAFtoALN", "wfmash", "bcftools", "samtools", "psmc",
             "fq2psmcfa", "splitfa"):
    path = Path(external[name]) if name in external and external[name] else profile / "bin" / name
    if not path.is_file():
        raise SystemExit(f"required executable absent (set VGP_{name.upper()}_BIN for companion builds): {name}")
    executables.append({"name": name, "path": str(path.resolve()), "sha256": sha(path)})

lock = json.loads(Path(os.environ["LOCK"]).read_text())
expected = {
    "sweepga": lock["source_identities"]["sweepga"]["accepted_executable_sha256"],
    "impg": lock["source_identities"]["impg"]["accepted_executable_sha256"],
    "gfaffix": lock["source_identities"]["impg"]["gfaffix_executable_sha256"],
}
expected.update(lock["source_identities"]["impg"]["companion_executable_sha256"])
for row in executables:
    if row["name"] in expected and row["sha256"] != expected[row["name"]]:
        raise SystemExit(f"accepted executable mismatch: {row['name']}")
result = {
    "schema_version": lock["schema_version"],
    "channel_commit": lock["channel_commit"],
    "resolved_channels_sha256": sha(tmp / "channels.resolved.scm"),
    "manifest_sha256": sha(Path(os.environ["LOCK"]).parent / "manifest.scm"),
    "profile": str(profile),
    "derivation": os.environ["DERIVATION"],
    "closure_sha256": sha(tmp / "closure.txt"),
    "closure_store_items": (tmp / "closure.txt").read_text().splitlines(),
    "executables": executables,
    "source_identities": lock["source_identities"],
    "reproducibility": {
        "psmc_guix_check": "passed",
        "psmc_output": os.environ["PSMC_OUTPUT"],
        "psmc_derivation": os.environ["PSMC_DERIVATION"],
        "method": "guix build --check of the identical source-pinned derivation"
    },
}
output = Path(os.environ["OUTPUT"])
output.parent.mkdir(parents=True, exist_ok=True)
partial = output.with_name(output.name + ".partial")
partial.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
partial.replace(output)
PY
