#!/usr/bin/env bash
set -euo pipefail

: "${VGP_ROOT:?VGP_ROOT is required}"
: "${VGP_REPO_ROOT:?VGP_REPO_ROOT is required}"
: "${VGP_SCALE_ROOT:?VGP_SCALE_ROOT is required}"
: "${VGP_ENVIRONMENT_CAPTURE:?VGP_ENVIRONMENT_CAPTURE is required}"
pair=${1:?selection id is required}
[[ $VGP_ROOT == /moosefs/erikg/vgp ]] || { echo "ERROR: noncanonical VGP_ROOT" >&2; exit 2; }
[[ $VGP_SCALE_ROOT == "$VGP_ROOT"/derived/scale-vgp-real-v1 ]] || {
    echo "ERROR: scale root does not derive from VGP_ROOT" >&2; exit 2;
}
[[ $pair == P04 || $pair == P07 ]] || { echo "ERROR: unsupported pair: $pair" >&2; exit 2; }
[[ -d /scratch && -w /scratch ]] || { echo "ERROR: /scratch unavailable" >&2; exit 2; }

scratch=$(mktemp -d "/scratch/vgp-scale-${pair}-${SLURM_JOB_ID:?}-XXXXXX")
cleanup() {
    [[ $scratch == /scratch/vgp-scale-"$pair"-"$SLURM_JOB_ID"-* ]] || exit 2
    rm -rf -- "$scratch"
}
trap cleanup EXIT INT TERM
export TMPDIR=$scratch TMP=$scratch TEMP=$scratch PYTHONPATH=$VGP_REPO_ROOT
cd "$scratch"
scratch_resolved=$(readlink -f -- "$scratch")
cwd_resolved=$(readlink -f -- /proc/$$/cwd)
[[ $cwd_resolved == "$scratch_resolved" ]] || { echo "ERROR: cwd escaped private /scratch" >&2; exit 2; }
for variable in TMPDIR TMP TEMP; do
    value=${!variable}
    [[ $(readlink -f -- "$value") == "$scratch_resolved" ]] || { echo "ERROR: $variable escaped private /scratch" >&2; exit 2; }
done

partial="$scratch/$pair.json"
profile=$(python3 - "$VGP_ENVIRONMENT_CAPTURE" <<'PY'
import json,sys
value=json.load(open(sys.argv[1]))
profile=value.get("profile","")
if not profile.startswith("/gnu/store/"):
    raise SystemExit("uncaptured Guix profile")
print(profile)
PY
)
[[ -x $profile/bin/python3 ]] || { echo "ERROR: realized Guix Python unavailable" >&2; exit 2; }
unset PYTHONHOME CONDA_PREFIX VIRTUAL_ENV LD_LIBRARY_PATH
export PATH="$profile/bin:/usr/bin:/bin"
"$profile/bin/python3" -m analysis.vgp_10_pilot verify-capture "$VGP_ENVIRONMENT_CAPTURE" >/dev/null
"$profile/bin/python3" "$VGP_REPO_ROOT/analysis/scale_vgp_real.py" verify-pair \
    --pair "$pair" --vgp-root "$VGP_ROOT" --output "$partial" --job-id "$SLURM_JOB_ID"

python3 - "$partial" "$VGP_ROOT" "$pair" <<'PY'
import json,sys
value=json.load(open(sys.argv[1]))
assert value["canonical_vgp_root"] == sys.argv[2]
assert value["selection_id"] == sys.argv[3]
assert value["status"] == "PASS"
PY
mkdir -p "$VGP_SCALE_ROOT/verification"
final="$VGP_SCALE_ROOT/verification/$pair.json"
promote="$VGP_SCALE_ROOT/verification/.$pair.${SLURM_JOB_ID}.partial"
cp --reflink=auto -- "$partial" "$promote"
sync "$promote"
mv -f -- "$promote" "$final"
cmp -- "$partial" "$final"
