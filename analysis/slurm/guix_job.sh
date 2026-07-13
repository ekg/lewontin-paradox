#!/usr/bin/env bash
set -euo pipefail

# The profile accepted below is realized and rooted only by prepare_guix.sh's
# committed `guix time-machine -C ... package -m ...` invocation.  Compute
# nodes cannot contact the daemon, so they verify that record and execute the
# identical shared-store profile rather than attempting a compute-side build.
SCRIPT_ROOT=$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)
TIER3_ROOT=${TIER3_ROOT:-$SCRIPT_ROOT}
export TIER3_ROOT

if (( $# < 2 )); then
    echo "usage: guix_job.sh ENVIRONMENT_RECORD COMMAND [ARG ...]" >&2
    exit 64
fi

ENVIRONMENT_RECORD=$1
shift

case "$ENVIRONMENT_RECORD" in
    /*) ;;
    *) echo "environment record must be an absolute path" >&2; exit 64 ;;
esac

readarray -t PROFILE_FIELDS < <(
    /usr/bin/python3 - "$ENVIRONMENT_RECORD" <<'PY'
import json
import sys
from pathlib import Path

record_path = Path(sys.argv[1])
record = json.loads(record_path.read_text(encoding="utf-8"))
for key in ("profile_store_path", "channels_file", "channels_sha256", "manifest_file", "manifest_sha256"):
    value = record.get(key)
    if not isinstance(value, str) or not value or "\n" in value:
        raise SystemExit("invalid environment record field: " + key)
    print(value)
PY
)

if (( ${#PROFILE_FIELDS[@]} != 5 )); then
    echo "environment record did not yield five provenance fields" >&2
    exit 65
fi
PROFILE=${PROFILE_FIELDS[0]}
CHANNELS=${PROFILE_FIELDS[1]}
CHANNELS_SHA256=${PROFILE_FIELDS[2]}
MANIFEST=${PROFILE_FIELDS[3]}
MANIFEST_SHA256=${PROFILE_FIELDS[4]}

case "$CHANNELS" in
    /*) ;;
    *) CHANNELS=${TIER3_ROOT:?TIER3_ROOT is required}/$CHANNELS ;;
esac
case "$MANIFEST" in
    /*) ;;
    *) MANIFEST=${TIER3_ROOT:?TIER3_ROOT is required}/$MANIFEST ;;
esac

case "$PROFILE" in
    /gnu/store/????????????????????????????????-*) ;;
    *) echo "recorded profile is not a Guix store path: $PROFILE" >&2; exit 65 ;;
esac
[[ -d "$PROFILE" ]] || { echo "shared Guix profile is not visible: $PROFILE" >&2; exit 69; }
[[ $(sha256sum "$CHANNELS" | awk '{print $1}') == "$CHANNELS_SHA256" ]] || {
    echo "committed channels lock differs from realization record" >&2
    exit 65
}
[[ $(sha256sum "$MANIFEST" | awk '{print $1}') == "$MANIFEST_SHA256" ]] || {
    echo "committed manifest differs from realization record" >&2
    exit 65
}

PYTHON_SITE=$(find "$PROFILE/lib" -maxdepth 2 -type d -path '*/python*/site-packages' -print -quit)
[[ -n "$PYTHON_SITE" ]] || { echo "realized profile has no Python site-packages" >&2; exit 69; }

JOB_KEY=${SLURM_JOB_ID:-login}-$(id -u)
BASE_TMP=${TIER3_SCRATCH_ROOT:-${SCRATCH:-/tmp}}
JOB_HOME=$BASE_TMP/tier3-job-home/$JOB_KEY
JOB_TMP=$BASE_TMP/tier3-job-tmp/$JOB_KEY
install -d -m 700 "$JOB_HOME" "$JOB_TMP"

ENVIRONMENT=(
    -i
    "PATH=$PROFILE/bin"
    "GUIX_PROFILE=$PROFILE"
    "GUIX_PYTHONPATH=$PYTHON_SITE"
    "PYTHONPATH=${TIER3_ROOT:-$PWD}"
    "HOME=$JOB_HOME"
    "TMPDIR=$JOB_TMP"
    "LANG=C"
    "LC_ALL=C"
    "PYTHONUTF8=1"
    "TZ=UTC"
)

while IFS='=' read -r NAME VALUE; do
    case "$NAME" in
        SLURM_*|TIER3_*|SCRATCH) ENVIRONMENT+=("$NAME=$VALUE") ;;
    esac
done < <(env)

env "${ENVIRONMENT[@]}" python3 "$TIER3_ROOT/analysis/run_tier3.py" \
    audit-environment "$ENVIRONMENT_RECORD" >/dev/null
exec env "${ENVIRONMENT[@]}" "$@"
