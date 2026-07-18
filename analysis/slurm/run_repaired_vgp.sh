#!/usr/bin/env bash
set -euo pipefail

# Dormant under the committed NO_GO.  The login-node runner supplies only
# gate-bound arguments and uses --export=NONE.  This worker never downloads.
RUN_ID=
GATE_SHA256=
AUTHORIZATION_TUPLE=
ACQUISITION_MANIFEST=
INVENTORY=
while (( $# )); do
    case "$1" in
        --run-id) RUN_ID=$2; shift 2 ;;
        --gate-sha256) GATE_SHA256=$2; shift 2 ;;
        --authorization-tuple) AUTHORIZATION_TUPLE=$2; shift 2 ;;
        --acquisition-manifest) ACQUISITION_MANIFEST=$2; shift 2 ;;
        --inventory) INVENTORY=$2; shift 2 ;;
        *) echo "unrecognized argument: $1" >&2; exit 64 ;;
    esac
done

fail() { echo "ERROR: $*" >&2; exit 65; }
[[ -n $RUN_ID && $GATE_SHA256 =~ ^[0-9a-f]{64}$ && $AUTHORIZATION_TUPLE =~ ^[0-9a-f]{64}$ ]] || fail "missing gate-bound identity"
[[ -n ${SLURM_JOB_ID:-} && -n ${SLURM_ARRAY_TASK_ID:-} ]] || fail "must run as a Slurm array element"
[[ -d ${SLURM_TMPDIR:-} ]] || fail "SLURM_TMPDIR is required"
[[ -r $ACQUISITION_MANIFEST && -r $INVENTORY ]] || fail "acquisition ledgers are unavailable"

# Explicit network prohibition.  Inputs must already be immutable promoted
# objects; no provider URL is dereferenced here.
unset http_proxy https_proxy HTTP_PROXY HTTPS_PROXY ALL_PROXY all_proxy NO_PROXY no_proxy
export CARGO_NET_OFFLINE=true
export GIT_CONFIG_NOSYSTEM=1
export GIT_TERMINAL_PROMPT=0

WORK=$SLURM_TMPDIR/$RUN_ID/${SLURM_ARRAY_TASK_ID}
install -d -m 700 "$WORK/inputs" "$WORK/output.partial" "$WORK/telemetry"

# Select one candidate deterministically and copy only its inventory objects to
# node-local scratch.  Python verifies ledger identity; sha256sum verifies the
# bytes after staging.  Absolute provider URLs are retained as evidence only.
/usr/bin/python3 - "$INVENTORY" "$SLURM_ARRAY_TASK_ID" "$AUTHORIZATION_TUPLE" "$WORK/input.tsv" <<'PY'
import csv
import sys
from pathlib import Path

inventory, index, expected_tuple, output = Path(sys.argv[1]), int(sys.argv[2]), sys.argv[3], Path(sys.argv[4])
with inventory.open(encoding="utf-8", newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
candidates = sorted({row["candidate_id"] for row in rows})
if index < 0 or index >= len(candidates):
    raise SystemExit("array index is outside exact promoted candidate set")
selected = [row for row in rows if row["candidate_id"] == candidates[index]]
if not selected or any(row["authorization_tuple_digest"] != expected_tuple for row in selected):
    raise SystemExit("inventory authorization tuple mismatch")
with output.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=selected[0].keys(), delimiter="\t", lineterminator="\n")
    writer.writeheader(); writer.writerows(selected)
PY

while IFS=$'\t' read -r object_path local_sha256 asset_role; do
    [[ -f $object_path && $local_sha256 =~ ^[0-9a-f]{64}$ ]] || fail "invalid promoted object row"
    destination=$WORK/inputs/${asset_role}
    cp --reflink=auto -- "$object_path" "$destination"
    [[ $(sha256sum "$destination" | awk '{print $1}') == "$local_sha256" ]] || fail "staged local payload SHA-256 mismatch"
    chmod 0440 "$destination"
done < <(/usr/bin/python3 - "$WORK/input.tsv" <<'PY'
import csv, sys
with open(sys.argv[1], encoding="utf-8", newline="") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        print(row["object_path"], row["local_sha256"], row["asset_role"], sep="\t")
PY
)

START_EPOCH=$(date +%s)
/usr/bin/time -v -o "$WORK/telemetry/time.txt" \
    "$PWD/analysis/slurm/guix_job.sh" "$PWD/analysis/tier3_environment.json" \
    python3 "$PWD/analysis/vgp_pilot_compute.py" \
        --input-ledger "$WORK/input.tsv" --staged-root "$WORK/inputs" \
        --output "$WORK/output.partial" --run-id "$RUN_ID" \
        --gate-sha256 "$GATE_SHA256" --authorization-tuple "$AUTHORIZATION_TUPLE"
END_EPOCH=$(date +%s)

# The compute program must have validated these scientific contracts before a
# result can be promoted.  Missing/false sentinels are terminal exclusions.
/usr/bin/python3 - "$WORK/output.partial/success.json" <<'PY'
import json, sys
from pathlib import Path
p = json.loads(Path(sys.argv[1]).read_text())
required = (
    "inputs_rehashed", "native_h1_annotation_exact", "sweepga_whole_haplotype_1to1",
    "impg_native_partitions", "normalized_vcf_tbi", "normalized_bcf_csi",
    "denominators_measured", "thresholds_enforced", "no_demographic_inference",
)
if any(p.get(key) is not True for key in required):
    raise SystemExit("compute success sentinel is incomplete")
PY

PROMOTION_ROOT=${VGP_RUN_PROMOTION_ROOT:?VGP_RUN_PROMOTION_ROOT must be set by the approved Slurm submission}/$RUN_ID
install -d -m 2750 "$PROMOTION_ROOT/.incoming"
TARGET=$PROMOTION_ROOT/${SLURM_ARRAY_TASK_ID}
PART=$PROMOTION_ROOT/.incoming/${SLURM_ARRAY_TASK_ID}.${SLURM_JOB_ID}.partial
cp -a -- "$WORK/output.partial" "$PART"
find "$PART" -type f -print0 | sort -z | xargs -0 sha256sum > "$PART/output.sha256"
sync "$PART"
[[ ! -e $TARGET ]] || fail "atomic promotion target already exists"
mv -T "$PART" "$TARGET"
printf '{"elapsed_seconds":%d,"slurm_job_id":"%s","array_task_id":"%s"}\n' \
    "$((END_EPOCH - START_EPOCH))" "$SLURM_JOB_ID" "$SLURM_ARRAY_TASK_ID" > "$TARGET/job_telemetry.json"
