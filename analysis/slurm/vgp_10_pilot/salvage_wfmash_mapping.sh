#!/usr/bin/env bash
# Salvage a completed wfmash mapping whose sbatch run died at the (pre-fix)
# fastga_scratch_guard finalize: the sweepga/wfmash mapping itself completed
# ("sweepga finished (guard_failed=0)", sweepga_status=0), the node-local
# scratch was retained as /scratch/vgp-map-<PAIR>-<JOBID>-*.failed-<JOBID>,
# and only the post-sweepga finalize stages never ran.
#
# Usage: salvage_wfmash_mapping.sh <PAIR> <JOBID> <NODE>
#
# Re-executes the exact post-sweepga stages of mapping_stage_wfmash.sh
# (enforce-paf, audit-paf, exclusion beds, mapping execution record,
# cross-filesystem promotion) inside the retained scratch on $NODE, marking
# every record salvage_recovered=true with the original job id and the guard
# failure reason. Idempotent: a FINAL with .complete.json short-circuits.

set -euo pipefail

PAIR=${1:?usage: salvage_wfmash_mapping.sh <PAIR> <JOBID> <NODE>}
JOBID=${2:?usage: salvage_wfmash_mapping.sh <PAIR> <JOBID> <NODE>}
NODE=${3:?usage: salvage_wfmash_mapping.sh <PAIR> <JOBID> <NODE>}

ROOT=/moosefs/erikg/lewontin-paradox
VGP_DATA_ROOT=/moosefs/erikg/vgp
RUN_ROOT=$VGP_DATA_ROOT/pilot/three-pair/vgp-wave-20260907-v1
INPUT_DIR=$VGP_DATA_ROOT/pilot/inputs/$PAIR
FINAL=$RUN_ROOT/$PAIR/mapping
VGP_TASK_ID=vgp-wave-20260907-v1
SUPPLEMENTAL_PROFILE=/gnu/store/8x4hx7d9hnv187yprjrzqyg0kxj2z32k-profile
GUARD_FAILURE_NOTE="fastga_scratch_guard finalize rejected the wfmash run (no live FastGA /proc snapshot); mapping itself completed successfully"

[[ -f $INPUT_DIR/input-manifest.json ]] || { echo "input manifest absent for $PAIR" >&2; exit 3; }
[[ -f $ROOT/analysis/sweepga_origin_main_build.json ]] || { echo "build provenance absent" >&2; exit 3; }
if [[ -f $FINAL/.complete.json ]]; then
    echo "RESUME: $PAIR mapping already complete: $FINAL"
    exit 0
fi

echo "[SALVAGE $PAIR job=$JOBID node=$NODE] locating retained scratch"
SCRATCH=$(ssh -o BatchMode=yes -o ConnectTimeout=10 "$NODE" \
    "ls -d /scratch/vgp-map-${PAIR}-${JOBID}-*.failed-${JOBID} \
           /mnt/sdb1/scratch/vgp-map-${PAIR}-${JOBID}-*.failed-${JOBID} 2>/dev/null | head -1")
[[ -n $SCRATCH ]] || { echo "retained scratch not found on $NODE for $PAIR/$JOBID" >&2; exit 4; }
echo "[SALVAGE $PAIR] scratch=$SCRATCH"

# ---- remote finalize: verbatim post-sweepga stages of mapping_stage_wfmash.sh
GUARD_NOTE_B64=$(printf %s "$GUARD_FAILURE_NOTE" | base64 -w0)
ssh -o BatchMode=yes -o ConnectTimeout=10 "$NODE" bash -s -- \
    "$PAIR" "$JOBID" "$SCRATCH" "$ROOT" \
    "$INPUT_DIR" "$FINAL" "$VGP_DATA_ROOT" \
    "$VGP_TASK_ID" "$SUPPLEMENTAL_PROFILE" "$GUARD_NOTE_B64" <<'REMOTE'
set -euo pipefail
PAIR=$1 JOBID=$2 SCRATCH=$3 ROOT=$4 INPUT_DIR=$5 FINAL=$6 DATA_ROOT=$7 TASK_ID=$8 PROFILE=$9
GUARD_NOTE=$(printf %s "${10}" | base64 -d)
export PATH="$PROFILE/bin:$PATH"
export PYTHONPATH="$ROOT"
export SLURM_JOB_ID="$JOBID" VGP_TASK_ID="$TASK_ID"
scratch=$SCRATCH
partial=$scratch/mapping.partial
native_paf=$partial/h2_to_h1.native.1to1.paf

[[ -s $native_paf ]] || { echo "native PAF missing/empty: $native_paf" >&2; exit 5; }
[[ -f $INPUT_DIR/input-manifest.json ]] || { echo "input manifest not visible on node" >&2; exit 5; }
[[ -d $partial ]] || { echo "mapping.partial absent" >&2; exit 5; }
[[ ! -e $FINAL ]] || { echo "FINAL already exists (without sentinel?) — refusing" >&2; exit 6; }

echo "[remote] enforce 1:1 multiplicity + audit PAF"
python3 -m analysis.vgp_10_pilot enforce-paf \
    "$native_paf" "$partial/h2_to_h1.1to1.paf" \
    >"$partial/exact_multiplicity_filter.json"
python3 -m analysis.vgp_10_pilot audit-paf \
    "$partial/h2_to_h1.1to1.paf" "$scratch/inputs/h1.fa" "$scratch/inputs/h2.fa" \
    >"$partial/multiplicity.json"

echo "[remote] derive exclusion beds"
python3 - "$partial/h2_to_h1.1to1.paf" "$scratch/inputs/h1.fa" \
    "$scratch/inputs/h2.fa" "$INPUT_DIR/h1_universe.bed" "$partial" <<'PY'
import sys
from pathlib import Path
from analysis.vgp_10_pilot import (
    low_complexity_intervals, non_acgt_intervals, paf_h1_intervals,
    parse_fasta, parse_paf, project_h2_non_acgt_to_h1, read_bed,
    subtract_intervals, write_bed,
)
paf, h1, h2, universe, out = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], Path(sys.argv[5])
records = parse_paf(paf)
one_to_one = paf_h1_intervals(records)
h1_sequences = parse_fasta(h1)
write_bed(out / "h1.1to1.bed", one_to_one)
write_bed(out / "not_1to1.bed", subtract_intervals(read_bed(universe), one_to_one))
write_bed(out / "h1_gap_or_N.bed", non_acgt_intervals(h1_sequences))
write_bed(out / "h2_gap_or_N.bed", project_h2_non_acgt_to_h1(records, parse_fasta(h2)))
write_bed(out / "repeat_or_low_complexity_primary.bed", low_complexity_intervals(h1_sequences))
PY

echo "[remote] write mapping execution record (salvage-annotated)"
python3 - "$partial/mapping_execution.json" "$DATA_ROOT" "$PAIR" \
    "$ROOT/analysis/sweepga_origin_main_build.json" "$JOBID" \
    "$INPUT_DIR/input-manifest.json" "$GUARD_NOTE" <<'PY'
import hashlib, json, os, sys
from pathlib import Path
from analysis.vgp_10_pilot import sha256_file

record_path, data_root, selection, build_json_path, job_id, manifest_path, guard_note = sys.argv[1:8]
build = json.load(open(build_json_path))["binary"]
sweepga_path = build["build_1_realpath"]
sweepga_sha = build["sha256_build_1"]
manifest = json.load(open(manifest_path))
Path(record_path).write_text(json.dumps({
    "schema_version": "vgp-wave-mapping-execution-v1",
    "task_id": os.environ.get("VGP_TASK_ID", "vgp-wave-20260907-v1"),
    "authorization_id": manifest.get("authorization_id"),
    "canonical_vgp_root": data_root,
    "selection_id": selection, "slurm_job_id": job_id,
    "mapping_backend": "wfmash",
    "aligner_provenance": "sweepga origin/main pinned build; wfmash backend, remap-lineage flag semantics (tier3a-origin-remap-20260716 command.txt), full-node exclusive threading (operator standardization 2026-09-19)",
    "required_option": "--num-mappings 1:1",
    "pinned_binary": {"path": sweepga_path, "sha256": sweepga_sha,
                      "build_provenance_json": {"path": build_json_path,
                                                "sha256": sha256_file(build_json_path)}},
    "query_target_orientation": "h2_query_h1_target",
    "input_manifest": {"path": manifest_path, "sha256": sha256_file(manifest_path)},
    "h1_accession_version": manifest.get("h1_accession_version"),
    "h2_accession_version": manifest.get("h2_accession_version"),
    "fallback_is_species_exclusion": False,
    "salvage_recovered": True,
    "salvage": {
        "original_slurm_job_id": job_id,
        "guard_failure": guard_note,
        "finalized_by": "salvage_wfmash_mapping.sh (verbatim post-sweepga stages)",
        "note": "sweepga/wfmash completed with status 0 before the incompatible fastga-scratch-guard finalize; post-sweepga stages re-executed verbatim from the retained scratch",
    },
}, sort_keys=True) + "\n")
PY

echo "[remote] record wfmash scratch contract (salvage form)"
python3 - "$partial/sweepga_scratch_contract.json" "$scratch" "$GUARD_NOTE" <<'PY'
import json, sys
from pathlib import Path
audit = Path(sys.argv[1]).parent / "sweepga_scratch_snapshots.jsonl"
n = sum(1 for line in audit.read_text().splitlines() if line.strip()) if audit.exists() else 0
Path(sys.argv[1]).write_text(json.dumps({
    "schema_version": "wfmash-scratch-contract-v1",
    "backend": "wfmash",
    "scratch": sys.argv[2],
    "snapshots_observed_by_fastga_guard_check": n,
    "escapes": 0,
    "verdict": "salvage_recovered_finalize_skipped",
    "note": sys.argv[3],
}, sort_keys=True) + "\n")
PY

echo "[remote] promote to $FINAL"
python3 - "$partial" "$FINAL" "$PAIR" "$JOBID" "$DATA_ROOT" "$GUARD_NOTE" <<'PY'
import json, sys
from pathlib import Path
from analysis.vgp_10_pilot import promote_stage_cross_filesystem
partial, final, selection, job_id, data_root, guard_note = sys.argv[1:7]
promote_stage_cross_filesystem(Path(partial), Path(final), {
    "selection_id": selection, "stage": "mapping", "atomic_promotion": True,
    "slurm_job_id": job_id, "canonical_vgp_root": data_root,
    "salvage_recovered": True, "salvage_guard_failure": guard_note,
}, job_id)
print(json.dumps({"promoted_to": str(final)}))
PY

echo "[remote] telemetry record"
mkdir -p "$(dirname "$FINAL")/telemetry-wfmash"
python3 - "$FINAL/../telemetry-wfmash/mapping.${JOBID}.json" "$PAIR" "$JOBID" <<'PY'
import json, sys, time
from pathlib import Path
Path(sys.argv[1]).write_text(json.dumps({
    "selection_id": sys.argv[2], "stage": "mapping",
    "job_id": sys.argv[3], "disposition": "salvage_recovered",
    "ended_epoch": int(time.time()),
    "note": "mapping completed on-node; post-sweepga finalize re-executed by salvage_wfmash_mapping.sh",
}, sort_keys=True) + "\n")
PY
echo "[remote] done"
REMOTE

echo "[SALVAGE $PAIR] verification"
PAF=$FINAL/h2_to_h1.1to1.paf
[[ -s $PAF ]] || { echo "salvage verification failed: $PAF missing" >&2; exit 7; }
echo "rows:   $(wc -l < "$PAF")"
echo "sha256: $(sha256sum "$PAF" | cut -d' ' -f1)"
ls "$FINAL"
[[ -f $FINAL/.complete.json ]] && echo "sentinel: present"
[[ -f $FINAL/mapping_execution.json ]] && grep -o '"salvage_recovered": true' "$FINAL/mapping_execution.json"
echo "[SALVAGE $PAIR] OK"
