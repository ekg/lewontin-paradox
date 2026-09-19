#!/usr/bin/env bash
#SBATCH --job-name=vgp-map-wfmash
#SBATCH --partition=workers
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --exclude=octopus11
# 2026-09-16: octopus11 empirically confirmed broken — job 2834221 assigned,
# RUNNING per slurm, but scontrol listpids empty, 0-byte logs, never executed.
# The legacy exclusion was justified; restored.
# 2026-09-17: --exclusive added — job 2837797 died when the SHARED node-local
# scratch (89 co-tenant dirs, /mnt/sdb1) transiently filled during FastGA's
# alignment phase ("IO write to .../sweepga-temp/_pair.*.N failed"). Exclusive
# allocation removes concurrent slurm scratch pressure on the node.
#SBATCH --cpus-per-task=48
#SBATCH --mem=200G
#SBATCH --time=24:00:00
#SBATCH --export=NONE
#SBATCH --output=/moosefs/erikg/vgp/logs/vgp-map-%x-%j.out
#SBATCH --error=/moosefs/erikg/vgp/logs/vgp-map-%x-%j.err

set -euo pipefail

fail() {
    echo "ERROR: $*" >&2
    exit 2
}

stage() {
    # Progress markers land in the slurm .out DURING the run (promoted-only
    # logging left operators blind on two failures).
    LAST_STAGE="$*"
    echo "[STAGE $(date -u +%FT%TZ)] $*"
}

# Whole-assembly mapping for the scale-out wave — WFMASH BACKEND, FULL NODE.
# 2026-09-19 operator standardization: wfmash becomes the wave backend at
# full-node threading (exclusive 48c). Evidence: fastga failed/hung on
# P09x2, P06, P10, P08 (IO write to _pair.* / 21h single-core hang; see
# /moosefs/erikg/vgp/derived/read-validation/runs/fastga-failure-evidence/)
# while succeeding only on P02/P07/P01/P05; wfmash completed every pair it
# was given (tier3a remap trio) and its one failure (P09 >48h at 32 threads)
# is addressed by full-node threading (~1.5x throughput => ~32h projection).
# Flag set = the remap lineage (command.txt of tier3a-origin-remap-20260716):
# --aligner wfmash --num-mappings 1:1 --scaffold-jump 0 --overlap 0.95
# --map-pct-identity 90 --min-aln-length 25k --scoring log-length-ani.
# Produces exactly the stage contracts analysis/slurm/run_vgp_bounded_pair.sh
# consumes (h2_to_h1.1to1.paf, h1.1to1.bed, exclusion beds, multiplicity
# audit, digest-complete .complete.json).
# Parameterized via VGP_SELECTION / VGP_INPUT_DIR / VGP_RUN_ROOT /
# VGP_SUPERSEDES / VGP_LINEAGE_NOTE; sbatch directives overridable per pair
# via CLI (-p / -c / --mem / -t / -J / --export).

if [[ ${1:-} != --inside-guix ]]; then
    ROOT=${SLURM_SUBMIT_DIR:-$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)}
    export TIER3_ROOT="$ROOT"
    # guix_job.sh scrubs the environment (only SLURM_*/TIER3_*/SCRATCH pass
    # through) — forward VGP_* parameterization knobs as TIER3_VGP_*.
    while IFS='=' read -r _name _value; do
        case $_name in VGP_*) export "TIER3_${_name}=${_value}" ;; esac
    done < <(env)
    exec "$ROOT/analysis/slurm/guix_job.sh" \
        "$ROOT/analysis/pilot_results/guix_environment.json" \
        bash "$0" --inside-guix
fi

ROOT=${SLURM_SUBMIT_DIR:-$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)}
export PYTHONPATH="$ROOT"
SUPPLEMENTAL_PROFILE=/gnu/store/8x4hx7d9hnv187yprjrzqyg0kxj2z32k-profile
[[ -d $SUPPLEMENTAL_PROFILE ]] || fail "pinned supplemental Guix profile is unavailable"
export PATH="$SUPPLEMENTAL_PROFILE/bin:$PATH"
# Parameterization knobs survive the guix re-exec via TIER3_VGP_* forwarding
# (guix_job.sh passes TIER3_* through); direct VGP_* reads kept for
# --inside-guix invocations (e.g. the instrumented repro).
SELECTION=${TIER3_VGP_SELECTION:-${VGP_SELECTION:-P09}}
VGP_DATA_ROOT=${TIER3_VGP_DATA_ROOT:-${VGP_DATA_ROOT:-/moosefs/erikg/vgp}}
INPUT_DIR=${TIER3_VGP_INPUT_DIR:-${VGP_INPUT_DIR:-$VGP_DATA_ROOT/pilot/inputs/$SELECTION}}
VGP_RUN_ID=${TIER3_VGP_RUN_ID:-${VGP_RUN_ID:-vgp-wave-20260907-v1}}
RUN_ROOT=${TIER3_VGP_RUN_ROOT:-${VGP_RUN_ROOT:-$VGP_DATA_ROOT/pilot/runs/$VGP_RUN_ID}}
VGP_TASK_ID=${TIER3_VGP_TASK_ID:-${VGP_TASK_ID:-vgp-wave-20260907-v1}}
VGP_NODE_LOCAL_BASE=${TIER3_VGP_NODE_LOCAL_BASE:-${VGP_NODE_LOCAL_BASE:-/scratch}}
BUILD_JSON=${TIER3_VGP_BUILD_JSON:-${VGP_BUILD_JSON:-$ROOT/analysis/sweepga_origin_main_build.json}}
PAIR_RUN="$RUN_ROOT/$SELECTION"
FINAL="$PAIR_RUN/mapping"
TELEMETRY="$PAIR_RUN/telemetry-wfmash"

LAST_STAGE=init
scratch=""
required_bytes=0
available_bytes=0

[[ $SLURM_JOB_ID ]] || fail "submit this mapping with Slurm"
[[ -d $VGP_NODE_LOCAL_BASE && -w $VGP_NODE_LOCAL_BASE ]] || fail "scratch unavailable"
stage "scratch filesystem guard ($VGP_NODE_LOCAL_BASE)"
case $(stat -f -c %T -- "$VGP_NODE_LOCAL_BASE") in
    nfs|nfs4|fuse*|lustre|gpfs|ceph|mfs|moosefs) fail "scratch is not node-local" ;;
esac

if [[ -f $FINAL/.complete.json ]]; then
    echo "RESUME: $SELECTION mapping already complete: $FINAL"
    exit 0
fi
[[ ! -e $FINAL ]] || fail "mapping final exists without completion sentinel: $FINAL"
[[ -f $INPUT_DIR/input-manifest.json ]] || fail "$SELECTION input manifest is absent"
[[ -f $BUILD_JSON ]] || fail "sweepga origin/main build provenance is absent"
mkdir -p "$TELEMETRY" "$PAIR_RUN"
started=$(date +%s)

stage "resolve pinned sweepga binary ($BUILD_JSON)"
readarray -t pinned < <(python3 - "$BUILD_JSON" <<'PY'
import hashlib, json, sys
from pathlib import Path
binary = json.load(open(sys.argv[1]))["binary"]
path = Path(binary["build_1_realpath"])
expected = binary["sha256_build_1"]
if not path.is_file():
    raise SystemExit(f"pinned sweepga binary is absent: {path}")
digest = hashlib.sha256(path.read_bytes()).hexdigest()
if digest != expected:
    raise SystemExit("pinned sweepga binary sha256 mismatch")
print(path)
print(digest)
print(hashlib.sha256(Path(sys.argv[1]).read_bytes()).hexdigest())
PY
)
sweepga=${pinned[0]}
sweepga_sha=${pinned[1]}
build_json_sha=${pinned[2]}
[[ -x $sweepga ]] || fail "pinned sweepga is not executable: $sweepga"
[[ -x $(dirname "$sweepga")/wfmash ]] || fail "pinned wfmash companion is absent"
export WFMASH_BIN_DIR="$(dirname "$sweepga")"
export PATH="$(dirname "$sweepga"):$PATH"
export TMPDIR="$VGP_NODE_LOCAL_BASE" TMP="$VGP_NODE_LOCAL_BASE" TEMP="$VGP_NODE_LOCAL_BASE"

stage "read input sizes ($INPUT_DIR)"
# Test seam: repro/small-pair runs may lower the staged-scratch floor;
# production default is 64 GiB.
scratch_floor_bytes=${TIER3_VGP_STAGED_SCRATCH_FLOOR_BYTES:-${VGP_STAGED_SCRATCH_FLOOR_BYTES:-$((64 * 2**30))}}
h1_bytes=$(stat -c %s -- "$INPUT_DIR/h1.fa")
h2_bytes=$(stat -c %s -- "$INPUT_DIR/h2.fa")
readarray -t plan < <(python3 - "$h1_bytes" "$h2_bytes" "$INPUT_DIR" "$scratch_floor_bytes" <<'PY'
import json, sys
resource = json.load(open(sys.argv[3] + "/resources.json"))["stages"]["mapping"]
staged = (int(sys.argv[1]) + int(sys.argv[2])) * 2
print(max(int(resource["scratch_bytes_high"]), staged + int(sys.argv[4])))
print(int(resource["cpus_per_task"]))
PY
)
required_bytes=${plan[0]}

available_bytes=$(df -PB1 -- "$VGP_NODE_LOCAL_BASE" | awk 'NR==2 {print $4}')
(( available_bytes >= required_bytes )) || \
    fail "insufficient measured scratch: $available_bytes < $required_bytes"
scratch=$(mktemp -d -- "$VGP_NODE_LOCAL_BASE/vgp-map-$SELECTION-${SLURM_JOB_ID}-XXXXXX")
partial="$scratch/mapping.partial"
mkdir -p "$partial" "$scratch/inputs"
export TMPDIR="$scratch" TMP="$scratch" TEMP="$scratch"

finish() {
    status=$?
    if (( status != 0 )); then
        stage "FAILURE (exit $status) after: $LAST_STAGE — writing failure.json + retaining scratch"
        python3 - "$TELEMETRY/mapping.${SLURM_JOB_ID}.json" "$partial" "$status" \
            "$started" "$required_bytes" "$available_bytes" "$VGP_DATA_ROOT" \
            "$SELECTION" "$LAST_STAGE" <<'PY' || true
import json, sys, time
from pathlib import Path
out, partial = Path(sys.argv[1]), Path(sys.argv[2])
diagnostics = {}
if partial.exists():
    for path in partial.rglob("*"):
        if path.is_file() and path.suffix in {".stderr", ".stdout", ".log", ".jsonl"}:
            diagnostics[str(path.relative_to(partial))] = \
                path.read_text(errors="replace")[-32768:]
out.write_text(json.dumps({
    "selection_id": sys.argv[8], "stage": "mapping",
    "job_id": out.stem.split(".")[-1], "disposition": "failure",
    "exit_status": int(sys.argv[3]), "started_epoch": int(sys.argv[4]),
    "ended_epoch": int(time.time()),
    "scratch_required_bytes": int(sys.argv[5]),
    "scratch_available_bytes_at_start": int(sys.argv[6]),
    "canonical_vgp_root": sys.argv[7], "last_stage": sys.argv[9],
    "diagnostic_tails": diagnostics,
}, sort_keys=True) + "\n")
PY
        # Failure evidence: retain the scratch workdir (rename = atomic, same
        # filesystem) instead of deleting it. Cleanup is a manual operator step.
        if [[ -n $scratch && -d $scratch ]]; then
            retained="${scratch}.failed-${SLURM_JOB_ID}"
            rm -rf -- "$retained" 2>/dev/null || true
            mv -- "$scratch" "$retained" 2>/dev/null || true
            echo "[STAGE $(date -u +%FT%TZ)] scratch retained at: $retained"
        fi
        exit "$status"
    fi
    [[ -z $scratch || $scratch == "$VGP_NODE_LOCAL_BASE/vgp-map-$SELECTION-${SLURM_JOB_ID}-"* ]] && \
        rm -rf -- "$scratch"
    exit "$status"
}
trap finish EXIT

stage "stage inputs into scratch ($INPUT_DIR -> $scratch/inputs)"
cp -- "$INPUT_DIR/h1.fa" "$scratch/inputs/h1.fa"
cp -- "$INPUT_DIR/h2.fa" "$scratch/inputs/h2.fa"
stage "digest audit vs input manifest"
python3 - "$INPUT_DIR/input-manifest.json" "$scratch/inputs/h1.fa" \
    "$scratch/inputs/h2.fa" "$partial/input_digest_audit.json" "$SELECTION" <<'PY'
import hashlib, json, sys
from pathlib import Path
manifest = json.load(open(sys.argv[1]))
audit = {"schema_version": "vgp-wave-mapping-input-digest-audit-v1",
         "selection_id": sys.argv[5], "roles": {}}
for role, staged in (("h1_fasta", Path(sys.argv[2])), ("h2_fasta", Path(sys.argv[3]))):
    asset = manifest["assets"][role]
    digest = hashlib.sha256(staged.read_bytes()).hexdigest()
    if digest != asset["sha256"]:
        raise SystemExit(f"staged {role} sha256 mismatch against input manifest")
    if staged.stat().st_size != asset["size_bytes"]:
        raise SystemExit(f"staged {role} size mismatch against input manifest")
    audit["roles"][role] = {"source_path": asset["path"],
                            "source_sha256": asset["sha256"],
                            "source_size_bytes": asset["size_bytes"],
                            "staged_sha256": digest,
                            "staged_size_bytes": staged.stat().st_size}
Path(sys.argv[4]).write_text(json.dumps(audit, sort_keys=True) + "\n")
PY

cd -- "$scratch"
scratch_resolved=$(readlink -f -- "$scratch")
cwd_resolved=$(readlink -f -- /proc/$$/cwd)
[[ $cwd_resolved == "$scratch_resolved" ]] || \
    fail "batch cwd outside private node-local scratch: $cwd_resolved"

threads=${SLURM_CPUS_PER_TASK:-32}
native_paf="$partial/h2_to_h1.native.1to1.paf"
stage "write command record + gate checks"
printf '%q ' "$sweepga" "$scratch/inputs/h2.fa" "$scratch/inputs/h1.fa" \
    --output-file "$native_paf" --aligner wfmash --num-mappings 1:1 --scaffold-jump 0 \
    --overlap 0.95 --map-pct-identity 90 --min-aln-length 25k \
    --scoring log-length-ani --threads "$threads" \
    --temp-dir "$scratch/sweepga-temp" > "$partial/command.txt"
printf '\n' >> "$partial/command.txt"
printf '%s  %s\n' "$sweepga_sha" "$sweepga" > "$partial/sweepga.sha256"
grep -F -- '--num-mappings 1:1' "$partial/command.txt" >/dev/null
grep -F -- '--aligner wfmash' "$partial/command.txt" >/dev/null
grep -F "$sweepga_sha" "$partial/sweepga.sha256" >/dev/null
mkdir -p "$scratch/sweepga-temp"

stage "launch sweepga/wfmash (pid written when backgrounded)"
"$sweepga" "$scratch/inputs/h2.fa" "$scratch/inputs/h1.fa" \
    --output-file "$native_paf" --aligner wfmash --num-mappings 1:1 --scaffold-jump 0 \
    --overlap 0.95 --map-pct-identity 90 --min-aln-length 25k \
    --scoring log-length-ani --threads "$threads" \
    --temp-dir "$scratch/sweepga-temp" \
    >"$partial/sweepga.stdout" 2>"$partial/sweepga.stderr" &
sweepga_pid=$!
echo "[STAGE $(date -u +%FT%TZ)] sweepga pid=$sweepga_pid scratch=$scratch"
guard_log="$partial/sweepga_scratch_snapshots.jsonl"
# wfmash-backend containment guard (equivalent rigor to fastga_scratch_guard):
# every poll inspects the whole sweepga process tree (sweepga + wfmash children)
# for scratch escape via cwd or writable fd targets, and samples RSS into the
# same jsonl. fastga_scratch_guard.py cannot be used here: its finalize
# requires a live FastGA /proc snapshot, which wfmash runs never produce
# (failure mode observed on jobs 2841278/2841279: mapping completed, finalize
# rejected). Same failure semantics: any escape -> guard_failed=1 -> kill.
guard_failed=0
wfmash_guard_check() {
    python3 - "$sweepga_pid" "$scratch" "$guard_log" <<'PY'
import json, os, sys, time
from pathlib import Path

sweepga_pid, scratch, audit = int(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
skip = {"/dev", "/proc", "pipe:", "socket:", "anon_inode:", "null"}

def tree_pids(root):
    pids, frontier = {root}, [root]
    while frontier:
        pid = frontier.pop()
        try:
            kids = [int(x) for x in Path(f"/proc/{pid}/task/{pid}/children").read_text().split()]
        except OSError:
            kids = []
        for kid in kids:
            if kid not in pids:
                pids.add(kid); frontier.append(kid)
    return sorted(pids)

def resolve(fd_or_cwd):
    try:
        return os.path.realpath(os.readlink(fd_or_cwd))
    except OSError:
        return None

samples, escapes = [], []
for pid in tree_pids(sweepga_pid):
    comm = "unknown"
    try: comm = Path(f"/proc/{pid}/comm").read_text().strip()
    except OSError: continue
    cwd = resolve(f"/proc/{pid}/cwd")
    try:
        rss_kb = int(next(
            line.split()[1] for line in Path(f"/proc/{pid}/status").read_text().splitlines()
            if line.startswith("VmRSS")))
    except (OSError, StopIteration):
        rss_kb = None
    cwd_ok = cwd is None or cwd.startswith(str(scratch))
    fd_escape = []
    try:
        for fd in Path(f"/proc/{pid}/fd").iterdir():
            target = resolve(str(fd))
            if target is None:
                continue
            if any(target.startswith(s) for s in skip):
                continue
            if not target.startswith(str(scratch)):
                fd_escape.append(target)
    except OSError:
        pass
    if not cwd_ok:
        escapes.append({"pid": pid, "kind": "cwd", "path": cwd})
    if fd_escape:
        escapes.append({"pid": pid, "kind": "fd", "paths": fd_escape[:5]})
    samples.append({"pid": pid, "comm": comm, "rss_kb": rss_kb,
                    "cwd": cwd, "cwd_ok": cwd_ok, "fd_escape": fd_escape})
record = {"ts": time.time(), "sweepga_pid": sweepga_pid,
          "backend": "wfmash", "samples": samples, "escapes": escapes}
with audit.open("a") as handle:
    handle.write(json.dumps(record, sort_keys=True) + "\n")
if escapes:
    print(json.dumps({"verdict": "escape", "escapes": escapes}), file=sys.stderr)
    raise SystemExit(1)
PY
}
while [[ $(ps -o stat= -p "$sweepga_pid" 2>/dev/null) != *Z* ]] && \
        kill -0 "$sweepga_pid" 2>/dev/null; do
    if ! wfmash_guard_check; then
        guard_failed=1
        pkill -TERM -P "$sweepga_pid" 2>/dev/null || true
        kill -TERM "$sweepga_pid" 2>/dev/null || true
        break
    fi
    sleep "${VGP_GUARD_INTERVAL_SECONDS:-5}"
done
stage "sweepga finished (guard_failed=$guard_failed)"
if wait "$sweepga_pid"; then sweepga_status=0; else sweepga_status=$?; fi
(( guard_failed == 0 )) || \
    fail "hard infrastructure error: aligner escaped private node-local scratch"
(( sweepga_status == 0 )) || exit "$sweepga_status"
# wfmash finalize: equivalent of fastga_scratch_guard finalize — requires at
# least one containment snapshot and zero escapes, then writes the contract.
python3 - "$scratch" "$guard_log" "$partial/sweepga_scratch_contract.json" <<'PY'
import json, sys
from pathlib import Path
scratch, audit, output = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
records = [json.loads(line) for line in audit.read_text().splitlines() if line.strip()]
if not records:
    raise SystemExit("wfmash scratch guard finalize: no containment snapshots observed")
escapes = [e for r in records for e in r.get("escapes", [])]
if escapes:
    raise SystemExit(f"wfmash scratch guard finalize: {len(escapes)} escape(s) recorded")
output.write_text(json.dumps({
    "schema_version": "wfmash-scratch-contract-v1",
    "backend": "wfmash",
    "scratch": str(scratch),
    "snapshots": len(records),
    "escapes": 0,
    "verdict": "contained",
}, sort_keys=True) + "\n")
PY

stage "enforce 1:1 multiplicity + audit PAF"
python3 -m analysis.vgp_10_pilot enforce-paf \
    "$native_paf" "$partial/h2_to_h1.1to1.paf" \
    >"$partial/exact_multiplicity_filter.json"
python3 -m analysis.vgp_10_pilot audit-paf \
    "$partial/h2_to_h1.1to1.paf" "$scratch/inputs/h1.fa" "$scratch/inputs/h2.fa" \
    >"$partial/multiplicity.json"
stage "derive exclusion beds"
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
stage "write mapping execution record"
python3 - "$partial/mapping_execution.json" "$VGP_DATA_ROOT" "$SELECTION" \
    "$BUILD_JSON" "$sweepga" "$sweepga_sha" "$build_json_sha" \
    "$INPUT_DIR/input-manifest.json" <<'PY'
import hashlib, json, os, sys
from pathlib import Path
from analysis.vgp_10_pilot import sha256_file
manifest = json.load(open(sys.argv[8]))
Path(sys.argv[1]).write_text(json.dumps({
    "schema_version": "vgp-wave-mapping-execution-v1",
    "task_id": os.environ.get("VGP_TASK_ID", "vgp-wave-20260907-v1"),
    "authorization_id": manifest.get("authorization_id"),
    "canonical_vgp_root": sys.argv[2],
    "selection_id": sys.argv[3], "slurm_job_id": os.environ["SLURM_JOB_ID"],
    "mapping_backend": "wfmash",
    "aligner_provenance": "sweepga origin/main pinned build; wfmash backend, remap-lineage flag semantics (tier3a-origin-remap-20260716 command.txt), full-node exclusive threading (operator standardization 2026-09-19)",
    "required_option": "--num-mappings 1:1",
    "pinned_binary": {"path": sys.argv[5], "sha256": sys.argv[6],
                      "build_provenance_json": {"path": sys.argv[4],
                                                "sha256": sys.argv[7]}},
    "query_target_orientation": "h2_query_h1_target",
    "input_manifest": {"path": sys.argv[8], "sha256": sha256_file(sys.argv[8])},
    "h1_accession_version": manifest.get("h1_accession_version"),
    "h2_accession_version": manifest.get("h2_accession_version"),
    "fallback_is_species_exclusion": False,
    **({"supersedes": os.environ["TIER3_VGP_SUPERSEDES"]} if os.environ.get("TIER3_VGP_SUPERSEDES") else {}),
    **({"lineage_note": os.environ["TIER3_VGP_LINEAGE_NOTE"]} if os.environ.get("TIER3_VGP_LINEAGE_NOTE") else {}),
}, sort_keys=True) + "\n")
PY
stage "promote to $FINAL"
python3 - "$partial" "$FINAL" "$SELECTION" "$SLURM_JOB_ID" "$VGP_DATA_ROOT" <<'PY'
import sys
from pathlib import Path
from analysis.vgp_10_pilot import promote_stage_cross_filesystem
promote_stage_cross_filesystem(Path(sys.argv[1]), Path(sys.argv[2]), {
    "selection_id": sys.argv[3], "stage": "mapping", "atomic_promotion": True,
    "slurm_job_id": sys.argv[4], "canonical_vgp_root": sys.argv[5],
}, sys.argv[4])
PY
stage "write telemetry"
python3 - "$TELEMETRY/mapping.${SLURM_JOB_ID}.json" "$started" "$required_bytes" \
    "$available_bytes" "$VGP_DATA_ROOT" "$SELECTION" <<'PY'
import json, os, sys, time
from pathlib import Path
Path(sys.argv[1]).write_text(json.dumps({
    "selection_id": sys.argv[6], "stage": "mapping",
    "job_id": os.environ["SLURM_JOB_ID"], "disposition": "success",
    "started_epoch": int(sys.argv[2]), "ended_epoch": int(time.time()),
    "scratch_required_bytes": int(sys.argv[3]),
    "scratch_available_bytes_at_start": int(sys.argv[4]),
    "canonical_vgp_root": sys.argv[5],
    "node_local_scratch_base": os.environ["VGP_NODE_LOCAL_BASE"],
}, sort_keys=True) + "\n")
PY
stage "COMPLETE: $FINAL"
printf '%s\n' "$FINAL"
