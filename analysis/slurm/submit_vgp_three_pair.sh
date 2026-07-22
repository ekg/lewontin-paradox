#!/usr/bin/env bash
set -euo pipefail

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
mode=${1:---dry-run}
[[ $mode == --dry-run || $mode == --submit ]] || { echo "usage: $0 [--dry-run|--submit]" >&2; exit 2; }
SELECTION="$ROOT/analysis/vgp_three_pair_selection_v1.json"
VGP_RUN_ID=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["run_id"])' "$SELECTION")
export VGP_RUN_ID VGP_DATA_ROOT=/moosefs/erikg/vgp
export VGP_ENVIRONMENT_CAPTURE="$ROOT/analysis/guix/vgp_10_pilot/realization.json"
export VGP_RESOURCE_PLAN="$ROOT/analysis/vgp_three_pair_resource_plan_v1.json"
export VGP_SELECTION_IDS="P02 P03" VGP_NODE_LOCAL_BASE=/scratch
export VGP_MAPPING_BACKEND=wfmash VGP_TASK_ID=run-vgp-three-pair

# P02's mapping completed and passed exact multiplicity before the graph reader
# incident.  Reuse only that closed stage; every graph/index/partition and all
# biological products are rebuilt in this run.
source_mapping="$VGP_DATA_ROOT/pilot/runs/vgp10-auth-20260718-v2-pilot-v1/P02/mapping"
pair_root="$VGP_DATA_ROOT/pilot/runs/$VGP_RUN_ID/P02"
target_mapping="$pair_root/mapping"
if [[ ! -f $target_mapping/.complete.json ]]; then
    [[ ! -e $target_mapping ]] || { echo "incomplete retained mapping target" >&2; exit 2; }
    python3 - "$source_mapping" <<'PY'
import json,sys
from pathlib import Path
from analysis.vgp_10_pilot import sha256_file
root=Path(sys.argv[1]); sentinel=json.loads((root/".complete.json").read_text())
for name,digest in sentinel["files"].items():
    path=root/name
    if not path.is_file() or sha256_file(path)!=digest: raise SystemExit(f"retained mapping digest mismatch: {name}")
if sentinel["selection_id"]!="P02" or sentinel["stage"]!="mapping": raise SystemExit("retained mapping identity mismatch")
PY
    if [[ $mode == --submit ]]; then
        mkdir -p "$pair_root"
        cp -a "$source_mapping" "$target_mapping"
        python3 - "$pair_root/retained_mapping_provenance.json" "$source_mapping" "$target_mapping" <<'PY'
import json,sys
from pathlib import Path
from analysis.vgp_10_pilot import sha256_file
Path(sys.argv[1]).write_text(json.dumps({
 "schema_version":"vgp-three-pair-retained-mapping-v1","selection_id":"P02",
 "source":sys.argv[2],"target":sys.argv[3],"source_job_id":"1782140",
 "reuse_scope":"closed exact SweepGA mapping only","graph_or_variant_products_reused":False,
 "exact_paf_sha256":sha256_file(Path(sys.argv[3])/"h2_to_h1.1to1.paf"),
},sort_keys=True)+"\n")
PY
    else
        echo "would copy validated P02 mapping: $source_mapping -> $target_mapping"
    fi
fi

if [[ $mode == --submit ]]; then
    mkdir -p "$VGP_DATA_ROOT/pilot/logs"
    overlap_job=$(sbatch --parsable \
        --output="$VGP_DATA_ROOT/pilot/logs/$VGP_RUN_ID-P03-overlap-%j.out" \
        --error="$VGP_DATA_ROOT/pilot/logs/$VGP_RUN_ID-P03-overlap-%j.err" \
        "$ROOT/analysis/slurm/run_vgp_three_pair_overlap.sh")
    printf 'controlled_overlap_job\t%s\n' "$overlap_job"
fi
exec "$ROOT/analysis/slurm/vgp_10_pilot/submit.sh" "$mode"
