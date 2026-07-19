#!/usr/bin/env bash
#SBATCH --job-name=vgp10-psmc-finalize
set -euo pipefail

source "$(dirname -- "$0")/common.sh"
VGP_SELECTION_ID=${1:?usage: psmc_finalize.sh SELECTION_ID}
VGP_STAGE_NAME=psmc/finalize
export VGP_SELECTION_ID VGP_STAGE_NAME
require_selection_id "$VGP_SELECTION_ID"
require_runtime
begin_stage "$VGP_STAGE_NAME"

input_dir="$VGP_DATA_ROOT/pilot/inputs/$VGP_SELECTION_ID"
python3 - "$VGP_PAIR_RUN" "$input_dir/psmc_scaling_scenarios.tsv" "$VGP_STAGE_PARTIAL" <<'PY'
import csv,json,math,sys
from pathlib import Path
from analysis.vgp_10_pilot import PilotError,parse_psmc_unscaled,scale_unscaled_trajectory,sha256_file
pair,scenario_path,out=Path(sys.argv[1]),Path(sys.argv[2]),Path(sys.argv[3])
primary=pair/"psmc/replicate-000/unscaled.psmc"
rows,theta=parse_psmc_unscaled(primary)
with (out/"unscaled_trajectory.tsv").open("w",newline="") as handle:
    writer=csv.DictWriter(handle,fieldnames=("interval","time_2N0","lambda"),
                          delimiter="\t",lineterminator="\n")
    writer.writeheader(); writer.writerows(rows)
if scenario_path.is_file():
    with scenario_path.open(newline="") as handle:
        scenarios=list(csv.DictReader(handle,delimiter="\t"))
else:
    # Absolute scaling is optional.  Its absence cannot veto the unscaled PSMC
    # trajectory or manufacture a zero-valued demographic estimate.
    scenarios=[]
for row in scenarios: row["theta_0"]=theta
_,scaled=scale_unscaled_trajectory(rows,scenarios)
fields=("scenario_id","interval","time_years","effective_size","mutation_rate_per_generation",
        "generation_time_years","psmc_bin_size_bp","mutation_rate_source","generation_time_source")
with (out/"scenario_scaled_trajectories.tsv").open("w",newline="") as handle:
    writer=csv.DictWriter(handle,fieldnames=fields,delimiter="\t",lineterminator="\n")
    writer.writeheader(); writer.writerows(scaled)
boot=[]
for replicate in range(1,201):
    path=pair/f"psmc/replicate-{replicate:03d}/bootstrap.unscaled.psmc"
    if not path.is_file(): raise SystemExit(f"missing bootstrap output: {replicate}")
    try:
        values,_=parse_psmc_unscaled(path)
        finite=all(math.isfinite(row["time_2N0"]) and math.isfinite(row["lambda"]) for row in values)
    except PilotError:
        finite=False
    boot.append({"replicate":replicate,"finite":str(finite).lower(),"sha256":sha256_file(path)})
with (out/"bootstrap_unscaled.tsv").open("w",newline="") as handle:
    writer=csv.DictWriter(handle,fieldnames=("replicate","finite","sha256"),
                          delimiter="\t",lineterminator="\n")
    writer.writeheader(); writer.writerows(boot)
finite=sum(row["finite"] == "true" for row in boot)
(out/"psmc_qc.json").write_text(json.dumps({
    "bootstrap_attempts":200,"finite_bootstraps":finite,"minimum_finite_required":190,
    "passed":finite>=190,"unscaled_primary_preserved":True,"scenario_scaling_separate":True,
    "absolute_scaling_status":"available" if scenarios else "not_available_nonblocking",
},sort_keys=True)+"\n")
if finite < 190: raise SystemExit("fewer than 190/200 finite bootstrap trajectories")
PY
promote_stage
