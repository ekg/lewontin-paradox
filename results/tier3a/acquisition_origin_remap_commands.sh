#!/usr/bin/env bash
# Submit and atomically publish the Tier 3A origin/main SweepGA correction.
set -euo pipefail

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
cd "$ROOT"
mkdir -p results/tier3a/logs

if [[ ${1:-} != finalize:* ]]; then
    JOB_ID=$(sbatch --parsable analysis/slurm/tier3a_origin_remap_array.sh)
    printf '%s\n' "$JOB_ID" > results/tier3a/logs/origin-remap.job-id
    printf 'submitted Tier 3A origin/main remap array %s\n' "$JOB_ID"
    printf 'after every array task completes, run:\n'
    printf '  %q %q\n' "$0" "finalize:$JOB_ID"
    exit 0
fi

JOB_ID=${1#finalize:}
if squeue -h -j "$JOB_ID" | grep -q .; then
    echo "Slurm array $JOB_ID is still active" >&2
    exit 75
fi
sacct -X -j "$JOB_ID" \
    --format=JobIDRaw,JobID,JobName,State,ExitCode,Elapsed,TotalCPU,MaxRSS,ReqMem,AllocCPUS,NodeList \
    --parsable2 > "results/tier3a/logs/origin-remap-${JOB_ID}.sacct"
if [[ $(awk -F'|' 'NR > 1 && $1 ~ /_[012]$/ && $4 == "COMPLETED" {n++} END {print n+0}' "results/tier3a/logs/origin-remap-${JOB_ID}.sacct") -ne 3 ]]; then
    echo "Not all three Slurm array tasks completed successfully" >&2
    exit 1
fi

python3 analysis/tier3a_origin_remap.py finalize \
    --base-manifest results/tier3a/acquisition_manifest.tsv \
    --work-root /moosefs/erikg/tier3data/tier3a-origin-remap-20260716 \
    --sacct "results/tier3a/logs/origin-remap-${JOB_ID}.sacct" \
    --mapping-job-id "$JOB_ID" \
    --output-manifest results/tier3a/acquisition_corrected_manifest.tsv \
    --supersession-ledger results/tier3a/acquisition_sweepga_supersession.tsv \
    --commands results/tier3a/acquisition_corrected_commands.tsv \
    --qc results/tier3a/acquisition_corrected_qc.md
python3 analysis/tier3a_origin_remap.py validate-manifest \
    --manifest results/tier3a/acquisition_corrected_manifest.tsv
