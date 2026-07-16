#!/usr/bin/env bash
set -euo pipefail

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
cd "$ROOT"
mkdir -p results/tier3a/logs

# Primary biological run: the three acquisition rows are array indices 0-2.
# Each task consumes the checksum-locked whole-H1/H2 SweepGA --num-mappings
# 1:1 PAF, then lets IMPG independently index, partition, and query it.
sbatch --array=0-2 analysis/slurm/tier3a_biological_array.sh

# After the array finishes, replace JOB_ID and capture measured Slurm resources.
JOB_ID=${1:-JOB_ID}
sacct -j "$JOB_ID" --allocations -P -n \
    -o JobIDRaw,JobName,State,ExitCode,Elapsed,MaxRSS,ReqMem,AllocCPUS,NodeList \
    > "results/tier3a/logs/tier3a-${JOB_ID}.sacct"

# The committed aggregation command is rerunnable after telemetry.tsv has been
# regenerated from sacct (the run manifest preserves the exact completed paths).
guix time-machine -C analysis/guix/channels.scm -- \
    shell -m analysis/guix/manifest.scm --pure -- \
    python3 analysis/tier3a_biological.py finalize \
        --manifest results/tier3a/acquisition_corrected_manifest.tsv \
        --work-root /moosefs/erikg/tier3data/tier3a-origin-biological-20260716 \
        --telemetry results/tier3a/logs/telemetry.tsv \
        --output-dir results/tier3a
