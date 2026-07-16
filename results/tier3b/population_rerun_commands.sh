#!/usr/bin/env bash
# Exact Tier 3B recovered-population rerun recipe.
# Scheduler history: application-failed array 1761582 / canceled dependency
# 1761583; performance-profile array 1761585 / canceled dependency 1761586;
# final pre-repair 1-Mb evidence array 1761588 / independent check 1761589;
# powered biological array 1761721 / independent check 1761722.
set -euo pipefail

ROOT=$(cd "$(dirname "$0")/../.." && pwd)
MANIFEST="$ROOT/results/tier3b/acquisition_manifest.tsv"
ENVIRONMENT_RECORD="$ROOT/analysis/pilot_results/guix_environment.json"
RAW="$ROOT/results/tier3b/population_raw"
SUBSET="$ROOT/results/tier3b/population_subset"
LOGS="$ROOT/results/tier3b/run_logs"
SCRATCH=/moosefs/erikg/tier3scratch/tier3b-recovery
export TIER3_ROOT="$ROOT" TIER3_ENVIRONMENT_RECORD="$ENVIRONMENT_RECORD" TIER3_SCRATCH_ROOT="$SCRATCH"

mkdir -p "$RAW" "$SUBSET" "$LOGS"

# Validate scheduler resource requests, then serialize the two full MooseFS
# scans with an array %1 throttle.  The 4-hour array request is based on the
# observed 8-minute scan time per 1 Mb and the repaired 21-Mb inputs, with a
# modest margin.  The independent implementation starts only after both
# biological tuple tasks succeed.
sbatch --test-only "$ROOT/analysis/slurm/tier3b_population_recovery_array.sh" \
  "$ROOT" "$MANIFEST" "$RAW" "$ENVIRONMENT_RECORD" "$SCRATCH"
sbatch --test-only "$ROOT/analysis/slurm/tier3b_population_independent_check.sh" \
  "$ROOT" "$MANIFEST" "$SUBSET" "$ENVIRONMENT_RECORD" "$SCRATCH"
ARRAY_JOB=$(sbatch --parsable "$ROOT/analysis/slurm/tier3b_population_recovery_array.sh" \
  "$ROOT" "$MANIFEST" "$RAW" "$ENVIRONMENT_RECORD" "$SCRATCH")
CHECK_JOB=$(sbatch --parsable --dependency="afterok:$ARRAY_JOB" \
  "$ROOT/analysis/slurm/tier3b_population_independent_check.sh" \
  "$ROOT" "$MANIFEST" "$SUBSET" "$ENVIRONMENT_RECORD" "$SCRATCH")
printf 'population array job: %s\nindependent check job: %s\n' "$ARRAY_JOB" "$CHECK_JOB"

# Wait without polling the shared inputs; only scheduler state is queried.
while squeue -h -j "$ARRAY_JOB,$CHECK_JOB" | grep -q .; do
  sleep 30
done

# After Slurm reports both jobs complete, capture sacct into the same columns
# used by the original run.  Use the array task job rows for state/exit and the
# corresponding .batch rows for peak RSS when available.
sacct -j "$ARRAY_JOB,$CHECK_JOB" -P \
  --format=JobID,JobIDRaw,JobName,State,Reason,Elapsed,TotalCPU,MaxRSS,ReqMem,AllocCPUS,ExitCode \
  > "$LOGS/population_sacct_raw.psv"

# Tuple IDs map to array indices in acquisition_manifest.tsv approved-row
# order; conversion rejects missing, failed, or nonzero-exit tasks.
"$ROOT/analysis/slurm/guix_job.sh" "$ENVIRONMENT_RECORD" \
  python3 "$ROOT/analysis/tier3b_population_recovery.py" telemetry \
  --manifest "$MANIFEST" --sacct "$LOGS/population_sacct_raw.psv" \
  --array-job "$ARRAY_JOB" --output "$LOGS/population_scheduler_telemetry.tsv"
"$ROOT/analysis/slurm/guix_job.sh" "$ENVIRONMENT_RECORD" \
  python3 "$ROOT/analysis/tier3b_population_recovery.py" collect \
  --manifest "$MANIFEST" --raw-dir "$RAW" --output-dir "$ROOT/results/tier3b" \
  --telemetry "$LOGS/population_scheduler_telemetry.tsv" \
  --environment "$ENVIRONMENT_RECORD"
"$ROOT/analysis/slurm/guix_job.sh" "$ENVIRONMENT_RECORD" \
  python3 "$ROOT/analysis/tier3b_population_recovery.py" validate \
  --manifest "$MANIFEST" --diversity "$ROOT/results/tier3b/population_diversity.tsv" \
  --independent "$ROOT/results/tier3b/population_independent_check.tsv" \
  2>&1 | tee "$LOGS/population_final_assertion.log"
"$ROOT/analysis/slurm/guix_job.sh" "$ENVIRONMENT_RECORD" \
  python3 -m pytest -q "$ROOT/analysis/tests" \
  2>&1 | tee "$LOGS/population_guix_tests.log"

{
  "$ROOT/analysis/slurm/guix_job.sh" "$ENVIRONMENT_RECORD" python3 --version
  "$ROOT/analysis/slurm/guix_job.sh" "$ENVIRONMENT_RECORD" bcftools --version
  "$ROOT/analysis/slurm/guix_job.sh" "$ENVIRONMENT_RECORD" samtools --version
  "$ROOT/analysis/slurm/guix_job.sh" "$ENVIRONMENT_RECORD" bgzip --version
  "$ROOT/analysis/slurm/guix_job.sh" "$ENVIRONMENT_RECORD" bedtools --version
  "$ROOT/analysis/slurm/guix_job.sh" "$ENVIRONMENT_RECORD" \
    python3 -m pytest -q "$ROOT/analysis/tests/test_tier3b.py" \
    "$ROOT/analysis/tests/test_tier3b_recovery.py"
} 2>&1 | tee "$LOGS/population_guix_smoke.log"
