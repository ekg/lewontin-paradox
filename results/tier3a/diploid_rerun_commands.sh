#!/usr/bin/env bash
set -euo pipefail

TIER3A_WORK_ROOT=/moosefs/erikg/tier3data/tier3a-origin-biological-rerun-20260716 sbatch --array=0 analysis/slurm/tier3a_biological_array.sh # menidia_menidia_fMenMen1
TIER3A_WORK_ROOT=/moosefs/erikg/tier3data/tier3a-origin-biological-rerun-20260716 sbatch --array=1 analysis/slurm/tier3a_biological_array.sh # spinachia_spinachia_SK-2024b
TIER3A_WORK_ROOT=/moosefs/erikg/tier3data/tier3a-origin-biological-rerun-20260716 sbatch --array=2 analysis/slurm/tier3a_biological_array.sh # tautogolabrus_adspersus_fTauAds1

# After the array completes, regenerate sacct telemetry and run:
python3 analysis/tier3a_biological.py finalize --manifest results/tier3a/acquisition_corrected_manifest.tsv --work-root /moosefs/erikg/tier3data/tier3a-origin-biological-rerun-20260716 --telemetry results/tier3a/logs/origin-rerun-telemetry.tsv --supersession-ledger results/tier3a/acquisition_sweepga_supersession.tsv --superseded-results results/tier3a/diploid_superseded_results.tsv --output-dir results/tier3a
