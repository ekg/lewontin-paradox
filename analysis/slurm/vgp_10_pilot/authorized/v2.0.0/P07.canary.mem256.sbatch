#!/usr/bin/env bash
#SBATCH --job-name=vgp10-P07-r256
#SBATCH --partition=highmem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=3-00:00:00
#SBATCH --signal=B:USR1@300
#SBATCH --output=/moosefs/erikg/vgp/pilot/logs/vgp10-auth-v2-P07-r256-%j.out
#SBATCH --error=/moosefs/erikg/vgp/pilot/logs/vgp10-auth-v2-P07-r256-%j.err
set -euo pipefail
export VGP_RETRY=1
SUBMIT_ROOT="${SLURM_SUBMIT_DIR:?submit from the repository root}"
CANARY_WORKER="$SUBMIT_ROOT/analysis/slurm/vgp_10_pilot/authorized/v2.0.0/run_canary.sh"
[[ -f $CANARY_WORKER ]] || { echo "canary worker absent: $CANARY_WORKER" >&2; exit 2; }
exec bash "$CANARY_WORKER"
