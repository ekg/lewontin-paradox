#!/usr/bin/env bash
#SBATCH --job-name=vgp10-P07-canary
#SBATCH --partition=highmem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=3-00:00:00
#SBATCH --signal=B:USR1@300
#SBATCH --output=/moosefs/erikg/vgp/pilot/logs/vgp10-auth-v2-P07-%j.out
#SBATCH --error=/moosefs/erikg/vgp/pilot/logs/vgp10-auth-v2-P07-%j.err
set -euo pipefail
export VGP_RETRY=0
exec bash "$(dirname -- "$0")/run_canary.sh"
