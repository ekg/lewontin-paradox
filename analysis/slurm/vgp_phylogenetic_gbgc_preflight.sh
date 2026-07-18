#!/usr/bin/env bash
#SBATCH --job-name=vgp-gbgc-preflight
#SBATCH --partition=lowmem
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:15:00
#SBATCH --output=analysis/vgp_phylogenetic_gbgc_preflight-%j.out

set -euo pipefail

project_root=${SLURM_SUBMIT_DIR:?Slurm submission directory is required}
cd "$project_root"
test -f analysis/guix/channels.scm
test -f analysis/run_vgp_phylogenetic_gbgc_pilot.py

channel_commit=44bbfc24e4bcc48d0e3343cd3d83452721af8c36
test -n "${SLURM_JOB_ID:-}"
profile=$(readlink -f "${VGP_GBGC_PROFILE:?realized Guix profile is required}")
case "$profile" in
  /gnu/store/????????????????????????????????-*) ;;
  *) echo "profile is not an immutable Guix store path: $profile" >&2; exit 65 ;;
esac
test -x "$profile/bin/python3"

job_home=${TMPDIR:-/tmp}/vgp-gbgc-preflight-${SLURM_JOB_ID}
install -d -m 700 "$job_home"

exec env -i \
  "PATH=$profile/bin" \
  "GUIX_PROFILE=$profile" \
  "GUIX_ENVIRONMENT=$profile" \
  "HOME=$job_home" \
  "TMPDIR=$job_home" \
  "LANG=C" \
  "LC_ALL=C" \
  "TZ=UTC" \
  "PYTHONUTF8=1" \
  "PYTHONPATH=$project_root" \
  "SLURM_JOB_ID=$SLURM_JOB_ID" \
  "VGP_GBGC_PINNED_GUIX=$channel_commit" \
  "$profile/bin/bash" -c \
  'python3 analysis/run_vgp_phylogenetic_gbgc_pilot.py --write --validate && pytest -q analysis/tests/test_vgp_phylogenetic_gbgc_pilot.py analysis/tests/test_gene_conversion_evidence_design.py'
