#!/usr/bin/env bash
#SBATCH --job-name=direct-gc-preflight
#SBATCH --partition=lowmem
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=analysis/direct_gene_conversion_preflight-%j.out

set -euo pipefail

project_root=${SLURM_SUBMIT_DIR:?Slurm submission directory is required}
cd "$project_root"
test -f analysis/guix/channels.scm
test -f analysis/run_direct_gene_conversion_pilot.py

channel_commit=44bbfc24e4bcc48d0e3343cd3d83452721af8c36
test -n "${SLURM_JOB_ID:-}"
profile=$(readlink -f "${DIRECT_GC_PROFILE:?realized Guix profile is required}")
expected_profile=/gnu/store/3c2mxm30rbzvnw7qsi235mrkk3m38fym-profile
test "$profile" = "$expected_profile"
case "$profile" in
  /gnu/store/????????????????????????????????-*) ;;
  *) echo "profile is not an immutable Guix store path: $profile" >&2; exit 65 ;;
esac
test -x "$profile/bin/python3"
test -x "$profile/bin/pytest"

job_home=${TMPDIR:-/tmp}/direct-gc-preflight-${SLURM_JOB_ID}
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
  "DIRECT_GC_PINNED_GUIX=$channel_commit" \
  "$profile/bin/bash" -c \
  'python3 analysis/run_direct_gene_conversion_pilot.py --write --validate && pytest -q analysis/tests/test_direct_gene_conversion_pilot.py analysis/tests/test_gene_conversion_evidence_design.py'
