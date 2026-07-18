#!/usr/bin/env bash
#SBATCH --job-name=vgp-gbgc-tests
#SBATCH --partition=lowmem
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=analysis/vgp_phylogenetic_gbgc_full_tests-%j.out

set -euo pipefail

project_root=${SLURM_SUBMIT_DIR:?Slurm submission directory is required}
cd "$project_root"
profile=$(readlink -f "${VGP_GBGC_FULL_TEST_PROFILE:?realized full-test Guix profile is required}")
case "$profile" in
  /gnu/store/????????????????????????????????-*) ;;
  *) echo "profile is not an immutable Guix store path: $profile" >&2; exit 65 ;;
esac
test -x "$profile/bin/pytest"

job_home=${TMPDIR:-/tmp}/vgp-gbgc-full-tests-${SLURM_JOB_ID:?}
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
  "PYTHONPATH=$project_root:$profile/lib/python3.10/site-packages" \
  "$profile/bin/pytest" -q analysis/tests
