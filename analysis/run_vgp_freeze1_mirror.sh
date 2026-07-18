#!/usr/bin/env bash
set -euo pipefail

project_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$project_root"

exec guix time-machine -C analysis/guix/channels.scm -- \
  shell -m analysis/guix/vgp-freeze1-manifest.scm --pure -- \
  python3 analysis/mirror_vgp_freeze1.py "$@"
