#!/usr/bin/env python3
"""Validate P07 re-alignment and emit its immutable contract packet."""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import subprocess
from pathlib import Path


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


p = argparse.ArgumentParser()
p.add_argument("--partial", type=Path, required=True)
p.add_argument("--input-manifest", type=Path, required=True)
p.add_argument("--prior-mapping", type=Path, required=True)
p.add_argument("--scratch", required=True)
p.add_argument("--job-id", required=True)
p.add_argument("--started", required=True)
p.add_argument("--capture-sha256", required=True)
p.add_argument("--amendment-sha256", required=True)
args = p.parse_args()
if not args.scratch.startswith("/scratch/vgp-scale-fastga-P07-"):
    raise SystemExit("invalid scratch root")
guard = json.loads((args.partial / "fastga_scratch_contract.json").read_text())
if not guard.get("contract_valid"):
    raise SystemExit("FastGA live scratch contract did not pass")
multiplicity = json.loads((args.partial / "multiplicity.json").read_text())
if (
    multiplicity.get("maximum_query_overlap_depth") != 1
    or multiplicity.get("maximum_target_overlap_depth") != 1
):
    raise SystemExit("absolute one-to-one multiplicity did not pass")
comparisons = {}
for name in ("h2_to_h1.native.1to1.paf", "h2_to_h1.1to1.paf"):
    observed = sha(args.partial / name)
    expected = sha(args.prior_mapping / name)
    comparisons[name] = {"revalidated_sha256": observed, "frozen_sha256": expected,
                         "exact_match": observed == expected}
    if observed != expected:
        raise SystemExit(f"revalidated mapping differs from frozen mapping: {name}")
manifest = json.loads(args.input_manifest.read_text())
repository_commit = subprocess.check_output(
    ["git", "-C", str(Path(__file__).resolve().parents[3]), "rev-parse", "HEAD"], text=True
).strip()
packet = {
    "schema_version": "vgp-real-scaleout-fastga-scratch-contract-v1",
    "canonical_vgp_root": "/moosefs/erikg/vgp",
    "selection_id": "P07",
    "job_id": args.job_id,
    "started_at_utc": args.started,
    "completed_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
    "private_working_directory": args.scratch,
    "resolved_node_local_scratch_root": guard["resolved_node_local_scratch_root"],
    "tmpdir": args.scratch,
    "fastga_live_contract_valid": True,
    "fastga_guard_snapshot_count": guard["snapshot_count"],
    "fastga_process_snapshot_count": guard["fastga_snapshot_count"],
    "observed_fastga_cwds": guard["observed_cwds"],
    "observed_fastga_temp_environments": guard["observed_temp_environments"],
    "observed_fastga_managed_paths": guard["observed_managed_open_paths"],
    "multiplicity_contract_valid": True,
    "command_contract": {"orientation": "H2_query_to_H1_reference", "num_mappings": "1:1",
                         "scaffold_jump": 0, "overlap": 0, "scoring": "log-length-ani"},
    "input_manifest": {"path": str(args.input_manifest), "sha256": sha(args.input_manifest),
                       "inputs": manifest["inputs"]},
    "environment_capture_sha256": args.capture_sha256,
    "fastga_amendment_sha256": args.amendment_sha256,
    "execution_repository_commit": repository_commit,
    "finalizer_sha256": sha(Path(__file__)),
    "frozen_mapping_exact_reproduction": comparisons,
    "downstream_binding": "same frozen P07 mapping; pi and PSMC remain a dependent pair",
    "contract_valid": True,
}
(args.partial / "contract.json").write_text(json.dumps(packet, indent=2, sort_keys=True) + "\n")
