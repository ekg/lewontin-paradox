#!/usr/bin/env python3
"""Validate the integrated VGP program without launching work or moving data."""

from __future__ import annotations

import csv
import hashlib
import json
import re
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = ROOT / "analysis"
INTEGRATION = ANALYSIS / "vgp_comprehensive_integration_manifest.json"
HISTORY = ANALYSIS / "vgp_historical_run_registry.json"
UTC_RE = re.compile(r"^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z$")
SHA_RE = re.compile(r"^[0-9a-f]{40}$")
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")


class IntegrationError(RuntimeError):
    """The promoted VGP tree violates the frozen integration contract."""


def load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise IntegrationError(f"JSON object required: {path}")
    return value


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def verify_integration_manifest(path: Path = INTEGRATION) -> dict[str, int]:
    manifest = load_json(path)
    if manifest.get("schema_version") != "vgp-comprehensive-integration-manifest-v1":
        raise IntegrationError("unexpected integration manifest schema")
    entries = manifest.get("integration_order")
    if not isinstance(entries, list) or len(entries) != 13:
        raise IntegrationError("exactly 13 promoted task entries are required")
    if [entry.get("order") for entry in entries] != list(range(1, 14)):
        raise IntegrationError("integration order must be contiguous and explicit")
    tasks = [entry.get("task") for entry in entries]
    if len(set(tasks)) != len(tasks):
        raise IntegrationError("integration tasks must be unique")
    position = {task: index for index, task in enumerate(tasks)}
    external = {"integrate-repaired-vgp"}
    for entry in entries:
        for key in ("integrated_commit", "integrated_parent", "patch_id"):
            if not SHA_RE.fullmatch(str(entry.get(key, ""))):
                raise IntegrationError(f"invalid {key} for {entry.get('task')}")
        commits = entry.get("source_commits")
        if not isinstance(commits, list) or not commits or not all(
            SHA_RE.fullmatch(str(commit)) for commit in commits
        ):
            raise IntegrationError(f"invalid source commits for {entry.get('task')}")
        for dependency in entry.get("dependencies", []):
            if dependency not in external and position.get(dependency, 10**9) >= position[entry["task"]]:
                raise IntegrationError(
                    f"dependency order violation: {entry['task']} before {dependency}"
                )
    if any(
        entries[index]["integrated_parent"] != entries[index - 1]["integrated_commit"]
        for index in range(1, len(entries))
    ):
        raise IntegrationError("promoted commits are not the recorded linear chain")
    source_count = sum(len(entry["source_commits"]) for entry in entries)
    if source_count != manifest.get("source_commit_count") or source_count != 17:
        raise IntegrationError("exactly 17 audited source commits are required")
    return {"tasks": len(entries), "source_commits": source_count}


def verify_historical_registry(path: Path = HISTORY, root: Path = ROOT) -> dict[str, int]:
    registry = load_json(path)
    if registry.get("schema_version") != "vgp-historical-run-registry-v1":
        raise IntegrationError("unexpected historical registry schema")
    policy = registry.get("authorization_policy", {})
    if policy.get("current_execution_authorization") is not None:
        raise IntegrationError("historical output must not be current execution authorization")
    if policy.get("historical_artifacts_are_authorization") is not False:
        raise IntegrationError("historical artifacts must be non-authorizing")
    if policy.get("historical_decision_inheritance_forbidden") is not True:
        raise IntegrationError("historical decision inheritance must be forbidden")
    runs = registry.get("runs")
    if not isinstance(runs, list) or len(runs) != 5:
        raise IntegrationError("five timestamped historical evidence packets are required")
    run_ids = [run.get("run_id") for run in runs]
    if len(set(run_ids)) != len(run_ids):
        raise IntegrationError("historical run IDs must be unique")
    artifacts = 0
    for run in runs:
        if run.get("authorization_binding") is not False:
            raise IntegrationError(f"historical run is authorization-binding: {run.get('run_id')}")
        if not str(run.get("classification", "")).startswith("historical_evidence_"):
            raise IntegrationError(f"historical classification absent: {run.get('run_id')}")
        if not UTC_RE.fullmatch(str(run.get("recorded_at_utc", ""))):
            raise IntegrationError(f"timestamp absent or invalid: {run.get('run_id')}")
        if not SHA_RE.fullmatch(str(run.get("source_integration_commit", ""))):
            raise IntegrationError(f"source commit absent or invalid: {run.get('run_id')}")
        for artifact in run.get("artifacts", []):
            relative = Path(str(artifact.get("path", "")))
            if relative.is_absolute() or ".." in relative.parts:
                raise IntegrationError(f"unsafe historical artifact path: {relative}")
            digest = str(artifact.get("sha256", ""))
            if not SHA256_RE.fullmatch(digest):
                raise IntegrationError(f"invalid artifact digest: {relative}")
            absolute = root / relative
            if not absolute.is_file() or sha256_file(absolute) != digest:
                raise IntegrationError(f"historical artifact drift: {relative}")
            artifacts += 1
    return {"historical_runs": len(runs), "historical_artifacts": artifacts}


def verify_rosters_and_tooling(root: Path = ROOT) -> dict[str, int]:
    primaries = read_tsv(root / "analysis/vgp_10_pair_manifest.tsv")
    alternates = read_tsv(root / "analysis/vgp_10_pair_alternates.tsv")
    if [row.get("selection_id") for row in primaries] != [f"P{index:02d}" for index in range(1, 11)]:
        raise IntegrationError("primary roster must be exactly P01 through P10")
    if [row.get("selection_id") for row in alternates] != [f"A{index:02d}" for index in range(1, 7)]:
        raise IntegrationError("alternate roster must be exactly A01 through A06")

    required = {
        "SweepGA": ("analysis/slurm/vgp_10_pilot/pair_stage.sh", "--num-mappings 1:1"),
        "IMPG index": ("analysis/slurm/vgp_10_pilot/pair_stage.sh", '"$impg" index'),
        "IMPG partition": ("analysis/slurm/vgp_10_pilot/pair_stage.sh", '"$impg" partition'),
        "IMPG query": ("analysis/slurm/vgp_10_pilot/impg_parallel_query.sh", '"$impg" query'),
        "IMPG lace": ("analysis/slurm/vgp_10_pilot/impg_hierarchical_lace.sh", '"$impg" lace'),
        "PSMC": ("analysis/slurm/vgp_10_pilot/psmc_array.sh", "verify_tool psmc"),
        "Slurm": ("analysis/slurm/vgp_10_pilot/submit.sh", "sbatch --parsable"),
        "mirror": ("analysis/mirror_vgp_freeze1.py", "def run_worker("),
        "atomic mirror promotion": ("analysis/mirror_vgp_freeze1.py", "def promote_verified("),
        "validation": ("analysis/validate_vgp_10_pilot_acquisition.py", "def validate("),
        "Guix channels": ("analysis/guix/vgp_10_pilot/channels.scm", "(commit \"44bbfc24"),
        "Guix manifest": ("analysis/guix/vgp_10_pilot/manifest.scm", "packages->manifest"),
    }
    for label, (relative, token) in required.items():
        path = root / relative
        if not path.is_file() or token not in path.read_text(encoding="utf-8"):
            raise IntegrationError(f"required {label} tooling absent: {relative}")
    return {
        "primary_pairs": len(primaries),
        "alternate_pairs": len(alternates),
        "tool_checks": len(required),
    }


def validate(root: Path = ROOT) -> dict[str, int]:
    result = {}
    result.update(verify_integration_manifest(root / "analysis/vgp_comprehensive_integration_manifest.json"))
    result.update(verify_historical_registry(root / "analysis/vgp_historical_run_registry.json", root))
    result.update(verify_rosters_and_tooling(root))
    return result


def main() -> int:
    print(json.dumps(validate(), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
