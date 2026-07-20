#!/usr/bin/env python3
"""Fail closed when a live FastGA process escapes its private node-local tree."""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
import re
from pathlib import Path
from typing import Iterable, Mapping


MANAGED_NAME = re.compile(
    r"(?:_tmp_|(?:^|[._-])pair(?:[._-]|$)|_algn\.|_uniq\.|"
    r"\.las(?:\s+\(deleted\))?$|\.1aln$|\.bps$|\.gix$|\.index$|\.fa(?:\.|$))"
)


class ScratchContractError(RuntimeError):
    """A FastGA working path violated the node-local scratch contract."""


def _resolved(path: os.PathLike[str] | str) -> Path:
    text = str(path)
    if text.endswith(" (deleted)"):
        text = text[: -len(" (deleted)")]
    return Path(text).resolve(strict=False)


def _require_under(path: os.PathLike[str] | str, root: Path, label: str) -> str:
    resolved = _resolved(path)
    try:
        resolved.relative_to(root)
    except ValueError as error:
        raise ScratchContractError(
            f"{label} outside private node-local scratch: {resolved} (root {root})"
        ) from error
    return str(resolved)


def validate_snapshot(
    *,
    scratch: os.PathLike[str] | str,
    pid: int,
    cwd: os.PathLike[str] | str,
    temp_env: Mapping[str, str],
    managed_open_paths: Iterable[os.PathLike[str] | str],
) -> dict[str, object]:
    """Validate one FastGA /proc snapshot and return its auditable form."""
    root = _resolved(scratch)
    cwd_resolved = _require_under(cwd, root, f"FastGA pid {pid} cwd")
    temps: dict[str, str] = {}
    for name in ("TMPDIR", "TMP", "TEMP"):
        value = temp_env.get(name)
        if not value:
            raise ScratchContractError(
                f"FastGA pid {pid} {name} outside private node-local scratch: unset"
            )
        temps[name] = _require_under(value, root, f"FastGA pid {pid} {name}")
    paths = sorted(
        {
            _require_under(path, root, f"FastGA pid {pid} managed open file")
            for path in managed_open_paths
        }
    )
    return {
        "pid": pid,
        "contract_valid": True,
        "scratch_resolved": str(root),
        "cwd_resolved": cwd_resolved,
        "temp_environment_resolved": temps,
        "managed_open_paths_resolved": paths,
    }


def _process_table() -> dict[int, int]:
    table: dict[int, int] = {}
    for entry in Path("/proc").iterdir():
        if not entry.name.isdigit():
            continue
        try:
            fields = (entry / "stat").read_text().split()
            table[int(entry.name)] = int(fields[3])
        except (FileNotFoundError, IndexError, PermissionError, ValueError):
            continue
    return table


def _descendants(parent: int) -> list[int]:
    table = _process_table()
    found: list[int] = []
    frontier = [parent]
    while frontier:
        current = frontier.pop()
        children = [pid for pid, ppid in table.items() if ppid == current]
        found.extend(children)
        frontier.extend(children)
    return found


def _environment(pid: int) -> dict[str, str]:
    raw = Path(f"/proc/{pid}/environ").read_bytes()
    values: dict[str, str] = {}
    for field in raw.split(b"\0"):
        if b"=" in field:
            name, value = field.split(b"=", 1)
            values[name.decode(errors="replace")] = value.decode(errors="replace")
    return values


def _managed_fds(pid: int) -> list[str]:
    values: list[str] = []
    for fd in Path(f"/proc/{pid}/fd").iterdir():
        try:
            target = os.readlink(fd)
        except (FileNotFoundError, PermissionError):
            continue
        if target.startswith("/") and MANAGED_NAME.search(Path(target).name):
            values.append(target)
    return values


def inspect_fastga_descendants(parent_pid: int, scratch: Path) -> dict[str, object]:
    snapshots: list[dict[str, object]] = []
    for pid in _descendants(parent_pid):
        try:
            if Path(f"/proc/{pid}/comm").read_text().strip() != "FastGA":
                continue
            snapshots.append(
                validate_snapshot(
                    scratch=scratch,
                    pid=pid,
                    cwd=os.readlink(f"/proc/{pid}/cwd"),
                    temp_env=_environment(pid),
                    managed_open_paths=_managed_fds(pid),
                )
            )
        except FileNotFoundError:
            # A process may exit between the process-table and fd snapshots.
            continue
    return {
        "observed_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "parent_pid": parent_pid,
        "fastga_processes": snapshots,
        "contract_valid": True,
    }


def _last_jsonl(path: Path) -> dict[str, object] | None:
    if not path.is_file() or path.stat().st_size == 0:
        return None
    with path.open("rb") as handle:
        handle.seek(0, os.SEEK_END)
        position = handle.tell() - 1
        if position >= 0:
            handle.seek(position)
            if handle.read(1) == b"\n":
                position -= 1
        while position >= 0:
            handle.seek(position)
            if handle.read(1) == b"\n":
                position += 1
                break
            position -= 1
        handle.seek(max(0, position))
        line = handle.readline().strip()
    return json.loads(line) if line else None


def _path_digest(paths: list[str]) -> str:
    return hashlib.sha256(("\0".join(paths) + "\0").encode()).hexdigest()


def _compact_repeated_paths(path: Path, row: dict[str, object]) -> None:
    previous = _last_jsonl(path)
    previous_by_pid = {
        snapshot["pid"]: snapshot
        for snapshot in (previous or {}).get("fastga_processes", [])
    }
    for snapshot in row.get("fastga_processes", []):
        paths = snapshot.get("managed_open_paths_resolved", [])
        digest = _path_digest(paths)
        snapshot["managed_open_path_count"] = len(paths)
        snapshot["managed_open_paths_sha256"] = digest
        prior = previous_by_pid.get(snapshot["pid"], {})
        prior_paths = prior.get("managed_open_paths_resolved", [])
        prior_digest = prior.get("managed_open_paths_sha256")
        if prior_digest is None and prior_paths:
            prior_digest = _path_digest(prior_paths)
        if prior_digest == digest:
            snapshot.pop("managed_open_paths_resolved", None)
            snapshot["managed_open_paths_unchanged"] = True


def _append_jsonl(path: Path, row: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    _compact_repeated_paths(path, row)
    with path.open("a") as handle:
        handle.write(json.dumps(row, sort_keys=True) + "\n")


def _check(args: argparse.Namespace) -> int:
    try:
        row = inspect_fastga_descendants(args.parent_pid, args.scratch)
    except ScratchContractError as error:
        row = {
            "observed_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
            "parent_pid": args.parent_pid,
            "contract_valid": False,
            "violation": str(error),
        }
        _append_jsonl(args.audit_jsonl, row)
        print(str(error), file=os.sys.stderr)
        return 70
    _append_jsonl(args.audit_jsonl, row)
    return 0


def _finalize(args: argparse.Namespace) -> int:
    rows = [json.loads(line) for line in args.audit_jsonl.read_text().splitlines() if line]
    violations = [row for row in rows if not row.get("contract_valid")]
    snapshots = [snapshot for row in rows for snapshot in row.get("fastga_processes", [])]
    if violations:
        raise ScratchContractError(f"live FastGA path violations recorded: {violations}")
    if not snapshots:
        raise ScratchContractError("no live FastGA /proc snapshot was observed")
    managed = sorted(
        {path for row in snapshots for path in row.get("managed_open_paths_resolved", [])}
    )
    result = {
        "schema_version": "vgp-fastga-node-local-scratch-contract-v1",
        "contract_valid": True,
        "canonical_requested_scratch_root": str(args.scratch),
        "resolved_node_local_scratch_root": str(_resolved(args.scratch)),
        "snapshot_count": len(rows),
        "fastga_snapshot_count": len(snapshots),
        "fastga_pids": sorted({row["pid"] for row in snapshots}),
        "observed_cwds": sorted({row["cwd_resolved"] for row in snapshots}),
        "observed_temp_environments": sorted(
            {json.dumps(row["temp_environment_resolved"], sort_keys=True) for row in snapshots}
        ),
        "observed_managed_open_paths": managed,
        "validated_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
    }
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)
    check = subparsers.add_parser("check")
    check.add_argument("--parent-pid", type=int, required=True)
    check.add_argument("--scratch", type=Path, required=True)
    check.add_argument("--audit-jsonl", type=Path, required=True)
    check.set_defaults(function=_check)
    finalize = subparsers.add_parser("finalize")
    finalize.add_argument("--scratch", type=Path, required=True)
    finalize.add_argument("--audit-jsonl", type=Path, required=True)
    finalize.add_argument("--output", type=Path, required=True)
    finalize.set_defaults(function=_finalize)
    args = parser.parse_args()
    try:
        return args.function(args)
    except ScratchContractError as error:
        parser.error(str(error))


if __name__ == "__main__":
    raise SystemExit(main())
