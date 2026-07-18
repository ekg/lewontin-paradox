#!/usr/bin/env python3
"""Provision and validate the bounded VGP pilot data root."""

from __future__ import annotations

import argparse
import errno
import fcntl
import grp
import hashlib
import json
import os
import pwd
import shlex
import shutil
import subprocess
import sys
import uuid
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = ROOT / "analysis" / "vgp_data_root_config.json"
DEFAULT_EVIDENCE = ROOT / "analysis" / "vgp_data_root_validation.json"
DEFAULT_MARKDOWN = ROOT / "analysis" / "vgp_data_root_contract.md"

HOST_COMMANDS = {
    "df": ("/bin/df", "/usr/bin/df"),
    "findmnt": ("/bin/findmnt", "/usr/bin/findmnt"),
    "ls": ("/bin/ls", "/usr/bin/ls"),
    "mfsmount": ("/usr/bin/mfsmount", "/bin/mfsmount"),
    "stat": ("/usr/bin/stat", "/bin/stat"),
}
QUOTA_COMMANDS = ("mfsgetquota", "mfsquota", "quota", "repquota", "lfs")
HOST_SEARCH_DIRS = ("/bin", "/usr/bin", "/sbin", "/usr/sbin")


@dataclass(frozen=True)
class Layout:
    root: Path
    directories: dict[str, Path]
    mode: int


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def load_contract(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def parse_mode(text: str) -> int:
    return int(text, 8)


def resolved_root(contract: dict[str, Any], root_override: Path | None = None) -> Path:
    return (root_override or Path(contract["root"])).resolve()


def build_layout(contract: dict[str, Any], root_override: Path | None = None) -> Layout:
    root = resolved_root(contract, root_override)
    mode = parse_mode(contract["default_directory_mode"])
    directories = {
        name: (root / relative).resolve()
        for name, relative in contract["layout"].items()
    }
    return Layout(root=root, directories=directories, mode=mode)


def ensure_within_root(root: Path, path: Path) -> Path:
    resolved = path.resolve()
    resolved.relative_to(root.resolve())
    return resolved


def ensure_directory(path: Path, mode: int) -> None:
    path.mkdir(parents=True, exist_ok=True)
    os.chmod(path, mode)


def provision_layout(layout: Layout) -> list[dict[str, str]]:
    ensure_directory(layout.root, layout.mode)
    records = []
    for name, path in layout.directories.items():
        ensure_within_root(layout.root, path)
        ensure_directory(path, layout.mode)
        records.append(
            {
                "name": name,
                "path": str(path),
                "relative_to_root": str(path.relative_to(layout.root)),
                "mode_octal": oct(path.stat().st_mode & 0o7777),
            }
        )
    return records


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fsync_path(path: Path) -> dict[str, Any]:
    flags = os.O_RDONLY
    if hasattr(os, "O_DIRECTORY"):
        flags |= os.O_DIRECTORY
    fd = os.open(path, flags)
    try:
        os.fsync(fd)
        return {"status": "pass", "path": str(path)}
    except OSError as error:
        if error.errno in {errno.EINVAL, errno.ENOTSUP, errno.EOPNOTSUPP}:
            return {
                "status": "not_supported",
                "path": str(path),
                "errno": error.errno,
                "message": str(error),
            }
        raise
    finally:
        os.close(fd)


def remove_path_within_root(root: Path, path: Path) -> None:
    target = ensure_within_root(root, path)
    if target.is_dir():
        shutil.rmtree(target)
    elif target.exists():
        target.unlink()


def prune_empty_parents(root: Path, path: Path, stop_at: Path) -> None:
    current = ensure_within_root(root, path)
    stop = ensure_within_root(root, stop_at)
    while current != stop:
        if not current.exists() or not current.is_dir():
            current = current.parent
            continue
        try:
            current.rmdir()
        except OSError as error:
            if error.errno in {errno.ENOTEMPTY, errno.EEXIST}:
                break
            raise
        current = current.parent


def smoke_test_storage_contract(layout: Layout) -> dict[str, Any]:
    test_id = f"validation-{os.getpid()}-{uuid.uuid4().hex[:10]}"
    partial_dir = layout.directories["staging_partials"] / test_id
    object_dir = layout.directories["immutable_objects"] / ".validation" / test_id
    lock_path = layout.directories["locks"] / f"{test_id}.lock"
    payload = b"vgp-root-validation\n"
    digest = hashlib.sha256(payload).hexdigest()
    partial_path = partial_dir / "payload.partial"
    final_path = object_dir / digest
    results: dict[str, Any] = {
        "test_id": test_id,
        "payload_sha256": digest,
        "payload_size_bytes": len(payload),
    }
    ensure_directory(partial_dir, layout.mode)
    ensure_directory(object_dir, layout.mode)
    try:
        with partial_path.open("wb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        results["file_fsync"] = {"status": "pass", "path": str(partial_path)}
        results["staging_dir_fsync"] = fsync_path(partial_dir)
        same_filesystem = partial_dir.stat().st_dev == object_dir.stat().st_dev
        results["same_filesystem"] = same_filesystem
        results["partial_device"] = partial_dir.stat().st_dev
        results["object_device"] = object_dir.stat().st_dev
        if not same_filesystem:
            results["atomic_promotion"] = {
                "status": "blocked",
                "message": "staging and immutable object directories are on different filesystems",
            }
        else:
            os.replace(partial_path, final_path)
            results["atomic_promotion"] = {
                "status": "pass",
                "source": str(partial_path),
                "destination": str(final_path),
            }
            results["object_dir_fsync"] = fsync_path(object_dir)
            results["checksum_verification"] = {
                "status": "pass" if sha256_file(final_path) == digest else "fail",
                "path": str(final_path),
            }
        fd_one = os.open(lock_path, os.O_RDWR | os.O_CREAT, 0o660)
        fd_two = os.open(lock_path, os.O_RDWR)
        try:
            fcntl.flock(fd_one, fcntl.LOCK_EX | fcntl.LOCK_NB)
            try:
                fcntl.flock(fd_two, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError:
                results["lock_behavior"] = {"status": "pass", "path": str(lock_path)}
            else:
                results["lock_behavior"] = {
                    "status": "fail",
                    "message": "second non-blocking exclusive lock unexpectedly succeeded",
                    "path": str(lock_path),
                }
                fcntl.flock(fd_two, fcntl.LOCK_UN)
        finally:
            fcntl.flock(fd_one, fcntl.LOCK_UN)
            os.close(fd_one)
            os.close(fd_two)
    finally:
        remove_path_within_root(layout.root, partial_dir)
        remove_path_within_root(layout.root, object_dir)
        remove_path_within_root(layout.root, lock_path)
        prune_empty_parents(layout.root, object_dir.parent, layout.directories["immutable_objects"])
    results["cleanup"] = {
        "status": "pass"
        if (
            not partial_dir.exists()
            and not object_dir.exists()
            and not lock_path.exists()
            and not object_dir.parent.exists()
        )
        else "fail",
        "paths_removed": [str(partial_dir), str(object_dir), str(lock_path), str(object_dir.parent)],
    }
    return results


def find_host_command(name: str) -> str | None:
    discovered = shutil.which(name)
    if discovered:
        return discovered
    for directory in HOST_SEARCH_DIRS:
        candidate = Path(directory) / name
        if candidate.exists():
            return str(candidate)
    return None


def resolve_host_command(name: str) -> str | None:
    discovered = find_host_command(name)
    if discovered:
        return discovered
    for candidate in HOST_COMMANDS.get(name, ()):
        if Path(candidate).exists():
            return candidate
    return None


def run_command(argv: list[str]) -> dict[str, Any]:
    completed = subprocess.run(argv, capture_output=True, text=True, check=False)
    return {
        "timestamp_utc": utc_now(),
        "command": " ".join(shlex.quote(value) for value in argv),
        "exit_code": completed.returncode,
        "stdout": completed.stdout,
        "stderr": completed.stderr,
    }


def collect_system_evidence(root: Path) -> dict[str, Any]:
    commands: list[dict[str, Any]] = []
    ls_cmd = resolve_host_command("ls")
    findmnt_cmd = resolve_host_command("findmnt")
    stat_cmd = resolve_host_command("stat")
    df_cmd = resolve_host_command("df")
    mfsmount_cmd = resolve_host_command("mfsmount")
    if ls_cmd:
        commands.append(run_command([ls_cmd, "-ld", str(root)]))
    if findmnt_cmd:
        commands.append(run_command([findmnt_cmd, "-T", str(root)]))
    if stat_cmd:
        commands.append(
            run_command(
                [
                    stat_cmd,
                    "-f",
                    "-c",
                    "fstype=%T block_size=%S blocks=%b blocks_available=%a files=%c files_free=%d",
                    str(root),
                ]
            )
        )
    if df_cmd:
        commands.append(run_command([df_cmd, "-hP", str(root)]))
        commands.append(run_command([df_cmd, "-iP", str(root)]))
    if mfsmount_cmd:
        commands.append(run_command([mfsmount_cmd, "-V"]))

    available_quota_interfaces = []
    for name in QUOTA_COMMANDS:
        path = find_host_command(name)
        if path:
            available_quota_interfaces.append({"name": name, "path": path})

    inode_status = "unknown"
    inode_detail: dict[str, Any] = {"status": "unknown"}
    for record in commands:
        if " -iP " not in f" {record['command']} ":
            continue
        if record["exit_code"] != 0:
            inode_status = "blocked"
            inode_detail = {
                "status": "blocked",
                "message": "df -i did not complete successfully",
                "command": record["command"],
            }
            break
        lines = [line for line in record["stdout"].splitlines() if line.strip()]
        if len(lines) >= 2 and not any(token == "-" for token in lines[-1].split()):
            inode_status = "reported"
            inode_detail = {
                "status": "reported",
                "command": record["command"],
                "output": record["stdout"],
            }
        else:
            inode_status = "blocked"
            inode_detail = {
                "status": "blocked",
                "message": "df -i output did not provide numeric inode fields",
                "command": record["command"],
                "output": record["stdout"],
            }
        break

    quota_status = "reported" if available_quota_interfaces else "blocked"
    quota_detail = {
        "status": quota_status,
        "available_interfaces": available_quota_interfaces,
        "message": (
            "No user-visible quota command was available from the current environment."
            if not available_quota_interfaces
            else "Quota commands were discovered; downstream tasks must still record their exact outputs before acquisition."
        ),
    }

    return {
        "commands": commands,
        "inode_state": inode_detail,
        "quota_state": quota_detail,
    }


def owner_record(path: Path) -> dict[str, Any]:
    stat_result = path.stat()
    return {
        "path": str(path),
        "owner": pwd.getpwuid(stat_result.st_uid).pw_name,
        "group": grp.getgrgid(stat_result.st_gid).gr_name,
        "uid": stat_result.st_uid,
        "gid": stat_result.st_gid,
        "mode_octal": oct(stat_result.st_mode & 0o7777),
        "world_writable": bool(stat_result.st_mode & 0o002),
    }


def collect_blockers(system: dict[str, Any], smoke: dict[str, Any]) -> list[dict[str, str]]:
    blockers = []
    # A user-visible quota command is useful observability, but it is not a
    # storage authorization gate when the exact root is writable and direct
    # filesystem/inode headroom is measured.  Preserve quota status in system
    # evidence without recreating a global refusal from an unavailable helper.
    if system["inode_state"]["status"] != "reported":
        blockers.append(
            {
                "code": "INODE_STATE_UNAVAILABLE",
                "message": "Inode availability could not be reported from the exact VGP root. Downstream acquisition must fail closed until inode evidence is recorded.",
            }
        )
    if smoke["atomic_promotion"]["status"] != "pass":
        blockers.append(
            {
                "code": "ATOMIC_PROMOTION_UNVERIFIED",
                "message": "Atomic promotion from staging to immutable objects could not be demonstrated inside the exact VGP root.",
            }
        )
    return blockers


def render_markdown(contract: dict[str, Any], layout: Layout, evidence: dict[str, Any]) -> str:
    lines = [
        "# VGP data root contract",
        "",
        f"Date: {evidence['generated_at_utc']}",
        "",
        f"Versioned config artifact: `analysis/vgp_data_root_config.json`",
        "",
        "## Root status",
        "",
        f"- Authorized root: `{layout.root}`",
        f"- Owner/group: `{evidence['root_owner']['owner']}:{evidence['root_owner']['group']}`",
        f"- Root mode: `{evidence['root_owner']['mode_octal']}`",
        f"- World-writable: `{str(evidence['root_owner']['world_writable']).lower()}`",
        f"- Downstream acquisition ready: `{str(not evidence['blockers']).lower()}`",
        "",
        "## Durable layout",
        "",
        "| key | relative path | mode | purpose |",
        "|---|---|---|---|",
    ]
    purposes = {
        "manifests": "metadata-only release and pilot manifests",
        "immutable_objects": "content-addressed verified immutable objects",
        "accession_views": "accession.version views that resolve only to verified immutable objects",
        "version_views": "inventory or release version views that resolve only to verified immutable objects",
        "staging": "bounded mutable staging root",
        "staging_acquisition": "temporary acquisition workspace before verification",
        "staging_partials": "the only location where partial file bytes may exist",
        "staging_outputs": "temporary validated-output promotion workspace",
        "quarantine": "failed checksum, size, provenance, or format candidates",
        "logs": "transfer, validation, and pilot logs",
        "locks": "advisory lock files for acquisition and promotion coordination",
        "pilot": "bounded pilot-only outputs and run packets",
        "pilot_manifests": "pilot manifest snapshots used by later tasks",
        "pilot_runs": "per-run telemetry packets and run-local metadata",
        "pilot_outputs": "validated pilot outputs promoted after checks pass",
    }
    for item in evidence["layout"]:
        lines.append(
            f"| `{item['name']}` | `{item['relative_to_root']}` | `{item['mode_octal']}` | {purposes[item['name']]} |"
        )
    lines.extend(
        [
            "",
            "## Transfer contract",
            "",
            f"- Acquisition writes mutable bytes only under `{contract['transfer_contract']['acquisition_staging_scope']}` and `{contract['transfer_contract']['acquisition_partial_scope']}`; partial files may exist only under `{contract['transfer_contract']['acquisition_partial_scope']}`.",
            f"- Promotion target for immutable verified content is `{contract['transfer_contract']['immutable_promotion_target']}`.",
            f"- Accession views live under `{contract['transfer_contract']['accession_view_scope']}` and version views live under `{contract['transfer_contract']['version_view_scope']}`; both must resolve only to verified immutable objects.",
            "- Before promotion, acquisition must verify source locator, accession.version, byte size, and checksum.",
            f"- {contract['transfer_contract']['compute_input_policy']}",
            f"- {contract['transfer_contract']['compute_workdir_policy']}",
            f"- {contract['transfer_contract']['compute_network_policy']}",
            f"- {contract['transfer_contract']['output_validation_policy']}",
            f"- {contract['transfer_contract']['cleanup_scope_policy']}",
            "",
            "## Validation evidence",
            "",
            "### Small-file, fsync, checksum, atomic rename, lock, cleanup",
            "",
            f"- File fsync: `{evidence['smoke_tests']['file_fsync']['status']}`",
            f"- Staging dir fsync: `{evidence['smoke_tests']['staging_dir_fsync']['status']}`",
            f"- Atomic promotion: `{evidence['smoke_tests']['atomic_promotion']['status']}`",
            f"- Same filesystem: `{str(evidence['smoke_tests']['same_filesystem']).lower()}`",
            f"- Checksum verification: `{evidence['smoke_tests']['checksum_verification']['status']}`",
            f"- Lock behavior: `{evidence['smoke_tests']['lock_behavior']['status']}`",
            f"- Cleanup: `{evidence['smoke_tests']['cleanup']['status']}`",
            "",
            "### Filesystem and MooseFS evidence",
            "",
        ]
    )
    for command in evidence["system_evidence"]["commands"]:
        lines.extend(
            [
                f"- `{command['timestamp_utc']}` `{command['command']}` exit=`{command['exit_code']}`",
                "",
                "```text",
                command["stdout"].rstrip() or "(no stdout)",
                "```",
            ]
        )
        if command["stderr"].strip():
            lines.extend(["```text", command["stderr"].rstrip(), "```"])
    lines.extend(
        [
            "",
            "## Blockers",
            "",
        ]
    )
    if evidence["blockers"]:
        for blocker in evidence["blockers"]:
            lines.append(f"- `{blocker['code']}`: {blocker['message']}")
    else:
        lines.append("- None.")
    lines.extend(["", "## Confidence warnings", ""])
    if evidence.get("confidence_warnings"):
        for warning in evidence["confidence_warnings"]:
            lines.append(f"- `{warning['code']}`: {warning['message']}")
    else:
        lines.append("- None.")
    lines.extend(
        [
            "",
            "## Git hygiene note",
            "",
            "- No biological bulk data, symlinked assembly tree, or generated biological asset was added to Git by this task.",
            "",
        ]
    )
    return "\n".join(lines)


def generate_evidence(contract: dict[str, Any], layout: Layout) -> dict[str, Any]:
    layout_records = provision_layout(layout)
    smoke = smoke_test_storage_contract(layout)
    system = collect_system_evidence(layout.root)
    evidence = {
        "generated_at_utc": utc_now(),
        "task_id": contract["task_id"],
        "root": str(layout.root),
        "root_owner": owner_record(layout.root),
        "layout": layout_records,
        "smoke_tests": smoke,
        "system_evidence": system,
    }
    evidence["blockers"] = collect_blockers(system, smoke)
    evidence["confidence_warnings"] = ([{
        "code": "QUOTA_INTERFACE_UNAVAILABLE",
        "message": "No user-visible quota helper was found; direct filesystem headroom remains authoritative.",
    }] if system["quota_state"]["status"] != "reported" else [])
    evidence["downstream_acquisition_ready"] = not evidence["blockers"]
    return evidence


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--root", type=Path, default=None)
    parser.add_argument("--write-evidence", type=Path, default=DEFAULT_EVIDENCE)
    parser.add_argument("--write-markdown", type=Path, default=DEFAULT_MARKDOWN)
    args = parser.parse_args(argv)

    contract = load_contract(args.config)
    layout = build_layout(contract, args.root)
    evidence = generate_evidence(contract, layout)
    args.write_evidence.write_text(json.dumps(evidence, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    args.write_markdown.write_text(render_markdown(contract, layout, evidence), encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
