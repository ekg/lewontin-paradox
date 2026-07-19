#!/usr/bin/env python3
"""Atomically promote a fully sentinel-verified VGP canary from checkpoint."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
from pathlib import Path
from typing import Sequence


class PromotionError(RuntimeError):
    """The canary is incomplete, has drifted, or cannot be promoted safely."""


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify_stage_sentinels(root: Path) -> dict[str, int]:
    sentinels = sorted(root.rglob(".complete.json"))
    if not sentinels:
        raise PromotionError(f"no completion sentinels under {root}")
    payload_files = total_bytes = 0
    for sentinel in sentinels:
        try:
            value = json.loads(sentinel.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as error:
            raise PromotionError(f"unreadable stage sentinel: {sentinel}") from error
        files = value.get("files")
        if not isinstance(files, dict):
            raise PromotionError(f"stage sentinel lacks file digest map: {sentinel}")
        for relative, expected in files.items():
            path = sentinel.parent / relative
            try:
                resolved = path.resolve(strict=True)
            except OSError as error:
                raise PromotionError(f"sentinel payload absent: {path}") from error
            if not resolved.is_relative_to(sentinel.parent.resolve()):
                raise PromotionError(f"sentinel payload escapes stage directory: {relative}")
            if not path.is_file():
                raise PromotionError(f"sentinel payload is not a file: {path}")
            observed = sha256_file(path)
            if observed != expected:
                raise PromotionError(f"stage payload digest mismatch: {path}")
            payload_files += 1
            total_bytes += path.stat().st_size
    return {
        "stage_sentinels": len(sentinels),
        "verified_payload_files": payload_files,
        "verified_payload_bytes": total_bytes,
    }


def required_completion(root: Path) -> None:
    required = [
        "preflight/.complete.json",
        "mapping/.complete.json",
        "impg/.complete.json",
        "variants/.complete.json",
        "consensus/.complete.json",
        "psmc/finalize/.complete.json",
        *[f"psmc/replicate-{replicate:03d}/.complete.json" for replicate in range(201)],
    ]
    absent = [relative for relative in required if not (root / relative).is_file()]
    if absent:
        raise PromotionError(f"canary lacks {len(absent)} required completion sentinels: {absent[:5]}")


def promote(source: Path, target: Path, canonical_root: Path, job_id: str) -> dict[str, object]:
    source, target, canonical_root = source.resolve(), target.resolve(), canonical_root.resolve()
    if not source.is_dir():
        raise PromotionError(f"checkpoint run is absent: {source}")
    if not source.is_relative_to(canonical_root) or not target.is_relative_to(canonical_root):
        raise PromotionError("source and target must remain under the canonical VGP root")
    required_completion(source)
    source_audit = verify_stage_sentinels(source)
    if target.exists():
        target_audit = verify_stage_sentinels(target)
        manifest = target / "promotion_manifest.json"
        if not manifest.is_file():
            raise PromotionError("existing target lacks promotion manifest")
        return {"status": "already_promoted", **target_audit, "target": str(target)}
    target.parent.mkdir(parents=True, exist_ok=True)
    partial = target.parent / f".{target.name}.{job_id}.partial"
    if partial.exists():
        raise PromotionError(f"refusing pre-existing promotion staging directory: {partial}")
    shutil.copytree(source, partial, copy_function=shutil.copy2)
    copied_audit = verify_stage_sentinels(partial)
    if copied_audit != source_audit:
        raise PromotionError("copied promotion audit differs from checkpoint audit")
    manifest = {
        "schema_version": "vgp-real-canary-promotion-v1",
        "task_id": "run-vgp-real-canary",
        "authorization_id": "vgp10-auth-20260718-v2",
        "selection_id": "P07",
        "slurm_job_id": job_id,
        "canonical_vgp_root": str(canonical_root),
        "checkpoint_source": str(source),
        "promotion_target": str(target),
        "source_preserved": True,
        "atomic_same_filesystem_rename": True,
        **source_audit,
    }
    manifest_path = partial / "promotion_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    manifest["promotion_manifest_sha256"] = sha256_file(manifest_path)
    os.replace(partial, target)
    final_audit = verify_stage_sentinels(target)
    if final_audit != source_audit:
        raise PromotionError("promoted target differs from verified checkpoint")
    return {"status": "promoted", "target": str(target), **manifest}


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--source", type=Path, required=True)
    result.add_argument("--target", type=Path, required=True)
    result.add_argument("--canonical-root", type=Path, default=Path("/moosefs/erikg/vgp"))
    result.add_argument("--slurm-job-id", required=True)
    result.add_argument("--output", type=Path)
    return result


def main(argv: Sequence[str] | None = None) -> int:
    args = parser().parse_args(argv)
    try:
        result = promote(args.source, args.target, args.canonical_root, args.slurm_job_id)
    except (PromotionError, OSError, shutil.Error) as error:
        print(f"ERROR: {error}")
        return 2
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        partial = args.output.with_suffix(args.output.suffix + ".partial")
        partial.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        partial.replace(args.output)
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
