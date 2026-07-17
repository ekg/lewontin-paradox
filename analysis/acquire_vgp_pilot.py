#!/usr/bin/env python3
"""Fail-closed acquisition entrypoint for the bounded VGP pilot."""

from __future__ import annotations

import argparse
import csv
import json
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Sequence

from analysis import gate_vgp_pilot as gate
from analysis.tier3_common import Tier3ValidationError


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_GATE = PROJECT_ROOT / "analysis" / "vgp_pilot_gate.json"
DEFAULT_MANIFEST = PROJECT_ROOT / "analysis" / "vgp_pilot_manifest.tsv"
DEFAULT_ROOT_CONFIG = PROJECT_ROOT / "analysis" / "vgp_data_root_config.json"
DEFAULT_OUTPUT_MANIFEST = PROJECT_ROOT / "analysis" / "vgp_pilot_acquisition_manifest.tsv"
DEFAULT_OUTPUT_REPORT = PROJECT_ROOT / "analysis" / "vgp_pilot_acquisition_report.md"

MANIFEST_FIELDS = [
    "run_id",
    "generated_at_utc",
    "record_type",
    "status",
    "candidate_id",
    "asset_role",
    "accession_version",
    "source_url",
    "expected_bytes",
    "observed_bytes",
    "expected_sha256",
    "observed_sha256",
    "provider_md5",
    "failure_code",
    "failure_source",
    "failure_message",
    "quarantine_path",
    "promoted_path",
]
SLURM_ENV_VARS = ("SLURM_JOB_ID", "SLURM_ARRAY_JOB_ID", "SLURM_ARRAY_TASK_ID", "SLURM_CLUSTER_NAME")


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def make_run_id() -> str:
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    return f"vgp-pilot-acquire-{stamp}"


def write_manifest(path: Path, rows: Sequence[Mapping[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def detected_slurm_environment() -> dict[str, str]:
    return {name: os.environ[name] for name in SLURM_ENV_VARS if name in os.environ}


def classify_failure_code(error: Exception) -> str:
    text = str(error)
    if "gate decision is NO_GO" in text:
        return "GATE_NO_GO"
    if "manifest digest mismatch" in text:
        return "MANIFEST_DIGEST_MISMATCH"
    if "root contract digest mismatch" in text:
        return "ROOT_CONTRACT_DIGEST_MISMATCH"
    if "cap vector digest mismatch" in text:
        return "CAP_VECTOR_DIGEST_MISMATCH"
    if "gate decision hash does not match" in text or "cap vector hash does not match the gate payload" in text:
        return "GATE_TAMPERED"
    return "PRECHECK_FAILED"


def refusal_rows(
    *,
    run_id: str,
    generated_at_utc: str,
    gate_path: Path,
    gate_payload: Mapping[str, Any] | None,
    error: Exception,
) -> list[dict[str, str]]:
    code = classify_failure_code(error)
    decision_sha = ""
    blockers = []
    if gate_payload is not None:
        decision_sha = str(gate_payload.get("decision_sha256", ""))
        blockers = list(gate_payload.get("blockers", []))

    rows = [
        {
            "run_id": run_id,
            "generated_at_utc": generated_at_utc,
            "record_type": "run_summary",
            "status": "refused_preflight",
            "candidate_id": "",
            "asset_role": "authorization_boundary",
            "accession_version": "",
            "source_url": str(gate_path),
            "expected_bytes": "0",
            "observed_bytes": "0",
            "expected_sha256": decision_sha,
            "observed_sha256": decision_sha,
            "provider_md5": "",
            "failure_code": code,
            "failure_source": str(gate_path),
            "failure_message": str(error),
            "quarantine_path": "",
            "promoted_path": "",
        }
    ]
    for blocker in blockers:
        rows.append(
            {
                "run_id": run_id,
                "generated_at_utc": generated_at_utc,
                "record_type": "gate_blocker",
                "status": "refused_preflight",
                "candidate_id": "",
                "asset_role": "authorization_boundary",
                "accession_version": "",
                "source_url": blocker["source"],
                "expected_bytes": "0",
                "observed_bytes": "0",
                "expected_sha256": "",
                "observed_sha256": "",
                "provider_md5": "",
                "failure_code": blocker["code"],
                "failure_source": blocker["source"],
                "failure_message": blocker["message"],
                "quarantine_path": "",
                "promoted_path": "",
            }
        )
    return rows


def render_report(
    *,
    run_id: str,
    generated_at_utc: str,
    gate_path: Path,
    gate_payload: Mapping[str, Any] | None,
    manifest_path: Path,
    root_config_path: Path,
    output_manifest_path: Path,
    refusal_error: Exception,
) -> str:
    slurm_env = detected_slurm_environment()
    decision_status = gate_payload["decision"]["status"] if gate_payload is not None else "UNREADABLE"
    manifest_digest = gate_payload["authorization_boundary"]["manifest_digest"] if gate_payload is not None else "UNAVAILABLE"
    root_digest = gate_payload["authorization_boundary"]["root_contract_digest"] if gate_payload is not None else "UNAVAILABLE"
    cap_digest = gate_payload["authorization_boundary"]["cap_vector_digest"] if gate_payload is not None else "UNAVAILABLE"
    lines = [
        "# VGP pilot acquisition report",
        "",
        f"- Run ID: `{run_id}`",
        f"- Generated at: `{generated_at_utc}`",
        f"- Gate path: `{gate_path}`",
        f"- Gate decision: `{decision_status}`",
        f"- Acquisition status: `refused_preflight`",
        f"- Refused before first biological byte: `true`",
        f"- Slurm environment detected: `{str(bool(slurm_env)).lower()}`",
        f"- Biological payload bytes transferred: `0`",
        f"- Verified immutable objects promoted: `0`",
        f"- Quarantine objects written: `0`",
        f"- Accessions promoted into views: `0`",
        f"- Output manifest: `{output_manifest_path}`",
        "",
        "## Authorization Boundary",
        "",
        f"- Manifest path: `{manifest_path}`",
        f"- Recorded manifest digest: `{manifest_digest}`",
        f"- Root contract path: `{root_config_path}`",
        f"- Recorded root contract digest: `{root_digest}`",
        f"- Recorded cap vector digest: `{cap_digest}`",
        "",
        "## Refusal Reason",
        "",
        f"- `{classify_failure_code(refusal_error)}`: {refusal_error}",
        "",
        "## Gate Blockers",
        "",
    ]
    blockers = list(gate_payload.get("blockers", [])) if gate_payload is not None else []
    if blockers:
        for blocker in blockers:
            lines.append(f"- `{blocker['code']}`: {blocker['message']} ({blocker['source']})")
    else:
        lines.append("- None recorded because the gate payload could not be read.")
    lines.extend(
        [
            "",
            "## Validation Notes",
            "",
            "- Live gate authorization was executed against the current manifest, root contract, and recomputed cap vector.",
            "- The current gate is not exactly `GO`, so no provider request, staging write, partial file, checksum verification, quarantine action, or promotion step was attempted.",
            "- Exact-reference/native-annotation linkage validation under pinned GNU Guix was not re-run because zero assets were authorized or acquired.",
            "- No Slurm command was submitted by this task.",
            "- No biological bulk file was added to Git by this task.",
        ]
    )
    if slurm_env:
        lines.extend(
            [
                "",
                "## Detected Slurm Environment",
                "",
                "```json",
                json.dumps(slurm_env, indent=2, sort_keys=True),
                "```",
            ]
        )
    lines.append("")
    return "\n".join(lines)


def run(
    *,
    gate_path: Path = DEFAULT_GATE,
    manifest_path: Path = DEFAULT_MANIFEST,
    root_config_path: Path = DEFAULT_ROOT_CONFIG,
    output_manifest_path: Path = DEFAULT_OUTPUT_MANIFEST,
    output_report_path: Path = DEFAULT_OUTPUT_REPORT,
) -> dict[str, Any]:
    run_id = make_run_id()
    generated_at_utc = utc_now()
    gate_payload = None
    refusal_error: Exception
    try:
        gate_payload = gate.load_gate(gate_path)
        gate.authorize_gate_action(gate_path, manifest_path, root_config_path, "acquire")
    except (Tier3ValidationError, FileNotFoundError, KeyError, ValueError) as error:
        refusal_error = error
    else:  # pragma: no cover - current task is intentionally fail-closed on the present NO_GO gate.
        raise Tier3ValidationError("authorized VGP pilot acquisition is not implemented in this repository state")

    rows = refusal_rows(
        run_id=run_id,
        generated_at_utc=generated_at_utc,
        gate_path=gate_path,
        gate_payload=gate_payload,
        error=refusal_error,
    )
    write_manifest(output_manifest_path, rows)
    report_text = render_report(
        run_id=run_id,
        generated_at_utc=generated_at_utc,
        gate_path=gate_path,
        gate_payload=gate_payload,
        manifest_path=manifest_path,
        root_config_path=root_config_path,
        output_manifest_path=output_manifest_path,
        refusal_error=refusal_error,
    )
    output_report_path.parent.mkdir(parents=True, exist_ok=True)
    output_report_path.write_text(report_text, encoding="utf-8")
    return {
        "run_id": run_id,
        "generated_at_utc": generated_at_utc,
        "status": "refused_preflight",
        "failure_code": classify_failure_code(refusal_error),
        "failure_message": str(refusal_error),
        "output_manifest": str(output_manifest_path),
        "output_report": str(output_report_path),
        "row_count": len(rows),
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gate", type=Path, default=DEFAULT_GATE)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--root-config", type=Path, default=DEFAULT_ROOT_CONFIG)
    parser.add_argument("--output-manifest", type=Path, default=DEFAULT_OUTPUT_MANIFEST)
    parser.add_argument("--output-report", type=Path, default=DEFAULT_OUTPUT_REPORT)
    args = parser.parse_args(argv)
    run(
        gate_path=args.gate,
        manifest_path=args.manifest,
        root_config_path=args.root_config,
        output_manifest_path=args.output_manifest,
        output_report_path=args.output_report,
    )
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
