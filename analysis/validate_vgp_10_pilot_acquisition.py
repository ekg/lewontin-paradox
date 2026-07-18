#!/usr/bin/env python3
"""Fail-closed validator for the committed VGP ten-pair acquisition handoff."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

from analysis import vgp_10_pilot_acquisition as acq


ROOT = Path(__file__).resolve().parents[1]
EXPECTED_PRIMARY_SHA = "bf6c9ff647aed332bfc002bf803e8307203b51432343f2eca6d95a6c80d82997"
EXPECTED_ALTERNATE_SHA = "55127b18f0f17f6673cc0367c60736207a7e2184198cd3855a7e6ea83f39c52e"


def _sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def validate(*, reverify_local: bool = False) -> list[str]:
    errors: list[str] = []
    required = [acq.PAIR_OUTPUT, acq.OBJECT_OUTPUT, acq.DIRECT_OUTPUT, acq.SUMMARY_OUTPUT, acq.HANDOFF_OUTPUT]
    for path in required:
        if not path.is_file():
            errors.append(f"missing deliverable: {path}")
    if errors:
        return errors
    if _sha(acq.PRIMARY) != EXPECTED_PRIMARY_SHA:
        errors.append("frozen primary roster digest drift")
    if _sha(acq.ALTERNATES) != EXPECTED_ALTERNATE_SHA:
        errors.append("frozen alternate roster digest drift")

    pairs = _rows(acq.PAIR_OUTPUT)
    objects = _rows(acq.OBJECT_OUTPUT)
    controls = _rows(acq.DIRECT_OUTPUT)
    summary = json.loads(acq.SUMMARY_OUTPUT.read_text(encoding="utf-8"))
    primary = [row for row in pairs if row["roster_type"] == "primary"]
    alternates = [row for row in pairs if row["roster_type"] == "alternate"]
    if [row["selection_id"] for row in primary] != [f"P{i:02d}" for i in range(1, 11)]:
        errors.append("primary closed world is not exactly P01..P10")
    if [row["selection_id"] for row in alternates] != [f"A{i:02d}" for i in range(1, 7)]:
        errors.append("alternate closed world is not exactly A01..A06")
    for row in primary:
        if row["activation_status"] != "active_primary" or row["core_identity_status"] != "pass" or row["core_acquisition_status"] != "verified":
            errors.append(f"{row['selection_id']}: primary identity/acquisition gate did not pass")
        if row["amendment_id"] != "none":
            errors.append(f"{row['selection_id']}: unexpected alternate amendment")
        if row["qv_status"] != "missing_exact_final_sequence_qv":
            errors.append(f"{row['selection_id']}: QV missingness was changed or imputed")
        if row["repeat_report_status"] != "missing_not_published":
            errors.append(f"{row['selection_id']}: repeat-report missingness was not explicit")
        for side in ("h1", "h2"):
            if not acq.HEX64.fullmatch(row[f"{side}_frozen_report_sha256"]):
                errors.append(f"{row['selection_id']}: invalid frozen {side} report SHA-256")
            for field in ("length_bp", "contigs", "contig_n50_bp", "scaffold_n50_bp"):
                if not row[f"{side}_{field}"].isdigit() or int(row[f"{side}_{field}"]) <= 0:
                    errors.append(f"{row['selection_id']}: invalid {side} {field}")
    for row in alternates:
        if row["activation_status"] != "standby_not_triggered" or row["core_acquisition_status"] != "not_activated":
            errors.append(f"{row['selection_id']}: alternate activated without amendment")

    if len({row["object_id"] for row in objects}) != len(objects):
        errors.append("duplicate VGP object IDs")
    allowed_status = {"verified", "reused", "missing", "superseded", "quarantined"}
    for row in objects:
        if row["status"] not in allowed_status:
            errors.append(f"{row['object_id']}: nonterminal status {row['status']}")
        if not row["source_url"].startswith("https://"):
            errors.append(f"{row['object_id']}: non-HTTPS source")
        try:
            if int(row["expected_bytes"]) < 0:
                raise ValueError
        except ValueError:
            errors.append(f"{row['object_id']}: invalid expected bytes")
        if row["object_role"] in acq.BULK_RAW_ROLES:
            errors.append(f"{row['object_id']}: forbidden raw-read role acquired")
        if row["status"] in {"verified", "reused"}:
            if not acq.HEX64.fullmatch(row["local_sha256"]):
                errors.append(f"{row['object_id']}: accepted object lacks local SHA-256")
            path = Path(row["local_path"])
            if not path.is_file():
                errors.append(f"{row['object_id']}: accepted CAS path is absent")
            elif reverify_local:
                observed_sha, _, observed_bytes = acq.hash_file(path)
                if observed_sha != row["local_sha256"] or observed_bytes != int(row["expected_bytes"]):
                    errors.append(f"{row['object_id']}: live local SHA-256/size revalidation failed")
    for selection in (f"P{i:02d}" for i in range(1, 11)):
        selected = [row for row in objects if row["selection_id"] == selection]
        roles = {(row["side"], row["object_role"]) for row in selected if row["status"] in {"verified", "reused"}}
        required_roles = {
            (side, role)
            for side in ("h1", "h2")
            for role in ("dataset_report", "checksum_catalog", "genome_fasta", "assembly_report", "assembly_stats")
        }
        if not required_roles.issubset(roles):
            errors.append(f"{selection}: accepted core object set is incomplete")
    repeat_rows = [row for row in objects if row["object_role"] == "repeat_report"]
    if len(repeat_rows) != 20 or any(row["status"] != "missing" for row in repeat_rows):
        errors.append("repeat-report closed-world missing ledger is not exactly 20 H1/H2 rows")
    if len([row for row in objects if row["roster_type"] == "alternate" and row["object_role"] == "genome_fasta" and row["status"] == "superseded"]) != 12:
        errors.append("alternate FASTA supersession ledger is not exactly 12 rows")

    if {row["control_id"] for row in controls} != {"D01"} or any(row["selection_status"] != "selected" for row in controls):
        errors.append("direct control is not exclusively selected D01")
    if len(controls) != 10:
        errors.append("D01 object/exclusion ledger is not exactly ten rows")
    relationship = controls[0].get("relationships", "") if controls else ""
    for tetrad in ("29", "38", "40", "58", "62", "34", "35", "39", "51", "53", "68", "69", "76"):
        if f"tetrad {tetrad}:" not in relationship:
            errors.append(f"D01 lacks complete relationship mapping for tetrad {tetrad}")
    excluded = [row for row in controls if row["object_role"] == "bulk_raw_archive_exclusion_aggregate"]
    if len(excluded) != 1 or excluded[0]["status"] != "superseded" or int(excluded[0]["expected_bytes"]) != 356_984_558_868:
        errors.append("D01 951-object/356984558868-byte raw archive exclusion drift")
    accepted_controls = [row for row in controls if row["status"] in {"verified", "reused"}]
    if len(accepted_controls) != 9:
        errors.append("D01 selected metadata/supplement/reference set is not exactly nine accepted objects")
    if reverify_local:
        for row in accepted_controls:
            path = Path(row["local_path"])
            if not path.is_file() or _sha(path) != row["local_sha256"] or path.stat().st_size != int(row["expected_bytes"]):
                errors.append(f"{row['object_id']}: live direct-control SHA-256/size revalidation failed")

    recomputed = acq.summarize_inventory([*objects, *controls])
    for key in ("planned", "transferred", "verified", "newly_promoted", "reused", "missing", "superseded", "quarantined", "terminal_disposition_reconciliation"):
        if summary["accounting"].get(key) != recomputed.get(key):
            errors.append(f"summary accounting drift for {key}")
    scope = summary.get("scope_proof", {})
    for key in ("unmanifested_bulk_raw_read_objects_acquired", "unmanifested_bulk_raw_read_bytes_acquired", "ena_fastq_payload_requests", "slurm_jobs_submitted"):
        if scope.get(key) != 0:
            errors.append(f"scope proof is nonzero for {key}")
    demo = summary.get("mechanism_demonstration", {})
    if demo.get("resume", {}).get("resume_from_bytes") != 5 or demo.get("mirror_reuse", {}).get("status") != "reused" or demo.get("quarantine", {}).get("status") != "quarantined":
        errors.append("resume/reuse/quarantine mechanism evidence is incomplete")
    return errors


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reverify-local", action="store_true")
    args = parser.parse_args()
    errors = validate(reverify_local=args.reverify_local)
    if errors:
        for error in errors:
            print("ERROR:", error)
        return 1
    print("VGP ten-pair acquisition validation passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
