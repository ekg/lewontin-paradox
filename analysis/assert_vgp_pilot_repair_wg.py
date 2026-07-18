#!/usr/bin/env python3
"""Assert the reviewed VGP repair task graph and definition safety gates.

This validator is intentionally read-only.  It inspects WG's graph.jsonl and
does not contact VGP/NCBI, transfer data, or submit compute work.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys
from typing import Iterable, Mapping


CONTROLLED_TASKS = (
    "quality-vgp-pilot-2",
    "repair-vgp-candidate",
    "regate-vgp-pilot",
    "acquire-repaired-vgp",
    "run-repaired-vgp",
    "audit-vgp-demography",
    "review-repaired-vgp",
    "synthesize-repaired-vgp",
)

# Exact scientific edges. WG-generated .assign-* lifecycle dependencies are
# checked for existence but excluded from this domain-edge comparison.
EXPECTED_AFTER = {
    "quality-vgp-pilot-2": {"integrate-vgp-pilot"},
    "repair-vgp-candidate": {"quality-vgp-pilot-2"},
    "regate-vgp-pilot": {"repair-vgp-candidate"},
    "audit-vgp-demography": {"repair-vgp-candidate"},
    "acquire-repaired-vgp": {"regate-vgp-pilot"},
    "run-repaired-vgp": {"acquire-repaired-vgp"},
    "review-repaired-vgp": {"run-repaired-vgp", "audit-vgp-demography"},
    "synthesize-repaired-vgp": {"review-repaired-vgp"},
}

REQUIRED_SNIPPETS = {
    "repair-vgp-candidate": (
        "717 physical catalog lines",
        "one header plus 716 data rows and 714 unique species",
        "batched",
        "immutable local response caching",
        "resumable checkpoints",
        "Remote SHA-256 or MD5 availability is not a pre-download eligibility requirement",
        "compute local SHA-256",
        "exact accession/version identity",
        "same-individual, correctly phased H1/H2 evidence",
        "post-alignment outputs, not circular pre-download prerequisites",
    ),
    "regate-vgp-pilot": (
        "Unknown quota is not evidence of adequate quota",
        "free space alone must not override a stricter integrated storage/root contract",
        "Remote SHA-256 or MD5 need not exist",
        "post-alignment outputs",
        "No downstream task may interpret any other state as GO",
    ),
    "acquire-repaired-vgp": (
        "Refuse NO_GO",
        "altered digests",
        "transfer zero bytes",
        "compute local SHA-256",
        "reverify it immediately before atomic read-only promotion",
        "same-individual/phased H1/H2 evidence",
    ),
    "run-repaired-vgp": (
        "refuse NO_GO",
        "altered local payload SHA-256",
        "NOT_SUBMITTED",
        "post-alignment outputs",
        "same-individual, correctly phased exact H1/H2 pairs",
        "must not be treated as population genotypes or as demographic input",
    ),
    "audit-vgp-demography": (
        "metadata and literature/source audit only",
        "PSMC:",
        "MSMC2:",
        "SMC++:",
        "VGP H1/H2 are assembly haplotypes",
        "not automatically a diploid heterozygosity dataset",
        "Excluded/circular estimates",
        "batched where possible",
        "cached immutably",
        "resumable",
    ),
    "review-repaired-vgp": (
        "reproduce refusal for NO_GO or altered bound digests",
        "measured post-alignment",
        "Keep PSMC, MSMC2, and SMC++ eligibility distinct",
        "Never assume VGP H1/H2 is a valid demographic genotype dataset",
        "same response pi as an independent predictor",
    ),
    "synthesize-repaired-vgp": (
        "Preserve NO_GO, zero-byte, NOT_SUBMITTED, and empty-result outcomes",
        "PSMC, MSMC2, and SMC++",
        "VGP H1/H2 must never be assumed to be a valid demographic genotype or population dataset",
        "genuinely independent Ne from estimates circular with the response pi",
        "require new explicit user authorization",
        "Do not create ready executable tasks",
    ),
}

CAP_TASKS = (
    "repair-vgp-candidate",
    "regate-vgp-pilot",
    "acquire-repaired-vgp",
    "run-repaired-vgp",
    "review-repaired-vgp",
    "synthesize-repaired-vgp",
)

CAP_SNIPPETS = (
    "6 species",
    "120 GiB compressed inputs",
    "750 GiB scratch",
    "1,500 core-hours",
    "2 concurrent",
    "256 GiB per job",
    "stricter integrated cap",
)


def default_graph_path() -> Path:
    project_root = os.environ.get("WG_PROJECT_ROOT")
    if project_root:
        return Path(project_root) / ".wg" / "graph.jsonl"
    return Path(".wg/graph.jsonl")


def load_tasks(path: Path) -> dict[str, dict[str, object]]:
    tasks: dict[str, dict[str, object]] = {}
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            try:
                record = json.loads(line)
            except json.JSONDecodeError as exc:
                raise ValueError(f"{path}:{line_number}: invalid JSON: {exc}") from exc
            if record.get("kind") != "task" or "id" not in record:
                continue
            task_id = str(record["id"])
            if task_id in tasks:
                raise ValueError(f"{path}:{line_number}: duplicate task id {task_id!r}")
            tasks[task_id] = record
    return tasks


def _as_strings(value: object) -> set[str]:
    if value is None:
        return set()
    if not isinstance(value, list) or not all(isinstance(item, str) for item in value):
        raise ValueError(f"expected a string list, got {value!r}")
    return set(value)


def _domain_dependencies(dependencies: Iterable[str]) -> set[str]:
    return {dependency for dependency in dependencies if not dependency.startswith(".")}


def validate_tasks(tasks: Mapping[str, Mapping[str, object]]) -> list[str]:
    errors: list[str] = []

    missing = sorted(set(CONTROLLED_TASKS) - set(tasks))
    if missing:
        errors.append(f"missing controlled tasks: {', '.join(missing)}")
        return errors

    # Every dependency on every controlled task must resolve, including WG's
    # dot-prefixed assignment lifecycle. This is the no-phantom assertion.
    for task_id in CONTROLLED_TASKS:
        try:
            dependencies = _as_strings(tasks[task_id].get("after"))
        except ValueError as exc:
            errors.append(f"{task_id}: invalid after list: {exc}")
            continue
        phantom = sorted(dependency for dependency in dependencies if dependency not in tasks)
        if phantom:
            errors.append(f"{task_id}: phantom dependencies: {', '.join(phantom)}")

        actual_domain = _domain_dependencies(dependencies)
        expected_domain = EXPECTED_AFTER[task_id]
        if actual_domain != expected_domain:
            errors.append(
                f"{task_id}: scientific after={sorted(actual_domain)!r}, "
                f"expected={sorted(expected_domain)!r}"
            )

    expected_children: dict[str, set[str]] = {task_id: set() for task_id in CONTROLLED_TASKS}
    for child, parents in EXPECTED_AFTER.items():
        for parent in parents:
            if parent in expected_children:
                expected_children[parent].add(child)

    # Assert the reverse representation too. This catches a malformed graph
    # record even if the forward 'after' side looks correct.
    for task_id, expected in expected_children.items():
        try:
            actual = _domain_dependencies(_as_strings(tasks[task_id].get("before")))
        except ValueError as exc:
            errors.append(f"{task_id}: invalid before list: {exc}")
            continue
        controlled_actual = actual.intersection(CONTROLLED_TASKS)
        if controlled_actual != expected:
            errors.append(
                f"{task_id}: controlled before={sorted(controlled_actual)!r}, "
                f"expected={sorted(expected)!r}"
            )

    # Explicit fork/join and fail-closed path checks make intent visible in
    # assertion output, instead of relying only on the general edge loop.
    if EXPECTED_AFTER["regate-vgp-pilot"] != {"repair-vgp-candidate"}:
        errors.append("internal assertion error: repair-to-regate fork changed")
    if EXPECTED_AFTER["audit-vgp-demography"] != {"repair-vgp-candidate"}:
        errors.append("internal assertion error: repair-to-audit fork changed")
    if EXPECTED_AFTER["review-repaired-vgp"] != {
        "run-repaired-vgp",
        "audit-vgp-demography",
    }:
        errors.append("internal assertion error: review join changed")
    if _domain_dependencies(_as_strings(tasks["acquire-repaired-vgp"].get("after"))) != {
        "regate-vgp-pilot"
    }:
        errors.append("acquisition bypasses or does not depend exclusively on regate GO path")
    if _domain_dependencies(_as_strings(tasks["run-repaired-vgp"].get("after"))) != {
        "acquire-repaired-vgp"
    }:
        errors.append("compute bypasses or does not depend exclusively on acquisition path")

    for task_id, snippets in REQUIRED_SNIPPETS.items():
        description = tasks[task_id].get("description")
        if not isinstance(description, str):
            errors.append(f"{task_id}: description is missing or not text")
            continue
        for snippet in snippets:
            if snippet not in description:
                errors.append(f"{task_id}: missing required definition text {snippet!r}")

    for task_id in CAP_TASKS:
        description = tasks[task_id].get("description")
        if not isinstance(description, str):
            continue
        for snippet in CAP_SNIPPETS:
            if snippet not in description:
                errors.append(f"{task_id}: missing cap text {snippet!r}")

    audit_description = str(tasks["audit-vgp-demography"].get("description", ""))
    for forbidden_authorization in (
        "download no bulk reads or VCFs",
        "submit no compute",
        "run no PSMC, MSMC2, SMC++, or other demographic inference",
    ):
        if forbidden_authorization not in audit_description:
            errors.append(
                "audit-vgp-demography: metadata-only refusal text missing "
                f"{forbidden_authorization!r}"
            )

    return errors


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--graph",
        type=Path,
        default=default_graph_path(),
        help="WG graph.jsonl path (default: $WG_PROJECT_ROOT/.wg/graph.jsonl)",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        tasks = load_tasks(args.graph)
        errors = validate_tasks(tasks)
    except (OSError, ValueError) as exc:
        print(f"VGP_PILOT_REPAIR_WG_ASSERTIONS_FAIL: {exc}", file=sys.stderr)
        return 2

    if errors:
        print("VGP_PILOT_REPAIR_WG_ASSERTIONS_FAIL", file=sys.stderr)
        for error in errors:
            print(f"- {error}", file=sys.stderr)
        return 1

    edge_count = sum(len(parents) for parents in EXPECTED_AFTER.values())
    print(
        "VGP_PILOT_REPAIR_WG_ASSERTIONS_OK "
        f"tasks={len(CONTROLLED_TASKS)} scientific_edges={edge_count} "
        "fork=repair:[regate,audit] join=review:[run,audit] go_bypass=none"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
