#!/usr/bin/env python3
"""Assert the quality-reviewed comprehensive VGP task graph.

The validator is deliberately read-only.  It parses WG's graph JSONL and does
not contact remote sources, transfer biological data, realize environments, or
submit scheduler work.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys
from typing import Iterable, Mapping


CONTROLLED_TASKS = (
    "design-vgp-comprehensive",
    "mirror-vgp-freeze1",
    "implement-vgp-10-pilot",
    "acquire-vgp-10-pilot",
    "design-gbgc-evidence",
    "run-vgp-10-pilot",
    "review-vgp-10-pilot",
    "scale-vgp-core",
    "pilot-vgp-phylo-gbgc",
    "pilot-pedigree-gbgc",
    "synthesize-vgp-program",
)

# Dot-prefixed WG lifecycle edges are validated for referential integrity but
# excluded from this exact scientific dependency comparison.
EXPECTED_AFTER = {
    "design-vgp-comprehensive": {"quality-vgp-psmc"},
    "mirror-vgp-freeze1": {"design-vgp-comprehensive"},
    "implement-vgp-10-pilot": {"design-vgp-comprehensive"},
    "acquire-vgp-10-pilot": {"design-vgp-comprehensive"},
    "design-gbgc-evidence": {"design-vgp-comprehensive"},
    "run-vgp-10-pilot": {
        "implement-vgp-10-pilot",
        "acquire-vgp-10-pilot",
    },
    "review-vgp-10-pilot": {"run-vgp-10-pilot"},
    "scale-vgp-core": {"review-vgp-10-pilot", "mirror-vgp-freeze1"},
    "pilot-vgp-phylo-gbgc": {"design-gbgc-evidence", "mirror-vgp-freeze1"},
    "pilot-pedigree-gbgc": {"design-gbgc-evidence", "acquire-vgp-10-pilot"},
    "synthesize-vgp-program": {
        "scale-vgp-core",
        "pilot-vgp-phylo-gbgc",
        "pilot-pedigree-gbgc",
    },
}

EXPECTED_TIMEOUTS = {
    "design-vgp-comprehensive": "1d",
    "mirror-vgp-freeze1": "14d",
    "implement-vgp-10-pilot": "2d",
    "acquire-vgp-10-pilot": "5d",
    "design-gbgc-evidence": "1d",
    "run-vgp-10-pilot": "7d",
    "review-vgp-10-pilot": "1d",
    "scale-vgp-core": "21d",
    "pilot-vgp-phylo-gbgc": "7d",
    "pilot-pedigree-gbgc": "7d",
    "synthesize-vgp-program": "2d",
}

REQUIRED_SNIPPETS = {
    "design-vgp-comprehensive": (
        "exactly ten primary VGP H1/H2 pairs plus at least five pre-ranked alternates",
        "Hi-C is supporting evidence for long-range phasing, not a universal core gate",
        "Annotation absence is never a genome-wide diversity or PSMC veto",
        "same-pair PSMC is not statistically independent of same-pair diversity",
        "H1/H2-only evidence must never be called direct conversion or biased transmission",
    ),
    "mirror-vgp-freeze1": (
        "dc1b2af5a7741b97d66fb10cb2bce97f41765cdf",
        "9c58420484a8b76a2d6175b7c26bf709e68bdc726a67fc7541b8c2b5a2fc13a4",
        "717 physical lines and 716 data rows",
        "rsync://hgdownload.soe.ucsc.edu/hubs/VGP/",
        "historical approximately 967 GB whole-collection and approximately 520 GB FASTA-only figures",
        "interruption-safe resume",
        "compute and record SHA-256",
        "no arbitrary laptop-scale global memory or total-byte ceiling",
    ),
    "implement-vgp-10-pilot": (
        "--num-mappings 1:1",
        "IMPG must consume that exact alignment",
        "native graph partitions",
        "reason-coded complement mask",
        "diploid H1-coordinate consensus",
        "at least 100 block bootstraps",
        "exact native annotation",
        "no arbitrary global memory or byte ceiling",
    ),
    "acquire-vgp-10-pilot": (
        "exactly the ten approved primary VGP pairs",
        "at least five alternates",
        "Hi-C absence is not by itself a core refusal",
        "Raw-read acquisition is selective",
        "compute/reverify local SHA-256",
        "Never silently substitute",
    ),
    "design-gbgc-evidence": (
        "Direct evidence uses complete pedigrees",
        "Population evidence uses multi-individual",
        "Historical phylogenetic evidence uses orthologous",
        "Non-allelic evidence tests conversion among paralogs",
        "cannot establish conversion direction, event transmission, or GC-biased transmission",
        "population and non-allelic branches remain design-only",
    ),
    "run-vgp-10-pilot": (
        "all ten immutable primary slots",
        "SweepGA --num-mappings 1:1",
        "IMPG index/partition/query/lace extraction",
        "at least 100 boundary-aware block bootstraps",
        "no silent substitutions",
        "use GNU Guix only",
    ),
    "review-vgp-10-pilot": (
        "at least 100 bootstrap attempts and at least 95% finite successful replicates",
        "at least 8 of 10 primary slots",
        "25% at the median and 50% at the 95th percentile",
        "NOT_RUN/DESIGN_ONLY",
        "cannot be presented as independent validation of pair-derived diversity",
    ),
    "scale-vgp-core": (
        "across every audited eligible same-individual",
        "Closed-world account for all 716 frozen catalog rows",
        "Hi-C absence is not an automatic veto",
        "Annotation absence is never a core veto",
        "at least 100 boundary-aware block bootstraps",
        "no arbitrary laptop-scale global memory or byte ceiling",
    ),
    "pilot-vgp-phylo-gbgc": (
        "Semi-complete genomes are admissible",
        "long-term substitution signatures",
        "does not observe direct conversion events or biased transmission",
        "does not substitute for multi-individual population-frequency evidence",
        "non-allelic estimand",
    ),
    "pilot-pedigree-gbgc": (
        "establish parent-of-origin, direction, and transmission",
        "static H1/H2 pair",
        "is not direct evidence",
        "paralogous/segmental-duplication candidates",
        "not cross-vertebrate rates",
    ),
    "synthesize-vgp-program": (
        "Do not claim that assembly-derived PSMC independently predicts or validates diversity",
        "Do not call H1/H2 heterozygous WS/SW states direct conversion or biased transmission",
        "four separate gene-conversion evidence rows",
        "population and non-allelic rows must remain NOT_RUN/DESIGN_ONLY",
        "same-pair diversity and PSMC",
    ),
}

DELIVERABLE_SNIPPETS = {
    "design-vgp-comprehensive": (
        "analysis/vgp_comprehensive_research_plan.md",
        "analysis/vgp_10_pair_manifest.tsv",
        "analysis/vgp_10_pair_alternates.tsv",
        "analysis/vgp_analysis_manifest.json",
    ),
    "mirror-vgp-freeze1": (
        "analysis/vgp_freeze1_source_inventory.tsv",
        "analysis/vgp_freeze1_mirror_manifest.tsv",
        "analysis/vgp_freeze1_mirror_summary.json",
        "analysis/vgp_freeze1_mirror_handoff.md",
    ),
    "implement-vgp-10-pilot": (
        "analysis/vgp_10_pilot_output_schema.json",
        "analysis/vgp_10_pilot_workflow_handoff.md",
    ),
    "acquire-vgp-10-pilot": (
        "analysis/vgp_10_pilot_acquisition_manifest.tsv",
        "analysis/vgp_10_pilot_object_inventory.tsv",
        "analysis/vgp_direct_control_acquisition_manifest.tsv",
        "analysis/vgp_10_pilot_acquisition_handoff.md",
    ),
    "design-gbgc-evidence": (
        "analysis/gene_conversion_evidence_plan.md",
        "analysis/gene_conversion_dataset_manifest.tsv",
        "analysis/gene_conversion_estimand_manifest.tsv",
        "analysis/gene_conversion_claim_matrix.tsv",
    ),
    "run-vgp-10-pilot": (
        "analysis/vgp_10_pilot_result_manifest.tsv",
        "analysis/vgp_10_pilot_qc.tsv",
        "analysis/vgp_10_pilot_resource_telemetry.tsv",
        "analysis/vgp_10_pilot_results.md",
    ),
    "review-vgp-10-pilot": (
        "analysis/vgp_10_pilot_review.md",
        "analysis/vgp_10_pilot_review_gates.tsv",
        "analysis/vgp_10_pilot_review_decision.json",
    ),
    "scale-vgp-core": (
        "analysis/vgp_core_scaleout_manifest.tsv",
        "analysis/vgp_core_scaleout_qc.tsv",
        "analysis/vgp_core_scaleout_telemetry.tsv",
        "analysis/vgp_core_scaleout_results.md",
    ),
    "pilot-vgp-phylo-gbgc": (
        "analysis/vgp_phylogenetic_gbgc_clade_manifest.tsv",
        "analysis/vgp_phylogenetic_gbgc_qc.tsv",
        "analysis/vgp_phylogenetic_gbgc_results.tsv",
        "analysis/vgp_phylogenetic_gbgc_pilot.md",
    ),
    "pilot-pedigree-gbgc": (
        "analysis/direct_gene_conversion_dataset_manifest.tsv",
        "analysis/direct_gene_conversion_tracts.tsv",
        "analysis/direct_gene_conversion_summary.tsv",
        "analysis/direct_gene_conversion_pilot.md",
    ),
    "synthesize-vgp-program": (
        "analysis/vgp_comprehensive_synthesis.md",
        "analysis/vgp_comprehensive_claim_ledger.tsv",
        "analysis/vgp_comprehensive_final_manifest.json",
    ),
}

# These were bounded execution limits in the repaired six-row annotation pilot.
# They are not scientifically valid global eligibility ceilings for this new
# catalog-wide program and must not reappear in any controlled definition.
FORBIDDEN_LEGACY_GLOBAL_CAPS = (
    "6 species",
    "120 GiB compressed inputs",
    "750 GiB scratch",
    "1,500 core-hours",
    "2 concurrent species",
    "256 GiB per job",
    "180 GB read",
    "200 GB write",
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


def _domain(dependencies: Iterable[str]) -> set[str]:
    return {item for item in dependencies if not item.startswith(".")}


def validate_tasks(tasks: Mapping[str, Mapping[str, object]]) -> list[str]:
    errors: list[str] = []
    missing = sorted(set(CONTROLLED_TASKS) - set(tasks))
    if missing:
        return [f"missing controlled tasks: {', '.join(missing)}"]

    for task_id in CONTROLLED_TASKS:
        task = tasks[task_id]
        try:
            after = _as_strings(task.get("after"))
        except ValueError as exc:
            errors.append(f"{task_id}: invalid after list: {exc}")
            continue
        phantom = sorted(item for item in after if item not in tasks)
        if phantom:
            errors.append(f"{task_id}: phantom dependencies: {', '.join(phantom)}")
        actual_domain = _domain(after)
        if actual_domain != EXPECTED_AFTER[task_id]:
            errors.append(
                f"{task_id}: scientific after={sorted(actual_domain)!r}, "
                f"expected={sorted(EXPECTED_AFTER[task_id])!r}"
            )
        expected_assignment = f".assign-{task_id}"
        if expected_assignment not in after:
            errors.append(f"{task_id}: missing WG lifecycle dependency {expected_assignment}")
        if task.get("timeout") != EXPECTED_TIMEOUTS[task_id]:
            errors.append(
                f"{task_id}: timeout={task.get('timeout')!r}, "
                f"expected={EXPECTED_TIMEOUTS[task_id]!r}"
            )

    expected_children: dict[str, set[str]] = {
        task_id: set() for task_id in CONTROLLED_TASKS
    }
    for child, parents in EXPECTED_AFTER.items():
        for parent in parents:
            if parent in expected_children:
                expected_children[parent].add(child)
    for task_id, expected in expected_children.items():
        try:
            before = _domain(_as_strings(tasks[task_id].get("before")))
        except ValueError as exc:
            errors.append(f"{task_id}: invalid before list: {exc}")
            continue
        actual_controlled = before.intersection(CONTROLLED_TASKS)
        if actual_controlled != expected:
            errors.append(
                f"{task_id}: controlled before={sorted(actual_controlled)!r}, "
                f"expected={sorted(expected)!r}"
            )

    for task_id in CONTROLLED_TASKS:
        description = tasks[task_id].get("description")
        if not isinstance(description, str):
            errors.append(f"{task_id}: description is missing or not text")
            continue
        for snippet in REQUIRED_SNIPPETS[task_id]:
            if snippet not in description:
                errors.append(f"{task_id}: missing required definition text {snippet!r}")
        for deliverable in DELIVERABLE_SNIPPETS[task_id]:
            if deliverable not in description:
                errors.append(f"{task_id}: missing deliverable {deliverable!r}")
        for legacy_cap in FORBIDDEN_LEGACY_GLOBAL_CAPS:
            if legacy_cap in description:
                errors.append(
                    f"{task_id}: contains forbidden legacy global cap {legacy_cap!r}"
                )

    # Critical fork/join assertions are repeated explicitly so a failure report
    # states the scientific bypass, not merely a generic edge mismatch.
    if _domain(_as_strings(tasks["run-vgp-10-pilot"].get("after"))) != {
        "implement-vgp-10-pilot",
        "acquire-vgp-10-pilot",
    }:
        errors.append("run-vgp-10-pilot bypasses the implementation/acquisition join")
    if _domain(_as_strings(tasks["scale-vgp-core"].get("after"))) != {
        "review-vgp-10-pilot",
        "mirror-vgp-freeze1",
    }:
        errors.append("scale-vgp-core bypasses independent review or full mirror completion")
    if _domain(_as_strings(tasks["synthesize-vgp-program"].get("after"))) != {
        "scale-vgp-core",
        "pilot-vgp-phylo-gbgc",
        "pilot-pedigree-gbgc",
    }:
        errors.append("synthesis bypasses a required executed evidence branch")

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
        errors = validate_tasks(load_tasks(args.graph))
    except (OSError, ValueError) as exc:
        print(f"VGP_COMPREHENSIVE_WG_ASSERTIONS_FAIL: {exc}", file=sys.stderr)
        return 2
    if errors:
        print("VGP_COMPREHENSIVE_WG_ASSERTIONS_FAIL", file=sys.stderr)
        for error in errors:
            print(f"- {error}", file=sys.stderr)
        return 1
    edge_count = sum(len(parents) for parents in EXPECTED_AFTER.values())
    print(
        "VGP_COMPREHENSIVE_WG_ASSERTIONS_OK "
        f"tasks={len(CONTROLLED_TASKS)} scientific_edges={edge_count} "
        "pilot=10 alternates>=5 mirror=full-freeze1 gbgc=4-disjoint-branches"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
