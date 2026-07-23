#!/usr/bin/env python3
"""Plan and audit bounded H1-coordinate IMPG work without a whole-genome graph.

The durable IMPG objects are an alignment index and its native partition table.
Every query is a complete, bounded group of consecutive H1-native partitions.
Native partition boundaries therefore also become output ownership boundaries:
no partition or callable base can be queried by two ranges, and no variant needs
post-hoc ownership guessing at an artificial tile edge.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Mapping, Sequence

from analysis.vgp_10_pilot import Interval, PilotError, canonical_json, sha256_file


PLAN_SCHEMA = "vgp-bounded-h1-range-plan-v1"


def read_fai(path: Path | str) -> list[tuple[str, int]]:
    rows: list[tuple[str, int]] = []
    seen: set[str] = set()
    with Path(path).open(encoding="ascii") as handle:
        for number, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                raise PilotError(f"{path}:{number}: malformed FASTA index")
            name, length = fields[0], int(fields[1])
            if not name or name in seen or length <= 0:
                raise PilotError(f"{path}:{number}: duplicate/empty FASTA index row")
            seen.add(name)
            rows.append((name, length))
    if not rows:
        raise PilotError("H1 FASTA index is empty")
    return rows


def read_h1_partitions(
    path: Path | str, dictionary: Mapping[str, int],
) -> list[tuple[Interval, str, int]]:
    """Read only H1-axis native rows and preserve their original line numbers."""
    rows: list[tuple[Interval, str, int]] = []
    with Path(path).open(encoding="ascii") as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                raise PilotError(f"{path}:{number}: native partition lacks an identifier")
            if fields[0] not in dictionary:
                continue
            interval = Interval(fields[0], int(fields[1]), int(fields[2]))
            if interval.end > dictionary[interval.contig]:
                raise PilotError(f"{path}:{number}: H1 partition exceeds staged dictionary")
            rows.append((interval, fields[3], number))
    if not rows:
        raise PilotError("native partition index contains no staged H1 rows")
    return rows


def _close_range(
    ranges: list[dict[str, object]], contig: str,
    rows: list[tuple[Interval, str, int]], ordinal: int,
) -> None:
    first, last = rows[0][0], rows[-1][0]
    ranges.append({
        "range_id": f"r{ordinal:06d}",
        "contig": contig,
        "start": first.start,
        "end": last.end,
        "length": last.end - first.start,
        "partition_count": len(rows),
        "first_partition_line": rows[0][2],
        "last_partition_line": rows[-1][2],
        "first_native_id": rows[0][1],
        "last_native_id": rows[-1][1],
        "query_required": True,
    })


def freeze_range_plan(
    fai_path: Path | str,
    partitions_path: Path | str,
    plan_path: Path | str,
    *,
    selection_id: str,
    target_bp: int = 5_000_000,
    hard_max_bp: int = 20_000_000,
) -> dict[str, object]:
    """Freeze complete native partitions into disjoint/exhaustive bounded ranges."""
    if target_bp <= 0 or hard_max_bp < target_bp:
        raise PilotError("range target/hard maximum is invalid")
    fai_rows = read_fai(fai_path)
    dictionary = dict(fai_rows)
    partitions = read_h1_partitions(partitions_path, dictionary)
    by_contig: dict[str, list[tuple[Interval, str, int]]] = defaultdict(list)
    for row in partitions:
        by_contig[row[0].contig].append(row)

    ranges: list[dict[str, object]] = []
    range_rows: dict[str, list[tuple[Interval, str, int]]] = {}
    ordinal = 0
    for contig, length in fai_rows:
        rows = sorted(by_contig.get(contig, []), key=lambda row: (row[0].start, row[0].end, row[2]))
        if not rows:
            # IMPG has no graph work for an unaligned H1 contig. It still owns
            # coordinate space in the genome-wide mask/PSMC population, where
            # it is represented as non-callable N rather than silently lost.
            for start in range(0, length, target_bp):
                end = min(start + target_bp, length)
                ranges.append({
                    "range_id": f"r{ordinal:06d}", "contig": contig,
                    "start": start, "end": end, "length": end - start,
                    "partition_count": 0, "first_partition_line": None,
                    "last_partition_line": None, "first_native_id": None,
                    "last_native_id": None, "query_required": False,
                })
                range_rows[str(ranges[-1]["range_id"])] = []
                ordinal += 1
            continue
        cursor = 0
        for interval, _, line_number in rows:
            if interval.start != cursor:
                kind = "overlap" if interval.start < cursor else "gap"
                raise PilotError(
                    f"H1 native partitions have a {kind} on {contig} before line {line_number}"
                )
            cursor = interval.end
        if cursor != length:
            raise PilotError(f"H1 native partitions do not reach the end of {contig}")

        current: list[tuple[Interval, str, int]] = []
        for row in rows:
            interval = row[0]
            if interval.length > hard_max_bp:
                raise PilotError(
                    f"single native partition exceeds hard bounded-range maximum: "
                    f"{contig}:{interval.start}-{interval.end}"
                )
            if current and interval.end - current[0][0].start > target_bp:
                _close_range(ranges, contig, current, ordinal)
                range_rows[ranges[-1]["range_id"]] = current  # type: ignore[index]
                ordinal += 1
                current = []
            current.append(row)
        if current:
            _close_range(ranges, contig, current, ordinal)
            range_rows[ranges[-1]["range_id"]] = current  # type: ignore[index]
            ordinal += 1

    maximum = max(int(row["length"]) for row in ranges)
    if maximum > hard_max_bp:
        raise PilotError("frozen range exceeds hard maximum")
    total_bp = sum(int(row["length"]) for row in ranges)
    expected_bp = sum(dictionary.values())
    total_partitions = sum(int(row["partition_count"]) for row in ranges)
    if total_bp != expected_bp or total_partitions != len(partitions):
        raise PilotError("range reduction is not exhaustive")

    result = {
        "schema_version": PLAN_SCHEMA,
        "selection_id": selection_id,
        "coordinate_system": "zero-based half-open H1 reference",
        "range_boundary_policy": "complete consecutive native H1 partition boundaries",
        "partition_owner_policy": "exactly one range per staged-H1 native partition row",
        "variant_owner_policy": "POS belongs to the unique half-open H1 range",
        "callable_owner_policy": "every H1 base belongs to exactly one half-open range",
        "target_bp": target_bp,
        "hard_max_bp": hard_max_bp,
        "range_count": len(ranges),
        "h1_contig_count": len(fai_rows),
        "h1_total_bp": expected_bp,
        "assigned_bp": total_bp,
        "partitioned_h1_bp": sum(
            int(row["length"]) for row in ranges if int(row["partition_count"]) > 0
        ),
        "nonquery_h1_bp": sum(
            int(row["length"]) for row in ranges if int(row["partition_count"]) == 0
        ),
        "native_h1_partition_rows": len(partitions),
        "assigned_native_partition_rows": total_partitions,
        "maximum_range_bp": maximum,
        "minimum_range_bp": min(int(row["length"]) for row in ranges),
        "all_genome_graph_materialized": False,
        "global_impg_lace_permitted": False,
        "global_impg_lace_created": False,
        "range_plan_disjoint": True,
        "range_plan_exhaustive": True,
        "native_partition_one_owner": True,
        "callable_base_one_owner_by_construction": True,
        "fai_path": str(fai_path),
        "fai_sha256": sha256_file(fai_path),
        "partitions_path": str(partitions_path),
        "partitions_sha256": sha256_file(partitions_path),
        "global_partition_assignment_ledger_materialized": False,
        "range_partition_loading": (
            "scan retained native partition index and materialize only the requested range BED"
        ),
        "ranges": ranges,
    }
    Path(plan_path).write_text(canonical_json(result), encoding="utf-8")
    return result


def emit_range_bed(
    plan_path: Path | str, partitions_path: Path | str, range_id: str,
    output_path: Path | str,
) -> dict[str, object]:
    """Load and materialize native rows for exactly one requested range."""
    plan = json.loads(Path(plan_path).read_text(encoding="utf-8"))
    if plan.get("schema_version") != PLAN_SCHEMA:
        raise PilotError("unknown bounded range plan schema")
    matches = [row for row in plan["ranges"] if str(row["range_id"]) == range_id]
    if len(matches) != 1:
        raise PilotError(f"requested range is absent or ambiguous: {range_id}")
    target = matches[0]
    contig, start, end = str(target["contig"]), int(target["start"]), int(target["end"])
    expected = int(target["partition_count"])
    rows: list[tuple[int, int, str, int]] = []
    with Path(partitions_path).open(encoding="ascii") as source:
        for number, line in enumerate(source, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4 or fields[0] != contig:
                continue
            row_start, row_end = int(fields[1]), int(fields[2])
            if row_start < end and row_end > start:
                if row_start < start or row_end > end:
                    raise PilotError("native partition crosses its frozen ownership boundary")
                rows.append((row_start, row_end, fields[3], number))
    rows.sort()
    if len(rows) != expected:
        raise PilotError(
            f"requested range partition census differs from plan: {len(rows)} != {expected}"
        )
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        "".join(
            f"{contig}\t{row_start}\t{row_end}\tquery{ordinal:09d}\t{native_id}\n"
            for ordinal, (row_start, row_end, native_id, _) in enumerate(rows)
        ),
        encoding="ascii",
    )
    return {
        "range_id": range_id, "contig": contig, "start": start, "end": end,
        "partition_count": len(rows), "query_required": bool(rows),
        "only_requested_range_materialized": True, "one_owner_census_passed": True,
    }


def choose_validation_ranges(plan: Mapping[str, object]) -> list[str]:
    ranges = list(plan.get("ranges", []))
    if not ranges:
        raise PilotError("cannot select strata from an empty range plan")
    indices = sorted({0, len(ranges) // 2, len(ranges) - 1})
    return [str(ranges[index]["range_id"]) for index in indices]


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    freeze = sub.add_parser("freeze-plan")
    freeze.add_argument("h1_fai")
    freeze.add_argument("partitions")
    freeze.add_argument("plan_json")
    freeze.add_argument("--selection-id", required=True)
    freeze.add_argument("--target-bp", type=int, default=5_000_000)
    freeze.add_argument("--hard-max-bp", type=int, default=20_000_000)
    emit = sub.add_parser("emit-range-bed")
    emit.add_argument("plan_json")
    emit.add_argument("partitions")
    emit.add_argument("range_id")
    emit.add_argument("output_bed")
    strata = sub.add_parser("validation-ranges")
    strata.add_argument("plan_json")
    args = parser.parse_args(argv)
    try:
        if args.command == "freeze-plan":
            result = freeze_range_plan(
                args.h1_fai, args.partitions, args.plan_json,
                selection_id=args.selection_id, target_bp=args.target_bp,
                hard_max_bp=args.hard_max_bp,
            )
        elif args.command == "emit-range-bed":
            result = emit_range_bed(
                args.plan_json, args.partitions, args.range_id, args.output_bed
            )
        else:
            plan = json.loads(Path(args.plan_json).read_text(encoding="utf-8"))
            result = {"range_ids": choose_validation_ranges(plan)}
    except (PilotError, OSError, ValueError, KeyError, json.JSONDecodeError) as error:
        print(f"ERROR: {error}", file=__import__("sys").stderr)
        return 2
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
