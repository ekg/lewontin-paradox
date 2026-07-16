#!/usr/bin/env python3
"""Select IMPG-native partition rows that intersect annotation target spans."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def read_bed(path: Path, require_name: bool) -> list[tuple[str, int, int, str]]:
    rows = []
    with path.open(encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < (4 if require_name else 3):
                raise SystemExit(f"BED row {line_number} lacks required columns: {path}")
            start, end = int(fields[1]), int(fields[2])
            if not 0 <= start < end:
                raise SystemExit(f"invalid BED interval at row {line_number}: {path}")
            rows.append((fields[0], start, end, fields[3] if len(fields) >= 4 else f"row{line_number}"))
    if not rows:
        raise SystemExit(f"empty BED: {path}")
    return rows


def select_partitions(
    partitions: list[tuple[str, int, int, str]],
    targets: list[tuple[str, int, int, str]],
) -> list[tuple[str, int, int, str, list[str]]]:
    """Intersect native partitions with targets using a per-contig interval sweep."""
    partitions_by_contig: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
    targets_by_contig: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
    for contig, start, end, name in partitions:
        partitions_by_contig[contig].append((start, end, name))
    for contig, start, end, name in targets:
        targets_by_contig[contig].append((start, end, name))

    selected: list[tuple[str, int, int, str, list[str]]] = []
    for contig in sorted(partitions_by_contig):
        contig_targets = sorted(targets_by_contig.get(contig, []))
        first_candidate = 0
        for start, end, partition_id in sorted(partitions_by_contig[contig]):
            while first_candidate < len(contig_targets) and contig_targets[first_candidate][1] <= start:
                first_candidate += 1
            hits = []
            target_index = first_candidate
            while target_index < len(contig_targets) and contig_targets[target_index][0] < end:
                left, right, name = contig_targets[target_index]
                if start < right:
                    hits.append(name)
                target_index += 1
            if hits:
                selected.append((contig, start, end, partition_id, sorted(hits)))
    return selected


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--partitions", type=Path, required=True, help="IMPG partition partitions.bed")
    parser.add_argument("--targets", type=Path, required=True, help="annotation execution spans BED")
    parser.add_argument("--focus-bed", type=Path, required=True)
    parser.add_argument("--mapping-tsv", type=Path, required=True)
    args = parser.parse_args()
    partitions = read_bed(args.partitions, True)
    targets = read_bed(args.targets, True)
    selected = select_partitions(partitions, targets)
    if not selected:
        raise SystemExit("no IMPG-native partition intersects an annotation target")
    args.focus_bed.parent.mkdir(parents=True, exist_ok=True)
    with args.focus_bed.open("w", encoding="utf-8") as handle:
        for contig, start, end, partition_id, _hits in selected:
            handle.write(f"{contig}\t{start}\t{end}\t{partition_id}\n")
    with args.mapping_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["impg_partition_id", "contig", "start_0based", "end_0based_exclusive", "annotation_execution_span_ids"])
        for contig, start, end, partition_id, hits in selected:
            writer.writerow([partition_id, contig, start, end, ",".join(hits)])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
