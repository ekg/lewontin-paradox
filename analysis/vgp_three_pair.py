#!/usr/bin/env python3
"""Exact staging and identifier audits for the three-pair VGP reliability pilot.

The P02 incident was caused by a graph reader failing to find sequences that
were demonstrably present in the source FASTA.  This module makes the boundary
explicit: immutable source bytes are verified, records are parsed into one
exact dictionary, uniformly reserialized, and every PAF/partition identifier
is censused against that staged dictionary.  Nothing is dropped to accommodate
a tool's identifier mismatch.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Iterator, Mapping, Sequence

from analysis.vgp_10_pilot import PilotError, canonical_json, sha256_file


def _fasta_records(path: Path) -> Iterator[tuple[str, str]]:
    name: str | None = None
    parts: list[str] = []
    seen: set[str] = set()
    with path.open(encoding="utf-8") as handle:
        for number, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(parts)
                name = line[1:].split()[0]
                if not name or name in seen:
                    raise PilotError(f"{path}:{number}: duplicate/empty FASTA name")
                seen.add(name)
                parts = []
            elif name is None:
                raise PilotError(f"{path}:{number}: sequence before FASTA header")
            else:
                sequence = line.upper()
                if set(sequence) - set("ACGTRYSWKMBDHVN.-"):
                    raise PilotError(f"{path}:{number}: unsupported sequence character")
                parts.append(sequence)
    if name is None:
        raise PilotError(f"empty FASTA: {path}")
    yield name, "".join(parts)


def _materialize_role(asset: Mapping[str, object], output: Path) -> list[dict[str, object]]:
    source = Path(str(asset.get("path", "")))
    expected_size = int(asset.get("size_bytes", 0))
    expected_sha = str(asset.get("sha256", ""))
    expected = asset.get("sequence_dictionary")
    if not source.is_file() or source.stat().st_size != expected_size:
        raise PilotError(f"source FASTA file/size mismatch: {source}")
    if len(expected_sha) != 64 or sha256_file(source) != expected_sha:
        raise PilotError(f"source FASTA SHA256 mismatch: {source}")
    if not isinstance(expected, list) or not expected:
        raise PilotError("source FASTA has no frozen sequence dictionary")

    observed: list[dict[str, object]] = []
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="ascii", newline="\n") as handle:
        for index, (name, sequence) in enumerate(_fasta_records(source)):
            md5 = hashlib.md5(sequence.encode("ascii")).hexdigest()
            row = {
                "name": name,
                "length": len(sequence),
                "md5": md5,
                "sequence_sha256": hashlib.sha256(sequence.encode("ascii")).hexdigest(),
            }
            manifest_row = {key: row[key] for key in ("name", "length", "md5")}
            if index >= len(expected) or manifest_row != expected[index]:
                raise PilotError(f"source dictionary mismatch at {name}")
            observed.append(row)
            handle.write(f">{name}\n")
            for offset in range(0, len(sequence), 80):
                handle.write(sequence[offset:offset + 80] + "\n")
    if len(observed) != len(expected):
        raise PilotError(
            f"source dictionary census mismatch: observed {len(observed)}, expected {len(expected)}"
        )
    return observed


def materialize_exact_staged_fastas(
    manifest_path: Path | str, output_dir: Path | str, audit_path: Path | str,
) -> dict[str, object]:
    """Verify and uniformly stage both haplotypes from their frozen dictionaries."""
    manifest_path, output_dir = Path(manifest_path), Path(output_dir)
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    assets = manifest.get("assets")
    if not isinstance(assets, Mapping):
        raise PilotError("input manifest assets are absent")
    roles: dict[str, object] = {}
    for role, filename in (("h1_fasta", "h1.fa"), ("h2_fasta", "h2.fa")):
        asset = assets.get(role)
        if not isinstance(asset, Mapping):
            raise PilotError(f"input manifest lacks {role}")
        staged = output_dir / filename
        records = _materialize_role(asset, staged)
        roles[role] = {
            "source_path": str(asset["path"]),
            "source_size_bytes": int(asset["size_bytes"]),
            "source_sha256": str(asset["sha256"]),
            "staged_path": str(staged),
            "staged_size_bytes": staged.stat().st_size,
            "staged_sha256": sha256_file(staged),
            "records": records,
            "record_count": len(records),
            "total_bp": sum(int(row["length"]) for row in records),
        }
    h1_rows = roles["h1_fasta"]["records"]  # type: ignore[index]
    universe = output_dir / "h1_universe.bed"
    universe.write_text(
        "".join(f"{row['name']}\t0\t{row['length']}\n" for row in h1_rows),
        encoding="utf-8",
    )
    result = {
        "schema_version": "vgp-three-pair-exact-staged-fasta-v1",
        "selection_id": manifest.get("selection_id"),
        "manifest_path": str(manifest_path),
        "manifest_sha256": sha256_file(manifest_path),
        "normalization": "uppercase sequence, exact first-token identifiers, 80-column wrapping",
        "logical_sequences_equal_to_frozen_source": True,
        "roles": roles,
        "h1_universe_path": str(universe),
        "h1_universe_sha256": sha256_file(universe),
    }
    Path(audit_path).write_text(canonical_json(result), encoding="utf-8")
    return result


def _aliases(
    alias_path: Path | str | None, canonical: Mapping[str, Mapping[str, object]],
) -> dict[str, str]:
    if alias_path is None:
        return {}
    value = json.loads(Path(alias_path).read_text(encoding="utf-8"))
    rows = value.get("aliases")
    if not isinstance(rows, list):
        raise PilotError("alias ledger lacks aliases list")
    result: dict[str, str] = {}
    for row in rows:
        if not isinstance(row, Mapping):
            raise PilotError("malformed alias ledger row")
        observed, target = str(row.get("observed_id", "")), str(row.get("canonical_id", ""))
        target_row = canonical.get(target)
        if not observed or observed in result or target_row is None:
            raise PilotError(f"alias does not resolve uniquely: {observed}")
        if int(row.get("length", -1)) != int(target_row["length"]):
            raise PilotError(f"alias length mismatch: {observed}")
        if row.get("sequence_sha256") != target_row["sequence_sha256"]:
            raise PilotError(f"alias sequence digest mismatch: {observed}")
        if len(str(row.get("source_fasta_sha256", ""))) != 64:
            raise PilotError(f"alias lacks source FASTA digest: {observed}")
        result[observed] = target
    return result


def audit_graph_identifiers(
    dictionary_path: Path | str, paf_path: Path | str, partitions_path: Path | str,
    output_path: Path | str, alias_path: Path | str | None = None,
) -> dict[str, object]:
    """Census graph IDs without omission; permit only sequence-digest aliases."""
    dictionary = json.loads(Path(dictionary_path).read_text(encoding="utf-8"))
    role_values = dictionary.get("roles", {})
    try:
        h1 = {row["name"]: row for row in role_values["h1_fasta"]["records"]}
        h2 = {row["name"]: row for row in role_values["h2_fasta"]["records"]}
    except (KeyError, TypeError) as error:
        raise PilotError("malformed staged FASTA dictionary") from error
    aliases = _aliases(alias_path, {**h1, **h2})

    resolved_aliases: set[str] = set()
    paf_rows = 0
    with Path(paf_path).open(encoding="utf-8") as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                raise PilotError(f"PAF row {number} has fewer than 12 columns")
            query, target = fields[0], fields[5]
            query_resolved, target_resolved = aliases.get(query, query), aliases.get(target, target)
            if query_resolved not in h2 or target_resolved not in h1:
                raise PilotError(f"unresolved graph ID at PAF row {number}: {query}, {target}")
            if int(fields[1]) != int(h2[query_resolved]["length"]):
                raise PilotError(f"PAF query length mismatch at row {number}")
            if int(fields[6]) != int(h1[target_resolved]["length"]):
                raise PilotError(f"PAF target length mismatch at row {number}")
            resolved_aliases.update({item for item in (query, target) if item in aliases})
            paf_rows += 1
    if paf_rows == 0:
        raise PilotError("graph PAF is empty")

    partition_rows = 0
    partition_contigs: set[str] = set()
    with Path(partitions_path).open(encoding="utf-8") as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                raise PilotError(f"partition row {number} has fewer than four columns")
            observed = fields[0]
            resolved = aliases.get(observed, observed)
            if resolved not in h1:
                raise PilotError(f"unresolved partition ID at row {number}: {observed}")
            start, end = int(fields[1]), int(fields[2])
            if start < 0 or end <= start or end > int(h1[resolved]["length"]):
                raise PilotError(f"partition coordinate outside staged FASTA at row {number}")
            if observed in aliases:
                resolved_aliases.add(observed)
            partition_contigs.add(resolved)
            partition_rows += 1
    if partition_rows == 0:
        raise PilotError("graph partitions are empty")
    result = {
        "schema_version": "vgp-three-pair-graph-id-audit-v1",
        "paf_rows": paf_rows,
        "partition_rows": partition_rows,
        "partition_contigs": len(partition_contigs),
        "staged_h1_records": len(h1),
        "staged_h2_records": len(h2),
        "unresolved_ids": 0,
        "digest_validated_aliases_available": len(aliases),
        "digest_validated_aliases_used": sorted(resolved_aliases),
        "silently_omitted_regions": 0,
    }
    Path(output_path).write_text(canonical_json(result), encoding="utf-8")
    return result


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    stage = sub.add_parser("stage-fastas")
    stage.add_argument("manifest")
    stage.add_argument("output_dir")
    stage.add_argument("audit_json")
    audit = sub.add_parser("audit-graph-ids")
    audit.add_argument("dictionary_json")
    audit.add_argument("paf")
    audit.add_argument("partitions")
    audit.add_argument("output_json")
    audit.add_argument("--aliases")
    args = parser.parse_args(argv)
    try:
        if args.command == "stage-fastas":
            result = materialize_exact_staged_fastas(args.manifest, args.output_dir, args.audit_json)
        else:
            result = audit_graph_identifiers(
                args.dictionary_json, args.paf, args.partitions, args.output_json, args.aliases
            )
    except (PilotError, OSError, json.JSONDecodeError, ValueError) as error:
        print(f"ERROR: {error}", file=__import__("sys").stderr)
        return 2
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
