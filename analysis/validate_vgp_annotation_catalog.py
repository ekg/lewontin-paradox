#!/usr/bin/env python3
"""Independent closed-world validation for the generated VGP annotation catalog."""

from __future__ import annotations

import collections
import csv
import gzip
import hashlib
import json
import os
import sys
from pathlib import Path
from typing import Any, Mapping, Sequence

import catalog_vgp_annotations as cataloger


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def read(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            block = handle.read(8 << 20)
            if not block:
                break
            value.update(block)
    return value.hexdigest()


def validate(repo_root: Path, vgp_root: Path) -> dict[str, Any]:
    analysis = repo_root / "analysis"
    catalog_path = analysis / "vgp_annotation_catalog.tsv"
    summary_path = analysis / "vgp_annotation_catalog_summary.json"
    binding_path = analysis / "vgp_annotation_assembly_bindings.tsv"
    pilot_path = analysis / "vgp_annotation_pilot_bindings.tsv"
    dictionary_path = analysis / "vgp_annotation_sequence_dictionary.tsv.gz"
    mirror_path = analysis / "vgp_freeze1_mirror_manifest.tsv"
    assembly_path = analysis / "vgp_freeze1_bgzf_manifest.tsv"
    metadata_path = vgp_root / "annotations/metadata/ncbi-datasets-v2.json"
    rows = read(catalog_path)
    bindings = read(binding_path)
    pilots = read(pilot_path)
    mirror = read(mirror_path)
    assemblies = read(assembly_path)
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))

    require(len(rows) == 1834, f"catalog row count changed: {len(rows)}")
    require(len(assemblies) == 581, f"assembly closed world changed: {len(assemblies)}")
    require(len(bindings) == 581, f"binding rows changed: {len(bindings)}")
    require(len(pilots) == 10, f"pilot rows changed: {len(pilots)}")
    require(len({row["assembly_accession_version"] for row in bindings}) == 581, "duplicate assembly binding")
    require({row["binding_status"] for row in bindings} == {"EXACT_DICTIONARY"}, "non-exact preferred assembly")
    require({row["binding_status"] for row in pilots} == {"EXACT_DICTIONARY"}, "non-exact pilot binding")

    mirrored_source = {
        row["source_relative_path"]: row
        for row in mirror
        if row.get("object_type") == "file" and cataloger.ANNOTATION_RE.search(row["source_relative_path"])
    }
    mirrored_catalog = {row["source_relative_path"]: row for row in rows if row["origin"] == "MIRRORED_FREEZE1"}
    require(mirrored_source.keys() == mirrored_catalog.keys(), "mirror/catalog annotation path set differs")
    require(len(mirrored_source) == 1481, f"mirror annotation count changed: {len(mirrored_source)}")
    require(
        sum(int(row["observed_bytes"]) for row in mirrored_source.values()) == 10740329039,
        "mirror annotation bytes changed",
    )

    fetched = [row for row in rows if row["origin"] == "FETCHED_NCBI_OFFICIAL"]
    reports = metadata["reports_by_accession"]
    advertised = {accession for accession, report in reports.items() if report.get("annotation_info")}
    require(len(advertised) == 353, f"advertised official annotation count changed: {len(advertised)}")
    require({row["assembly_accession_version"] for row in fetched} == advertised, "official fetch set differs")
    require(len(fetched) == 353, f"fetched official row count changed: {len(fetched)}")
    require(sum(int(row["size_bytes"]) for row in fetched) == 5138794094, "fetched bytes changed")
    require(not summary["closed_world"]["official_fetch_failures"], "official fetch failures present")

    hard_error_fields = [
        "column_errors",
        "coordinate_errors",
        "cds_phase_errors",
        "missing_parent_references",
        "decompression_errors",
    ]
    for row in rows:
        physical = Path(row["physical_path"])
        require(physical.is_file(), f"missing physical annotation {physical}")
        require(physical.stat().st_size == int(row["size_bytes"]), f"size drift for {physical}")
        require(row["digest_verified"] == "true", f"unverified digest for {physical}")
        require(row["parse_validation_status"] == "PASS", f"parse failure for {physical}")
        require(row["species"], f"species missing for {physical}")
        require(row["source_url"], f"source URL missing for {physical}")
        require(row["annotation_provider"], f"provider missing for {physical}")
        for field in hard_error_fields:
            require(int(row[field]) == 0, f"{field}={row[field]} for {physical}")
        if row["origin"] == "FETCHED_NCBI_OFFICIAL":
            require(physical.name == row["actual_sha256"], f"CAS filename/digest mismatch for {physical}")
        for view_text in json.loads(row["view_paths"]):
            view = Path(view_text)
            require(view.is_symlink(), f"annotation view is not a symlink: {view}")
            require(view.resolve(strict=True) == physical.resolve(strict=True), f"view target mismatch: {view}")

    rejected = [row for row in rows if row["accepted"] != "true"]
    require(len(rejected) == 1, f"expected one dictionary exception, observed {len(rejected)}")
    require(rejected[0]["binding_status"] == "DICTIONARY_MISMATCH", "unexpected rejection class")
    require(rejected[0]["assembly_accession_version"] == "GCF_014108235.1", "unexpected rejected assembly")
    require(int(rejected[0]["unresolved_sequence_count"]) == 90, "unexpected unresolved sequence count")
    require(sum(row["accepted"] == "true" for row in rows) == 1833, "accepted count differs")
    require(
        collections.Counter(row["binding_status"] for row in rows)
        == {"EXACT_DICTIONARY": 1828, "VALIDATED_ALIAS": 5, "DICTIONARY_MISMATCH": 1},
        "annotation binding classification changed",
    )

    dictionary_counts: collections.Counter[tuple[str, str]] = collections.Counter()
    with gzip.open(dictionary_path, "rt", encoding="utf-8", newline="") as handle:
        for dictionary_row in csv.DictReader(handle, delimiter="\t"):
            dictionary_counts[(dictionary_row["annotation_object_id"], dictionary_row["physical_path"])] += 1
            if dictionary_row["name_binding"] != "UNRESOLVED":
                require(dictionary_row["assembly_length"], "bound sequence lacks assembly length")
                require(
                    int(dictionary_row["observed_max_end"] or 0)
                    <= 2 * int(dictionary_row["assembly_length"]),
                    "sequence coordinate exceeds circular-aware hard bound",
                )
    for row in rows:
        key = (row["object_id"], row["physical_path"])
        require(dictionary_counts[key] == int(row["observed_sequence_count"]), f"dictionary row loss for {key}")

    for binding in bindings:
        require(Path(binding["annotation_path"]).is_file(), f"preferred annotation missing: {binding}")
        require(Path(binding["assembly_bgzf_path"]).is_file(), f"BGZF missing: {binding}")
        require(Path(binding["assembly_fai_path"]).is_file(), f"FAI missing: {binding}")
    for pilot in pilots:
        require(Path(pilot["annotation_path"]).is_file(), f"pilot annotation missing: {pilot}")
        require(Path(pilot["assembly_bgzf_path"]).is_file(), f"pilot BGZF missing: {pilot}")
        require(Path(pilot["assembly_fai_path"]).is_file(), f"pilot FAI missing: {pilot}")

    shared_manifest_dir = vgp_root / "annotations/manifests"
    expected_shared = [
        catalog_path,
        analysis / "vgp_annotation_catalog.json",
        binding_path,
        pilot_path,
        dictionary_path,
        summary_path,
        analysis / "vgp_annotation_catalog_handoff.md",
    ]
    for source in expected_shared:
        shared = shared_manifest_dir / source.name
        require(shared.is_file(), f"shared manifest missing: {shared}")
        require(digest(shared) == digest(source), f"shared/repository manifest differs: {source.name}")

    return {
        "schema_version": "vgp-annotation-catalog-validation-v1",
        "validated_at_utc": cataloger.utc_now(),
        "status": "PASS",
        "catalog_sha256": digest(catalog_path),
        "catalog_json_sha256": digest(analysis / "vgp_annotation_catalog.json"),
        "sequence_dictionary_sha256": digest(dictionary_path),
        "physical_annotations": len(rows),
        "physical_annotation_bytes": sum(int(row["size_bytes"]) for row in rows),
        "accepted_annotations": sum(row["accepted"] == "true" for row in rows),
        "parsed_annotations": sum(row["parse_validation_status"] == "PASS" for row in rows),
        "mirrored_annotations": len(mirrored_catalog),
        "fetched_official_annotations": len(fetched),
        "assembly_bindings": len(bindings),
        "pilot_bindings": len(pilots),
        "binding_status_counts": dict(sorted(collections.Counter(row["binding_status"] for row in rows).items())),
        "preferred_assembly_status_counts": dict(
            sorted(collections.Counter(row["binding_status"] for row in bindings).items())
        ),
        "pilot_status_counts": dict(sorted(collections.Counter(row["binding_status"] for row in pilots).items())),
        "circular_origin_spanning_sequences": summary["validation_errors"]["circular_coordinate_wrap_count"],
        "rejected_exceptions": [
            {
                "assembly_accession_version": row["assembly_accession_version"],
                "physical_path": row["physical_path"],
                "binding_status": row["binding_status"],
                "binding_reason": row["binding_reason"],
            }
            for row in rejected
        ],
    }


def main(argv: Sequence[str] | None = None) -> int:
    repo_root = Path(__file__).resolve().parent.parent
    vgp_root = Path("/moosefs/erikg/vgp")
    if argv:
        if len(argv) > 1:
            raise SystemExit("usage: validate_vgp_annotation_catalog.py [VGP_ROOT]")
        vgp_root = Path(argv[0])
    report = validate(repo_root, vgp_root)
    output = repo_root / "analysis/vgp_annotation_catalog_validation.json"
    cataloger.json_dump(output, report)
    shared = vgp_root / "annotations/manifests" / output.name
    shared.parent.mkdir(parents=True, exist_ok=True)
    shared.write_bytes(output.read_bytes())
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
