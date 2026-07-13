#!/usr/bin/env python3
"""Generate, acquire, and validate the checksum-locked Tier 3 manifest."""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import ssl
import tempfile
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence

from analysis.tier3_common import (
    Tier3ValidationError,
    collect_fourfold_sites,
    fasta_dictionary,
    parse_gff,
    read_fasta,
    resolve_contig_aliases,
    verify_file,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SCHEMA = Path(__file__).with_name("schemas") / "tier3_manifest.schema.json"

BUFFALO_COMMIT = "b8f91d5c34675733db8cae8dcab625dcbb55c30a"
BUFFALO_URL = (
    "https://raw.githubusercontent.com/vsbuffalo/paradox_variation/"
    f"{BUFFALO_COMMIT}/data/combined_data.tsv"
)
BUFFALO_SHA256 = "df559451dad94b53ba8675e09811708107a57eeb6ffe8f72b944bcbbf3a1f2eb"
BUFFALO_SIZE_BYTES = 171863
BUFFALO_CORE_COUNT = 173
BUFFALO_CORE_RULE = "diversity_not_NA_and_pred_log10_N_not_NA"

_CREDENTIAL_ASSIGNMENT = re.compile(
    r"(?i)(?:credential|password|passwd|secret|token|access[_-]?token|api[_-]?key)\s*[=:]\s*\S+"
)


def _jsonschema_validate(instance: Mapping[str, Any], schema_path: Path) -> None:
    try:
        import jsonschema
    except ImportError as error:  # pragma: no cover - the Guix manifest supplies it
        raise Tier3ValidationError(
            "python-jsonschema is required; enter the pinned Guix shell"
        ) from error
    schema = json.loads(schema_path.read_text(encoding="utf-8"))
    validator = jsonschema.Draft202012Validator(schema, format_checker=jsonschema.FormatChecker())
    errors = sorted(validator.iter_errors(instance), key=lambda item: list(item.absolute_path))
    if errors:
        error = errors[0]
        location = ".".join(str(part) for part in error.absolute_path) or "<root>"
        raise Tier3ValidationError(f"manifest schema violation at {location}: {error.message}")


def _walk(value: Any) -> Iterable[Any]:
    yield value
    if isinstance(value, Mapping):
        for child in value.values():
            yield from _walk(child)
    elif isinstance(value, list):
        for child in value:
            yield from _walk(child)


def _artifacts(value: Any) -> Iterable[Mapping[str, Any]]:
    for candidate in _walk(value):
        if isinstance(candidate, Mapping) and {"logical_name", "uri", "sha256", "size_bytes"} <= set(candidate):
            yield candidate


def _verify_artifact(artifact: Mapping[str, Any]) -> None:
    uri = urllib.parse.urlparse(artifact["uri"])
    if uri.scheme != "file":
        return
    if uri.netloc not in {"", "localhost"}:
        raise Tier3ValidationError(f"non-local authority in file URI: {artifact['uri']}")
    path = Path(urllib.request.url2pathname(uri.path))
    verify_file(path, artifact["sha256"], artifact["size_bytes"])


def _local_artifact_path(artifact: Optional[Mapping[str, Any]]) -> Optional[Path]:
    if not artifact:
        return None
    uri = urllib.parse.urlparse(artifact["uri"])
    if uri.scheme != "file" or uri.netloc not in {"", "localhost"}:
        return None
    return Path(urllib.request.url2pathname(uri.path))


def _read_contig_mapping(path: Path) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if fields[:2] in (["annotation_contig", "fasta_contig"], ["gff_contig", "fasta_contig"]):
                continue
            if len(fields) != 2 or not all(fields):
                raise Tier3ValidationError(f"invalid contig mapping line {line_number} in {path}")
            if fields[0] in mapping:
                raise Tier3ValidationError(f"duplicate annotation contig {fields[0]!r} in {path}")
            mapping[fields[0]] = fields[1]
    return mapping


def _content_coordinate_audit(dataset: Mapping[str, Any]) -> None:
    """Recompute exact-reference claims whenever the tuple is locally staged."""

    reference, annotation = dataset.get("reference"), dataset.get("annotation")
    if not reference:
        return
    fasta_path = _local_artifact_path(reference["fasta"])
    if fasta_path is None:
        return
    fasta = read_fasta(fasta_path)
    if annotation:
        gff_path = _local_artifact_path(annotation["file"])
        mapping_path = _local_artifact_path(annotation["contig_mapping"])
        if gff_path is not None and mapping_path is not None:
            parsed = parse_gff(gff_path)
            aliases = _read_contig_mapping(mapping_path)
            resolved = resolve_contig_aliases(
                fasta_dictionary(fasta), parsed.sequence_regions, aliases
            )
            annotation_fasta = dict(fasta)
            for annotation_name, fasta_name in resolved.items():
                annotation_fasta[annotation_name] = fasta[fasta_name]
            # Reconstruct every retained canonical CDS, not merely a sample.
            collect_fourfold_sites(annotation_fasta, parsed, annotation["genetic_code"])

    variation = dataset.get("variation")
    if variation and variation["kind"] == "deposited_exact_reference_variants":
        bcf_path = _local_artifact_path(variation["normalized_bcf"])
        if bcf_path is not None:
            try:
                import pysam
            except ImportError as error:  # pragma: no cover - supplied by Guix
                raise Tier3ValidationError("pysam is required for the local REF-allele audit") from error
            with pysam.VariantFile(str(bcf_path)) as records:
                for record in records:
                    if record.contig not in fasta:
                        raise Tier3ValidationError(
                            f"dataset {dataset['dataset_id']}: VCF contig {record.contig!r} absent from FASTA"
                        )
                    observed = fasta[record.contig][record.start : record.start + len(record.ref)].upper()
                    if observed != record.ref.upper():
                        raise Tier3ValidationError(
                            f"dataset {dataset['dataset_id']}: VCF REF mismatch at "
                            f"{record.contig}:{record.pos} ({record.ref} != {observed})"
                        )

    denominator = dataset.get("denominator")
    if denominator:
        bed_path = _local_artifact_path(denominator["mask"])
        if bed_path is not None:
            with bed_path.open("r", encoding="utf-8") as handle:
                for line_number, raw in enumerate(handle, 1):
                    if not raw.strip() or raw.startswith("#"):
                        continue
                    fields = raw.rstrip("\n").split("\t")
                    if len(fields) < 3 or fields[0] not in fasta:
                        raise Tier3ValidationError(f"invalid callable BED contig at line {line_number}")
                    start, end = int(fields[1]), int(fields[2])
                    if not 0 <= start < end <= len(fasta[fields[0]]):
                        raise Tier3ValidationError(f"callable BED out of FASTA bounds at line {line_number}")


def _validate_dataset(dataset: Mapping[str, Any], verify_local_files: bool) -> None:
    dataset_id = dataset["dataset_id"]
    reference = dataset.get("reference")
    annotation = dataset.get("annotation")
    variation = dataset.get("variation")
    denominator = dataset.get("denominator")
    samples = dataset.get("samples")
    observables = dataset["observables"]
    annotation_derived = observables["gc3"] == "eligible" or observables["pi_S_over_pi_W"] == "eligible"
    diversity_eligible = dataset["diversity_eligibility"] == "eligible"

    if annotation_derived or diversity_eligible:
        if not reference or not annotation:
            raise Tier3ValidationError(f"dataset {dataset_id}: reference and annotation are required")
        if annotation["assembly_accession"] != reference["assembly_accession"]:
            raise Tier3ValidationError(
                f"dataset {dataset_id}: annotation assembly {annotation['assembly_accession']} does not match "
                f"reference assembly {reference['assembly_accession']}"
            )
        if annotation["status"] != "native":
            raise Tier3ValidationError(
                f"dataset {dataset_id}: primary annotation-derived results require a native exact-assembly annotation"
            )

    if diversity_eligible:
        if not variation:
            raise Tier3ValidationError(f"dataset {dataset_id}: eligible diversity lacks a variant/alignment source")
        if not denominator:
            raise Tier3ValidationError(f"dataset {dataset_id}: eligible diversity lacks an explicit denominator")
        if not samples:
            raise Tier3ValidationError(f"dataset {dataset_id}: eligible diversity lacks sample/ploidy policy")
        if denominator["sample_list_sha256"] != samples["sample_list"]["sha256"]:
            raise Tier3ValidationError(f"dataset {dataset_id}: denominator and sample-list checksums disagree")
        if variation["kind"] == "deposited_exact_reference_variants" and denominator["kind"] not in {
            "all_sites_vcf", "gvcf", "cohort_callable_mask"
        }:
            raise Tier3ValidationError(
                f"dataset {dataset_id}: deposited variant source lacks explicit genotype callability"
            )
        if variation["kind"] == "direct_wfmash_h2_to_h1" and denominator["kind"] != "h1_reference_alignable_mask":
            raise Tier3ValidationError(
                f"dataset {dataset_id}: direct WFMASH requires an H1-reference alignable denominator"
            )

    if reference and annotation:
        fasta_accession = reference["assembly_accession"]
        gff_accession = annotation["assembly_accession"]
        if fasta_accession != gff_accession:
            raise Tier3ValidationError(f"dataset {dataset_id}: FASTA/GFF assembly-coordinate mismatch")

    if verify_local_files:
        for artifact in _artifacts(dataset):
            _verify_artifact(artifact)
        _content_coordinate_audit(dataset)


def validate_manifest(
    manifest: Mapping[str, Any],
    *,
    schema_path: Path = DEFAULT_SCHEMA,
    verify_local_files: bool = True,
) -> Mapping[str, Any]:
    """Validate JSON Schema plus invariants that JSON Schema cannot express."""

    _jsonschema_validate(manifest, schema_path)
    for value in _walk(manifest):
        if isinstance(value, str) and _CREDENTIAL_ASSIGNMENT.search(value):
            raise Tier3ValidationError("manifest contains a credential-like assignment")
    dataset_ids = [dataset["dataset_id"] for dataset in manifest["datasets"]]
    duplicates = sorted({item for item in dataset_ids if dataset_ids.count(item) > 1})
    if duplicates:
        raise Tier3ValidationError(f"duplicate dataset_id values: {duplicates!r}")
    for dataset in manifest["datasets"]:
        _validate_dataset(dataset, verify_local_files)
    if verify_local_files:
        for artifact in _artifacts(manifest.get("guix", {})):
            _verify_artifact(artifact)
    return manifest


def load_and_validate_manifest(
    path: Path, *, schema_path: Path = DEFAULT_SCHEMA, verify_local_files: bool = True
) -> Mapping[str, Any]:
    try:
        manifest = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise Tier3ValidationError(f"cannot read manifest {path}: {error}") from error
    return validate_manifest(manifest, schema_path=schema_path, verify_local_files=verify_local_files)


def buffalo_core_rows(path: Path) -> List[Dict[str, str]]:
    """Read and verify Buffalo's 173-species source cohort.

    Presence, rather than numeric validity, is the published core criterion.
    This preserves the one row whose non-positive value is excluded only later
    from logarithmic model fitting.
    """

    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"species", "diversity", "pred_log10_N"}
        if reader.fieldnames is None or not required <= set(reader.fieldnames):
            raise Tier3ValidationError(f"Buffalo table lacks columns {sorted(required)!r}")
        core = [
            dict(row)
            for row in reader
            if row["diversity"] not in {"", "NA", "NaN"}
            and row["pred_log10_N"] not in {"", "NA", "NaN"}
        ]
    names = [row["species"].strip() for row in core]
    if len(core) != BUFFALO_CORE_COUNT:
        raise Tier3ValidationError(
            f"Buffalo core expected {BUFFALO_CORE_COUNT} rows under {BUFFALO_CORE_RULE}, observed {len(core)}"
        )
    if len(names) != len(set(names)):
        raise Tier3ValidationError("Buffalo core scientific names are not unique")
    return sorted(core, key=lambda row: row["species"].encode("utf-8"))


def _https_context() -> ssl.SSLContext:
    """Return a strict TLS context rooted in the active Guix environment."""

    configured = os.environ.get("SSL_CERT_DIR")
    if configured:
        return ssl.create_default_context(capath=configured)
    environment = os.environ.get("GUIX_ENVIRONMENT")
    if environment:
        cert_dir = Path(environment) / "etc" / "ssl" / "certs"
        if not cert_dir.is_dir():
            raise Tier3ValidationError(
                "active Guix environment lacks etc/ssl/certs; add nss-certs to the manifest"
            )
        return ssl.create_default_context(capath=str(cert_dir))
    return ssl.create_default_context()


def acquire_buffalo(destination: Path) -> List[Dict[str, str]]:
    """Atomically fetch the immutable Buffalo table and verify it before use."""

    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
    fd, temporary_name = tempfile.mkstemp(prefix=f".{destination.name}.", dir=str(destination.parent))
    os.close(fd)
    temporary = Path(temporary_name)
    try:
        request = urllib.request.Request(BUFFALO_URL, headers={"User-Agent": "tier3-foundations/1"})
        with urllib.request.urlopen(
            request, timeout=120, context=_https_context()
        ) as response, temporary.open("wb") as output:
            while True:
                chunk = response.read(1024 * 1024)
                if not chunk:
                    break
                output.write(chunk)
        verify_file(temporary, BUFFALO_SHA256, BUFFALO_SIZE_BYTES)
        rows = buffalo_core_rows(temporary)
        os.replace(str(temporary), str(destination))
        return rows
    finally:
        if temporary.exists():
            temporary.unlink()


def buffalo_provenance(rows: Sequence[Mapping[str, str]]) -> Dict[str, Any]:
    if len(rows) != BUFFALO_CORE_COUNT:
        raise Tier3ValidationError(f"Buffalo provenance requires all {BUFFALO_CORE_COUNT} core species")
    return {
        "source": {
            "repository": "https://github.com/vsbuffalo/paradox_variation.git",
            "commit": BUFFALO_COMMIT,
            "path": "data/combined_data.tsv",
            "url": BUFFALO_URL,
            "sha256": BUFFALO_SHA256,
            "size_bytes": BUFFALO_SIZE_BYTES,
        },
        "core_rule": BUFFALO_CORE_RULE,
        "core_count": BUFFALO_CORE_COUNT,
        "species": [
            {
                "scientific_name": row["species"],
                "buffalo_diversity": row["diversity"],
                "buffalo_pred_log10_N": row["pred_log10_N"],
            }
            for row in rows
        ],
    }


def write_json_atomic(value: Mapping[str, Any], destination: Path) -> None:
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    text = json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False) + "\n"
    fd, temporary_name = tempfile.mkstemp(prefix=f".{destination.name}.", dir=str(destination.parent), text=True)
    os.close(fd)
    temporary = Path(temporary_name)
    try:
        temporary.write_text(text, encoding="utf-8")
        os.replace(str(temporary), str(destination))
    finally:
        if temporary.exists():
            temporary.unlink()


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subcommands = parser.add_subparsers(dest="command", required=True)
    validate = subcommands.add_parser("validate", help="validate a Tier 3 manifest")
    validate.add_argument("manifest", type=Path)
    validate.add_argument("--no-files", action="store_true", help="skip staged file checksum verification")
    acquire = subcommands.add_parser("acquire-buffalo", help="fetch pinned Buffalo data and emit provenance")
    acquire.add_argument("destination", type=Path, help="raw TSV destination outside git")
    acquire.add_argument("--provenance", required=True, type=Path, help="small provenance JSON output")
    emit = subcommands.add_parser("emit", help="canonicalize and validate an assembled manifest")
    emit.add_argument("source", type=Path)
    emit.add_argument("destination", type=Path)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _parser().parse_args(argv)
    if args.command == "validate":
        load_and_validate_manifest(args.manifest, verify_local_files=not args.no_files)
        print(f"valid: {args.manifest}")
    elif args.command == "acquire-buffalo":
        rows = acquire_buffalo(args.destination)
        write_json_atomic(buffalo_provenance(rows), args.provenance)
        print(f"verified {len(rows)} Buffalo core species: {args.destination}")
    elif args.command == "emit":
        manifest = load_and_validate_manifest(args.source, verify_local_files=True)
        write_json_atomic(manifest, args.destination)
        print(f"wrote validated canonical manifest: {args.destination}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
