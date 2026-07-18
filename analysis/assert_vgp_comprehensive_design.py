#!/usr/bin/env python3
"""Fail-closed assertions for the frozen comprehensive VGP design artifacts.

This validator is read-only and uses only committed design metadata.  If the
external frozen catalog is present, it also verifies its digest, row count,
and the selected catalog-row accessions.  It performs no network request,
payload transfer, environment realization, or scheduler operation.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
PRIMARY = ROOT / "analysis/vgp_10_pair_manifest.tsv"
ALTERNATES = ROOT / "analysis/vgp_10_pair_alternates.tsv"
MANIFEST = ROOT / "analysis/vgp_analysis_manifest.json"
SCHEMA = ROOT / "analysis/schemas/vgp_analysis_manifest.schema.json"
PLAN = ROOT / "analysis/vgp_comprehensive_research_plan.md"
FROZEN_CATALOG = Path(
    "/moosefs/erikg/vgp/manifests/"
    "VGPPhase1-freeze-1.0.commit-dc1b2af5a7741b97d66fb10cb2bce97f41765cdf.tsv"
)

CATALOG_COMMIT = "dc1b2af5a7741b97d66fb10cb2bce97f41765cdf"
CATALOG_SHA256 = "9c58420484a8b76a2d6175b7c26bf709e68bdc726a67fc7541b8c2b5a2fc13a4"
PRIMARY_SHA256 = "bf6c9ff647aed332bfc002bf803e8307203b51432343f2eca6d95a6c80d82997"
ALTERNATE_SHA256 = "55127b18f0f17f6673cc0367c60736207a7e2184198cd3855a7e6ea83f39c52e"

EXPECTED_PRIMARY = {
    "P01": (101, "Camelus dromedarius", "GCA_036321535.1", "GCA_036321565.1"),
    "P02": (66, "Pseudorca crassidens", "GCA_039906515.1", "GCA_039906525.1"),
    "P03": (263, "Colius striatus", "GCA_028858725.2", "GCA_028858625.2"),
    "P04": (254, "Falco naumanni", "GCA_017639655.1", "GCA_017639645.1"),
    "P05": (378, "Candoia aspera", "GCA_035149785.1", "GCA_035125265.1"),
    "P06": (434, "Dendropsophus ebraccatus", "GCA_027789765.1", "GCA_027789725.1"),
    "P07": (670, "Spinachia spinachia", "GCA_048126635.1", "GCA_048127205.1"),
    "P08": (595, "Menidia menidia", "GCA_048628825.1", "GCA_048544195.1"),
    "P09": (676, "Heterodontus francisci", "GCA_036365525.1", "GCA_036365495.1"),
    "P10": (678, "Hemiscyllium ocellatum", "GCA_020745735.1", "GCA_020745765.1"),
}

EXPECTED_ALTERNATES = {
    "A01": (73, "Inia geoffrensis", "GCA_036417435.1", "GCA_036417475.1"),
    "A02": (197, "Lonchura striata domestica", "GCA_046129695.1", "GCA_046129705.1"),
    "A03": (387, "Anolis sagrei", "GCA_037176765.1", "GCA_037176775.1"),
    "A04": (429, "Xenopus petersii", "GCA_038501925.1", "GCA_038501915.1"),
    "A05": (553, "Syngnathus typhle", "GCA_048301445.1", "GCA_048301605.1"),
    "A06": (675, "Hydrolagus colliei", "GCA_035084275.1", "GCA_035084065.1"),
}

REQUIRED_ROW_FIELDS = {
    "selection_id",
    "rank",
    "catalog_row",
    "catalog_commit",
    "catalog_sha256",
    "catalog_taxid",
    "resolved_taxid",
    "species",
    "clade",
    "assembly_generation",
    "genome_size_stratum",
    "expected_diversity_stratum",
    "expected_diversity_evidence",
    "biosample",
    "individual_or_isolate",
    "h1_accession_version",
    "h1_label",
    "h1_role",
    "h1_date",
    "h1_technology",
    "h1_length_bp",
    "h1_contigs",
    "h1_contig_n50_bp",
    "h1_scaffold_n50_bp",
    "h1_report_sha256",
    "h2_accession_version",
    "h2_label",
    "h2_role",
    "h2_date",
    "h2_technology",
    "h2_length_bp",
    "h2_contigs",
    "h2_contig_n50_bp",
    "h2_scaffold_n50_bp",
    "h2_report_sha256",
    "reciprocal_pair_evidence",
    "read_hifi_provenance",
    "long_range_phasing_evidence",
    "qv_evidence",
    "completeness_evidence",
    "duplication_collapse_evidence",
    "native_annotation_accession",
    "annotation_status",
    "source_urls",
    "checksum_contract",
    "license_or_reuse_terms",
    "core_confidence_tier",
    "alternate_replacement_rule",
}

REQUIRED_PLAN_PHRASES = (
    "annotation is not a core gate",
    "not independent validation of the diversity estimate",
    "--num-mappings 1:1",
    "retained query-overlap depth above 1",
    "masked bases become `N`",
    "exactly 200 primary block-bootstrap attempts",
    "Direct pedigree/gamete",
    "Population frequency-spectrum gBGC",
    "Historical phylogenetic",
    "Non-allelic/paralog",
    "NOT_RUN_DESIGN_ONLY",
    "downloaded no biological payload",
    "submitted no local, batch, or Slurm analysis job",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{path}: missing TSV header")
        rows = list(reader)
        return reader.fieldnames, rows


def _validate_rows(
    path: Path,
    rows: list[dict[str, str]],
    expected: dict[str, tuple[int, str, str, str]],
    errors: list[str],
) -> None:
    if len(rows) != len(expected):
        errors.append(f"{path}: observed {len(rows)} rows, expected {len(expected)}")
    by_id = {row.get("selection_id", ""): row for row in rows}
    if set(by_id) != set(expected):
        errors.append(f"{path}: selection IDs {sorted(by_id)!r}, expected {sorted(expected)!r}")
    ranks = [row.get("rank") for row in rows]
    if ranks != [str(i) for i in range(1, len(rows) + 1)]:
        errors.append(f"{path}: ranks are not exactly 1..{len(rows)} in row order")

    accessions: set[str] = set()
    for selection_id, row in by_id.items():
        missing = sorted(field for field in REQUIRED_ROW_FIELDS if not row.get(field, "").strip())
        if missing:
            errors.append(f"{path}:{selection_id}: empty required fields: {', '.join(missing)}")
            continue
        if row["catalog_commit"] != CATALOG_COMMIT or row["catalog_sha256"] != CATALOG_SHA256:
            errors.append(f"{path}:{selection_id}: catalog identity drift")
        if row["catalog_taxid"] != row["resolved_taxid"]:
            errors.append(f"{path}:{selection_id}: unresolved catalog/current TaxId discrepancy")
        if not row["biosample"].startswith(("SAMN", "SAMEA")):
            errors.append(f"{path}:{selection_id}: invalid BioSample {row['biosample']!r}")
        if "each report links the other accession" not in row["reciprocal_pair_evidence"]:
            errors.append(f"{path}:{selection_id}: reciprocal evidence is not explicit")
        if "shared BioSample and isolate" not in row["reciprocal_pair_evidence"]:
            errors.append(f"{path}:{selection_id}: same-individual evidence is not explicit")
        if "before" not in row["checksum_contract"] or "SHA256" not in row["checksum_contract"]:
            errors.append(f"{path}:{selection_id}: incomplete checksum/promotion contract")
        for side in ("h1", "h2"):
            accession = row[f"{side}_accession_version"]
            if not accession.startswith("GCA_") or "." not in accession:
                errors.append(f"{path}:{selection_id}: invalid exact {side} accession {accession!r}")
            if accession in accessions:
                errors.append(f"{path}:{selection_id}: reused assembly accession {accession}")
            accessions.add(accession)
            report_digest = row[f"{side}_report_sha256"]
            if len(report_digest) != 64 or any(c not in "0123456789abcdef" for c in report_digest):
                errors.append(f"{path}:{selection_id}: invalid {side} report SHA-256")
        h1_len = int(row["h1_length_bp"])
        h2_len = int(row["h2_length_bp"])
        ratio = min(h1_len, h2_len) / max(h1_len, h2_len)
        if ratio < 0.8:
            errors.append(f"{path}:{selection_id}: preflight length ratio {ratio:.3f} < 0.8")
        if min(int(row["h1_contig_n50_bp"]), int(row["h2_contig_n50_bp"])) < 1_000_000:
            errors.append(f"{path}:{selection_id}: contig N50 below frozen 1 Mb gate")

        observed = (
            int(row["catalog_row"]),
            row["species"],
            row["h1_accession_version"],
            row["h2_accession_version"],
        )
        if selection_id in expected and observed != expected[selection_id]:
            errors.append(f"{path}:{selection_id}: identity {observed!r}, expected {expected[selection_id]!r}")


def validate_design(root: Path = ROOT, *, validate_external_catalog: bool = True) -> list[str]:
    errors: list[str] = []
    required_paths = [PRIMARY, ALTERNATES, MANIFEST, SCHEMA, PLAN]
    for path in required_paths:
        if not path.is_file():
            errors.append(f"missing required artifact: {path}")
    if errors:
        return errors

    if sha256(PRIMARY) != PRIMARY_SHA256:
        errors.append("primary manifest SHA-256 drift")
    if sha256(ALTERNATES) != ALTERNATE_SHA256:
        errors.append("alternate manifest SHA-256 drift")

    primary_header, primaries = read_tsv(PRIMARY)
    alternate_header, alternates = read_tsv(ALTERNATES)
    if set(primary_header) != REQUIRED_ROW_FIELDS:
        errors.append("primary TSV schema differs from frozen field set")
    if set(alternate_header) != REQUIRED_ROW_FIELDS | {"replaces_clade"}:
        errors.append("alternate TSV schema differs from frozen field set")
    _validate_rows(PRIMARY, primaries, EXPECTED_PRIMARY, errors)
    _validate_rows(ALTERNATES, alternates, EXPECTED_ALTERNATES, errors)

    primary_clades = {row["clade"] for row in primaries}
    required_clades = {"Mammalia", "Aves", "Reptilia", "Amphibia", "Actinopterygii", "Chondrichthyes"}
    if primary_clades != required_clades:
        errors.append(f"primary clade coverage {sorted(primary_clades)!r}, expected {sorted(required_clades)!r}")
    generations = {"early" if row["assembly_generation"].startswith("early") else "later" for row in primaries}
    if generations != {"early", "later"}:
        errors.append("both early and later assembly generations are not represented")
    for field, expected_values in (
        ("genome_size_stratum", {"small", "medium", "large"}),
        ("expected_diversity_stratum", {"low", "medium", "high"}),
    ):
        observed = {row[field] for row in primaries}
        if observed != expected_values:
            errors.append(f"{field} coverage {sorted(observed)!r}, expected {sorted(expected_values)!r}")
    if {row["replaces_clade"] for row in alternates} != required_clades:
        errors.append("alternates do not provide exactly one pre-ranked replacement per required clade")

    manifest: dict[str, Any] = json.loads(MANIFEST.read_text(encoding="utf-8"))
    schema: dict[str, Any] = json.loads(SCHEMA.read_text(encoding="utf-8"))
    try:
        import jsonschema

        jsonschema.Draft202012Validator(schema, format_checker=jsonschema.FormatChecker()).validate(manifest)
    except ImportError:
        # The standard-library validator remains useful on login nodes.  The
        # pinned full suite has jsonschema and exercises the actual JSON Schema.
        pass
    except Exception as exc:  # jsonschema exceptions vary by installed version
        errors.append(f"analysis manifest schema validation failed: {exc}")

    selection = manifest.get("selection", {})
    if selection.get("primary_manifest", {}).get("sha256") != PRIMARY_SHA256:
        errors.append("JSON does not bind the primary manifest digest")
    if selection.get("alternate_manifest", {}).get("sha256") != ALTERNATE_SHA256:
        errors.append("JSON does not bind the alternate manifest digest")
    zero = manifest.get("gates", {}).get("zero_tolerance", {})
    if any(value != 0 for value in zero.values()):
        errors.append("one or more zero-tolerance gates permit a nonzero violation")
    branches = manifest.get("evidence_branches", {})
    if branches.get("population", {}).get("status") != "NOT_RUN_DESIGN_ONLY":
        errors.append("population branch is not frozen design-only")
    if branches.get("non_allelic", {}).get("status") != "NOT_RUN_DESIGN_ONLY":
        errors.append("non-allelic branch is not frozen design-only")
    if branches.get("h1_h2_only_direct_conversion_claim_allowed") is not False:
        errors.append("H1/H2-only direct-conversion claim is not forbidden")
    if branches.get("h1_h2_only_biased_transmission_claim_allowed") is not False:
        errors.append("H1/H2-only biased-transmission claim is not forbidden")

    plan_text = PLAN.read_text(encoding="utf-8")
    for phrase in REQUIRED_PLAN_PHRASES:
        if phrase.casefold() not in plan_text.casefold():
            errors.append(f"research plan missing required phrase {phrase!r}")

    if validate_external_catalog and FROZEN_CATALOG.is_file():
        if sha256(FROZEN_CATALOG) != CATALOG_SHA256:
            errors.append("external frozen catalog SHA-256 drift")
        lines = FROZEN_CATALOG.read_text(encoding="utf-8").splitlines()
        if len(lines) != 717:
            errors.append(f"external frozen catalog has {len(lines)} lines, expected 717")
        catalog_rows = list(csv.DictReader(lines, delimiter="\t"))
        for selection_id, (row_number, species, h1, h2) in {**EXPECTED_PRIMARY, **EXPECTED_ALTERNATES}.items():
            row = catalog_rows[row_number - 1]
            catalog_other = row["Accession #s other high-quality haplotypes"].strip()
            if not catalog_other:
                catalog_other = row["Accession #s alternate haplotypes"].strip()
            if row["Scientific Name"].strip() != species:
                errors.append(f"catalog row {row_number}/{selection_id}: species drift")
            if row["Accession # for main haplotype"].strip() != h1:
                errors.append(f"catalog row {row_number}/{selection_id}: H1 accession drift")
            if h2 not in {item.strip() for item in catalog_other.split(",")}:
                errors.append(f"catalog row {row_number}/{selection_id}: H2 accession drift")

    return errors


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--no-external-catalog",
        action="store_true",
        help="skip the optional in-place frozen catalog cross-check",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    errors = validate_design(validate_external_catalog=not args.no_external_catalog)
    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        return 1
    print(
        "VGP_COMPREHENSIVE_DESIGN_ASSERTIONS_OK "
        "primaries=10 alternates=6 clades=6 generations=2 sizes=3 diversity=3 "
        "catalog=717/716 zero_tolerance=all branches=7"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
