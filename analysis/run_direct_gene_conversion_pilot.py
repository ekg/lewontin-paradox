#!/usr/bin/env python3
"""Fail-closed direct meiotic gene-conversion pilot for D01.

The acquired D01 control contains the complete tetrad/run inventory, exact
TAIR10 reference, and the publication's marker/event/validation supplements.
It deliberately does not contain any FASTQ/BAM payload.  This program audits
and independently reconciles the published directional candidates, but never
promotes them to validated direct events or rates because raw read-backed,
per-tetrad callable, mapping/paralogy, and empirical error evidence is absent.

Repository results may only be written by an authorized Slurm job using the
frozen Guix channel/profile.  Validation is read-only and can run elsewhere.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import re
import sys
from typing import Iterable, Iterator
from xml.etree import ElementTree as ET
from zipfile import BadZipFile, ZipFile


ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = ROOT / "analysis"
DESIGN_DATASETS = ANALYSIS / "gene_conversion_dataset_manifest.tsv"
ACQUISITION = ANALYSIS / "vgp_direct_control_acquisition_manifest.tsv"
ACQUISITION_SUMMARY = ANALYSIS / "vgp_10_pilot_acquisition_summary.json"
GUIX_REALIZATION = ANALYSIS / "guix/vgp_10_pilot/realization.json"

DATASET_OUTPUT = "direct_gene_conversion_dataset_manifest.tsv"
TRACT_OUTPUT = "direct_gene_conversion_tracts.tsv"
SUMMARY_OUTPUT = "direct_gene_conversion_summary.tsv"
REPORT_OUTPUT = "direct_gene_conversion_pilot.md"

MANIFEST_VERSION = "direct-gene-conversion-pilot-v1.0.0"
ACQUISITION_VERSION = "vgp-10-pilot-acquisition-v1.0.0"
CHANNEL_COMMIT = "44bbfc24e4bcc48d0e3343cd3d83452721af8c36"
NOT_ESTIMABLE = "NOT_ESTIMABLE_RAW_EVIDENCE_AND_CALLABILITY_GATE"
PUBLISHED_ONLY = "PUBLISHED_DIRECTIONAL_CANDIDATE_NOT_RAW_REVALIDATED"
TETRADS = ("29", "34", "35", "38", "39", "40", "51", "53", "58", "62", "68", "69", "76")
HIGH_DEPTH_TETRADS = ("29", "38", "40", "58", "62")
PARENTS = ("Col-0", "Ler_WUR", "Cvi_NEW_correct")

M = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
R = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{path}: missing header")
        return reader.fieldnames, list(reader)


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    if not rows:
        raise ValueError(f"refusing empty output: {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _column_index(cell_reference: str) -> int:
    match = re.match(r"[A-Z]+", cell_reference)
    if match is None:
        raise ValueError(f"invalid XLSX cell reference: {cell_reference}")
    value = 0
    for character in match.group(0):
        value = value * 26 + ord(character) - 64
    return value - 1


def xlsx_rows(path: Path, sheet_name: str) -> Iterator[list[str]]:
    """Stream one XLSX sheet as dense rows, using only the Python stdlib."""
    with ZipFile(path) as archive:
        workbook = ET.fromstring(archive.read("xl/workbook.xml"))
        relationships = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        targets = {node.attrib["Id"]: node.attrib["Target"] for node in relationships}
        shared: list[str] = []
        if "xl/sharedStrings.xml" in archive.namelist():
            with archive.open("xl/sharedStrings.xml") as handle:
                for _, node in ET.iterparse(handle, events=("end",)):
                    if node.tag == M + "si":
                        shared.append("".join(part.text or "" for part in node.iter(M + "t")))
                        node.clear()
        sheets = workbook.find(M + "sheets")
        if sheets is None:
            raise ValueError(f"{path}: workbook has no sheets")
        sheet = next((node for node in sheets if node.attrib.get("name") == sheet_name), None)
        if sheet is None:
            raise ValueError(f"{path}: missing sheet {sheet_name}")
        target = "xl/" + targets[sheet.attrib[R + "id"]].lstrip("/")
        with archive.open(target) as handle:
            for _, node in ET.iterparse(handle, events=("end",)):
                if node.tag != M + "row":
                    continue
                values: dict[int, str] = {}
                for cell in node.findall(M + "c"):
                    value_node = cell.find(M + "v")
                    value = "" if value_node is None else value_node.text or ""
                    if cell.attrib.get("t") == "s" and value:
                        value = shared[int(value)]
                    elif cell.attrib.get("t") == "inlineStr":
                        value = "".join(part.text or "" for part in cell.iter(M + "t"))
                    values[_column_index(cell.attrib["r"])] = value
                if values:
                    yield [values.get(index, "") for index in range(max(values) + 1)]
                node.clear()


def table_rows(path: Path, sheet: str, header_row: int) -> list[dict[str, str]]:
    rows = list(xlsx_rows(path, sheet))
    if len(rows) <= header_row:
        raise ValueError(f"{path} sheet {sheet}: no data after header")
    header = rows[header_row]
    if not all(header):
        raise ValueError(f"{path} sheet {sheet}: blank header cell")
    result: list[dict[str, str]] = []
    for values in rows[header_row + 1 :]:
        padded = values + [""] * (len(header) - len(values))
        if any(padded[: len(header)]):
            result.append(dict(zip(header, padded[: len(header)])))
    return result


def _split_semicolon_count(value: str) -> int:
    return len([item for item in value.split(";") if item])


def _split_semicolon_sum(value: str) -> int:
    return sum(int(item) for item in value.split(";") if item)


def _is_sha256(value: str) -> bool:
    return re.fullmatch(r"[0-9a-f]{64}", value) is not None


def _simple_base(value: str) -> bool:
    return value in {"A", "C", "G", "T"}


def _strength(value: str) -> str:
    if value in {"A", "T"}:
        return "W"
    if value in {"C", "G"}:
        return "S"
    return "NA"


def audit_inputs() -> dict[str, object]:
    for path in (DESIGN_DATASETS, ACQUISITION, ACQUISITION_SUMMARY, GUIX_REALIZATION):
        if not path.is_file():
            raise ValueError(f"missing upstream artifact: {path}")

    _, design_rows = read_tsv(DESIGN_DATASETS)
    design = {row["dataset_id"]: row for row in design_rows}
    if design.get("D01", {}).get("execution_state") != "AUTHORIZED_DOWNSTREAM_PILOT":
        raise ValueError("D01 is not authorized by the frozen design")
    if design.get("D02", {}).get("pilot_role") != "alternate":
        raise ValueError("D02 is not the frozen alternate")

    _, acquired = read_tsv(ACQUISITION)
    if len(acquired) != 10 or {row["control_id"] for row in acquired} != {"D01"}:
        raise ValueError("direct-control manifest must contain exactly ten D01 rows")
    if {row["manifest_version"] for row in acquired} != {ACQUISITION_VERSION}:
        raise ValueError("direct-control acquisition version drift")
    by_role = {row["object_role"]: row for row in acquired}
    required_roles = {
        "ena_exact_run_object_manifest",
        "article_version_metadata",
        "published_sample_and_event_tables",
        "published_marker_filter_and_validation_tables",
        "published_methods_supplement",
        "reference_dataset_report",
        "tair10_reference_fasta",
        "tair10_refseq_annotation",
        "tair10_assembly_report",
        "bulk_raw_archive_exclusion_aggregate",
    }
    if set(by_role) != required_roles:
        raise ValueError("direct-control object-role inventory drift")

    for role, row in by_role.items():
        if role == "bulk_raw_archive_exclusion_aggregate":
            if row["status"] != "superseded" or row["local_path"] or row["local_sha256"]:
                raise ValueError("raw archive must remain explicitly superseded and absent")
            continue
        path = Path(row["local_path"])
        if row["status"] != "reused" or not path.is_file():
            raise ValueError(f"{role}: accepted object is not locally present/reused")
        if path.stat().st_size != int(row["observed_bytes"]):
            raise ValueError(f"{role}: byte-size mismatch")
        if not _is_sha256(row["local_sha256"]) or sha256_file(path) != row["local_sha256"]:
            raise ValueError(f"{role}: local SHA-256 mismatch")

    summary = json.loads(ACQUISITION_SUMMARY.read_text(encoding="utf-8"))
    if summary["manifest_version"] != ACQUISITION_VERSION:
        raise ValueError("acquisition summary version drift")
    if summary["direct_control"] != {
        "control_id": "D01",
        "missing_objects": 0,
        "raw_archive_disposition": "superseded_not_selected",
        "selection_status": "selected",
        "verified_or_reused_objects": 9,
    }:
        raise ValueError("direct-control acquisition summary drift")
    if summary["scope_proof"]["ena_fastq_payload_requests"] != 0:
        raise ValueError("unexpected FASTQ payload acquisition")
    if summary["scope_proof"]["slurm_jobs_submitted"] != 0:
        raise ValueError("upstream acquisition unexpectedly submitted Slurm jobs")

    realization = json.loads(GUIX_REALIZATION.read_text(encoding="utf-8"))
    profile = realization.get("profile", "")
    if realization.get("channel_commit") != CHANNEL_COMMIT or not re.fullmatch(r"/gnu/store/[0-9a-z]{32}-profile", profile):
        raise ValueError("frozen Guix realization identity drift")

    runs_path = Path(by_role["ena_exact_run_object_manifest"]["local_path"])
    _, runs = read_tsv(runs_path)
    if len(runs) != 472 or {row["study_accession"] for row in runs} != {"PRJEB4500"}:
        raise ValueError("ENA run inventory drift")
    if {row["secondary_study_accession"] for row in runs} != {"ERP003793"}:
        raise ValueError("ENA secondary study drift")
    if any(row["instrument_model"] != "Illumina HiSeq 2000" for row in runs):
        raise ValueError("read-platform drift")
    fastq_objects = sum(_split_semicolon_count(row["fastq_ftp"]) for row in runs)
    fastq_bytes = sum(_split_semicolon_sum(row["fastq_bytes"]) for row in runs)
    if fastq_objects != 951 or fastq_bytes != 356_984_558_868:
        raise ValueError("ENA FASTQ object accounting drift")
    for row in runs:
        if not (
            _split_semicolon_count(row["fastq_ftp"])
            == _split_semicolon_count(row["fastq_bytes"])
            == _split_semicolon_count(row["fastq_md5"])
        ):
            raise ValueError(f"{row['run_accession']}: FASTQ URL/size/MD5 cardinality mismatch")
        if any(re.fullmatch(r"[0-9a-f]{32}", value) is None for value in row["fastq_md5"].split(";")):
            raise ValueError(f"{row['run_accession']}: invalid provider MD5")

    biological_product_libraries = {f"{tetrad}_{product}" for tetrad in TETRADS for product in range(1, 5)}
    # 29_4_batch4 is an explicitly declared technical library for product 29_4.
    # It belongs in raw-object accounting but is not a 53rd biological product.
    bounded_libraries = biological_product_libraries | set(PARENTS) | {"29_4_batch4"}
    bounded_runs = [row for row in runs if row["library_name"] in bounded_libraries]
    bounded_objects = sum(_split_semicolon_count(row["fastq_ftp"]) for row in bounded_runs)
    bounded_bytes = sum(_split_semicolon_sum(row["fastq_bytes"]) for row in bounded_runs)
    if bounded_objects != 697 or bounded_bytes != 265_809_761_864:
        raise ValueError("bounded 13-tetrad plus parent FASTQ accounting drift")

    event_path = Path(by_role["published_sample_and_event_tables"]["local_path"])
    filter_path = Path(by_role["published_marker_filter_and_validation_tables"]["local_path"])
    sample_rows = table_rows(event_path, "A", 1)
    tetrad_samples = [
        row for row in sample_rows
        if row["Sample type"] in {"Tetrad", "Tetrad (shallow sequenced)"}
    ]
    if len(tetrad_samples) != 52 or {row["Sample"] for row in tetrad_samples} != biological_product_libraries:
        raise ValueError("published sample table does not contain the 52 expected products")
    if any(float(row["Avg. nuclear genome coverage"]) <= 0 for row in tetrad_samples):
        raise ValueError("non-positive published tetrad coverage")
    high_depth = {row["Sample"].split("_")[0] for row in tetrad_samples if row["Sample type"] == "Tetrad"}
    if high_depth != set(HIGH_DEPTH_TETRADS):
        raise ValueError("published high-depth tetrad stratum drift")

    marker_rows = table_rows(filter_path, "B", 1)
    if len(marker_rows) != 137_339:
        raise ValueError("published filtered-marker count drift")
    if any(row["Col-0"] == row["Ler"] for row in marker_rows):
        raise ValueError("filtered parental marker list contains a non-informative Col/Ler state")
    simple_snv_markers = sum(
        _simple_base(row["Col-0"]) and _simple_base(row["Ler"])
        for row in marker_rows
    )
    ws_markers = sum(
        _simple_base(row["Col-0"])
        and _simple_base(row["Ler"])
        and _strength(row["Col-0"]) != _strength(row["Ler"])
        for row in marker_rows
    )

    nuclear_lengths: dict[str, int] = {}
    current = ""
    fasta_path = Path(by_role["tair10_reference_fasta"]["local_path"])
    with gzip.open(fasta_path, "rt", encoding="ascii") as handle:
        for line in handle:
            if line.startswith(">"):
                current = line[1:].split()[0]
                nuclear_lengths[current] = 0
            else:
                nuclear_lengths[current] += len(line.strip())
    nuclear_accessions = {"NC_003070.9", "NC_003071.7", "NC_003074.8", "NC_003075.7", "NC_003076.8"}
    if not nuclear_accessions <= set(nuclear_lengths):
        raise ValueError("TAIR10 nuclear sequence dictionary drift")
    nuclear_bases = sum(nuclear_lengths[name] for name in nuclear_accessions)
    if nuclear_bases != 119_146_348:
        raise ValueError("TAIR10 nuclear span drift")

    co_rows = table_rows(event_path, "F", 2)
    co_marker_rows = table_rows(event_path, "G", 2)
    nco_marker_rows = table_rows(event_path, "D", 2)
    if len(co_rows) != 71:
        raise ValueError("published exact-CO table drift")
    if len({row["CO-ID"] for row in co_rows}) != 71:
        raise ValueError("duplicate exact CO identifier")
    co_candidates = [row for row in co_rows if int(float(row["Minimal length"])) > 0]
    if len(co_candidates) != 44:
        raise ValueError("published CO-associated conversion count drift")
    co_by_id = {row["CO-ID"]: row for row in co_candidates}
    if {row["CO-ID"] for row in co_marker_rows} != set(co_by_id):
        raise ValueError("CO converted-marker/event identifier reconciliation failed")
    for marker in co_marker_rows:
        event = co_by_id[marker["CO-ID"]]
        if {marker["Gamete1"], marker["Gamete2"]} != {event["Gamete 1"], event["Gamete 2"]}:
            raise ValueError(f"CO {marker['CO-ID']}: converted-marker product identity mismatch")
        marker_left, marker_right = sorted((int(marker["Begin"]), int(marker["End"])))
        tract_left, tract_right = sorted((int(event["COCT begin"]), int(event["COCT end"])))
        if marker_left < tract_left or marker_right > tract_right:
            raise ValueError(f"CO {marker['CO-ID']}: converted marker falls outside published COCT")
    nco_ids = sorted({int(row["NCO ID"]) for row in nco_marker_rows})
    if nco_ids != [1, 2, 3, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16]:
        raise ValueError("published NCO identifier set drift")
    primer_rows = table_rows(filter_path, "D", 1)
    if len(primer_rows) != 27:
        raise ValueError("published tetrad NCO primer table drift")
    if not set(nco_ids) <= {int(row["NCO-ID"]) for row in primer_rows}:
        raise ValueError("published accepted NCO set lacks primer-inventory coverage")
    for marker in nco_marker_rows:
        if not int(marker["left end"]) <= int(marker["Locus"]) <= int(marker["right end"]):
            raise ValueError(f"NCO {marker['NCO ID']}: converted marker falls outside published bounds")

    return {
        "design": design,
        "acquired": acquired,
        "by_role": by_role,
        "summary": summary,
        "realization": realization,
        "runs": runs,
        "fastq_objects": fastq_objects,
        "fastq_bytes": fastq_bytes,
        "bounded_runs": len(bounded_runs),
        "bounded_objects": bounded_objects,
        "bounded_bytes": bounded_bytes,
        "sample_rows": sample_rows,
        "marker_rows": marker_rows,
        "simple_snv_markers": simple_snv_markers,
        "ws_markers": ws_markers,
        "nuclear_bases": nuclear_bases,
        "co_rows": co_rows,
        "co_marker_rows": co_marker_rows,
        "nco_marker_rows": nco_marker_rows,
        "primer_rows": primer_rows,
        "input_sha256": {
            "design": sha256_file(DESIGN_DATASETS),
            "acquisition": sha256_file(ACQUISITION),
            "acquisition_summary": sha256_file(ACQUISITION_SUMMARY),
            "run_manifest": sha256_file(runs_path),
            "event_supplement": sha256_file(event_path),
            "filter_supplement": sha256_file(filter_path),
            "reference_fasta": sha256_file(fasta_path),
        },
    }


def dataset_rows(audit: dict[str, object], job_id: str) -> list[dict[str, str]]:
    digest = audit["input_sha256"]
    assert isinstance(digest, dict)
    common = {
        "manifest_version": MANIFEST_VERSION,
        "analysis_slurm_job_id": job_id,
        "guix_channel_commit": CHANNEL_COMMIT,
    }
    return [
        {
            **common,
            "dataset_id": "D01",
            "pilot_role": "primary_executable",
            "activation_status": "SELECTED_BLOCKED_RAW_EVIDENCE_AND_CALLABILITY",
            "dataset": "Arabidopsis recombinant tetrads and doubled haploids",
            "species": "Arabidopsis thaliana",
            "sex_or_gamete": "male meiosis; qrt1 pollen tetrads; Cvi receptor cross",
            "relationship_and_meiosis_structure": "13 complete four-product Col/Ler tetrads; 52 products; 13 independent male meioses; 3 biological parents; 10 DH products excluded from the direct tetrad estimand",
            "parent_of_origin_direction": "ESTABLISHED_PUBLISHED; Col/Ler parental alleles and all four products polarize 2:2 versus 3:1/1:3 segregation",
            "phasing_provenance": "published all-four-product reciprocal haplotype/CO and conversion-marker tables; no independent read-backed phase reconstruction in this run",
            "exact_accessions": "PRJEB4500;ERP003793;doi:10.7554/eLife.01426.v1",
            "assembly_and_annotation": "TAIR10.1;GCF_000001735.4;RefSeq GFF on identical assembly",
            "reference_nuclear_bases": str(audit["nuclear_bases"]),
            "read_provenance": f"ENA exact manifest: 472 runs; {audit['fastq_objects']} FASTQ objects; Illumina HiSeq 2000 WGS; per-object provider MD5 present",
            "read_payload_status": f"NOT_ACQUIRED; all {audit['fastq_objects']} FASTQs/{audit['fastq_bytes']} bytes superseded; bounded tetrad+parent subset {audit['bounded_objects']} FASTQs/{audit['bounded_bytes']} bytes also absent",
            "checksums": ";".join(f"{key}={value}" for key, value in sorted(digest.items())),
            "consent_access_terms": "non-human plant material; human consent not applicable; public ENA study; eLife supplements CC BY 3.0; ENA/EMBL-EBI terms and study citation; third-party rights not transferred",
            "ascertainment": "five high-depth tetrads have exact CO/COCT and NCO tables; eight shallow tetrads have approximate COs and no equivalent NCO discovery/callability table; 137339 common filtered parental markers are not a per-product callable mask",
            "callable_opportunity_status": "NOT_ESTIMABLE; no per-tetrad intersection of parental and four-product depth/genotype/mapping/copy/phase masks",
            "raw_validation_status": "FAIL; no FASTQ/BAM pileups, independent genotypes, sample-swap check, de-novo adjudication, or blinded raw-evidence review",
            "assembly_paralog_status": "FAIL; no candidate-level unique-mapping, repeat, rearrangement, copy-number, segmental-duplication, or paralog mask",
            "error_model_status": "NOT_RUN; empirical genotyping/phasing/mapping errors unavailable; no defensible FDR/FNR spike-in surface",
            "decision": "PUBLISHED_DIRECTIONAL_CANDIDATES_INVENTORIED_EXCLUDED_FROM_DIRECT_RATE",
            "alternate_activation_rule": "D02 not activated: D01 direction is valid and public terms permit analysis, but the versioned acquisition intentionally omitted raw payload; activating D02 requires a pre-result manifest amendment and separately authorized acquisition",
        },
        {
            **common,
            "dataset_id": "D02",
            "pilot_role": "alternate",
            "activation_status": "NOT_ACTIVATED_NO_VERSIONED_AMENDMENT",
            "dataset": "Illumina Platinum Genomes CEPH pedigree 1463",
            "species": "Homo sapiens",
            "sex_or_gamete": "parental transmissions; sex-specific origin would require frozen pedigree/phasing audit",
            "relationship_and_meiosis_structure": "design-only 17-member three-generation pedigree; no local sample or transmission payload",
            "parent_of_origin_direction": "NOT_AUDITED_NOT_ACTIVATED",
            "phasing_provenance": "NOT_ACQUIRED_NOT_USED",
            "exact_accessions": "phs001224;PRJEB3381;PRJEB3246",
            "assembly_and_annotation": "GRCh38 GCA_000001405.15 required if activated; not acquired",
            "reference_nuclear_bases": "NOT_MEASURED",
            "read_provenance": "NOT_ACQUIRED_NOT_USED",
            "read_payload_status": "NOT_ACQUIRED_NOT_USED",
            "checksums": "NOT_APPLICABLE_NOT_ACTIVATED",
            "consent_access_terms": "NOT_REAUDITED; dbGaP/public-access split remains an activation gate",
            "ascertainment": "NOT_MEASURED_NOT_ACTIVATED",
            "callable_opportunity_status": "NOT_MEASURED_NOT_ACTIVATED",
            "raw_validation_status": "NOT_MEASURED_NOT_ACTIVATED",
            "assembly_paralog_status": "NOT_MEASURED_NOT_ACTIVATED",
            "error_model_status": "NOT_MEASURED_NOT_ACTIVATED",
            "decision": "RETAINED_PREAPPROVED_ALTERNATE_ONLY",
            "alternate_activation_rule": "requires D01 intrinsic access/identity/direction/callability failure plus append-only versioned amendment before result unblinding; never an H1/H2 substitute",
        },
    ]


def _group_by(rows: Iterable[dict[str, str]], key: str) -> dict[str, list[dict[str, str]]]:
    grouped: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        grouped.setdefault(row[key], []).append(row)
    return grouped


def _ws_counts(parent_a: str, parent_b: str, resolved: str) -> tuple[int, int, int]:
    if not (_simple_base(parent_a) and _simple_base(parent_b)) or _strength(parent_a) == _strength(parent_b):
        return (0, 0, 0)
    if resolved not in {parent_a, parent_b}:
        return (0, 0, 1)
    return (1 if _strength(resolved) == "S" else 0, 1 if _strength(resolved) == "W" else 0, 0)


def tract_rows(audit: dict[str, object], job_id: str) -> list[dict[str, str]]:
    co_rows = audit["co_rows"]
    co_marker_rows = audit["co_marker_rows"]
    nco_marker_rows = audit["nco_marker_rows"]
    assert isinstance(co_rows, list) and isinstance(co_marker_rows, list) and isinstance(nco_marker_rows, list)
    co_markers = _group_by(co_marker_rows, "CO-ID")
    nco_markers = _group_by(nco_marker_rows, "NCO ID")
    rows: list[dict[str, str]] = []

    for event in co_rows:
        if int(float(event["Minimal length"])) <= 0:
            continue
        markers = co_markers.get(event["CO-ID"], [])
        s_count = w_count = ambiguous = 0
        for marker in markers:
            if marker["Gamete1(Allele)"] != marker["Gamete2(Allele)"]:
                ambiguous += 1
                continue
            if marker["Gamete1(Allele)"] == "Col-0":
                resolved = marker["Col-0-allele"]
            elif marker["Gamete1(Allele)"] == "Ler":
                resolved = marker["Ler-allele"]
            else:
                ambiguous += 1
                continue
            s, w, a = _ws_counts(marker["Col-0-allele"], marker["Ler-allele"], resolved)
            s_count += s
            w_count += w
            ambiguous += a
        tetrad = event["Gamete 1"].split("_")[0]
        rows.append({
            "manifest_version": MANIFEST_VERSION,
            "candidate_id": f"D01_CO_{int(event['CO-ID']):03d}",
            "source_event_id": f"CO-{event['CO-ID']}",
            "tetrad_id": tetrad,
            "meiotic_product_ids": f"{event['Gamete 1']};{event['Gamete 2']}",
            "chromosome": event["Chr"],
            "association_class": "CROSSOVER_ASSOCIATED",
            "direction_evidence": "published reciprocal tetrad haplotypes and same-allele converted marker table",
            "inner_tract_bp": str(int(float(event["Minimal length"]))),
            "outer_tract_bp": str(int(float(event["Maximal length"]))),
            "left_outer_bound": event["Flanking marker upstream"],
            "right_outer_bound": event["Flanking marker downstream"],
            "converted_marker_count": str(len(markers)),
            "informative_ws_marker_count": str(s_count + w_count),
            "s_resolved_marker_count": str(s_count),
            "w_resolved_marker_count": str(w_count),
            "ambiguous_or_non_snv_marker_count": str(ambiguous),
            "published_validation_evidence": "supplement_1F_COCT_bounds_and_1G_converted_markers",
            "independent_computational_review": "PASS_TABLE_RECONCILIATION_ONLY; marker rows contained within published tract bounds and product/tetrad IDs reconciled",
            "raw_read_validation": "NOT_RUN_FASTQ_BAM_ABSENT",
            "mendelian_transmission_validation": "PUBLISHED_ONLY_NOT_RECALLED_FROM_ALL_FOUR_PRODUCTS",
            "paralog_segmental_duplication_mask": "NOT_RUN_EXCLUDED_FROM_ALLELIC_ESTIMATE",
            "mapping_rearrangement_copy_qc": "NOT_RUN_EXCLUDED_FROM_ALLELIC_ESTIMATE",
            "de_novo_and_phase_error_qc": "NOT_RUN_EXCLUDED_FROM_ALLELIC_ESTIMATE",
            "manual_blinded_review": "NOT_RUN_NO_RAW_EVIDENCE_PACKET",
            "analysis_status": PUBLISHED_ONLY,
            "direct_rate_inclusion": "EXCLUDED",
            "exclusion_reason": "NO_RAW_READ_REVALIDATION;NO_PER_TETRAD_CALLABLE_MASK;NO_CANDIDATE_PARALOG_MAPPING_COPY_QC;NO_EMPIRICAL_ERROR_MODEL",
            "analysis_slurm_job_id": job_id,
        })

    for event_id, markers in sorted(nco_markers.items(), key=lambda item: int(item[0])):
        first = markers[0]
        if any(
            marker[field] != first[field]
            for marker in markers[1:]
            for field in ("Tetrads", "Chr", "NCO", "Background", "left end", "right end", "MinNCOsize", "MaxNCOsize")
        ):
            raise ValueError(f"NCO {event_id}: inconsistent tract rows")
        s_count = w_count = ambiguous = 0
        for marker in markers:
            s, w, a = _ws_counts(marker["Background allele"], marker["NCO_allele"], marker["NCO_allele"])
            s_count += s
            w_count += w
            ambiguous += a
        rows.append({
            "manifest_version": MANIFEST_VERSION,
            "candidate_id": f"D01_NCO_{int(event_id):03d}",
            "source_event_id": f"NCO-{event_id}",
            "tetrad_id": first["Tetrads"].split("_")[0],
            "meiotic_product_ids": first["Tetrads"],
            "chromosome": first["Chr"],
            "association_class": "NON_CROSSOVER",
            "direction_evidence": f"published background={first['Background']} and converted={first['NCO']} parental allele in one tetrad product",
            "inner_tract_bp": str(int(float(first["MinNCOsize"]))),
            "outer_tract_bp": str(int(float(first["MaxNCOsize"]))),
            "left_outer_bound": first["left end"],
            "right_outer_bound": first["right end"],
            "converted_marker_count": str(len(markers)),
            "informative_ws_marker_count": str(s_count + w_count),
            "s_resolved_marker_count": str(s_count),
            "w_resolved_marker_count": str(w_count),
            "ambiguous_or_non_snv_marker_count": str(ambiguous),
            "published_validation_evidence": "supplement_1D_NCO_bounds; supplement_2D_candidate-specific_PCR_primer_inventory",
            "independent_computational_review": "PASS_TABLE_RECONCILIATION_ONLY; duplicated-marker rows share product, parental direction, and bounds; primer event IDs inventoried",
            "raw_read_validation": "NOT_RUN_FASTQ_BAM_ABSENT",
            "mendelian_transmission_validation": "PUBLISHED_ONLY_NOT_RECALLED_FROM_ALL_FOUR_PRODUCTS",
            "paralog_segmental_duplication_mask": "NOT_RUN_EXCLUDED_FROM_ALLELIC_ESTIMATE",
            "mapping_rearrangement_copy_qc": "NOT_RUN_EXCLUDED_FROM_ALLELIC_ESTIMATE",
            "de_novo_and_phase_error_qc": "NOT_RUN_EXCLUDED_FROM_ALLELIC_ESTIMATE",
            "manual_blinded_review": "NOT_RUN_NO_RAW_EVIDENCE_PACKET",
            "analysis_status": PUBLISHED_ONLY,
            "direct_rate_inclusion": "EXCLUDED",
            "exclusion_reason": "NO_RAW_READ_REVALIDATION;NO_PER_TETRAD_CALLABLE_MASK;NO_CANDIDATE_PARALOG_MAPPING_COPY_QC;NO_EMPIRICAL_ERROR_MODEL",
            "analysis_slurm_job_id": job_id,
        })
    if len(rows) != 58 or sum(row["association_class"] == "CROSSOVER_ASSOCIATED" for row in rows) != 44:
        raise ValueError("candidate tract accounting drift")
    return rows


def summary_rows(audit: dict[str, object], tracts: list[dict[str, str]], job_id: str) -> list[dict[str, str]]:
    s_count = sum(int(row["s_resolved_marker_count"]) for row in tracts)
    w_count = sum(int(row["w_resolved_marker_count"]) for row in tracts)
    ambiguous = sum(int(row["ambiguous_or_non_snv_marker_count"]) for row in tracts)
    published_marker_opportunities = len(audit["marker_rows"]) * len(TETRADS)  # type: ignore[arg-type]
    rows: list[dict[str, str]] = []

    def add(estimand: str, stratum: str, status: str, numerator: str, denominator: str,
            estimate: str, lower: str, upper: str, unit: str, uncertainty: str, reason: str) -> None:
        rows.append({
            "manifest_version": MANIFEST_VERSION,
            "estimand": estimand,
            "stratum": stratum,
            "status": status,
            "numerator": numerator,
            "denominator": denominator,
            "estimate": estimate,
            "ci_lower": lower,
            "ci_upper": upper,
            "unit": unit,
            "uncertainty_method": uncertainty,
            "reason_or_interpretation": reason,
            "analysis_slurm_job_id": job_id,
        })

    add("PILOT_DISPOSITION", "D01", "NOT_EXECUTED_INPUT_GATE", "0 validated events", "0 callable meioses/bases established", NOT_ESTIMABLE, NOT_ESTIMABLE, NOT_ESTIMABLE, "decision", "not applicable", "direction and complete relationships pass, but raw evidence, callable masks, paralog/mapping/copy QC, and error calibration fail")
    add("RELATIONSHIP_MEIOSES_AUDITED", "all tetrads", "AUDITED_METADATA", "13 complete tetrads;52 products", "13 independent male meioses", "13", "13", "13", "meioses", "exact metadata count", "all four products present for every frozen tetrad; 5 high-depth and 8 shallow")
    add("REFERENCE_NUCLEAR_SPAN", "TAIR10.1", "AUDITED_REFERENCE_NOT_CALLABLE", str(audit["nuclear_bases"]), "five nuclear chromosomes", str(audit["nuclear_bases"]), str(audit["nuclear_bases"]), str(audit["nuclear_bases"]), "bp", "exact sequence count", "reference span is not a callable meiotic denominator")
    add("PUBLISHED_FILTERED_MARKERS", "common study list", "AUDITED_NOT_CALLABLE", str(len(audit["marker_rows"])), "published supplement 2B", str(len(audit["marker_rows"])), str(len(audit["marker_rows"])), str(len(audit["marker_rows"])), "parental marker sites", "exact table count", "no per-product depth/genotype/mapping state; cannot substitute for callable opportunity")
    add("PUBLISHED_WS_MARKERS", "simple SNVs in common study list", "AUDITED_NOT_CALLABLE", str(audit["ws_markers"]), str(audit["simple_snv_markers"]), str(audit["ws_markers"]), str(audit["ws_markers"]), str(audit["ws_markers"]), "W/S parental marker sites", "exact table count", "directionally informative simple-SNV design sites only; indels/complex markers excluded; not converted or callable sites")
    add("DESIGN_MARKER_MEIOSIS_OPPORTUNITIES", "13 tetrads x common marker list", "NOT_A_CALLABLE_DENOMINATOR", str(published_marker_opportunities), "137339 markers x 13 tetrads", str(published_marker_opportunities), str(published_marker_opportunities), str(published_marker_opportunities), "marker-meiosis design cells", "arithmetic audit", "reported solely to prevent accidental use as callable opportunity")
    add("PUBLISHED_DIRECTIONAL_CANDIDATE_TRACTS", "all", "PUBLISHED_ONLY_EXCLUDED", "58", "44 CO-associated + 14 NCO", "58", "58", "58", "published candidate tracts", "exact table reconciliation", "direct direction is present, but no candidate passes independent raw/paralog/error validation")
    add("PUBLISHED_DIRECTIONAL_CANDIDATE_TRACTS", "crossover-associated", "PUBLISHED_ONLY_EXCLUDED", "44", "71 exactly localized high-depth crossovers", "44", "44", "44", "published candidate tracts", "exact table reconciliation", "conversion-marker-containing COCT rows; 71/131 COs are exact and only five high-depth tetrads contribute")
    add("PUBLISHED_DIRECTIONAL_CANDIDATE_TRACTS", "non-crossover", "PUBLISHED_ONLY_EXCLUDED", "14", "five high-depth tetrads", "14", "14", "14", "published candidate tracts", "exact unique event-ID count", "PCR primers are an assay inventory, not raw sequence validation in this run; IDs 4 and 9 are absent from accepted published NCO table")
    add("VALIDATED_ALLELIC_EVENTS", "all", "NOT_ESTIMABLE_INPUT_GATE", "0 admitted", "callable opportunity unavailable", NOT_ESTIMABLE, NOT_ESTIMABLE, NOT_ESTIMABLE, "validated events", "not applicable", "zero is an admitted-analysis count, not evidence that biological events are absent")
    for stratum in ("all", "crossover-associated", "non-crossover"):
        add("D_EVT_RATE_PER_MEIOSIS", stratum, "NOT_ESTIMABLE_INPUT_GATE", NOT_ESTIMABLE, "callable complete meioses unavailable", NOT_ESTIMABLE, NOT_ESTIMABLE, NOT_ESTIMABLE, "validated events/callable meiosis", "exact/cluster interval not computable", "five-versus-eight depth ascertainment and absent raw callability forbid a rate")
        add("D_EVT_RATE_PER_BASE", stratum, "NOT_ESTIMABLE_INPUT_GATE", NOT_ESTIMABLE, "callable base-meiosis unavailable", NOT_ESTIMABLE, NOT_ESTIMABLE, NOT_ESTIMABLE, "validated events/callable bp/meiosis", "Poisson interval not computable", "TAIR10 span and marker list are not callable base-meiosis denominators")
    add("D_TRACT_LENGTH", "all", "NOT_ESTIMABLE_INPUT_GATE", "58 published interval-censored candidates excluded", "0 validated tracts", NOT_ESTIMABLE, NOT_ESTIMABLE, NOT_ESTIMABLE, "bp", "Turnbull/interval-censored distribution not fit", "published inner/outer bounds retained row-wise; validation gate prevents promotion")
    add("D_CO_ASSOCIATION", "near/overlap CO", "NOT_ESTIMABLE_INPUT_GATE", "44 published CO-associated candidates excluded", "matched callable near-CO and non-CO opportunity unavailable", NOT_ESTIMABLE, NOT_ESTIMABLE, NOT_ESTIMABLE, "rate ratio/odds ratio", "permutation/model not run", "no callable matched opportunity or CO uncertainty model across all 13 tetrads")
    add("D_GCBIAS", "published simple W/S candidate markers", "NOT_ESTIMABLE_POWER_AND_INPUT_GATE", f"S={s_count};W={w_count};ambiguous_or_nonSNV={ambiguous}", "0 validated informative event clusters", NOT_ESTIMABLE, NOT_ESTIMABLE, NOT_ESTIMABLE, "Pr(S resolution)", "event-cluster exact/bootstrap interval not computable", "published marker tally is an audit only; linked markers are not independent and reciprocal WS/SW detection is uncalibrated; threshold 100 sites/50 meioses not met")
    add("FALSE_DISCOVERY_RATE", "passing direct-event class", "NOT_ESTIMABLE_INPUT_GATE", "no empirical null simulation", "no raw-derived error/callability map", NOT_ESTIMABLE, NOT_ESTIMABLE, NOT_ESTIMABLE, "FDR", "Mendelian tetrad simulation not run", "genotyping, phasing, mapping, de novo, CO, rearrangement, copy, and paralog errors cannot be parameterized")
    add("FALSE_NEGATIVE_RATE", "tract length x marker density", "NOT_ESTIMABLE_INPUT_GATE", "no recovered spike-ins", "no injected events on empirical raw callability map", NOT_ESTIMABLE, NOT_ESTIMABLE, NOT_ESTIMABLE, "FNR", "spike-in sensitivity not run", "reciprocal WS/SW, tract-length, marker-density, and depth sensitivities unavailable")
    add("SENSITIVITY_ANALYSES", "phasing/mapping/tract/paralog/ascertainment", "NOT_RUN_UPSTREAM_GATE", "0", "five required sensitivity families", NOT_ESTIMABLE, NOT_ESTIMABLE, NOT_ESTIMABLE, "analysis families", "not applicable", "each remains an explicit gate; the five-high-depth-only exact tract ascertainment is documented")
    add("INDEPENDENT_REVIEW", "all 58 published candidates", "TABLE_RECONCILIATION_ONLY", "58 structurally reconciled", "58 published candidates", "58", "58", "58", "candidate rows", "independent parser/table consistency checks", "not blinded/manual raw-evidence review; all remain excluded")
    add("POPULATION_FREQUENCY_GBGC", "out of scope", "NOT_MEASURED", "NOT_MEASURED", "NOT_MEASURED", "NOT_MEASURED", "NOT_MEASURED", "NOT_MEASURED", "population B/SFS", "not applicable", "population-frequency summaries are not direct events")
    add("HISTORICAL_PHYLOGENETIC_GBGC", "out of scope", "NOT_MEASURED", "NOT_MEASURED", "NOT_MEASURED", "NOT_MEASURED", "NOT_MEASURED", "NOT_MEASURED", "branch substitution bias", "not applicable", "phylogenetic substitutions are not direct events")
    add("NON_ALLELIC_CONVERSION", "paralogs/segmental duplications", "NOT_MEASURED", "NOT_MEASURED", "NOT_MEASURED", "NOT_MEASURED", "NOT_MEASURED", "NOT_MEASURED", "non-allelic events", "not applicable", "candidates lacking copy/paralog exclusion cannot enter the allelic estimate")
    add("CROSS_VERTEBRATE_TRANSFER", "VGP sensitivity prior", "NOT_MEASURED", "NOT_MEASURED", "NOT_MEASURED", "NOT_MEASURED", "NOT_MEASURED", "NOT_MEASURED", "bounded sensitivity prior", "not applicable", "no Arabidopsis direct estimate was validated; no cross-species transfer is supplied")
    return rows


def report_text(audit: dict[str, object], tracts: list[dict[str, str]], job_id: str) -> str:
    s_count = sum(int(row["s_resolved_marker_count"]) for row in tracts)
    w_count = sum(int(row["w_resolved_marker_count"]) for row in tracts)
    ambiguous = sum(int(row["ambiguous_or_non_snv_marker_count"]) for row in tracts)
    return f"""# Direct meiotic gene-conversion pilot: fail-closed execution report

**Disposition:** `NOT_EXECUTED_RAW_EVIDENCE_AND_CALLABILITY_GATE`<br>
**Dataset:** D01, *Arabidopsis thaliana*, PRJEB4500 / ERP003793 / eLife 01426 v1<br>
**Execution:** pinned GNU Guix channel `{CHANNEL_COMMIT}`, immutable profile `{audit['realization']['profile']}`, authorized Slurm job `{job_id}`<br>
**Manifest:** `{MANIFEST_VERSION}`

## Outcome

D01 passes the relationship and direction gates: the frozen study contains 13 complete male
Col/Ler four-product tetrads (52 products), exact Col/Ler/Cvi parent libraries, and published
reciprocal haplotype, 2:2, 3:1/1:3, crossover, and converted-marker tables. This is direct
transmission structure, not a trio label, static H1/H2 pair, unpolarized genotype pair, population
frequency, or phylogenetic substitution comparison.

The pilot does **not** produce a validated allelic event rate. The upstream versioned acquisition
deliberately superseded all {audit['fastq_objects']} ENA FASTQ objects ({audit['fastq_bytes']} bytes),
including the bounded 13-tetrad-plus-parent subset of {audit['bounded_objects']} objects
({audit['bounded_bytes']} bytes). It acquired no BAMs or per-product genotype/callability masks.
Consequently no candidate can be re-called from all four products, checked for sample exchange or
de novo mutation, inspected against pileups, or passed through candidate-specific mapping,
rearrangement, copy-number, repeat, segmental-duplication, and paralog filters. No empirical error
surface exists for a defensible Mendelian-null simulation or reciprocal spike-in sensitivity.

The Slurm job therefore performs an independent structural reconciliation of the immutable
published tables and emits a fail-closed ledger. It inventories 58 published directional
candidates—44 crossover-associated and 14 non-crossover—but marks every candidate
`{PUBLISHED_ONLY}` and excludes all 58 from direct-rate, tract-distribution, crossover-association,
and GC-bias estimates. An admitted count of zero is not a biological zero.

## Relationship, sex, direction, and ascertainment audit

- The biological units are 13 male meioses recovered as qrt1 pollen tetrads after a Col/Ler donor
  parent was crossed to a Cvi receptor. Every tetrad ID ({', '.join(TETRADS)}) has products 1–4.
  The read table contains 52 tetrad products, three biological parents, and 10 doubled haploids;
  separate Col/Ler DH and tetrad libraries are technical/library provenance, not extra parents.
- All four products establish reciprocal Col/Ler haplotypes and distinguish 2:2 segregation from
  3:1/1:3 conversion. The published NCO table identifies background and converted parental allele;
  the CO marker table identifies the resolved Col or Ler allele in reciprocal products.
- Exact tract-resolution ascertainment is restricted to five high-depth tetrads
  ({', '.join(HIGH_DEPTH_TETRADS)}). The eight shallow tetrads have approximate crossover rows and no
  equivalent NCO discovery/callability table. Discarding those eight or pretending equal
  sensitivity would bias a per-meiosis rate.
- The ENA snapshot contains 472 runs and {audit['fastq_objects']} provider-MD5-bound FASTQ objects,
  all Illumina HiSeq 2000 genomic WGS. URLs, byte counts, and MD5 cardinalities reconcile, but the
  manifest is provenance, not local raw evidence.
- The material is non-human, so human-participant consent is not applicable. ENA is public and the
  eLife supplements are CC BY 3.0; ENA/EMBL-EBI terms and study citation apply, and public access
  does not transfer third-party rights.

## Assembly and callable opportunity

The exact reference is TAIR10.1 / GCF_000001735.4 with matching RefSeq annotation. The five nuclear
chromosomes contain {audit['nuclear_bases']} reference bases. Supplement 2B has
{len(audit['marker_rows'])} filtered Col/Ler markers, of which {audit['simple_snv_markers']} are
simple SNVs and {audit['ws_markers']} are W/S parental differences. The arithmetic design grid is
{len(audit['marker_rows']) * len(TETRADS)} marker–meiosis cells.

None of those quantities is a callable denominator. Callability requires, for each tetrad, the
intersection of confident biallelic parental state; adequate genotype likelihood/depth in all four
products; separable Cvi contribution; continuous phase; unique diploid/copy-normal mapping; and
repeat, rearrangement, and paralog exclusion. Those per-product states are absent. The nuclear span
cannot be used as callable bases, and the common marker list cannot be multiplied by 13 and used as
callable marker opportunities. Events per meiosis, per base-meiosis, and per callable marker are
therefore `{NOT_ESTIMABLE}` with no confidence interval.

## Candidate tract reconciliation and validation boundary

Supplement 1F contains 71 exactly localized crossovers in the five high-depth tetrads; 44 have a
positive inner converted-tract bound and a corresponding converted-marker inventory. Supplement 1D
contains 14 unique NCO IDs (18 marker rows because four events contain multiple converted markers).
Supplement 2D contains the published tetrad-NCO PCR primer inventory. The parser checks event IDs,
tetrad/product identity, parental direction, inner/outer bounds, marker containment, and duplicated
NCO-row consistency. Inner and outer bounds remain interval-censored; midpoints are not analyzed as
observed lengths.

This is independent computational **table reconciliation only**. It is neither a blinded manual
review nor raw read-backed validation. Candidate-level pileups, all-four-product Mendelian recalls,
copy/paralog masks, and rearrangement/mapping evidence are absent, so the predeclared manual/raw
subset has size zero and all candidates remain excluded. The tract table preserves this status and
reason codes row by row so published candidates cannot inflate the allelic estimate.

## Error model, simulations, and sensitivity

False-discovery and false-negative bounds are not estimable. Genotyping error, depth, allele
balance, phase switches, mapping ambiguity, de novo mutation, crossover uncertainty, structural
rearrangements, copy state, and marker deserts cannot be measured from the selected objects. A
simulation with invented error parameters would not validate this dataset. Mendelian-null
simulation, tract spike-ins, reciprocal W→S/S→W spike-ins, and length/marker-density recovery are
therefore explicit `NOT_RUN_UPSTREAM_GATE` results, as are sensitivities to phasing, mapping, tract
definition, paralog mask, and informative-site ascertainment. The documented five-high-depth versus
eight-shallow split is an observed ascertainment failure, not a corrected sensitivity analysis.

## Direct estimands and GC transmission

No candidate passes the raw-evidence gate, so validated event counts, event rates, interval-censored
tract distributions, and near-CO enrichment are not estimated. Across published candidate rows the
parser can structurally tally simple W/S markers as S={s_count}, W={w_count}, with {ambiguous}
ambiguous/non-SNV marker rows, but this is not a GC-transmission estimate: linked markers are not
independent event clusters, raw genotype support is absent, reciprocal detection is uncalibrated,
and the preregistered threshold of 100 informative converted mismatches across 50 meioses cannot be
met by 13 meioses. `D_GCBIAS` is `NOT_ESTIMABLE_POWER_AND_INPUT_GATE`, never zero bias.

Population-frequency gBGC, historical phylogenetic substitution bias, and non-allelic conversion
are `NOT_MEASURED`. No H1/H2 evidence is used. No Arabidopsis number is transferred as a human,
vertebrate, or VGP rate, and this blocked pilot supplies no bounded cross-species sensitivity prior.

## Alternate and safe continuation

D02 remains the pre-approved alternate but is not activated. D01 is public and directionally valid;
the present failure is the deliberate omission of raw payload and a per-tetrad callable product, not
an intrinsic access or transmission failure. Activating D02 now would require a versioned pre-result
manifest amendment plus a fresh access/consent and acquisition audit. It must never be replaced by a
VGP H1/H2 pair.

A valid D01 continuation must first authorize and checksum-bind the bounded raw objects (or an exact
study BAM/pileup release), construct per-tetrad callable and exclusion masks on GCF_000001735.4,
re-call all four products and parents, apply assembly-aware copy/paralog/rearrangement filters,
predeclare a blinded candidate subset, and submit a separately fingerprinted pinned-Guix Slurm
analysis with empirical null and reciprocal spike-ins. The blocked-output generator must refuse to
run if such inputs later appear; biological execution then requires a new manifest version.
"""


def production_environment_errors() -> list[str]:
    errors: list[str] = []
    if os.environ.get("DIRECT_GC_PINNED_GUIX") != CHANNEL_COMMIT:
        errors.append(f"write requires pinned channel {CHANNEL_COMMIT}")
    if not os.environ.get("SLURM_JOB_ID"):
        errors.append("write and biological preflight require Slurm")
    profile = os.environ.get("GUIX_PROFILE", "")
    if re.fullmatch(r"/gnu/store/[0-9a-z]{32}-profile", profile) is None:
        errors.append("GUIX_PROFILE must be an immutable Guix store profile")
    else:
        realization = json.loads(GUIX_REALIZATION.read_text(encoding="utf-8"))
        if profile != realization.get("profile"):
            errors.append("GUIX_PROFILE differs from the frozen realized profile")
    return errors


def validate_artifacts(directory: Path = ANALYSIS) -> list[str]:
    errors: list[str] = []
    paths = {name: directory / name for name in (DATASET_OUTPUT, TRACT_OUTPUT, SUMMARY_OUTPUT, REPORT_OUTPUT)}
    for name, path in paths.items():
        if not path.is_file() or path.stat().st_size == 0:
            errors.append(f"missing or empty artifact: {name}")
    if errors:
        return errors
    try:
        _, datasets = read_tsv(paths[DATASET_OUTPUT])
        _, tracts = read_tsv(paths[TRACT_OUTPUT])
        _, summary = read_tsv(paths[SUMMARY_OUTPUT])
    except Exception as exc:  # pragma: no cover - defensive diagnostic
        return [str(exc)]

    if len(datasets) != 2 or {row["dataset_id"] for row in datasets} != {"D01", "D02"}:
        errors.append("dataset manifest must contain exactly D01 and D02")
    d01 = next((row for row in datasets if row["dataset_id"] == "D01"), {})
    d02 = next((row for row in datasets if row["dataset_id"] == "D02"), {})
    if d01.get("activation_status") != "SELECTED_BLOCKED_RAW_EVIDENCE_AND_CALLABILITY":
        errors.append("D01 blocked status drift")
    if d02.get("activation_status") != "NOT_ACTIVATED_NO_VERSIONED_AMENDMENT":
        errors.append("D02 was activated without an amendment")
    if len(tracts) != 58:
        errors.append("tract ledger must contain 58 published candidates")
    if sum(row.get("association_class") == "CROSSOVER_ASSOCIATED" for row in tracts) != 44:
        errors.append("CO-associated candidate count drift")
    if sum(row.get("association_class") == "NON_CROSSOVER" for row in tracts) != 14:
        errors.append("NCO candidate count drift")
    if any(row.get("analysis_status") != PUBLISHED_ONLY or row.get("direct_rate_inclusion") != "EXCLUDED" for row in tracts):
        errors.append("a published candidate was promoted to the direct estimate")
    if len({row.get("candidate_id") for row in tracts}) != len(tracts):
        errors.append("duplicate candidate ID")
    if any(int(row["inner_tract_bp"]) <= 0 or int(row["outer_tract_bp"]) < int(row["inner_tract_bp"]) for row in tracts):
        errors.append("invalid interval-censored tract bounds")
    forbidden_measured = {"D_EVT_RATE_PER_MEIOSIS", "D_EVT_RATE_PER_BASE", "D_TRACT_LENGTH", "D_CO_ASSOCIATION", "D_GCBIAS", "FALSE_DISCOVERY_RATE", "FALSE_NEGATIVE_RATE"}
    for row in summary:
        if row.get("estimand") in forbidden_measured and row.get("estimate") not in {NOT_ESTIMABLE, "NOT_MEASURED"}:
            errors.append(f"{row.get('estimand')}: blocked estimand has a numeric estimate")
    required = {
        "PILOT_DISPOSITION", "RELATIONSHIP_MEIOSES_AUDITED", "REFERENCE_NUCLEAR_SPAN",
        "PUBLISHED_FILTERED_MARKERS", "DESIGN_MARKER_MEIOSIS_OPPORTUNITIES",
        "PUBLISHED_DIRECTIONAL_CANDIDATE_TRACTS", "VALIDATED_ALLELIC_EVENTS",
        "D_EVT_RATE_PER_MEIOSIS", "D_EVT_RATE_PER_BASE", "D_TRACT_LENGTH", "D_CO_ASSOCIATION",
        "D_GCBIAS", "FALSE_DISCOVERY_RATE", "FALSE_NEGATIVE_RATE", "SENSITIVITY_ANALYSES",
        "INDEPENDENT_REVIEW", "POPULATION_FREQUENCY_GBGC", "HISTORICAL_PHYLOGENETIC_GBGC",
        "NON_ALLELIC_CONVERSION", "CROSS_VERTEBRATE_TRANSFER",
    }
    missing = required - {row.get("estimand", "") for row in summary}
    if missing:
        errors.append(f"missing summary estimands: {sorted(missing)}")
    for estimand in ("POPULATION_FREQUENCY_GBGC", "HISTORICAL_PHYLOGENETIC_GBGC", "NON_ALLELIC_CONVERSION", "CROSS_VERTEBRATE_TRANSFER"):
        rows = [row for row in summary if row.get("estimand") == estimand]
        if len(rows) != 1 or rows[0].get("status") != "NOT_MEASURED":
            errors.append(f"{estimand}: must be exactly one NOT_MEASURED row")
    report = paths[REPORT_OUTPUT].read_text(encoding="utf-8")
    phrases = (
        "NOT_EXECUTED_RAW_EVIDENCE_AND_CALLABILITY_GATE",
        "13 complete male",
        "44 crossover-associated",
        "14 non-crossover",
        "callable denominator",
        "D02 remains the pre-approved alternate but is not activated",
        "Population-frequency gBGC, historical phylogenetic substitution bias, and non-allelic conversion",
        "No H1/H2 evidence is used",
    )
    for phrase in phrases:
        if phrase not in report:
            errors.append(f"report missing required phrase: {phrase}")
    if "event rate is" in report.lower() or "vertebrate rate" in report.lower() and "No Arabidopsis number" not in report:
        errors.append("report contains an unsafe rate claim")
    return errors


def write_outputs() -> None:
    environment_errors = production_environment_errors()
    if environment_errors:
        raise ValueError("; ".join(environment_errors))
    audit = audit_inputs()
    job_id = os.environ["SLURM_JOB_ID"]
    tracts = tract_rows(audit, job_id)
    write_tsv(ANALYSIS / DATASET_OUTPUT, dataset_rows(audit, job_id))
    write_tsv(ANALYSIS / TRACT_OUTPUT, tracts)
    write_tsv(ANALYSIS / SUMMARY_OUTPUT, summary_rows(audit, tracts, job_id))
    (ANALYSIS / REPORT_OUTPUT).write_text(report_text(audit, tracts, job_id), encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true", help="write the four fail-closed artifacts")
    parser.add_argument("--validate", action="store_true", help="validate the committed/current artifacts")
    args = parser.parse_args(argv)
    try:
        if args.write:
            write_outputs()
        if args.validate:
            errors = validate_artifacts()
            if errors:
                for error in errors:
                    print(f"ERROR: {error}", file=sys.stderr)
                return 1
        if not args.write and not args.validate:
            audit_inputs()
    except (BadZipFile, KeyError, OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
