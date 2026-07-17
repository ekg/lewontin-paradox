#!/usr/bin/env python3
"""Freeze the VGP Phase 1 manifest and resolve bounded pilot metadata.

This script is intentionally metadata-only.  It freezes the exact raw TSV
under the external VGP data root, reconciles the synchronized guidance counts,
queries official NCBI metadata for the apparent Tier 3A/Tier 3C candidate
rows, and emits bounded pilot manifests for later review/gating tasks.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import math
import re
import struct
import time
import urllib.parse
import urllib.request
import urllib.error
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, MutableMapping, Optional, Sequence, Tuple


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_ROOT_CONFIG = PROJECT_ROOT / "analysis" / "vgp_data_root_config.json"
DEFAULT_PROVENANCE = PROJECT_ROOT / "analysis" / "vgp_phase1_freeze_provenance.json"
DEFAULT_MANIFEST = PROJECT_ROOT / "analysis" / "vgp_pilot_manifest.tsv"
DEFAULT_REJECTIONS = PROJECT_ROOT / "analysis" / "vgp_pilot_rejections.tsv"
DEFAULT_BUDGET = PROJECT_ROOT / "analysis" / "vgp_pilot_size_budget.tsv"

VGP_REPO_API = "https://api.github.com/repos/VGP/vgp-phase1"
VGP_RAW_TEMPLATE = "https://raw.githubusercontent.com/VGP/vgp-phase1/{commit}/VGPPhase1-freeze-1.0.tsv"
VGP_COMMIT_URL = (
    "https://api.github.com/repos/VGP/vgp-phase1/commits"
    "?path=VGPPhase1-freeze-1.0.tsv&per_page=1"
)
NCBI_POLICY_URL = "https://www.ncbi.nlm.nih.gov/home/about/policies/"
NCBI_DATASET_REPORT = (
    "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{accession}/dataset_report"
)
NCBI_DOWNLOAD_SUMMARY = (
    "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{accession}/download_summary"
    "?include_annotation_type=GENOME_GFF"
)
NCBI_ESEARCH_ASSEMBLY = (
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    "?db=assembly&term={term}&retmode=json"
)
NCBI_ESUMMARY_ASSEMBLY = (
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    "?db=assembly&id={uid}&retmode=json"
)
NCBI_EFETCH_TAXONOMY = (
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    "?db=taxonomy&id={taxid}&retmode=xml"
)

EXPECTED_COUNTS = {
    "unique_species": 714,
    "completed": 223,
    "completed_annotated": 120,
    "triple_eligible": 40,
    "triple_eligible_fish": 13,
    "completed_refseq_fish": 46,
}

ANNOTATION_OK = {"Completed NCBI", "Ready"}

USER_AGENT = "lewontin-paradox-freeze-vgp-manifest/1.0"

SEED_MODALITY_TIER3A = "tier3a_seed"
SEED_MODALITY_TIER3C = "tier3c_seed"

GROUP_MAP = {
    "actinopterygii": "Fishes",
    "sarcopterygii": "Fishes",
    "mammalia": "Mammals",
    "aves": "Birds",
    "reptilia": "Reptiles",
    "amphibia": "Amphibians",
    "chondrichthyes": "Other",
    "hyperoartia": "Other",
    "myxini": "Other",
}

CANDIDATE_SCHEMA_COLUMNS = [
    "inventory_release",
    "record_role",
    "candidate_id",
    "catalog_row_number",
    "scientific_name_source",
    "ncbi_taxid",
    "ncbi_current_name",
    "class",
    "order",
    "source_catalog_url",
    "source_catalog_revision",
    "source_catalog_sha256",
    "source_retrieved_at_utc",
    "catalog_evidence_status",
    "h1_accession_version",
    "h1_release_version_or_date",
    "h1_assembly_name",
    "h1_haplotype_role",
    "h1_fasta_url",
    "h1_provider_md5",
    "h1_fasta_sha256",
    "h1_sequence_set_sha256",
    "h1_length_bp",
    "h1_contig_count",
    "h1_contig_n50_bp",
    "h2_accession_version",
    "h2_release_version_or_date",
    "h2_assembly_name_or_label",
    "h2_haplotype_role",
    "h2_fasta_url",
    "h2_provider_md5",
    "h2_fasta_sha256",
    "h2_sequence_set_sha256",
    "h2_length_bp",
    "h2_contig_count",
    "h2_contig_n50_bp",
    "biosample_accession",
    "individual_or_isolate_id",
    "h1_h2_relationship",
    "haplotype_contig_map_sha256",
    "haplotype_contig_relationship_audit",
    "pair_evidence_url",
    "pair_evidence_retrieved_at_utc",
    "annotation_accession_version",
    "annotation_reference_accession_version",
    "annotation_release_version_or_date",
    "annotation_provider_and_pipeline",
    "annotation_gff_url",
    "annotation_provider_md5",
    "annotation_gff_sha256",
    "annotation_native_status",
    "annotation_contig_map_sha256",
    "annotation_contig_audit",
    "cds_reconstruction_audit",
    "variant_resource_accession",
    "variant_reference_accession_version",
    "variant_url",
    "variant_sha256",
    "callability_resource_accession",
    "callability_reference_accession_version",
    "callability_url",
    "callability_sha256",
    "callable_bases",
    "callable_fraction",
    "queryable_gene_count",
    "queryable_gene_bases",
    "resource_retrieved_at_utc",
    "license_or_reuse_terms",
    "license_evidence_url",
    "evidence_summary",
    "uncertainty_status",
    "assembly_composition_eligible",
    "assembly_diversity_eligible",
    "population_genomic_eligible",
    "demographic_eligible",
    "acceptance_status",
    "explicit_acceptance_or_rejection_reason",
    "blocking_requirement_ids",
]

EXTRA_MANIFEST_COLUMNS = [
    "seed_modalities",
    "lineage_group",
    "manifest_order_code",
    "manifest_family_name",
    "manifest_annotation_status",
    "manifest_refseq_accession",
    "linked_h2_accessions_ncbi",
    "h1_fasta_compressed_bytes",
    "h1_fasta_uncompressed_bytes",
    "h2_fasta_compressed_bytes",
    "h2_fasta_uncompressed_bytes",
    "annotation_gff_compressed_bytes",
    "annotation_gff_uncompressed_bytes",
    "same_individual_evidence",
    "same_individual_status",
    "annotation_file_status",
    "pilot_selection_group",
    "pilot_genome_stratum",
    "pilot_contiguity_stratum",
    "pilot_annotation_stratum",
    "pilot_resource_stratum",
    "pilot_selected",
    "predicted_download_bytes_exact",
    "predicted_persistent_storage_bytes_exact",
    "predicted_core_hours_base",
    "predicted_core_hours_high",
    "predicted_peak_memory_gib_base",
    "predicted_peak_memory_gib_high",
    "predicted_wall_hours_base",
    "predicted_wall_hours_high",
    "predicted_scratch_gb_base",
    "predicted_scratch_gb_high",
    "predicted_inode_count_base",
    "predicted_inode_count_high",
    "predicted_moosefs_read_gb_base",
    "predicted_moosefs_read_gb_high",
    "predicted_moosefs_write_gb_base",
    "predicted_moosefs_write_gb_high",
    "predicted_metadata_operations_base",
    "predicted_metadata_operations_high",
]


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def sha256_bytes(data: bytes) -> str:
    import hashlib

    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    import hashlib

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def slugify(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", value.lower()).strip("_")


def parse_int(value: str) -> Optional[int]:
    text = value.strip().replace(",", "")
    if not text:
        return None
    match = re.search(r"\d+", text)
    if not match:
        return None
    return int(match.group(0))


def parse_float(value: str) -> Optional[float]:
    text = value.strip().replace(",", "")
    if not text:
        return None
    return float(text)


def split_multi(value: str) -> List[str]:
    if not value.strip():
        return []
    parts = [piece.strip() for piece in re.split(r"[;,]", value) if piece.strip()]
    return parts


def number_or_zero(value: Any) -> float:
    if value is None:
        return 0.0
    if isinstance(value, (int, float)):
        return float(value)
    text = str(value).strip()
    if not text or text.startswith("UNRESOLVED"):
        return 0.0
    return float(text.replace(",", ""))


def ftp_to_https(url: str) -> str:
    if url.startswith("ftp://"):
        return "https://" + url[len("ftp://") :]
    return url


def request_url(url: str, *, method: str = "GET", headers: Optional[Mapping[str, str]] = None) -> bytes:
    merged_headers = {"User-Agent": USER_AGENT}
    if headers:
        merged_headers.update(headers)
    last_error: Optional[Exception] = None
    for attempt in range(6):
        request = urllib.request.Request(url, headers=merged_headers, method=method)
        try:
            with urllib.request.urlopen(request, timeout=60) as response:
                return response.read()
        except urllib.error.HTTPError as error:
            last_error = error
            if error.code not in {429, 500, 502, 503, 504}:
                raise
            retry_after = error.headers.get("Retry-After")
            delay = float(retry_after) if retry_after else min(60.0, 1.5 * (attempt + 1))
            time.sleep(delay)
        except Exception as error:
            last_error = error
            if attempt == 5:
                raise
            time.sleep(min(10.0, 1.5 * (attempt + 1)))
    assert last_error is not None
    raise last_error


def json_get(url: str) -> Mapping[str, Any]:
    return json.loads(request_url(url).decode("utf-8"))


def text_get(url: str) -> str:
    return request_url(url).decode("utf-8")


def head_headers(url: str) -> Mapping[str, str]:
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT}, method="HEAD")
    with urllib.request.urlopen(request, timeout=60) as response:
        return dict(response.headers.items())


def gzip_uncompressed_size(url: str) -> int:
    trailer = request_url(url, headers={"Range": "bytes=-4"})
    return struct.unpack("<I", trailer)[0]


def md5_catalog(url: str) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    for raw in text_get(url).splitlines():
        line = raw.strip()
        if not line:
            continue
        checksum, relative = line.split(None, 1)
        mapping[relative.strip().lstrip("./")] = checksum
    return mapping


def taxonomy_lineage(url: str) -> List[Tuple[str, str]]:
    document = ET.fromstring(request_url(url))
    lineage: List[Tuple[str, str]] = []
    for taxon in document.findall(".//LineageEx/Taxon"):
        name = taxon.findtext("ScientificName") or ""
        rank = (taxon.findtext("Rank") or "").lower()
        lineage.append((rank, name))
    return lineage


def lineage_rank(lineage: Sequence[Tuple[str, str]], rank: str) -> Optional[str]:
    for observed_rank, name in lineage:
        if observed_rank == rank:
            return name
    return None


@dataclass(frozen=True)
class SeedRow:
    record_role: str
    candidate_id: str
    catalog_row_number: int
    scientific_name: str
    manifest_lineage: str
    manifest_order_code: str
    manifest_family: str
    manifest_annotation_status: str
    ncbi_taxid_seed: Optional[int]
    h1_accession: str
    h2_accessions: Tuple[str, ...]
    refseq_annotation_accession: str
    seed_modalities: Tuple[str, ...]


@dataclass
class FileMetadata:
    url: str
    provider_md5: Optional[str]
    compressed_bytes: Optional[int]
    uncompressed_bytes: Optional[int]
    last_modified: Optional[str]
    status: str


@dataclass
class AssemblyMetadata:
    accession: str
    uid: str
    report_url: str
    esummary_url: str
    download_summary_url: str
    taxid: int
    scientific_name: str
    current_name: str
    biosample_accession: Optional[str]
    isolate: Optional[str]
    tolid: Optional[str]
    yggdrasil_individual: Optional[str]
    release_date: Optional[str]
    assembly_name: Optional[str]
    diploid_role: Optional[str]
    linked_assemblies: Tuple[str, ...]
    paired_accession: Optional[str]
    refseq_genbank_are_different: Optional[bool]
    paired_differences: Optional[str]
    assembly_type: Optional[str]
    total_sequence_length: Optional[int]
    contig_count: Optional[int]
    contig_n50: Optional[int]
    ftp_path: Optional[str]
    genbank_synonym: Optional[str]
    refseq_synonym: Optional[str]
    synonym_similarity: Optional[str]
    annotation_name: Optional[str]
    annotation_release_date: Optional[str]
    annotation_provider_pipeline: Optional[str]
    annotation_status: Optional[str]
    fasta: FileMetadata
    gff: Optional[FileMetadata]
    taxonomy_lineage: List[Tuple[str, str]]


class Resolver:
    def __init__(self) -> None:
        self._assembly_cache: Dict[Tuple[str, str], AssemblyMetadata] = {}

    def resolve_assembly(self, accession: str, *, mode: str = "primary") -> AssemblyMetadata:
        accession = accession.strip()
        cache_key = (accession, mode)
        if cache_key in self._assembly_cache:
            return self._assembly_cache[cache_key]

        want_report = mode in {"primary", "pair"}
        want_taxonomy = mode == "primary"

        search_term = urllib.parse.quote(f"{accession}[Assembly Accession]")
        esearch_url = NCBI_ESEARCH_ASSEMBLY.format(term=search_term)
        esearch = json_get(esearch_url)
        idlist = esearch["esearchresult"]["idlist"]
        if not idlist:
            raise RuntimeError(f"assembly UID not found for {accession}")
        uid = idlist[0]
        esummary_url = NCBI_ESUMMARY_ASSEMBLY.format(uid=uid)
        summary = json_get(esummary_url)["result"][uid]
        report_url = NCBI_DATASET_REPORT.format(accession=accession)
        report = None
        if want_report:
            report_payload = json_get(report_url)
            report = report_payload["reports"][0] if report_payload.get("reports") else None
        taxid = int((report or {}).get("organism", {}).get("tax_id") or summary.get("taxid") or summary.get("speciestaxid"))
        taxonomy_url = NCBI_EFETCH_TAXONOMY.format(taxid=taxid)
        lineage = taxonomy_lineage(taxonomy_url) if want_taxonomy and taxid else []

        assembly_info = (report or {}).get("assembly_info") or {}
        assembly_stats = (report or {}).get("assembly_stats") or self._summary_stats(summary)
        biosample = assembly_info.get("biosample") or {}
        attr_map = {
            item.get("name"): item.get("value")
            for item in biosample.get("attributes", [])
            if item.get("name") and item.get("value")
        }
        ftp_path = (
            summary.get("ftppath_refseq")
            if accession.startswith("GCF_")
            else summary.get("ftppath_genbank")
        ) or summary.get("ftppath_genbank") or summary.get("ftppath_refseq") or ""
        ftp_https = ftp_to_https(ftp_path) if ftp_path else ""
        basename = ""
        if ftp_path:
            basename = Path(urllib.parse.urlparse(ftp_https).path).name
        md5_map: Dict[str, str] = {}
        if ftp_https:
            md5_map = md5_catalog(f"{ftp_https}/md5checksums.txt")

        fasta_name = f"{basename}_genomic.fna.gz" if basename else ""
        fasta = self._file_metadata(
            f"{ftp_https}/{fasta_name}" if basename else "",
            md5_map.get(fasta_name),
            basename != "",
        )

        gff: Optional[FileMetadata] = None
        annotation_info = (report or {}).get("annotation_info") or {}
        if basename and (accession.startswith("GCF_") or annotation_info):
            gff_name = f"{basename}_genomic.gff.gz"
            gff = self._file_metadata(
                f"{ftp_https}/{gff_name}",
                md5_map.get(gff_name),
                True,
            )

        paired_assembly = assembly_info.get("paired_assembly") or {}
        metadata = AssemblyMetadata(
            accession=accession,
            uid=uid,
            report_url=report_url,
            esummary_url=esummary_url,
            download_summary_url=NCBI_DOWNLOAD_SUMMARY.format(accession=accession) if mode == "primary" else "not_requested",
            taxid=taxid,
            scientific_name=(report or {}).get("organism", {}).get("organism_name") or summary.get("speciesname"),
            current_name=(report or {}).get("organism", {}).get("organism_name") or summary.get("speciesname"),
            biosample_accession=biosample.get("accession") or summary.get("biosampleaccn"),
            isolate=attr_map.get("isolate"),
            tolid=attr_map.get("ToLID"),
            yggdrasil_individual=attr_map.get("yggdrasil_individual"),
            release_date=assembly_info.get("release_date") or summary.get("seqreleasedate") or summary.get("submissiondate"),
            assembly_name=assembly_info.get("assembly_name") or summary.get("assemblyname"),
            diploid_role=assembly_info.get("diploid_role"),
            linked_assemblies=tuple(
                entry["linked_assembly"]
                for entry in assembly_info.get("linked_assemblies", [])
                if entry.get("linked_assembly")
            ),
            paired_accession=paired_assembly.get("accession"),
            refseq_genbank_are_different=paired_assembly.get("refseq_genbank_are_different"),
            paired_differences=paired_assembly.get("differences"),
            assembly_type=assembly_info.get("assembly_type"),
            total_sequence_length=parse_int(str(assembly_stats.get("total_sequence_length", ""))),
            contig_count=parse_int(str(assembly_stats.get("number_of_contigs", ""))),
            contig_n50=parse_int(str(assembly_stats.get("contig_n50", ""))),
            ftp_path=ftp_https or None,
            genbank_synonym=(summary.get("synonym") or {}).get("genbank") or None,
            refseq_synonym=(summary.get("synonym") or {}).get("refseq") or None,
            synonym_similarity=(summary.get("synonym") or {}).get("similarity") or None,
            annotation_name=annotation_info.get("name"),
            annotation_release_date=annotation_info.get("release_date"),
            annotation_provider_pipeline=_annotation_provider_pipeline(annotation_info),
            annotation_status=annotation_info.get("status"),
            fasta=fasta,
            gff=gff,
            taxonomy_lineage=lineage,
        )
        self._assembly_cache[cache_key] = metadata
        return metadata

    @staticmethod
    def _summary_stats(summary: Mapping[str, Any]) -> Mapping[str, Any]:
        meta = summary.get("meta") or ""
        return {
            "total_sequence_length": _meta_stat(meta, "total_length"),
            "number_of_contigs": _meta_stat(meta, "contig_count"),
            "contig_n50": _meta_stat(meta, "contig_n50"),
        }

    def _file_metadata(self, url: str, provider_md5: Optional[str], required: bool) -> FileMetadata:
        if not url:
            return FileMetadata(
                url="UNRESOLVED",
                provider_md5=provider_md5,
                compressed_bytes=None,
                uncompressed_bytes=None,
                last_modified=None,
                status="missing_url",
            )
        try:
            headers = head_headers(url)
            compressed = parse_int(headers.get("Content-Length", "") or "")
            uncompressed = gzip_uncompressed_size(url)
            return FileMetadata(
                url=url,
                provider_md5=provider_md5,
                compressed_bytes=compressed,
                uncompressed_bytes=uncompressed,
                last_modified=headers.get("Last-Modified"),
                status="ok" if compressed is not None and provider_md5 else "missing_checksum_or_size",
            )
        except Exception as error:  # pragma: no cover - exercised in live run
            status = "required_probe_failed" if required else "optional_probe_failed"
            return FileMetadata(
                url=url,
                provider_md5=provider_md5,
                compressed_bytes=None,
                uncompressed_bytes=None,
                last_modified=None,
                status=f"{status}:{error.__class__.__name__}",
            )


def _annotation_provider_pipeline(annotation_info: Mapping[str, Any]) -> Optional[str]:
    if not annotation_info:
        return None
    parts = [
        annotation_info.get("provider"),
        annotation_info.get("pipeline"),
        annotation_info.get("software_version"),
        annotation_info.get("method"),
    ]
    return ";".join(part for part in parts if part)


def _meta_stat(meta: str, category: str) -> Optional[int]:
    match = re.search(
        rf'category="{re.escape(category)}"[^>]*>(\d+)</Stat>',
        meta,
    )
    if not match:
        return None
    return int(match.group(1))


def load_root_config(path: Path) -> Mapping[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def resolve_latest_vgp_commit() -> Mapping[str, str]:
    payload = json_get(VGP_COMMIT_URL)
    commit = payload[0]
    return {
        "commit": commit["sha"],
        "committed_at": commit["commit"]["committer"]["date"],
        "html_url": commit["html_url"],
    }


def resolve_vgp_repo_license() -> Mapping[str, Any]:
    repo = json_get(VGP_REPO_API)
    return {
        "repo_url": repo["html_url"],
        "default_branch": repo["default_branch"],
        "license_spdx_id": (repo.get("license") or {}).get("spdx_id"),
        "license_name": (repo.get("license") or {}).get("name"),
    }


def freeze_raw_manifest(
    *,
    root_config: Mapping[str, Any],
    source_commit: Optional[str] = None,
) -> Mapping[str, Any]:
    resolved = resolve_latest_vgp_commit()
    commit = source_commit or resolved["commit"]
    raw_url = VGP_RAW_TEMPLATE.format(commit=commit)
    retrieved_at = utc_now()
    payload = request_url(raw_url)
    sha256 = sha256_bytes(payload)
    line_count = len(payload.decode("utf-8").splitlines())
    manifests_dir = Path(root_config["root"]) / root_config["layout"]["manifests"]
    manifests_dir.mkdir(parents=True, exist_ok=True)
    frozen_path = manifests_dir / f"VGPPhase1-freeze-1.0.commit-{commit}.tsv"
    frozen_path.write_bytes(payload)
    return {
        "path": str(frozen_path),
        "sha256": sha256,
        "size_bytes": len(payload),
        "line_count": line_count,
        "retrieved_at_utc": retrieved_at,
        "source_commit": commit,
        "source_commit_committed_at_utc": resolved["committed_at"],
        "source_commit_html_url": resolved["html_url"],
        "source_url": raw_url,
        "license": resolve_vgp_repo_license(),
        "text": payload.decode("utf-8"),
    }


def parse_vgp_rows(text: str) -> List[Mapping[str, str]]:
    return list(csv.DictReader(io.StringIO(text), delimiter="\t"))


def build_seed_rows(rows: Sequence[Mapping[str, str]]) -> List[SeedRow]:
    seeds: Dict[str, SeedRow] = {}
    for index, row in enumerate(rows, start=2):
        completed = row["Status"].strip() == "4"
        if not completed:
            continue
        scientific_name = row["Scientific Name"].strip()
        h1_accession = row["Accession # for main haplotype"].strip()
        if not scientific_name or not h1_accession:
            continue
        seed_modalities: List[str] = []
        annotation_status = row["Annotation status"].strip()
        if annotation_status in ANNOTATION_OK and row["Assembly IDs other high-quality haplotypes"].strip():
            seed_modalities.append(SEED_MODALITY_TIER3A)
        if row["Lineage"].strip() == "Fishes" and row["RefSeq annotation main haplotype"].strip():
            seed_modalities.append(SEED_MODALITY_TIER3C)
        if not seed_modalities:
            continue
        key = h1_accession
        merged_modalities = set(seed_modalities)
        if key in seeds:
            merged_modalities.update(seeds[key].seed_modalities)
        seeds[key] = SeedRow(
            record_role="candidate",
            candidate_id=f"{slugify(scientific_name)}_{slugify(h1_accession)}",
            catalog_row_number=index,
            scientific_name=scientific_name,
            manifest_lineage=row["Lineage"].strip(),
            manifest_order_code=row["Order"].strip(),
            manifest_family=row["Family Scientific Name"].strip(),
            manifest_annotation_status=annotation_status,
            ncbi_taxid_seed=parse_int(row["NCBI taxon ID"]),
            h1_accession=h1_accession,
            h2_accessions=tuple(split_multi(row["Accession #s other high-quality haplotypes"])),
            refseq_annotation_accession=row["RefSeq annotation main haplotype"].strip(),
            seed_modalities=tuple(sorted(merged_modalities)),
        )
    return sorted(seeds.values(), key=lambda item: item.catalog_row_number)


def compute_counts(rows: Sequence[Mapping[str, str]]) -> Mapping[str, Any]:
    completed = [row for row in rows if row["Status"].strip() == "4"]
    completed_annotated = [
        row for row in completed if row["Annotation status"].strip() in ANNOTATION_OK
    ]
    triple = [
        row
        for row in completed_annotated
        if row["Assembly IDs other high-quality haplotypes"].strip()
    ]
    fish_triple = [row for row in triple if row["Lineage"].strip() == "Fishes"]
    fish_completed_refseq = [
        row
        for row in completed
        if row["Lineage"].strip() == "Fishes" and row["RefSeq annotation main haplotype"].strip()
    ]
    observed = {
        "unique_species": len({row["Scientific Name"].strip() for row in rows if row["Scientific Name"].strip()}),
        "completed": len(completed),
        "completed_annotated": len(completed_annotated),
        "triple_eligible": len(triple),
        "triple_eligible_fish": len(fish_triple),
        "completed_refseq_fish": len(fish_completed_refseq),
    }
    discrepancies = []
    for key, expected in EXPECTED_COUNTS.items():
        if observed[key] != expected:
            discrepancies.append(
                {
                    "metric": key,
                    "expected": expected,
                    "observed": observed[key],
                    "delta": observed[key] - expected,
                }
            )
    return {
        "observed": observed,
        "discrepancies": discrepancies,
        "current_fish_triple_species": sorted({row["Scientific Name"] for row in fish_triple}),
        "current_fish_completed_refseq_species": sorted(
            {row["Scientific Name"] for row in fish_completed_refseq}
        ),
    }


def resolve_pair(
    seed: SeedRow,
    h1: AssemblyMetadata,
    resolver: Resolver,
) -> Tuple[Optional[AssemblyMetadata], str, str]:
    if SEED_MODALITY_TIER3A not in seed.seed_modalities:
        return None, "not_applicable", "not_applicable"
    linked = set(h1.linked_assemblies)
    seed_h2 = [accession for accession in seed.h2_accessions if accession]
    linked_match = [accession for accession in seed_h2 if accession in linked]
    chosen: Optional[str] = None
    audit = ""
    if len(linked_match) == 1:
        chosen = linked_match[0]
        audit = "manifest_h2_matches_ncbi_linked_assembly"
    elif len(linked) == 1 and len(seed_h2) == 1 and next(iter(linked)) == seed_h2[0]:
        chosen = seed_h2[0]
        audit = "single_manifest_h2_matches_single_ncbi_link"
    elif len(seed_h2) == 1 and not linked:
        chosen = seed_h2[0]
        audit = "manifest_single_h2_without_ncbi_link"
    elif len(seed_h2) == 1 and len(linked) == 1 and next(iter(linked)) != seed_h2[0]:
        return None, "reject_pair_mismatch", "manifest_h2_and_ncbi_linked_h2_disagree"
    else:
        return None, "reject_pair_unresolved", "unable_to_resolve_single_h2_accession"
    h2 = resolver.resolve_assembly(chosen, mode="pair")
    evidence = []
    if h1.biosample_accession and h2.biosample_accession and h1.biosample_accession == h2.biosample_accession:
        evidence.append(f"shared_biosample:{h1.biosample_accession}")
    if h1.isolate and h2.isolate and h1.isolate == h2.isolate:
        evidence.append(f"shared_isolate:{h1.isolate}")
    if h1.tolid and h2.tolid and h1.tolid == h2.tolid:
        evidence.append(f"shared_tolid:{h1.tolid}")
    if (
        h1.yggdrasil_individual
        and h2.yggdrasil_individual
        and h1.yggdrasil_individual == h2.yggdrasil_individual
    ):
        evidence.append(f"shared_yggdrasil_individual:{h1.yggdrasil_individual}")
    if h1.accession in h2.linked_assemblies:
        evidence.append("reciprocal_linked_assemblies")
    if not evidence:
        return h2, "reject_same_individual_unresolved", audit
    return h2, "ok", audit + ";" + ";".join(evidence)


def annotation_bundle(
    seed: SeedRow,
    h1: AssemblyMetadata,
    resolver: Resolver,
) -> Tuple[Optional[AssemblyMetadata], str]:
    if h1.gff and h1.annotation_name:
        return h1, "native_exact_h1_full_annotation"
    refseq = seed.refseq_annotation_accession.strip()
    if not refseq.startswith("GCF_"):
        return None, "missing_refseq_annotation_accession"
    annotation = resolver.resolve_assembly(refseq, mode="annotation")
    if not annotation.gff:
        return None, "refseq_annotation_gff_unavailable"
    if annotation.paired_accession and annotation.paired_accession != h1.accession:
        return None, "refseq_paired_accession_mismatch"
    if annotation.genbank_synonym and annotation.genbank_synonym != h1.accession:
        return None, "refseq_genbank_synonym_mismatch"
    if annotation.refseq_genbank_are_different:
        return None, f"refseq_genbank_different:{annotation.paired_differences or 'unspecified'}"
    if annotation.synonym_similarity and annotation.synonym_similarity != "identical":
        return None, f"refseq_similarity_{annotation.synonym_similarity.lower()}"
    if (
        annotation.total_sequence_length != h1.total_sequence_length
        or annotation.contig_count != h1.contig_count
    ):
        return None, "refseq_h1_assembly_stats_mismatch"
    return annotation, "refseq_identical_to_h1_pending_contig_audit"


def normalize_group(class_name: Optional[str], manifest_lineage: str) -> str:
    if class_name:
        normalized = GROUP_MAP.get(class_name.lower())
        if normalized:
            return normalized
    if manifest_lineage == "Fishes":
        return "Fishes"
    if manifest_lineage in {"Mammals", "Birds", "Reptiles", "Amphibians"}:
        return manifest_lineage
    return "Other"


def selection_categories(rows: Sequence[MutableMapping[str, Any]], key: str) -> None:
    numeric_rows = [row for row in rows if row.get(key) is not None]
    values = sorted(float(row[key]) for row in numeric_rows)
    if not values:
        return
    first = values[max(0, len(values) // 3 - 1)]
    second = values[max(0, 2 * len(values) // 3 - 1)]
    for row in rows:
        value = row.get(key)
        if value is None:
            row[f"{key}_stratum"] = "unknown"
        elif float(value) <= first:
            row[f"{key}_stratum"] = "small"
        elif float(value) <= second:
            row[f"{key}_stratum"] = "medium"
        else:
            row[f"{key}_stratum"] = "large"


def predicted_scale(row: Mapping[str, Any], *, include_h2: bool) -> float:
    h1_length = number_or_zero(row.get("h1_length_bp"))
    h2_length = number_or_zero(row.get("h2_length_bp")) if include_h2 else 0.0
    contigs = max(number_or_zero(row.get("h1_contig_count")), number_or_zero(row.get("h2_contig_count")), 1.0)
    annotation_gb = (number_or_zero(row.get("annotation_gff_uncompressed_bytes")) / 1_000_000_000) or 0.05
    genome_factor = max(0.35, (h1_length + h2_length) / 1_200_000_000)
    contig_factor = max(0.5, math.sqrt(contigs / 500.0))
    annotation_factor = max(0.35, annotation_gb / 0.60)
    return max(0.35, 0.55 * genome_factor + 0.20 * contig_factor + 0.25 * annotation_factor)


def predicted_resources(row: MutableMapping[str, Any]) -> None:
    composition_scale = predicted_scale(row, include_h2=False)
    diversity_scale = predicted_scale(row, include_h2=True)
    has_diversity = row["assembly_diversity_eligible"] == "yes"
    exact_download = sum(
        int(value)
        for value in [
            row.get("h1_fasta_compressed_bytes"),
            row.get("h2_fasta_compressed_bytes"),
            row.get("annotation_gff_compressed_bytes"),
        ]
        if value is not None and not str(value).startswith("UNRESOLVED")
    )
    row["predicted_download_bytes_exact"] = exact_download
    row["predicted_persistent_storage_bytes_exact"] = exact_download
    if has_diversity:
        base_core = round((0.10 + 0.05 + 1.40 + 0.40 + 0.15) * diversity_scale, 4)
        high_core = round((0.50 + 0.20 + 25.0 + 6.0 + 1.5) * diversity_scale, 4)
        row["predicted_core_hours_base"] = base_core
        row["predicted_core_hours_high"] = high_core
        row["predicted_peak_memory_gib_base"] = round(max(4, 64 * diversity_scale), 2)
        row["predicted_peak_memory_gib_high"] = round(min(192.0, max(8, 96 * diversity_scale)), 2)
        row["predicted_wall_hours_base"] = round((0.10 + 0.05 + 0.175 + 0.05 + 0.075) * diversity_scale, 4)
        row["predicted_wall_hours_high"] = round((0.50 + 0.20 + 3.125 + 0.75 + 0.75) * diversity_scale, 4)
        row["predicted_scratch_gb_base"] = round((7.5 + 0.5 + 7.5 + 3.0 + 0.5) * diversity_scale, 4)
        row["predicted_scratch_gb_high"] = round((37.5 + 2.0 + 37.5 + 15.0 + 2.0) * diversity_scale, 4)
        row["predicted_inode_count_base"] = int(round((150 + 300 + 500 + 1500 + 1000) * diversity_scale))
        row["predicted_inode_count_high"] = int(round((500 + 1000 + 3000 + 10000 + 5000) * diversity_scale))
        row["predicted_moosefs_read_gb_base"] = round((4.0 + 4.0 + 2.5 + 1.0 + 0.5) * diversity_scale, 4)
        row["predicted_moosefs_read_gb_high"] = round((20.0 + 20.0 + 12.5 + 5.0 + 2.0) * diversity_scale, 4)
        row["predicted_moosefs_write_gb_base"] = round((0.02 + 0.02 + 0.6 + 0.4 + 0.08) * diversity_scale, 4)
        row["predicted_moosefs_write_gb_high"] = round((0.05 + 0.10 + 6.0 + 2.0 + 0.3) * diversity_scale, 4)
        row["predicted_metadata_operations_base"] = int(round((1000 + 3000 + 5000 + 15000 + 5000) * diversity_scale))
        row["predicted_metadata_operations_high"] = int(round((5000 + 15000 + 50000 + 100000 + 30000) * diversity_scale))
    else:
        base_core = round((0.10 + 0.05 + 0.135) * composition_scale, 4)
        high_core = round((0.50 + 0.20 + 0.532) * composition_scale, 4)
        row["predicted_core_hours_base"] = base_core
        row["predicted_core_hours_high"] = high_core
        row["predicted_peak_memory_gib_base"] = round(max(2, 4 * composition_scale), 2)
        row["predicted_peak_memory_gib_high"] = round(min(32.0, max(4, 16 * composition_scale)), 2)
        row["predicted_wall_hours_base"] = round((0.10 + 0.05 + 0.0675) * composition_scale, 4)
        row["predicted_wall_hours_high"] = round((0.50 + 0.20 + 0.266) * composition_scale, 4)
        row["predicted_scratch_gb_base"] = round((7.5 + 0.5 + 0.3) * composition_scale, 4)
        row["predicted_scratch_gb_high"] = round((37.5 + 2.0 + 1.0) * composition_scale, 4)
        row["predicted_inode_count_base"] = int(round((150 + 300 + 300) * composition_scale))
        row["predicted_inode_count_high"] = int(round((500 + 1000 + 1000) * composition_scale))
        row["predicted_moosefs_read_gb_base"] = round((4.0 + 4.0 + 1.0) * composition_scale, 4)
        row["predicted_moosefs_read_gb_high"] = round((20.0 + 20.0 + 5.0) * composition_scale, 4)
        row["predicted_moosefs_write_gb_base"] = round((0.02 + 0.02 + 0.03) * composition_scale, 4)
        row["predicted_moosefs_write_gb_high"] = round((0.05 + 0.10 + 0.2) * composition_scale, 4)
        row["predicted_metadata_operations_base"] = int(round((1000 + 3000 + 2000) * composition_scale))
        row["predicted_metadata_operations_high"] = int(round((5000 + 15000 + 10000) * composition_scale))


def select_pilot(rows: List[MutableMapping[str, Any]], *, limit: int = 6) -> None:
    eligible = [
        row
        for row in rows
        if row["assembly_composition_eligible"] == "yes" or row["assembly_diversity_eligible"] == "yes"
    ]
    selection_categories(eligible, "combined_genome_bp")
    selection_categories(eligible, "combined_contig_count")
    selection_categories(eligible, "annotation_gff_uncompressed_bytes")
    selection_categories(eligible, "predicted_core_hours_high")
    for row in eligible:
        row["pilot_selection_group"] = normalize_group(row.get("class"), row["lineage_group"])
        row["pilot_genome_stratum"] = row.get("combined_genome_bp_stratum", "unknown")
        row["pilot_contiguity_stratum"] = row.get("combined_contig_count_stratum", "unknown")
        row["pilot_annotation_stratum"] = row.get("annotation_gff_uncompressed_bytes_stratum", "unknown")
        row["pilot_resource_stratum"] = row.get("predicted_core_hours_high_stratum", "unknown")

    covered: Dict[str, set] = {
        "group": set(),
        "genome": set(),
        "contiguity": set(),
        "annotation": set(),
        "resource": set(),
    }
    selected: List[MutableMapping[str, Any]] = []
    remaining = eligible[:]
    while remaining and len(selected) < limit:
        def score(row: Mapping[str, Any]) -> Tuple[float, float, float]:
            novelty = 0
            novelty += 1 if row["pilot_selection_group"] not in covered["group"] else 0
            novelty += 1 if row["pilot_genome_stratum"] not in covered["genome"] else 0
            novelty += 1 if row["pilot_contiguity_stratum"] not in covered["contiguity"] else 0
            novelty += 1 if row["pilot_annotation_stratum"] not in covered["annotation"] else 0
            novelty += 1 if row["pilot_resource_stratum"] not in covered["resource"] else 0
            diversity_bonus = 0.25 if row["assembly_diversity_eligible"] == "yes" else 0.0
            size_penalty = -float(row["predicted_download_bytes_exact"] or 0)
            return novelty + diversity_bonus, size_penalty, -row["catalog_row_number"]

        chosen = max(remaining, key=score)
        selected.append(chosen)
        covered["group"].add(chosen["pilot_selection_group"])
        covered["genome"].add(chosen["pilot_genome_stratum"])
        covered["contiguity"].add(chosen["pilot_contiguity_stratum"])
        covered["annotation"].add(chosen["pilot_annotation_stratum"])
        covered["resource"].add(chosen["pilot_resource_stratum"])
        remaining.remove(chosen)

    for row in rows:
        row["pilot_selected"] = "yes" if row in selected else "no"
        if row["pilot_selected"] == "yes":
            row["acceptance_status"] = "selected_pilot"
        elif row["assembly_composition_eligible"] == "yes" or row["assembly_diversity_eligible"] == "yes":
            row["acceptance_status"] = "eligible_not_selected"


def build_row(seed: SeedRow, frozen: Mapping[str, Any], resolver: Resolver) -> MutableMapping[str, Any]:
    h1 = resolver.resolve_assembly(seed.h1_accession, mode="primary")
    h2, pair_status, pair_audit = resolve_pair(seed, h1, resolver)
    annotation_source, annotation_status = annotation_bundle(seed, h1, resolver)
    class_name = lineage_rank(h1.taxonomy_lineage, "class")
    order_name = lineage_rank(h1.taxonomy_lineage, "order")
    lineage_group = normalize_group(class_name, seed.manifest_lineage)
    evidence_retrieved_at = utc_now()

    h1_ok = h1.fasta.provider_md5 and h1.fasta.compressed_bytes and h1.fasta.uncompressed_bytes
    annotation_ok = (
        annotation_source
        and annotation_source.gff
        and annotation_source.gff.provider_md5
        and annotation_source.gff.compressed_bytes
        and annotation_source.gff.uncompressed_bytes
    )
    h2_ok = (
        h2
        and h2.fasta.provider_md5
        and h2.fasta.compressed_bytes
        and h2.fasta.uncompressed_bytes
        and pair_status == "ok"
    )
    composition_eligible = "yes" if annotation_ok and h1_ok else "no"
    diversity_eligible = "yes" if composition_eligible == "yes" and h2_ok else "no"

    reasons = []
    blocking = []
    if composition_eligible == "no":
        reasons.append("missing_exact_h1_annotation_or_file_size_checksum_evidence")
        blocking.append("B03")
    if SEED_MODALITY_TIER3A in seed.seed_modalities and diversity_eligible == "no":
        reasons.append(pair_status)
        blocking.append("B02")
    if not h1_ok:
        blocking.append("B01")
    acceptance_status = "rejected_unresolved"
    if composition_eligible == "yes" or diversity_eligible == "yes":
        acceptance_status = "eligible_pending_selection"

    row: MutableMapping[str, Any] = {
        "inventory_release": "phase1-freeze-1.0-2026-07-17",
        "record_role": seed.record_role,
        "candidate_id": seed.candidate_id,
        "catalog_row_number": seed.catalog_row_number,
        "scientific_name_source": seed.scientific_name,
        "ncbi_taxid": h1.taxid,
        "ncbi_current_name": h1.current_name,
        "class": class_name or "UNRESOLVED",
        "order": order_name or "UNRESOLVED",
        "source_catalog_url": frozen["source_url"],
        "source_catalog_revision": frozen["source_commit"],
        "source_catalog_sha256": frozen["sha256"],
        "source_retrieved_at_utc": frozen["retrieved_at_utc"],
        "catalog_evidence_status": "current_frozen_catalog_row",
        "h1_accession_version": h1.accession,
        "h1_release_version_or_date": h1.release_date or "UNRESOLVED",
        "h1_assembly_name": h1.assembly_name or "UNRESOLVED",
        "h1_haplotype_role": h1.diploid_role or h1.assembly_type or "UNRESOLVED",
        "h1_fasta_url": h1.fasta.url,
        "h1_provider_md5": h1.fasta.provider_md5 or "UNRESOLVED",
        "h1_fasta_sha256": "UNRESOLVED_NOT_DOWNLOADED",
        "h1_sequence_set_sha256": "UNRESOLVED_NOT_DOWNLOADED",
        "h1_length_bp": h1.total_sequence_length or "UNRESOLVED",
        "h1_contig_count": h1.contig_count or "UNRESOLVED",
        "h1_contig_n50_bp": h1.contig_n50 or "UNRESOLVED",
        "h2_accession_version": h2.accession if h2 else "UNRESOLVED",
        "h2_release_version_or_date": h2.release_date if h2 and h2.release_date else "UNRESOLVED",
        "h2_assembly_name_or_label": h2.assembly_name if h2 and h2.assembly_name else ",".join(seed.h2_accessions) or "UNRESOLVED",
        "h2_haplotype_role": h2.diploid_role if h2 and h2.diploid_role else "UNRESOLVED",
        "h2_fasta_url": h2.fasta.url if h2 else "UNRESOLVED",
        "h2_provider_md5": h2.fasta.provider_md5 if h2 and h2.fasta.provider_md5 else "UNRESOLVED",
        "h2_fasta_sha256": "UNRESOLVED_NOT_DOWNLOADED" if h2 else "UNRESOLVED",
        "h2_sequence_set_sha256": "UNRESOLVED_NOT_DOWNLOADED" if h2 else "UNRESOLVED",
        "h2_length_bp": h2.total_sequence_length if h2 and h2.total_sequence_length else "UNRESOLVED",
        "h2_contig_count": h2.contig_count if h2 and h2.contig_count else "UNRESOLVED",
        "h2_contig_n50_bp": h2.contig_n50 if h2 and h2.contig_n50 else "UNRESOLVED",
        "biosample_accession": h1.biosample_accession or "UNRESOLVED",
        "individual_or_isolate_id": h1.tolid or h1.isolate or h1.yggdrasil_individual or "UNRESOLVED",
        "h1_h2_relationship": "same_individual_linked_pair" if pair_status == "ok" else pair_status,
        "haplotype_contig_map_sha256": "UNRESOLVED_PRECOMPUTE",
        "haplotype_contig_relationship_audit": pair_audit,
        "pair_evidence_url": h1.report_url,
        "pair_evidence_retrieved_at_utc": evidence_retrieved_at,
        "annotation_accession_version": (
            annotation_source.annotation_name or annotation_source.accession
            if annotation_source
            else "UNRESOLVED"
        ),
        "annotation_reference_accession_version": h1.accession,
        "annotation_release_version_or_date": (
            annotation_source.annotation_release_date
            if annotation_source and annotation_source.annotation_release_date
            else "UNRESOLVED"
        ),
        "annotation_provider_and_pipeline": (
            annotation_source.annotation_provider_pipeline
            if annotation_source and annotation_source.annotation_provider_pipeline
            else "UNRESOLVED"
        ),
        "annotation_gff_url": (
            annotation_source.gff.url
            if annotation_source and annotation_source.gff
            else "UNRESOLVED"
        ),
        "annotation_provider_md5": (
            annotation_source.gff.provider_md5
            if annotation_source and annotation_source.gff and annotation_source.gff.provider_md5
            else "UNRESOLVED"
        ),
        "annotation_gff_sha256": (
            "UNRESOLVED_NOT_DOWNLOADED" if annotation_source and annotation_source.gff else "UNRESOLVED"
        ),
        "annotation_native_status": annotation_status if annotation_ok else annotation_status,
        "annotation_contig_map_sha256": "UNRESOLVED_PRECOMPUTE",
        "annotation_contig_audit": "metadata_only_pending_payload_audit" if annotation_ok else annotation_status,
        "cds_reconstruction_audit": "metadata_only_pending_payload_audit" if annotation_ok else annotation_status,
        "variant_resource_accession": "UNRESOLVED_PRECOMPUTE",
        "variant_reference_accession_version": h1.accession,
        "variant_url": "UNRESOLVED_PRECOMPUTE",
        "variant_sha256": "UNRESOLVED_PRECOMPUTE",
        "callability_resource_accession": "UNRESOLVED_PRECOMPUTE",
        "callability_reference_accession_version": h1.accession,
        "callability_url": "UNRESOLVED_PRECOMPUTE",
        "callability_sha256": "UNRESOLVED_PRECOMPUTE",
        "callable_bases": "UNRESOLVED_PRECOMPUTE",
        "callable_fraction": "UNRESOLVED_PRECOMPUTE",
        "queryable_gene_count": "UNRESOLVED_PRECOMPUTE",
        "queryable_gene_bases": "UNRESOLVED_PRECOMPUTE",
        "resource_retrieved_at_utc": evidence_retrieved_at,
        "license_or_reuse_terms": "NCBI_molecular_data_no_NCBI_restriction_submitter_rights_not_transferred",
        "license_evidence_url": NCBI_POLICY_URL,
        "evidence_summary": (
            f"seed_modalities={','.join(seed.seed_modalities)}; "
            f"h1_report={h1.report_url}; "
            f"h1_esummary={h1.esummary_url}; "
            f"h1_download_summary={h1.download_summary_url}; "
            f"taxonomy={NCBI_EFETCH_TAXONOMY.format(taxid=h1.taxid)}; "
            f"pair_status={pair_status}; "
            f"pair_audit={pair_audit}; "
            f"manifest_refseq={seed.refseq_annotation_accession or 'none'}; "
            f"annotation_status={annotation_status}; "
            f"annotation_report={annotation_source.report_url if annotation_source else 'none'}"
        ),
        "uncertainty_status": (
            "low"
            if diversity_eligible == "yes"
            else "medium"
            if composition_eligible == "yes"
            else "high"
        ),
        "assembly_composition_eligible": composition_eligible,
        "assembly_diversity_eligible": diversity_eligible,
        "population_genomic_eligible": "no",
        "demographic_eligible": "no",
        "acceptance_status": acceptance_status,
        "explicit_acceptance_or_rejection_reason": (
            "resolved_metadata_with_exact_sizes_and_checksums"
            if acceptance_status != "rejected_unresolved"
            else ";".join(reasons) or "unresolved_metadata"
        ),
        "blocking_requirement_ids": ";".join(sorted(set(blocking))) or "",
        "seed_modalities": ",".join(seed.seed_modalities),
        "lineage_group": seed.manifest_lineage,
        "manifest_order_code": seed.manifest_order_code,
        "manifest_family_name": seed.manifest_family,
        "manifest_annotation_status": seed.manifest_annotation_status,
        "manifest_refseq_accession": seed.refseq_annotation_accession or "UNRESOLVED",
        "linked_h2_accessions_ncbi": ",".join(h1.linked_assemblies),
        "h1_fasta_compressed_bytes": h1.fasta.compressed_bytes,
        "h1_fasta_uncompressed_bytes": h1.fasta.uncompressed_bytes,
        "h2_fasta_compressed_bytes": h2.fasta.compressed_bytes if h2 else None,
        "h2_fasta_uncompressed_bytes": h2.fasta.uncompressed_bytes if h2 else None,
        "annotation_gff_compressed_bytes": (
            annotation_source.gff.compressed_bytes if annotation_source and annotation_source.gff else None
        ),
        "annotation_gff_uncompressed_bytes": (
            annotation_source.gff.uncompressed_bytes if annotation_source and annotation_source.gff else None
        ),
        "same_individual_evidence": pair_audit,
        "same_individual_status": "yes" if pair_status == "ok" else "no",
        "annotation_file_status": annotation_source.gff.status if annotation_source and annotation_source.gff else "missing",
        "pilot_selection_group": "UNASSIGNED",
        "pilot_genome_stratum": "UNASSIGNED",
        "pilot_contiguity_stratum": "UNASSIGNED",
        "pilot_annotation_stratum": "UNASSIGNED",
        "pilot_resource_stratum": "UNASSIGNED",
        "pilot_selected": "no",
        "combined_genome_bp": (h1.total_sequence_length or 0) + (h2.total_sequence_length or 0 if h2 else 0),
        "combined_contig_count": max(h1.contig_count or 0, h2.contig_count or 0 if h2 else 0),
    }
    predicted_resources(row)
    return row


def write_tsv(path: Path, rows: Sequence[Mapping[str, Any]], columns: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(columns), delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    key: (
                        ""
                        if row.get(key) is None
                        else row[key]
                    )
                    for key in columns
                }
            )


def budget_rows(
    rows: Sequence[Mapping[str, Any]],
    *,
    root_config: Mapping[str, Any],
) -> List[Mapping[str, Any]]:
    selected = [row for row in rows if row["pilot_selected"] == "yes"]
    output: List[Mapping[str, Any]] = []
    for row in rows:
        output.append(
            {
                "row_type": "candidate",
                "candidate_id": row["candidate_id"],
                "scientific_name": row["scientific_name_source"],
                "pilot_selected": row["pilot_selected"],
                "pilot_selection_group": row["pilot_selection_group"],
                "assembly_composition_eligible": row["assembly_composition_eligible"],
                "assembly_diversity_eligible": row["assembly_diversity_eligible"],
                "download_bytes_exact": row["predicted_download_bytes_exact"],
                "persistent_storage_bytes_exact": row["predicted_persistent_storage_bytes_exact"],
                "core_hours_base": row["predicted_core_hours_base"],
                "core_hours_high": row["predicted_core_hours_high"],
                "peak_memory_gib_base": row["predicted_peak_memory_gib_base"],
                "peak_memory_gib_high": row["predicted_peak_memory_gib_high"],
                "wall_hours_base": row["predicted_wall_hours_base"],
                "wall_hours_high": row["predicted_wall_hours_high"],
                "scratch_gb_base": row["predicted_scratch_gb_base"],
                "scratch_gb_high": row["predicted_scratch_gb_high"],
                "inode_count_base": row["predicted_inode_count_base"],
                "inode_count_high": row["predicted_inode_count_high"],
                "moosefs_read_gb_base": row["predicted_moosefs_read_gb_base"],
                "moosefs_read_gb_high": row["predicted_moosefs_read_gb_high"],
                "moosefs_write_gb_base": row["predicted_moosefs_write_gb_base"],
                "moosefs_write_gb_high": row["predicted_moosefs_write_gb_high"],
                "metadata_operations_base": row["predicted_metadata_operations_base"],
                "metadata_operations_high": row["predicted_metadata_operations_high"],
                "formula_version": "vgp_pilot_budget_v1",
                "quota_state": "unavailable_fail_closed",
                "gate_readiness": "NO_GO_pending_quota_and_human_gate",
            }
        )

    def aggregate(name: str, subset: Sequence[Mapping[str, Any]]) -> Mapping[str, Any]:
        return {
            "row_type": name,
            "candidate_id": name,
            "scientific_name": name,
            "pilot_selected": "n/a",
            "pilot_selection_group": "n/a",
            "assembly_composition_eligible": "n/a",
            "assembly_diversity_eligible": "n/a",
            "download_bytes_exact": sum(int(row["predicted_download_bytes_exact"]) for row in subset),
            "persistent_storage_bytes_exact": sum(
                int(row["predicted_persistent_storage_bytes_exact"]) for row in subset
            ),
            "core_hours_base": round(sum(float(row["predicted_core_hours_base"]) for row in subset), 4),
            "core_hours_high": round(sum(float(row["predicted_core_hours_high"]) for row in subset), 4),
            "peak_memory_gib_base": round(
                max((float(row["predicted_peak_memory_gib_base"]) for row in subset), default=0.0), 4
            ),
            "peak_memory_gib_high": round(
                max((float(row["predicted_peak_memory_gib_high"]) for row in subset), default=0.0), 4
            ),
            "wall_hours_base": round(sum(float(row["predicted_wall_hours_base"]) for row in subset), 4),
            "wall_hours_high": round(sum(float(row["predicted_wall_hours_high"]) for row in subset), 4),
            "scratch_gb_base": round(
                max((float(row["predicted_scratch_gb_base"]) for row in subset), default=0.0), 4
            ),
            "scratch_gb_high": round(
                max((float(row["predicted_scratch_gb_high"]) for row in subset), default=0.0), 4
            ),
            "inode_count_base": sum(int(row["predicted_inode_count_base"]) for row in subset),
            "inode_count_high": sum(int(row["predicted_inode_count_high"]) for row in subset),
            "moosefs_read_gb_base": round(
                sum(float(row["predicted_moosefs_read_gb_base"]) for row in subset), 4
            ),
            "moosefs_read_gb_high": round(
                sum(float(row["predicted_moosefs_read_gb_high"]) for row in subset), 4
            ),
            "moosefs_write_gb_base": round(
                sum(float(row["predicted_moosefs_write_gb_base"]) for row in subset), 4
            ),
            "moosefs_write_gb_high": round(
                sum(float(row["predicted_moosefs_write_gb_high"]) for row in subset), 4
            ),
            "metadata_operations_base": sum(int(row["predicted_metadata_operations_base"]) for row in subset),
            "metadata_operations_high": sum(int(row["predicted_metadata_operations_high"]) for row in subset),
            "formula_version": "vgp_pilot_budget_v1",
            "quota_state": "unavailable_fail_closed",
            "gate_readiness": "NO_GO_pending_quota_and_human_gate",
        }

    output.append(aggregate("aggregate_selected", selected))
    output.append(aggregate("aggregate_eligible", [row for row in rows if row["acceptance_status"] != "rejected_unresolved"]))
    output.append(
        {
            "row_type": "filesystem_context",
            "candidate_id": "filesystem_context",
            "scientific_name": "filesystem_context",
            "pilot_selected": "n/a",
            "pilot_selection_group": "n/a",
            "assembly_composition_eligible": "n/a",
            "assembly_diversity_eligible": "n/a",
            "download_bytes_exact": "",
            "persistent_storage_bytes_exact": "",
            "core_hours_base": "",
            "core_hours_high": "",
            "peak_memory_gib_base": "",
            "peak_memory_gib_high": "",
            "wall_hours_base": "",
            "wall_hours_high": "",
            "scratch_gb_base": "",
            "scratch_gb_high": "",
            "inode_count_base": "",
            "inode_count_high": "",
            "moosefs_read_gb_base": "",
            "moosefs_read_gb_high": "",
            "moosefs_write_gb_base": "",
            "moosefs_write_gb_high": "",
            "metadata_operations_base": "",
            "metadata_operations_high": "",
            "formula_version": "vgp_pilot_budget_v1",
            "quota_state": "unavailable_fail_closed",
            "gate_readiness": (
                f"NO_GO_pending_quota_interface_under_{root_config['root']}; "
                "filesystem free space exists but user quota evidence is unavailable"
            ),
        }
    )
    return output


def rejection_rows(rows: Sequence[Mapping[str, Any]]) -> List[Mapping[str, Any]]:
    rejects = []
    for row in rows:
        if row["acceptance_status"] == "rejected_unresolved":
            rejects.append(
                {
                    "candidate_id": row["candidate_id"],
                    "scientific_name": row["scientific_name_source"],
                    "catalog_row_number": row["catalog_row_number"],
                    "seed_modalities": row["seed_modalities"],
                    "h1_accession_version": row["h1_accession_version"],
                    "h2_accession_version": row["h2_accession_version"],
                    "ncbi_taxid": row["ncbi_taxid"],
                    "class": row["class"],
                    "order": row["order"],
                    "acceptance_status": row["acceptance_status"],
                    "explicit_rejection_reason": row["explicit_acceptance_or_rejection_reason"],
                    "blocking_requirement_ids": row["blocking_requirement_ids"],
                    "pair_evidence_url": row["pair_evidence_url"],
                    "resource_retrieved_at_utc": row["resource_retrieved_at_utc"],
                    "annotation_native_status": row["annotation_native_status"],
                    "annotation_file_status": row["annotation_file_status"],
                    "same_individual_status": row["same_individual_status"],
                    "same_individual_evidence": row["same_individual_evidence"],
                }
            )
    return rejects


def validate_rows(rows: Sequence[Mapping[str, Any]], frozen: Mapping[str, Any]) -> None:
    if frozen["line_count"] != 717:
        raise RuntimeError(f"expected 717 lines in pinned TSV, observed {frozen['line_count']}")
    if sum(1 for row in rows if row["pilot_selected"] == "yes") > 6:
        raise RuntimeError("selected pilot exceeds six-species ceiling")
    seen = set()
    for row in rows:
        candidate_id = row["candidate_id"]
        if candidate_id in seen:
            raise RuntimeError(f"duplicate candidate_id {candidate_id}")
        seen.add(candidate_id)
        if row["pilot_selected"] == "yes":
            required = [
                "h1_accession_version",
                "annotation_accession_version",
                "annotation_gff_url",
                "h1_fasta_url",
                "h1_provider_md5",
            ]
            if row["assembly_diversity_eligible"] == "yes":
                required.extend(["h2_accession_version", "h2_fasta_url", "h2_provider_md5"])
            for key in required:
                if not row.get(key) or str(row[key]).startswith("UNRESOLVED"):
                    raise RuntimeError(f"selected row {candidate_id} missing required field {key}")
            if not row.get("predicted_download_bytes_exact"):
                raise RuntimeError(f"selected row {candidate_id} lacks exact download bytes")


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root-config", type=Path, default=DEFAULT_ROOT_CONFIG)
    parser.add_argument("--provenance-out", type=Path, default=DEFAULT_PROVENANCE)
    parser.add_argument("--manifest-out", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--rejections-out", type=Path, default=DEFAULT_REJECTIONS)
    parser.add_argument("--budget-out", type=Path, default=DEFAULT_BUDGET)
    parser.add_argument("--source-commit", default=None)
    parser.add_argument("--sleep-seconds", type=float, default=0.0)
    args = parser.parse_args(argv)

    root_config = load_root_config(args.root_config)
    frozen = freeze_raw_manifest(root_config=root_config, source_commit=args.source_commit)
    rows = parse_vgp_rows(frozen["text"])
    counts = compute_counts(rows)
    seeds = build_seed_rows(rows)

    resolver = Resolver()
    resolved_rows: List[MutableMapping[str, Any]] = []
    for seed in seeds:
        resolved_rows.append(build_row(seed, frozen, resolver))
        if args.sleep_seconds:
            time.sleep(args.sleep_seconds)
    select_pilot(resolved_rows, limit=6)
    validate_rows(resolved_rows, frozen)

    manifest_columns = CANDIDATE_SCHEMA_COLUMNS + EXTRA_MANIFEST_COLUMNS
    write_tsv(args.manifest_out, resolved_rows, manifest_columns)
    write_tsv(
        args.rejections_out,
        rejection_rows(resolved_rows),
        [
            "candidate_id",
            "scientific_name",
            "catalog_row_number",
            "seed_modalities",
            "h1_accession_version",
            "h2_accession_version",
            "ncbi_taxid",
            "class",
            "order",
            "acceptance_status",
            "explicit_rejection_reason",
            "blocking_requirement_ids",
            "pair_evidence_url",
            "resource_retrieved_at_utc",
            "annotation_native_status",
            "annotation_file_status",
            "same_individual_status",
            "same_individual_evidence",
        ],
    )
    write_tsv(
        args.budget_out,
        budget_rows(resolved_rows, root_config=root_config),
        [
            "row_type",
            "candidate_id",
            "scientific_name",
            "pilot_selected",
            "pilot_selection_group",
            "assembly_composition_eligible",
            "assembly_diversity_eligible",
            "download_bytes_exact",
            "persistent_storage_bytes_exact",
            "core_hours_base",
            "core_hours_high",
            "peak_memory_gib_base",
            "peak_memory_gib_high",
            "wall_hours_base",
            "wall_hours_high",
            "scratch_gb_base",
            "scratch_gb_high",
            "inode_count_base",
            "inode_count_high",
            "moosefs_read_gb_base",
            "moosefs_read_gb_high",
            "moosefs_write_gb_base",
            "moosefs_write_gb_high",
            "metadata_operations_base",
            "metadata_operations_high",
            "formula_version",
            "quota_state",
            "gate_readiness",
        ],
    )

    provenance = {
        "generated_at_utc": utc_now(),
        "task_id": "freeze-vgp-manifest",
        "source_catalog": {
            key: value
            for key, value in frozen.items()
            if key not in {"text"}
        },
        "counts": counts,
        "candidate_summary": {
            "seed_count": len(seeds),
            "resolved_row_count": len(resolved_rows),
            "selected_count": sum(1 for row in resolved_rows if row["pilot_selected"] == "yes"),
            "rejected_count": sum(
                1 for row in resolved_rows if row["acceptance_status"] == "rejected_unresolved"
            ),
        },
        "outputs": {
            "manifest_tsv": str(args.manifest_out),
            "rejections_tsv": str(args.rejections_out),
            "budget_tsv": str(args.budget_out),
        },
        "quota_gate": {
            "status": "NO_GO_pending_quota_interface",
            "source": str(args.root_config),
        },
    }
    args.provenance_out.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
