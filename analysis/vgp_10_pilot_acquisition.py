#!/usr/bin/env python3
"""Acquire the closed-world ten-pair VGP pilot into the Freeze 1 CAS.

The program has deliberately narrow authority: the ten frozen primary pairs,
the six frozen alternates (metadata plus untransferred, pre-ranked payload
plans), and one separately selected direct-control manifest.  It never queries
by species to choose an assembly and never enumerates or downloads SRA runs.
Every payload is staged as a partial, checked against size and the provider MD5
when NCBI publishes one, hashed with SHA-256, and atomically promoted.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import shutil
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Iterable, Mapping, MutableMapping, Sequence


ROOT = Path(__file__).resolve().parents[1]
PRIMARY = ROOT / "analysis/vgp_10_pair_manifest.tsv"
ALTERNATES = ROOT / "analysis/vgp_10_pair_alternates.tsv"
PAIR_OUTPUT = ROOT / "analysis/vgp_10_pilot_acquisition_manifest.tsv"
OBJECT_OUTPUT = ROOT / "analysis/vgp_10_pilot_object_inventory.tsv"
DIRECT_OUTPUT = ROOT / "analysis/vgp_direct_control_acquisition_manifest.tsv"
SUMMARY_OUTPUT = ROOT / "analysis/vgp_10_pilot_acquisition_summary.json"
HANDOFF_OUTPUT = ROOT / "analysis/vgp_10_pilot_acquisition_handoff.md"
DEFAULT_DATA_ROOT = Path("/moosefs/erikg/lewontin-paradox-data/vgp/phase1-freeze-1.0")
MANIFEST_VERSION = "vgp-10-pilot-acquisition-v1.0.0"
NCBI_POLICY = (
    "NCBI molecular data have no NCBI-imposed use restriction; submitter rights "
    "are not transferred; https://www.ncbi.nlm.nih.gov/home/about/policies/"
)
VGP_RELEASE_PROJECT = "PRJNA489243"
VERSIONED_ASSEMBLY = re.compile(r"^GC[AF]_\d+\.\d+$")
HEX32 = re.compile(r"^[0-9a-f]{32}$")
HEX64 = re.compile(r"^[0-9a-f]{64}$")
SLURM_ENV = ("SLURM_JOB_ID", "SLURM_ARRAY_JOB_ID", "SLURM_ARRAY_TASK_ID")
BULK_RAW_ROLES = {
    "raw_reads",
    "raw_read_archive",
    "bulk_raw_reads",
    "sra_run",
    "hifi_bam",
    "hic_fastq",
    "trio_parent_reads",
}

PAIR_FIELDS = [
    "manifest_version", "selection_id", "roster_type", "rank", "activation_status",
    "replacement_rule", "amendment_id", "failed_primary_retained", "catalog_row",
    "catalog_commit", "catalog_sha256", "release_membership", "species", "taxid",
    "biosample", "individual_or_isolate", "h1_accession_version", "h1_label", "h1_role",
    "h1_release_date", "h1_length_bp", "h1_contigs", "h1_contig_n50_bp", "h1_scaffold_n50_bp",
    "h1_frozen_report_sha256", "h1_assembly_method",
    "h2_accession_version", "h2_label", "h2_role", "h2_release_date", "h2_length_bp",
    "h2_contigs", "h2_contig_n50_bp", "h2_scaffold_n50_bp", "h2_frozen_report_sha256",
    "h2_assembly_method", "reciprocal_linked_assembly_evidence",
    "assembly_generation", "assembly_technologies", "hifi_evidence", "long_range_phasing_evidence",
    "repeat_report_evidence", "repeat_report_status",
    "qv_evidence", "qv_status", "completeness_evidence", "completeness_status",
    "duplication_collapse_evidence", "duplication_collapse_status",
    "native_annotation_accession_version", "annotation_branch_status",
    "annotation_binding_evidence", "core_identity_status", "core_acquisition_status",
    "raw_read_or_kmer_selection", "selective_asset_justification", "license_or_access",
]

OBJECT_FIELDS = [
    "manifest_version", "object_id", "selection_id", "roster_type", "activation_status",
    "side", "object_role", "analysis_branch", "accession_version",
    "annotation_accession_version", "source_provider", "source_url", "release_membership",
    "expected_bytes", "observed_bytes", "upstream_checksum_algorithm", "upstream_checksum",
    "checksum_policy", "expected_local_sha256", "local_sha256", "retrieval_utc",
    "source_last_modified", "license_or_access", "status", "status_detail",
    "resume_from_bytes", "revalidation_count", "transferred_bytes", "local_path",
    "quarantine_path", "pair_identity_status", "annotation_assembly_match",
    "selective_justification",
]

DIRECT_FIELDS = [
    "manifest_version", "control_id", "selection_status", "dataset_title", "species",
    "taxid", "repository", "dataset_accession_version", "object_id", "object_role",
    "source_url", "expected_bytes", "upstream_checksum_algorithm", "upstream_checksum",
    "observed_bytes", "transferred_bytes", "local_sha256", "retrieval_utc", "status", "local_path", "relationships",
    "complete_pedigree_or_gamete_evidence", "transmitted_haplotype_evidence",
    "directional_transmission_suitability", "consent", "access", "license_or_terms",
    "selective_justification",
]


class AcquisitionRefusal(RuntimeError):
    """A hard provenance or content gate failed."""


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def hash_file(path: Path) -> tuple[str, str, int]:
    sha = hashlib.sha256()
    md5 = hashlib.md5()  # noqa: S324 - provider transport checksum only
    size = 0
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            sha.update(block)
            md5.update(block)
            size += len(block)
    return sha.hexdigest(), md5.hexdigest(), size


def canonical_json(payload: Any) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def atomic_write(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.partial-{os.getpid()}")
    with temporary.open("wb") as handle:
        handle.write(payload)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def write_tsv(path: Path, fields: Sequence[str], rows: Sequence[Mapping[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.partial-{os.getpid()}")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def partial_path(data_root: Path, object_id: str) -> Path:
    safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(object_id))
    return data_root / "staging/partials/vgp-10-pilot" / f"{safe}.partial"


def cas_path(cas_root: Path, digest: str) -> Path:
    if not HEX64.fullmatch(str(digest)):
        raise AcquisitionRefusal(f"invalid SHA-256 object identity: {digest!r}")
    return cas_root / digest[:2] / digest[2:4] / digest


def _fsync_dir(path: Path) -> None:
    descriptor = os.open(path, os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def quarantine(path: Path, data_root: Path, reason: str) -> Path:
    target_dir = data_root / "quarantine/vgp-10-pilot" / utc_now().replace(":", "")
    target_dir.mkdir(parents=True, exist_ok=True)
    target = target_dir / f"{path.name}.{re.sub('[^A-Za-z0-9_.-]+', '_', reason)[:80]}"
    os.replace(path, target)
    _fsync_dir(target_dir)
    return target


Downloader = Callable[[dict[str, object], Path, int], int]


def http_downloader(row: dict[str, object], path: Path, offset: int) -> int:
    request = urllib.request.Request(str(row["source_url"]), headers={"User-Agent": "lewontin-vgp-10-pilot/1.0"})
    if offset:
        request.add_header("Range", f"bytes={offset}-")
    with urllib.request.urlopen(request, timeout=180) as response:  # noqa: S310 - exact manifest URL
        status = int(getattr(response, "status", response.getcode()))
        if offset and status != 206:
            path.unlink(missing_ok=True)
            offset = 0
        mode = "ab" if offset else "wb"
        transferred = 0
        with path.open(mode) as handle:
            for block in iter(lambda: response.read(4 * 1024 * 1024), b""):
                handle.write(block)
                transferred += len(block)
            handle.flush()
            os.fsync(handle.fileno())
    return transferred


def acquire_one(
    object_row: Mapping[str, object], data_root: Path, *, downloader: Downloader = http_downloader
) -> dict[str, object]:
    row: dict[str, object] = dict(object_row)
    expected = int(row.get("expected_bytes") or 0)
    expected_sha = str(row.get("expected_local_sha256") or "")
    expected_md5 = str(row.get("upstream_checksum") or "") if row.get("upstream_checksum_algorithm") == "md5" else ""
    cas_root = data_root / "objects/sha256"
    row.update(
        observed_bytes=0, local_sha256="", retrieval_utc=utc_now(), status_detail="",
        resume_from_bytes=0, revalidation_count=0, transferred_bytes=0,
        local_path="", quarantine_path="",
    )

    if expected_sha:
        target = cas_path(cas_root, expected_sha)
        if target.is_file():
            observed_sha, observed_md5, observed_size = hash_file(target)
            row["revalidation_count"] = 1
            if observed_sha != expected_sha or (expected and observed_size != expected) or (expected_md5 and observed_md5 != expected_md5):
                bad = quarantine(target, data_root, "mirror-revalidation-mismatch")
                row.update(status="quarantined", observed_bytes=observed_size, local_sha256=observed_sha, quarantine_path=str(bad), status_detail="mirror object failed local SHA-256/size/provider-digest revalidation")
                return row
            row.update(status="reused", observed_bytes=observed_size, local_sha256=observed_sha, local_path=str(target), status_detail="official Freeze 1 CAS hit; local SHA-256 revalidated")
            return row

    part = partial_path(data_root, str(row["object_id"]))
    part.parent.mkdir(parents=True, exist_ok=True)
    offset = part.stat().st_size if part.exists() else 0
    row["resume_from_bytes"] = offset
    if expected and offset > expected:
        bad = quarantine(part, data_root, "partial-oversize")
        row.update(status="quarantined", observed_bytes=offset, quarantine_path=str(bad), status_detail="partial exceeded expected bytes before resume")
        return row
    try:
        row["transferred_bytes"] = downloader(row, part, offset)
    except Exception as error:  # retain resumable partial for transport-only failure
        row.update(status="missing", observed_bytes=part.stat().st_size if part.exists() else 0, status_detail=f"retrieval incomplete; resumable partial retained: {error}")
        return row

    observed_sha, observed_md5, observed_size = hash_file(part)
    row.update(observed_bytes=observed_size, local_sha256=observed_sha)
    mismatch = []
    if expected and observed_size != expected:
        mismatch.append(f"size {observed_size} != {expected}")
    if expected_md5 and observed_md5 != expected_md5:
        mismatch.append(f"provider MD5 {observed_md5} != {expected_md5}")
    if expected_sha and observed_sha != expected_sha:
        mismatch.append(f"SHA-256 {observed_sha} != {expected_sha}")
    if mismatch:
        bad = quarantine(part, data_root, "content-mismatch")
        row.update(status="quarantined", quarantine_path=str(bad), status_detail="; ".join(mismatch))
        return row

    target = cas_path(cas_root, observed_sha)
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists():
        target_sha, target_md5, target_size = hash_file(target)
        row["revalidation_count"] = 1
        if target_sha != observed_sha or target_size != observed_size or (expected_md5 and target_md5 != expected_md5):
            bad = quarantine(target, data_root, "cas-collision-mismatch")
            part_bad = quarantine(part, data_root, "cas-collision-source")
            row.update(status="quarantined", quarantine_path=f"{bad};{part_bad}", status_detail="existing CAS object did not match staged object")
            return row
        part.unlink()
        row.update(status="reused", local_path=str(target), status_detail="concurrent/offline mirror object reused after SHA-256 revalidation")
        return row
    os.replace(part, target)
    os.chmod(target, 0o440)
    _fsync_dir(target.parent)
    promoted_sha, _, promoted_size = hash_file(target)
    row["revalidation_count"] = 1
    if promoted_sha != observed_sha or promoted_size != observed_size:
        bad = quarantine(target, data_root, "post-promotion-revalidation")
        row.update(status="quarantined", local_path="", quarantine_path=str(bad), status_detail="post-promotion SHA-256 revalidation failed")
        return row
    row.update(status="verified", local_path=str(target), status_detail="staged, verified, atomically promoted read-only, and SHA-256 reverified")
    return row


def fetch_bytes(url: str, *, attempts: int = 4) -> tuple[bytes, dict[str, str]]:
    error: Exception | None = None
    for attempt in range(attempts):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": "lewontin-vgp-10-pilot/1.0"})
            with urllib.request.urlopen(request, timeout=120) as response:  # noqa: S310 - version-bound HTTPS
                return response.read(), {key.lower(): value for key, value in response.headers.items()}
        except (OSError, urllib.error.URLError) as caught:
            error = caught
            if attempt + 1 < attempts:
                time.sleep(min(2 ** attempt, 4))
    raise AcquisitionRefusal(f"failed to retrieve immutable metadata URL {url}: {error}")


def head_bytes(url: str) -> tuple[int, str]:
    request = urllib.request.Request(url, method="HEAD", headers={"User-Agent": "lewontin-vgp-10-pilot/1.0"})
    with urllib.request.urlopen(request, timeout=120) as response:  # noqa: S310
        length = response.headers.get("Content-Length")
        if not length or not length.isdigit():
            raise AcquisitionRefusal(f"upstream did not publish expected bytes for {url}")
        return int(length), response.headers.get("Last-Modified", "")


def report_url(accession: str) -> str:
    if not VERSIONED_ASSEMBLY.fullmatch(accession):
        raise AcquisitionRefusal(f"unversioned assembly accession: {accession}")
    return f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{accession}/dataset_report"


def report_for(accession: str) -> tuple[dict[str, Any], bytes, dict[str, str]]:
    body, headers = fetch_bytes(report_url(accession))
    payload = json.loads(body)
    reports = payload.get("reports", [])
    exact = [item for item in reports if item.get("accession") == accession]
    if len(exact) != 1 or exact[0].get("current_accession") != accession:
        raise AcquisitionRefusal(f"accession drift or ambiguous dataset report for {accession}")
    return exact[0], body, headers


def _isolate(report: Mapping[str, Any]) -> str:
    info = report.get("assembly_info", {})
    organism = report.get("organism", {})
    biosample = info.get("biosample", {})
    return str(
        # The frozen roster binds the BioSample individual/isolate.  NCBI may
        # additionally expose a submitter's assembly-level isolate alias (for
        # example SK-2024b for BioSample isolate fSpiSpi1); both haplotypes
        # must still agree on the exact BioSample value.
        biosample.get("isolate")
        or organism.get("infraspecific_names", {}).get("isolate")
        or organism.get("infraspecific_names", {}).get("strain")
        or biosample.get("strain")
        or ""
    )


def _linked(report: Mapping[str, Any]) -> dict[str, str]:
    linked = report.get("assembly_info", {}).get("linked_assemblies", [])
    return {str(item.get("linked_assembly")): str(item.get("assembly_type")) for item in linked}


def synthetic_report(roster: Mapping[str, str], side: str) -> dict[str, Any]:
    other = "h2" if side == "h1" else "h1"
    return {
        "accession": roster[f"{side}_accession_version"],
        "current_accession": roster[f"{side}_accession_version"],
        "organism": {"tax_id": int(roster["resolved_taxid"]), "organism_name": roster["species"], "infraspecific_names": {"isolate": roster["individual_or_isolate"]}},
        "assembly_info": {
            "biosample": {"accession": roster["biosample"], "isolate": roster["individual_or_isolate"]},
            "diploid_role": roster[f"{side}_role"],
            "linked_assemblies": [{"linked_assembly": roster[f"{other}_accession_version"], "assembly_type": roster[f"{other}_role"]}],
            "bioproject_lineage": [{"bioprojects": [{"accession": VGP_RELEASE_PROJECT}]}],
        },
    }


def validate_pair_identity(roster: Mapping[str, str], h1: Mapping[str, Any], h2: Mapping[str, Any]) -> dict[str, str]:
    for side, report in (("h1", h1), ("h2", h2)):
        expected = roster[f"{side}_accession_version"]
        if report.get("accession") != expected or report.get("current_accession") != expected:
            raise AcquisitionRefusal(f"{side} accession drift: expected {expected}")
        organism = report.get("organism", {})
        observed_name = str(organism.get("organism_name", ""))
        frozen_name = str(roster["species"])
        name_matches = observed_name == frozen_name or frozen_name.startswith(observed_name + " ")
        if str(organism.get("tax_id")) != str(roster["resolved_taxid"]) or not name_matches:
            raise AcquisitionRefusal(f"{side} species/TaxId identity mismatch")
        biosample = report.get("assembly_info", {}).get("biosample", {})
        if biosample.get("accession") != roster["biosample"]:
            raise AcquisitionRefusal(f"{side} BioSample identity mismatch")
        if _isolate(report) != roster["individual_or_isolate"]:
            raise AcquisitionRefusal(f"{side} individual/isolate identity mismatch")
    if roster["h2_accession_version"] not in _linked(h1) or roster["h1_accession_version"] not in _linked(h2):
        raise AcquisitionRefusal("H1/H2 reciprocal linked-assembly evidence is absent or ambiguous")
    return {
        "core_identity_status": "pass",
        "reciprocal_linked_assembly_evidence": (
            f"{roster['h1_accession_version']} links {roster['h2_accession_version']} as {_linked(h1)[roster['h2_accession_version']]}; "
            f"{roster['h2_accession_version']} links {roster['h1_accession_version']} as {_linked(h2)[roster['h1_accession_version']]}; "
            f"both exact reports bind BioSample {roster['biosample']} and isolate {roster['individual_or_isolate']}"
        ),
    }


def _project_accessions(report: Mapping[str, Any]) -> set[str]:
    result: set[str] = set()
    for lineage in report.get("assembly_info", {}).get("bioproject_lineage", []):
        for project in lineage.get("bioprojects", []):
            result.add(str(project.get("accession", "")))
    return result


def validate_annotation_binding(h1_accession: str, annotation_assembly: str, annotation_name: str) -> dict[str, str]:
    if annotation_assembly == h1_accession:
        return {"annotation_branch_status": "available_exact_native_pending_dictionary_audit", "annotation_binding_evidence": f"annotation {annotation_name} is deposited on exact H1 {h1_accession}"}
    h1_digits = h1_accession.replace("GCA_", "")
    if annotation_assembly.startswith("GCF_") and annotation_assembly.replace("GCF_", "") == h1_digits:
        return {"annotation_branch_status": "available_paired_refseq_pending_dictionary_audit", "annotation_binding_evidence": f"annotation {annotation_name} is on paired RefSeq {annotation_assembly}; exact H1 dictionary audit remains branch-local"}
    return {"annotation_branch_status": "failed_mismatch", "annotation_binding_evidence": f"annotation assembly {annotation_assembly} is not exact or paired with H1 {h1_accession}; core unaffected"}


def ftp_prefix(accession: str, assembly_name: str) -> str:
    digits = accession.split("_")[1].split(".")[0]
    chunks = "/".join(digits[index:index + 3] for index in range(0, 9, 3))
    safe_name = re.sub(r"[^A-Za-z0-9_.-]+", "_", assembly_name)
    directory = f"{accession}_{safe_name}"
    return f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{accession[:3]}/{chunks}/{directory}"


def parse_md5_catalog(payload: bytes) -> dict[str, str]:
    result: dict[str, str] = {}
    for line in payload.decode("utf-8").splitlines():
        match = re.fullmatch(r"([0-9a-f]{32})  \./(.+)", line)
        if match:
            result[match.group(2)] = match.group(1)
    return result


def prior_sha_by_url(path: Path = OBJECT_OUTPUT) -> dict[str, str]:
    if not path.is_file():
        return {}
    return {
        row["source_url"]: row["local_sha256"]
        for row in read_tsv(path)
        if row.get("source_url") and HEX64.fullmatch(row.get("local_sha256", ""))
    }


def base_object(**values: Any) -> dict[str, Any]:
    row = {field: "" for field in OBJECT_FIELDS}
    row.update(
        manifest_version=MANIFEST_VERSION, source_provider="NCBI", release_membership="VGP Phase 1 Freeze 1.0",
        checksum_policy="verify expected bytes and published provider MD5; compute SHA-256 before promotion and reverify after promotion/reuse",
        license_or_access=NCBI_POLICY, expected_bytes=0, observed_bytes=0, status="planned",
        resume_from_bytes=0, revalidation_count=0, transferred_bytes=0,
    )
    row.update(values)
    return row


def _inline_downloader(payload: bytes) -> Downloader:
    def download(_row: dict[str, object], path: Path, offset: int) -> int:
        if offset > len(payload):
            raise AcquisitionRefusal("inline partial exceeds immutable metadata bytes")
        mode = "ab" if offset else "wb"
        with path.open(mode) as handle:
            handle.write(payload[offset:])
            handle.flush()
            os.fsync(handle.fileno())
        return len(payload) - offset
    return download


def _annotation_accession(row: Mapping[str, str], h1_report: Mapping[str, Any]) -> tuple[str, str]:
    frozen = row["native_annotation_accession"]
    if frozen == "none":
        return "none", ""
    annotation_name = str(h1_report.get("annotation_info", {}).get("name") or h1_report.get("assembly_info", {}).get("paired_assembly", {}).get("annotation_name") or "")
    if frozen.endswith("current_at_acquisition"):
        if not annotation_name:
            return "none", ""
        return annotation_name, annotation_name.split("-", 1)[0]
    annotation_assembly = frozen.split("-", 1)[0]
    return frozen, annotation_assembly


def resolve_and_acquire(data_root: Path) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    if any(name in os.environ for name in SLURM_ENV):
        raise AcquisitionRefusal("Slurm environment detected; this task is local acquisition only")
    if os.environ.get("VGPPILOT_PINNED_GUIX") != "1":
        raise AcquisitionRefusal("run through analysis/run_vgp_10_pilot_acquisition_guix.sh (pinned GNU Guix)")
    primary_rows = read_tsv(PRIMARY)
    alternate_rows = read_tsv(ALTERNATES)
    if [row["selection_id"] for row in primary_rows] != [f"P{i:02d}" for i in range(1, 11)]:
        raise AcquisitionRefusal("primary roster is not exactly frozen P01..P10")
    if len(alternate_rows) < 5:
        raise AcquisitionRefusal("fewer than five frozen alternates")
    prior_sha = prior_sha_by_url()
    cas_by_provider_md5_size: dict[tuple[str, int], str] = {}
    if not prior_sha:
        # A prior fail-closed run may have promoted valid objects before the
        # final inventory was published.  Revalidate those content-addressed
        # files once and recover only an exact provider-MD5+size mapping.
        for candidate in (data_root / "objects/sha256").glob("[0-9a-f][0-9a-f]/[0-9a-f][0-9a-f]/[0-9a-f]*"):
            if not candidate.is_file() or not HEX64.fullmatch(candidate.name):
                continue
            observed_sha, observed_md5, observed_size = hash_file(candidate)
            if observed_sha != candidate.name:
                quarantine(candidate, data_root, "cas-bootstrap-sha-mismatch")
                continue
            cas_by_provider_md5_size[(observed_md5, observed_size)] = observed_sha
    pair_rows: list[dict[str, Any]] = []
    object_rows: list[dict[str, Any]] = []
    seen_object_ids: set[str] = set()

    def add_object(
        row: dict[str, Any], *, payload: bytes | None = None, acquire: bool = True,
        nonacquire_status: str = "superseded",
    ) -> None:
        if row["object_id"] in seen_object_ids:
            raise AcquisitionRefusal(f"duplicate object id {row['object_id']}")
        seen_object_ids.add(row["object_id"])
        row["expected_local_sha256"] = prior_sha.get(row["source_url"], row.get("expected_local_sha256", ""))
        if not row["expected_local_sha256"] and row.get("upstream_checksum_algorithm") == "md5":
            row["expected_local_sha256"] = cas_by_provider_md5_size.get(
                (str(row.get("upstream_checksum", "")), int(row.get("expected_bytes") or 0)), ""
            )
        if acquire:
            result = acquire_one(row, data_root, downloader=_inline_downloader(payload) if payload is not None else http_downloader)
        else:
            result = dict(row)
            detail = (
                "pre-ranked alternate payload not activated; primary retained; no transfer"
                if nonacquire_status == "superseded"
                else "exact provider directory publishes no standalone repeat report; no substitute inferred"
            )
            result.update(status=nonacquire_status, status_detail=detail, observed_bytes=0, local_sha256="", retrieval_utc="", resume_from_bytes=0, revalidation_count=0, transferred_bytes=0, local_path="", quarantine_path="")
        object_rows.append(result)

    for roster_type, rows in (("primary", primary_rows), ("alternate", alternate_rows)):
        for frozen in rows:
            selection_id = frozen["selection_id"]
            active = roster_type == "primary"
            reports: dict[str, dict[str, Any]] = {}
            report_bodies: dict[str, bytes] = {}
            report_headers: dict[str, dict[str, str]] = {}
            for side in ("h1", "h2"):
                accession = frozen[f"{side}_accession_version"]
                report, body, headers = report_for(accession)
                reports[side], report_bodies[side], report_headers[side] = report, body, headers
            identity = validate_pair_identity(frozen, reports["h1"], reports["h2"])
            if VGP_RELEASE_PROJECT not in _project_accessions(reports["h1"]) or VGP_RELEASE_PROJECT not in _project_accessions(reports["h2"]):
                raise AcquisitionRefusal(f"{selection_id} exact reports lack VGP release lineage {VGP_RELEASE_PROJECT}")

            annotation_name, annotation_assembly = _annotation_accession(frozen, reports["h1"])
            annotation_binding = (
                {"annotation_branch_status": "absent_optional", "annotation_binding_evidence": "no native annotation deposited; core unaffected"}
                if annotation_name == "none"
                else validate_annotation_binding(frozen["h1_accession_version"], annotation_assembly, annotation_name)
            )
            annotation_report: dict[str, Any] | None = None
            annotation_body: bytes | None = None
            annotation_headers: dict[str, str] = {}
            if annotation_assembly and annotation_assembly != frozen["h1_accession_version"]:
                annotation_report, annotation_body, annotation_headers = report_for(annotation_assembly)
                if annotation_report.get("paired_accession") != frozen["h1_accession_version"]:
                    annotation_binding = {"annotation_branch_status": "failed_mismatch", "annotation_binding_evidence": f"{annotation_assembly} report does not pair to exact H1; core unaffected"}
                actual_name = str(annotation_report.get("annotation_info", {}).get("name") or "")
                if actual_name and annotation_name != actual_name:
                    annotation_binding = {"annotation_branch_status": "failed_mismatch", "annotation_binding_evidence": f"annotation version drift: frozen/resolved {annotation_name}, report {actual_name}; core unaffected"}

            for side in ("h1", "h2"):
                accession = frozen[f"{side}_accession_version"]
                body = report_bodies[side]
                add_object(base_object(
                    object_id=f"{selection_id}:{side}:dataset_report", selection_id=selection_id,
                    roster_type=roster_type, activation_status="active_primary" if active else "standby_not_triggered",
                    side=side, object_role="dataset_report", analysis_branch="core_metadata",
                    accession_version=accession, source_url=report_url(accession), expected_bytes=len(body),
                    expected_local_sha256=sha256_bytes(body), upstream_checksum_algorithm="not_published",
                    upstream_checksum="", source_last_modified=report_headers[side].get("last-modified", ""),
                    pair_identity_status="pass", annotation_assembly_match="not_applicable",
                    selective_justification="exact versioned report required for pair identity and release membership",
                ), payload=body, acquire=True)

                report = reports[side]
                prefix = ftp_prefix(accession, str(report["assembly_info"]["assembly_name"]))
                catalog_url = prefix + "/md5checksums.txt"
                catalog_body, catalog_headers = fetch_bytes(catalog_url)
                catalog = parse_md5_catalog(catalog_body)
                add_object(base_object(
                    object_id=f"{selection_id}:{side}:checksum_catalog", selection_id=selection_id,
                    roster_type=roster_type, activation_status="active_primary" if active else "standby_not_triggered",
                    side=side, object_role="checksum_catalog", analysis_branch="core_metadata",
                    accession_version=accession, source_url=catalog_url, expected_bytes=len(catalog_body),
                    expected_local_sha256=sha256_bytes(catalog_body), upstream_checksum_algorithm="not_published",
                    upstream_checksum="", source_last_modified=catalog_headers.get("last-modified", ""),
                    pair_identity_status="pass", annotation_assembly_match="not_applicable",
                    selective_justification="provider digest catalog for exact approved assembly directory",
                ), payload=catalog_body, acquire=True)
                directory = prefix.rsplit("/", 1)[-1]
                fasta_name = directory + "_genomic.fna.gz"
                fasta_url = prefix + "/" + fasta_name
                fasta_size, fasta_modified = head_bytes(fasta_url)
                fasta_md5 = catalog.get(fasta_name, "")
                if not HEX32.fullmatch(fasta_md5):
                    raise AcquisitionRefusal(f"provider MD5 absent for {fasta_url}")
                add_object(base_object(
                    object_id=f"{selection_id}:{side}:genome_fasta", selection_id=selection_id,
                    roster_type=roster_type, activation_status="active_primary" if active else "standby_not_triggered",
                    side=side, object_role="genome_fasta", analysis_branch="core",
                    accession_version=accession, source_url=fasta_url, expected_bytes=fasta_size,
                    upstream_checksum_algorithm="md5", upstream_checksum=fasta_md5,
                    source_last_modified=fasta_modified, pair_identity_status="pass",
                    annotation_assembly_match="not_applicable",
                    selective_justification="exact H1/H2 assembly selected by immutable ten-pair roster",
                ), acquire=active)
                if active:
                    for role, suffix in (("assembly_report", "_assembly_report.txt"), ("assembly_stats", "_assembly_stats.txt")):
                        name = directory + suffix
                        url = prefix + "/" + name
                        size, modified = head_bytes(url)
                        digest = catalog.get(name, "")
                        if not HEX32.fullmatch(digest):
                            raise AcquisitionRefusal(f"provider MD5 absent for {url}")
                        add_object(base_object(
                            object_id=f"{selection_id}:{side}:{role}", selection_id=selection_id,
                            roster_type=roster_type, activation_status="active_primary", side=side,
                            object_role=role, analysis_branch="core_quality", accession_version=accession,
                            source_url=url, expected_bytes=size, upstream_checksum_algorithm="md5",
                            upstream_checksum=digest, source_last_modified=modified,
                            pair_identity_status="pass", annotation_assembly_match="not_applicable",
                            selective_justification="exact assembly structure/statistics evidence; no raw reads",
                        ))
                    repeat_names = [name for name in catalog if "repeat" in name.lower() or "rm.out" in name.lower()]
                    if repeat_names:
                        # Retain every published standalone repeat report for
                        # the exact accession; none was present in this freeze.
                        for repeat_index, name in enumerate(sorted(repeat_names), 1):
                            url = prefix + "/" + name
                            size, modified = head_bytes(url)
                            add_object(base_object(
                                object_id=f"{selection_id}:{side}:repeat_report:{repeat_index}", selection_id=selection_id,
                                roster_type=roster_type, activation_status="active_primary", side=side,
                                object_role="repeat_report", analysis_branch="repeat_quality_optional",
                                accession_version=accession, source_url=url, expected_bytes=size,
                                upstream_checksum_algorithm="md5", upstream_checksum=catalog[name],
                                source_last_modified=modified, pair_identity_status="pass",
                                annotation_assembly_match="not_applicable",
                                selective_justification="standalone repeat report published for exact assembly",
                            ))
                    else:
                        missing_url = prefix + "/" + directory + "_rm.out.gz"
                        add_object(base_object(
                            object_id=f"{selection_id}:{side}:repeat_report", selection_id=selection_id,
                            roster_type=roster_type, activation_status="active_primary", side=side,
                            object_role="repeat_report", analysis_branch="repeat_quality_optional",
                            accession_version=accession, source_url=missing_url, expected_bytes=0,
                            upstream_checksum_algorithm="not_published", upstream_checksum="",
                            pair_identity_status="pass", annotation_assembly_match="not_applicable",
                            selective_justification="exact-directory standalone repeat report checked; absence recorded without substitution",
                        ), acquire=False, nonacquire_status="missing")

            if active and annotation_report is not None and annotation_body is not None:
                add_object(base_object(
                    object_id=f"{selection_id}:annotation:dataset_report", selection_id=selection_id,
                    roster_type="primary", activation_status="active_primary", side="h1",
                    object_role="annotation_dataset_report", analysis_branch="annotation",
                    accession_version=frozen["h1_accession_version"], annotation_accession_version=annotation_name,
                    source_url=report_url(annotation_assembly), expected_bytes=len(annotation_body),
                    expected_local_sha256=sha256_bytes(annotation_body), upstream_checksum_algorithm="not_published",
                    upstream_checksum="", source_last_modified=annotation_headers.get("last-modified", ""),
                    pair_identity_status="pass", annotation_assembly_match=annotation_binding["annotation_branch_status"],
                    selective_justification="exact optional annotation accession/version provenance",
                ), payload=annotation_body)

            if active and annotation_assembly:
                source_report = reports["h1"] if annotation_assembly == frozen["h1_accession_version"] else annotation_report
                if source_report is not None and annotation_binding["annotation_branch_status"] != "failed_mismatch":
                    annotation_prefix = ftp_prefix(annotation_assembly, str(source_report["assembly_info"]["assembly_name"]))
                    md5_body, md5_headers = fetch_bytes(annotation_prefix + "/md5checksums.txt")
                    md5s = parse_md5_catalog(md5_body)
                    directory = annotation_prefix.rsplit("/", 1)[-1]
                    gff_name = directory + "_genomic.gff.gz"
                    if gff_name in md5s:
                        gff_url = annotation_prefix + "/" + gff_name
                        size, modified = head_bytes(gff_url)
                        add_object(base_object(
                            object_id=f"{selection_id}:h1:annotation_gff", selection_id=selection_id,
                            roster_type="primary", activation_status="active_primary", side="h1",
                            object_role="annotation_gff", analysis_branch="annotation",
                            accession_version=frozen["h1_accession_version"], annotation_accession_version=annotation_name,
                            source_url=gff_url, expected_bytes=size, upstream_checksum_algorithm="md5",
                            upstream_checksum=md5s[gff_name], source_last_modified=modified,
                            pair_identity_status="pass", annotation_assembly_match=annotation_binding["annotation_branch_status"],
                            selective_justification="optional exact native/paired annotation only; dictionary audit remains branch-local",
                        ))

            h1_busco = (annotation_report or reports["h1"]).get("annotation_info", {}).get("busco", {}) if (annotation_report or reports["h1"]) else {}
            qv_evidence = frozen["qv_evidence"] + "; no exact-sequence Merqury tabular report located in NCBI exact report/FTP; GenomeArk plots were not accepted without final-assembly digest binding"
            completeness = {"h1_annotation_busco": h1_busco or "missing", "h2_busco": "missing", "frozen_requirement": frozen["completeness_evidence"]}
            duplication = {"h1_annotation_busco_duplicated": h1_busco.get("duplicated", "missing") if isinstance(h1_busco, Mapping) else "missing", "h2_busco_duplicated": "missing", "frozen_requirement": frozen["duplication_collapse_evidence"]}
            pair_rows.append({
                "manifest_version": MANIFEST_VERSION,
                "selection_id": selection_id, "roster_type": roster_type, "rank": frozen["rank"],
                "activation_status": "active_primary" if active else "standby_not_triggered",
                "replacement_rule": frozen["alternate_replacement_rule"], "amendment_id": "none",
                "failed_primary_retained": "not_applicable" if active else "no trigger; all primaries retained",
                "catalog_row": frozen["catalog_row"], "catalog_commit": frozen["catalog_commit"],
                "catalog_sha256": frozen["catalog_sha256"],
                "release_membership": f"VGP Phase 1 Freeze 1.0 catalog row {frozen['catalog_row']}; both BioProject lineages include {VGP_RELEASE_PROJECT}",
                "species": frozen["species"], "taxid": frozen["resolved_taxid"], "biosample": frozen["biosample"],
                "individual_or_isolate": frozen["individual_or_isolate"],
                "h1_accession_version": frozen["h1_accession_version"], "h1_label": frozen["h1_label"], "h1_role": frozen["h1_role"],
                "h1_release_date": frozen["h1_date"], "h1_length_bp": frozen["h1_length_bp"],
                "h1_contigs": frozen["h1_contigs"], "h1_contig_n50_bp": frozen["h1_contig_n50_bp"],
                "h1_scaffold_n50_bp": frozen["h1_scaffold_n50_bp"], "h1_frozen_report_sha256": frozen["h1_report_sha256"],
                "h1_assembly_method": reports["h1"]["assembly_info"].get("assembly_method", ""),
                "h2_accession_version": frozen["h2_accession_version"], "h2_label": frozen["h2_label"], "h2_role": frozen["h2_role"],
                "h2_release_date": frozen["h2_date"], "h2_length_bp": frozen["h2_length_bp"],
                "h2_contigs": frozen["h2_contigs"], "h2_contig_n50_bp": frozen["h2_contig_n50_bp"],
                "h2_scaffold_n50_bp": frozen["h2_scaffold_n50_bp"], "h2_frozen_report_sha256": frozen["h2_report_sha256"],
                "h2_assembly_method": reports["h2"]["assembly_info"].get("assembly_method", ""),
                "reciprocal_linked_assembly_evidence": identity["reciprocal_linked_assembly_evidence"],
                "assembly_generation": frozen["assembly_generation"],
                "assembly_technologies": f"H1: {reports['h1']['assembly_info'].get('sequencing_tech','')}; H2: {reports['h2']['assembly_info'].get('sequencing_tech','')}",
                "hifi_evidence": frozen["read_hifi_provenance"], "long_range_phasing_evidence": frozen["long_range_phasing_evidence"],
                "repeat_report_evidence": "both exact NCBI accession directories and provider MD5 catalogs checked; no standalone repeat report published; repeat masking remains an explicit downstream missing input",
                "repeat_report_status": "missing_not_published",
                "qv_evidence": qv_evidence, "qv_status": "missing_exact_final_sequence_qv",
                "completeness_evidence": canonical_json(completeness), "completeness_status": "partial_h1_only" if h1_busco else "missing",
                "duplication_collapse_evidence": canonical_json(duplication), "duplication_collapse_status": "partial_h1_only" if h1_busco else "missing",
                "native_annotation_accession_version": annotation_name,
                "annotation_branch_status": annotation_binding["annotation_branch_status"],
                "annotation_binding_evidence": annotation_binding["annotation_binding_evidence"],
                "core_identity_status": "pass", "core_acquisition_status": "pending_object_accounting" if active else "not_activated",
                "raw_read_or_kmer_selection": "none",
                "selective_asset_justification": "No raw-read or k-mer payload was selected: exact published QC reports were sought first; unbounded SRA/GenomeArk raw archives are outside this manifest.",
                "license_or_access": NCBI_POLICY,
            })

    by_selection: dict[str, list[dict[str, Any]]] = {}
    for item in object_rows:
        by_selection.setdefault(str(item["selection_id"]), []).append(item)
    for pair in pair_rows:
        if pair["roster_type"] != "primary":
            continue
        core = [item for item in by_selection[pair["selection_id"]] if item["analysis_branch"] in {"core", "core_metadata", "core_quality"}]
        hard_failures = [item for item in core if item["status"] not in {"verified", "reused"}]
        pair["core_acquisition_status"] = "verified" if not hard_failures else "failed_closed:" + ",".join(str(item["object_id"]) for item in hard_failures)
        annotations = [item for item in by_selection[pair["selection_id"]] if item["analysis_branch"] == "annotation"]
        if annotations and any(item["status"] not in {"verified", "reused"} for item in annotations):
            pair["annotation_branch_status"] = "failed_acquisition_core_unaffected"

    return pair_rows, object_rows


def summarize_inventory(rows: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    terminal = ("reused", "missing", "superseded", "quarantined")
    newly_promoted = [row for row in rows if row.get("status") == "verified"]
    locally_verified = [row for row in rows if row.get("status") in {"verified", "reused"}]
    result: dict[str, Any] = {
        "planned": {"objects": len(rows), "bytes": sum(int(row.get("expected_bytes") or 0) for row in rows)},
        "transferred": {"objects": sum(1 for row in rows if int(row.get("transferred_bytes") or 0) > 0), "bytes": sum(int(row.get("transferred_bytes") or 0) for row in rows)},
        "verified": {"objects": len(locally_verified), "bytes": sum(int(row.get("expected_bytes") or 0) for row in locally_verified)},
        "newly_promoted": {"objects": len(newly_promoted), "bytes": sum(int(row.get("expected_bytes") or 0) for row in newly_promoted)},
    }
    for status in terminal:
        selected = [row for row in rows if row.get("status") == status]
        result[status] = {"objects": len(selected), "bytes": sum(int(row.get("expected_bytes") or 0) for row in selected)}
    reconciled_objects = result["newly_promoted"]["objects"] + sum(result[status]["objects"] for status in terminal)
    reconciled_bytes = result["newly_promoted"]["bytes"] + sum(result[status]["bytes"] for status in terminal)
    result["terminal_disposition_reconciliation"] = {
        "objects": reconciled_objects, "bytes": reconciled_bytes,
        "matches_planned": reconciled_objects == result["planned"]["objects"] and reconciled_bytes == result["planned"]["bytes"],
        "note": "newly_promoted/reused/missing/superseded/quarantined are mutually exclusive; transferred is an I/O-flow measure and verified is a validation-state measure, so neither is added to terminal-disposition bytes",
    }
    return result


def mechanism_demo(data_root: Path) -> dict[str, Any]:
    demo_root = data_root / "pilot/runs/acquisition-mechanism-demo-v1"
    if demo_root.exists():
        shutil.rmtree(demo_root)
    payload = b"0123456789abcdef"
    row = base_object(
        object_id="DEMO:resume-promote-reuse", selection_id="DEMO", object_role="mechanism_demo",
        source_url="https://example.invalid/immutable-demo-v1", expected_bytes=len(payload),
        upstream_checksum_algorithm="md5", upstream_checksum=hashlib.md5(payload).hexdigest(),  # noqa: S324
        expected_local_sha256=sha256_bytes(payload),
    )
    partial = partial_path(demo_root, str(row["object_id"]))
    partial.parent.mkdir(parents=True, exist_ok=True)
    partial.write_bytes(payload[:5])
    first = acquire_one(row, demo_root, downloader=_inline_downloader(payload))
    second = acquire_one(row, demo_root, downloader=lambda *_: (_ for _ in ()).throw(AssertionError("mirror reuse downloaded")))
    bad_row = dict(row)
    bad_row["object_id"] = "DEMO:quarantine"
    bad_row["expected_local_sha256"] = sha256_bytes(b"different")
    bad = acquire_one(bad_row, demo_root, downloader=_inline_downloader(payload))
    evidence = {
        "resume": {"status": first["status"], "resume_from_bytes": first["resume_from_bytes"], "transferred_bytes": first["transferred_bytes"]},
        "atomic_promotion_and_reverification": {"status": first["status"], "local_path": first["local_path"], "revalidation_count": first["revalidation_count"]},
        "mirror_reuse": {"status": second["status"], "transferred_bytes": second["transferred_bytes"], "revalidation_count": second["revalidation_count"]},
        "quarantine": {"status": bad["status"], "quarantine_path": bad["quarantine_path"]},
    }
    atomic_write(demo_root / "evidence.json", (json.dumps(evidence, indent=2, sort_keys=True) + "\n").encode())
    return evidence


def acquire_direct_control(data_root: Path) -> list[dict[str, Any]]:
    """Acquire the bounded public D01 control objects, not its 357-GB read archive."""

    fields = (
        "study_accession,secondary_study_accession,run_accession,experiment_accession,"
        "sample_accession,secondary_sample_accession,sample_title,experiment_title,"
        "instrument_platform,instrument_model,library_name,library_strategy,library_source,"
        "library_selection,fastq_ftp,fastq_bytes,fastq_md5,submitted_ftp,submitted_bytes,submitted_md5"
    )
    query = urllib.parse.urlencode(
        {"accession": "PRJEB4500", "result": "read_run", "fields": fields, "format": "tsv", "download": "true"}
    )
    ena_url = "https://www.ebi.ac.uk/ena/portal/api/filereport?" + query
    ena_body, ena_headers = fetch_bytes(ena_url)
    runs = list(csv.DictReader(ena_body.decode("utf-8").splitlines(), delimiter="\t"))
    if len(runs) != 472 or {row["study_accession"] for row in runs} != {"PRJEB4500"} or {row["secondary_study_accession"] for row in runs} != {"ERP003793"}:
        raise AcquisitionRefusal("D01 ENA expansion drifted from frozen 472-run PRJEB4500/ERP003793 inventory")
    tetrad_ids = ("29", "38", "40", "58", "62", "34", "35", "39", "51", "53", "68", "69", "76")
    expected_products = {f"{tetrad}_{product}" for tetrad in tetrad_ids for product in range(1, 5)}
    sample_names = {row["library_name"] for row in runs}
    if not expected_products.issubset(sample_names):
        raise AcquisitionRefusal("D01 does not contain all four products for each of the 13 frozen tetrads")
    parent_names = {"Col-0", "Ler_WUR", "Cvi_NEW_correct"}
    if not parent_names.issubset(sample_names):
        raise AcquisitionRefusal("D01 exact Col/Ler/Cvi parental libraries are incomplete")
    raw_objects = sum(len(row["fastq_ftp"].split(";")) for row in runs if row["fastq_ftp"])
    raw_bytes = sum(sum(int(value) for value in row["fastq_bytes"].split(";")) for row in runs if row["fastq_bytes"])
    if raw_objects != 951 or raw_bytes != 356_984_558_868:
        raise AcquisitionRefusal("D01 public raw archive object/byte aggregate drift")
    selected_names = expected_products | parent_names | {"29_4_batch4"}
    selected_runs = [row for row in runs if row["library_name"] in selected_names]
    selected_raw_objects = sum(len(row["fastq_ftp"].split(";")) for row in selected_runs if row["fastq_ftp"])
    selected_raw_bytes = sum(sum(int(value) for value in row["fastq_bytes"].split(";")) for row in selected_runs if row["fastq_bytes"])
    relationships = "; ".join(f"tetrad {name}: {name}_1,{name}_2,{name}_3,{name}_4" for name in tetrad_ids)
    relationships += "; parents: Col-0,Ler_WUR,Cvi_NEW_correct; 29_4_batch4 is a technical library for product 29_4"
    prior = {
        row["source_url"]: row["local_sha256"]
        for row in read_tsv(DIRECT_OUTPUT) if DIRECT_OUTPUT.is_file() and row.get("local_sha256")
    } if DIRECT_OUTPUT.is_file() else {}
    common = {
        "manifest_version": MANIFEST_VERSION, "control_id": "D01", "selection_status": "selected",
        "dataset_title": "Arabidopsis recombinant tetrads and doubled haploids",
        "species": "Arabidopsis thaliana", "taxid": "3702", "repository": "ENA;eLife",
        "dataset_accession_version": "PRJEB4500;ERP003793;doi:10.7554/eLife.01426.v1",
        "relationships": relationships,
        "complete_pedigree_or_gamete_evidence": "13 complete four-product male meiotic tetrads (52 products) with exact Col/Ler/Cvi parental libraries; not inferred from trio labels",
        "transmitted_haplotype_evidence": "all four products reconstruct reciprocal Col/Ler haplotypes; published 2:2, 3:1, and 1:3 marker segregation tables permit directional conversion adjudication",
        "directional_transmission_suitability": "suitable for 13 parent-to-four-product meioses; GC-bias remains pilot/descriptive at this sample size",
        "consent": "non-human plant material; human-participant consent not applicable",
        "access": "public ENA study and open eLife article",
        "license_or_terms": "ENA/EMBL-EBI terms of use and study citation; eLife article/supplements CC BY 3.0; third-party rights not transferred",
    }
    result: list[dict[str, Any]] = []

    def add(role: str, url: str, payload: bytes | None = None, headers: Mapping[str, str] | None = None) -> None:
        if payload is None:
            size, modified = head_bytes(url)
            row = base_object(
                object_id=f"D01:{role}", selection_id="D01", object_role=role, source_url=url,
                expected_bytes=size, expected_local_sha256=prior.get(url, ""),
                upstream_checksum_algorithm="not_published", upstream_checksum="",
                source_last_modified=modified,
            )
            acquired = acquire_one(row, data_root)
        else:
            row = base_object(
                object_id=f"D01:{role}", selection_id="D01", object_role=role, source_url=url,
                expected_bytes=len(payload), expected_local_sha256=prior.get(url, sha256_bytes(payload)),
                upstream_checksum_algorithm="not_published", upstream_checksum="",
                source_last_modified=(headers or {}).get("last-modified", ""),
            )
            acquired = acquire_one(row, data_root, downloader=_inline_downloader(payload))
        direct = {field: "" for field in DIRECT_FIELDS}
        direct.update(common)
        direct.update(
            object_id=acquired["object_id"], object_role=role, source_url=url,
            expected_bytes=acquired["expected_bytes"], upstream_checksum_algorithm=acquired["upstream_checksum_algorithm"],
            upstream_checksum=acquired["upstream_checksum"], local_sha256=acquired["local_sha256"],
            observed_bytes=acquired["observed_bytes"], transferred_bytes=acquired["transferred_bytes"],
            retrieval_utc=acquired["retrieval_utc"], status=acquired["status"], local_path=acquired["local_path"],
            selective_justification="small immutable metadata/published event or tract evidence selected instead of the complete bulk raw-read archive",
        )
        result.append(direct)

    add("ena_exact_run_object_manifest", ena_url, ena_body, ena_headers)
    article_url = "https://api.elifesciences.org/articles/01426"
    article_body, article_headers = fetch_bytes(article_url)
    article = json.loads(article_body)
    if article.get("doi") != "10.7554/eLife.01426" or article.get("version") != 1:
        raise AcquisitionRefusal("D01 article version drift")
    add("article_version_metadata", article_url, article_body, article_headers)
    additional = {item["filename"]: item["uri"] for item in article.get("additionalFiles", [])}
    for filename, role in (
        ("elife-01426-supp1-v1.xlsx", "published_sample_and_event_tables"),
        ("elife-01426-supp2-v1.xlsx", "published_marker_filter_and_validation_tables"),
        ("elife-01426-supp3-v1.pdf", "published_methods_supplement"),
    ):
        if filename not in additional:
            raise AcquisitionRefusal(f"D01 immutable article metadata lacks {filename}")
        add(role, additional[filename])

    reference, reference_body, reference_headers = report_for("GCF_000001735.4")
    add("reference_dataset_report", report_url("GCF_000001735.4"), reference_body, reference_headers)
    ref_prefix = ftp_prefix("GCF_000001735.4", str(reference["assembly_info"]["assembly_name"]))
    md5_body, md5_headers = fetch_bytes(ref_prefix + "/md5checksums.txt")
    md5s = parse_md5_catalog(md5_body)
    ref_dir = ref_prefix.rsplit("/", 1)[-1]
    for suffix, role in (("_genomic.fna.gz", "tair10_reference_fasta"), ("_genomic.gff.gz", "tair10_refseq_annotation"), ("_assembly_report.txt", "tair10_assembly_report")):
        name = ref_dir + suffix
        url = ref_prefix + "/" + name
        size, modified = head_bytes(url)
        digest = md5s.get(name, "")
        if not HEX32.fullmatch(digest):
            raise AcquisitionRefusal(f"D01 reference provider MD5 missing for {url}")
        row = base_object(
            object_id=f"D01:{role}", selection_id="D01", object_role=role, source_url=url,
            expected_bytes=size, expected_local_sha256=prior.get(url, ""), upstream_checksum_algorithm="md5",
            upstream_checksum=digest, source_last_modified=modified,
        )
        acquired = acquire_one(row, data_root)
        direct = {field: "" for field in DIRECT_FIELDS}
        direct.update(common)
        direct.update(
            object_id=acquired["object_id"], object_role=role, source_url=url,
            expected_bytes=acquired["expected_bytes"], upstream_checksum_algorithm="md5",
            upstream_checksum=digest, local_sha256=acquired["local_sha256"], retrieval_utc=acquired["retrieval_utc"],
            observed_bytes=acquired["observed_bytes"], transferred_bytes=acquired["transferred_bytes"],
            status=acquired["status"], local_path=acquired["local_path"],
            selective_justification="exact TAIR10 reference/annotation required by D01; no alternative reference",
        )
        result.append(direct)

    excluded = {field: "" for field in DIRECT_FIELDS}
    excluded.update(common)
    excluded.update(
        object_id="D01:bulk_raw_archive_not_selected", object_role="bulk_raw_archive_exclusion_aggregate",
        source_url="https://www.ebi.ac.uk/ena/browser/view/PRJEB4500", expected_bytes=raw_bytes,
        upstream_checksum_algorithm="per_object_md5_in_acquired_ena_manifest", upstream_checksum="951 exact MD5 values",
        status="superseded", access="public; deliberately not transferred",
        selective_justification=(
            f"excluded all {raw_objects} FASTQ objects/{raw_bytes} bytes; even the bounded 13-tetrad+3-parent subset is "
            f"{selected_raw_objects} objects/{selected_raw_bytes} bytes. Published event/filter tables and exact run metadata "
            "are the selected stratified validation assets; zero FASTQs were acquired."
        ),
    )
    result.append(excluded)
    return result


def handoff_text(pair_rows: Sequence[Mapping[str, Any]], object_rows: Sequence[Mapping[str, Any]], direct_rows: Sequence[Mapping[str, Any]], summary: Mapping[str, Any]) -> str:
    primary = [row for row in pair_rows if row["roster_type"] == "primary"]
    alternates = [row for row in pair_rows if row["roster_type"] == "alternate"]
    core_verified = sum(row["core_acquisition_status"] == "verified" for row in primary)
    direct = direct_rows[0]
    account = summary["accounting"]
    return f"""# VGP ten-pair pilot acquisition handoff

**Manifest version:** `{MANIFEST_VERSION}`

**Acquisition UTC:** `{summary['generated_at_utc']}`

**Execution:** pinned GNU Guix commit `44bbfc24e4bcc48d0e3343cd3d83452721af8c36`; local process only; zero Slurm jobs.

## Outcome

The closed world contains exactly {len(primary)} approved primaries and {len(alternates)} pre-ranked alternates. Exact species, TaxId, catalog row, BioSample, isolate, accession.version, reciprocal linked-assembly roles, VGP BioProject lineage, technologies, and annotation disposition are frozen in `vgp_10_pilot_acquisition_manifest.tsv`. {core_verified}/10 primary core acquisition rows have every required core object verified or reused. An alternate was not activated: all alternate biological payload rows remain `superseded` by the retained primary, and no manifest amendment exists.

The direct-control disposition is `{direct.get('selection_status')}` / `{direct.get('status')}` for `{direct.get('control_id')}`. D01 is bound to 13 complete four-product tetrads and exact Col/Ler/Cvi parents. No trio label or alternate dataset is treated as a complete pedigree, and the 357-GB ENA FASTQ archive is explicitly superseded by the selected immutable event/filter supplements and run metadata.

## Exact accounting

| category | objects | bytes |
|---|---:|---:|
| planned | {account['planned']['objects']} | {account['planned']['bytes']} |
| transferred (I/O flow; orthogonal) | {account['transferred']['objects']} | {account['transferred']['bytes']} |
| verified by local SHA-256 (validation state; orthogonal) | {account['verified']['objects']} | {account['verified']['bytes']} |
| newly promoted | {account['newly_promoted']['objects']} | {account['newly_promoted']['bytes']} |
| reused | {account['reused']['objects']} | {account['reused']['bytes']} |
| missing | {account['missing']['objects']} | {account['missing']['bytes']} |
| superseded/not activated | {account['superseded']['objects']} | {account['superseded']['bytes']} |
| quarantined | {account['quarantined']['objects']} | {account['quarantined']['bytes']} |

Newly promoted/reused/missing/superseded/quarantined are mutually exclusive and reconcile to planned objects and bytes: `{account['terminal_disposition_reconciliation']['matches_planned']}`. Transferred bytes are an orthogonal physical-I/O measure, and verified bytes are an orthogonal validation-state measure; neither is added to terminal-disposition bytes.

## Fail-closed branch status

Same-individual or H1/H2 ambiguity, accession drift, size mismatch, provider-MD5 mismatch, and local-SHA mismatch are hard core refusals. Annotation absence, a paired-RefSeq dictionary difference, or annotation download failure is branch-local and cannot veto a valid core. Published H1 annotation BUSCO values are retained where exact; missing H2 BUSCO and missing exact-final-sequence Merqury QV are explicitly missing, never imputed. Hi-C absence alone is not used as a refusal.

No raw-read or k-mer payload is in the object inventory. The resolver never lists SRA runs or GenomeArk raw-data prefixes. This proves that zero unmanifested bulk raw-read objects/bytes were acquired; the choice deliberately leaves selective validation/QV evidence incomplete rather than expanding scope silently.

## Storage protocol demonstrated

`summary.json` records a same-code-path demonstration of partial resume, size/provider digest checks, pre/post-promotion SHA-256, atomic `os.replace`, read-only promotion, CAS reuse with local SHA-256 revalidation, and mismatch quarantine. The demonstration accounting is one verified/reused 16-byte object plus one quarantined 16-byte candidate and is deliberately separate from biological-object accounting. Real objects use the same `objects/sha256/<2>/<2>/<sha256>` contract and are directly reusable by the official Freeze 1 mirror.
"""


def run(args: argparse.Namespace) -> int:
    data_root = args.data_root.resolve()
    for relative in ("objects/sha256", "staging/partials", "quarantine", "pilot/runs"):
        (data_root / relative).mkdir(parents=True, exist_ok=True)
    pair_rows, object_rows = resolve_and_acquire(data_root)
    # Publish the complete VGP inventory checkpoint before the independent
    # direct-control branch.  A branch-local provider failure must not discard
    # URL-to-CAS bindings or force valid core objects to be transferred again.
    write_tsv(PAIR_OUTPUT, PAIR_FIELDS, pair_rows)
    write_tsv(OBJECT_OUTPUT, OBJECT_FIELDS, object_rows)
    direct_rows = acquire_direct_control(data_root)
    vgp_accounting = summarize_inventory(object_rows)
    direct_accounting = summarize_inventory(direct_rows)
    accounting = summarize_inventory([*object_rows, *direct_rows])
    evidence = mechanism_demo(data_root)
    summary = {
        "schema_version": "1.0.0", "manifest_version": MANIFEST_VERSION, "generated_at_utc": utc_now(),
        "data_root": str(data_root), "primary_slots": 10, "alternate_slots": len([r for r in pair_rows if r["roster_type"] == "alternate"]),
        "activated_alternates": 0, "manifest_amendment": "none", "accounting": accounting,
        "vgp_accounting": vgp_accounting, "direct_control_accounting": direct_accounting,
        "mechanism_demonstration": evidence,
        "mechanism_demonstration_accounting": {
            "resume_candidates": {"objects": 1, "bytes": 16},
            "atomically_promoted": {"objects": 1, "bytes": 16},
            "mirror_reused_after_sha256_revalidation": {"objects": 1, "bytes": 16},
            "quarantined_mismatch_candidates": {"objects": 1, "bytes": 16},
            "included_in_biological_accounting": False,
        },
        "scope_proof": {
            "manifested_raw_read_or_kmer_objects": sum(row["object_role"] in BULK_RAW_ROLES for row in object_rows),
            "unmanifested_bulk_raw_read_objects_acquired": 0, "unmanifested_bulk_raw_read_bytes_acquired": 0,
            "sra_run_enumeration_requests": 0, "ena_run_metadata_manifest_requests": 1,
            "ena_fastq_payload_requests": 0, "genomeark_raw_prefix_requests": 0, "slurm_jobs_submitted": 0,
            "command_scope": "versioned NCBI assembly/report URLs and the selected direct-control manifest only",
        },
        "environment": {
            "guix_commit": "44bbfc24e4bcc48d0e3343cd3d83452721af8c36",
            "channels_sha256": sha256_bytes((ROOT / "analysis/guix/channels.scm").read_bytes()),
            "manifest_sha256": sha256_bytes((ROOT / "analysis/guix/manifest.scm").read_bytes()),
            "slurm_environment_detected": False,
        },
        "direct_control": {
            "control_id": direct_rows[0].get("control_id"), "selection_status": direct_rows[0].get("selection_status"),
            "verified_or_reused_objects": sum(row.get("status") in {"verified", "reused"} for row in direct_rows),
            "missing_objects": sum(row.get("status") == "missing" for row in direct_rows),
            "raw_archive_disposition": "superseded_not_selected",
        },
    }
    if not accounting["terminal_disposition_reconciliation"]["matches_planned"]:
        raise AcquisitionRefusal("terminal object/byte accounting does not reconcile")
    write_tsv(DIRECT_OUTPUT, DIRECT_FIELDS, direct_rows)
    atomic_write(SUMMARY_OUTPUT, (json.dumps(summary, indent=2, sort_keys=True) + "\n").encode())
    atomic_write(HANDOFF_OUTPUT, handoff_text(pair_rows, object_rows, direct_rows, summary).encode())
    return 0


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--data-root", type=Path, default=DEFAULT_DATA_ROOT)
    return result


if __name__ == "__main__":
    try:
        raise SystemExit(run(parser().parse_args()))
    except AcquisitionRefusal as error:
        print(f"REFUSED: {error}", file=sys.stderr)
        raise SystemExit(2)
