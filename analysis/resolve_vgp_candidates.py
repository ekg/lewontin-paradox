#!/usr/bin/env python3
"""Repair the bounded VGP pilot resolver from immutable official metadata.

The default operation is offline and deterministic: it reads the integrated
frozen VGP catalog plus response objects already present in
``analysis/vgp_resolution_cache`` and regenerates the three pilot TSVs and the
cache index.  ``--refresh-cache`` performs metadata-only NCBI requests.  It
never downloads a FASTA, GFF, population payload, or other biological object.

The resolver deliberately distinguishes metadata eligibility from execution
authorization.  A row may have an exact NCBI assembly and native annotation
while remaining unselected because a stricter integrated cap is unresolved.
"""

from __future__ import annotations

import argparse
import csv
import email.utils
import hashlib
import json
import os
import platform
import random
import re
import shutil
import sys
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List, Mapping, MutableMapping, Optional, Sequence, Tuple

_IMPORT_ROOT = Path(__file__).resolve().parents[1]
if str(_IMPORT_ROOT) not in sys.path:
    sys.path.insert(0, str(_IMPORT_ROOT))

from analysis import freeze_vgp_manifest as freeze


PROJECT_ROOT = Path(__file__).resolve().parents[1]
CACHE_ROOT = PROJECT_ROOT / "analysis" / "vgp_resolution_cache"
PRIOR_ROOT = CACHE_ROOT / "prior_refusal"
DEFAULT_CATALOG = Path(
    "/moosefs/erikg/vgp/manifests/"
    "VGPPhase1-freeze-1.0.commit-dc1b2af5a7741b97d66fb10cb2bce97f41765cdf.tsv"
)
DEFAULT_MANIFEST = PROJECT_ROOT / "analysis" / "vgp_pilot_manifest.tsv"
DEFAULT_REJECTIONS = PROJECT_ROOT / "analysis" / "vgp_pilot_rejections.tsv"
DEFAULT_BUDGET = PROJECT_ROOT / "analysis" / "vgp_pilot_size_budget.tsv"
DEFAULT_ROOT_VALIDATION = PROJECT_ROOT / "analysis" / "vgp_data_root_validation.json"
DEFAULT_INDEX = CACHE_ROOT / "index.json"
DEFAULT_CHECKPOINT = CACHE_ROOT / "checkpoint.json"

CATALOG_COMMIT = "dc1b2af5a7741b97d66fb10cb2bce97f41765cdf"
CATALOG_SHA256 = "9c58420484a8b76a2d6175b7c26bf709e68bdc726a67fc7541b8c2b5a2fc13a4"
NCBI_DATASETS_VERSION = "18.33.1"
SOFTWARE_VERSION = "vgp-candidate-resolver/2.0"
USER_AGENT = "lewontin-paradox-vgp-candidate-resolver/2.0 (metadata-only)"

PRIORITIZED_CANDIDATE_IDS = (
    "camelus_dromedarius_gca_036321535_1",
    "colius_striatus_gca_028858725_2",
    "candoia_aspera_gca_035149785_1",
    "dendropsophus_ebraccatus_gca_027789765_1",
    "lepisosteus_oculatus_gca_040954835_1",
    "heterodontus_francisci_gca_036365525_1",
)

REPAIRED_COLUMNS = (
    "resolution_priority",
    "resolution_status",
    "resolution_reason_codes",
    "metadata_cache_request_keys",
    "h1_exact_version_status",
    "taxon_identity_status",
    "annotation_sequence_region_linkage_status",
    "phase_evidence_status",
    "resolved_modality",
    "acquisition_obligations",
    "post_alignment_measurement_contract_id",
    "post_alignment_result_disposition",
)

BUDGET_CAPS = {
    "species": 6,
    "compressed_inputs_gib": 120.0,
    "scratch_gib": 750.0,
    "core_hours": 1500.0,
    "concurrent_species": 2,
    "memory_per_job_gib": 256.0,
}

POST_ALIGNMENT_MEASUREMENT_CONTRACT = {
    "contract_id": "vgp_post_alignment_denominators_v1",
    "phase": "post_alignment_pre_result_acceptance",
    "executable_command": (
        "python3 analysis/resolve_vgp_candidates.py measure "
        "--metrics-json METRICS.json --output-json ACCEPTANCE.json"
    ),
    "inputs": {
        "callable_bases": "integer bases passing the pinned alignment/callability mask",
        "callable_fraction": "callable_bases divided by exact H1 assembly length",
        "queryable_gene_count": "protein-coding genes with at least one queryable CDS base",
        "queryable_gene_bases": "union of queryable protein-coding CDS bases",
    },
    "minimum_thresholds": {
        "callable_bases": 10_000_000,
        "callable_fraction": 0.50,
        "queryable_gene_count": 1_000,
        "queryable_gene_bases": 1_000_000,
    },
    "failure_disposition": "exclude_downstream_result",
    "pre_download_prerequisite": False,
}


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def canonical_json_bytes(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("utf-8")


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _atomic_write(path: Path, value: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix=path.name + ".", dir=str(path.parent))
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(value)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, str(path))
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def catalog_statistics(text: str) -> Dict[str, Any]:
    physical_lines = len(text.splitlines())
    rows = list(csv.DictReader(text.splitlines(), delimiter="\t"))
    counts = Counter(
        row.get("Scientific Name", "").strip()
        for row in rows
        if row.get("Scientific Name", "").strip()
    )
    duplicates = [
        {"scientific_name": name, "multiplicity": multiplicity}
        for name, multiplicity in sorted(counts.items())
        if multiplicity > 1
    ]
    return {
        "physical_lines": physical_lines,
        "header_lines": 1 if physical_lines else 0,
        "data_rows": len(rows),
        "unique_species": len(counts),
        "data_row_excess_over_unique_species": len(rows) - len(counts),
        "duplicated_species": duplicates,
    }


def _normalize(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {str(key): _normalize(value[key]) for key in sorted(value)}
    if isinstance(value, (list, tuple)):
        normalized = [_normalize(item) for item in value]
        if all(isinstance(item, (str, int, float, bool)) or item is None for item in normalized):
            return sorted(normalized, key=lambda item: json.dumps(item, sort_keys=True))
        return normalized
    return value


@dataclass(frozen=True)
class NormalizedRequest:
    source: str
    source_version: str
    method: str
    endpoint: str
    parameters: Mapping[str, Any]

    def __post_init__(self) -> None:
        object.__setattr__(self, "method", self.method.upper())
        object.__setattr__(self, "endpoint", self.endpoint.rstrip("/"))
        object.__setattr__(self, "parameters", _normalize(dict(self.parameters)))

    def as_dict(self) -> Dict[str, Any]:
        return {
            "source": self.source,
            "source_version": self.source_version,
            "method": self.method,
            "endpoint": self.endpoint,
            "parameters": self.parameters,
        }

    @property
    def cache_key(self) -> str:
        return sha256_bytes(canonical_json_bytes(self.as_dict()))


class CacheConflictError(RuntimeError):
    pass


class ImmutableResponseCache:
    def __init__(self, root: Path):
        self.root = Path(root)

    def _paths(self, request: NormalizedRequest) -> Tuple[Path, Path]:
        return (
            self.root / "responses" / (request.cache_key + ".bin"),
            self.root / "entries" / (request.cache_key + ".json"),
        )

    def load(self, request: NormalizedRequest) -> Optional[Tuple[bytes, Dict[str, Any]]]:
        response_path, entry_path = self._paths(request)
        if not response_path.is_file() or not entry_path.is_file():
            return None
        body = response_path.read_bytes()
        entry = json.loads(entry_path.read_text(encoding="utf-8"))
        if entry.get("request_key") != request.cache_key:
            raise CacheConflictError("cache request key does not match entry")
        if entry.get("request") != request.as_dict():
            raise CacheConflictError("cache normalized request does not match entry")
        if sha256_bytes(body) != entry.get("response_sha256"):
            raise CacheConflictError("cache response digest does not match bytes")
        return body, entry

    def store(
        self,
        request: NormalizedRequest,
        body: bytes,
        *,
        response_headers: Mapping[str, str],
        retrieved_at_utc: str,
        software_environment: Mapping[str, Any],
    ) -> Dict[str, Any]:
        response_path, entry_path = self._paths(request)
        existing = self.load(request)
        if existing is not None:
            if existing[0] != body:
                raise CacheConflictError("immutable cache conflict for " + request.cache_key)
            return existing[1]
        if response_path.exists() or entry_path.exists():
            raise CacheConflictError(
                "immutable cache conflict: incomplete object already exists for "
                + request.cache_key
            )
        try:
            relative_response_path = str(response_path.relative_to(PROJECT_ROOT))
        except ValueError:
            relative_response_path = str(response_path)
        entry = {
            "request_key": request.cache_key,
            "request": request.as_dict(),
            "endpoint": request.endpoint,
            "retrieved_at_utc": retrieved_at_utc,
            "response_sha256": sha256_bytes(body),
            "response_size_bytes": len(body),
            "response_headers": {
                str(key).lower(): str(value) for key, value in sorted(response_headers.items())
            },
            "software_environment": dict(software_environment),
            "response_path": relative_response_path,
        }
        _atomic_write(response_path, body)
        _atomic_write(entry_path, json.dumps(entry, indent=2, sort_keys=True).encode("utf-8") + b"\n")
        return entry

    def entries(self) -> List[Dict[str, Any]]:
        if not (self.root / "entries").is_dir():
            return []
        return [
            json.loads(path.read_text(encoding="utf-8"))
            for path in sorted((self.root / "entries").glob("*.json"))
        ]


@dataclass
class TransportResponse:
    status: int
    body: bytes
    headers: Mapping[str, str]


def _request_url(request: NormalizedRequest) -> str:
    query = request.parameters.get("_query", {})
    if not query:
        return request.endpoint
    return request.endpoint + "?" + urllib.parse.urlencode(query)


def urllib_transport(request: NormalizedRequest) -> TransportResponse:
    headers = {"User-Agent": USER_AGENT, "Accept": "application/json"}
    wire = urllib.request.Request(_request_url(request), headers=headers, method=request.method)
    try:
        with urllib.request.urlopen(wire, timeout=90) as response:
            return TransportResponse(response.status, response.read(), dict(response.headers.items()))
    except urllib.error.HTTPError as error:
        return TransportResponse(error.code, error.read(), dict(error.headers.items()))


class BatchedCachedClient:
    RETRYABLE = frozenset((429, 500, 502, 503, 504))

    def __init__(
        self,
        *,
        cache: ImmutableResponseCache,
        checkpoint_path: Path,
        transport: Callable[[NormalizedRequest], TransportResponse],
        sleep: Callable[[float], None] = time.sleep,
        jitter: Callable[[int], float] = lambda attempt: random.uniform(0.0, min(1.0, 0.1 * attempt)),
        minimum_interval_seconds: float = 0.34,
        maximum_attempts: int = 6,
        maximum_backoff_seconds: float = 30.0,
        software_environment: Optional[Mapping[str, Any]] = None,
    ):
        self.cache = cache
        self.checkpoint_path = Path(checkpoint_path)
        self.transport = transport
        self.sleep = sleep
        self.jitter = jitter
        self.minimum_interval_seconds = minimum_interval_seconds
        self.maximum_attempts = maximum_attempts
        self.maximum_backoff_seconds = maximum_backoff_seconds
        self.software_environment = dict(software_environment or environment_identity())
        self._last_request_finished: Optional[float] = None

    def _checkpoint(self) -> Dict[str, Any]:
        if self.checkpoint_path.is_file():
            return json.loads(self.checkpoint_path.read_text(encoding="utf-8"))
        return {"schema_version": 1, "completed_request_keys": [], "updated_at_utc": None}

    def _mark_complete(self, request_key: str) -> None:
        checkpoint = self._checkpoint()
        keys = set(checkpoint.get("completed_request_keys", []))
        keys.add(request_key)
        checkpoint["completed_request_keys"] = sorted(keys)
        checkpoint["updated_at_utc"] = utc_now()
        _atomic_write(
            self.checkpoint_path,
            json.dumps(checkpoint, indent=2, sort_keys=True).encode("utf-8") + b"\n",
        )

    @staticmethod
    def _retry_after(value: Optional[str]) -> Optional[float]:
        if not value:
            return None
        try:
            return max(0.0, float(value))
        except ValueError:
            parsed = email.utils.parsedate_to_datetime(value)
            now = datetime.now(parsed.tzinfo or timezone.utc)
            return max(0.0, (parsed - now).total_seconds())

    def fetch(self, request: NormalizedRequest) -> Tuple[bytes, Dict[str, Any]]:
        cached = self.cache.load(request)
        if cached is not None:
            self._mark_complete(request.cache_key)
            return cached
        last_status = None
        for attempt in range(1, self.maximum_attempts + 1):
            if self._last_request_finished is not None and self.minimum_interval_seconds > 0:
                elapsed = time.monotonic() - self._last_request_finished
                if elapsed < self.minimum_interval_seconds:
                    self.sleep(self.minimum_interval_seconds - elapsed)
            response = self.transport(request)
            self._last_request_finished = time.monotonic()
            last_status = response.status
            if 200 <= response.status < 300:
                entry = self.cache.store(
                    request,
                    response.body,
                    response_headers=response.headers,
                    retrieved_at_utc=utc_now(),
                    software_environment=self.software_environment,
                )
                self._mark_complete(request.cache_key)
                return response.body, entry
            if response.status not in self.RETRYABLE or attempt == self.maximum_attempts:
                raise RuntimeError(
                    "metadata request failed: status={0} endpoint={1}".format(
                        response.status, request.endpoint
                    )
                )
            retry_after = self._retry_after(
                response.headers.get("Retry-After") or response.headers.get("retry-after")
            )
            delay = retry_after if retry_after is not None else min(
                self.maximum_backoff_seconds, 2.0 ** (attempt - 1)
            )
            self.sleep(min(self.maximum_backoff_seconds, delay) + self.jitter(attempt))
        raise RuntimeError("metadata request exhausted retries: status={0}".format(last_status))

    def fetch_accession_batches(
        self,
        *,
        accessions: Sequence[str],
        batch_size: int,
        request_factory: Callable[[Tuple[str, ...]], NormalizedRequest],
    ) -> List[Tuple[bytes, Dict[str, Any]]]:
        normalized = sorted(set(accession.strip() for accession in accessions if accession.strip()))
        results = []
        for offset in range(0, len(normalized), batch_size):
            batch = tuple(normalized[offset : offset + batch_size])
            results.append(self.fetch(request_factory(batch)))
        return results


def environment_identity() -> Dict[str, Any]:
    return {
        "software": SOFTWARE_VERSION,
        "python": platform.python_version(),
        "implementation": platform.python_implementation(),
        "platform": platform.platform(),
        "guix_channels_sha256": sha256_file(PROJECT_ROOT / "analysis/guix/channels.scm"),
        "guix_manifest_sha256": sha256_file(PROJECT_ROOT / "analysis/guix/manifest.scm"),
    }


def payload_acquisition_obligations(
    *,
    accession_version: str,
    url: str,
    expected_size_bytes: Optional[int],
    official_checksum_algorithm: Optional[str],
    official_checksum: Optional[str],
) -> Dict[str, Any]:
    exact_accession = bool(re.match(r"^GC[AF]_\d+\.\d+$", accession_version or ""))
    deterministic_url = bool(re.match(r"^https://\S+$", url or ""))
    finite_size = isinstance(expected_size_bytes, int) and expected_size_bytes > 0
    checksum_present = bool(official_checksum_algorithm and official_checksum)
    steps = []
    if exact_accession and deterministic_url and finite_size:
        steps = [
            "stage_full_payload",
            "verify_expected_size_bytes:{0}".format(expected_size_bytes),
            "compute_local_sha256",
            "reverify_local_sha256_before_promotion",
        ]
        if checksum_present:
            steps.append(
                "verify_official_{0}:{1}".format(
                    official_checksum_algorithm.lower(), official_checksum.lower()
                )
            )
        steps.append("atomic_promote_read_only")
    return {
        "accession_version": accession_version,
        "url": url,
        "expected_size_bytes": expected_size_bytes,
        "pre_download_eligible": exact_accession and deterministic_url and finite_size,
        "remote_checksum_required_for_pre_download": False,
        "official_checksum_verification_required": checksum_present,
        "steps": steps,
    }


def pre_download_eligibility(row: Mapping[str, Any], *, require_diversity: bool) -> Dict[str, Any]:
    blockers = []
    h1 = str(row.get("h1_accession_version", ""))
    if not re.match(r"^GC[AF]_\d+\.\d+$", h1):
        blockers.append("H1_EXACT_VERSION_MISSING")
    try:
        if int(str(row.get("ncbi_taxid", ""))) <= 0:
            raise ValueError
    except ValueError:
        blockers.append("TAXON_IDENTITY_NOT_PROVEN")
    annotation = str(row.get("annotation_accession_version", ""))
    if not annotation or annotation.startswith("UNRESOLVED"):
        blockers.append("ANNOTATION_EXACT_VERSION_MISSING")
    if row.get("annotation_reference_accession_version") != h1:
        blockers.append("ANNOTATION_H1_REFERENCE_MISMATCH")
    if row.get("annotation_sequence_region_linkage_status") != "proven_official_exact_h1":
        blockers.append("ANNOTATION_SEQUENCE_REGION_LINKAGE_NOT_PROVEN")
    if require_diversity:
        if not re.match(r"^GC[AF]_\d+\.\d+$", str(row.get("h2_accession_version", ""))):
            blockers.append("TIER3A_H2_EXACT_VERSION_MISSING")
        if row.get("same_individual_status") != "yes":
            blockers.append("TIER3A_SAME_INDIVIDUAL_NOT_AFFIRMATIVE")
        if row.get("phase_evidence_status") != "affirmative_correctly_phased":
            blockers.append("TIER3A_PHASE_EVIDENCE_NOT_AFFIRMATIVE")
        if not str(row.get("pair_evidence_url", "")).startswith("https://"):
            blockers.append("TIER3A_PAIR_EVIDENCE_URL_MISSING")
    return {"eligible": not blockers, "blockers": blockers}


def evaluate_post_alignment_measurements(metrics: Mapping[str, Any]) -> Dict[str, Any]:
    failures = []
    thresholds = POST_ALIGNMENT_MEASUREMENT_CONTRACT["minimum_thresholds"]
    for field in ("callable_bases", "callable_fraction", "queryable_gene_count", "queryable_gene_bases"):
        if metrics.get(field) is None:
            failures.append(field + "_missing")
            continue
        try:
            value = float(metrics[field])
        except (TypeError, ValueError):
            failures.append(field + "_invalid")
            continue
        if value < float(thresholds[field]):
            failures.append(field + "_below_minimum")
    return {
        "accepted": not failures,
        "failed_thresholds": failures,
        "result_disposition": "include" if not failures else "exclude",
    }


def _load_tsv(path: Path) -> Tuple[List[str], List[Dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def _write_tsv(path: Path, fieldnames: Sequence[str], rows: Sequence[Mapping[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix=path.name + ".", dir=str(path.parent))
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for row in rows:
                writer.writerow({key: row.get(key, "") for key in fieldnames})
        os.replace(temporary, str(path))
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def _candidate_accessions(rows: Sequence[Mapping[str, str]]) -> Dict[str, str]:
    by_id = {row["candidate_id"]: row for row in rows}
    missing = [candidate_id for candidate_id in PRIORITIZED_CANDIDATE_IDS if candidate_id not in by_id]
    if missing:
        raise RuntimeError("prioritized candidates absent from frozen baseline: " + ",".join(missing))
    return {
        candidate_id: by_id[candidate_id]["manifest_refseq_accession"]
        for candidate_id in PRIORITIZED_CANDIDATE_IDS
    }


def _report_request(accessions: Sequence[str]) -> NormalizedRequest:
    ordered = tuple(sorted(accessions))
    return NormalizedRequest(
        source="NCBI Datasets",
        source_version=NCBI_DATASETS_VERSION,
        method="GET",
        endpoint=(
            "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/"
            + ",".join(ordered)
            + "/dataset_report"
        ),
        parameters={"accessions": ordered, "report_type": "dataset_report"},
    )


def _download_summary_request(accessions: Sequence[str]) -> NormalizedRequest:
    ordered = tuple(sorted(accessions))
    return NormalizedRequest(
        source="NCBI Datasets",
        source_version=NCBI_DATASETS_VERSION,
        method="GET",
        endpoint=(
            "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/"
            + ",".join(ordered)
            + "/download_summary"
        ),
        parameters={
            "accessions": ordered,
            "report_type": "download_summary",
            "_query": {"include_annotation_type": "GENOME_GFF"},
        },
    )


def _sequence_request(accession: str) -> NormalizedRequest:
    return NormalizedRequest(
        source="NCBI Datasets",
        source_version=NCBI_DATASETS_VERSION,
        method="GET",
        endpoint=(
            "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/"
            + accession
            + "/sequence_reports"
        ),
        parameters={
            "accessions": [accession],
            "report_type": "sequence_reports",
            "batch_endpoint_available": False,
            "_query": {"page_size": 10000},
        },
    )


def _ftp_base(accession: str, assembly_name: str) -> str:
    digits = accession.split("_")[1].split(".")[0]
    padded = digits.zfill(9)
    groups = "/".join((padded[0:3], padded[3:6], padded[6:9]))
    safe_name = re.sub(r"[^A-Za-z0-9._-]+", "_", assembly_name)
    return "https://ftp.ncbi.nlm.nih.gov/genomes/all/{0}/{1}/{2}_{3}".format(
        accession[:3], groups, accession, safe_name
    )


def _ftp_request(url: str, *, method: str, accession: str, object_type: str) -> NormalizedRequest:
    return NormalizedRequest(
        source="NCBI genomes FTP",
        source_version=accession,
        method=method,
        endpoint=url,
        parameters={"accession": accession, "object_type": object_type},
    )


def _software_client() -> BatchedCachedClient:
    return BatchedCachedClient(
        cache=ImmutableResponseCache(CACHE_ROOT),
        checkpoint_path=DEFAULT_CHECKPOINT,
        transport=urllib_transport,
        software_environment=environment_identity(),
    )


def snapshot_prior_refusal() -> None:
    PRIOR_ROOT.mkdir(parents=True, exist_ok=True)
    for source in (DEFAULT_MANIFEST, DEFAULT_REJECTIONS, DEFAULT_BUDGET):
        target = PRIOR_ROOT / source.name
        if target.exists():
            continue
        shutil.copyfile(str(source), str(target))


def refresh_cache() -> None:
    snapshot_prior_refusal()
    _columns, rows = _load_tsv(PRIOR_ROOT / "vgp_pilot_manifest.tsv")
    accessions_by_id = _candidate_accessions(rows)
    accessions = list(accessions_by_id.values())
    client = _software_client()
    report_body, _report_entry = client.fetch(_report_request(accessions))
    client.fetch(_download_summary_request(accessions))
    reports = json.loads(report_body.decode("utf-8")).get("reports", [])
    reports_by_accession = {report["accession"]: report for report in reports}
    for accession in sorted(accessions):
        if accession not in reports_by_accession:
            raise RuntimeError("NCBI batch report omitted exact accession " + accession)
        report = reports_by_accession[accession]
        client.fetch(_sequence_request(accession))
        assembly_name = report.get("assembly_info", {}).get("assembly_name", "")
        if not assembly_name:
            raise RuntimeError("NCBI report lacks assembly name for " + accession)
        base = _ftp_base(accession, assembly_name)
        basename = base.rsplit("/", 1)[-1]
        client.fetch(_ftp_request(base + "/md5checksums.txt", method="GET", accession=accession, object_type="md5_catalog"))
        client.fetch(
            _ftp_request(
                base + "/" + basename + "_genomic.fna.gz",
                method="HEAD",
                accession=accession,
                object_type="genomic_fasta",
            )
        )
        client.fetch(
            _ftp_request(
                base + "/" + basename + "_genomic.gff.gz",
                method="HEAD",
                accession=accession,
                object_type="genome_gff",
            )
        )


def _require_cached(request: NormalizedRequest) -> Tuple[bytes, Dict[str, Any]]:
    cached = ImmutableResponseCache(CACHE_ROOT).load(request)
    if cached is None:
        raise RuntimeError("required immutable response is not cached: " + request.cache_key)
    return cached


def _parse_md5(body: bytes) -> Dict[str, str]:
    parsed = {}
    for raw in body.decode("utf-8").splitlines():
        parts = raw.strip().split(None, 1)
        if len(parts) == 2:
            parsed[parts[1].lstrip("./")] = parts[0].lower()
    return parsed


def _content_length(entry: Mapping[str, Any]) -> int:
    headers = entry.get("response_headers", {})
    value = headers.get("content-length")
    if not value or not str(value).isdigit() or int(value) <= 0:
        raise RuntimeError("official HEAD response lacks finite Content-Length")
    return int(value)


def _native_linkage(accession: str, report: Mapping[str, Any], sequence_body: bytes) -> bool:
    annotation = report.get("annotation_info") or {}
    if not annotation.get("name") or not annotation.get("release_date"):
        return False
    sequence = json.loads(sequence_body.decode("utf-8"))
    reports = sequence.get("reports", [])
    total_count = sequence.get("total_count")
    if not reports or total_count != len(reports):
        return False
    return all(item.get("assembly_accession") == accession for item in reports)


def _stringify_obligations(obligations: Sequence[Mapping[str, Any]]) -> str:
    return json.dumps(list(obligations), sort_keys=True, separators=(",", ":"))


def _prepare_resolved_rows() -> Tuple[List[str], List[Dict[str, Any]], Dict[str, Any]]:
    baseline_path = PRIOR_ROOT / "vgp_pilot_manifest.tsv"
    columns, baseline = _load_tsv(baseline_path)
    accessions_by_id = _candidate_accessions(baseline)
    accessions = list(accessions_by_id.values())
    report_body, report_entry = _require_cached(_report_request(accessions))
    download_body, download_entry = _require_cached(_download_summary_request(accessions))
    reports = json.loads(report_body.decode("utf-8")).get("reports", [])
    reports_by_accession = {report["accession"]: report for report in reports}
    download_summary = json.loads(download_body.decode("utf-8"))
    if download_summary.get("record_count") != len(accessions):
        raise RuntimeError("batched download summary does not cover prioritized accessions")

    rank = {candidate_id: index + 1 for index, candidate_id in enumerate(PRIORITIZED_CANDIDATE_IDS)}
    output = []
    metadata_eligible = []
    for baseline_row in baseline:
        row: Dict[str, Any] = dict(baseline_row)
        row.update({column: "" for column in REPAIRED_COLUMNS})
        row["pilot_selected"] = "no"
        row["post_alignment_measurement_contract_id"] = POST_ALIGNMENT_MEASUREMENT_CONTRACT["contract_id"]
        row["post_alignment_result_disposition"] = "not_measured_exclude_until_adequate"
        candidate_id = row["candidate_id"]
        if candidate_id not in rank:
            row["resolution_status"] = "rejected_not_prioritized_small_pilot"
            row["resolution_reason_codes"] = "NOT_IN_PRIORITIZED_CLADE_STRATIFIED_SET"
            row["acceptance_status"] = "rejected_not_prioritized"
            row["explicit_acceptance_or_rejection_reason"] = (
                "outside_six_species_prioritized_clade_stratified_metadata_resolution"
            )
            row["blocking_requirement_ids"] = "PILOT_LIMIT"
            output.append(row)
            continue

        accession = accessions_by_id[candidate_id]
        report = reports_by_accession.get(accession)
        if report is None:
            row["resolution_status"] = "rejected_unresolved"
            row["resolution_reason_codes"] = "OFFICIAL_BATCH_REPORT_MISSING"
            row["acceptance_status"] = "rejected_unresolved"
            output.append(row)
            continue
        sequence_body, sequence_entry = _require_cached(_sequence_request(accession))
        assembly_info = report.get("assembly_info", {})
        assembly_name = assembly_info.get("assembly_name", "")
        base = _ftp_base(accession, assembly_name)
        basename = base.rsplit("/", 1)[-1]
        fasta_url = base + "/" + basename + "_genomic.fna.gz"
        gff_url = base + "/" + basename + "_genomic.gff.gz"
        md5_body, md5_entry = _require_cached(
            _ftp_request(base + "/md5checksums.txt", method="GET", accession=accession, object_type="md5_catalog")
        )
        _empty, fasta_entry = _require_cached(
            _ftp_request(fasta_url, method="HEAD", accession=accession, object_type="genomic_fasta")
        )
        _empty, gff_entry = _require_cached(
            _ftp_request(gff_url, method="HEAD", accession=accession, object_type="genome_gff")
        )
        md5s = _parse_md5(md5_body)
        fasta_size = _content_length(fasta_entry)
        gff_size = _content_length(gff_entry)
        fasta_name = fasta_url.rsplit("/", 1)[-1]
        gff_name = gff_url.rsplit("/", 1)[-1]
        annotation = report.get("annotation_info") or {}
        assembly_stats = report.get("assembly_stats") or {}
        native = _native_linkage(accession, report, sequence_body)
        taxid = int(report.get("organism", {}).get("tax_id", 0))
        catalog_taxid = int(row["ncbi_taxid"])
        row.update(
            {
                "resolution_priority": rank[candidate_id],
                "resolution_status": "metadata_eligible_execution_cap_blocked",
                "resolution_reason_codes": "QUOTA_UNAVAILABLE",
                "metadata_cache_request_keys": ";".join(
                    entry["request_key"]
                    for entry in (
                        report_entry,
                        download_entry,
                        sequence_entry,
                        md5_entry,
                        fasta_entry,
                        gff_entry,
                    )
                ),
                "h1_exact_version_status": "official_current_exact_version",
                "taxon_identity_status": "exact_taxid_match" if taxid == catalog_taxid else "taxid_mismatch",
                "annotation_sequence_region_linkage_status": (
                    "proven_official_exact_h1" if native else "not_proven"
                ),
                "phase_evidence_status": "not_applicable_composition_only",
                "resolved_modality": "tier3c_composition",
                "h1_accession_version": accession,
                "h1_release_version_or_date": assembly_info.get("release_date", ""),
                "h1_assembly_name": assembly_name,
                "h1_haplotype_role": assembly_info.get("assembly_type", "haploid_reference"),
                "h1_fasta_url": fasta_url,
                "h1_provider_md5": md5s.get(fasta_name, ""),
                "h1_fasta_sha256": "ACQUISITION_TIME_LOCAL_SHA256_REQUIRED",
                "h1_sequence_set_sha256": "ACQUISITION_TIME_AFTER_STAGING_REQUIRED",
                "h1_length_bp": assembly_stats.get("total_sequence_length", ""),
                "h1_contig_count": assembly_stats.get("number_of_contigs", ""),
                "h1_contig_n50_bp": assembly_stats.get("contig_n50", ""),
                "h2_accession_version": "",
                "h2_release_version_or_date": "",
                "h2_assembly_name_or_label": "",
                "h2_haplotype_role": "",
                "h2_fasta_url": "",
                "h2_provider_md5": "",
                "h2_fasta_sha256": "",
                "h2_sequence_set_sha256": "",
                "h2_length_bp": "",
                "h2_contig_count": "",
                "h2_contig_n50_bp": "",
                "ncbi_taxid": taxid,
                "ncbi_current_name": report.get("organism", {}).get("organism_name", ""),
                "annotation_accession_version": annotation.get("name", ""),
                "annotation_reference_accession_version": accession,
                "annotation_release_version_or_date": annotation.get("release_date", ""),
                "annotation_provider_and_pipeline": ";".join(
                    str(annotation.get(key, ""))
                    for key in ("provider", "pipeline", "software_version", "method")
                    if annotation.get(key)
                ),
                "annotation_gff_url": gff_url,
                "annotation_provider_md5": md5s.get(gff_name, ""),
                "annotation_gff_sha256": "ACQUISITION_TIME_LOCAL_SHA256_REQUIRED",
                "annotation_native_status": "official_ncbi_native_exact_h1",
                "annotation_contig_audit": "official_complete_sequence_report_exact_h1_scope",
                "annotation_file_status": "official_exact_location_size_known_not_downloaded",
                "callability_reference_accession_version": accession,
                "callable_bases": "POST_ALIGNMENT_REQUIRED",
                "callable_fraction": "POST_ALIGNMENT_REQUIRED",
                "queryable_gene_count": "POST_ALIGNMENT_REQUIRED",
                "queryable_gene_bases": "POST_ALIGNMENT_REQUIRED",
                "resource_retrieved_at_utc": report_entry["retrieved_at_utc"],
                "h1_fasta_compressed_bytes": fasta_size,
                "h1_fasta_uncompressed_bytes": "UNRESOLVED_ACQUISITION_TIME_MEASURE_REQUIRED",
                "h2_fasta_compressed_bytes": "UNRESOLVED_NOT_REQUIRED_COMPOSITION_ONLY",
                "h2_fasta_uncompressed_bytes": "UNRESOLVED_NOT_REQUIRED_COMPOSITION_ONLY",
                "annotation_gff_compressed_bytes": gff_size,
                "annotation_gff_uncompressed_bytes": "UNRESOLVED_ACQUISITION_TIME_MEASURE_REQUIRED",
                "same_individual_status": "not_applicable_composition_only",
                "same_individual_evidence": "not_applicable_composition_only",
                "pair_evidence_url": "",
                "assembly_composition_eligible": "yes",
                "assembly_diversity_eligible": "no",
                "population_genomic_eligible": "no",
                "demographic_eligible": "no",
                "acceptance_status": "blocked_stricter_cap",
                "explicit_acceptance_or_rejection_reason": (
                    "metadata_eligible_but_not_selected_user_quota_unavailable_fail_closed"
                ),
                "blocking_requirement_ids": "QUOTA_UNAVAILABLE",
                "pilot_selection_group": row["lineage_group"],
            }
        )
        obligations = [
            payload_acquisition_obligations(
                accession_version=accession,
                url=fasta_url,
                expected_size_bytes=fasta_size,
                official_checksum_algorithm="md5" if md5s.get(fasta_name) else None,
                official_checksum=md5s.get(fasta_name),
            ),
            payload_acquisition_obligations(
                accession_version=accession,
                url=gff_url,
                expected_size_bytes=gff_size,
                official_checksum_algorithm="md5" if md5s.get(gff_name) else None,
                official_checksum=md5s.get(gff_name),
            ),
        ]
        row["acquisition_obligations"] = _stringify_obligations(obligations)
        eligibility = pre_download_eligibility(row, require_diversity=False)
        if taxid != catalog_taxid:
            eligibility["blockers"].append("TAXON_IDENTITY_MISMATCH")
            eligibility["eligible"] = False
        if not all(item["pre_download_eligible"] for item in obligations):
            eligibility["blockers"].append("PAYLOAD_LOCATION_OR_SIZE_UNRESOLVED")
            eligibility["eligible"] = False
        if not eligibility["eligible"]:
            row["assembly_composition_eligible"] = "no"
            row["resolution_status"] = "rejected_exact_pre_download_requirement"
            row["resolution_reason_codes"] = ";".join(eligibility["blockers"])
            row["acceptance_status"] = "rejected_unresolved"
            row["blocking_requirement_ids"] = row["resolution_reason_codes"]
        else:
            metadata_eligible.append(candidate_id)
        freeze.predicted_resources(row)
        output.append(row)

    catalog_text = DEFAULT_CATALOG.read_text(encoding="utf-8")
    stats = catalog_statistics(catalog_text)
    context = {
        "catalog_statistics": stats,
        "metadata_eligible_candidate_ids": metadata_eligible,
        "prioritized_candidate_ids": list(PRIORITIZED_CANDIDATE_IDS),
    }
    return columns + [column for column in REPAIRED_COLUMNS if column not in columns], output, context


def _rejection_rows(rows: Sequence[Mapping[str, Any]]) -> Tuple[List[str], List[Dict[str, Any]]]:
    fields = [
        "candidate_id",
        "scientific_name",
        "catalog_row_number",
        "resolution_priority",
        "seed_modalities",
        "resolved_modality",
        "h1_accession_version",
        "h2_accession_version",
        "ncbi_taxid",
        "class",
        "order",
        "acceptance_status",
        "explicit_rejection_reason",
        "blocking_requirement_ids",
        "annotation_reference_accession_version",
        "annotation_sequence_region_linkage_status",
        "same_individual_status",
        "phase_evidence_status",
        "metadata_cache_request_keys",
        "resolution_reason_codes",
    ]
    rejected = []
    for row in rows:
        if row.get("pilot_selected") == "yes":
            continue
        rejected.append(
            {
                "candidate_id": row.get("candidate_id", ""),
                "scientific_name": row.get("scientific_name_source", ""),
                "catalog_row_number": row.get("catalog_row_number", ""),
                "resolution_priority": row.get("resolution_priority", ""),
                "seed_modalities": row.get("seed_modalities", ""),
                "resolved_modality": row.get("resolved_modality", ""),
                "h1_accession_version": row.get("h1_accession_version", ""),
                "h2_accession_version": row.get("h2_accession_version", ""),
                "ncbi_taxid": row.get("ncbi_taxid", ""),
                "class": row.get("class", ""),
                "order": row.get("order", ""),
                "acceptance_status": row.get("acceptance_status", ""),
                "explicit_rejection_reason": row.get("explicit_acceptance_or_rejection_reason", ""),
                "resolution_reason_codes": row.get("resolution_reason_codes", ""),
                "blocking_requirement_ids": row.get("blocking_requirement_ids", ""),
                "annotation_reference_accession_version": row.get("annotation_reference_accession_version", ""),
                "annotation_sequence_region_linkage_status": row.get("annotation_sequence_region_linkage_status", ""),
                "same_individual_status": row.get("same_individual_status", ""),
                "phase_evidence_status": row.get("phase_evidence_status", ""),
                "metadata_cache_request_keys": row.get("metadata_cache_request_keys", ""),
            }
        )
    return fields, rejected


def _budget_rows(rows: Sequence[Mapping[str, Any]], quota_status: str) -> Tuple[List[str], List[Dict[str, Any]], Dict[str, Any]]:
    fields = [
        "row_type",
        "candidate_id",
        "scientific_name",
        "pilot_selected",
        "metadata_eligible",
        "download_bytes_exact",
        "compressed_inputs_gib",
        "core_hours_high",
        "peak_memory_gib_high",
        "scratch_gib_high",
        "acquisition_obligations",
        "post_alignment_measurement_contract_id",
        "cap_status",
        "blocking_cap",
    ]
    budget = []
    proposed = []
    selected = []
    for row in rows:
        if not row.get("resolution_priority"):
            continue
        eligible = row.get("assembly_composition_eligible") == "yes"
        item = {
            "row_type": "candidate",
            "candidate_id": row["candidate_id"],
            "scientific_name": row["scientific_name_source"],
            "pilot_selected": row["pilot_selected"],
            "metadata_eligible": "yes" if eligible else "no",
            "download_bytes_exact": row["predicted_download_bytes_exact"],
            "compressed_inputs_gib": round(int(row["predicted_download_bytes_exact"]) / (1024.0 ** 3), 6),
            "core_hours_high": row["predicted_core_hours_high"],
            "peak_memory_gib_high": row["predicted_peak_memory_gib_high"],
            "scratch_gib_high": row["predicted_scratch_gb_high"],
            "acquisition_obligations": row["acquisition_obligations"],
            "post_alignment_measurement_contract_id": POST_ALIGNMENT_MEASUREMENT_CONTRACT["contract_id"],
            "cap_status": "blocked" if quota_status != "available" else "within_numeric_caps",
            "blocking_cap": "QUOTA_UNAVAILABLE" if quota_status != "available" else "",
        }
        budget.append(item)
        if eligible:
            proposed.append(item)
        if row["pilot_selected"] == "yes":
            selected.append(item)

    def aggregate(label: str, subset: Sequence[Mapping[str, Any]]) -> Dict[str, Any]:
        scratch_values = sorted((float(item["scratch_gib_high"]) for item in subset), reverse=True)
        return {
            "row_type": label,
            "candidate_id": label,
            "scientific_name": label,
            "pilot_selected": "n/a",
            "metadata_eligible": "n/a",
            "download_bytes_exact": sum(int(item["download_bytes_exact"]) for item in subset),
            "compressed_inputs_gib": round(sum(float(item["compressed_inputs_gib"]) for item in subset), 6),
            "core_hours_high": round(sum(float(item["core_hours_high"]) for item in subset), 4),
            "peak_memory_gib_high": max((float(item["peak_memory_gib_high"]) for item in subset), default=0.0),
            "scratch_gib_high": round(sum(scratch_values[: BUDGET_CAPS["concurrent_species"]]), 4),
            "acquisition_obligations": "see_candidate_rows",
            "post_alignment_measurement_contract_id": POST_ALIGNMENT_MEASUREMENT_CONTRACT["contract_id"],
            "cap_status": "blocked" if quota_status != "available" else "within_numeric_caps",
            "blocking_cap": "QUOTA_UNAVAILABLE" if quota_status != "available" else "",
        }

    proposed_aggregate = aggregate("aggregate_proposed_metadata_eligible", proposed)
    selected_aggregate = aggregate("aggregate_selected", selected)
    budget.extend((proposed_aggregate, selected_aggregate))
    numeric_pass = (
        len(proposed) <= BUDGET_CAPS["species"]
        and float(proposed_aggregate["compressed_inputs_gib"]) <= BUDGET_CAPS["compressed_inputs_gib"]
        and float(proposed_aggregate["scratch_gib_high"]) <= BUDGET_CAPS["scratch_gib"]
        and float(proposed_aggregate["core_hours_high"]) <= BUDGET_CAPS["core_hours"]
        and float(proposed_aggregate["peak_memory_gib_high"]) <= BUDGET_CAPS["memory_per_job_gib"]
    )
    summary = {
        "caps": BUDGET_CAPS,
        "proposed_species_count": len(proposed),
        "selected_species_count": len(selected),
        "numeric_caps_pass": numeric_pass,
        "stricter_integrated_quota_status": quota_status,
        "overall_cap_vector_pass": numeric_pass and quota_status == "available",
        "scratch_aggregation": "sum_of_two_largest_per_species_high_estimates_for_concurrency_2",
    }
    return fields, budget, summary


def render_outputs(index_path: Path = DEFAULT_INDEX) -> Dict[str, Any]:
    if not (PRIOR_ROOT / "vgp_pilot_manifest.tsv").is_file():
        raise RuntimeError("prior refusal snapshot missing; run --refresh-cache once")
    columns, rows, context = _prepare_resolved_rows()
    root_validation = json.loads(DEFAULT_ROOT_VALIDATION.read_text(encoding="utf-8"))
    quota = root_validation.get("system_evidence", {}).get("quota_state", {})
    quota_status = "available" if quota.get("status") == "available" else "unavailable"
    manifest_rows = [row for row in rows if row.get("resolution_priority")]
    rejection_fields, rejects = _rejection_rows(rows)
    budget_fields, budget, budget_summary = _budget_rows(rows, quota_status)

    # Unknown quota is a stricter cap.  Metadata-eligible rows remain proposed,
    # never selected, until an independent gate records exact quota evidence.
    if not budget_summary["overall_cap_vector_pass"]:
        for row in rows:
            row["pilot_selected"] = "no"
    _write_tsv(DEFAULT_MANIFEST, columns, manifest_rows)
    _write_tsv(DEFAULT_REJECTIONS, rejection_fields, rejects)
    _write_tsv(DEFAULT_BUDGET, budget_fields, budget)

    old_manifest_fields, old_manifest = _load_tsv(PRIOR_ROOT / "vgp_pilot_manifest.tsv")
    _old_rejection_fields, old_rejections = _load_tsv(PRIOR_ROOT / "vgp_pilot_rejections.tsv")
    cache = ImmutableResponseCache(CACHE_ROOT)
    selected = [row for row in rows if row.get("pilot_selected") == "yes"]
    index = {
        "schema_version": 2,
        "task_id": "repair-vgp-candidate",
        "generated_at_utc": utc_now(),
        "resolver_software": environment_identity(),
        "catalog_source": {
            "path": str(DEFAULT_CATALOG),
            "source_commit": CATALOG_COMMIT,
            "source_version": "VGPPhase1-freeze-1.0",
            "sha256": sha256_file(DEFAULT_CATALOG),
        },
        "catalog_statistics": context["catalog_statistics"],
        "lookup_policy": {
            "ncbi_batch_endpoints_used": [
                _report_request(list(_candidate_accessions(old_manifest).values())).endpoint,
                _download_summary_request(list(_candidate_accessions(old_manifest).values())).endpoint,
            ],
            "sequence_report_batch_endpoint_available": False,
            "rate_limit_minimum_interval_seconds": 0.34,
            "retry_after_honored": True,
            "maximum_attempts": 6,
            "bounded_exponential_backoff_seconds": 30.0,
            "jitter": "uniform_bounded",
            "checkpoint_path": str(DEFAULT_CHECKPOINT.relative_to(PROJECT_ROOT)),
            "immutable_cache_key": "sha256(canonical normalized request + source/version)",
        },
        "responses": cache.entries(),
        "prioritized_resolution": {
            "candidate_ids": context["prioritized_candidate_ids"],
            "metadata_eligible_candidate_ids": context["metadata_eligible_candidate_ids"],
            "clade_strategy": "one mammal, bird, reptile, amphibian, actinopterygian, and chondrichthyan",
            "selected_candidate_ids": [row["candidate_id"] for row in selected],
        },
        "selection_comparison": {
            "old_selected_count": sum(row.get("pilot_selected") == "yes" for row in old_manifest),
            "old_rejected_count": len(old_rejections),
            "repaired_selected_count": len(selected),
            "repaired_rejected_count": len(rejects),
            "zero_selection_irreducible_blockers": ["QUOTA_UNAVAILABLE"] if not selected else [],
            "interpretation": (
                "metadata repair resolved six exact native-H1 composition candidates, but the "
                "integrated per-user quota interface remains unavailable; free space does not "
                "override that stricter cap, so zero rows are selected"
            ),
        },
        "resource_budget": budget_summary,
        "post_alignment_measurement_contract": POST_ALIGNMENT_MEASUREMENT_CONTRACT,
        "remaining_obligations": {
            "acquisition": (
                "For any later independently gated row: stage each full exact payload, verify "
                "advertised size, verify official MD5 when present, compute local SHA-256, "
                "reverify it before atomic read-only promotion. No acquisition is authorized now."
            ),
            "run": (
                "Measure callable bases/fraction and queryable protein-coding gene count/bases "
                "after alignment; exclude absent or sub-threshold results. No run is authorized now."
            ),
        },
        "outputs": {
            "manifest": {"path": "analysis/vgp_pilot_manifest.tsv", "sha256": sha256_file(DEFAULT_MANIFEST)},
            "rejections": {"path": "analysis/vgp_pilot_rejections.tsv", "sha256": sha256_file(DEFAULT_REJECTIONS)},
            "size_budget": {"path": "analysis/vgp_pilot_size_budget.tsv", "sha256": sha256_file(DEFAULT_BUDGET)},
            "resolution_report": {
                "path": "analysis/vgp_pilot_resolution_report.md",
                "sha256": sha256_file(PROJECT_ROOT / "analysis/vgp_pilot_resolution_report.md"),
            },
            "prior_refusal_manifest_sha256": sha256_file(PRIOR_ROOT / "vgp_pilot_manifest.tsv"),
            "prior_refusal_rejections_sha256": sha256_file(PRIOR_ROOT / "vgp_pilot_rejections.tsv"),
            "prior_refusal_budget_sha256": sha256_file(PRIOR_ROOT / "vgp_pilot_size_budget.tsv"),
        },
        "authorization_boundary": {
            "biological_payloads_downloaded": 0,
            "jobs_submitted": 0,
            "demographic_inference_performed": False,
            "full_catalog_acquisition_authorized": False,
            "raw_population_bulk_download_authorized": False,
            "decision": "NO_GO_QUOTA_UNAVAILABLE",
        },
    }
    if index["catalog_source"]["sha256"] != CATALOG_SHA256:
        raise RuntimeError("frozen VGP catalog digest changed")
    if index["catalog_statistics"] != {
        "physical_lines": 717,
        "header_lines": 1,
        "data_rows": 716,
        "unique_species": 714,
        "data_row_excess_over_unique_species": 2,
        "duplicated_species": [
            {"scientific_name": "Lophostoma evotis", "multiplicity": 2},
            {"scientific_name": "Micronycteris microtis", "multiplicity": 2},
        ],
    }:
        raise RuntimeError("frozen VGP catalog units or duplicate multiplicities changed")
    _atomic_write(index_path, json.dumps(index, indent=2, sort_keys=True).encode("utf-8") + b"\n")
    return index


def validate_resolution_cache_index(path: Path) -> Dict[str, Any]:
    index = json.loads(Path(path).read_text(encoding="utf-8"))
    if index.get("schema_version") != 2:
        raise RuntimeError("unsupported resolution cache index schema")
    for entry in index.get("responses", []):
        request_data = entry["request"]
        request = NormalizedRequest(
            source=request_data["source"],
            source_version=request_data["source_version"],
            method=request_data["method"],
            endpoint=request_data["endpoint"],
            parameters=request_data["parameters"],
        )
        cached = ImmutableResponseCache(CACHE_ROOT).load(request)
        if cached is None or cached[1]["response_sha256"] != entry["response_sha256"]:
            raise RuntimeError("cache index response missing or invalid: " + entry["request_key"])
        if request.source == "NCBI Datasets":
            observed_version = entry.get("response_headers", {}).get("x-datasets-version")
            if observed_version != request.source_version or observed_version != NCBI_DATASETS_VERSION:
                raise RuntimeError(
                    "NCBI Datasets response version differs from pinned request: "
                    + entry["request_key"]
                )
    for output in ("manifest", "rejections", "size_budget"):
        record = index["outputs"][output]
        if sha256_file(PROJECT_ROOT / record["path"]) != record["sha256"]:
            raise RuntimeError("repaired output digest mismatch: " + output)
    return index


def _measure(args: argparse.Namespace) -> int:
    metrics = json.loads(args.metrics_json.read_text(encoding="utf-8"))
    result = evaluate_post_alignment_measurements(metrics)
    result["contract_id"] = POST_ALIGNMENT_MEASUREMENT_CONTRACT["contract_id"]
    _atomic_write(args.output_json, json.dumps(result, indent=2, sort_keys=True).encode("utf-8") + b"\n")
    return 0 if result["accepted"] else 2


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--refresh-cache", action="store_true", help="perform metadata-only official requests")
    parser.add_argument("--validate-only", action="store_true")
    subparsers = parser.add_subparsers(dest="command")
    measure = subparsers.add_parser("measure", help="evaluate post-alignment denominator metrics")
    measure.add_argument("--metrics-json", type=Path, required=True)
    measure.add_argument("--output-json", type=Path, required=True)
    args = parser.parse_args(argv)
    if args.command == "measure":
        return _measure(args)
    if args.refresh_cache:
        refresh_cache()
    if args.validate_only:
        validate_resolution_cache_index(DEFAULT_INDEX)
        return 0
    render_outputs()
    validate_resolution_cache_index(DEFAULT_INDEX)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
