#!/usr/bin/env python3
"""Cache bounded NCBI metadata for the repaired six-species demography audit.

The default mode is offline and only verifies/reads immutable cache objects.
``--refresh-cache`` performs four small metadata requests (two batched searches
followed by two batched summaries).  It never requests sequence, alignment,
VCF, FASTA, FASTQ, BAM, CRAM, or another biological payload.
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
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "analysis" / "vgp_pilot_manifest.tsv"
CACHE = ROOT / "analysis" / "vgp_demography_cache"
INDEX = CACHE / "index.json"
CHECKPOINT = CACHE / "checkpoint.json"
SOFTWARE = "vgp-demography-metadata-audit/1.0"
USER_AGENT = SOFTWARE + " (metadata-only; contact project maintainers)"
MIN_INTERVAL_SECONDS = 0.34
MAX_ATTEMPTS = 5

PRIORITIZED_IDS = (
    "camelus_dromedarius_gca_036321535_1",
    "colius_striatus_gca_028858725_2",
    "candoia_aspera_gca_035149785_1",
    "dendropsophus_ebraccatus_gca_027789765_1",
    "lepisosteus_oculatus_gca_040954835_1",
    "heterodontus_francisci_gca_036365525_1",
)


def now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def canonical(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode()


def digest(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def file_digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def atomic_write(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=path.name + ".", dir=str(path.parent))
    try:
        with os.fdopen(fd, "wb") as handle:
            handle.write(data)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def load_rows() -> list[dict[str, str]]:
    with MANIFEST.open(newline="", encoding="utf-8") as handle:
        by_id = {row["candidate_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    if set(by_id) != set(PRIORITIZED_IDS):
        raise RuntimeError("repaired manifest does not contain exactly the bounded six candidates")
    return [by_id[candidate_id] for candidate_id in PRIORITIZED_IDS]


def normalized_request(endpoint: str, parameters: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "source": "NCBI Entrez E-utilities",
        "source_version": "live service; response Date header recorded",
        "method": "GET",
        "endpoint": endpoint,
        "parameters": {key: parameters[key] for key in sorted(parameters)},
    }


def request_key(request: Mapping[str, Any]) -> str:
    return digest(canonical(request))


class Cache:
    def paths(self, key: str) -> tuple[Path, Path]:
        return CACHE / "responses" / (key + ".json"), CACHE / "entries" / (key + ".json")

    def load(self, request: Mapping[str, Any]) -> tuple[bytes, dict[str, Any]] | None:
        key = request_key(request)
        response_path, entry_path = self.paths(key)
        if not response_path.is_file() or not entry_path.is_file():
            return None
        body = response_path.read_bytes()
        entry = json.loads(entry_path.read_text(encoding="utf-8"))
        if entry["request_key"] != key or entry["request"] != request:
            raise RuntimeError("immutable cache request mismatch: " + key)
        if entry["response_sha256"] != digest(body):
            raise RuntimeError("immutable cache response digest mismatch: " + key)
        json.loads(body)
        return body, entry

    def store(
        self,
        request: Mapping[str, Any],
        body: bytes,
        headers: Mapping[str, str],
        retrieved_at: str,
    ) -> dict[str, Any]:
        key = request_key(request)
        existing = self.load(request)
        if existing is not None:
            if existing[0] != body:
                raise RuntimeError("immutable cache conflict: " + key)
            return existing[1]
        json.loads(body)
        response_path, entry_path = self.paths(key)
        relative = response_path.relative_to(ROOT).as_posix()
        entry = {
            "request_key": key,
            "request": request,
            "response_path": relative,
            "response_sha256": digest(body),
            "response_size_bytes": len(body),
            "retrieved_at_utc": retrieved_at,
            "response_headers": {str(k).lower(): str(v) for k, v in headers.items()},
            "software_environment": {
                "software": SOFTWARE,
                "python": platform.python_version(),
                "implementation": platform.python_implementation(),
                "platform": platform.platform(),
                "guix_channels_sha256": file_digest(ROOT / "analysis/guix/channels.scm"),
                "guix_manifest_sha256": file_digest(ROOT / "analysis/guix/manifest.scm"),
            },
        }
        atomic_write(response_path, body)
        atomic_write(entry_path, json.dumps(entry, indent=2, sort_keys=True).encode() + b"\n")
        return entry


class Client:
    def __init__(self) -> None:
        self.cache = Cache()
        self.last_request = 0.0
        self.retry_after_events: list[dict[str, Any]] = []

    def fetch(self, endpoint: str, parameters: Mapping[str, Any], refresh: bool) -> tuple[dict[str, Any], dict[str, Any]]:
        request = normalized_request(endpoint, parameters)
        cached = self.cache.load(request)
        if cached is not None:
            return json.loads(cached[0]), cached[1]
        if not refresh:
            raise RuntimeError("missing immutable cache object; rerun once with --refresh-cache: " + request_key(request))
        query = urllib.parse.urlencode(parameters)
        url = endpoint + "?" + query
        for attempt in range(1, MAX_ATTEMPTS + 1):
            delay = MIN_INTERVAL_SECONDS - (time.monotonic() - self.last_request)
            if delay > 0:
                time.sleep(delay)
            try:
                req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT, "Accept": "application/json"})
                with urllib.request.urlopen(req, timeout=60) as response:
                    self.last_request = time.monotonic()
                    body = response.read()
                    entry = self.cache.store(request, body, dict(response.headers), now())
                    return json.loads(body), entry
            except urllib.error.HTTPError as error:
                self.last_request = time.monotonic()
                if error.code not in {429, 500, 502, 503, 504} or attempt == MAX_ATTEMPTS:
                    raise
                header = error.headers.get("Retry-After")
                if header:
                    try:
                        sleep_seconds = float(header)
                    except ValueError:
                        parsed = email.utils.parsedate_to_datetime(header)
                        sleep_seconds = max(0.0, parsed.timestamp() - time.time())
                    self.retry_after_events.append({"attempt": attempt, "value": header, "seconds": sleep_seconds})
                else:
                    sleep_seconds = min(8.0, 0.5 * (2 ** (attempt - 1))) + random.uniform(0, 0.1)
                time.sleep(sleep_seconds)
        raise AssertionError("unreachable")


def write_checkpoint(entries: list[dict[str, Any]], state: str) -> None:
    payload = {
        "schema_version": 1,
        "task_id": "audit-vgp-demography",
        "state": state,
        "completed_request_keys": [entry["request_key"] for entry in entries],
        "pending_request_keys": [],
        "updated_at_utc": now(),
    }
    atomic_write(CHECKPOINT, json.dumps(payload, indent=2, sort_keys=True).encode() + b"\n")


def run(refresh: bool) -> dict[str, Any]:
    rows = load_rows()
    biosamples = [row["biosample_accession"] for row in rows]
    taxids = [row["ncbi_taxid"] for row in rows]
    endpoint_search = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    endpoint_summary = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    client = Client()
    entries: list[dict[str, Any]] = []

    taxonomy_query = " OR ".join(f"txid{value}[Organism:exp]" for value in taxids)
    taxonomy, entry = client.fetch(endpoint_search, {"db": "taxonomy", "term": taxonomy_query, "retmax": "100", "retmode": "json"}, refresh)
    entries.append(entry)
    taxonomy_ids = taxonomy["esearchresult"]["idlist"]
    _, entry = client.fetch(endpoint_summary, {"db": "taxonomy", "id": ",".join(taxonomy_ids), "retmode": "json"}, refresh)
    entries.append(entry)

    sra_query = " OR ".join(f"{value}[BioSample]" for value in biosamples)
    sra, entry = client.fetch(endpoint_search, {"db": "sra", "term": sra_query, "retmax": "10000", "retmode": "json"}, refresh)
    entries.append(entry)
    sra_ids = sra["esearchresult"]["idlist"]
    _, entry = client.fetch(endpoint_summary, {"db": "sra", "id": ",".join(sra_ids), "retmode": "json"}, refresh)
    entries.append(entry)

    write_checkpoint(entries, "complete")
    index = {
        "schema_version": 1,
        "task_id": "audit-vgp-demography",
        "generated_at_utc": now(),
        "manifest": {"path": MANIFEST.relative_to(ROOT).as_posix(), "sha256": file_digest(MANIFEST)},
        "audit_denominator": {
            "interpretation": "six repaired prioritized metadata-eligible candidates; execution remains quota-blocked",
            "candidate_ids": list(PRIORITIZED_IDS),
            "biosample_accessions": biosamples,
            "ncbi_taxids": taxids,
        },
        "lookup_policy": {
            "metadata_only": True,
            "batching": "one OR search per database followed by one comma-ID batch summary",
            "ncbi_requests": 4,
            "rate_limit_minimum_interval_seconds": MIN_INTERVAL_SECONDS,
            "maximum_attempts": MAX_ATTEMPTS,
            "retry_after_aware": True,
            "retry_after_events_observed": client.retry_after_events,
            "resume_checkpoint": CHECKPOINT.relative_to(ROOT).as_posix(),
            "immutable_cache_key": "sha256(canonical normalized request)",
        },
        "authorization_boundary": {
            "biological_payloads_downloaded": 0,
            "raw_population_bulk_download_authorized": False,
            "demographic_inference_performed": False,
            "jobs_submitted": 0,
        },
        "responses": entries,
    }
    atomic_write(INDEX, json.dumps(index, indent=2, sort_keys=True).encode() + b"\n")
    return index


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--refresh-cache", action="store_true")
    args = parser.parse_args()
    result = run(args.refresh_cache)
    print("VGP_DEMOGRAPHY_METADATA_OK responses=%d" % len(result["responses"]))


if __name__ == "__main__":
    main()
