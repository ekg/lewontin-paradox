#!/usr/bin/env python3
"""Build and validate the project-neutral VGP Freeze 1 annotation catalog.

The catalog is intentionally closed over three namespaces:

* source-relative objects in the immutable Freeze 1 mirror;
* content-addressed objects already reachable through an annotation view; and
* exact official NCBI GFF3 files advertised by the current assembly metadata.

Large payloads are never copied into an accession view.  A fetched payload is
written once to the canonical VGP SHA-256 CAS and views are relative symlinks.
The immutable Freeze 1 source objects are only opened for reading.
"""

from __future__ import annotations

import argparse
import collections
import concurrent.futures
import csv
import dataclasses
import datetime as dt
import errno
import gzip
import hashlib
import io
import json
import os
import re
import shutil
import sys
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request
import zlib
from pathlib import Path
from typing import Any, Iterable, Iterator, Mapping, Sequence


SCHEMA_VERSION = "vgp-annotation-catalog-v1"
PARSER_VERSION = 2
ANNOTATION_RE = re.compile(r"\.(gff3?|gtf)(?:\.(gz|bgz|bgzf))?$", re.I)
ACCESSION_RE = re.compile(r"\b(GC[AF]_\d{9}\.\d+)\b")
SEQUENCE_REGION_RE = re.compile(br"^##sequence-region\s+(\S+)\s+(\d+)\s+(\d+)\s*$")
HEADER_ACCESSION_RE = re.compile(br"^#!genome-build-accession\s+\S+:(GC[AF]_\d{9}\.\d+)\s*$")
GFF_ID_RE = re.compile(br"(?:^|;)ID=([^;]+)")
GFF_PARENT_RE = re.compile(br"(?:^|;)Parent=([^;]+)")
GTF_GENE_RE = re.compile(br'(?:^|;\s*)gene_id\s+"([^"]+)"')
GTF_TRANSCRIPT_RE = re.compile(br'(?:^|;\s*)transcript_id\s+"([^"]+)"')
GZIP_SUFFIXES = {"gz", "bgz", "bgzf"}
GZIP_ERROR = getattr(gzip, "BadGzipFile", OSError)
HARD_PARSE_FIELDS = (
    "decompression_errors",
    "column_errors",
    "coordinate_errors",
    "cds_phase_errors",
    "missing_parent_references",
)


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def json_dump(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")
    os.replace(temporary, path)


def write_tsv(path: Path, rows: Iterable[Mapping[str, Any]], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    opener = gzip.open if path.suffix == ".gz" else open
    kwargs: dict[str, Any] = {"mode": "wt", "encoding": "utf-8", "newline": ""}
    with opener(temporary, **kwargs) as handle:  # type: ignore[arg-type]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: scalar(row.get(field, "")) for field in fields})
    os.replace(temporary, path)


def scalar(value: Any) -> Any:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (dict, list, tuple)):
        return json.dumps(value, sort_keys=True, separators=(",", ":"))
    return value


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256_file(path: Path, chunk_size: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        chunk = handle.read(chunk_size)
        while chunk:
            digest.update(chunk)
            chunk = handle.read(chunk_size)
    return digest.hexdigest()


def md5_file(path: Path, chunk_size: int = 8 << 20) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        chunk = handle.read(chunk_size)
        while chunk:
            digest.update(chunk)
            chunk = handle.read(chunk_size)
    return digest.hexdigest()


def detect_format(path_or_name: str) -> tuple[str, str]:
    match = ANNOTATION_RE.search(path_or_name)
    if not match:
        raise ValueError(f"not an annotation suffix: {path_or_name}")
    suffix_format = match.group(1).lower()
    annotation_format = "GFF3" if suffix_format in {"gff", "gff3"} else "GTF"
    suffix_compression = (match.group(2) or "").lower()
    return annotation_format, suffix_compression


def detect_compression(path: Path, suffix_compression: str = "") -> str:
    with path.open("rb") as handle:
        header = handle.read(18)
    if not header.startswith(b"\x1f\x8b"):
        return "NONE"
    # BGZF stores a BC subfield in the gzip extra field.  Eagerly read only the
    # short fixed header: all valid BGZF blocks expose BC within these bytes.
    if len(header) >= 18 and header[3] & 0x04 and header[12:14] == b"BC":
        return "BGZF"
    return "GZIP"


def annotation_basename_class(name: str) -> tuple[str, str]:
    lower = name.lower()
    if ".ncbirefseq." in lower:
        return "NCBI RefSeq/UCSC genePred export", "ncbiRefSeq"
    if ".ncbigene." in lower:
        return "NCBI Gene/UCSC genePred export", "ncbiGene"
    if ".xenorefgene." in lower:
        return "UCSC xenoRefGene", "xenoRefGene"
    if ".augustus." in lower:
        return "AUGUSTUS", "AUGUSTUS"
    if "catliftoffgenes" in lower:
        match = re.search(r"\.catliftoffgenes([^.]*)", name, re.I)
        return "Comparative Annotation Toolkit/LiftOff", match.group(1) if match else "CAT-LiftOff"
    if "veupathgenes" in lower:
        return "VEuPathDB/UCSC genePred export", "VEuPathDB-freeze1"
    if ".ensgene." in lower:
        match = re.search(r"\.ensGene\.([^.]+)\.gtf", name, re.I)
        return "Ensembl/UCSC genePred export", match.group(1) if match else "ensGene"
    return "unknown", "freeze1-snapshot"


def safe_component(value: str) -> str:
    value = value.strip().replace("/", "_")
    return re.sub(r"[^A-Za-z0-9._+-]+", "_", value) or "unknown"


def cas_path(vgp_root: Path, digest: str) -> Path:
    return vgp_root / "objects" / "sha256" / digest[:2] / digest[2:4] / digest


def relative_symlink(target: Path, link: Path) -> None:
    link.parent.mkdir(parents=True, exist_ok=True)
    expected = os.path.relpath(target, link.parent)
    if link.is_symlink():
        if os.readlink(link) == expected and link.resolve() == target.resolve():
            return
        raise RuntimeError(f"refusing to replace divergent symlink {link}")
    if link.exists():
        if link.samefile(target):
            return
        raise RuntimeError(f"refusing to replace existing view object {link}")
    temporary = link.with_name(f".{link.name}.{os.getpid()}.tmp")
    os.symlink(expected, temporary)
    os.replace(temporary, link)


def request_json(url: str, retries: int = 6) -> dict[str, Any]:
    for attempt in range(retries):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": "VGP-annotation-catalog/1.0"})
            with urllib.request.urlopen(request, timeout=180) as response:
                return json.load(response)
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError):
            if attempt + 1 == retries:
                raise
            time.sleep(min(30, 2**attempt))
    raise AssertionError("unreachable")


def fetch_ncbi_metadata(accessions: Sequence[str], cache_path: Path, refresh: bool) -> dict[str, dict[str, Any]]:
    if cache_path.exists() and not refresh:
        cached = json.loads(cache_path.read_text(encoding="utf-8"))
        reports = cached.get("reports_by_accession", {})
        if set(accessions).issubset(reports) or cached.get("closed_world_accessions") == list(accessions):
            return reports

    reports: dict[str, dict[str, Any]] = {}
    failures: dict[str, str] = {}
    # The endpoint defaults to a 20-report page.  Keep batches at 20 so no
    # requested accession is silently truncated.
    for start in range(0, len(accessions), 20):
        batch = list(accessions[start : start + 20])
        url = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/" + ",".join(batch) + "/dataset_report"
        payload = request_json(url)
        for report in payload.get("reports", []):
            reports[report["accession"]] = report
        time.sleep(0.2)

    # Some suppressed or recently versioned records are omitted from a batch
    # response.  Resolve them individually and retain a closed-world failure
    # reason rather than pretending that they were never requested.
    for accession in sorted(set(accessions) - set(reports)):
        url = f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{accession}/dataset_report"
        try:
            payload = request_json(url)
            found = payload.get("reports", [])
            if found:
                reports[found[0]["accession"]] = found[0]
            else:
                failures[accession] = "EMPTY_DATASET_REPORT"
        except urllib.error.HTTPError as error:
            failures[accession] = f"HTTP_{error.code}"
        except Exception as error:  # persisted and surfaced in the final summary
            failures[accession] = f"{type(error).__name__}:{error}"
        time.sleep(0.2)

    json_dump(
        cache_path,
        {
            "schema_version": "ncbi-datasets-v2-assembly-snapshot-v1",
            "retrieved_at_utc": utc_now(),
            "endpoint": "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{accessions}/dataset_report",
            "closed_world_accessions": list(accessions),
            "reports_by_accession": reports,
            "unresolved_accessions": failures,
        },
    )
    return reports


def official_gff_url(report: Mapping[str, Any]) -> tuple[str, str, str]:
    accession = str(report["accession"])
    assembly_name = str(report["assembly_info"]["assembly_name"]).replace(" ", "_").replace("/", "_")
    stem = f"{accession}_{assembly_name}"
    route = f"{accession[:3]}/{accession[4:7]}/{accession[7:10]}/{accession[10:13]}/{stem}"
    base = f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{route}"
    return f"{base}/{stem}_genomic.gff.gz", f"{base}/md5checksums.txt", f"{stem}_genomic.gff.gz"


def discover_official_gff_url(report: Mapping[str, Any]) -> tuple[str, str, str]:
    """Resolve exceptional NCBI FTP stems from the authoritative directory.

    Most FTP stems are a mechanical accession + assembly-name join.  A few
    historical records sanitize '+' or retain an older assembly-name spelling,
    so the dataset-report spelling is insufficient.  Directory discovery is a
    narrow fallback and still selects only the exact accession version.
    """

    accession = str(report["accession"])
    parent = (
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/"
        f"{accession[:3]}/{accession[4:7]}/{accession[7:10]}/{accession[10:13]}/"
    )
    listing = fetch_text(parent)
    matches = sorted(
        set(
            re.findall(
                r'href="(' + re.escape(accession) + r'_[^"/]+/?)"',
                listing,
            )
        )
    )
    if len(matches) != 1:
        raise RuntimeError(f"expected one exact FTP directory for {accession}, observed {matches}")
    stem = matches[0].rstrip("/")
    base = parent + urllib.parse.quote(stem, safe="._-")
    basename = f"{stem}_genomic.gff.gz"
    return f"{base}/{basename}", f"{base}/md5checksums.txt", basename


def fetch_text(url: str, retries: int = 6) -> str:
    for attempt in range(retries):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": "VGP-annotation-catalog/1.0"})
            with urllib.request.urlopen(request, timeout=180) as response:
                return response.read().decode("utf-8")
        except (urllib.error.URLError, TimeoutError):
            if attempt + 1 == retries:
                raise
            time.sleep(min(30, 2**attempt))
    raise AssertionError("unreachable")


def expected_md5(md5_text: str, basename: str) -> str:
    for line in md5_text.splitlines():
        pieces = line.split(None, 1)
        if len(pieces) == 2 and pieces[1].lstrip("*./") == basename:
            return pieces[0].lower()
    raise RuntimeError(f"{basename} absent from official md5checksums.txt")


def download_one_official(
    report: Mapping[str, Any], vgp_root: Path, annotation_root: Path, refresh: bool
) -> dict[str, Any]:
    accession = str(report["accession"])
    info = dict(report.get("annotation_info", {}))
    annotation_name = str(info.get("name", f"{accession}-annotation"))
    source_url, checksum_url, basename = official_gff_url(report)
    provenance_dir = annotation_root / "provenance"
    provenance_dir.mkdir(parents=True, exist_ok=True)

    existing: dict[str, Any] | None = None
    for provenance_path in provenance_dir.glob("*.json"):
        try:
            candidate = json.loads(provenance_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        if candidate.get("source_url") == source_url and candidate.get("state") == "VERIFIED":
            physical = Path(candidate["cas_path"])
            if physical.exists() and physical.stat().st_size == candidate.get("size_bytes"):
                existing = candidate
                break
    if existing and not refresh:
        view = annotation_root / "by-accession" / accession / safe_component(annotation_name) / basename
        relative_symlink(Path(existing["cas_path"]), view)
        existing = dict(existing)
        existing["view_path"] = str(view)
        return existing

    try:
        md5_text = fetch_text(checksum_url)
    except urllib.error.HTTPError as error:
        if error.code != 404:
            raise
        source_url, checksum_url, basename = discover_official_gff_url(report)
        md5_text = fetch_text(checksum_url)
    advertised_md5 = expected_md5(md5_text, basename)
    staging_dir = annotation_root / "staging"
    staging_dir.mkdir(parents=True, exist_ok=True)
    fd, temporary_name = tempfile.mkstemp(prefix=f".{accession}.", suffix=".part", dir=staging_dir)
    os.close(fd)
    temporary = Path(temporary_name)
    sha256 = hashlib.sha256()
    md5 = hashlib.md5()
    size = 0
    started = utc_now()
    try:
        request = urllib.request.Request(source_url, headers={"User-Agent": "VGP-annotation-catalog/1.0"})
        with urllib.request.urlopen(request, timeout=600) as response, temporary.open("wb") as output:
            chunk = response.read(8 << 20)
            while chunk:
                output.write(chunk)
                sha256.update(chunk)
                md5.update(chunk)
                size += len(chunk)
                chunk = response.read(8 << 20)
            output.flush()
            os.fsync(output.fileno())
        if md5.hexdigest() != advertised_md5:
            raise RuntimeError(
                f"official MD5 mismatch for {source_url}: {md5.hexdigest()} != {advertised_md5}"
            )
        digest = sha256.hexdigest()
        destination = cas_path(vgp_root, digest)
        destination.parent.mkdir(parents=True, exist_ok=True)
        if destination.exists():
            if destination.stat().st_size != size or sha256_file(destination) != digest:
                raise RuntimeError(f"divergent object already occupies CAS path {destination}")
            temporary.unlink()
        else:
            os.replace(temporary, destination)
            # Persist the containing directory entry before publishing views.
            directory_fd = os.open(destination.parent, os.O_RDONLY)
            try:
                os.fsync(directory_fd)
            finally:
                os.close(directory_fd)
        view = annotation_root / "by-accession" / accession / safe_component(annotation_name) / basename
        relative_symlink(destination, view)
        provenance = {
            "schema_version": "vgp-official-annotation-provenance-v1",
            "state": "VERIFIED",
            "assembly_accession_version": accession,
            "annotation_accession_version": annotation_name,
            "species": report.get("organism", {}).get("organism_name", ""),
            "annotation_provider": info.get("provider", ""),
            "annotation_pipeline": info.get("pipeline", ""),
            "annotation_software_version": info.get("software_version", ""),
            "annotation_release_date": info.get("release_date", ""),
            "source_url": source_url,
            "checksum_url": checksum_url,
            "advertised_md5": advertised_md5,
            "observed_md5": md5.hexdigest(),
            "sha256": digest,
            "size_bytes": size,
            "cas_path": str(destination),
            "view_path": str(view),
            "retrieval_started_utc": started,
            "retrieval_completed_utc": utc_now(),
        }
        json_dump(provenance_dir / f"{digest}.json", provenance)
        return provenance
    except Exception:
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass
        raise


def fetch_official_annotations(
    reports: Mapping[str, Mapping[str, Any]], vgp_root: Path, annotation_root: Path, workers: int, refresh: bool
) -> tuple[list[dict[str, Any]], dict[str, str]]:
    advertised = [report for report in reports.values() if report.get("annotation_info")]
    fetched: list[dict[str, Any]] = []
    failures: dict[str, str] = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, workers)) as executor:
        future_to_accession = {
            executor.submit(download_one_official, report, vgp_root, annotation_root, refresh): str(report["accession"])
            for report in advertised
        }
        for future in concurrent.futures.as_completed(future_to_accession):
            accession = future_to_accession[future]
            try:
                fetched.append(future.result())
            except Exception as error:
                failures[accession] = f"{type(error).__name__}:{error}"
                print(f"fetch failed {accession}: {error}", file=sys.stderr, flush=True)
    fetched.sort(key=lambda row: (row["assembly_accession_version"], row["source_url"]))
    return fetched, failures


class HashingRawReader(io.RawIOBase):
    def __init__(self, raw: io.BufferedReader):
        self.raw = raw
        self.sha256 = hashlib.sha256()
        self.bytes_read = 0

    def readable(self) -> bool:
        return True

    def readinto(self, buffer: bytearray | memoryview) -> int:
        chunk = self.raw.read(len(buffer))
        if not chunk:
            return 0
        count = len(chunk)
        buffer[:count] = chunk
        self.sha256.update(chunk)
        self.bytes_read += count
        return count

    def close(self) -> None:
        try:
            self.raw.close()
        finally:
            super().close()


def open_annotation_stream(path: Path, compression: str) -> tuple[Iterator[bytes], HashingRawReader, Any]:
    raw_file = path.open("rb")
    hashing = HashingRawReader(raw_file)
    buffered = io.BufferedReader(hashing, buffer_size=1 << 20)
    if compression in {"GZIP", "BGZF"}:
        stream: Any = gzip.GzipFile(fileobj=buffered, mode="rb")
    else:
        stream = buffered
    return iter(stream), hashing, stream


def parse_annotation_task(task: Mapping[str, Any]) -> dict[str, Any]:
    path = Path(task["physical_path"])
    annotation_format = str(task["format"])
    compression = str(task["compression"])
    feature_counts: collections.Counter[str] = collections.Counter()
    source_counts: collections.Counter[str] = collections.Counter()
    observed: dict[str, int] = {}
    declared: dict[str, int] = {}
    gff_ids: set[bytes] = set()
    gff_parents: set[bytes] = set()
    embedded_accession = ""
    feature_rows = 0
    comment_rows = 0
    blank_rows = 0
    fasta_rows = 0
    column_errors = 0
    coordinate_errors = 0
    cds_phase_errors = 0
    gtf_parent_attribute_errors = 0
    duplicate_declared_sequence_regions = 0
    circular_sequences: set[str] = set()
    utf8_errors = 0
    decompression_errors = 0
    uncompressed_bytes = 0
    hashing: HashingRawReader | None = None
    stream: Any = None
    try:
        lines, hashing, stream = open_annotation_stream(path, compression)
        in_fasta = False
        for raw_line in lines:
            uncompressed_bytes += len(raw_line)
            line = raw_line.rstrip(b"\r\n")
            if in_fasta:
                fasta_rows += 1
                continue
            if line == b"##FASTA":
                in_fasta = True
                comment_rows += 1
                continue
            if not line:
                blank_rows += 1
                continue
            if line.startswith(b"#"):
                comment_rows += 1
                region = SEQUENCE_REGION_RE.match(line)
                if region:
                    name = region.group(1).decode("utf-8", "replace")
                    if "\ufffd" in name:
                        utf8_errors += 1
                    start, end = int(region.group(2)), int(region.group(3))
                    if name in declared and declared[name] != end:
                        duplicate_declared_sequence_regions += 1
                    declared[name] = end
                    if start != 1 or end < start:
                        coordinate_errors += 1
                accession_match = HEADER_ACCESSION_RE.match(line)
                if accession_match:
                    embedded_accession = accession_match.group(1).decode("ascii")
                continue
            columns = line.split(b"\t")
            if len(columns) != 9:
                column_errors += 1
                continue
            feature_rows += 1
            try:
                seqid = columns[0].decode("utf-8")
                source = columns[1].decode("utf-8")
                feature = columns[2].decode("utf-8")
            except UnicodeDecodeError:
                utf8_errors += 1
                seqid = columns[0].decode("utf-8", "replace")
                source = columns[1].decode("utf-8", "replace")
                feature = columns[2].decode("utf-8", "replace")
            feature_counts[feature] += 1
            source_counts[source] += 1
            try:
                start, end = int(columns[3]), int(columns[4])
            except ValueError:
                coordinate_errors += 1
                continue
            if start < 1 or end < start:
                coordinate_errors += 1
            observed[seqid] = max(observed.get(seqid, 0), end)
            if feature == "CDS" and columns[7] not in {b"0", b"1", b"2"}:
                cds_phase_errors += 1
            attributes = columns[8]
            if feature.lower() == "region" and (
                b"Is_circular=true" in attributes
                or b"is_circular=true" in attributes
                or b"circular=true" in attributes
            ):
                circular_sequences.add(seqid)
            if annotation_format == "GFF3":
                id_match = GFF_ID_RE.search(attributes)
                if id_match:
                    gff_ids.add(id_match.group(1))
                parent_match = GFF_PARENT_RE.search(attributes)
                if parent_match:
                    gff_parents.update(parent_match.group(1).split(b","))
            else:
                gene_match = GTF_GENE_RE.search(attributes)
                transcript_match = GTF_TRANSCRIPT_RE.search(attributes)
                lower_feature = feature.lower()
                if lower_feature in {"cds", "exon", "start_codon", "stop_codon", "utr", "five_prime_utr", "three_prime_utr"}:
                    if not gene_match or not transcript_match:
                        gtf_parent_attribute_errors += 1
                elif lower_feature in {"transcript", "mrna"} and not gene_match:
                    gtf_parent_attribute_errors += 1
        stream.close()
    except (OSError, EOFError, GZIP_ERROR, zlib.error) as error:
        decompression_errors += 1
        error_message = f"{type(error).__name__}:{error}"
    else:
        error_message = ""
    finally:
        if stream is not None:
            try:
                stream.close()
            except Exception:
                pass
    missing_parents = len(gff_parents - gff_ids) if annotation_format == "GFF3" else 0
    result = dict(task)
    result.update(
        {
            "actual_sha256": hashing.sha256.hexdigest() if hashing else "",
            "compressed_bytes_read": hashing.bytes_read if hashing else 0,
            "uncompressed_bytes": uncompressed_bytes,
            "feature_rows": feature_rows,
            "comment_rows": comment_rows,
            "blank_rows": blank_rows,
            "embedded_fasta_rows": fasta_rows,
            "feature_counts": dict(sorted(feature_counts.items())),
            "feature_source_counts": dict(sorted(source_counts.items())),
            "observed_sequence_max_end": observed,
            "declared_sequence_dictionary": declared,
            "circular_sequences": sorted(circular_sequences),
            "observed_sequence_count": len(observed),
            "embedded_sequence_count": len(declared),
            "embedded_assembly_accession_version": embedded_accession,
            "column_errors": column_errors,
            "coordinate_errors": coordinate_errors,
            "cds_phase_errors": cds_phase_errors,
            "gtf_parent_attribute_errors": gtf_parent_attribute_errors,
            "parent_reference_count": len(gff_parents),
            "defined_id_count": len(gff_ids),
            "missing_parent_references": missing_parents,
            "duplicate_declared_sequence_regions": duplicate_declared_sequence_regions,
            "utf8_errors": utf8_errors,
            "decompression_errors": decompression_errors,
            "parse_error_message": error_message,
        }
    )
    return result


def load_fai(path: Path) -> dict[str, int]:
    dictionary: dict[str, int] = {}
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            columns = line.rstrip("\n").split("\t")
            if len(columns) < 2:
                raise RuntimeError(f"malformed FAI {path}:{line_number}")
            name, length_text = columns[:2]
            length = int(length_text)
            if name in dictionary:
                raise RuntimeError(f"duplicate sequence name {name!r} in {path}")
            dictionary[name] = length
    return dictionary


def parse_assembly_report(path: Path, fai: Mapping[str, int]) -> tuple[dict[str, tuple[str, int]], list[str]]:
    aliases: dict[str, set[tuple[str, int]]] = collections.defaultdict(set)
    errors: list[str] = []
    if not path.exists():
        return {}, ["ASSEMBLY_REPORT_MISSING"]
    header: list[str] | None = None
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("# Sequence-Name"):
                header = line[2:].rstrip("\n").split("\t")
                continue
            if line.startswith("#") or not line.strip():
                continue
            columns = line.rstrip("\n").split("\t")
            if header is None:
                if len(columns) >= 10:
                    header = [
                        "Sequence-Name",
                        "Sequence-Role",
                        "Assigned-Molecule",
                        "Assigned-Molecule-Location/Type",
                        "GenBank-Accn",
                        "Relationship",
                        "RefSeq-Accn",
                        "Assembly-Unit",
                        "Sequence-Length",
                        "UCSC-style-name",
                    ]
                else:
                    errors.append("ASSEMBLY_REPORT_HEADER_MISSING")
                    break
            record = dict(zip(header, columns))
            try:
                length = int(record["Sequence-Length"])
            except (KeyError, ValueError):
                errors.append("ASSEMBLY_REPORT_LENGTH_INVALID")
                continue
            names = {
                record.get("Sequence-Name", ""),
                record.get("GenBank-Accn", ""),
                record.get("RefSeq-Accn", ""),
                record.get("UCSC-style-name", ""),
            } - {"", "na", "NA"}
            canonical_candidates = [(name, fai[name]) for name in names if name in fai and fai[name] == length]
            if len(canonical_candidates) != 1:
                continue
            canonical, canonical_length = canonical_candidates[0]
            for name in names:
                aliases[name].add((canonical, canonical_length))
    return {name: next(iter(values)) for name, values in aliases.items() if len(values) == 1}, errors


def bind_annotation(
    annotation: Mapping[str, Any], assembly: Mapping[str, Any], fai: Mapping[str, int], aliases: Mapping[str, tuple[str, int]]
) -> dict[str, Any]:
    accession = str(assembly["accession_version"])
    embedded_accession = str(annotation.get("embedded_assembly_accession_version", ""))
    result: dict[str, Any] = {
        "target_assembly_accession_version": accession,
        "assembly_bgzf_path": assembly["derived_path"],
        "assembly_fai_path": assembly["fai_path"],
        "assembly_bgzf_sha256": assembly["derived_bgzf_sha256"],
        "assembly_fai_sha256": assembly["fai_sha256"],
        "assembly_sequence_count": len(fai),
        "binding_status": "",
        "binding_reason": "",
        "exact_name_count": 0,
        "validated_alias_count": 0,
        "unresolved_sequence_count": 0,
        "length_mismatch_count": 0,
        "coordinate_overrun_count": 0,
        "circular_coordinate_wrap_count": 0,
        "alias_map": {},
    }
    if embedded_accession and embedded_accession != accession:
        result["binding_status"] = "ACCESSION_MISMATCH"
        result["binding_reason"] = f"embedded accession {embedded_accession} != target {accession}"
        return result
    observed = dict(annotation.get("observed_sequence_max_end", {}))
    declared = dict(annotation.get("declared_sequence_dictionary", {}))
    names = set(observed) | set(declared)
    if not names:
        result["binding_status"] = "DICTIONARY_MISMATCH"
        result["binding_reason"] = "annotation has no feature or embedded dictionary sequence names"
        return result
    alias_map: dict[str, str] = {}
    unresolved: list[str] = []
    length_mismatches: list[str] = []
    overruns: list[str] = []
    circular_wraps: list[str] = []
    circular_sequences = set(annotation.get("circular_sequences", []))
    exact = 0
    alias_count = 0
    for name in sorted(names):
        if name in fai:
            canonical, length = name, fai[name]
            exact += 1
        elif name in aliases:
            canonical, length = aliases[name]
            # Equal length is explicit in the assembly report.  If the GFF
            # embeds a sequence-region length, demand equality a second time.
            if name in declared and int(declared[name]) != length:
                length_mismatches.append(name)
                continue
            alias_map[name] = canonical
            alias_count += 1
        else:
            unresolved.append(name)
            continue
        if name in declared and int(declared[name]) != length:
            length_mismatches.append(name)
        observed_end = int(observed.get(name, 0))
        if observed_end > length:
            if name in circular_sequences and observed_end <= 2 * length:
                circular_wraps.append(name)
            else:
                overruns.append(name)
    result.update(
        {
            "exact_name_count": exact,
            "validated_alias_count": alias_count,
            "unresolved_sequence_count": len(unresolved),
            "length_mismatch_count": len(length_mismatches),
            "coordinate_overrun_count": len(overruns),
            "circular_coordinate_wrap_count": len(circular_wraps),
            "alias_map": alias_map,
        }
    )
    if unresolved or length_mismatches or overruns:
        result["binding_status"] = "DICTIONARY_MISMATCH"
        reasons = []
        if unresolved:
            reasons.append(f"{len(unresolved)} unresolved names")
        if length_mismatches:
            reasons.append(f"{len(length_mismatches)} length mismatches")
        if overruns:
            reasons.append(f"{len(overruns)} coordinate overruns")
        result["binding_reason"] = "; ".join(reasons)
    elif alias_count:
        result["binding_status"] = "VALIDATED_ALIAS"
        result["binding_reason"] = f"{alias_count} one-to-one equal-length assembly-report aliases"
    else:
        result["binding_status"] = "EXACT_DICTIONARY"
        result["binding_reason"] = f"all {exact} annotation sequence names and lengths match the target FAI"
    return result


def source_annotation_rows(
    mirror_rows: Sequence[Mapping[str, str]],
    reports: Mapping[str, Mapping[str, Any]],
    freeze_commit: str,
    fallback_species: Mapping[str, str],
) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []
    for mirror in mirror_rows:
        relative = mirror.get("source_relative_path", "")
        if mirror.get("object_type") != "file" or not ANNOTATION_RE.search(relative):
            continue
        accession = mirror["accession_version"]
        physical = Path(mirror["durable_path"])
        annotation_format, suffix_compression = detect_format(relative)
        provider, filename_version = annotation_basename_class(physical.name)
        info = reports.get(accession, {}).get("annotation_info", {})
        version = filename_version
        if filename_version in {"ncbiRefSeq", "ncbiGene"} and info.get("name"):
            version = str(info["name"])
        elif filename_version == "freeze1-snapshot":
            version = freeze_commit
        output.append(
            {
                "object_id": mirror.get("local_sha256", ""),
                "origin": "MIRRORED_FREEZE1",
                "physical_path": str(physical),
                "source_relative_path": relative,
                "source_url": mirror.get("source_endpoint", "")
                or f"rsync://hgdownload.soe.ucsc.edu/hubs/{relative}",
                "size_bytes": int(mirror.get("observed_bytes") or mirror.get("size_bytes") or 0),
                "expected_sha256": mirror.get("local_sha256", ""),
                "upstream_checksum_algorithm": mirror.get("upstream_checksum_algorithm", ""),
                "upstream_checksum": mirror.get("upstream_checksum", ""),
                "cas_path": mirror.get("cas_path", ""),
                "view_paths": [],
                "assembly_accession_version": accession,
                "species": reports.get(accession, {}).get("organism", {}).get("organism_name", "")
                or fallback_species.get(accession, ""),
                "annotation_accession_version": version,
                "annotation_provider": provider,
                "annotation_pipeline": filename_version,
                "annotation_release_date": "",
                "format": annotation_format,
                "compression": detect_compression(physical, suffix_compression),
            }
        )
    return output


def fetched_annotation_rows(fetched: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []
    for provenance in fetched:
        physical = Path(provenance["cas_path"])
        output.append(
            {
                "object_id": provenance["sha256"],
                "origin": "FETCHED_NCBI_OFFICIAL",
                "physical_path": str(physical),
                "source_relative_path": "",
                "source_url": provenance["source_url"],
                "size_bytes": int(provenance["size_bytes"]),
                "expected_sha256": provenance["sha256"],
                "upstream_checksum_algorithm": "md5",
                "upstream_checksum": provenance["advertised_md5"],
                "cas_path": str(physical),
                "view_paths": [provenance["view_path"]],
                "assembly_accession_version": provenance["assembly_accession_version"],
                "species": provenance.get("species", ""),
                "annotation_accession_version": provenance.get("annotation_accession_version", ""),
                "annotation_provider": provenance.get("annotation_provider", ""),
                "annotation_pipeline": provenance.get("annotation_pipeline", ""),
                "annotation_release_date": provenance.get("annotation_release_date", ""),
                "format": "GFF3",
                "compression": detect_compression(physical, "gz"),
            }
        )
    return output


def attach_view_paths(rows: list[dict[str, Any]], roots: Sequence[Path]) -> list[dict[str, Any]]:
    by_physical: dict[Path, list[dict[str, Any]]] = collections.defaultdict(list)
    for row in rows:
        try:
            by_physical[Path(row["physical_path"]).resolve()].append(row)
        except OSError:
            pass
    for root in roots:
        if not root.exists():
            continue
        for directory, _, filenames in os.walk(root, followlinks=False):
            for filename in filenames:
                path = Path(directory) / filename
                if not (path.is_symlink() and ANNOTATION_RE.search(filename)):
                    continue
                try:
                    target = path.resolve(strict=True)
                except OSError:
                    continue
                for row in by_physical.get(target, []):
                    views = set(row.get("view_paths", []))
                    views.add(str(path))
                    row["view_paths"] = sorted(views)
    return rows


def consolidate_logical_objects(rows: Sequence[dict[str, Any]]) -> list[dict[str, Any]]:
    # The same payload can legitimately be present in multiple source-relative
    # paths.  Preserve a row for every physical source object; add an explicit
    # payload multiplicity rather than collapsing the closed-world inventory.
    counts = collections.Counter(row["object_id"] for row in rows)
    for row in rows:
        row["payload_path_multiplicity"] = counts[row["object_id"]]
    return list(rows)


def parser_cache_key(row: Mapping[str, Any]) -> str:
    token = f"parser-v{PARSER_VERSION}\0{row['expected_sha256']}\0{row['physical_path']}\0{row['size_bytes']}"
    return hashlib.sha256(token.encode()).hexdigest()


def parse_all_annotations(
    rows: Sequence[dict[str, Any]], cache_dir: Path, workers: int, refresh: bool
) -> list[dict[str, Any]]:
    cache_dir.mkdir(parents=True, exist_ok=True)
    results: dict[str, dict[str, Any]] = {}
    tasks: list[dict[str, Any]] = []
    for row in rows:
        key = parser_cache_key(row)
        cache = cache_dir / f"{key}.json"
        if cache.exists() and not refresh:
            try:
                cached = json.loads(cache.read_text(encoding="utf-8"))
                if cached.get("expected_sha256") == row.get("expected_sha256"):
                    # Parser caches contain expensive content-derived fields,
                    # but catalog metadata (species, URLs, view paths) can be
                    # enriched on a later run.  Always overlay the current
                    # inventory row so a cache never freezes stale metadata.
                    cached.update(row)
                    results[key] = cached
                    continue
            except (OSError, json.JSONDecodeError):
                pass
        task = dict(row)
        task["parser_cache_key"] = key
        tasks.append(task)
    if tasks:
        with concurrent.futures.ProcessPoolExecutor(max_workers=max(1, workers)) as executor:
            future_to_task = {executor.submit(parse_annotation_task, task): task for task in tasks}
            completed = 0
            for future in concurrent.futures.as_completed(future_to_task):
                task = future_to_task[future]
                result = future.result()
                key = task["parser_cache_key"]
                json_dump(cache_dir / f"{key}.json", result)
                results[key] = result
                completed += 1
                if completed % 25 == 0 or completed == len(tasks):
                    print(f"parsed {completed}/{len(tasks)} uncached annotations", flush=True)
    return [results[parser_cache_key(row)] for row in rows]


def find_assembly_report(mirror_rows: Sequence[Mapping[str, str]], accession: str) -> Path | None:
    candidates = [
        Path(row["durable_path"])
        for row in mirror_rows
        if row.get("accession_version") == accession
        and row.get("object_type") == "file"
        and row.get("source_relative_path", "").endswith("_assembly_report.txt")
        and Path(row.get("durable_path", "")).exists()
    ]
    return candidates[0] if candidates else None


def fallback_species_by_accession(
    mirror_rows: Sequence[Mapping[str, str]], accessions: Sequence[str]
) -> dict[str, str]:
    """Recover scientific names from frozen hub/report metadata.

    NCBI Datasets omits a small set of withdrawn-but-frozen accessions from
    dataset_report responses.  The Freeze 1 release still carries its own
    immutable UCSC hub and NCBI assembly-report metadata, so those names remain
    recoverable without substituting a newer assembly.
    """

    output: dict[str, str] = {}
    rows_by_accession: dict[str, list[Mapping[str, str]]] = collections.defaultdict(list)
    for row in mirror_rows:
        rows_by_accession[row.get("accession_version", "")].append(row)
    for accession in accessions:
        rows = rows_by_accession.get(accession, [])
        hub_candidates = [
            Path(row["durable_path"])
            for row in rows
            if row.get("source_relative_path", "").endswith("/hub.txt")
            and Path(row.get("durable_path", "")).exists()
        ]
        for hub in hub_candidates:
            with hub.open(encoding="utf-8", errors="replace") as handle:
                for line in handle:
                    if line.startswith("scientificName "):
                        output[accession] = line.split(None, 1)[1].strip()
                        break
            if accession in output:
                break
        if accession in output:
            continue
        report = find_assembly_report(mirror_rows, accession)
        if report:
            with report.open(encoding="utf-8", errors="replace") as handle:
                for line in handle:
                    if line.startswith("# Organism name:"):
                        value = line.split(":", 1)[1].strip()
                        output[accession] = re.sub(r"\s+\([^()]*(?:\([^()]*\)[^()]*)*\)\s*$", "", value)
                        break
        if accession in output:
            continue
        description_candidates = [
            Path(row["durable_path"])
            for row in rows
            if row.get("source_relative_path", "").endswith(".description.html")
            and Path(row.get("durable_path", "")).exists()
        ]
        for description in description_candidates:
            text = description.read_text(encoding="utf-8", errors="replace")
            match = re.search(r"Taxonomic name:\s*([^,<]+),\s*taxonomy ID:", text, re.I)
            if match:
                output[accession] = re.sub(r"<[^>]+>", "", match.group(1)).strip()
                break
    return output


def binding_priority(row: Mapping[str, Any]) -> tuple[int, int, str]:
    origin = str(row.get("origin", ""))
    provider = str(row.get("annotation_provider", "")).lower()
    pipeline = str(row.get("annotation_pipeline", "")).lower()
    if origin == "FETCHED_NCBI_OFFICIAL":
        score = 100
    elif "ncbirefseq" in pipeline or "ncbi refseq" in provider:
        score = 85
    elif "ncbigene" in pipeline:
        score = 75
    elif "cat" in pipeline or "liftoff" in provider:
        score = 65
    elif "ensembl" in provider:
        score = 60
    elif "augustus" in pipeline:
        score = 40
    elif "xenoref" in pipeline:
        score = 20
    else:
        score = 10
    return (score, int(row.get("feature_rows", 0)), str(row.get("physical_path", "")))


def validate_and_bind(
    parsed: list[dict[str, Any]], assembly_rows: Sequence[Mapping[str, str]], mirror_rows: Sequence[Mapping[str, str]]
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, dict[str, int]]]:
    assemblies = {row["accession_version"]: row for row in assembly_rows}
    fai_by_accession: dict[str, dict[str, int]] = {}
    aliases_by_accession: dict[str, dict[str, tuple[str, int]]] = {}
    for accession, assembly in assemblies.items():
        fai = load_fai(Path(assembly["fai_path"]))
        fai_by_accession[accession] = fai
        report_path = find_assembly_report(mirror_rows, accession)
        aliases, _ = parse_assembly_report(report_path, fai) if report_path else ({}, ["ASSEMBLY_REPORT_MISSING"])
        aliases_by_accession[accession] = aliases

    for row in parsed:
        hard_errors = sum(int(row.get(field, 0)) for field in HARD_PARSE_FIELDS)
        digest_matches = row.get("actual_sha256") == row.get("expected_sha256")
        row["digest_verified"] = digest_matches
        row["parse_validation_status"] = "PASS" if hard_errors == 0 and digest_matches else "FAIL"
        accession = row["assembly_accession_version"]
        if accession not in assemblies:
            row.update(
                {
                    "binding_status": "MISSING",
                    "binding_reason": f"assembly {accession} absent from 581-row BGZF manifest",
                    "target_assembly_accession_version": accession,
                }
            )
        else:
            row.update(
                bind_annotation(row, assemblies[accession], fai_by_accession[accession], aliases_by_accession[accession])
            )
        row["accepted"] = row["parse_validation_status"] == "PASS" and row["binding_status"] in {
            "EXACT_DICTIONARY",
            "VALIDATED_ALIAS",
        }

    assembly_bindings: list[dict[str, Any]] = []
    candidates: dict[str, list[dict[str, Any]]] = collections.defaultdict(list)
    for row in parsed:
        if row.get("accepted"):
            candidates[row["target_assembly_accession_version"]].append(row)
    for accession, assembly in sorted(assemblies.items()):
        available = sorted(candidates.get(accession, []), key=binding_priority, reverse=True)
        if available:
            chosen = available[0]
            assembly_bindings.append(
                {
                    "assembly_accession_version": accession,
                    "binding_status": chosen["binding_status"],
                    "annotation_object_id": chosen["object_id"],
                    "annotation_path": chosen["physical_path"],
                    "annotation_view_paths": chosen.get("view_paths", []),
                    "annotation_accession_version": chosen.get("annotation_accession_version", ""),
                    "annotation_provider": chosen.get("annotation_provider", ""),
                    "annotation_format": chosen.get("format", ""),
                    "annotation_feature_rows": chosen.get("feature_rows", 0),
                    "assembly_bgzf_path": assembly["derived_path"],
                    "assembly_fai_path": assembly["fai_path"],
                    "accepted_candidate_count": len(available),
                    "reason": chosen["binding_reason"],
                }
            )
        else:
            rejected = [row for row in parsed if row.get("assembly_accession_version") == accession]
            statuses = collections.Counter(row.get("binding_status", "") for row in rejected)
            assembly_bindings.append(
                {
                    "assembly_accession_version": accession,
                    "binding_status": "MISSING",
                    "annotation_object_id": "",
                    "annotation_path": "",
                    "annotation_view_paths": [],
                    "annotation_accession_version": "",
                    "annotation_provider": "",
                    "annotation_format": "",
                    "annotation_feature_rows": 0,
                    "assembly_bgzf_path": assembly["derived_path"],
                    "assembly_fai_path": assembly["fai_path"],
                    "accepted_candidate_count": 0,
                    "reason": f"no accepted annotation; rejected binding statuses {dict(statuses)}",
                }
            )
    return parsed, assembly_bindings, fai_by_accession


def build_pilot_bindings(
    pilot_rows: Sequence[Mapping[str, str]], assembly_bindings: Sequence[Mapping[str, Any]]
) -> list[dict[str, Any]]:
    by_accession = {row["assembly_accession_version"]: row for row in assembly_bindings}
    output: list[dict[str, Any]] = []
    for pilot in pilot_rows:
        annotation_label = pilot.get("native_annotation_accession", "")
        annotation_reference = annotation_label.split("-")[0] if annotation_label and annotation_label != "none" else ""
        chosen = by_accession.get(annotation_reference)
        if chosen:
            status = chosen["binding_status"]
            reason = f"native annotation reference {annotation_reference} bound to its exact Freeze 1 assembly"
        elif not annotation_reference:
            status = "MISSING"
            reason = "pilot roster explicitly records no native annotation"
        else:
            status = "MISSING"
            reason = f"native annotation reference {annotation_reference} absent from Freeze 1 assembly bindings"
        output.append(
            {
                "selection_id": pilot.get("selection_id", ""),
                "species": pilot.get("species", ""),
                "pair_h1_accession_version": pilot.get("h1_accession_version", ""),
                "pair_h2_accession_version": pilot.get("h2_accession_version", ""),
                "native_annotation_accession_version": annotation_label,
                "annotation_reference_accession_version": annotation_reference,
                "binding_status": status,
                "annotation_path": chosen.get("annotation_path", "") if chosen else "",
                "assembly_bgzf_path": chosen.get("assembly_bgzf_path", "") if chosen else "",
                "assembly_fai_path": chosen.get("assembly_fai_path", "") if chosen else "",
                "reason": reason,
            }
        )
    return output


CATALOG_FIELDS = [
    "object_id",
    "origin",
    "physical_path",
    "source_relative_path",
    "source_url",
    "size_bytes",
    "actual_sha256",
    "digest_verified",
    "upstream_checksum_algorithm",
    "upstream_checksum",
    "cas_path",
    "view_paths",
    "payload_path_multiplicity",
    "format",
    "compression",
    "assembly_accession_version",
    "embedded_assembly_accession_version",
    "target_assembly_accession_version",
    "species",
    "annotation_accession_version",
    "annotation_provider",
    "annotation_pipeline",
    "annotation_release_date",
    "feature_rows",
    "feature_counts",
    "feature_source_counts",
    "uncompressed_bytes",
    "observed_sequence_count",
    "embedded_sequence_count",
    "circular_sequences",
    "parse_validation_status",
    "accepted",
    "column_errors",
    "coordinate_errors",
    "cds_phase_errors",
    "gtf_parent_attribute_errors",
    "parent_reference_count",
    "defined_id_count",
    "missing_parent_references",
    "duplicate_declared_sequence_regions",
    "utf8_errors",
    "decompression_errors",
    "binding_status",
    "binding_reason",
    "exact_name_count",
    "validated_alias_count",
    "unresolved_sequence_count",
    "length_mismatch_count",
    "coordinate_overrun_count",
    "circular_coordinate_wrap_count",
    "alias_map",
    "assembly_bgzf_path",
    "assembly_fai_path",
    "assembly_bgzf_sha256",
    "assembly_fai_sha256",
]

ASSEMBLY_BINDING_FIELDS = [
    "assembly_accession_version",
    "binding_status",
    "annotation_object_id",
    "annotation_path",
    "annotation_view_paths",
    "annotation_accession_version",
    "annotation_provider",
    "annotation_format",
    "annotation_feature_rows",
    "assembly_bgzf_path",
    "assembly_fai_path",
    "accepted_candidate_count",
    "reason",
]

PILOT_FIELDS = [
    "selection_id",
    "species",
    "pair_h1_accession_version",
    "pair_h2_accession_version",
    "native_annotation_accession_version",
    "annotation_reference_accession_version",
    "binding_status",
    "annotation_path",
    "assembly_bgzf_path",
    "assembly_fai_path",
    "reason",
]

DICTIONARY_FIELDS = [
    "annotation_object_id",
    "physical_path",
    "assembly_accession_version",
    "sequence_name",
    "declared_length",
    "observed_max_end",
    "bound_sequence_name",
    "assembly_length",
    "name_binding",
]


def sequence_dictionary_rows(
    catalog: Sequence[Mapping[str, Any]], fai_by_accession: Mapping[str, Mapping[str, int]]
) -> Iterator[dict[str, Any]]:
    for row in catalog:
        accession = str(row["assembly_accession_version"])
        fai = fai_by_accession.get(accession, {})
        observed = row.get("observed_sequence_max_end", {})
        declared = row.get("declared_sequence_dictionary", {})
        aliases = row.get("alias_map", {})
        for name in sorted(set(observed) | set(declared)):
            bound = name if name in fai else aliases.get(name, "")
            yield {
                "annotation_object_id": row["object_id"],
                "physical_path": row["physical_path"],
                "assembly_accession_version": accession,
                "sequence_name": name,
                "declared_length": declared.get(name, ""),
                "observed_max_end": observed.get(name, ""),
                "bound_sequence_name": bound,
                "assembly_length": fai.get(bound, ""),
                "name_binding": "EXACT" if bound == name else ("VALIDATED_ALIAS" if bound else "UNRESOLVED"),
            }


def summary_document(
    catalog: Sequence[Mapping[str, Any]],
    assembly_bindings: Sequence[Mapping[str, Any]],
    pilot_bindings: Sequence[Mapping[str, Any]],
    fetched_failures: Mapping[str, str],
    metadata_count: int,
    metadata_requested: int,
) -> dict[str, Any]:
    discovered = [row for row in catalog if row["origin"] == "MIRRORED_FREEZE1"]
    fetched = [row for row in catalog if row["origin"] == "FETCHED_NCBI_OFFICIAL"]
    unique_payloads = {row["object_id"]: row for row in catalog}
    return {
        "schema_version": SCHEMA_VERSION,
        "generated_at_utc": utc_now(),
        "closed_world": {
            "mirror_annotation_objects": len(discovered),
            "mirror_annotation_bytes": sum(int(row["size_bytes"]) for row in discovered),
            "fetched_official_annotation_objects": len(fetched),
            "fetched_official_annotation_bytes": sum(int(row["size_bytes"]) for row in fetched),
            "all_physical_annotation_objects": len(catalog),
            "all_physical_annotation_bytes": sum(int(row["size_bytes"]) for row in catalog),
            "unique_payloads": len(unique_payloads),
            "unique_payload_bytes": sum(int(row["size_bytes"]) for row in unique_payloads.values()),
            "ncbi_metadata_reports": metadata_count,
            "ncbi_metadata_unresolved": metadata_requested - metadata_count,
            "ncbi_metadata_closed_world_accessions": metadata_requested,
            "official_fetch_failures": dict(sorted(fetched_failures.items())),
        },
        "formats": dict(sorted(collections.Counter(row["format"] for row in catalog).items())),
        "compression": dict(sorted(collections.Counter(row["compression"] for row in catalog).items())),
        "origins": dict(sorted(collections.Counter(row["origin"] for row in catalog).items())),
        "parse_validation": dict(
            sorted(collections.Counter(row["parse_validation_status"] for row in catalog).items())
        ),
        "annotation_binding_status": dict(sorted(collections.Counter(row["binding_status"] for row in catalog).items())),
        "accepted_annotations": sum(bool(row["accepted"]) for row in catalog),
        "assembly_accounting": {
            "total": len(assembly_bindings),
            "status_counts": dict(
                sorted(collections.Counter(row["binding_status"] for row in assembly_bindings).items())
            ),
        },
        "pilot_accounting": {
            "total": len(pilot_bindings),
            "status_counts": dict(sorted(collections.Counter(row["binding_status"] for row in pilot_bindings).items())),
        },
        "validation_errors": {
            field: sum(int(row.get(field, 0)) for row in catalog)
            for field in (
                "column_errors",
                "coordinate_errors",
                "cds_phase_errors",
                "gtf_parent_attribute_errors",
                "missing_parent_references",
                "duplicate_declared_sequence_regions",
                "utf8_errors",
                "decompression_errors",
                "length_mismatch_count",
                "coordinate_overrun_count",
                "circular_coordinate_wrap_count",
                "unresolved_sequence_count",
            )
        },
    }


def handoff_markdown(
    summary: Mapping[str, Any],
    catalog: Sequence[Mapping[str, Any]],
    assembly_bindings: Sequence[Mapping[str, Any]],
    pilot_bindings: Sequence[Mapping[str, Any]],
) -> str:
    closed = summary["closed_world"]
    lines = [
        "# VGP annotation catalog handoff",
        "",
        f"Generated `{summary['generated_at_utc']}` under schema `{summary['schema_version']}`.",
        "",
        "## Closed-world result",
        "",
        f"- Mirrored source objects: **{closed['mirror_annotation_objects']:,}** files, **{closed['mirror_annotation_bytes']:,}** bytes.",
        f"- Fetched official objects: **{closed['fetched_official_annotation_objects']:,}** files, **{closed['fetched_official_annotation_bytes']:,}** bytes.",
        f"- All catalog physical objects: **{closed['all_physical_annotation_objects']:,}** files, **{closed['all_physical_annotation_bytes']:,}** bytes.",
        f"- Unique payloads: **{closed['unique_payloads']:,}**, **{closed['unique_payload_bytes']:,}** bytes.",
        f"- Accepted, fully parsed annotations: **{summary['accepted_annotations']:,}**.",
        f"- Assembly accounting: `{json.dumps(summary['assembly_accounting']['status_counts'], sort_keys=True)}` over **{summary['assembly_accounting']['total']}** FASTAs.",
        f"- Pilot accounting: `{json.dumps(summary['pilot_accounting']['status_counts'], sort_keys=True)}` over **{summary['pilot_accounting']['total']}** pairs.",
        "",
        "The TSV catalog is the exhaustive physical-object handoff. The table below lists the single preferred accepted annotation selected for every assembly.",
        "",
        "## Ten pilot pairs",
        "",
        "| Pilot | Species | Native annotation | Status | Exact annotation path |",
        "|---|---|---|---|---|",
    ]
    for row in pilot_bindings:
        lines.append(
            f"| {row['selection_id']} | {row['species']} | `{row['native_annotation_accession_version']}` | {row['binding_status']} | `{row['annotation_path'] or 'MISSING'}` |"
        )
    lines.extend(
        [
            "",
            "## Preferred annotation path for every Freeze 1 assembly",
            "",
            "| Assembly | Status | Annotation | Physical path | BGZF assembly | FAI |",
            "|---|---|---|---|---|---|",
        ]
    )
    for row in assembly_bindings:
        lines.append(
            f"| `{row['assembly_accession_version']}` | {row['binding_status']} | `{row['annotation_accession_version'] or 'MISSING'}` | `{row['annotation_path'] or 'MISSING'}` | `{row['assembly_bgzf_path']}` | `{row['assembly_fai_path']}` |"
        )
    lines.extend(
        [
            "",
            "## Every accepted exact or alias-validated annotation location",
            "",
            "This is the human-readable path handoff for every accepted annotation, not only the preferred per-assembly selection.",
            "",
            "| Assembly | Binding | Provider/version | Format | Physical path | Accession views |",
            "|---|---|---|---|---|---|",
        ]
    )
    for row in catalog:
        if not row.get("accepted"):
            continue
        label = f"{row.get('annotation_provider', '')} / {row.get('annotation_accession_version', '')}"
        views = "; ".join(f"`{path}`" for path in row.get("view_paths", [])) or ""
        lines.append(
            f"| `{row['assembly_accession_version']}` | {row['binding_status']} | {label} | {row['format']} | `{row['physical_path']}` | {views} |"
        )
    rejected = [row for row in catalog if not row.get("accepted")]
    if rejected:
        lines.extend(
            [
                "",
                "## Parsed but rejected annotation exceptions",
                "",
                "| Assembly | Classification | Reason | Physical path |",
                "|---|---|---|---|",
            ]
        )
        for row in rejected:
            lines.append(
                f"| `{row['assembly_accession_version']}` | {row['binding_status']} | {row['binding_reason']} | `{row['physical_path']}` |"
            )
    lines.append("")
    return "\n".join(lines)


def copy_shared_manifests(repo_outputs: Sequence[Path], shared_dir: Path) -> None:
    shared_dir.mkdir(parents=True, exist_ok=True)
    for source in repo_outputs:
        destination = shared_dir / source.name
        temporary = destination.with_name(f".{destination.name}.{os.getpid()}.tmp")
        shutil.copyfile(source, temporary)
        os.replace(temporary, destination)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    repo_root = Path(__file__).resolve().parent.parent
    parser.add_argument("--vgp-root", type=Path, default=Path("/moosefs/erikg/vgp"))
    parser.add_argument("--mirror-manifest", type=Path, default=repo_root / "analysis/vgp_freeze1_mirror_manifest.tsv")
    parser.add_argument("--assembly-manifest", type=Path, default=repo_root / "analysis/vgp_freeze1_bgzf_manifest.tsv")
    parser.add_argument("--pilot-manifest", type=Path, default=repo_root / "analysis/vgp_10_pair_manifest.tsv")
    parser.add_argument("--output-dir", type=Path, default=repo_root / "analysis")
    parser.add_argument("--workers", type=int, default=min(12, os.cpu_count() or 1))
    parser.add_argument("--download-workers", type=int, default=6)
    parser.add_argument("--no-fetch", action="store_true", help="catalog existing objects without retrieving advertised official GFFs")
    parser.add_argument("--refresh-metadata", action="store_true")
    parser.add_argument("--refresh-downloads", action="store_true")
    parser.add_argument("--refresh-parse", action="store_true")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    annotation_root = args.vgp_root / "annotations"
    annotation_root.mkdir(parents=True, exist_ok=True)
    mirror_rows = read_tsv(args.mirror_manifest)
    assembly_rows = read_tsv(args.assembly_manifest)
    pilot_rows = read_tsv(args.pilot_manifest)
    if len(assembly_rows) != 581:
        raise RuntimeError(f"closed-world assembly manifest has {len(assembly_rows)} rows, expected 581")
    if len(pilot_rows) != 10:
        raise RuntimeError(f"pilot manifest has {len(pilot_rows)} rows, expected 10")
    accessions = [row["accession_version"] for row in assembly_rows]
    metadata_path = annotation_root / "metadata" / "ncbi-datasets-v2.json"
    reports = fetch_ncbi_metadata(accessions, metadata_path, args.refresh_metadata)

    fetched: list[dict[str, Any]] = []
    fetch_failures: dict[str, str] = {}
    if not args.no_fetch:
        fetched, fetch_failures = fetch_official_annotations(
            reports, args.vgp_root, annotation_root, args.download_workers, args.refresh_downloads
        )
    else:
        for path in (annotation_root / "provenance").glob("*.json") if (annotation_root / "provenance").exists() else []:
            provenance = json.loads(path.read_text(encoding="utf-8"))
            if provenance.get("state") == "VERIFIED" and Path(provenance.get("cas_path", "")).exists():
                fetched.append(provenance)

    freeze_commit = mirror_rows[0].get("catalog_commit", "") if mirror_rows else ""
    fallback_species = fallback_species_by_accession(mirror_rows, accessions)
    rows = source_annotation_rows(mirror_rows, reports, freeze_commit, fallback_species)
    rows.extend(fetched_annotation_rows(fetched))
    attach_view_paths(rows, [args.vgp_root / "views", annotation_root / "by-accession"])
    rows = consolidate_logical_objects(rows)
    expected_mirror_count = sum(
        1
        for row in mirror_rows
        if row.get("object_type") == "file" and ANNOTATION_RE.search(row.get("source_relative_path", ""))
    )
    observed_mirror_count = sum(row["origin"] == "MIRRORED_FREEZE1" for row in rows)
    if observed_mirror_count != expected_mirror_count:
        raise RuntimeError(f"mirror inventory lost annotations: {observed_mirror_count} != {expected_mirror_count}")

    parsed = parse_all_annotations(rows, annotation_root / "validation" / "objects", args.workers, args.refresh_parse)
    parsed, assembly_bindings, fai_by_accession = validate_and_bind(parsed, assembly_rows, mirror_rows)
    pilot_bindings = build_pilot_bindings(pilot_rows, assembly_bindings)
    summary = summary_document(parsed, assembly_bindings, pilot_bindings, fetch_failures, len(reports), len(accessions))

    output_dir: Path = args.output_dir
    catalog_tsv = output_dir / "vgp_annotation_catalog.tsv"
    catalog_json = output_dir / "vgp_annotation_catalog.json"
    binding_tsv = output_dir / "vgp_annotation_assembly_bindings.tsv"
    pilot_tsv = output_dir / "vgp_annotation_pilot_bindings.tsv"
    dictionary_tsv = output_dir / "vgp_annotation_sequence_dictionary.tsv.gz"
    summary_json = output_dir / "vgp_annotation_catalog_summary.json"
    handoff_md = output_dir / "vgp_annotation_catalog_handoff.md"
    write_tsv(catalog_tsv, parsed, CATALOG_FIELDS)
    write_tsv(binding_tsv, assembly_bindings, ASSEMBLY_BINDING_FIELDS)
    write_tsv(pilot_tsv, pilot_bindings, PILOT_FIELDS)
    write_tsv(dictionary_tsv, sequence_dictionary_rows(parsed, fai_by_accession), DICTIONARY_FIELDS)
    json_dump(
        catalog_json,
        {
            "schema_version": SCHEMA_VERSION,
            "generated_at_utc": summary["generated_at_utc"],
            "summary": summary,
            "annotations": [{field: row.get(field, "") for field in CATALOG_FIELDS} for row in parsed],
            "assembly_bindings": assembly_bindings,
            "pilot_bindings": pilot_bindings,
            "sequence_dictionary_manifest": str(dictionary_tsv),
        },
    )
    json_dump(summary_json, summary)
    handoff_md.write_text(handoff_markdown(summary, parsed, assembly_bindings, pilot_bindings), encoding="utf-8")
    outputs = [catalog_tsv, catalog_json, binding_tsv, pilot_tsv, dictionary_tsv, summary_json, handoff_md]
    copy_shared_manifests(outputs, annotation_root / "manifests")
    print(json.dumps(summary, indent=2, sort_keys=True))

    if summary["closed_world"]["official_fetch_failures"]:
        return 2
    if summary["assembly_accounting"]["total"] != 581 or summary["pilot_accounting"]["total"] != 10:
        return 3
    if summary["parse_validation"].get("FAIL", 0):
        return 4
    if summary["assembly_accounting"]["status_counts"].get("MISSING", 0):
        return 5
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
