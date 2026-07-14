#!/usr/bin/env python3
"""Exact-reference NCBI assembly composition runner for Tier 3c.

The runner deliberately separates discovery from analysis.  Discovery can
rank current NCBI assembly records, but an analysis is performed only from a
checksum-locked manifest row naming one accession *including its version*.
Remote objects are downloaded to a sibling temporary file, verified with the
manifest SHA-256 and byte count, and atomically published before they are read.

Primary GC3 is emitted only for an annotation declared native to the exact
FASTA accession.  An absent native annotation is structured missingness;
present-but-mismatched or corrupt FASTA/GFF inputs are fatal and no result is
published.  Whole-genome GC remains independently available from an eligible
exact nuclear FASTA.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
import platform
import re
import sys
import tempfile
import urllib.parse
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple, Union

from analysis.tier3_common import (
    DNA,
    GFFAnnotation,
    Tier3ValidationError,
    fasta_dictionary,
    parse_gff,
    read_fasta,
    reconstruct_cds,
    resolve_contig_aliases,
    sha256_file,
    verify_file,
)
from analysis.tier3_manifest import load_and_validate_manifest


NCBI_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
USER_AGENT = "lewontin-paradox-tier3c/1.0 (checksum-locked research pipeline)"
ACCESSION_RE = re.compile(r"^GC[AF]_[0-9]+\.[1-9][0-9]*$")
STORE_PATH_RE = re.compile(r"^(/gnu/store/[0-9a-z]{32}-[^/]+)")
STOP_CODONS_1 = frozenset({"TAA", "TAG", "TGA"})

# Declared before pilot execution and preserved as historical checks.  The
# 2026-07-14 independent audit showed that the uncited dm6/hg38 GC3 anchors do
# not define a comparable statistic and that both upper bounds exclude valid
# published gene-level values.  Do not widen or use these bands for promotion;
# analysis/tier3c_control_audit.py applies the disclosed post-hoc exact
# cross-implementation gate.  Assembly-GC bands remain useful separate checks.
PILOT_RANGES: Mapping[str, Mapping[str, Any]] = {
    "Drosophila melanogaster": {
        "assembly_accession": "GCF_000001215.4",
        "gc3": (0.50, 0.60),
        "whole_genome_gc": (0.40, 0.44),
        "source": "historical uncited dm6 GC3 assertion (about 0.55); failed and superseded for promotion by the 2026-07-14 audited-control gate",
    },
    "Homo sapiens": {
        "assembly_accession": "GCF_000001405.40",
        "gc3": (0.47, 0.57),
        "whole_genome_gc": (0.39, 0.43),
        "source": "historical uncited hg38 GC3 assertion (about 0.52); failed and superseded for promotion by the 2026-07-14 audited-control gate",
    },
}


@dataclass(frozen=True)
class AssemblyCandidate:
    """The NCBI fields used by the frozen assembly-ranking policy."""

    accession: str
    taxid: int
    species_taxid: int
    refseq_category: str
    assembly_level: str
    release_date: str
    ftp_path: str
    assembly_name: str = ""
    suppressed: bool = False

    def __post_init__(self) -> None:
        if not ACCESSION_RE.fullmatch(self.accession):
            raise Tier3ValidationError(
                f"assembly accession must include GCA/GCF prefix and version: {self.accession!r}"
            )
        if self.taxid <= 0 or self.species_taxid <= 0:
            raise Tier3ValidationError(f"invalid NCBI taxon identifiers for {self.accession}")


def _release_ordinal(text: str) -> int:
    if not text:
        return 0
    # NCBI uses this sentinel for legacy assemblies with no usable release
    # date.  It must rank behind dated records, not make the entire species
    # undiscoverable.
    if text.startswith("1/01/01"):
        return 0
    normalized = text[:10].replace("/", "-")
    try:
        return dt.date.fromisoformat(normalized).toordinal()
    except ValueError:
        for pattern in ("%Y/%m/%d", "%Y-%m-%d", "%b %d, %Y"):
            try:
                return dt.datetime.strptime(text, pattern).date().toordinal()
            except ValueError:
                pass
    raise Tier3ValidationError(f"unparseable NCBI assembly release date {text!r}")


def _assembly_rank(candidate: AssemblyCandidate, taxon_id: Optional[int]) -> Tuple[Any, ...]:
    taxon_rank = 0
    if taxon_id is not None:
        taxon_rank = 0 if candidate.taxid == taxon_id else 1 if candidate.species_taxid == taxon_id else 2
    category = candidate.refseq_category.strip().lower()
    category_rank = {"reference genome": 0, "representative genome": 1, "na": 2, "": 2}.get(category, 3)
    level = re.sub(r"[^a-z]", "", candidate.assembly_level.lower())
    level_rank = {"completegenome": 0, "chromosome": 1, "scaffold": 2, "contig": 3}.get(level, 4)
    prefix_rank = 0 if candidate.accession.startswith("GCF_") else 1
    return (
        1 if candidate.suppressed else 0,
        taxon_rank,
        category_rank,
        prefix_rank,
        level_rank,
        -_release_ordinal(candidate.release_date),
        candidate.accession.encode("utf-8"),
    )


def rank_assemblies(
    candidates: Iterable[AssemblyCandidate], *, taxon_id: Optional[int] = None
) -> List[AssemblyCandidate]:
    """Return candidates in a deterministic, documented best-first order."""

    return sorted(candidates, key=lambda item: _assembly_rank(item, taxon_id))


def select_assembly(
    candidates: Iterable[AssemblyCandidate],
    *,
    taxon_id: Optional[int] = None,
    exact_accession: Optional[str] = None,
) -> AssemblyCandidate:
    """Choose one assembly, requiring literal accession+version when pinned."""

    candidate_list = list(candidates)
    if exact_accession is not None:
        if not ACCESSION_RE.fullmatch(exact_accession):
            raise Tier3ValidationError(
                f"exact assembly accession must include GCA/GCF prefix and version: {exact_accession!r}"
            )
        matches = [item for item in candidate_list if item.accession == exact_accession and not item.suppressed]
        if len(matches) != 1:
            raise Tier3ValidationError(
                f"expected one current exact assembly accession {exact_accession}, observed {len(matches)}"
            )
        selected = matches[0]
        if taxon_id is not None and taxon_id not in {selected.taxid, selected.species_taxid}:
            raise Tier3ValidationError(
                f"exact assembly accession {exact_accession} is not assigned to requested taxon {taxon_id}"
            )
        return selected
    ranked = [item for item in rank_assemblies(candidate_list, taxon_id=taxon_id) if not item.suppressed]
    if not ranked:
        raise Tier3ValidationError("NCBI discovery returned no current assemblies")
    return ranked[0]


def parse_ncbi_esummary(payload: Mapping[str, Any]) -> List[AssemblyCandidate]:
    """Parse the stable subset of the NCBI Assembly ESummary JSON response."""

    result = payload.get("result")
    if not isinstance(result, Mapping):
        raise Tier3ValidationError("NCBI ESummary response lacks result object")
    uids = result.get("uids", [])
    if not isinstance(uids, list):
        raise Tier3ValidationError("NCBI ESummary result.uids is not a list")
    candidates: List[AssemblyCandidate] = []
    for uid in uids:
        record = result.get(str(uid))
        if not isinstance(record, Mapping):
            raise Tier3ValidationError(f"NCBI ESummary lacks record for UID {uid}")
        accession = str(record.get("assemblyaccession") or record.get("assembly_accession") or "")
        refseq_path = str(record.get("ftppath_refseq") or record.get("ftp_path_refseq") or "")
        genbank_path = str(record.get("ftppath_genbank") or record.get("ftp_path_genbank") or "")
        # The path must belong to the accession represented by this row.  GCA
        # and GCF partners are related assemblies, not interchangeable exact
        # accession strings.
        ftp_path = (refseq_path or genbank_path) if accession.startswith("GCF_") else (genbank_path or refseq_path)
        candidates.append(
            AssemblyCandidate(
                accession=accession,
                taxid=int(record.get("taxid") or 0),
                species_taxid=int(record.get("species_taxid") or record.get("taxid") or 0),
                refseq_category=str(record.get("refseq_category") or ""),
                assembly_level=str(record.get("assemblystatus") or record.get("assembly_level") or ""),
                release_date=str(record.get("asmreleasedate_genbank") or record.get("seqreleasedate") or record.get("release_date") or ""),
                ftp_path=_https_url(ftp_path) if ftp_path else "",
                assembly_name=str(record.get("assemblyname") or record.get("assembly_name") or ""),
                suppressed=str(record.get("assemblystatus") or "").lower() in {"suppressed", "replaced"},
            )
        )
    return candidates


def _read_json_response(response: Any) -> Mapping[str, Any]:
    try:
        data = response.read()
        value = json.loads(data.decode("utf-8") if isinstance(data, bytes) else data)
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise Tier3ValidationError(f"invalid NCBI JSON response: {error}") from error
    if not isinstance(value, Mapping):
        raise Tier3ValidationError("NCBI response is not a JSON object")
    return value


def discover_assemblies(
    scientific_name: str,
    *,
    taxon_id: Optional[int] = None,
    opener: Callable[..., Any] = urllib.request.urlopen,
) -> List[AssemblyCandidate]:
    """Discover current NCBI assemblies without silently choosing one."""

    if not scientific_name.strip():
        raise Tier3ValidationError("scientific name is empty")
    organism_term = f"txid{taxon_id}[Organism:exp]" if taxon_id else f'"{scientific_name}"[Organism]'

    def search(term: str) -> List[Any]:
        search_url = f"{NCBI_EUTILS}/esearch.fcgi?" + urllib.parse.urlencode(
            {"db": "assembly", "term": term, "retmode": "json", "retmax": "500"}
        )
        request = urllib.request.Request(search_url, headers={"User-Agent": USER_AGENT})
        with opener(request, timeout=120) as response:
            search_payload = _read_json_response(response)
        values = search_payload.get("esearchresult", {}).get("idlist", [])
        if not isinstance(values, list):
            raise Tier3ValidationError("NCBI ESearch result.idlist is not a list")
        return values

    try:
        # Ask first for the small, policy-preferred RefSeq category set.  A
        # broad model-organism query can exceed 3,000 assemblies, causing the
        # reference assembly to fall outside ESearch's 500-row retrieval cap.
        preferred_term = (
            f"{organism_term} AND "
            '("reference genome"[RCAT] OR "representative genome"[RCAT])'
        )
        ids = search(preferred_term)
        if not ids:
            ids = search(organism_term)
        if not ids:
            return []
        summary_url = f"{NCBI_EUTILS}/esummary.fcgi"
        # Large model-organism searches can return 500 assembly UIDs.  POST
        # avoids an HTTP 414 while preserving the complete candidate set.
        body = urllib.parse.urlencode(
            {"db": "assembly", "id": ",".join(str(item) for item in ids), "retmode": "json"}
        ).encode("ascii")
        request = urllib.request.Request(
            summary_url,
            data=body,
            headers={"User-Agent": USER_AGENT, "Content-Type": "application/x-www-form-urlencoded"},
        )
        with opener(request, timeout=120) as response:
            return parse_ncbi_esummary(_read_json_response(response))
    except Tier3ValidationError:
        raise
    except Exception as error:
        raise Tier3ValidationError(f"NCBI assembly discovery failed: {error}") from error


def ncbi_artifact_urls(candidate: AssemblyCandidate) -> Mapping[str, str]:
    """Return canonical provider URLs for a discovered NCBI assembly."""

    if not candidate.ftp_path:
        raise Tier3ValidationError(f"assembly {candidate.accession} lacks an NCBI file path")
    root = _https_url(candidate.ftp_path).rstrip("/")
    basename = Path(urllib.parse.urlparse(root).path).name
    if not basename.startswith(candidate.accession + "_"):
        raise Tier3ValidationError(
            f"NCBI path basename {basename!r} does not encode exact accession {candidate.accession}"
        )
    return {
        "fasta": f"{root}/{basename}_genomic.fna.gz",
        "gff": f"{root}/{basename}_genomic.gff.gz",
        "assembly_report": f"{root}/{basename}_assembly_report.txt",
        "provider_checksums": f"{root}/md5checksums.txt",
    }


def _https_url(url: str) -> str:
    parsed = urllib.parse.urlparse(url)
    if parsed.scheme == "ftp":
        return urllib.parse.urlunparse(parsed._replace(scheme="https"))
    if parsed.scheme != "https":
        raise Tier3ValidationError(f"remote acquisition requires HTTPS, observed {parsed.scheme or 'no'} scheme")
    return url


def _fsync_directory(path: Path) -> None:
    try:
        descriptor = os.open(str(path), os.O_RDONLY)
    except OSError:
        return
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def acquire_verified(
    url: str,
    destination: Union[str, Path],
    *,
    expected_sha256: str,
    expected_size: int,
    opener: Callable[..., Any] = urllib.request.urlopen,
) -> bool:
    """Acquire one HTTPS object and publish it only after exact verification.

    Returns ``True`` when a new object was published and ``False`` when the
    already staged destination was valid.  A corrupt old destination is never
    trusted; it is replaced only after a complete verified retry succeeds.
    """

    if not re.fullmatch(r"[0-9a-f]{64}", expected_sha256):
        raise Tier3ValidationError("download requires a lowercase SHA-256 checksum")
    if expected_size < 0:
        raise Tier3ValidationError("download requires a non-negative byte count")
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
    if destination.is_file():
        try:
            verify_file(destination, expected_sha256, expected_size)
            return False
        except Tier3ValidationError:
            pass
    fd, temporary_name = tempfile.mkstemp(prefix=f".{destination.name}.", dir=str(destination.parent))
    temporary = Path(temporary_name)
    try:
        request = urllib.request.Request(_https_url(url), headers={"User-Agent": USER_AGENT})
        try:
            with opener(request, timeout=300) as response, os.fdopen(fd, "wb") as output:
                fd = -1
                while True:
                    chunk = response.read(1024 * 1024)
                    if not chunk:
                        break
                    output.write(chunk)
                output.flush()
                os.fsync(output.fileno())
        except Exception as error:
            raise Tier3ValidationError(f"download failed for {url}: {error}") from error
        verify_file(temporary, expected_sha256, expected_size)
        os.replace(str(temporary), str(destination))
        _fsync_directory(destination.parent)
        return True
    finally:
        if fd >= 0:
            os.close(fd)
        if temporary.exists():
            temporary.unlink()


def _artifact_local_path(artifact: Mapping[str, Any], cache_directory: Path) -> Path:
    required = {"logical_name", "uri", "sha256", "size_bytes"}
    missing = required - set(artifact)
    if missing:
        raise Tier3ValidationError(f"artifact lacks required fields: {sorted(missing)!r}")
    if not re.fullmatch(r"[0-9a-f]{64}", str(artifact["sha256"])):
        raise Tier3ValidationError("artifact SHA-256 must be 64 lowercase hexadecimal characters")
    if not isinstance(artifact["size_bytes"], int) or artifact["size_bytes"] < 0:
        raise Tier3ValidationError("artifact size_bytes must be a non-negative integer")
    parsed = urllib.parse.urlparse(str(artifact["uri"]))
    if parsed.scheme == "file":
        if parsed.netloc not in {"", "localhost"}:
            raise Tier3ValidationError(f"file URI has a non-local authority: {artifact['uri']}")
        path = Path(urllib.request.url2pathname(parsed.path))
        verify_file(path, str(artifact["sha256"]), int(artifact["size_bytes"]))
        return path
    logical_name = str(artifact["logical_name"])
    if not logical_name or Path(logical_name).name != logical_name:
        raise Tier3ValidationError(f"unsafe artifact logical_name {logical_name!r}")
    path = cache_directory / logical_name
    acquire_verified(
        str(artifact["uri"]),
        path,
        expected_sha256=str(artifact["sha256"]),
        expected_size=int(artifact["size_bytes"]),
    )
    return path


def _read_contig_mapping(path: Path) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if fields == ["annotation_contig", "fasta_contig"] or fields == ["gff_contig", "fasta_contig"]:
                continue
            if len(fields) != 2 or not fields[0] or not fields[1]:
                raise Tier3ValidationError(f"invalid contig mapping at line {line_number} in {path}")
            if fields[0] in mapping:
                raise Tier3ValidationError(f"duplicate annotation contig {fields[0]!r} in {path}")
            mapping[fields[0]] = fields[1]
    if not mapping:
        raise Tier3ValidationError("annotation contig mapping is empty")
    return mapping


def _canonical_digest(value: Any) -> str:
    data = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode("utf-8")
    return hashlib.sha256(data).hexdigest()


def _whole_genome_gc(fasta: Mapping[str, str]) -> Dict[str, Any]:
    total = sum(len(sequence) for sequence in fasta.values())
    callable_bases = sum(base in DNA for sequence in fasta.values() for base in sequence.upper())
    gc_bases = sum(base in {"G", "C"} for sequence in fasta.values() for base in sequence.upper())
    if callable_bases == 0:
        raise Tier3ValidationError("exact nuclear FASTA has no callable A/C/G/T bases")
    return {
        "status": "available",
        "value": gc_bases / callable_bases,
        "gc_bases": gc_bases,
        "callable_bases": callable_bases,
        "total_bases": total,
        "ambiguous_bases_excluded": total - callable_bases,
        "contigs": len(fasta),
    }


def _reject_duplicate_cds_coordinates(annotation: GFFAnnotation) -> None:
    """Reject repeated biological segments even when their GFF line differs."""

    for transcript in annotation.transcripts.values():
        observed = set()
        for segment in transcript.segments:
            identity = (segment.contig, segment.start, segment.end, segment.strand, segment.phase)
            if identity in observed:
                raise Tier3ValidationError(
                    f"duplicate CDS segment in transcript {transcript.transcript_id!r} at line {segment.line_number}"
                )
            observed.add(identity)


def _translated_cds(fasta: Mapping[str, str], transcript: Any) -> Tuple[str, bool]:
    sequence = reconstruct_cds(fasta, transcript)
    codons = [sequence[offset : offset + 3] for offset in range(0, len(sequence), 3)]
    if any(codon in STOP_CODONS_1 for codon in codons[:-1]):
        raise Tier3ValidationError(f"internal stop in transcript {transcript.transcript_id!r}")
    terminal_stop = bool(codons and codons[-1] in STOP_CODONS_1)
    return sequence[:-3] if terminal_stop else sequence, terminal_stop


def _validated_canonical(
    fasta: Mapping[str, str], annotation: GFFAnnotation
) -> Tuple[Dict[str, Tuple[Any, str, bool]], Dict[str, int]]:
    genes: Dict[str, List[Any]] = {}
    for transcript in annotation.transcripts.values():
        if transcript.segments:
            genes.setdefault(transcript.gene_id, []).append(transcript)
    selected: Dict[str, Tuple[Any, str, bool]] = {}
    exclusions: Dict[str, int] = {}
    for gene_id, candidates in genes.items():
        exceptional = [item for item in candidates if item.annotation_exception]
        if exceptional:
            exclusions["annotated_translation_exception"] = (
                exclusions.get("annotated_translation_exception", 0) + len(exceptional)
            )
        candidates = [item for item in candidates if not item.annotation_exception]
        if not candidates:
            exclusions["gene_without_valid_cds"] = exclusions.get("gene_without_valid_cds", 0) + 1
            continue
        provider = [item for item in candidates if item.provider_canonical]
        if len(provider) > 1:  # parse_gff normally catches this as well.
            raise Tier3ValidationError(f"gene {gene_id!r} has multiple provider-canonical transcripts")
        pool = provider if provider else candidates
        valid: List[Tuple[Any, str, bool]] = []
        for transcript in pool:
            try:
                sequence, terminal_stop = _translated_cds(fasta, transcript)
            except Tier3ValidationError as error:
                reason = "invalid_provider_canonical" if provider else "invalid_candidate_transcript"
                exclusions[reason] = exclusions.get(reason, 0) + 1
                if provider:
                    raise Tier3ValidationError(
                        f"provider-canonical transcript {transcript.transcript_id!r} is invalid: {error}"
                    ) from error
                continue
            valid.append((transcript, sequence, terminal_stop))
        if not valid:
            exclusions["gene_without_valid_cds"] = exclusions.get("gene_without_valid_cds", 0) + 1
            continue
        chosen = sorted(
            valid,
            key=lambda item: (-len(item[1]), item[0].transcript_id.encode("utf-8")),
        )[0]
        selected[gene_id] = chosen
    if not selected:
        raise Tier3ValidationError("annotation has no valid canonical protein-coding CDS")
    return selected, exclusions


def _gc3(
    fasta: Mapping[str, str], annotation: GFFAnnotation
) -> Tuple[Dict[str, Any], Dict[str, Tuple[Any, str, bool]], Dict[str, int]]:
    selected, exclusions = _validated_canonical(fasta, annotation)
    gc_bases = callable_thirds = terminal_stops = 0
    for _gene_id, (_transcript, sequence, terminal_stop) in selected.items():
        third_bases = sequence[2::3]
        gc_bases += sum(base in {"G", "C"} for base in third_bases)
        callable_thirds += len(third_bases)
        terminal_stops += int(terminal_stop)
    if callable_thirds == 0:
        raise Tier3ValidationError("canonical CDS contain no callable third positions")
    return (
        {
            "status": "available",
            "value": gc_bases / callable_thirds,
            "gc_bases": gc_bases,
            "callable_third_positions": callable_thirds,
            "genes": len(selected),
            "transcripts": len(selected),
            "terminal_stop_codons_excluded": terminal_stops,
        },
        selected,
        exclusions,
    )


def _provider_cds_audit(
    dataset_id: str,
    selected: Mapping[str, Tuple[Any, str, bool]],
    provider_sequences: Optional[Mapping[str, str]],
) -> Dict[str, Any]:
    ordered = sorted(
        (item[0].transcript_id for item in selected.values()),
        key=lambda transcript_id: hashlib.sha256((dataset_id + transcript_id).encode("utf-8")).digest(),
    )[:100]
    audit = {
        "sample_rule": "first_100_by_sha256_dataset_id_plus_transcript_id_or_all_if_fewer",
        "sampled_cds_count": len(ordered),
        "sampled_cds_mismatches": 0,
        "provider_comparison": "provider_sequences_not_deposited",
    }
    if provider_sequences is None:
        return audit
    by_transcript = {item[0].transcript_id: item for item in selected.values()}
    for transcript_id in ordered:
        if transcript_id not in provider_sequences:
            raise Tier3ValidationError(f"provider CDS FASTA lacks sampled transcript {transcript_id!r}")
        reconstructed, _terminal = _translated_cds_from_selected(by_transcript[transcript_id])
        provider = provider_sequences[transcript_id].upper()
        # Providers differ on whether a terminal stop is retained.  Normalize
        # only that explicit terminal codon, never any internal bases.
        provider_trimmed = provider[:-3] if provider[-3:] in STOP_CODONS_1 else provider
        if reconstructed != provider_trimmed:
            raise Tier3ValidationError(f"sampled CDS mismatch for transcript {transcript_id!r}")
    audit["provider_comparison"] = "exact_match"
    return audit


def _translated_cds_from_selected(selected: Tuple[Any, str, bool]) -> Tuple[str, bool]:
    _transcript, sequence, terminal = selected
    return sequence, terminal


def _store_path(executable: Union[str, Path]) -> Optional[str]:
    resolved = str(Path(executable).resolve())
    match = STORE_PATH_RE.match(resolved)
    return match.group(1) if match else None


def collect_environment_provenance(
    manifest: Optional[Mapping[str, Any]] = None, *, require_guix: bool = True
) -> Dict[str, Any]:
    """Record the realized Guix environment and every executable used."""

    environment = os.environ.get("GUIX_ENVIRONMENT")
    python_executable = str(Path(sys.executable).resolve())
    python_store = _store_path(python_executable)
    if require_guix and (not environment or not str(Path(environment).resolve()).startswith("/gnu/store/")):
        raise Tier3ValidationError("Tier 3c must run inside the pinned pure GNU Guix environment")
    if require_guix and python_store is None:
        raise Tier3ValidationError("active Python executable is not from the Guix store")
    channel_commit = None
    if manifest:
        channel_commit = manifest.get("guix", {}).get("channel", {}).get("commit")
    return {
        "manager": "gnu-guix",
        "guix_environment": str(Path(environment).resolve()) if environment else None,
        "channel_commit": channel_commit,
        "platform": platform.platform(),
        "tools": {
            "python3": {
                "version": platform.python_version(),
                "executable": python_executable,
                "store_path": python_store,
            }
        },
    }


def _validate_environment_record(environment: Mapping[str, Any]) -> None:
    if environment.get("manager") != "gnu-guix":
        raise Tier3ValidationError("output environment manager must be gnu-guix")
    if not str(environment.get("guix_environment", "")).startswith("/gnu/store/"):
        raise Tier3ValidationError("output must record the realized Guix store environment")
    tools = environment.get("tools")
    if not isinstance(tools, Mapping) or not tools:
        raise Tier3ValidationError("output must record tool versions and Guix store paths")
    for name, record in tools.items():
        if not isinstance(record, Mapping) or not record.get("version"):
            raise Tier3ValidationError(f"tool {name!r} lacks a recorded version")
        if not str(record.get("store_path", "")).startswith("/gnu/store/"):
            raise Tier3ValidationError(f"tool {name!r} lacks a Guix store path")


def atomic_write_json(value: Mapping[str, Any], destination: Union[str, Path]) -> bool:
    """Write canonical JSON durably and atomically; leave identical output alone."""

    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    data = (json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False) + "\n").encode("utf-8")
    if destination.is_file() and destination.read_bytes() == data:
        return False
    fd, temporary_name = tempfile.mkstemp(prefix=f".{destination.name}.", dir=str(destination.parent))
    temporary = Path(temporary_name)
    try:
        with os.fdopen(fd, "wb") as handle:
            fd = -1
            handle.write(data)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(str(temporary), str(destination))
        _fsync_directory(destination.parent)
        return True
    finally:
        if fd >= 0:
            os.close(fd)
        if temporary.exists():
            temporary.unlink()


def _annotation_provenance_base(annotation: Mapping[str, Any]) -> Dict[str, Any]:
    return {
        "provider": annotation.get("provider"),
        "release": annotation.get("release"),
        "assembly_accession": annotation.get("assembly_accession"),
        "status": annotation.get("status"),
        "native_vs_projected": annotation.get("native_vs_projected", annotation.get("status")),
        "genetic_code": annotation.get("genetic_code"),
    }


def analyze_dataset(
    dataset: Mapping[str, Any],
    output: Union[str, Path],
    *,
    cache_directory: Optional[Union[str, Path]] = None,
    environment: Optional[Mapping[str, Any]] = None,
    manifest: Optional[Mapping[str, Any]] = None,
) -> Dict[str, Any]:
    """Analyze one exact-manifest dataset and atomically emit its JSON record."""

    output = Path(output)
    cache = Path(cache_directory) if cache_directory else output.parent / ".tier3c-cache"
    dataset_id = dataset.get("dataset_id")
    if not isinstance(dataset_id, str) or not dataset_id:
        raise Tier3ValidationError("Tier 3c dataset lacks a stable dataset_id")
    species = dataset.get("species")
    if not isinstance(species, Mapping) or not species.get("scientific_name"):
        raise Tier3ValidationError("Tier 3c dataset lacks a scientific species name")
    reference = dataset.get("reference")
    if not isinstance(reference, Mapping):
        raise Tier3ValidationError("Tier 3c dataset lacks an exact reference record")
    accession = str(reference.get("assembly_accession") or "")
    if not ACCESSION_RE.fullmatch(accession):
        raise Tier3ValidationError("reference assembly accession must include exact GCA/GCF version")
    if not reference.get("provider") or not reference.get("release"):
        raise Tier3ValidationError("exact reference must record provider and release")
    fasta_artifact = reference.get("fasta")
    if not isinstance(fasta_artifact, Mapping):
        raise Tier3ValidationError("Tier 3c reference lacks a checksum-locked FASTA artifact")

    annotation = dataset.get("annotation")
    if annotation is not None:
        if not isinstance(annotation, Mapping):
            raise Tier3ValidationError("annotation record is not an object")
        annotation_accession = str(annotation.get("assembly_accession") or "")
        if annotation_accession != accession:
            raise Tier3ValidationError(
                f"FASTA/GFF assembly accession mismatch: {accession} != {annotation_accession}"
            )
        if annotation.get("exact_reference_assertion") is not True:
            raise Tier3ValidationError("annotation lacks an exact-reference assertion")
        if not annotation.get("provider") or not annotation.get("release"):
            raise Tier3ValidationError("annotation must record provider and release")

    fasta_path = _artifact_local_path(fasta_artifact, cache)
    fasta = read_fasta(fasta_path)
    if reference.get("nuclear_contigs_only") is not True:
        raise Tier3ValidationError(
            "whole-genome GC requires a manifest assertion that the FASTA contains nuclear contigs only"
        )
    genome_gc = _whole_genome_gc(fasta)
    fasta_dict = fasta_dictionary(fasta)
    reference_record: Dict[str, Any] = {
        "accession": accession,
        "provider": reference.get("provider"),
        "release": reference.get("release"),
        "fasta_uri": fasta_artifact["uri"],
        "fasta_sha256": fasta_artifact["sha256"],
        "fasta_size_bytes": fasta_artifact["size_bytes"],
        "contig_dictionary_sha256": _canonical_digest(fasta_dict),
        "contig_dictionary": fasta_dict,
        "nuclear_contigs_only": True,
    }

    gc3: Dict[str, Any]
    annotation_record: Optional[Dict[str, Any]] = None
    if annotation is None:
        gc3 = {"status": "unavailable", "reason": "native_annotation_absent"}
    else:
        annotation_record = _annotation_provenance_base(annotation)
        if annotation.get("status") != "native":
            gc3 = {
                "status": "unavailable",
                "reason": "annotation_not_native",
                "annotation_status": annotation.get("status"),
            }
        else:
            gff_artifact = annotation.get("file")
            mapping_artifact = annotation.get("contig_mapping")
            if not isinstance(gff_artifact, Mapping) or not isinstance(mapping_artifact, Mapping):
                raise Tier3ValidationError("native annotation requires checksum-locked GFF and contig mapping")
            gff_path = _artifact_local_path(gff_artifact, cache)
            mapping_path = _artifact_local_path(mapping_artifact, cache)
            # Native annotations without any CDS still support exact-reference
            # whole-genome GC.  Preserve their dictionaries/provenance and
            # report annotation-derived GC3 as structured unavailability.
            parsed = parse_gff(gff_path, require_cds=False)
            _reject_duplicate_cds_coordinates(parsed)
            aliases = _read_contig_mapping(mapping_path)
            resolved = resolve_contig_aliases(fasta_dict, parsed.sequence_regions, aliases)
            annotation_fasta = dict(fasta)
            annotation_fasta.update(
                {annotation_name: fasta[fasta_name] for annotation_name, fasta_name in resolved.items()}
            )
            sequence_regions_digest = _canonical_digest(parsed.sequence_regions)
            sequence_regions_artifact = annotation.get("sequence_regions")
            if sequence_regions_artifact is not None:
                if not isinstance(sequence_regions_artifact, Mapping):
                    raise Tier3ValidationError("annotation sequence_regions artifact is invalid")
                _artifact_local_path(sequence_regions_artifact, cache)
                sequence_regions_digest = str(sequence_regions_artifact["sha256"])
            annotation_record.update(
                {
                    "gff_uri": gff_artifact["uri"],
                    "gff_sha256": gff_artifact["sha256"],
                    "gff_size_bytes": gff_artifact["size_bytes"],
                    "fasta_sha256": fasta_artifact["sha256"],
                    "sequence_regions_sha256": sequence_regions_digest,
                    "sequence_region_count": len(parsed.sequence_regions),
                    "sequence_regions": parsed.sequence_regions,
                    "contig_mapping_sha256": mapping_artifact["sha256"],
                    "contig_mapping_uri": mapping_artifact["uri"],
                    "contig_mapping_size_bytes": mapping_artifact["size_bytes"],
                    "contig_mapping": resolved,
                    "contig_dictionary_validated": True,
                    "all_retained_cds_validated": False,
                }
            )
            genetic_code = annotation.get("genetic_code")
            has_cds = any(transcript.segments for transcript in parsed.transcripts.values())
            if not has_cds:
                gc3 = {
                    "status": "unavailable",
                    "reason": "native_annotation_has_no_cds",
                }
                annotation_record["cds_audit"] = {
                    "all_retained_cds_validated": False,
                    "reason": "native_annotation_has_no_cds",
                    "retained_genes": 0,
                    "retained_transcripts": 0,
                }
            elif genetic_code != 1:
                gc3 = {
                    "status": "unavailable",
                    "reason": "unsupported_nuclear_genetic_code",
                    "genetic_code": genetic_code,
                }
                annotation_record["cds_audit"] = {
                    "all_retained_cds_validated": False,
                    "reason": "unsupported_nuclear_genetic_code",
                }
            else:
                gc3, selected, exclusions = _gc3(annotation_fasta, parsed)
                provider_sequences = None
                provider_artifact = annotation.get("provider_cds") or annotation.get("cds_fasta")
                if provider_artifact is not None:
                    if not isinstance(provider_artifact, Mapping):
                        raise Tier3ValidationError("provider CDS artifact is invalid")
                    provider_sequences = read_fasta(_artifact_local_path(provider_artifact, cache))
                    annotation_record["provider_cds_sha256"] = provider_artifact["sha256"]
                    annotation_record["provider_cds_size_bytes"] = provider_artifact["size_bytes"]
                sampled_audit = _provider_cds_audit(
                    dataset_id, selected, provider_sequences
                )
                annotation_record["all_retained_cds_validated"] = True
                annotation_record["cds_audit"] = {
                    **sampled_audit,
                    "retained_genes": len(selected),
                    "retained_transcripts": len(selected),
                    "exclusions": exclusions,
                    "all_retained_cds_phase_translation_passed": True,
                }

    environment_record = dict(environment) if environment is not None else collect_environment_provenance(manifest)
    _validate_environment_record(environment_record)
    notes = dataset.get("notes")
    if isinstance(notes, str):
        note_list = [notes] if notes else []
    elif isinstance(notes, list) and all(isinstance(item, str) for item in notes):
        note_list = list(notes)
    else:
        note_list = []
    if gc3["status"] != "available":
        note_list.append(f"gc3_unavailable:{gc3['reason']}")
    result: Dict[str, Any] = {
        "schema_version": "tier3c-composition-v1",
        "dataset_id": dataset_id,
        "species": species,
        "reference": reference_record,
        "annotation_provenance": annotation_record,
        "whole_genome_gc": genome_gc,
        "gc3": gc3,
        "environment": environment_record,
        "notes": note_list,
    }
    atomic_write_json(result, output)
    return result


def validate_pilot(species: str, result: Mapping[str, Any]) -> List[str]:
    """Return explicit deviations from a pilot's predeclared validation bands."""

    if species not in PILOT_RANGES:
        raise Tier3ValidationError(f"no predeclared Tier 3c pilot for {species!r}")
    expected = PILOT_RANGES[species]
    failures: List[str] = []
    reference = result.get("reference")
    observed_accession = reference.get("accession") if isinstance(reference, Mapping) else None
    if observed_accession != expected["assembly_accession"]:
        failures.append(
            f"reference accession {observed_accession!r} != frozen {expected['assembly_accession']}"
        )
    annotation = result.get("annotation_provenance")
    if not isinstance(annotation, Mapping) or annotation.get("status") != "native":
        failures.append("pilot lacks native exact-reference annotation provenance")
    elif annotation.get("assembly_accession") != expected["assembly_accession"]:
        failures.append("pilot annotation accession does not match frozen reference")
    for key in ("gc3", "whole_genome_gc"):
        observed_record = result.get(key)
        if not isinstance(observed_record, Mapping) or observed_record.get("status") != "available":
            failures.append(f"{key} unavailable")
            continue
        observed = observed_record.get("value")
        lower, upper = expected[key]
        if not isinstance(observed, (int, float)) or not lower <= float(observed) <= upper:
            failures.append(f"{key}={observed!r} outside predeclared [{lower}, {upper}]")
    return failures


def _candidate_json(candidate: AssemblyCandidate) -> Dict[str, Any]:
    value = dict(candidate.__dict__)
    value["artifact_urls"] = ncbi_artifact_urls(candidate) if candidate.ftp_path else None
    return value


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    discover = commands.add_parser("discover", help="rank current NCBI assembly records")
    discover.add_argument("scientific_name")
    discover.add_argument("--taxon-id", type=int)
    discover.add_argument("--exact-accession")
    discover.add_argument("--output", type=Path)
    run = commands.add_parser("run", help="analyze one checksum-locked manifest dataset")
    run.add_argument("manifest", type=Path)
    run.add_argument("dataset_id")
    run.add_argument("output", type=Path)
    run.add_argument("--cache-directory", type=Path)
    pilot = commands.add_parser("pilot-check", help="validate a completed pilot against frozen ranges")
    pilot.add_argument("species", choices=sorted(PILOT_RANGES))
    pilot.add_argument("result", type=Path)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _parser().parse_args(argv)
    if args.command == "discover":
        candidates = discover_assemblies(args.scientific_name, taxon_id=args.taxon_id)
        if args.exact_accession:
            selected = select_assembly(
                candidates, taxon_id=args.taxon_id, exact_accession=args.exact_accession
            )
            payload: Any = _candidate_json(selected)
        else:
            payload = [_candidate_json(item) for item in rank_assemblies(candidates, taxon_id=args.taxon_id)]
        if args.output:
            atomic_write_json({"assemblies": payload} if isinstance(payload, list) else payload, args.output)
        else:
            print(json.dumps(payload, indent=2, sort_keys=True))
        return 0
    if args.command == "run":
        manifest = load_and_validate_manifest(args.manifest, verify_local_files=False)
        matches = [item for item in manifest["datasets"] if item["dataset_id"] == args.dataset_id]
        if len(matches) != 1:
            raise Tier3ValidationError(
                f"manifest must contain exactly one dataset_id {args.dataset_id!r}, observed {len(matches)}"
            )
        analyze_dataset(
            matches[0],
            args.output,
            cache_directory=args.cache_directory,
            manifest=manifest,
        )
        print(f"wrote exact-reference Tier 3c composition: {args.output}")
        return 0
    result = json.loads(args.result.read_text(encoding="utf-8"))
    failures = validate_pilot(args.species, result)
    if failures:
        raise Tier3ValidationError("; ".join(failures))
    print(f"pilot passed predeclared ranges: {args.species}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Tier3ValidationError as error:
        print(f"tier3c: {error}", file=sys.stderr)
        raise SystemExit(2)
