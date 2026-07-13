#!/usr/bin/env python3
"""Collect fail-closed Tier 3a VGP individual-diversity inputs.

The deposited route preserves an exact-reference VCF/BCF plus an explicit
sample-callable mask.  The assembly fallback treats H1 as the reference and
streams an extended-CIGAR H2-to-H1 WFMASH PAF into a callable BED and a
biallelic SNV BCF.  Sparse variant records are never used as a denominator.

IMPG support in this module is deliberately limited to construction and audit
helpers.  The frozen Guix manifest currently ships source only; callers must
not execute IMPG until its pinned truth test has explicitly passed.
"""

from __future__ import annotations

import argparse
import bisect
import collections
import json
import re
import shutil
import subprocess
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple, Union

try:
    from analysis.tier3_common import (
        DNA,
        Tier3ValidationError,
        collect_fourfold_sites,
        fasta_dictionary,
        merge_intervals,
        parse_extended_cigar,
        parse_gff,
        read_fasta,
        resolve_contig_aliases,
        reverse_complement,
        sha256_file,
        traverse_paf,
    )
except ModuleNotFoundError:  # direct ``python analysis/<script>.py`` execution
    from tier3_common import (  # type: ignore[no-redef]
        DNA,
        Tier3ValidationError,
        collect_fourfold_sites,
        fasta_dictionary,
        merge_intervals,
        parse_extended_cigar,
        parse_gff,
        read_fasta,
        resolve_contig_aliases,
        reverse_complement,
        sha256_file,
        traverse_paf,
    )


WFMASH_COMMIT = "e040aa10e87cab44ed5a4db005e784be62b0bd21"
IMPG_COMMIT = "101df81eb28a809c8fac97d297acd9fcfbbfa048"
WFMASH_POLICY_ARGUMENTS = ("-p", "90", "-s", "5k", "-l", "25k", "-4")
DEFAULT_EDGE_EXCLUSION_BP = 100
DEFAULT_INDEL_FLANK_BP = 10
DEFAULT_IMPG_CORE_BP = 1_000_000
DEFAULT_IMPG_PADDING_BP = 10_000
_STORE_PATH = re.compile(r"^/gnu/store/[0-9a-z]{32}-[^/]+$")


@dataclass(frozen=True)
class ImpgCoreWindow:
    contig: str
    core_start: int
    core_end: int
    query_start: int
    query_end: int


def _path(path: Union[str, Path]) -> Path:
    return Path(path).expanduser().resolve()


def _validate_store_root(path: Union[str, Path], tool: str, *, executable_required: bool = True) -> Path:
    root = _path(path)
    if not _STORE_PATH.match(str(root)):
        raise Tier3ValidationError(f"{tool} must be recorded as a pinned Guix store path: {root}")
    executable = root / "bin" / tool
    if executable_required and not executable.is_file():
        raise Tier3ValidationError(f"pinned {tool} executable is absent: {executable}")
    return root


def _artifact(path: Path) -> Dict[str, Any]:
    return {
        "path": str(path.resolve()),
        "sha256": sha256_file(path),
        "size_bytes": path.stat().st_size,
    }


def _write_json(path: Path, value: Mapping[str, Any]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(path)


def _read_bed(
    path: Union[str, Path], reference: Mapping[str, str]
) -> Dict[str, List[Tuple[int, int]]]:
    intervals: Dict[str, List[Tuple[int, int]]] = collections.defaultdict(list)
    with Path(path).open("r", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 3 or fields[0] not in reference:
                raise Tier3ValidationError(f"invalid callable BED contig at line {line_number}")
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError as error:
                raise Tier3ValidationError(f"non-integer callable BED coordinate at line {line_number}") from error
            if not 0 <= start < end <= len(reference[fields[0]]):
                raise Tier3ValidationError(f"callable BED out of FASTA bounds at line {line_number}")
            intervals[fields[0]].append((start, end))
    if not intervals:
        raise Tier3ValidationError("callable BED has no intervals")
    return {contig: merge_intervals(values) for contig, values in intervals.items()}


def _position_in_intervals(position: int, intervals: Sequence[Tuple[int, int]]) -> bool:
    starts = [start for start, _end in intervals]
    index = bisect.bisect_right(starts, position) - 1
    return index >= 0 and position < intervals[index][1]


def _positions_to_bed(path: Path, positions: Set[Tuple[str, int]]) -> int:
    count = 0
    with path.open("w", encoding="utf-8") as handle:
        current_contig: Optional[str] = None
        current_start = previous = -1
        for contig, position in sorted(positions):
            count += 1
            if contig == current_contig and position == previous + 1:
                previous = position
                continue
            if current_contig is not None:
                handle.write(f"{current_contig}\t{current_start}\t{previous + 1}\n")
            current_contig, current_start, previous = contig, position, position
        if current_contig is not None:
            handle.write(f"{current_contig}\t{current_start}\t{previous + 1}\n")
    return count


def _paf_snv_alleles(
    paf_lines: Iterable[str],
    target_fasta: Mapping[str, str],
    query_fasta: Mapping[str, str],
    callable_positions: Set[Tuple[str, int]],
) -> List[Dict[str, Any]]:
    """Recover H2 alleles only after ``traverse_paf`` accepted each base."""

    projections: Dict[Tuple[str, int], List[str]] = collections.defaultdict(list)
    for line_number, raw in enumerate(paf_lines, 1):
        if not raw.strip() or raw.startswith("#"):
            continue
        fields = raw.rstrip("\n").split("\t")
        qname, qstart, qend, strand = fields[0], int(fields[2]), int(fields[3]), fields[4]
        tname, tstart = fields[5], int(fields[7])
        cigar_tags = [field[5:] for field in fields[12:] if field.startswith("cg:Z:")]
        if len(cigar_tags) != 1:
            raise Tier3ValidationError(f"PAF line {line_number} must have exactly one cg:Z tag")
        query_piece = query_fasta[qname][qstart:qend].upper()
        if strand == "-":
            query_piece = reverse_complement(query_piece)
        query_offset, target_position = 0, tstart
        for length, operation in parse_extended_cigar(cigar_tags[0]):
            if operation in {"=", "X"}:
                if operation == "X":
                    for delta in range(length):
                        projections[(tname, target_position + delta)].append(
                            query_piece[query_offset + delta]
                        )
                query_offset += length
                target_position += length
            elif operation == "I":
                query_offset += length
            else:
                target_position += length

    variants: List[Dict[str, Any]] = []
    for key in sorted(callable_positions):
        alleles = projections.get(key, [])
        if not alleles:
            continue
        if len(alleles) != 1:
            raise Tier3ValidationError(f"accepted SNV {key!r} unexpectedly has multiple H2 alleles")
        ref = target_fasta[key[0]][key[1]].upper()
        alt = alleles[0].upper()
        if ref not in DNA or alt not in DNA or ref == alt:
            raise Tier3ValidationError(f"invalid accepted X allele at {key[0]}:{key[1] + 1}")
        variants.append(
            {"contig": key[0], "position_1based": key[1] + 1, "ref": ref, "alt": alt}
        )
    return variants


def _nonunique_query_targets(paf_lines: Iterable[str]) -> Set[Tuple[str, int]]:
    """Return H1 bases reached by an H2 base used in multiple mappings."""

    query_to_targets: Dict[Tuple[str, int], List[Tuple[str, int]]] = collections.defaultdict(list)
    for line_number, raw in enumerate(paf_lines, 1):
        if not raw.strip() or raw.startswith("#"):
            continue
        fields = raw.rstrip("\n").split("\t")
        if any(field == "tp:A:S" for field in fields[12:]):
            raise Tier3ValidationError(f"secondary PAF mapping is forbidden at line {line_number}")
        qname, qstart, qend, strand = fields[0], int(fields[2]), int(fields[3]), fields[4]
        tname, target_position = fields[5], int(fields[7])
        cigar_tags = [field[5:] for field in fields[12:] if field.startswith("cg:Z:")]
        if len(cigar_tags) != 1:
            raise Tier3ValidationError(f"PAF line {line_number} must have exactly one cg:Z tag")
        query_offset = 0
        for length, operation in parse_extended_cigar(cigar_tags[0]):
            if operation in {"=", "X"}:
                for delta in range(length):
                    oriented_offset = query_offset + delta
                    query_position = (
                        qstart + oriented_offset
                        if strand == "+"
                        else qend - 1 - oriented_offset
                    )
                    query_to_targets[(qname, query_position)].append(
                        (tname, target_position + delta)
                    )
                query_offset += length
                target_position += length
            elif operation == "I":
                query_offset += length
            else:
                target_position += length
    return {
        target
        for targets in query_to_targets.values()
        if len(targets) != 1
        for target in targets
    }


def _write_snv_bcf(
    path: Path,
    index_path: Path,
    reference: Mapping[str, str],
    variants: Sequence[Mapping[str, Any]],
    sample: str,
    bcftools_root: Path,
) -> None:
    try:
        import pysam
    except ImportError as error:  # pragma: no cover - supplied by the Guix manifest
        raise Tier3ValidationError("pysam is required; enter the pinned pure Guix environment") from error

    header = pysam.VariantHeader()
    header.add_meta("fileformat", value="VCFv4.3")
    header.add_meta(
        "source",
        value=f"tier3a-direct-wfmash-{WFMASH_COMMIT}",
    )
    header.add_meta(
        "INFO",
        items=[("ID", "T3SRC"), ("Number", "1"), ("Type", "String"),
               ("Description", "Tier 3 source modality")],
    )
    header.add_meta(
        "FORMAT",
        items=[("ID", "GT"), ("Number", "1"), ("Type", "String"),
               ("Description", "H1 reference and H2 alternate haplotype")],
    )
    for contig, sequence in reference.items():
        header.contigs.add(contig, length=len(sequence))
    header.add_sample(sample)
    with pysam.VariantFile(str(path), "wb", header=header) as output:
        for item in variants:
            start = int(item["position_1based"]) - 1
            record = output.new_record(
                contig=str(item["contig"]),
                start=start,
                stop=start + 1,
                alleles=(str(item["ref"]), str(item["alt"])),
                qual=60,
                filter="PASS",
            )
            record.info["T3SRC"] = "direct_wfmash_extended_cigar"
            record.samples[sample]["GT"] = (0, 1)
            record.samples[sample].phased = True
            output.write(record)
    completed = subprocess.run(
        [str(bcftools_root / "bin" / "bcftools"), "index", "--force", "--csi", str(path)],
        text=True,
        capture_output=True,
        check=False,
    )
    if completed.returncode:
        raise Tier3ValidationError(f"bcftools index failed: {completed.stderr.strip()}")
    generated = Path(str(path) + ".csi")
    if generated != index_path:
        generated.replace(index_path)
    if not index_path.is_file():
        raise Tier3ValidationError("bcftools did not create the required CSI index")


def run_pinned_wfmash(
    h1_fasta: Union[str, Path],
    h2_fasta: Union[str, Path],
    output_paf: Union[str, Path],
    wfmash_store_path: Union[str, Path],
    *,
    threads: int = 2,
) -> List[str]:
    """Run the frozen H2-query/H1-target WFMASH command without a shell."""

    if threads < 1:
        raise Tier3ValidationError("WFMASH threads must be positive")
    root = _validate_store_root(wfmash_store_path, "wfmash")
    command = [
        str(root / "bin" / "wfmash"),
        str(_path(h1_fasta)),
        str(_path(h2_fasta)),
        *WFMASH_POLICY_ARGUMENTS,
        "-t",
        str(threads),
    ]
    output_paf = Path(output_paf)
    with output_paf.open("w", encoding="utf-8") as handle:
        completed = subprocess.run(command, stdout=handle, stderr=subprocess.PIPE, text=True, check=False)
    if completed.returncode:
        raise Tier3ValidationError(f"pinned WFMASH failed: {completed.stderr.strip()}")
    if not output_paf.is_file() or output_paf.stat().st_size == 0:
        raise Tier3ValidationError("pinned WFMASH produced no PAF records")
    return command


def _audit_variant_reference(
    variants_path: Union[str, Path],
    reference: Mapping[str, str],
    *,
    selected_sample: Optional[str] = None,
) -> Dict[str, Any]:
    try:
        import pysam
    except ImportError as error:  # pragma: no cover - supplied by the Guix manifest
        raise Tier3ValidationError("pysam is required; enter the pinned pure Guix environment") from error

    records = ref_bases = 0
    with pysam.VariantFile(str(variants_path)) as variants:
        observed_dictionary = {
            name: variants.header.contigs[name].length for name in variants.header.contigs
        }
        if set(observed_dictionary) != set(reference):
            raise Tier3ValidationError("VCF/BCF contig dictionary does not exactly match reference FASTA")
        for contig, sequence in reference.items():
            if observed_dictionary[contig] != len(sequence):
                raise Tier3ValidationError(f"VCF/BCF contig length mismatch for {contig!r}")
        samples = list(variants.header.samples)
        if selected_sample is not None and selected_sample not in samples:
            raise Tier3ValidationError(f"sample {selected_sample!r} is absent from VCF/BCF")
        for record in variants:
            if record.contig not in reference:
                raise Tier3ValidationError(f"VCF contig {record.contig!r} absent from reference FASTA")
            observed = reference[record.contig][record.start : record.start + len(record.ref)].upper()
            if observed != record.ref.upper():
                raise Tier3ValidationError(
                    f"VCF REF mismatch at {record.contig}:{record.pos}: {record.ref} != {observed}"
                )
            records += 1
            ref_bases += len(record.ref)
    return {
        "records_audited": records,
        "ref_bases_audited": ref_bases,
        "contig_dictionary_passed": True,
        "ref_allele_audit_passed": True,
        "samples": samples,
    }


def collect_deposited_variants(
    reference_fasta: Union[str, Path],
    deposited_variants: Union[str, Path],
    callable_bed: Optional[Union[str, Path]],
    output_dir: Union[str, Path],
    sample: str,
    bcftools_store_path: Union[str, Path],
    *,
    reference_assembly_accession: str,
    variant_reference_accession: str,
    annotation_provenance: Optional[Mapping[str, Any]] = None,
    annotation_gff: Optional[Union[str, Path]] = None,
) -> Dict[str, Any]:
    """Preserve and normalize a deposited exact-reference individual tuple.

    Variant-only input is returned as explicitly unavailable.  It is never
    normalized into an apparently eligible row with a fabricated invariant
    denominator.
    """

    if reference_assembly_accession != variant_reference_accession:
        raise Tier3ValidationError(
            "deposited variant reference accession does not match the exact FASTA assembly"
        )
    reference_path = _path(reference_fasta)
    variants_path = _path(deposited_variants)
    reference = read_fasta(reference_path)
    if annotation_provenance:
        _audit_annotation_provenance(
            annotation_provenance,
            reference_path,
            reference,
            reference_assembly_accession,
            annotation_gff,
        )
    original_qc = _audit_variant_reference(variants_path, reference, selected_sample=sample)
    if callable_bed is None:
        return {
            "status": "unavailable_missing_callable_denominator",
            "modality": "deposited_exact_reference_variants_plus_mask",
            "statistic_label": "individual_snv_heterozygosity",
            "population_pi": None,
            "callable_bases": None,
            "variant_only_reference_length_assumption_used": False,
            "original_variant_qc": original_qc,
        }
    callable_path = _path(callable_bed)
    intervals = _read_bed(callable_path, reference)
    callable_bases = sum(
        end - start for values in intervals.values() for start, end in values
    )
    if callable_bases == 0:
        raise Tier3ValidationError("deposited callable BED has no invariant denominator")
    bcftools_root = _validate_store_root(bcftools_store_path, "bcftools")
    output = _path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    reference_copy = output / "reference.fa"
    original_copy = output / ("deposited.original" + "".join(variants_path.suffixes))
    mask_copy = output / "deposited.callable.bed"
    if reference_path != reference_copy:
        shutil.copyfile(reference_path, reference_copy)
    if variants_path != original_copy:
        shutil.copyfile(variants_path, original_copy)
    if callable_path != mask_copy:
        shutil.copyfile(callable_path, mask_copy)
    normalized = output / "deposited.normalized.bcf"
    normalized_index = output / "deposited.normalized.bcf.csi"
    command = [
        str(bcftools_root / "bin" / "bcftools"),
        "norm",
        "--fasta-ref",
        str(reference_copy),
        "--multiallelics",
        "-any",
        "--output-type",
        "b",
        "--output",
        str(normalized),
        str(variants_path),
    ]
    completed = subprocess.run(command, text=True, capture_output=True, check=False)
    if completed.returncode:
        raise Tier3ValidationError(f"bcftools norm failed: {completed.stderr.strip()}")
    completed = subprocess.run(
        [
            str(bcftools_root / "bin" / "bcftools"),
            "index",
            "--force",
            "--csi",
            str(normalized),
        ],
        text=True,
        capture_output=True,
        check=False,
    )
    if completed.returncode:
        raise Tier3ValidationError(f"bcftools index failed: {completed.stderr.strip()}")
    generated_index = Path(str(normalized) + ".csi")
    if generated_index != normalized_index:
        generated_index.replace(normalized_index)
    normalized_qc = _audit_variant_reference(normalized, reference, selected_sample=sample)
    provenance_path = output / "provenance.json"
    record: Dict[str, Any] = {
        "status": "eligible",
        "modality": "deposited_exact_reference_variants_plus_mask",
        "statistic_label": "individual_snv_heterozygosity",
        "population_pi": None,
        "reference_assembly_accession": reference_assembly_accession,
        "sample": sample,
        "callable_bases_from_bed_intervals": callable_bases,
        "variant_only_reference_length_assumption_used": False,
        "annotation_provenance": dict(annotation_provenance) if annotation_provenance else None,
        "tools": {"bcftools_store_path": str(bcftools_root)},
        "normalization_command": command,
        "original_variant_qc": original_qc,
        "normalized_variant_qc": normalized_qc,
        "outputs": {
            "reference_fasta": str(reference_copy),
            "original_variants": str(original_copy),
            "callable_bed": str(mask_copy),
            "normalized_bcf": str(normalized),
            "normalized_bcf_index": str(normalized_index),
        },
        "artifacts": {
            "reference_fasta": _artifact(reference_copy),
            "original_variants": _artifact(original_copy),
            "callable_bed": _artifact(mask_copy),
            "normalized_bcf": _artifact(normalized),
            "normalized_bcf_index": _artifact(normalized_index),
        },
    }
    _write_json(provenance_path, record)
    record["provenance_path"] = str(provenance_path)
    return record


def _audit_annotation_provenance(
    provenance: Mapping[str, Any],
    h1_path: Path,
    h1_reference: Mapping[str, str],
    h1_assembly_accession: Optional[str],
    annotation_gff: Optional[Union[str, Path]],
) -> None:
    required = {
        "provider", "release", "assembly_accession_version", "status",
        "genetic_code", "fasta_sha256", "gff_sha256", "contig_mapping",
    }
    missing = sorted(required - set(provenance))
    if missing:
        raise Tier3ValidationError(f"annotation provenance lacks fields: {missing!r}")
    if provenance["status"] != "native":
        raise Tier3ValidationError("H1 primary annotation must be native, never projected")
    if provenance["assembly_accession_version"] != h1_assembly_accession:
        raise Tier3ValidationError("annotation assembly does not match the H1 target assembly")
    if provenance["fasta_sha256"] != sha256_file(h1_path):
        raise Tier3ValidationError("annotation provenance FASTA checksum does not match H1")
    if int(provenance["genetic_code"]) != 1:
        raise Tier3ValidationError("only frozen nuclear genetic code 1 is supported")
    if provenance.get("audit_passed") is False or provenance.get("sampled_cds_mismatches", 0) != 0:
        raise Tier3ValidationError("annotation CDS reconstruction audit failed")
    if annotation_gff is None:
        raise Tier3ValidationError("annotation provenance requires the native H1 GFF payload")
    gff_path = _path(annotation_gff)
    if provenance["gff_sha256"] != sha256_file(gff_path):
        raise Tier3ValidationError("annotation provenance GFF checksum mismatch")
    annotation = parse_gff(gff_path)
    aliases = provenance["contig_mapping"]
    if not isinstance(aliases, Mapping):
        raise Tier3ValidationError("annotation contig_mapping provenance must be an object")
    resolved = resolve_contig_aliases(
        fasta_dictionary(h1_reference), annotation.sequence_regions, aliases
    )
    annotation_fasta = dict(h1_reference)
    for annotation_contig, fasta_contig in resolved.items():
        annotation_fasta[annotation_contig] = h1_reference[fasta_contig]
    collect_fourfold_sites(annotation_fasta, annotation, int(provenance["genetic_code"]))


def collect_direct_alignment(
    h1_fasta: Union[str, Path],
    h2_fasta: Union[str, Path],
    paf_path: Optional[Union[str, Path]],
    output_dir: Union[str, Path],
    sample: str,
    wfmash_store_path: Union[str, Path],
    bcftools_store_path: Union[str, Path],
    phase_qc_passed: bool,
    collapse_qc_passed: bool,
    edge_exclusion_bp: int = DEFAULT_EDGE_EXCLUSION_BP,
    indel_flank_bp: int = DEFAULT_INDEL_FLANK_BP,
    *,
    h1_accessibility_bed: Optional[Union[str, Path]] = None,
    h1_assembly_accession: Optional[str] = None,
    h2_assembly_accession: Optional[str] = None,
    annotation_assembly_accession: Optional[str] = None,
    annotation_provenance: Optional[Mapping[str, Any]] = None,
    annotation_gff: Optional[Union[str, Path]] = None,
    threads: int = 2,
) -> Dict[str, Any]:
    """Create the complete direct-alignment artifact tuple.

    ``phase_qc_passed`` and ``collapse_qc_passed`` are explicit upstream QC
    decisions, not inferred from a good-looking alignment.  A false value is a
    hard eligibility failure.
    """

    if not sample or any(character.isspace() for character in sample):
        raise Tier3ValidationError("sample must be a stable non-whitespace identifier")
    if not phase_qc_passed:
        raise Tier3ValidationError("phase/identity QC failed for the H1/H2 individual")
    if not collapse_qc_passed:
        raise Tier3ValidationError("collapse/duplication QC failed for the H1/H2 individual")
    if annotation_assembly_accession and annotation_assembly_accession != h1_assembly_accession:
        raise Tier3ValidationError("annotation assembly does not match the H1 target assembly")
    h1_path, h2_path = _path(h1_fasta), _path(h2_fasta)
    if h1_path == h2_path:
        raise Tier3ValidationError("self alignment is forbidden: H1 and H2 paths are identical")
    target, query = read_fasta(h1_path), read_fasta(h2_path)
    if annotation_provenance:
        _audit_annotation_provenance(
            annotation_provenance,
            h1_path,
            target,
            h1_assembly_accession,
            annotation_gff,
        )
    wfmash_root = _validate_store_root(wfmash_store_path, "wfmash")
    bcftools_root = _validate_store_root(bcftools_store_path, "bcftools")
    output = _path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    raw_paf = output / "raw.wfmash.paf"
    command: Optional[List[str]] = None
    if paf_path is None:
        command = run_pinned_wfmash(h1_path, h2_path, raw_paf, wfmash_root, threads=threads)
    else:
        source_paf = _path(paf_path)
        if not source_paf.is_file():
            raise Tier3ValidationError(f"PAF does not exist: {source_paf}")
        if source_paf != raw_paf:
            shutil.copyfile(source_paf, raw_paf)

    with raw_paf.open("r", encoding="utf-8") as paf_handle:
        result = traverse_paf(
            paf_handle,
            target,
            query,
            edge_exclusion_bp=edge_exclusion_bp,
            indel_flank_bp=indel_flank_bp,
        )
    callable_positions = set(result.callable_positions)
    exclusion_counts = collections.Counter(result.exclusion_counts)
    with raw_paf.open("r", encoding="utf-8") as paf_handle:
        nonunique_query_targets = _nonunique_query_targets(paf_handle)
    if nonunique_query_targets:
        callable_positions.difference_update(nonunique_query_targets)
        exclusion_counts["multiple_query_projection"] += len(nonunique_query_targets)
    if h1_accessibility_bed is not None:
        accessible = _read_bed(h1_accessibility_bed, target)
        retained: Set[Tuple[str, int]] = set()
        for contig, position in callable_positions:
            if _position_in_intervals(position, accessible.get(contig, [])):
                retained.add((contig, position))
            else:
                exclusion_counts["h1_accessibility_mask"] += 1
        callable_positions = retained

    with raw_paf.open("r", encoding="utf-8") as paf_handle:
        variants = _paf_snv_alleles(paf_handle, target, query, callable_positions)
    if len(variants) != len(result.snv_positions & callable_positions):
        raise Tier3ValidationError("accepted X count does not agree with reconstructed H1/H2 SNV alleles")

    callable_bed = output / "h1.callable.bed"
    callable_bases = _positions_to_bed(callable_bed, callable_positions)
    normalized_bcf = output / "h1_vs_h2.normalized.snvs.bcf"
    normalized_index = output / "h1_vs_h2.normalized.snvs.bcf.csi"
    _write_snv_bcf(normalized_bcf, normalized_index, target, variants, sample, bcftools_root)

    with raw_paf.open("r", encoding="utf-8") as paf_handle:
        mapping_records = sum(
            1 for line in paf_handle if line.strip() and not line.startswith("#")
        )
    qc_path = output / "accepted_mapping_qc.json"
    qc = {
        "status": "eligible" if callable_bases else "unavailable_no_callable_alignment",
        "policy_version": "tier3-decisions-v1",
        "h1_role": "target_reference",
        "h2_role": "query_alternate_haplotype",
        "mapping_records": mapping_records,
        "phase_identity_audit_passed": True,
        "collapse_qc_passed": True,
        "edge_exclusion_bp": edge_exclusion_bp,
        "indel_flank_bp": indel_flank_bp,
        "operation_counts": result.operation_counts,
        "exclusion_counts": dict(sorted(exclusion_counts.items())),
        "callable_bases": callable_bases,
        "heterozygous_snvs": len(variants),
        "raw_paf_sha256": sha256_file(raw_paf),
        "wfmash_command": command,
    }
    _write_json(qc_path, qc)
    provenance_path = output / "provenance.json"
    provenance: Dict[str, Any] = {
        "status": qc["status"],
        "modality": "direct_wfmash_extended_cigar",
        "statistic_label": "individual_snv_heterozygosity",
        "population_pi": None,
        "assemblies": {
            "h1_target": h1_assembly_accession,
            "h2_query": h2_assembly_accession,
            "annotation_h1": (
                annotation_assembly_accession
                or (annotation_provenance or {}).get("assembly_accession_version")
            ),
        },
        "inputs": {
            "h1_fasta": _artifact(h1_path),
            "h2_fasta": _artifact(h2_path),
            "raw_wfmash_paf": _artifact(raw_paf),
        },
        "annotation_provenance": dict(annotation_provenance) if annotation_provenance else None,
        "tools": {
            "wfmash_commit": WFMASH_COMMIT,
            "wfmash_store_path": str(wfmash_root),
            "bcftools_store_path": str(bcftools_root),
        },
        "outputs": {
            "raw_paf": str(raw_paf),
            "accepted_mapping_qc": str(qc_path),
            "callable_bed": str(callable_bed),
            "normalized_snv_bcf": str(normalized_bcf),
            "normalized_snv_bcf_index": str(normalized_index),
        },
    }
    _write_json(provenance_path, provenance)
    record = {
        **provenance,
        **qc,
        "variant_alleles": variants,
        "provenance_path": str(provenance_path),
        "artifacts": {
            "raw_paf": _artifact(raw_paf),
            "accepted_mapping_qc": _artifact(qc_path),
            "callable_bed": _artifact(callable_bed),
            "normalized_snv_bcf": _artifact(normalized_bcf),
            "normalized_snv_bcf_index": _artifact(normalized_index),
        },
    }
    return record


def impg_core_windows(
    contig_lengths: Mapping[str, int],
    *,
    core_size_bp: int = DEFAULT_IMPG_CORE_BP,
    padding_bp: int = DEFAULT_IMPG_PADDING_BP,
) -> List[ImpgCoreWindow]:
    """Return disjoint ownership cores and clipped padded query windows."""

    if core_size_bp <= 0 or padding_bp < 0:
        raise Tier3ValidationError("IMPG core size must be positive and padding non-negative")
    windows: List[ImpgCoreWindow] = []
    for contig, length in contig_lengths.items():
        if length <= 0:
            raise Tier3ValidationError(f"invalid contig length for IMPG windows: {contig}={length}")
        for core_start in range(0, length, core_size_bp):
            core_end = min(length, core_start + core_size_bp)
            windows.append(
                ImpgCoreWindow(
                    contig,
                    core_start,
                    core_end,
                    max(0, core_start - padding_bp),
                    min(length, core_end + padding_bp),
                )
            )
    return windows


def own_normalized_records(
    records: Iterable[Mapping[str, Any]], windows: Sequence[ImpgCoreWindow]
) -> List[Dict[str, Any]]:
    """Deduplicate normalized alleles and assign POS-1 to one ownership core."""

    by_contig: Dict[str, List[ImpgCoreWindow]] = collections.defaultdict(list)
    for window in windows:
        by_contig[window.contig].append(window)
    for values in by_contig.values():
        values.sort(key=lambda item: item.core_start)
        for left, right in zip(values, values[1:]):
            if left.core_end != right.core_start:
                raise Tier3ValidationError("IMPG ownership cores are not disjoint and contiguous")

    unique: Dict[Tuple[str, int, str, Tuple[str, ...]], Mapping[str, Any]] = {}
    for record in records:
        contig, pos, ref = str(record["contig"]), int(record["pos"]), str(record["ref"])
        alt_value = record.get("alt", record.get("alts"))
        alts = (str(alt_value),) if isinstance(alt_value, str) else tuple(str(item) for item in alt_value)
        if pos < 1 or not ref or not alts:
            raise Tier3ValidationError("invalid normalized IMPG record")
        unique.setdefault((contig, pos, ref, alts), record)

    owned: List[Dict[str, Any]] = []
    for key, record in sorted(unique.items()):
        contig, pos, _ref, _alts = key
        anchor = pos - 1
        owners = [
            window for window in by_contig.get(contig, [])
            if window.core_start <= anchor < window.core_end
        ]
        if len(owners) != 1:
            raise Tier3ValidationError(
                f"normalized IMPG allele {contig}:{pos} has {len(owners)} ownership cores"
            )
        value = dict(record)
        value["owner_core_start"] = owners[0].core_start
        value["owner_core_end"] = owners[0].core_end
        owned.append(value)
    return owned


def audit_phase_orientation(
    expected: Mapping[Tuple[str, int, str, str], Sequence[Optional[int]]],
    observed: Mapping[Tuple[str, int, str, str], Sequence[Optional[int]]],
) -> Dict[str, Any]:
    """Audit IMPG phase orientation rather than trusting PanSN path order."""

    if not expected:
        return {
            "result": "failed",
            "audited_sites": 0,
            "matching_sites": 0,
            "inverted_sites": 0,
            "phase_sensitive_eligible": False,
        }
    matching = inverted = missing = 0
    for key, expected_gt in expected.items():
        observed_gt = observed.get(key)
        expected_tuple = tuple(expected_gt)
        if observed_gt is None:
            missing += 1
        elif tuple(observed_gt) == expected_tuple:
            matching += 1
        elif tuple(observed_gt) == tuple(reversed(expected_tuple)):
            inverted += 1
        else:
            missing += 1
    if matching == len(expected):
        result = "orientation_matches"
    elif inverted == len(expected):
        result = "orientation_inverted"
    else:
        result = "failed"
    # Frozen v1 forbids phase-sensitive IMPG use even when an audit happens to
    # match; the validated build inverted orientation in the truth test.
    return {
        "result": result,
        "audited_sites": len(expected),
        "matching_sites": matching,
        "inverted_sites": inverted,
        "missing_or_discordant_sites": missing,
        "phase_sensitive_eligible": False,
    }


def build_impg_query_plan(
    contig_lengths: Mapping[str, int],
    *,
    impg_store_path: Union[str, Path],
    callable_bed: Union[str, Path],
    truth_test_passed: bool,
    lace_threads: int = 2,
) -> Dict[str, Any]:
    """Build, but do not execute, the frozen optional IMPG concordance plan."""

    if not truth_test_passed:
        raise Tier3ValidationError("IMPG remains disabled until its pinned Guix truth test passes")
    if lace_threads < 2:
        raise Tier3ValidationError("IMPG lace requires at least two threads")
    root = _validate_store_root(impg_store_path, "impg")
    if not Path(callable_bed).is_file():
        raise Tier3ValidationError("IMPG cannot run without the direct PAF callable BED")
    windows = impg_core_windows(contig_lengths)
    return {
        "role": "optional_orthogonal_concordance",
        "impg_commit": IMPG_COMMIT,
        "impg_store_path": str(root),
        "windows": [asdict(window) for window in windows],
        "query_output": "vcf:poa",
        "merge_distance_bp": 0,
        "lace_threads": lace_threads,
        "postprocess": [
            "lace",
            "global_fasta_based_bcftools_norm",
            "exact_allele_deduplication",
            "bcf_conversion",
            "csi_index",
            "normalized_pos_minus_1_core_ownership",
            "phase_orientation_audit",
        ],
        "command_requirements": {
            "query": "impg query --merge-distance 0 -o vcf:poa <padded-region>",
            "lace": f"impg lace -t {lace_threads}",
            "normalization_scope": "global_after_all_padded_windows_before_core_ownership",
            "index": "bcftools index --csi",
        },
        "callable_denominator": str(Path(callable_bed).resolve()),
        "impg_supplies_denominator": False,
        "phase_sensitive_use": False,
    }


def _cli() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    alignment = subparsers.add_parser("alignment", help="collect direct H2-to-H1 WFMASH artifacts")
    alignment.add_argument("--h1", type=Path, required=True)
    alignment.add_argument("--h2", type=Path, required=True)
    alignment.add_argument("--paf", type=Path)
    alignment.add_argument("--output-dir", type=Path, required=True)
    alignment.add_argument("--sample", required=True)
    alignment.add_argument("--wfmash-store-path", required=True)
    alignment.add_argument("--bcftools-store-path", required=True)
    alignment.add_argument("--phase-qc-passed", action="store_true")
    alignment.add_argument("--collapse-qc-passed", action="store_true")
    alignment.add_argument("--h1-assembly-accession")
    alignment.add_argument("--h2-assembly-accession")
    alignment.add_argument("--annotation-assembly-accession")
    alignment.add_argument("--annotation-gff", type=Path)
    alignment.add_argument("--annotation-provenance-json", type=Path)
    alignment.add_argument("--h1-accessibility-bed", type=Path)
    alignment.add_argument("--threads", type=int, default=2)
    deposited = subparsers.add_parser("deposited", help="collect deposited exact-reference variants")
    deposited.add_argument("--reference", type=Path, required=True)
    deposited.add_argument("--variants", type=Path, required=True)
    deposited.add_argument("--callable-bed", type=Path)
    deposited.add_argument("--output-dir", type=Path, required=True)
    deposited.add_argument("--sample", required=True)
    deposited.add_argument("--bcftools-store-path", required=True)
    deposited.add_argument("--reference-assembly-accession", required=True)
    deposited.add_argument("--variant-reference-accession", required=True)
    deposited.add_argument("--annotation-gff", type=Path)
    deposited.add_argument("--annotation-provenance-json", type=Path)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _cli().parse_args(argv)
    if args.command == "alignment":
        annotation_provenance = (
            json.loads(args.annotation_provenance_json.read_text(encoding="utf-8"))
            if args.annotation_provenance_json
            else None
        )
        result = collect_direct_alignment(
            args.h1,
            args.h2,
            args.paf,
            args.output_dir,
            args.sample,
            args.wfmash_store_path,
            args.bcftools_store_path,
            args.phase_qc_passed,
            args.collapse_qc_passed,
            h1_accessibility_bed=args.h1_accessibility_bed,
            h1_assembly_accession=args.h1_assembly_accession,
            h2_assembly_accession=args.h2_assembly_accession,
            annotation_assembly_accession=args.annotation_assembly_accession,
            annotation_provenance=annotation_provenance,
            annotation_gff=args.annotation_gff,
            threads=args.threads,
        )
        print(json.dumps(result, indent=2, sort_keys=True))
    elif args.command == "deposited":
        annotation_provenance = (
            json.loads(args.annotation_provenance_json.read_text(encoding="utf-8"))
            if args.annotation_provenance_json
            else None
        )
        result = collect_deposited_variants(
            args.reference,
            args.variants,
            args.callable_bed,
            args.output_dir,
            args.sample,
            args.bcftools_store_path,
            reference_assembly_accession=args.reference_assembly_accession,
            variant_reference_accession=args.variant_reference_accession,
            annotation_provenance=annotation_provenance,
            annotation_gff=args.annotation_gff,
        )
        print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
