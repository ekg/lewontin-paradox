#!/usr/bin/env bash
# Reproduce the Tier 3A public-data acquisition in the frozen GNU Guix profile.
# Raw biological inputs are deliberately staged outside Git on MooseFS.
#
# SweepGA was built from clean code at $SWEEPGA_COMMIT (only documentation was
# modified in the source checkout) inside the pinned Guix shell described by
# acquisition_toolchain_manifest.scm.  The recorded build used the Guix GCC,
# Clang, GSL, jemalloc, HTSlib, and compression libraries and an explicit
# Rust/Cargo 1.94.1 front-end wrapper.  The complete literal build environment
# was:
#   guix shell -m results/tier3a/acquisition_toolchain_manifest.scm --pure \
#     --preserve=HOME --preserve=CARGO_HOME -- bash -c \
#     'PATH=/moosefs/erikg/impg/.guix/bin:$PATH CC=gcc CXX=g++ \
#      CARGO_TARGET_X86_64_UNKNOWN_LINUX_GNU_LINKER=gcc \
#      LIBCLANG_PATH=$GUIX_ENVIRONMENT/lib cargo build --release \
#      --manifest-path /moosefs/erikg/sweepga/Cargo.toml'
# IMPG/gfaffix were built from the exact staged source commit under the same
# Guix development manifest; their immutable SHA-256 values are in the TSV.
set -euo pipefail

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
ROOT=$(CDPATH= cd -- "$SCRIPT_DIR/../.." && pwd)
ENV_RECORD="$ROOT/analysis/pilot_results/guix_environment.json"
GUIX_JOB="$ROOT/analysis/slurm/guix_job.sh"
DATA_ROOT=${TIER3A_DATA_ROOT:-/moosefs/erikg/tier3data/tier3a-acquisition-20260716}
SCRATCH_ROOT=${TIER3_SCRATCH_ROOT:-/moosefs/erikg/tier3scratch/acquire-tier3a}
PROFILE=/gnu/store/z9v2f6faha9cwjz0sm5iphhlzisgi077-profile
WFMASH_ROOT=/gnu/store/w9x6axr2w0hhvzzm10gzlp06jg07806d-wfmash-tier3-0.24.2-12.e040aa1
BCFTOOLS_ROOT=/gnu/store/rlb2gljax8lzmhhidbvbzp3al1ad1mww-bcftools-1.14
SWEEPGA_COMMIT=018e4ce49d2c125820e0ac50dc5feaa02d423683
IMPG_COMMIT=101df81eb28a809c8fac97d297acd9fcfbbfa048
HANDOFF_TOOLS="$DATA_ROOT/tools/sweepga-impg-guix-handoff-018e4ce-101df81"
SWEEPGA_BIN="$HANDOFF_TOOLS/sweepga"
IMPG_BIN="$HANDOFF_TOOLS/impg"
export TIER3_ROOT="$ROOT" TIER3_SCRATCH_ROOT="$SCRATCH_ROOT"

run_guix() {
    "$GUIX_JOB" "$ENV_RECORD" "$@"
}

record_versions() {
    {
        echo "guix=$(guix --version | sed -n '1p')"
        run_guix python3 --version
        run_guix samtools --version | head -n2
        run_guix bcftools --version | head -n2
        echo "wfmash_commit=e040aa10e87cab44ed5a4db005e784be62b0bd21"
        echo "wfmash_store_path=$WFMASH_ROOT"
        run_guix "$SWEEPGA_BIN" --version
        run_guix "$IMPG_BIN" --version
        sha256sum "$SWEEPGA_BIN" "$HANDOFF_TOOLS/wfmash" "$IMPG_BIN" "$HANDOFF_TOOLS/gfaffix"
        echo "sweepga_commit=$SWEEPGA_COMMIT"
        echo "impg_commit=$IMPG_COMMIT"
        echo "impg_embedded_sweepga_commit=ddd31d39b6a68fc972025b048076032341b66835"
        echo "profile_store_path=$PROFILE"
    } > "$DATA_ROOT/guix_tool_versions.txt"
}

install -d -m 700 "$DATA_ROOT" "$SCRATCH_ROOT"

case "${1:-all}" in
download|all)
run_guix python3 - "$DATA_ROOT" <<'PY'
import gzip
import hashlib
import json
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

root = Path(sys.argv[1]).resolve()
inventory_url = (
    "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/bioproject/"
    "PRJNA489243/dataset_report?page_size=1000"
)
candidates = [
    {
        "id": "spinachia_spinachia_SK-2024b",
        "h1": "GCA_048126635.1", "h1_name": "SpiSpi1_v1.hap1",
        "h2": "GCA_048127205.1", "h2_name": "SpiSpi1_v1.hap2",
    },
    {
        "id": "menidia_menidia_fMenMen1",
        "h1": "GCA_048628825.1", "h1_name": "ASM4862882v1",
        "h2": "GCA_048544195.1", "h2_name": "ASM4854419v1",
    },
    {
        "id": "tautogolabrus_adspersus_fTauAds1",
        "h1": "GCA_020745685.1", "h1_name": "fTauAds1.pri.cur",
        "h2": "GCA_020745675.1", "h2_name": "fTauAds1.alt.cur",
    },
]

def fetch(url: str, destination: Path) -> None:
    """Retry/resume an immutable public object and atomically promote it."""
    if destination.is_file() and destination.stat().st_size:
        return
    destination.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
    temporary = destination.with_suffix(destination.suffix + ".part")
    for attempt in range(1, 7):
        try:
            offset = temporary.stat().st_size if temporary.exists() else 0
            request = urllib.request.Request(url)
            if offset:
                request.add_header("Range", f"bytes={offset}-")
            with urllib.request.urlopen(request, timeout=180) as response:
                append = offset and response.status == 206
                mode = "ab" if append else "wb"
                with temporary.open(mode) as output:
                    while True:
                        block = response.read(8 * 1024 * 1024)
                        if not block:
                            break
                        output.write(block)
            if not temporary.stat().st_size:
                raise OSError("empty response")
            temporary.replace(destination)
            return
        except (OSError, urllib.error.URLError) as error:
            if attempt == 6:
                raise
            print(f"retry {attempt}/6 for {url}: {error}", file=sys.stderr)
            time.sleep(2 ** attempt)

def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()

def md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()

def gunzip(source: Path, destination: Path) -> None:
    if destination.is_file() and destination.stat().st_size:
        return
    temporary = destination.with_suffix(destination.suffix + ".part")
    with gzip.open(source, "rb") as inp, temporary.open("wb") as out:
        while True:
            block = inp.read(8 * 1024 * 1024)
            if not block:
                break
            out.write(block)
    temporary.replace(destination)

inventory_path = root / "authoritative_vgp_PRJNA489243_dataset_report.json"
inventory = json.loads(inventory_path.read_text()) if inventory_path.is_file() else {}
if len(inventory.get("reports", [])) < int(inventory.get("total_count", 1)):
    reports = []
    page_token = None
    page_urls = []
    total_count = None
    while True:
        parameters = {"page_size": 1000}
        if page_token:
            parameters["page_token"] = page_token
        url = (
            "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/bioproject/PRJNA489243/"
            "dataset_report?" + urllib.parse.urlencode(parameters)
        )
        page_urls.append(url)
        page = json.load(urllib.request.urlopen(url, timeout=180))
        reports.extend(page.get("reports", []))
        total_count = int(page.get("total_count", total_count or len(reports)))
        page_token = page.get("next_page_token")
        if not page_token:
            break
    if len(reports) != total_count:
        raise SystemExit(f"incomplete VGP inventory: {len(reports)}/{total_count}")
    temporary = inventory_path.with_suffix(".json.part")
    temporary.write_text(
        json.dumps(
            {"reports": reports, "total_count": total_count, "source_page_urls": page_urls},
            separators=(",", ":"),
        ) + "\n",
        encoding="utf-8",
    )
    temporary.replace(inventory_path)

for candidate in candidates:
    candidate_dir = root / candidate["id"]
    metadata_dir = candidate_dir / "metadata"
    metadata_dir.mkdir(parents=True, exist_ok=True, mode=0o700)
    for role in ("h1", "h2"):
        accession = candidate[role]
        assembly_name = candidate[f"{role}_name"]
        digits = accession.split("_")[1].split(".")[0]
        prefix = f"{accession}_{assembly_name}"
        ftp_root = (
            "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/"
            f"{digits[:3]}/{digits[3:6]}/{digits[6:9]}/{prefix}/"
        )
        fetch(
            f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/"
            f"{accession}/dataset_report",
            metadata_dir / f"{role}.dataset_report.json",
        )
        wanted = [
            "md5checksums.txt", "uncompressed_checksums.txt",
            f"{prefix}_assembly_report.txt", f"{prefix}_assembly_stats.txt",
            f"{prefix}_genomic.fna.gz",
        ]
        if role == "h1":
            wanted.append(f"{prefix}_genomic.gff.gz")
        for filename in wanted:
            try:
                fetch(ftp_root + filename, metadata_dir / filename if filename.endswith(".txt") else candidate_dir / filename)
            except urllib.error.HTTPError as error:
                # Older assemblies need not publish the newer uncompressed checksum ledger.
                if filename == "uncompressed_checksums.txt" and error.code == 404:
                    continue
                raise

        checksum_file = metadata_dir / "md5checksums.txt"
        # Both roles have their own file with the same basename; preserve after verification.
        role_checksum_file = metadata_dir / f"{role}.md5checksums.txt"
        if checksum_file.exists():
            checksum_file.replace(role_checksum_file)
        checksums = {}
        for line in role_checksum_file.read_text(encoding="utf-8").splitlines():
            fields = line.split(maxsplit=1)
            if len(fields) == 2:
                checksums[fields[1].removeprefix("./")] = fields[0]
        payloads = [candidate_dir / f"{prefix}_genomic.fna.gz"]
        if role == "h1":
            payloads.append(candidate_dir / f"{prefix}_genomic.gff.gz")
        for payload in payloads:
            expected = checksums.get(payload.name)
            if not expected or md5(payload) != expected:
                raise SystemExit(f"NCBI MD5 verification failed for {payload}")

        fasta_gz = candidate_dir / f"{prefix}_genomic.fna.gz"
        fasta = candidate_dir / f"{role}.fna"
        gunzip(fasta_gz, fasta)
        if role == "h1":
            gunzip(candidate_dir / f"{prefix}_genomic.gff.gz", candidate_dir / "h1.native.gff")

    (metadata_dir / "selection.json").write_text(
        json.dumps(candidate, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
PY

for candidate_dir in "$DATA_ROOT"/*_*_*; do
    [[ -d "$candidate_dir" ]] || continue
    run_guix samtools faidx "$candidate_dir/h1.fna"
    run_guix samtools faidx "$candidate_dir/h2.fna"
    run_guix python3 - "$candidate_dir" <<'PY'
import hashlib
import json
import re
import sys
from pathlib import Path

root = Path(sys.argv[1])
fai = {}
for line in (root / "h1.fna.fai").read_text(encoding="utf-8").splitlines():
    name, length, *_ = line.split("\t")
    fai[name] = int(length)
regions = {}
with (root / "h1.native.gff").open(encoding="utf-8") as handle:
    for line in handle:
        match = re.match(r"##sequence-region\s+(\S+)\s+1\s+(\d+)$", line.rstrip())
        if match:
            regions[match.group(1)] = int(match.group(2))
if not regions:
    raise SystemExit("native GFF has no ##sequence-region dictionary")
if regions != fai:
    raise SystemExit("native GFF sequence-region dictionary does not exactly match H1 FASTA")
mapping = root / "h1.annotation_contig_map.tsv"
with mapping.open("w", encoding="utf-8") as out:
    out.write("annotation_contig\tfasta_contig\tlength\n")
    for name in sorted(regions):
        out.write(f"{name}\t{name}\t{regions[name]}\n")
PY
done
;;
esac

case "${1:-all}" in
smoke|all)
for candidate_dir in "$DATA_ROOT"/*_*_*; do
    [[ -f "$candidate_dir/h1.fna.fai" ]] || continue
    smoke="$candidate_dir/wfmash_smoke"
    install -d -m 700 "$smoke"
    h1_contig=$(awk -F '\t' '$2 > n {n=$2; name=$1} END {print name}' "$candidate_dir/h1.fna.fai")
    h2_contig=$(awk -F '\t' '$2 > n {n=$2; name=$1} END {print name}' "$candidate_dir/h2.fna.fai")
    run_guix samtools faidx "$candidate_dir/h1.fna" "${h1_contig}:1-250000" > "$smoke/h1.250kb.fna"
    run_guix samtools faidx "$candidate_dir/h2.fna" "${h2_contig}:1-250000" > "$smoke/h2.250kb.fna"
    run_guix samtools faidx "$smoke/h1.250kb.fna"
    run_guix samtools faidx "$smoke/h2.250kb.fna"
    # The older fTauAds1 alternate is scaffold-level, so its longest scaffold
    # is not necessarily homologous to the longest H1 chromosome. Locate an
    # H2 segment with the same pinned WFMASH policy, then freeze that region.
    if [[ $(basename "$candidate_dir") == tautogolabrus_adspersus_fTauAds1 ]]; then
        run_guix "$WFMASH_ROOT/bin/wfmash" "$smoke/h1.250kb.fna" "$candidate_dir/h2.fna" \
            -p 90 -w 5k -l 25k -o -4 -t 2 > "$smoke/h1_region_vs_full_h2.paf"
        read -r qname qstart qend tstart tend < <(
            run_guix python3 - "$smoke/h1_region_vs_full_h2.paf" <<'PY'
import sys
rows = [line.rstrip().split("\t") for line in open(sys.argv[1]) if line.strip()]
if not rows:
    raise SystemExit("no homologous H2 smoke region located")
row = max(rows, key=lambda x: (int(x[11]), int(x[3]) - int(x[2])))
print(row[0], row[2], row[3], row[7], row[8])
PY
        )
        run_guix samtools faidx "$candidate_dir/h1.fna" \
            "${h1_contig}:$((tstart + 1))-${tend}" > "$smoke/h1.250kb.fna"
        run_guix samtools faidx "$candidate_dir/h2.fna" \
            "${qname}:$((qstart + 1))-${qend}" > "$smoke/h2.250kb.fna"
        run_guix python3 -c "from pathlib import Path; [p.unlink(missing_ok=True) for p in map(Path, ['$smoke/h1.250kb.fna.fai', '$smoke/h2.250kb.fna.fai'])]"
    fi
    run_guix samtools faidx "$smoke/h1.250kb.fna"
    run_guix samtools faidx "$smoke/h2.250kb.fna"
    accession_h1=$(run_guix python3 -c "import json; print(json.load(open('$candidate_dir/metadata/selection.json'))['h1'])")
    accession_h2=$(run_guix python3 -c "import json; print(json.load(open('$candidate_dir/metadata/selection.json'))['h2'])")
    run_guix python3 "$ROOT/analysis/tier3a_vgp_collect.py" alignment \
        --h1 "$smoke/h1.250kb.fna" --h2 "$smoke/h2.250kb.fna" \
        --output-dir "$smoke/parser" --sample "$(basename "$candidate_dir")" \
        --wfmash-store-path "$WFMASH_ROOT" --bcftools-store-path "$BCFTOOLS_ROOT" \
        --phase-qc-passed --collapse-qc-passed \
        --h1-assembly-accession "$accession_h1" --h2-assembly-accession "$accession_h2" \
        --threads 2 > "$smoke/parser_result.json"
    run_guix python3 "$ROOT/results/tier3a/acquisition_build_mask.py" \
        --h1 "$smoke/h1.250kb.fna" --h2 "$smoke/h2.250kb.fna" \
        --paf "$smoke/parser/raw.wfmash.paf" --output-dir "$smoke/interval_parser" \
        --sample "$(basename "$candidate_dir")" --bcftools "$BCFTOOLS_ROOT/bin/bcftools" \
        > "$smoke/interval_parser_result.json"
    run_guix python3 - "$smoke/parser/accepted_mapping_qc.json" "$smoke/interval_parser/accepted_mapping_qc.json" <<'PY'
import json
import sys
left, right = (json.load(open(path)) for path in sys.argv[1:])
for field in ("callable_bases", "heterozygous_snvs", "operation_counts"):
    if left[field] != right[field]:
        raise SystemExit(f"interval/per-base smoke parser disagreement for {field}: {left[field]} != {right[field]}")
PY
done
;;
esac

case "${1:-all}" in
align|all)
for candidate_dir in "$DATA_ROOT"/*_*_*; do
    [[ -f "$candidate_dir/h1.fna.fai" ]] || continue
    alignment="$candidate_dir/alignment"
    install -d -m 700 "$alignment"
    if [[ ! -s "$alignment/raw.wfmash.paf" ]]; then
        temporary="$alignment/raw.wfmash.paf.part"
        run_guix "$WFMASH_ROOT/bin/wfmash" "$candidate_dir/h1.fna" "$candidate_dir/h2.fna" \
            -p 90 -w 5k -l 25k -o -4 -t "${TIER3A_THREADS:-8}" \
            > "$temporary" 2> "$alignment/wfmash.stderr.log"
        [[ -s "$temporary" ]] || { echo "empty full WFMASH PAF for $candidate_dir" >&2; exit 65; }
        mv "$temporary" "$alignment/raw.wfmash.paf"
    fi
done
;;
esac

case "${1:-all}" in
mask|all)
for candidate_dir in "$DATA_ROOT"/*_*_*; do
    [[ -s "$candidate_dir/alignment/raw.wfmash.paf" ]] || continue
    if [[ ! -s "$candidate_dir/alignment/accepted_mapping_qc.json" ]]; then
        run_guix python3 "$ROOT/results/tier3a/acquisition_build_mask.py" \
            --h1 "$candidate_dir/h1.fna" --h2 "$candidate_dir/h2.fna" \
            --paf "$candidate_dir/alignment/raw.wfmash.paf" \
            --output-dir "$candidate_dir/alignment" --sample "$(basename "$candidate_dir")" \
            --bcftools "$BCFTOOLS_ROOT/bin/bcftools" \
            > "$candidate_dir/alignment/parser_result.json"
    fi
done
;;
esac

case "${1:-all}" in
bounded|all)
# Retained sensitivity experiment: SweepGA filters the earlier standalone
# WFMASH PAF at 1:1, 5:5, and 10:10.  This is not the production mapping path.
# "M:N" sets query/target retained-hit cardinalities under the explicit 0.95
# rejection-overlap policy; raw interval concurrency is not that thresholded
# count. Sequence keys are coordinate domains, not annotation partitions.
[[ -x "$SWEEPGA_BIN" ]] || { echo "missing Guix-built SweepGA: $SWEEPGA_BIN" >&2; exit 66; }
for candidate_dir in "$DATA_ROOT"/*_*_*; do
    [[ -s "$candidate_dir/alignment/raw.wfmash.paf" ]] || continue
    bounded="$candidate_dir/sweepga"
    install -d -m 700 "$bounded"
    for cap in 1 5 10; do
        run_guix "$SWEEPGA_BIN" "$candidate_dir/alignment/raw.wfmash.paf" \
            --output-file "$bounded/cap${cap}.paf" \
            --num-mappings "${cap}:${cap}" --scaffold-jump 0 \
            --overlap 0.95 --scoring log-length-ani --threads 2 \
            > "$bounded/cap${cap}.stdout.log" 2> "$bounded/cap${cap}.stderr.log"
        [[ -s "$bounded/cap${cap}.paf" ]] || { echo "empty cap ${cap}:${cap} PAF" >&2; exit 65; }
    done
done
run_guix python3 "$ROOT/results/tier3a/acquisition_enrich_manifest.py" audit-caps \
    --data-root "$DATA_ROOT"
for candidate_dir in "$DATA_ROOT"/*_*_*; do
    [[ -s "$candidate_dir/sweepga/selected_cap.txt" ]] || continue
    selected=$(cut -d: -f1 "$candidate_dir/sweepga/selected_cap.txt")
    primary="$candidate_dir/sweepga/primary"
    rm -rf "$primary"
    run_guix python3 "$ROOT/results/tier3a/acquisition_build_mask.py" \
        --h1 "$candidate_dir/h1.fna" --h2 "$candidate_dir/h2.fna" \
        --paf "$candidate_dir/sweepga/cap${selected}.paf" \
        --output-dir "$primary" --sample "$(basename "$candidate_dir")" \
        --bcftools "$BCFTOOLS_ROOT/bin/bcftools" \
        > "$candidate_dir/sweepga/primary_result.json"
done
;;
esac

case "${1:-}" in
legacy-impg-truth)
# Historical region-direct sensitivity check. It is deliberately excluded from
# `all` and is not the production partition/query workflow documented below.
# IMPG's graph VCF uses region-local contig names, so truth checks are kept in
# two non-overlapping 10 kb H1 ownership windows.  The region FASTA headers are
# rewritten to those exact 0-based half-open names before normalization.
candidate_dir="$DATA_ROOT/spinachia_spinachia_SK-2024b"
truth="$candidate_dir/impg_truth"
install -d -m 700 "$truth/regions" "$truth/tmp"
printf '%s\n' 'CM106590.1:50000-60000' 'CM106590.1:105000-115000' > "$truth/regions.txt"
: > "$truth/h1.regions.fna"
i=0
while IFS=: read -r contig span; do
    i=$((i + 1)); start=${span%-*}; end=${span#*-}; one=$((start + 1))
    name="${contig}:${start}-${end}"
    run_guix samtools faidx "$candidate_dir/h1.fna" "${contig}:${one}-${end}" | \
        run_guix python3 -c \
        'import sys; x=sys.stdin.readlines(); x[0]=">"+sys.argv[1]+"\n"; sys.stdout.writelines(x)' \
        "$name" >> "$truth/h1.regions.fna"
    run_guix "$IMPG_BIN" query -a "$candidate_dir/sweepga/cap1.paf" \
        -r "$name" --no-merge -o vcf:poa \
        --sequence-files "$candidate_dir/h1.fna" "$candidate_dir/h2.fna" \
        -t 2 > "$truth/regions/region${i}.vcf" \
        2> "$truth/regions/region${i}.stderr.log"
    [[ $(grep -vc '^#' "$truth/regions/region${i}.vcf") -gt 0 ]] || \
        { echo "header-only IMPG truth VCF for $name" >&2; exit 65; }
done < "$truth/regions.txt"
run_guix samtools faidx "$truth/h1.regions.fna"
printf '%s\n' "$truth/regions/region1.vcf" "$truth/regions/region2.vcf" > "$truth/vcf.list"
run_guix "$IMPG_BIN" lace --file-list "$truth/vcf.list" --format vcf \
    --reference "$candidate_dir/h1.fna" --output "$truth/laced.vcf" \
    --compress none --temp-dir "$truth/tmp" -t 2 \
    > "$truth/lace.stdout.log" 2> "$truth/lace.stderr.log"
run_guix bcftools norm -f "$candidate_dir/h1.fna" -m -any "$truth/laced.vcf" -Ou | \
    run_guix bcftools sort -Ou -T "$truth/bcftools-sort" | \
    run_guix bcftools norm -d exact -Ob -o "$truth/laced.normalized.bcf"
run_guix bcftools index -f "$truth/laced.normalized.bcf"
[[ $(run_guix bcftools view -H "$truth/laced.normalized.bcf" | wc -l) -gt 0 ]] || \
    { echo "header-only normalized IMPG truth BCF" >&2; exit 65; }
;;
esac

case "${1:-all}" in
direct-sweepga|all)
# Production mapping. SweepGA reads both complete FASTAs and invokes its
# documented WFMASH backend because the ranked FastGA attempt failed at
# GIXmake's whole-assembly index-size limit. SweepGA applies the 1:1 cap in
# the same command. There is no annotation-driven preprocessing partition.
for candidate_dir in "$DATA_ROOT"/*_*_*; do
    [[ -s "$candidate_dir/h1.fna" && -s "$candidate_dir/h2.fna" ]] || continue
    output="$candidate_dir/sweepga/direct"
    scratch="$SCRATCH_ROOT/direct-$(basename "$candidate_dir")"
    install -d -m 700 "$output" "$scratch"
    [[ -s "$output/cap1.paf" ]] && continue
    temporary="$output/cap1.paf.part"
    run_guix bash -c \
      'set -euo pipefail; export PATH="$1:$PATH" WFMASH_BIN_DIR="$1"; exec "$1/sweepga" "$2" "$3" --output-file "$4" --aligner wfmash --map-pct-identity 90 --min-aln-length 25k --num-mappings 1:1 --scaffold-jump 0 --overlap 0.95 --scoring log-length-ani --threads "$5" --temp-dir "$6"' \
      bash "$HANDOFF_TOOLS" "$candidate_dir/h1.fna" "$candidate_dir/h2.fna" \
      "$temporary" "${TIER3A_THREADS:-8}" "$scratch" \
      > "$output/cap1.stdout.log" 2> "$output/cap1.stderr.log"
    [[ -s "$temporary" ]] || { echo "empty direct SweepGA PAF for $candidate_dir" >&2; exit 65; }
    mv "$temporary" "$output/cap1.paf"
done
;;
esac

case "${1:-all}" in
cap-recheck|all)
# Independently reapply the exact pinned SweepGA 1:1/0.95 policy to its own
# output. The 12 mandatory PAF columns must be an ordered fixed point. SweepGA
# appends an st tag when reading PAF, so optional-column byte identity is not
# the correct invariant. Raw interval concurrency is reported separately and
# is not the cap definition when the documented overlap threshold is 0.95.
for candidate_dir in "$DATA_ROOT"/*_*_*; do
    input="$candidate_dir/sweepga/direct/cap1.paf"
    [[ -s "$input" ]] || continue
    output="$candidate_dir/sweepga/direct/cap1.recheck.paf"
    temporary="$output.part"
    run_guix "$SWEEPGA_BIN" "$input" --output-file "$temporary" \
      --num-mappings 1:1 --scaffold-jump 0 --overlap 0.95 \
      --scoring log-length-ani --threads 2 \
      > "$candidate_dir/sweepga/direct/cap1.recheck.stdout.log" \
      2> "$candidate_dir/sweepga/direct/cap1.recheck.stderr.log"
    run_guix python3 - "$input" "$temporary" <<'PY'
import sys
def core(path):
    with open(path, encoding="utf-8") as handle:
        return [tuple(line.rstrip().split("\t")[:12]) for line in handle if line.strip()]
left, right = map(core, sys.argv[1:])
if not left or left != right:
    raise SystemExit("SweepGA 1:1 exact-policy mandatory PAF fields are not a fixed point")
PY
    mv "$temporary" "$output"
done
;;
esac

case "${1:-all}" in
queries|all)
# Compile provider-native H1 coding features only after whole-haplotype mapping.
# Overlap-merging is an execution optimization; feature identity and CDS phase
# remain in the deterministic query manifest and span-to-feature map.
for candidate_dir in "$DATA_ROOT"/*_*_*; do
    [[ -s "$candidate_dir/sweepga/direct/cap1.paf" ]] || continue
    release=$(run_guix python3 - "$candidate_dir/metadata/h1.dataset_report.json" <<'PY'
import json, sys
print(json.load(open(sys.argv[1]))["reports"][0]["annotation_info"]["name"])
PY
)
    run_guix python3 "$ROOT/results/tier3a/acquisition_build_queries.py" \
      --dataset-id "$(basename "$candidate_dir")" --gff "$candidate_dir/h1.native.gff" \
      --fai "$candidate_dir/h1.fna.fai" --h2-fai "$candidate_dir/h2.fna.fai" \
      --bounded-paf "$candidate_dir/sweepga/direct/cap1.paf" \
      --output-dir "$candidate_dir/annotation_query" --annotation-release "$release" \
      > "$candidate_dir/annotation_query.build.json"
done
;;
esac

case "${1:-}" in
impg-production)
# Run only on downstream compute: IMPG, not SweepGA, partitions the implicit
# graph and focuses it on annotation-defined regions. `query` emits VCF;
# bcftools performs normalization, target trimming, BCF encoding, and indexing.
for candidate_dir in "$DATA_ROOT"/*_*_*; do
    [[ -s "$candidate_dir/annotation_query/impg_execution_spans.bed" ]] || continue
    work="$candidate_dir/impg_production"
    install -d -m 700 "$work/partitions" "$work/calls" "$work/tmp"
    run_guix "$IMPG_BIN" index -a "$candidate_dir/sweepga/direct/cap1.paf" \
      -i "$work/h1_h2.impg" -t "${TIER3A_THREADS:-8}"
    run_guix "$IMPG_BIN" partition -a "$candidate_dir/sweepga/direct/cap1.paf" \
      -i "$work/h1_h2.impg" -w 2000 -d 0 --min-missing-size 1 \
      --min-boundary-distance 0 -o bed --output-folder "$work/partitions" \
      -t "${TIER3A_THREADS:-8}"
    run_guix python3 "$ROOT/results/tier3a/acquisition_select_impg_partitions.py" \
      --partitions "$work/partitions/partitions.bed" \
      --targets "$candidate_dir/annotation_query/impg_execution_spans.bed" \
      --focus-bed "$work/focus.native_partitions.bed" \
      --mapping-tsv "$work/partition_annotation_map.tsv"
    run_guix "$IMPG_BIN" query -a "$candidate_dir/sweepga/direct/cap1.paf" \
      -i "$work/h1_h2.impg" -b "$work/focus.native_partitions.bed" -d 0 \
      -o vcf:poa --sequence-files "$candidate_dir/h1.fna" "$candidate_dir/h2.fna" \
      -O "$work/calls" -t "${TIER3A_THREADS:-8}"
    find "$work/calls" -name '*.vcf' -type f | sort > "$work/vcf.list"
    [[ -s "$work/vcf.list" ]] || { echo "IMPG emitted no regional VCFs" >&2; exit 65; }
    run_guix "$IMPG_BIN" lace -l "$work/vcf.list" --format vcf -o "$work/laced.vcf" \
      --reference "$candidate_dir/h1.fna" --compress none --temp-dir "$work/tmp" \
      -t "${TIER3A_THREADS:-8}"
    run_guix bcftools norm -f "$candidate_dir/h1.fna" -m -any "$work/laced.vcf" -Ou | \
      run_guix bcftools sort -Ou -T "$work/bcftools-sort" | \
      run_guix bcftools view -R "$candidate_dir/annotation_query/impg_execution_spans.bed" -Ou | \
      run_guix bcftools norm -d exact -Ob -o "$work/normalized.targeted.bcf"
    run_guix bcftools index -f "$work/normalized.targeted.bcf"
    run_guix bcftools view -Oz -o "$work/normalized.targeted.vcf.gz" "$work/normalized.targeted.bcf"
    run_guix bcftools index -f -t "$work/normalized.targeted.vcf.gz"
done
;;
esac

case "${1:-all}" in
manifest|all)
record_versions
run_guix python3 - "$DATA_ROOT" "$ROOT/results/tier3a/acquisition_manifest.tsv" "$ROOT" <<'PY'
import csv
import hashlib
import json
import sys
from pathlib import Path

data_root, output_path, repo_root = map(Path, sys.argv[1:])
channel_commit = "44bbfc24e4bcc48d0e3343cd3d83452721af8c36"
profile = "/gnu/store/z9v2f6faha9cwjz0sm5iphhlzisgi077-profile"
wfmash_root = "/gnu/store/w9x6axr2w0hhvzzm10gzlp06jg07806d-wfmash-tier3-0.24.2-12.e040aa1"

def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()

def report(candidate, role):
    return json.loads((candidate / "metadata" / f"{role}.dataset_report.json").read_text())["reports"][0]

def provider_md5(candidate, role, basename):
    for line in (candidate / "metadata" / f"{role}.md5checksums.txt").read_text().splitlines():
        fields = line.split(maxsplit=1)
        if len(fields) == 2 and fields[1].removeprefix("./") == basename:
            return fields[0]
    raise SystemExit(f"provider MD5 absent for {basename}")

fields = [
    "dataset_id", "eligibility_status", "scientific_name", "taxon_id", "biosample_accession",
    "individual_id", "sex", "source_project", "assembly_release", "h1_accession_version",
    "h1_assembly_name", "h1_haplotype_label", "h1_fasta_path", "h1_fasta_sha256",
    "h1_fasta_size_bytes", "h1_compressed_path", "h1_compressed_sha256", "h1_provider_md5",
    "h1_fai_path", "h1_fai_sha256", "h1_contig_count", "h1_total_length", "h1_contig_n50",
    "h1_scaffold_n50", "h2_accession_version", "h2_assembly_name", "h2_haplotype_label",
    "h2_fasta_path", "h2_fasta_sha256", "h2_fasta_size_bytes", "h2_compressed_path",
    "h2_compressed_sha256", "h2_provider_md5", "h2_fai_path", "h2_fai_sha256",
    "h2_contig_count", "h2_total_length", "h2_contig_n50", "h2_scaffold_n50",
    "haplotype_span_ratio_h2_over_h1", "phase_identity_qc", "collapse_duplication_qc",
    "assembly_method", "sequencing_technology", "annotation_accession_version",
    "annotation_provider", "annotation_release_date", "annotation_software",
    "annotation_status", "annotation_genetic_code", "annotation_gff_path",
    "annotation_gff_sha256", "annotation_gff_size_bytes", "annotation_compressed_path",
    "annotation_compressed_sha256", "annotation_provider_md5", "annotation_provenance_path",
    "annotation_provenance_sha256", "annotation_contig_map_path", "annotation_contig_map_sha256",
    "annotation_contig_map_rows", "annotation_contig_audit", "modality",
    "deposited_variant_status", "raw_wfmash_paf_path", "raw_wfmash_paf_sha256",
    "callable_bed_path", "callable_bed_sha256", "callable_bases", "heterozygous_snvs",
    "normalized_bcf_path", "normalized_bcf_sha256", "normalized_bcf_csi_path",
    "normalized_bcf_csi_sha256", "callable_construction", "callable_exclusions_json",
    "smoke_qc_path", "smoke_qc_sha256", "smoke_callable_bases", "smoke_parser_status",
    "guix_manager", "guix_channel_commit", "guix_channels_path", "guix_channels_sha256",
    "guix_manifest_path", "guix_manifest_sha256", "guix_profile_store_path",
    "wfmash_store_path", "wfmash_commit", "bcftools_store_path", "commands_path",
    "authoritative_inventory_path", "authoritative_inventory_sha256",
]
rows = []
for candidate in sorted(path for path in data_root.iterdir() if path.is_dir() and (path / "metadata/selection.json").is_file()):
    selection = json.loads((candidate / "metadata/selection.json").read_text())
    h1_report, h2_report = report(candidate, "h1"), report(candidate, "h2")
    h1_info, h2_info = h1_report["assembly_info"], h2_report["assembly_info"]
    h1_stats, h2_stats = h1_report["assembly_stats"], h2_report["assembly_stats"]
    annotation = h1_report["annotation_info"]
    h1_prefix = f"{selection['h1']}_{selection['h1_name']}"
    h2_prefix = f"{selection['h2']}_{selection['h2_name']}"
    h1, h2 = candidate / "h1.fna", candidate / "h2.fna"
    h1_gz = candidate / f"{h1_prefix}_genomic.fna.gz"
    h2_gz = candidate / f"{h2_prefix}_genomic.fna.gz"
    gff, gff_gz = candidate / "h1.native.gff", candidate / f"{h1_prefix}_genomic.gff.gz"
    mapping = candidate / "h1.annotation_contig_map.tsv"
    provenance = {
        "provider": annotation["provider"],
        "release": annotation["name"],
        "release_date": annotation["release_date"],
        "assembly_accession_version": selection["h1"],
        "status": "native",
        "native_vs_projected": "native_exact_assembly_submitted_annotation",
        "genetic_code": 1,
        "fasta_sha256": sha256(h1),
        "gff_sha256": sha256(gff),
        "contig_mapping": {},
        "contig_mapping_path": str(mapping.resolve()),
        "contig_mapping_sha256": sha256(mapping),
        "contig_dictionary_audit_passed": True,
        "projected_lifted_congener_or_agent_denovo_used": False,
    }
    provenance_path = candidate / "h1.annotation_provenance.json"
    provenance_path.write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")
    qc_path = candidate / "alignment/accepted_mapping_qc.json"
    if not qc_path.is_file():
        raise SystemExit(f"whole-assembly QC missing: {qc_path}")
    qc = json.loads(qc_path.read_text())
    if qc.get("status") != "eligible" or int(qc.get("callable_bases", 0)) <= 0:
        raise SystemExit(f"noneligible whole-assembly QC: {qc_path}")
    smoke_path = candidate / "wfmash_smoke/parser/accepted_mapping_qc.json"
    smoke = json.loads(smoke_path.read_text())
    isolate = h1_report["organism"].get("infraspecific_names", {}).get("isolate", "")
    sex = h1_report["organism"].get("infraspecific_names", {}).get("sex", "")
    h1_role = h1_info["diploid_role"]
    h2_role = h2_info["diploid_role"]
    linked_h1 = {item["linked_assembly"] for item in h1_info.get("linked_assemblies", [])}
    linked_h2 = {item["linked_assembly"] for item in h2_info.get("linked_assemblies", [])}
    if selection["h2"] not in linked_h1 or selection["h1"] not in linked_h2:
        raise SystemExit(f"nonreciprocal H1/H2 links for {candidate.name}")
    h1_contigs = sum(1 for line in (candidate / "h1.fna.fai").read_text().splitlines() if line)
    h2_contigs = sum(1 for line in (candidate / "h2.fna.fai").read_text().splitlines() if line)
    map_rows = sum(1 for line in mapping.read_text().splitlines()[1:] if line)
    if h1_contigs != map_rows:
        raise SystemExit(f"annotation map is not an exact H1 dictionary for {candidate.name}")
    bcf = candidate / "alignment/h1_vs_h2.normalized.snvs.bcf"
    csi = Path(str(bcf) + ".csi")
    row = {
        "dataset_id": candidate.name, "eligibility_status": "eligible_biological",
        "scientific_name": h1_report["organism"]["organism_name"], "taxon_id": h1_report["organism"]["tax_id"],
        "biosample_accession": h1_info["biosample"]["accession"], "individual_id": isolate, "sex": sex,
        "source_project": "VGP PRJNA489243 lineage", "assembly_release": f"{h1_info['release_date']} / {h2_info['release_date']}",
        "h1_accession_version": selection["h1"], "h1_assembly_name": h1_info["assembly_name"], "h1_haplotype_label": h1_role,
        "h1_fasta_path": str(h1.resolve()), "h1_fasta_sha256": sha256(h1), "h1_fasta_size_bytes": h1.stat().st_size,
        "h1_compressed_path": str(h1_gz.resolve()), "h1_compressed_sha256": sha256(h1_gz), "h1_provider_md5": provider_md5(candidate, "h1", h1_gz.name),
        "h1_fai_path": str((candidate / "h1.fna.fai").resolve()), "h1_fai_sha256": sha256(candidate / "h1.fna.fai"), "h1_contig_count": h1_contigs,
        "h1_total_length": h1_stats["total_sequence_length"], "h1_contig_n50": h1_stats["contig_n50"], "h1_scaffold_n50": h1_stats["scaffold_n50"],
        "h2_accession_version": selection["h2"], "h2_assembly_name": h2_info["assembly_name"], "h2_haplotype_label": h2_role,
        "h2_fasta_path": str(h2.resolve()), "h2_fasta_sha256": sha256(h2), "h2_fasta_size_bytes": h2.stat().st_size,
        "h2_compressed_path": str(h2_gz.resolve()), "h2_compressed_sha256": sha256(h2_gz), "h2_provider_md5": provider_md5(candidate, "h2", h2_gz.name),
        "h2_fai_path": str((candidate / "h2.fna.fai").resolve()), "h2_fai_sha256": sha256(candidate / "h2.fna.fai"), "h2_contig_count": h2_contigs,
        "h2_total_length": h2_stats["total_sequence_length"], "h2_contig_n50": h2_stats["contig_n50"], "h2_scaffold_n50": h2_stats["scaffold_n50"],
        "haplotype_span_ratio_h2_over_h1": f"{int(h2_stats['total_sequence_length']) / int(h1_stats['total_sequence_length']):.9f}",
        "phase_identity_qc": "passed_reciprocal_NCBI_link_and_same_isolate",
        "collapse_duplication_qc": "passed_provider_purge_dups_and_comparable_haplotype_span;independent_kmer_switch_error_not_deposited",
        "assembly_method": h1_info["assembly_method"], "sequencing_technology": h1_info["sequencing_tech"],
        "annotation_accession_version": annotation["name"], "annotation_provider": annotation["provider"],
        "annotation_release_date": annotation["release_date"], "annotation_software": f"{annotation['pipeline']} {annotation['software_version']}; {annotation['method']}",
        "annotation_status": "native_exact_H1_full_annotation", "annotation_genetic_code": 1,
        "annotation_gff_path": str(gff.resolve()), "annotation_gff_sha256": sha256(gff), "annotation_gff_size_bytes": gff.stat().st_size,
        "annotation_compressed_path": str(gff_gz.resolve()), "annotation_compressed_sha256": sha256(gff_gz), "annotation_provider_md5": provider_md5(candidate, "h1", gff_gz.name),
        "annotation_provenance_path": str(provenance_path.resolve()), "annotation_provenance_sha256": sha256(provenance_path),
        "annotation_contig_map_path": str(mapping.resolve()), "annotation_contig_map_sha256": sha256(mapping), "annotation_contig_map_rows": map_rows,
        "annotation_contig_audit": "passed_exact_sequence_region_name_length_bijection",
        "modality": "base_manifest_pending_SweepGA_IMPG_enrichment", "deposited_variant_status": "unavailable_no_exact_H1_calls_plus_mask_in_authoritative_assembly_release",
        "raw_wfmash_paf_path": qc["artifacts"]["raw_paf"]["path"], "raw_wfmash_paf_sha256": qc["artifacts"]["raw_paf"]["sha256"],
        "callable_bed_path": qc["artifacts"]["callable_bed"]["path"], "callable_bed_sha256": qc["artifacts"]["callable_bed"]["sha256"],
        "callable_bases": qc["callable_bases"], "heterozygous_snvs": qc["heterozygous_snvs"],
        "normalized_bcf_path": str(bcf.resolve()), "normalized_bcf_sha256": sha256(bcf), "normalized_bcf_csi_path": str(csi.resolve()), "normalized_bcf_csi_sha256": sha256(csi),
        "callable_construction": "H2-to-H1 pinned WFMASH; unique query/target projection; MAPQ 1-254; extended =/X; 100bp edges; 10bp indel flanks; both alleles ACGT",
        "callable_exclusions_json": json.dumps(qc["exclusion_counts"], sort_keys=True, separators=(",", ":")),
        "smoke_qc_path": str(smoke_path.resolve()), "smoke_qc_sha256": sha256(smoke_path), "smoke_callable_bases": smoke["callable_bases"], "smoke_parser_status": smoke["status"],
        "guix_manager": "GNU Guix", "guix_channel_commit": channel_commit,
        "guix_channels_path": "analysis/guix/channels.scm", "guix_channels_sha256": sha256(repo_root / "analysis/guix/channels.scm"),
        "guix_manifest_path": "analysis/guix/manifest.scm", "guix_manifest_sha256": sha256(repo_root / "analysis/guix/manifest.scm"),
        "guix_profile_store_path": profile, "wfmash_store_path": wfmash_root,
        "wfmash_commit": "e040aa10e87cab44ed5a4db005e784be62b0bd21", "bcftools_store_path": "/gnu/store/rlb2gljax8lzmhhidbvbzp3al1ad1mww-bcftools-1.14",
        "commands_path": "results/tier3a/acquisition_commands.sh",
        "authoritative_inventory_path": str((data_root / "authoritative_vgp_PRJNA489243_dataset_report.json").resolve()),
        "authoritative_inventory_sha256": sha256(data_root / "authoritative_vgp_PRJNA489243_dataset_report.json"),
    }
    rows.append(row)
if len(rows) < 3:
    raise SystemExit("fewer than three eligible biological tuples")
output_path.parent.mkdir(parents=True, exist_ok=True)
with output_path.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
PY
run_guix python3 "$ROOT/results/tier3a/acquisition_enrich_manifest.py" enrich-manifest \
    --data-root "$DATA_ROOT" --manifest "$ROOT/results/tier3a/acquisition_manifest.tsv" \
    --sweepga-binary "$SWEEPGA_BIN" --impg-binary "$IMPG_BIN"
;;
esac

case "${1:-all}" in
validate|all)
run_guix python3 - "$ROOT/results/tier3a/acquisition_manifest.tsv" "$ROOT/results/tier3a/acquisition_sources.tsv" <<'PY'
import csv
import collections
import hashlib
import json
import sys
from pathlib import Path

import pysam

manifest, sources = map(Path, sys.argv[1:])
def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()

def paf_core(path):
    with Path(path).open(encoding="utf-8") as handle:
        return [tuple(line.rstrip().split("\t")[:12]) for line in handle if line.strip()]

with manifest.open(newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
if len(rows) < 3 or any(row["eligibility_status"] != "eligible_biological" for row in rows):
    raise SystemExit("manifest lacks three eligible biological tuples")
if len({row["dataset_id"] for row in rows}) != len(rows):
    raise SystemExit("duplicate manifest dataset ID")
handoff = json.loads(Path(rows[0]["handoff_validation_path"]).read_text(encoding="utf-8"))
biological = handoff.get("biological", {})
if not handoff.get("pass"):
    raise SystemExit("SweepGA/IMPG handoff proof did not pass")
if handoff.get("short_n_probe", {}).get("accepted") is not False:
    raise SystemExit("pinned SweepGA short -n semantics were not verified")
if [handoff.get("multiplicity_probe", {}).get(str(cap), {}).get("records") for cap in (1, 5, 10)] != [1, 5, 10]:
    raise SystemExit("SweepGA 1/5/10 two-axis multiplicity probe failed")
if biological.get("observed_record_count", 0) <= 0 or biological.get("callable_denominator_bp", 0) <= 0:
    raise SystemExit("biological IMPG proof lacks variants or callable denominator")
if biological.get("native_partition_count", 0) <= 0 or not biological.get("normalized_bcf_indexed") or not biological.get("normalized_vcf_indexed"):
    raise SystemExit("biological IMPG partition/query/index proof failed")
if biological.get("annotation", {}).get("cds_phases_observed") != [0, 1, 2]:
    raise SystemExit("biological IMPG proof lacks all CDS phases")
sequence_qc = biological.get("direct_sequence_validation", {})
if not sequence_qc.get("all_ref_alleles_match_h1") or not sequence_qc.get("h1_with_observed_alts_equals_aligned_h2"):
    raise SystemExit("biological IMPG sequence-level variant validation failed")
path_hash_pairs = [
    ("h1_fasta_path", "h1_fasta_sha256"), ("h1_compressed_path", "h1_compressed_sha256"),
    ("h1_fai_path", "h1_fai_sha256"), ("h2_fasta_path", "h2_fasta_sha256"),
    ("h2_compressed_path", "h2_compressed_sha256"), ("h2_fai_path", "h2_fai_sha256"),
    ("annotation_gff_path", "annotation_gff_sha256"),
    ("annotation_compressed_path", "annotation_compressed_sha256"),
    ("annotation_provenance_path", "annotation_provenance_sha256"),
    ("annotation_contig_map_path", "annotation_contig_map_sha256"),
    ("raw_wfmash_paf_path", "raw_wfmash_paf_sha256"),
    ("sweepga_bounded_paf_path", "sweepga_bounded_paf_sha256"),
    ("sweepga_direct_cap_recheck_path", "sweepga_direct_cap_recheck_sha256"),
    ("sweepga_cap_sensitivity_path", "sweepga_cap_sensitivity_sha256"),
    ("sweepga_binary_path", "sweepga_binary_sha256"),
    ("sweepga_fastga_attempt_report_path", "sweepga_fastga_attempt_report_sha256"),
    ("impg_binary_path", "impg_binary_sha256"),
    ("impg_gfaffix_path", "impg_gfaffix_sha256"),
    ("handoff_validation_path", "handoff_validation_sha256"),
    ("wfmash_sensitivity_bcf_path", "wfmash_sensitivity_bcf_sha256"),
    ("annotation_gene_manifest_path", "annotation_gene_manifest_sha256"),
    ("annotation_query_manifest_path", "annotation_query_manifest_sha256"),
    ("annotation_execution_spans_path", "annotation_execution_spans_sha256"),
    ("annotation_span_feature_map_path", "annotation_span_feature_map_sha256"),
    ("haplotype_contig_map_path", "haplotype_contig_map_sha256"),
    ("annotation_query_qc_path", "annotation_query_qc_sha256"),
    ("guix_toolchain_manifest_path", "guix_toolchain_manifest_sha256"),
    ("guix_tool_versions_path", "guix_tool_versions_sha256"),
    ("staged_object_inventory_path", "staged_object_inventory_sha256"),
    ("callable_bed_path", "callable_bed_sha256"), ("normalized_bcf_path", "normalized_bcf_sha256"),
    ("normalized_bcf_csi_path", "normalized_bcf_csi_sha256"), ("smoke_qc_path", "smoke_qc_sha256"),
    ("authoritative_inventory_path", "authoritative_inventory_sha256"),
]
for row in rows:
    if row["primary_mapping_engine"].split()[0] != "SweepGA":
        raise SystemExit(f"non-SweepGA primary modality for {row['dataset_id']}")
    if row["sweepga_selected_cap_query"] != "1" or row["sweepga_selected_cap_target"] != "1":
        raise SystemExit(f"unexpected promoted cap for {row['dataset_id']}")
    if min(float(row["sweepga_cap1_query_coverage"]), float(row["sweepga_cap1_target_coverage"])) < 0.80:
        raise SystemExit(f"selected cap coverage below floor for {row['dataset_id']}")
    if min(float(row["sweepga_direct_query_coverage"]), float(row["sweepga_direct_target_coverage"])) < 0.80:
        raise SystemExit(f"direct whole-haplotype SweepGA coverage below floor for {row['dataset_id']}")
    if not row["sweepga_direct_cap_fixed_point_status"].startswith("passed_identical_"):
        raise SystemExit(f"direct SweepGA 1:1 fixed-point validation failed for {row['dataset_id']}")
    if row["production_variant_extractor"].split()[0] != "IMPG":
        raise SystemExit(f"non-IMPG production variant extractor for {row['dataset_id']}")
    if min(int(row["targeted_gene_count"]), int(row["queryable_gene_count"]), int(row["queryable_gene_union_bases"]), int(row["queryable_CDS_rows"])) <= 0:
        raise SystemExit(f"empty annotation-derived query inventory for {row['dataset_id']}")
    if int(row["haplotype_contig_map_rows"]) <= 0:
        raise SystemExit(f"empty H1/H2 contig map for {row['dataset_id']}")
    phases = __import__("json").loads(row["queryable_CDS_phase_counts_json"])
    if set(phases) != {"0", "1", "2"} or min(map(int, phases.values())) <= 0:
        raise SystemExit(f"incomplete CDS phase verification for {row['dataset_id']}")
    for path_field, hash_field in path_hash_pairs:
        path = Path(row[path_field])
        if not path.is_file() or not path.stat().st_size or sha256(path) != row[hash_field]:
            raise SystemExit(f"missing/hash-mismatched {path_field} for {row['dataset_id']}")
    if paf_core(row["sweepga_bounded_paf_path"]) != paf_core(row["sweepga_direct_cap_recheck_path"]):
        raise SystemExit(f"direct SweepGA cap fixed-point core records disagree for {row['dataset_id']}")
    if int(row["callable_bases"]) <= 0 or int(row["smoke_callable_bases"]) <= 0:
        raise SystemExit(f"zero callable denominator for {row['dataset_id']}")
    reference = pysam.FastaFile(row["h1_fasta_path"])
    lengths = dict(zip(reference.references, reference.lengths))
    with open(row["annotation_query_manifest_path"], newline="") as query_handle:
        query_rows = list(csv.DictReader(query_handle, delimiter="\t"))
    feature_ids = [item["feature_row_id"] for item in query_rows]
    if len(feature_ids) != len(set(feature_ids)) or not feature_ids:
        raise SystemExit(f"nonunique/empty query feature IDs for {row['dataset_id']}")
    queryable_rows = [item for item in query_rows if item["queryable"] == "yes"]
    observed_phases = collections.Counter(item["phase"] for item in queryable_rows)
    if dict(sorted(observed_phases.items())) != {key: int(value) for key, value in sorted(phases.items())}:
        raise SystemExit(f"query manifest/phase-count disagreement for {row['dataset_id']}")
    if len(queryable_rows) != int(row["queryable_CDS_rows"]):
        raise SystemExit(f"query manifest/queryable-CDS disagreement for {row['dataset_id']}")
    for item in queryable_rows:
        start, end = int(item["start_0based"]), int(item["end_0based_exclusive"])
        if not item["gene_ids"] or not item["transcript_ids"] or not item["cds_id"]:
            raise SystemExit(f"queryable CDS identity missing for {row['dataset_id']}")
        if item["contig"] not in lengths or not 0 <= start < end <= lengths[item["contig"]]:
            raise SystemExit(f"queryable CDS coordinate mismatch for {row['dataset_id']}")
    execution_span_count = execution_span_bases = 0
    execution_previous = {}
    with open(row["annotation_execution_spans_path"], encoding="utf-8") as execution_bed:
        for line_number, line in enumerate(execution_bed, 1):
            contig, start, end, _span_id = line.rstrip().split("\t")[:4]
            start, end = int(start), int(end)
            if contig not in lengths or not 0 <= start < end <= lengths[contig]:
                raise SystemExit(f"invalid annotation execution span {line_number} for {row['dataset_id']}")
            if start <= execution_previous.get(contig, -1):
                raise SystemExit(f"overlapping/unsorted annotation execution spans for {row['dataset_id']}")
            execution_previous[contig] = end
            execution_span_count += 1
            execution_span_bases += end - start
    if execution_span_count != int(row["annotation_execution_span_count"]) or execution_span_bases != int(row["annotation_execution_span_union_bases"]):
        raise SystemExit(f"annotation execution span totals disagree for {row['dataset_id']}")
    with open(row["haplotype_contig_map_path"], newline="") as contig_map_handle:
        contig_map_rows = list(csv.DictReader(contig_map_handle, delimiter="\t"))
    if len(contig_map_rows) != int(row["haplotype_contig_map_rows"]):
        raise SystemExit(f"H1/H2 contig-map row total disagrees for {row['dataset_id']}")
    if {item["h1_contig"] for item in contig_map_rows if item["h1_contig"]} != set(lengths):
        raise SystemExit(f"H1/H2 contig map does not account for every H1 contig in {row['dataset_id']}")
    callable_bases = 0
    previous = {}
    with open(row["callable_bed_path"], encoding="utf-8") as bed:
        for line_number, line in enumerate(bed, 1):
            contig, start, end = line.rstrip().split("\t")[:3]
            start, end = int(start), int(end)
            if contig not in lengths or not 0 <= start < end <= lengths[contig]:
                raise SystemExit(f"invalid callable BED line {line_number} for {row['dataset_id']}")
            if start < previous.get(contig, 0):
                raise SystemExit(f"overlapping/unsorted callable BED for {row['dataset_id']}")
            previous[contig] = end
            callable_bases += end - start
    if callable_bases != int(row["callable_bases"]):
        raise SystemExit(f"callable BED/manifest denominator disagreement for {row['dataset_id']}")
    with pysam.VariantFile(row["normalized_bcf_path"]) as variants:
        dictionary = {name: variants.header.contigs[name].length for name in variants.header.contigs}
        if dictionary != lengths:
            raise SystemExit(f"BCF/H1 dictionary disagreement for {row['dataset_id']}")
        records = sum(1 for _ in variants)
    if records != int(row["heterozygous_snvs"]):
        raise SystemExit(f"BCF/manifest SNV disagreement for {row['dataset_id']}")
    if not row["impg_truth_status"].startswith("passed_native_H1_annotation_biological_"):
        raise SystemExit(f"missing biological IMPG handoff validation for {row['dataset_id']}")
inventory_path = Path(rows[0]["staged_object_inventory_path"])
if any(Path(row["staged_object_inventory_path"]) != inventory_path for row in rows):
    raise SystemExit("manifest rows disagree on staged object inventory")
with inventory_path.open(newline="") as handle:
    inventory_rows = list(csv.DictReader(handle, delimiter="\t"))
listed = set()
for item in inventory_rows:
    path = Path(item["absolute_path"])
    if not path.is_file() or path.stat().st_size != int(item["size_bytes"]) or sha256(path) != item["sha256"]:
        raise SystemExit(f"staged object inventory mismatch: {path}")
    listed.add(path.resolve())
actual = {
    path.resolve()
    for path in inventory_path.parent.rglob("*")
    if path.is_file() and path != inventory_path
}
if listed != actual:
    raise SystemExit(f"staged object inventory coverage mismatch: listed={len(listed)} actual={len(actual)}")
if not sources.is_file() or sum(1 for _ in sources.open()) <= 1:
    raise SystemExit("source ledger is empty")
print(f"validated {len(rows)} eligible biological tuples")
PY
;;
esac

case "${1:-all}" in
versions|all)
record_versions
;;
esac

echo "$DATA_ROOT"
