#!/usr/bin/env bash
set -euo pipefail

ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
CHANNELS="$ROOT/analysis/guix/channels.scm"
MANIFEST="$ROOT/analysis/guix/sweepga_impg_smoke_manifest.scm"
RUN_SUPPLEMENTAL_CONTROLLED=${RUN_SUPPLEMENTAL_CONTROLLED:-0}

if [[ ${1:-} != --inside-guix ]]; then
    exec guix time-machine -C "$CHANNELS" -- shell -m "$MANIFEST" --pure -- \
        env RUN_SUPPLEMENTAL_CONTROLLED="$RUN_SUPPLEMENTAL_CONTROLLED" \
        bash "$0" --inside-guix "${@:1}"
fi
shift

PUBLISH_RESULT=
if [[ $# -gt 0 ]]; then
    OUT=$1
else
    OUT=$(mktemp -d /tmp/sweepga-impg-handoff.XXXXXX)
    PUBLISH_RESULT="$ROOT/analysis/sweepga_impg_observed.json"
    trap 'rm -rf "$OUT"' EXIT
fi
SWEEPGA_SOURCE=/moosefs/erikg/sweepga
IMPG_SOURCE=/moosefs/erikg/impg
IMPG_SYNG_COMMIT=dd00f52b688c0fb78cb7f25336ef9ac9f6a3e109
IMPG_GFAFFIX_COMMIT=460e0dd798a9da7d12aef4f9181419d71489da95
BUILD_TOOLS=/moosefs/erikg/tier3data/tier3a-acquisition-20260716/tools/sweepga-impg-guix-handoff-018e4ce-101df81
SWEEPGA="$BUILD_TOOLS/sweepga"
IMPG="$BUILD_TOOLS/impg"
IMPG_TOOLS="$BUILD_TOOLS"
BIO=/moosefs/erikg/tier3data/tier3a-acquisition-20260716/spinachia_spinachia_SK-2024b
BIO_GFF="$BIO/h1.native.gff"
BIO_CONTIG_MAP="$BIO/h1.annotation_contig_map.tsv"
BIO_PROVENANCE="$BIO/h1.annotation_provenance.json"
EXPECTED="$ROOT/analysis/tests/fixtures/sweepga_impg_expected.tsv"

SWEEPGA_COMMIT=018e4ce49d2c125820e0ac50dc5feaa02d423683
IMPG_COMMIT=101df81eb28a809c8fac97d297acd9fcfbbfa048
SWEEPGA_SHA=1a5440529f5eff91cb7d82876a83a5282df66fb5e2c4b1a9c6caa0bdb83de7b1
IMPG_SHA=c587dc2326cd24f887b1fcb3938404229ad0f0a27ef0773e90c287b1ade160d4
GFAFFIX_SHA=4bc1c5e236a8fe6aa1dbcff6e6cf515e8a70c808549990a515c3c5212776a627

rm -rf "$OUT"
mkdir -p "$OUT"/{controlled,biological,multiplicity,help}
export PATH="$IMPG_TOOLS:$PATH"
export SWEEPGA IMPG
export TMPDIR="$OUT/tmp"
mkdir -p "$TMPDIR"

fail() {
    echo "ERROR: $*" >&2
    exit 1
}

sha_of() {
    sha256sum "$1" | awk '{print $1}'
}

[[ $(git -C "$SWEEPGA_SOURCE" rev-parse HEAD) == "$SWEEPGA_COMMIT" ]] || fail "unexpected SweepGA commit"
[[ $(git -C "$IMPG_SOURCE" rev-parse HEAD) == "$IMPG_COMMIT" ]] || fail "unexpected IMPG commit"
git -C "$SWEEPGA_SOURCE" diff --quiet HEAD -- Cargo.toml Cargo.lock build.rs src || fail "modified SweepGA build input"
git -C "$IMPG_SOURCE" diff --quiet HEAD -- Cargo.toml Cargo.lock build.rs src vendor/gfaffix || fail "modified IMPG build input"
[[ $(git -C "$IMPG_SOURCE/vendor/syng" rev-parse HEAD) == "$IMPG_SYNG_COMMIT" ]] || fail "unexpected IMPG syng submodule commit"
[[ $(git -C "$IMPG_SOURCE/vendor/gfaffix" rev-parse HEAD) == "$IMPG_GFAFFIX_COMMIT" ]] || fail "unexpected IMPG gfaffix submodule commit"
[[ $(sha_of "$SWEEPGA") == "$SWEEPGA_SHA" ]] || fail "unexpected SweepGA executable hash"
[[ $(sha_of "$IMPG") == "$IMPG_SHA" ]] || fail "unexpected IMPG executable hash"
[[ $(sha_of "$IMPG_TOOLS/gfaffix") == "$GFAFFIX_SHA" ]] || fail "unexpected gfaffix executable hash"
while read -r expected_hash tool_name; do
    [[ $(sha_of "$BUILD_TOOLS/$tool_name") == "$expected_hash" ]] || fail "unexpected $tool_name companion hash"
done <<'EOF'
dc5446ff411723862219f378892c94d57abeea882d66a47e96d62b3b27e80b57 ALNtoPAF
8b602f0aab9871a2b666dc6a0a2540cb81eac4d1b9f7461acc01a5c998ea558a FAtoGDB
88fcc7f0eb076a36727677c863557325c2a8ab41a71419c2e50bdcd9a4d3bb4e FastGA
17bed139553e965aa7a1ccb18ff419e936c098c88257a1eabfd4e648d30cf9c0 GIXmake
410ea649af9d018358dda5ca038f3a2e6cd7ca0530f971376d17ea0d4f693400 GIXrm
77d7dd475859c1fae7de095a0a9e67d819c9e2209dc120dc790d1938967c3278 ONEview
f2d181ef7bbb9acf7b113c174b584761424dc6c227020e525dca045d65c3ca6d PAFtoALN
0d8a3a72cfda75a30c38e81b90320c5d212d24b8c312ad22fe97d67e553fc0f6 wfmash
EOF

for binary in "$SWEEPGA" "$IMPG" "$IMPG_TOOLS/gfaffix"; do
    readelf -l "$binary" | grep -q 'Requesting program interpreter: /gnu/store/' ||
        fail "non-Guix ELF interpreter for $binary"
    ldd "$binary" | awk '
        /=> \// && $1 !~ "^/gnu/store/" && $3 !~ "^/gnu/store/" { exit 1 }
    ' || fail "non-Guix runtime dependency for $binary"
done

"$SWEEPGA" --help > "$OUT/help/sweepga.txt"
"$IMPG" index --help > "$OUT/help/impg_index.txt"
"$IMPG" partition --help > "$OUT/help/impg_partition.txt"
"$IMPG" query --help > "$OUT/help/impg_query.txt"
"$IMPG" lace --help > "$OUT/help/impg_lace.txt"

# The assigned task called this option "-n".  Probe it literally and retain
# the negative result: this SweepGA revision only defines --num-mappings.
set +e
"$SWEEPGA" -n 1:1 --version > "$OUT/help/sweepga_short_n.stdout" 2> "$OUT/help/sweepga_short_n.stderr"
SHORT_N_STATUS=$?
set -e
[[ $SHORT_N_STATUS -eq 2 ]] || fail "SweepGA -n probe unexpectedly succeeded"
grep -q "unexpected argument '-n'" "$OUT/help/sweepga_short_n.stderr" || fail "missing -n rejection"

# A dense, completely overlapping PAF fixture demonstrates that M:N is an
# active-overlap cap on both coordinate axes.  The first two records are an
# exact score/start tie; input index is the final deterministic tie-break.
python3 - "$OUT/multiplicity/raw.paf" <<'PY'
import sys
from pathlib import Path

rows = []
for i in range(12):
    length = 1000 if i < 2 else 1000 - i
    label = "tie_first" if i == 0 else ("tie_second" if i == 1 else f"rank_{i + 1}")
    rows.append(
        f"Q#1#chrR\t2000\t0\t{length}\t+\tT#2#chrR\t2000\t0\t{length}"
        f"\t{length}\t{length}\t60\tcg:Z:{length}=\tid:Z:{label}"
    )
Path(sys.argv[1]).write_text("\n".join(rows) + "\n")
PY

for cap in 1 5 10; do
    "$SWEEPGA" "$OUT/multiplicity/raw.paf" \
        --output-file "$OUT/multiplicity/cap${cap}.paf" \
        --num-mappings "${cap}:${cap}" --scaffold-jump 0 \
        --overlap 0.95 --scoring log-length-ani --threads 2 \
        > "$OUT/multiplicity/cap${cap}.stdout" \
        2> "$OUT/multiplicity/cap${cap}.stderr"
done

python3 - "$OUT" <<'PY'
import json
import sys
from pathlib import Path

root = Path(sys.argv[1]) / "multiplicity"
observed = {}
for cap in (1, 5, 10):
    lines = (root / f"cap{cap}.paf").read_text().splitlines()
    observed[str(cap)] = {
        "accepted_form": f"{cap}:{cap}",
        "records": len(lines),
        "labels": [next(x[5:] for x in line.split("\t") if x.startswith("id:Z:")) for line in lines],
    }
    if len(lines) != cap:
        raise SystemExit(f"cap {cap} retained {len(lines)}, expected {cap}")
if observed["1"]["labels"] != ["tie_first"]:
    raise SystemExit("deterministic input-index tie-break changed")
(root / "observed.json").write_text(json.dumps(observed, indent=2, sort_keys=True) + "\n")
PY

# Optional controlled interface regression: one SNP, a 3-bp insertion, and a
# 2-bp deletion.  This is deliberately not a completion gate; the biological
# native-annotation proof below is the scientific acceptance path.
if [[ $RUN_SUPPLEMENTAL_CONTROLLED == 1 ]]; then
python3 - "$OUT/controlled" <<'PY'
import random
import sys
from pathlib import Path

out = Path(sys.argv[1])
rng = random.Random(20260716)
h1 = "".join(rng.choice("ACGT") for _ in range(20000))
alt = {"A": "C", "C": "G", "G": "T", "T": "A"}[h1[5300]]
assert (h1[5300], alt, h1[5499], h1[5699:5702]) == ("G", "T", "G", "GAT")
h2 = list(h1)
h2[5300] = alt
h2 = "".join(h2)
h2 = h2[:5500] + "TTA" + h2[5500:]
h2 = h2[:5703] + h2[5705:]
for name, seq, path in (
    ("H1#1#chrS", h1, out / "h1.fa"),
    ("H2#2#chrS", h2, out / "h2.fa"),
):
    path.write_text(f">{name}\n" + "\n".join(seq[i:i + 80] for i in range(0, len(seq), 80)) + "\n")
PY

"$SWEEPGA" "$OUT/controlled/h1.fa" "$OUT/controlled/h2.fa" \
    --output-file "$OUT/controlled/h1_h2.1to1.paf" \
    --num-mappings 1:1 --scaffold-jump 0 --overlap 0.95 \
    --scoring log-length-ani --threads 2 \
    > "$OUT/controlled/sweepga.stdout" 2> "$OUT/controlled/sweepga.stderr"

"$IMPG" index -a "$OUT/controlled/h1_h2.1to1.paf" \
    -i "$OUT/controlled/h1_h2.impg" -t 2 \
    > "$OUT/controlled/index.stdout" 2> "$OUT/controlled/index.stderr"

mkdir -p "$OUT/controlled/partitions" "$OUT/controlled/calls"
"$IMPG" partition -a "$OUT/controlled/h1_h2.1to1.paf" \
    -i "$OUT/controlled/h1_h2.impg" -w 1000 -d 0 \
    --min-missing-size 1 --min-boundary-distance 0 \
    -o bed --output-folder "$OUT/controlled/partitions" -t 2 \
    > "$OUT/controlled/partition.stdout" 2> "$OUT/controlled/partition.stderr"

python3 - "$OUT/controlled/partitions/partitions.bed" "$OUT/controlled/focus.bed" <<'PY'
import sys
from pathlib import Path

rows = Path(sys.argv[1]).read_text().splitlines()
candidates = [x for x in rows if x.startswith("H1#1#chrS\t")]
row = next(x for x in candidates if int(x.split("\t")[1]) <= 5300 and int(x.split("\t")[2]) > 5702)
fields = row.split("\t")
Path(sys.argv[2]).write_text("\t".join(fields[:3] + ["controlled_truth"]) + "\n")
PY

"$IMPG" query -a "$OUT/controlled/h1_h2.1to1.paf" \
    -i "$OUT/controlled/h1_h2.impg" -b "$OUT/controlled/focus.bed" \
    -d 0 -o vcf:poa --sequence-files "$OUT/controlled/h1.fa" "$OUT/controlled/h2.fa" \
    -O "$OUT/controlled/calls" -t 2 \
    > "$OUT/controlled/query.stdout" 2> "$OUT/controlled/query.stderr"

"$IMPG" lace -f "$OUT/controlled/calls/controlled_truth.vcf" --format vcf \
    -o "$OUT/controlled/laced.vcf" --reference "$OUT/controlled/h1.fa" \
    --compress none -t 2 \
    > "$OUT/controlled/lace.stdout" 2> "$OUT/controlled/lace.stderr"
bcftools norm -f "$OUT/controlled/h1.fa" -m -any -Ob \
    -o "$OUT/controlled/normalized.bcf" "$OUT/controlled/laced.vcf" \
    > "$OUT/controlled/norm.stdout" 2> "$OUT/controlled/norm.stderr"
bcftools index -f "$OUT/controlled/normalized.bcf"
bcftools norm -f "$OUT/controlled/h1.fa" -m -any -Oz \
    -o "$OUT/controlled/normalized.vcf.gz" "$OUT/controlled/laced.vcf" \
    > "$OUT/controlled/norm_vcf.stdout" 2> "$OUT/controlled/norm_vcf.stderr"
bcftools index -f -t "$OUT/controlled/normalized.vcf.gz"
fi

# Biological proof: only extract a 20-kb representative region.  The staged
# whole biological FASTAs remain read-only and SweepGA, not WFMASH, starts the
# mapping proof on these two excerpts.
samtools faidx "$BIO/h1.fna" 'CM106590.1:50001-70000' > "$OUT/biological/h1.raw.fa"
samtools faidx "$BIO/h2.fna" 'CM106672.1:45001-65000' > "$OUT/biological/h2.raw.fa"
python3 - "$OUT/biological" <<'PY'
import sys
from pathlib import Path

out = Path(sys.argv[1])
for raw, dest, name in (
    (out / "h1.raw.fa", out / "h1.fa", "BIOH1#1#CM106590.1"),
    (out / "h2.raw.fa", out / "h2.fa", "BIOH2#2#CM106672.1"),
):
    lines = raw.read_text().splitlines()
    dest.write_text(f">{name}\n" + "\n".join(lines[1:]) + "\n")
PY

# Select targets from the original H1-native GFF, not from variant density.
# Genes must be protein coding and fully contained in the 20-kb excerpt.  CDS
# rows retain their original phase and feature identity for the provenance
# audit.  The one overlapping but boundary-truncated coding gene is recorded
# as excluded rather than silently relabelled as callable.
python3 - "$BIO_GFF" "$BIO_CONTIG_MAP" "$BIO_PROVENANCE" \
    "$OUT/biological/annotation_selection.json" "$OUT/biological/target_genes.bed" \
    "$OUT/biological/target_cds.tsv" <<'PY'
import csv
import json
import sys
from pathlib import Path

gff, contig_map, provenance, result_path, bed_path, cds_path = map(Path, sys.argv[1:])
contig = "CM106590.1"
excerpt_start, excerpt_end = 50001, 70000  # 1-based, closed source coordinates

def attrs(text):
    return dict(item.split("=", 1) for item in text.split(";") if "=" in item)

genes, cds = [], []
with gff.open() as handle:
    for line in handle:
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) != 9 or fields[0] != contig:
            continue
        start, end = int(fields[3]), int(fields[4])
        a = attrs(fields[8])
        if fields[2] == "gene" and a.get("gene_biotype") == "protein_coding" and start <= excerpt_end and end >= excerpt_start:
            genes.append({
                "id": a["ID"], "locus_tag": a["locus_tag"], "start": start,
                "end": end, "strand": fields[6], "name": a.get("Name", a["locus_tag"]),
            })
        elif fields[2] == "CDS" and start <= excerpt_end and end >= excerpt_start:
            cds.append({
                "id": a["ID"], "parent": a["Parent"], "locus_tag": a["locus_tag"],
                "start": start, "end": end, "strand": fields[6], "phase": fields[7],
            })

with contig_map.open() as handle:
    dictionary = {r["annotation_contig"]: r for r in csv.DictReader(handle, delimiter="\t")}
if contig not in dictionary or dictionary[contig]["fasta_contig"] != contig:
    raise SystemExit("annotation contig does not map exactly to the H1 FASTA")
prov = json.loads(provenance.read_text())
if not (prov.get("status") == "native" and prov.get("contig_dictionary_audit_passed")):
    raise SystemExit("annotation provenance is not native/audited")

targeted = [g for g in genes if g["start"] >= excerpt_start and g["end"] <= excerpt_end]
excluded = [g | {"reason": "partial_excerpt_boundary"} for g in genes if g not in targeted]
target_loci = {g["locus_tag"] for g in targeted}
target_cds = [c for c in cds if c["locus_tag"] in target_loci]
if len(targeted) != 3 or len(excluded) != 2 or not target_cds:
    raise SystemExit("native annotation target fixture changed")

def local_interval(row):
    return row["start"] - excerpt_start, row["end"] - excerpt_start + 1

bed_rows = []
for g in targeted:
    start, end = local_interval(g)
    bed_rows.append((start, end, g["id"]))
Path(bed_path).write_text("".join(f"BIOH1#1#{contig}\t{s}\t{e}\t{i}\n" for s, e, i in bed_rows))
with cds_path.open("w") as handle:
    handle.write("contig\tstart_1based\tend_1based\tstrand\tphase\tfeature_id\tparent\tlocus_tag\n")
    for c in target_cds:
        handle.write(f"{contig}\t{c['start']}\t{c['end']}\t{c['strand']}\t{c['phase']}\t{c['id']}\t{c['parent']}\t{c['locus_tag']}\n")

target_bp = sum(e - s for s, e, _ in bed_rows)
excluded_bp_in_excerpt = sum(min(g["end"], excerpt_end) - max(g["start"], excerpt_start) + 1 for g in excluded)
result = {
    "source_gff": str(gff),
    "source_gff_sha256": prov["gff_sha256"],
    "annotation_release": prov["release"],
    "annotation_status": prov["status"],
    "native_vs_projected": prov["native_vs_projected"],
    "contig_dictionary_audit_passed": True,
    "contig": contig,
    "fasta_contig": dictionary[contig]["fasta_contig"],
    "fasta_contig_length": int(dictionary[contig]["length"]),
    "excerpt_source_coordinates_1based_closed": [excerpt_start, excerpt_end],
    "candidate_protein_coding_genes": len(genes),
    "candidate_gene_bp_in_excerpt": target_bp + excluded_bp_in_excerpt,
    "targeted_genes": targeted,
    "targeted_gene_count": len(targeted),
    "targeted_gene_bp": target_bp,
    "excluded_genes": excluded,
    "excluded_gene_count": len(excluded),
    "excluded_gene_bp_in_excerpt": excluded_bp_in_excerpt,
    "cds_feature_rows": len(target_cds),
    "cds_feature_ids": sorted({c["id"] for c in target_cds}),
    "cds_phases_observed": sorted({int(c["phase"]) for c in target_cds}),
}
Path(result_path).write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
PY

"$SWEEPGA" "$OUT/biological/h1.fa" "$OUT/biological/h2.fa" \
    --output-file "$OUT/biological/h1_h2.1to1.paf" \
    --num-mappings 1:1 --scaffold-jump 0 --overlap 0.95 \
    --scoring log-length-ani --threads 2 \
    > "$OUT/biological/sweepga.stdout" 2> "$OUT/biological/sweepga.stderr"

"$IMPG" index -a "$OUT/biological/h1_h2.1to1.paf" \
    -i "$OUT/biological/h1_h2.impg" -t 2 \
    > "$OUT/biological/index.stdout" 2> "$OUT/biological/index.stderr"
mkdir -p "$OUT/biological/partitions" "$OUT/biological/calls"
"$IMPG" partition -a "$OUT/biological/h1_h2.1to1.paf" \
    -i "$OUT/biological/h1_h2.impg" -w 2000 -d 0 \
    --min-missing-size 1 --min-boundary-distance 0 \
    -o bed --output-folder "$OUT/biological/partitions" -t 2 \
    > "$OUT/biological/partition.stdout" 2> "$OUT/biological/partition.stderr"

python3 - "$OUT/biological/partitions/partitions.bed" \
    "$OUT/biological/target_genes.bed" "$OUT/biological/h1_h2.1to1.paf" \
    "$OUT/biological/focus.bed" "$OUT/biological/annotation_query_qc.json" <<'PY'
import json
import sys
from pathlib import Path

partition_path, target_path, paf_path, focus_path, qc_path = map(Path, sys.argv[1:])
targets = [line.split("\t") for line in target_path.read_text().splitlines() if line]
paf = [line.split("\t") for line in paf_path.read_text().splitlines() if line]

def covered(start, end):
    return any(row[0] == "BIOH1#1#CM106590.1" and int(row[2]) <= start and int(row[3]) >= end for row in paf)

queryable = [row for row in targets if covered(int(row[1]), int(row[2]))]
mapping_excluded = [row for row in targets if row not in queryable]
rows, partition_ids = [], set()
for line in partition_path.read_text().splitlines():
    fields = line.split("\t")
    if fields[0] != "BIOH1#1#CM106590.1":
        continue
    pstart, pend = int(fields[1]), int(fields[2])
    if any(pstart < int(t[2]) and int(t[1]) < pend for t in queryable):
        partition_ids.add(fields[3])
        rows.append("\t".join(fields[:3] + [f"bio_partition_{fields[3]}"]))
if not rows:
    raise SystemExit("no native IMPG partitions intersect annotation targets")
focus_path.write_text("\n".join(rows) + "\n")
qc = {
    "queryable_gene_count": len(queryable),
    "queryable_gene_bp": sum(int(r[2]) - int(r[1]) for r in queryable),
    "one_to_one_excluded_gene_count": len(mapping_excluded),
    "one_to_one_excluded_gene_bp": sum(int(r[2]) - int(r[1]) for r in mapping_excluded),
    "one_to_one_excluded_gene_ids": [r[3] for r in mapping_excluded],
    "native_partitions_selected": len(partition_ids),
    "native_focus_bp": sum(int(r.split("\t")[2]) - int(r.split("\t")[1]) for r in rows),
    "selection_rule": "native IMPG partitions intersecting exact queryable gene intervals",
    "additional_preprocessing_padding_bp": 0,
}
qc_path.write_text(json.dumps(qc, indent=2, sort_keys=True) + "\n")
PY

"$IMPG" query -a "$OUT/biological/h1_h2.1to1.paf" \
    -i "$OUT/biological/h1_h2.impg" -b "$OUT/biological/focus.bed" \
    -d 0 -o vcf:poa --sequence-files "$OUT/biological/h1.fa" "$OUT/biological/h2.fa" \
    -O "$OUT/biological/calls" -t 2 \
    > "$OUT/biological/query.stdout" 2> "$OUT/biological/query.stderr"

find "$OUT/biological/calls" -name '*.vcf' -type f | sort > "$OUT/biological/vcf.list"
[[ -s "$OUT/biological/vcf.list" ]] || fail "biological IMPG query emitted no VCFs"
"$IMPG" lace -l "$OUT/biological/vcf.list" --format vcf \
    -o "$OUT/biological/laced.vcf" --reference "$OUT/biological/h1.fa" \
    --compress none -t 2 \
    > "$OUT/biological/lace.stdout" 2> "$OUT/biological/lace.stderr"
bcftools norm -f "$OUT/biological/h1.fa" -m -any -Ob \
    -o "$OUT/biological/normalized.untrimmed.bcf" "$OUT/biological/laced.vcf" \
    > "$OUT/biological/norm.stdout" 2> "$OUT/biological/norm.stderr"
bcftools index -f "$OUT/biological/normalized.untrimmed.bcf"
bcftools view -R "$OUT/biological/target_genes.bed" -Ou \
    "$OUT/biological/normalized.untrimmed.bcf" | \
    bcftools norm -d exact -Ob -o "$OUT/biological/normalized.bcf"
bcftools index -f "$OUT/biological/normalized.bcf"
bcftools view -Oz -o "$OUT/biological/normalized.vcf.gz" \
    "$OUT/biological/normalized.bcf" \
    > "$OUT/biological/norm_vcf.stdout" 2> "$OUT/biological/norm_vcf.stderr"
bcftools index -f -t "$OUT/biological/normalized.vcf.gz"

python3 - "$OUT" "$EXPECTED" "$SHORT_N_STATUS" "$CHANNELS" "$MANIFEST" <<'PY'
import csv
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

out = Path(sys.argv[1])
expected_path = Path(sys.argv[2])
short_status = int(sys.argv[3])
channels = sys.argv[4]
manifest = sys.argv[5]

def paf_metrics(path):
    rows = [line.split("\t") for line in Path(path).read_text().splitlines() if line]
    def max_depth(axis):
        name_col, start_col, end_col = ((0, 2, 3) if axis == "query" else (5, 7, 8))
        by_name = {}
        for row in rows:
            events = by_name.setdefault(row[name_col], [])
            events.extend(((int(row[start_col]), 1), (int(row[end_col]), -1)))
        best = 0
        for events in by_name.values():
            depth = 0
            for _, delta in sorted(events, key=lambda x: (x[0], x[1])):
                depth += delta
                best = max(best, depth)
        return best
    return {"records": len(rows), "max_query_overlap_depth": max_depth("query"), "max_target_overlap_depth": max_depth("target")}

def query_records(path):
    text = subprocess.check_output(["bcftools", "query", "-f", "%CHROM\\t%POS\\t%REF\\t%ALT\\n", str(path)], text=True)
    return [line.split("\t") for line in text.splitlines() if line]

with expected_path.open() as handle:
    expected = [row for row in csv.DictReader(handle, delimiter="\t") if row["proof"] == "controlled"]
expected_records = [[r["chrom"], r["pos"], r["ref"], r["alt"]] for r in expected]
controlled_executed = (out / "controlled/normalized.bcf").is_file()
controlled_records = query_records(out / "controlled/normalized.bcf") if controlled_executed else []
bio_records = query_records(out / "biological/normalized.bcf")
bio_untrimmed_records = query_records(out / "biological/normalized.untrimmed.bcf")
bio_mapping = paf_metrics(out / "biological/h1_h2.1to1.paf")
focus_rows = [line.split("\t") for line in (out / "biological/focus.bed").read_text().splitlines() if line]
callable_bp = sum(int(row[2]) - int(row[1]) for row in focus_rows)
partitions = {row[3] for row in focus_rows}
multiplicity = json.loads((out / "multiplicity/observed.json").read_text())
annotation = json.loads((out / "biological/annotation_selection.json").read_text())
annotation_qc = json.loads((out / "biological/annotation_query_qc.json").read_text())
target_rows = [line.split("\t") for line in (out / "biological/target_genes.bed").read_text().splitlines() if line]

def fasta_sequence(path):
    return "".join(line for line in Path(path).read_text().splitlines() if not line.startswith(">" )).upper()

bio_paf = next(line.split("\t") for line in (out / "biological/h1_h2.1to1.paf").read_text().splitlines() if line)
qstart, qend, strand, tstart, tend = int(bio_paf[2]), int(bio_paf[3]), bio_paf[4], int(bio_paf[7]), int(bio_paf[8])
if strand != "+":
    raise SystemExit("biological direct sequence verifier currently requires the observed plus-strand mapping")
h1 = fasta_sequence(out / "biological/h1.fa")
h2 = fasta_sequence(out / "biological/h2.fa")
pieces, last, allele_checks = [], qstart, []
for chrom, pos_text, ref, alt in bio_records:
    pos = int(pos_text)
    zero = pos - 1
    ref_matches = h1[zero:zero + len(ref)] == ref.upper()
    gene_ids = [r[3] for r in target_rows if int(r[1]) < zero + len(ref) and zero < int(r[2])]
    allele_checks.append({
        "chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
        "kind": "SNP" if len(ref) == len(alt) == 1 else ("insertion" if len(alt) > len(ref) else "deletion"),
        "h1_ref_matches": ref_matches, "target_gene_ids": gene_ids,
    })
    if not ref_matches or not gene_ids or zero < last or zero + len(ref) > qend:
        raise SystemExit("biological allele failed H1/annotation/aligned-span validation")
    pieces.extend((h1[last:zero], alt.upper()))
    last = zero + len(ref)
pieces.append(h1[last:qend])
reconstructed_alt = "".join(pieces)
alt_haplotype_matches_h2 = reconstructed_alt == h2[tstart:tend]
variant_kinds = {r["kind"] for r in allele_checks}

controlled_vcf_indexed = (out / "controlled/normalized.vcf.gz.tbi").is_file()
controlled_bcf_indexed = (out / "controlled/normalized.bcf.csi").is_file()
bio_vcf_indexed = (out / "biological/normalized.vcf.gz.tbi").is_file()
bio_bcf_indexed = (out / "biological/normalized.bcf.csi").is_file()

result = {
    "schema_version": "sweepga-impg-handoff-v1",
    "pass": (
        bool(bio_records)
        and bio_mapping["records"] > 0
        and bio_mapping["max_query_overlap_depth"] <= 1
        and bio_mapping["max_target_overlap_depth"] <= 1
        and callable_bp > 0
        and bio_vcf_indexed
        and bio_bcf_indexed
        and annotation["annotation_status"] == "native"
        and annotation["cds_phases_observed"] == [0, 1, 2]
        and annotation_qc["queryable_gene_count"] > 0
        and annotation_qc["queryable_gene_bp"] > 0
        and alt_haplotype_matches_h2
        and "SNP" in variant_kinds
        and bool({"insertion", "deletion"} & variant_kinds)
    ),
    "environment": {
        "guix_channel_commit": "44bbfc24e4bcc48d0e3343cd3d83452721af8c36",
        "channels_file": channels,
        "manifest_file": manifest,
        "guix_profile": os.environ["GUIX_ENVIRONMENT"],
        "rust_toolchain_source_sha256": "294b3d81fa72e62581276290c60c81eb8b58498d333d422ca1dfc432877d0c40",
        "patchelf_source_sha256": "1952b2a782ba576279c211ee942e341748fdb44997f704dd53def46cd055470b",
        "sweepga_commit": "018e4ce49d2c125820e0ac50dc5feaa02d423683",
        "sweepga_version": subprocess.check_output([os.environ["SWEEPGA"], "--version"], text=True).strip(),
        "sweepga_sha256": hashlib.sha256(Path(os.environ["SWEEPGA"]).read_bytes()).hexdigest(),
        "impg_commit": "101df81eb28a809c8fac97d297acd9fcfbbfa048",
        "impg_syng_submodule_commit": "dd00f52b688c0fb78cb7f25336ef9ac9f6a3e109",
        "impg_gfaffix_submodule_commit": "460e0dd798a9da7d12aef4f9181419d71489da95",
        "impg_version": subprocess.check_output([os.environ["IMPG"], "--version"], text=True).strip(),
        "impg_sha256": hashlib.sha256(Path(os.environ["IMPG"]).read_bytes()).hexdigest(),
        "gfaffix_sha256": hashlib.sha256((Path(os.environ["IMPG"]).parent / "gfaffix").read_bytes()).hexdigest(),
        "sweepga_companion_sha256": {
            name: hashlib.sha256((Path(os.environ["SWEEPGA"]).parent / name).read_bytes()).hexdigest()
            for name in ("ALNtoPAF", "FAtoGDB", "FastGA", "GIXmake", "GIXrm", "ONEview", "PAFtoALN", "wfmash")
        },
        "cargo_version": "cargo 1.94.1 (29ea6fb6a 2026-03-24)",
        "rustc_version": "rustc 1.94.1 (e408947bf 2026-03-25)",
        "bcftools_version": subprocess.check_output(["bcftools", "--version"], text=True).splitlines()[0],
        "samtools_version": subprocess.check_output(["samtools", "--version"], text=True).splitlines()[0],
    },
    "short_n_probe": {
        "command": "sweepga -n 1:1 --version",
        "exit_status": short_status,
        "accepted": False,
        "supported_option": "--num-mappings 1:1",
    },
    "multiplicity_probe": multiplicity,
    "controlled": {
        "completion_gate": False,
        "executed": controlled_executed,
        "sweepga_mapping": paf_metrics(out / "controlled/h1_h2.1to1.paf") if controlled_executed else None,
        "native_partition_count": len({line.split("\t")[3] for line in (out / "controlled/partitions/partitions.bed").read_text().splitlines()}) if controlled_executed else 0,
        "focused_callable_bp": sum(int(x.split("\t")[2]) - int(x.split("\t")[1]) for x in (out / "controlled/focus.bed").read_text().splitlines()) if controlled_executed else 0,
        "expected_records": expected_records,
        "observed_records": controlled_records,
        "exact_match": controlled_executed and controlled_records == expected_records,
        "normalized_vcf_indexed": controlled_vcf_indexed,
        "normalized_bcf_indexed": controlled_bcf_indexed,
    },
    "biological": {
        "pair": "Spinachia spinachia SK-2024b H1/H2",
        "source_h1": "/moosefs/erikg/tier3data/tier3a-acquisition-20260716/spinachia_spinachia_SK-2024b/h1.fna:CM106590.1:50001-70000",
        "source_h2": "/moosefs/erikg/tier3data/tier3a-acquisition-20260716/spinachia_spinachia_SK-2024b/h2.fna:CM106672.1:45001-65000",
        "sweepga_mapping": bio_mapping,
        "native_partition_count": len(partitions),
        "native_partition_focus_bp": callable_bp,
        "callable_denominator_bp": annotation_qc["queryable_gene_bp"],
        "callable_scope": "queryable H1-native protein-coding gene spans; not genome-wide",
        "observed_record_count": len(bio_records),
        "observed_records": bio_records,
        "annotation": annotation,
        "annotation_query_qc": annotation_qc | {
            "native_partition_context_bp_beyond_targets": callable_bp - annotation_qc["queryable_gene_bp"],
            "untrimmed_record_count": len(bio_untrimmed_records),
            "trimmed_deduplicated_record_count": len(bio_records),
            "records_removed_by_target_trim_or_exact_dedup": len(bio_untrimmed_records) - len(bio_records),
        },
        "direct_sequence_validation": {
            "aligned_h1_local_half_open": [qstart, qend],
            "aligned_h2_local_half_open": [tstart, tend],
            "all_ref_alleles_match_h1": all(r["h1_ref_matches"] for r in allele_checks),
            "h1_with_observed_alts_equals_aligned_h2": alt_haplotype_matches_h2,
            "alleles": allele_checks,
        },
        "nonzero_mapping": paf_metrics(out / "biological/h1_h2.1to1.paf")["records"] > 0,
        "nonzero_callable_denominator": annotation_qc["queryable_gene_bp"] > 0,
        "normalized_vcf_indexed": bio_vcf_indexed,
        "normalized_bcf_indexed": bio_bcf_indexed,
    },
}
(out / "observed_vs_expected.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
if not result["pass"]:
    raise SystemExit("observed-versus-expected assertions failed")
print(json.dumps({"pass": True, "result": str(out / 'observed_vs_expected.json')}))
PY

if [[ -n $PUBLISH_RESULT ]]; then
    cp "$OUT/observed_vs_expected.json" "$PUBLISH_RESULT"
    echo "published machine-readable result: $PUBLISH_RESULT"
fi
