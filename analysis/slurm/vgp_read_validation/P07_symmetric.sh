#!/usr/bin/env bash
#SBATCH --job-name=vgp-symread-P07
#SBATCH --partition=workers
#SBATCH --cpus-per-task=32
#SBATCH --mem=220G
#SBATCH --time=48:00:00
#SBATCH --output=/moosefs/erikg/vgp/logs/vgp-symmetric-read-P07-%j.out
#SBATCH --error=/moosefs/erikg/vgp/logs/vgp-symmetric-read-P07-%j.err
set -euo pipefail

# Symmetric read-vs-assembly validation for VGP P07 (Spinachia spinachia).
# Addresses review defect I-1: reads are mapped to BOTH haplotypes; the
# assembly callable denominator and SNP sites are lifted H1<->H2 through the
# corrected 1:1 SweepGA PAF; every statistic is computed per frame so no
# single reference is a truth oracle.  Follows analysis/slurm/vgp_read_validation/
# P07_validate.sh conventions (pinned profile, node-local scratch, durable
# promotion) but is additive: no existing script is modified.

readonly PROFILE=/moosefs/erikg/vgp/derived/read-validation/environment/profile
readonly NODE_LOCAL_BASE=/scratch
readonly REQUIRED_SCRATCH_BYTES=450000000000

readonly READS_ROOT=/moosefs/erikg/vgp/views/accession/P07
readonly HIFI_VIEW=$READS_ROOT/SRR25606782/SRR25606782_subreads.fastq.gz
readonly R1_VIEW=$READS_ROOT/SRR30200290/SRR30200290_1.fastq.gz
readonly R2_VIEW=$READS_ROOT/SRR30200290/SRR30200290_2.fastq.gz

readonly ACQ=/moosefs/erikg/tier3data/tier3a-acquisition-20260716/spinachia_spinachia_SK-2024b
readonly H1_FASTA=$ACQ/h1.fna
readonly H2_FASTA=$ACQ/h2.fna
readonly H1_EXPECTED_SHA256=438faaebe34180e563b700d911eb80973589f8ad4d5a70861747067621aaf6ba
readonly H2_EXPECTED_SHA256=0bbd50ea5954e53e47e4bd80f6e01a22e0b80bfef50dd4781b4acaf5f5fb9418

readonly BP=/moosefs/erikg/vgp/pilot/three-pair/vgp-three-pair-20260722-v1/P07/bounded-production
readonly ASSEMBLY_BCF=$BP/variants/normalized.bcf
readonly ASSEMBLY_CALLABLE=$BP/consensus/masks/callable.bed
readonly ASSEMBLY_JOIN_QC=$BP/consensus/join_qc.json
readonly ASSEMBLY_CALLABLE_BP=270531638

readonly PAF=/moosefs/erikg/tier3data/tier3a-origin-remap-20260716/spinachia_spinachia_SK-2024b/mapping/production.1to1.paf

# Pinned input digests (review should-fix #3): every assembly-side input is
# verified after staging; read files are verified against the acquisition
# manifest (vgp_validation_reads_manifest_v1.json local_sha256/observed_bytes).
readonly PAF_EXPECTED_SHA256=bfd2b8c3d23db5078d181b1a49a816fe57b8309d48839c1a1e54265389beb962
readonly BCF_EXPECTED_SHA256=8e96cec90df57e4190d928798b86f8b9623f943f1ad2c597ee2c06c4f5047c74
readonly BCF_CSI_EXPECTED_SHA256=20d2e26c3cc5e240cc07a11c0259c6854a3fa0590b7a8f53f5aef195d9969f58
readonly CALLABLE_EXPECTED_SHA256=5933fd23bcfcf1d9e0ce1b5ef06c5983af0dcbc720a8c48cb27bcfb5367f7384
readonly JOINQC_EXPECTED_SHA256=57f19e7887b4957c3ae286fd813f58d80261caefd96a0c91a7bc3359b6f962aa

# Optional BAM-reuse resume mode: point SYMREAD_REUSE_BAMS at a promoted run
# directory, a bams/ directory, or a bare directory containing
# <platform>.<frame>.bam(+.bai). Each of the four expected BAMs resolves
# independently: a present file must pass samtools quickcheck (hard fail on
# present-but-corrupt); a missing file is mapped fresh (read staging kicks
# back in automatically for mixed partial sets). Reused BAM sha256 digests
# are recorded in input_manifest.json; the reused-vs-mapped decision in
# execution.json. Set SYMREAD_REUSE_DRY_RUN=1 to print the resolution table
# and exit before staging or compute (configuration check; needs no Slurm).
readonly REUSE_BAMS_DIR=${SYMREAD_REUSE_BAMS:-}
readonly DRY_RUN=${SYMREAD_REUSE_DRY_RUN:-}
readonly READS_MANIFEST_STAGE=analysis/vgp_validation_reads_manifest_v1.json

if [[ -n ${SLURM_JOB_ID:-} ]]; then
    readonly REPOSITORY_ROOT=${SLURM_SUBMIT_DIR:?submit from the repository root}
else
    REPOSITORY_ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/../../.." && pwd)
    readonly REPOSITORY_ROOT
fi
cd "$REPOSITORY_ROOT"
readonly READS_MANIFEST="$REPOSITORY_ROOT/$READS_MANIFEST_STAGE"

[[ -x $PROFILE/bin/python3 ]] || { echo "missing pinned profile: $PROFILE" >&2; exit 2; }
export GUIX_PROFILE="$PROFILE"
export GUIX_PYTHONPATH=${GUIX_PYTHONPATH:-}
# shellcheck disable=SC1091
source "$GUIX_PROFILE/etc/profile"
export PYTHONPATH="$REPOSITORY_ROOT"
export LC_ALL=C
export LANG=C

for tool in python3 minimap2 samtools bcftools bedtools sha256sum gzip time; do
    command -v "$tool" >/dev/null || { echo "missing pinned executable: $tool" >&2; exit 2; }
done

# --------------------- BAM reuse resolution (before staging/scratch checks)
declare -A REUSE_SRC=() REUSE_SHA=()
resolve_reuse_bam() {  # echo first existing candidate path for a BAM base name
    local base=$1 candidate
    for candidate in "$REUSE_BAMS_DIR/output/bams/$base" \
                     "$REUSE_BAMS_DIR/bams/$base" \
                     "$REUSE_BAMS_DIR/$base"; do
        [[ -f $candidate ]] && { printf '%s\n' "$candidate"; return 0; }
    done
    return 1
}
if [[ -n $REUSE_BAMS_DIR ]]; then
    [[ -d $REUSE_BAMS_DIR ]] || {
        echo "SYMREAD_REUSE_BAMS is not a directory: $REUSE_BAMS_DIR" >&2; exit 2; }
    for ref in h1 h2; do
        for platform in illumina hifi; do
            base="${platform}.${ref}.bam"
            if src=$(resolve_reuse_bam "$base"); then
                [[ -f ${src}.bai ]] || {
                    echo "reuse BAM present but .bai missing: $src" >&2; exit 2; }
                samtools quickcheck -v "$src" || {
                    echo "reuse BAM failed quickcheck (present-but-corrupt): $src" >&2; exit 2; }
                REUSE_SRC[$base]=$src
                REUSE_SHA[$base]=$(sha256sum "$src" | awk '{print $1}')
            fi
        done
    done
    (( ${#REUSE_SRC[@]} > 0 )) || {
        echo "SYMREAD_REUSE_BAMS resolved none of the 4 expected BAMs under: $REUSE_BAMS_DIR" >&2; exit 2; }
    echo "BAM reuse resolution: ${#REUSE_SRC[@]}/4 reusable from $REUSE_BAMS_DIR"
    for ref in h1 h2; do
        for platform in illumina hifi; do
            base="${platform}.${ref}.bam"
            if [[ -n ${REUSE_SRC[$base]:-} ]]; then
                echo "  REUSE  $base  ${REUSE_SHA[$base]}  ${REUSE_SRC[$base]}"
            else
                echo "  MAP    $base  (fresh mapping; reads will be staged)"
            fi
        done
    done
    if [[ -n $DRY_RUN ]]; then
        echo "SYMREAD_REUSE_DRY_RUN: resolution table above; exiting before staging/compute"
        exit 0
    fi
fi
[[ -n ${SLURM_JOB_ID:-} && -d $NODE_LOCAL_BASE && -w $NODE_LOCAL_BASE ]] || {
    echo "this computation requires Slurm-managed node-local scratch" >&2
    exit 2
}
case $(stat -f -c %T -- "$NODE_LOCAL_BASE") in
    nfs|nfs4|fuse*|lustre|gpfs|ceph) echo "scratch is not node-local" >&2; exit 2 ;;
esac
available_scratch=$(df -PB1 -- "$NODE_LOCAL_BASE" | awk 'NR==2 {print $4}')
(( available_scratch >= REQUIRED_SCRATCH_BYTES )) || {
    echo "insufficient node-local scratch: $available_scratch < $REQUIRED_SCRATCH_BYTES" >&2
    exit 2
}

readonly STARTED_UTC=$(date -u +%Y-%m-%dT%H:%M:%SZ)
readonly WORK=$(mktemp -d "$NODE_LOCAL_BASE/vgp-symread-P07-${SLURM_JOB_ID}-XXXXXX")
readonly OUT="$WORK/output"
mkdir -p "$OUT/transfer" "$OUT/bams" "$OUT/h1" "$OUT/h2" "$OUT/telemetry"
cleanup() {
    status=$?
    if (( status == 0 )); then
        rm -rf -- "$WORK"
    else
        echo "FAILED_WORKDIR=$WORK" >&2
        echo "retaining node-local failure evidence for diagnosis/resume" >&2
    fi
}
trap cleanup EXIT

exec > >(tee "$OUT/telemetry/worker.stdout") 2> >(tee "$OUT/telemetry/worker.stderr" >&2)
cp -- "$0" "$OUT/executed_worker.sh"
cp -- analysis/vgp_symmetric_read_test.py "$OUT/executed_symmetric_module.py"
cp -- analysis/vgp_read_validation.py "$OUT/executed_validation_module.py"

# Module isolation (postmortem for job 2825969): every python -m analysis.*
# invocation must resolve against byte-frozen staged copies, never live
# repo files, so mid-run repository edits cannot alter a running job.
mkdir -p "$WORK/modules/analysis"
cp -- "$REPOSITORY_ROOT/analysis/vgp_symmetric_read_test.py" "$WORK/modules/analysis/"
cp -- "$REPOSITORY_ROOT/analysis/vgp_read_validation.py" "$WORK/modules/analysis/"
cd "$WORK/modules"
export PYTHONPATH=""

run_timed() {
    local label=$1
    shift
    echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] START $label: $*"
    command time -v -o "$OUT/telemetry/${label}.time.txt" -- "$@"
    echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] DONE $label"
}

for input in "$H1_FASTA" "$H2_FASTA" "$ASSEMBLY_BCF" "$ASSEMBLY_BCF.csi" \
    "$ASSEMBLY_CALLABLE" "$ASSEMBLY_JOIN_QC" "$PAF"; do
    [[ -r $input ]] || { echo "missing input: $input" >&2; exit 2; }
done
needs_reads=1
if [[ -n $REUSE_BAMS_DIR ]] && (( ${#REUSE_SRC[@]} == 4 )); then needs_reads=0; fi
if (( needs_reads )); then
    for input in "$HIFI_VIEW" "$R1_VIEW" "$R2_VIEW"; do
        [[ -r $input ]] || { echo "missing input: $input" >&2; exit 2; }
    done
fi

# ---------------------------------------------------------------- stage inputs
if (( needs_reads )); then
    run_timed stage_hifi cp --reflink=auto -- "$HIFI_VIEW" "$WORK/hifi.fastq.gz"
    run_timed stage_r1   cp --reflink=auto -- "$R1_VIEW" "$WORK/r1.fastq.gz"
    run_timed stage_r2   cp --reflink=auto -- "$R2_VIEW" "$WORK/r2.fastq.gz"
fi
run_timed stage_h1   cp -- "$H1_FASTA" "$WORK/h1.fa"
run_timed stage_h2   cp -- "$H2_FASTA" "$WORK/h2.fa"
samtools faidx "$WORK/h1.fa"
samtools faidx "$WORK/h2.fa"
cp -- "$ASSEMBLY_BCF" "$WORK/assembly.bcf"
cp -- "$ASSEMBLY_BCF.csi" "$WORK/assembly.bcf.csi"
cp -- "$ASSEMBLY_CALLABLE" "$WORK/assembly.callable.bed"
cp -- "$ASSEMBLY_JOIN_QC" "$WORK/assembly.join_qc.json"
cp -- "$PAF" "$WORK/production.1to1.paf"

for pair in "h1.fa:$H1_EXPECTED_SHA256" "h2.fa:$H2_EXPECTED_SHA256" \
    "assembly.bcf:$BCF_EXPECTED_SHA256" "assembly.bcf.csi:$BCF_CSI_EXPECTED_SHA256" \
    "assembly.callable.bed:$CALLABLE_EXPECTED_SHA256" \
    "assembly.join_qc.json:$JOINQC_EXPECTED_SHA256" \
    "production.1to1.paf:$PAF_EXPECTED_SHA256"; do
    file=${pair%%:*}
    expected=${pair##*:}
    observed=$(sha256sum "$WORK/$file" | awk '{print $1}')
    [[ $observed == "$expected" ]] || {
        echo "pinned digest mismatch for $file: $observed != $expected" >&2
        exit 2
    }
done
echo "verified staged reference, assembly-callset, and PAF digests"

if [[ -n $REUSE_BAMS_DIR ]]; then
    : > "$WORK/reuse_bams.tsv"
    for base in "${!REUSE_SRC[@]}"; do
        printf '%s\t%s\t%s\n' "$base" "${REUSE_SHA[$base]}" "${REUSE_SRC[$base]}" >> "$WORK/reuse_bams.tsv"
    done
    sort -o "$WORK/reuse_bams.tsv" "$WORK/reuse_bams.tsv"
fi

# ------------------------------- input manifest (every consumed input digest)
export SYMREAD_WORK="$WORK" H1_FASTA H2_FASTA ASSEMBLY_BCF ASSEMBLY_CALLABLE \
    ASSEMBLY_JOIN_QC PAF HIFI_VIEW R1_VIEW R2_VIEW SYMREAD_REUSE_BAMS=${REUSE_BAMS_DIR} \
    REPOSITORY_ROOT SYMREAD_STAGED_MODULES="$WORK/modules/analysis" \
    SYMREAD_REUSE_TSV="$WORK/reuse_bams.tsv"
python3 - "$READS_MANIFEST" "$OUT/input_manifest.json" <<'PY'
import hashlib, json, os, sys
from pathlib import Path

manifest_path, out_path = sys.argv[1], sys.argv[2]
work = Path(os.environ["SYMREAD_WORK"])
reads_manifest = json.loads(Path(manifest_path).read_text())
by_view = {o["accession_view_path"]: o for o in reads_manifest["objects"]}

record = {"schema_version": "vgp-symmetric-read-input-manifest-v1"}

def sha256_of(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()

records = []
for staged, source in (("h1.fa", os.environ["H1_FASTA"]),
                       ("h2.fa", os.environ["H2_FASTA"]),
                       ("assembly.bcf", os.environ["ASSEMBLY_BCF"]),
                       ("assembly.bcf.csi", os.environ["ASSEMBLY_BCF"] + ".csi"),
                       ("assembly.callable.bed", os.environ["ASSEMBLY_CALLABLE"]),
                       ("assembly.join_qc.json", os.environ["ASSEMBLY_JOIN_QC"]),
                       ("production.1to1.paf", os.environ["PAF"])):
    digest = sha256_of(work / staged)
    records.append({"staged": staged, "source_path": source, "sha256": digest,
                    "size_bytes": (work / staged).stat().st_size})

reuse_dir = os.environ.get("SYMREAD_REUSE_BAMS", "")
reuse_tsv = os.environ.get("SYMREAD_REUSE_TSV", "")
bam_records = []
if reuse_tsv and os.path.exists(reuse_tsv):
    for line in Path(reuse_tsv).read_text().splitlines():
        base, digest, source = line.split("\t")
        bam_records.append({"bam": base, "sha256": digest, "source_path": source,
                            "verification": "samtools-quickcheck+sha256-recorded"})

modules_dir = Path(os.environ["SYMREAD_STAGED_MODULES"])
module_records = []
for name in ("vgp_symmetric_read_test.py", "vgp_read_validation.py"):
    staged_mod = modules_dir / name
    module_records.append({"staged": f"modules/analysis/{name}",
                           "source_path": str(Path(os.environ["REPOSITORY_ROOT"]) / "analysis" / name),
                           "sha256": sha256_of(staged_mod),
                           "size_bytes": staged_mod.stat().st_size})

read_records = []
if len(bam_records) < 4:
    for staged, view in (("hifi.fastq.gz", os.environ["HIFI_VIEW"]),
                         ("r1.fastq.gz", os.environ["R1_VIEW"]),
                         ("r2.fastq.gz", os.environ["R2_VIEW"])):
        entry = by_view.get(view)
        if entry is None:
            raise SystemExit(f"reads manifest has no object for view path {view}")
        path = work / staged
        digest = sha256_of(path)
        size = path.stat().st_size
        if digest != entry["local_sha256"]:
            raise SystemExit(f"staged {staged} sha256 {digest} != manifest {entry['local_sha256']}")
        if size != entry["observed_bytes"]:
            raise SystemExit(f"staged {staged} size {size} != manifest {entry['observed_bytes']}")
        read_records.append({"staged": staged, "source_path": view,
                             "sha256": digest, "size_bytes": size,
                             "verified_against": "vgp_validation_reads_manifest_v1.json"})
    record["reads"] = read_records
if reuse_dir:
    record["bam_reuse"] = {"source_dir": reuse_dir, "bams": bam_records}
record["python_modules"] = module_records
record["inputs"] = records
Path(out_path).write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
print(f"wrote {out_path}: {len(records)} pinned inputs + {len(module_records)} staged modules"
      + (f" + {len(bam_records)} reused BAMs" if bam_records else "")
      + (f" + {len(read_records)} verified read files" if read_records else ""))
PY

# ------------------------------------------------- reconstruct H1 denominator
cut -f1,2 "$WORK/h1.fa.fai" > "$WORK/h1.genome"
bcftools query -i 'TYPE!="snp"' -f '%CHROM\t%POS0\t%END\n' "$WORK/assembly.bcf" | \
    bedtools intersect -a stdin -b "$WORK/assembly.callable.bed" -f 1.0 -u | \
    bedtools slop -b 10 -g "$WORK/h1.genome" | \
    bedtools sort -g "$WORK/h1.genome" | bedtools merge > "$WORK/assembly.indel_flanks.bed"
bedtools subtract -a "$WORK/assembly.callable.bed" -b "$WORK/assembly.indel_flanks.bed" \
    > "$WORK/assembly.pi.callable.bed"
python3 - "$WORK/assembly.pi.callable.bed" "$WORK/assembly.join_qc.json" "$ASSEMBLY_CALLABLE_BP" <<'PY'
import json,sys
observed=sum(int(row.split()[2])-int(row.split()[1]) for row in open(sys.argv[1]) if row.strip())
expected=json.load(open(sys.argv[2]))["consensus"]["consensus_callable_bp"]
assert expected == int(sys.argv[3]), "join_qc drifted from frozen expectation"
if observed != expected:
    raise SystemExit(f"final assembly callable denominator mismatch: {observed} != {expected}")
print(f"reproduced inherited final callable denominator: {observed}")
PY

printf 'chrom\tposition_1based\tref\talt\n' > "$WORK/assembly.snps.tsv"
bcftools query -i 'TYPE="snp" && N_ALT=1' -f '%CHROM\t%POS\t%REF\t%ALT\n' \
    "$WORK/assembly.bcf" >> "$WORK/assembly.snps.tsv"

# ------------------------------------------------------ symmetric PAF transfer
run_timed transfer python3 -m analysis.vgp_symmetric_read_test transfer \
    --paf "$WORK/production.1to1.paf" \
    --bed "$WORK/assembly.pi.callable.bed" \
    --sites "$WORK/assembly.snps.tsv" \
    --source-label h1 --target-label h2 \
    --output-dir "$OUT/transfer"
cp -- "$OUT/transfer/transfer_stats.json" "$OUT/transfer_stats.json"

# ------------------------------------------------------------------ mapping
map_one() {  # ref label preset reads...
    local ref=$1 label=$2 preset=$3
    shift 3
    OUT_BAM="$WORK/${label}.${ref}.bam" REF_FA="$WORK/${ref}.fa" RG_ID="$label" PRESET="$preset" \
    run_timed "map_${label}_${ref}" bash -c \
        'set -o pipefail; minimap2 -t 24 -I 8G -ax "$PRESET" --secondary=no -R "@RG\\tID:${RG_ID}\\tSM:P07" "$REF_FA" "$@" | samtools sort -@ 8 -m 4G -o "$OUT_BAM" -' \
        bash "$@"
    samtools index -@ 16 "$WORK/${label}.${ref}.bam"
    samtools quickcheck -v "$WORK/${label}.${ref}.bam"
    samtools flagstat -@ 16 "$WORK/${label}.${ref}.bam" > "$OUT/${ref}/${label}.flagstat.txt"
    samtools stats -@ 16 "$WORK/${label}.${ref}.bam" > "$OUT/${ref}/${label}.stats.txt"
    samtools coverage "$WORK/${label}.${ref}.bam" > "$OUT/${ref}/${label}.coverage.tsv"
    cp -- "$WORK/${label}.${ref}.bam" "$WORK/${label}.${ref}.bam.bai" "$OUT/bams/"
}

reuse_bam() {  # base platform ref — copy a quickcheck+digest-recorded BAM
    local base=$1 platform=$2 ref=$3
    local src=${REUSE_SRC[$base]}
    cp -- "$src" "$src.bai" "$WORK/"
    samtools quickcheck -v "$WORK/$base"
    samtools flagstat -@ 16 "$WORK/$base" > "$OUT/${ref}/${platform}.flagstat.txt"
    samtools stats -@ 16 "$WORK/$base" > "$OUT/${ref}/${platform}.stats.txt"
    samtools coverage "$WORK/$base" > "$OUT/${ref}/${platform}.coverage.tsv"
    cp -- "$WORK/$base" "$WORK/$base.bai" "$OUT/bams/"
}

for ref in h1 h2; do
    ref_fa="$WORK/${ref}.fa"
    case $ref in
        h1) frame_bed="$OUT/transfer/h1.symmetric.pi.callable.bed"
            frame_sites="$OUT/transfer/h1.symmetric.assembly.snps.tsv" ;;
        h2) frame_bed="$OUT/transfer/h2.pi.callable.bed"
            frame_sites="$OUT/transfer/h2.assembly.snps.tsv" ;;
    esac
    cut -f1,2 "$WORK/${ref}.fa.fai" > "$WORK/${ref}.genome"

    for platform in illumina hifi; do
        base="${platform}.${ref}.bam"
        if [[ -n ${REUSE_SRC[$base]:-} ]]; then
            reuse_bam "$base" "$platform" "$ref"
        elif [[ $platform == illumina ]]; then
            map_one "$ref" illumina sr "$WORK/r1.fastq.gz" "$WORK/r2.fastq.gz"
        else
            map_one "$ref" hifi map-hifi "$WORK/hifi.fastq.gz"
        fi
    done

    # ------------------------------------------------ depth masks (dp suite)
    run_timed "depth_masks_${ref}" bash -c \
        'set -o pipefail; samtools depth -aa -q 20 -Q 20 -b "$1" "$2" | python3 -m analysis.vgp_read_validation depth-masks --output-dir "$3" --mask dp5_100:5:100 --mask dp10_60:10:60 --mask dp10_80:10:80 --mask dp15_80:15:80 --mask dp20_80:20:80 --mask dp10_100:10:100' \
        bash "$frame_bed" "$WORK/illumina.${ref}.bam" "$OUT/${ref}/masks"

    # -------------------------------------- calls on primary mask (dp10_80)
    # Masks are derived from ILLUMINA depth per frame (inherited convention);
    # both platforms' calls are restricted to the same per-frame dp10_80 mask
    # so π denominators are identical across platforms within a frame.
    for platform in illumina hifi; do
        run_timed "call_${platform}_${ref}" bash -c \
            'set -o pipefail; bcftools mpileup --threads 8 -Ou -f "$1" -R "$2" -q 20 -Q 20 -a FORMAT/DP,FORMAT/AD "$3" | bcftools call --threads 8 -m -Ob -o "$4"' \
            bash "$ref_fa" "$OUT/${ref}/masks/dp10_80.bed" "$WORK/${platform}.${ref}.bam" "$WORK/read.${platform}.${ref}.raw.bcf"
        bcftools norm --threads 16 -f "$ref_fa" -m -any -Ob -o "$WORK/read.${platform}.${ref}.norm.bcf" "$WORK/read.${platform}.${ref}.raw.bcf"
        bcftools index --threads 16 "$WORK/read.${platform}.${ref}.norm.bcf"
        printf 'chrom\tposition_1based\tref\talt\tquality\tgenotype\tdepth\tallelic_depths\n' \
            > "$OUT/${ref}/${platform}.read.snps.tsv"
        bcftools view --regions-overlap 0 -R "$OUT/${ref}/masks/dp10_80.bed" \
            -i 'QUAL>=30 && TYPE="snp" && N_ALT=1' -Oz -o "$WORK/read.${platform}.${ref}.mask.vcf.gz" "$WORK/read.${platform}.${ref}.norm.bcf"
        bcftools index --threads 8 "$WORK/read.${platform}.${ref}.mask.vcf.gz"
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL[\t%GT\t%DP\t%AD]\n' \
            "$WORK/read.${platform}.${ref}.mask.vcf.gz" >> "$OUT/${ref}/${platform}.read.snps.tsv"
    done

    # ---------------------------- pileup evidence at assembly SNVs, per frame
    awk 'BEGIN{OFS="\t"} NR>1 {print $1,$2-1,$2}' "$frame_sites" > "$WORK/${ref}.assembly.snps.bed"
    for platform in illumina hifi; do
        run_timed "pileup_${platform}_${ref}" samtools mpileup -aa -q 20 -Q 20 \
            -l "$WORK/${ref}.assembly.snps.bed" -f "$ref_fa" \
            -o "$WORK/${platform}.${ref}.assembly-sites.pileup" "$WORK/${platform}.${ref}.bam"
        max_depth=$([[ $platform == illumina ]] && echo 80 || echo 120)
        python3 -m analysis.vgp_read_validation assembly-evidence \
            --assembly-sites "$frame_sites" \
            --pileup "$WORK/${platform}.${ref}.assembly-sites.pileup" \
            --minimum-depth 10 --maximum-depth "$max_depth" \
            --output "$OUT/${ref}/${platform}.assembly_evidence.tsv" \
            --summary "$OUT/${ref}/${platform}.assembly_evidence.json"
        mask_bp=$(python3 - "$OUT/${ref}/masks/depth_mask_summary.json" <<PY
import json,sys
print(json.load(open(sys.argv[1]))["masks"]["dp10_80"]["callable_bp"])
PY
)
        run_timed "metrics_${platform}_${ref}" python3 -m analysis.vgp_symmetric_read_test metrics \
            --frame "$ref" --platform "$platform" \
            --assembly-sites "$frame_sites" \
            --frame-callable-bed "$frame_bed" \
            --read-variants "$OUT/${ref}/${platform}.read.snps.tsv" \
            --assembly-evidence "$OUT/${ref}/${platform}.assembly_evidence.tsv" \
            --callable-bp "$mask_bp" \
            --output "$OUT/${ref}/${platform}.metrics.json"
    done
done

# ------------------------------------------------------------------- report
python3 -m analysis.vgp_symmetric_read_test report \
    --metrics-json "$OUT/h1/illumina.metrics.json" \
    --metrics-json "$OUT/h1/hifi.metrics.json" \
    --metrics-json "$OUT/h2/illumina.metrics.json" \
    --metrics-json "$OUT/h2/hifi.metrics.json" \
    --output "$OUT/symmetric_report.md"

# -------------------------------------------------------------- execution rec
python3 - "$OUT" "$STARTED_UTC" "${SLURM_JOB_ID:-local}" <<'PY'
import hashlib, json, os, platform, shutil, sys, time
from pathlib import Path
out, started, job = Path(sys.argv[1]), sys.argv[2], sys.argv[3]
tools = {}
for name in ("python3", "minimap2", "samtools", "bcftools", "bedtools"):
    resolved = shutil.which(name)
    if resolved is None:
        raise SystemExit(f"missing executable during capture: {name}")
    path = Path(resolved).resolve()
    tools[name] = {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
result = {
    "schema_version": "vgp-symmetric-read-validation-execution-v1",
    "selection_id": "P07",
    "symmetric_design": "reads mapped to H1 and H2 independently; assembly denominator and SNP sites lifted through 1:1 PAF; per-frame metrics",
    "slurm_job_id": job,
    "started_utc": started,
    "completed_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    "node": platform.node(),
    "profile": "/moosefs/erikg/vgp/derived/read-validation/environment/profile",
    "tools": tools,
    "repository_commit": os.popen(f"git -C {os.environ['REPOSITORY_ROOT']} rev-parse HEAD").read().strip(),
    "symmetric_module_sha256": hashlib.sha256((out / "executed_symmetric_module.py").read_bytes()).hexdigest(),
    "input_manifest_sha256": hashlib.sha256((out / "input_manifest.json").read_bytes()).hexdigest(),
    "bam_reuse_source": os.environ.get("SYMREAD_REUSE_BAMS", ""),
    "assembly_pi": {"callable_bp": 270531638, "heterozygous_snps": 316631,
                    "pi": 0.0011704028495181033},
}
reuse_tsv = os.environ.get("SYMREAD_REUSE_TSV", "")
if reuse_tsv and Path(reuse_tsv).exists():
    rows = [line.split("\t") for line in Path(reuse_tsv).read_text().splitlines()]
    reused_names = {r[0] for r in rows}
    mapped = [f"{platform}.{frame}.bam" for frame in ("h1", "h2")
              for platform in ("illumina", "hifi")
              if f"{platform}.{frame}.bam" not in reused_names]
    result["bam_reuse"] = {
        "source_dir": os.environ.get("SYMREAD_REUSE_BAMS", ""),
        "reused": sorted(({"bam": r[0], "sha256": r[1], "source_path": r[2]} for r in rows),
                          key=lambda r: r["bam"]),
        "mapped": sorted(mapped),
    }
(out / "execution.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
PY

# ------------------------------------------------------------ durable output
find "$OUT" -type f ! -name output_manifest.tsv ! -name worker.stdout ! -name worker.stderr -print0 | sort -z | \
    xargs -0 sha256sum | sed "s#  $OUT/#\t#" | awk 'BEGIN{OFS="\t";print "sha256","relative_path"}{print $1,$2}' \
    > "$OUT/output_manifest.tsv"

readonly CANONICAL_PARTIAL="/moosefs/erikg/vgp/staging/outputs/vgp-symmetric-read-P07-${SLURM_JOB_ID}.partial"
readonly CANONICAL_FINAL="/moosefs/erikg/vgp/derived/read-validation/runs/P07-symmetric/slurm-${SLURM_JOB_ID}"
[[ ! -e $CANONICAL_PARTIAL && ! -e $CANONICAL_FINAL ]] || {
    echo "refusing to overwrite validation promotion target" >&2; exit 2;
}
mkdir -p "${CANONICAL_PARTIAL%/*}" "${CANONICAL_FINAL%/*}"
cp -a -- "$OUT" "$CANONICAL_PARTIAL"
(cd "$CANONICAL_PARTIAL" && tail -n +2 output_manifest.tsv | while IFS=$'\t' read -r digest relative; do
    printf '%s  %s\n' "$digest" "$relative"
done | sha256sum -c -)
mv -- "$CANONICAL_PARTIAL" "$CANONICAL_FINAL"
echo "PROMOTED=$CANONICAL_FINAL"
