#!/usr/bin/env bash
#SBATCH --job-name=vgp-readval-P07
#SBATCH --partition=workers
#SBATCH --cpus-per-task=32
#SBATCH --mem=220G
#SBATCH --time=2-00:00:00
#SBATCH --output=/moosefs/erikg/vgp/logs/vgp-read-validation-P07-%j.out
#SBATCH --error=/moosefs/erikg/vgp/logs/vgp-read-validation-P07-%j.err
set -euo pipefail

readonly TASK_ID=validate-vgp-pilot-reads
readonly SELECTION_ID=P07
readonly VGP_ROOT_CONFIG=analysis/vgp_data_root_config.json
readonly PROFILE=/moosefs/erikg/vgp/derived/read-validation/environment/profile
readonly NODE_LOCAL_BASE=/scratch
readonly REQUIRED_SCRATCH_BYTES=400000000000

if [[ -n ${SLURM_JOB_ID:-} ]]; then
    readonly REPOSITORY_ROOT=${SLURM_SUBMIT_DIR:?submit from the repository root}
else
    REPOSITORY_ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/../../.." && pwd)
    readonly REPOSITORY_ROOT
fi
cd "$REPOSITORY_ROOT"

readonly VGP_ROOT=$(
    python3 - "$VGP_ROOT_CONFIG" <<'PY'
import json,sys
value=json.load(open(sys.argv[1]))["root"]
if value != "/moosefs/erikg/vgp":
    raise SystemExit(f"refusing noncanonical VGP root: {value}")
print(value)
PY
)
[[ -x $PROFILE/bin/python3 ]] || { echo "missing pinned profile: $PROFILE" >&2; exit 2; }
export GUIX_PROFILE="$PROFILE"
# Guix's generated profile uses `${GUIX_PYTHONPATH:+...}`; initialize it
# before sourcing under `set -u` on clean Slurm nodes.
export GUIX_PYTHONPATH=${GUIX_PYTHONPATH:-}
# shellcheck disable=SC1091
source "$GUIX_PROFILE/etc/profile"
export PYTHONPATH="$REPOSITORY_ROOT"
export LC_ALL=C
export LANG=C

for tool in python3 minimap2 samtools bcftools bedtools jellyfish psmc fq2psmcfa sha256sum gzip time; do
    command -v "$tool" >/dev/null || { echo "missing pinned executable: $tool" >&2; exit 2; }
done
[[ -n ${SLURM_JOB_ID:-} && -d $NODE_LOCAL_BASE && -w $NODE_LOCAL_BASE ]] || {
    echo "this biological computation requires Slurm-managed node-local scratch" >&2
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
readonly WORK=$(mktemp -d "$NODE_LOCAL_BASE/vgp-readval-P07-${SLURM_JOB_ID}-XXXXXX")
readonly OUT="$WORK/output"
readonly TELEMETRY="$OUT/telemetry"
mkdir -p "$OUT" "$TELEMETRY" "$OUT/masks" "$OUT/kmer" "$OUT/sites" "$OUT/psmc"
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

exec > >(tee "$TELEMETRY/worker.stdout") 2> >(tee "$TELEMETRY/worker.stderr" >&2)
cp -- "$0" "$OUT/executed_worker.sh"
cp -- analysis/vgp_read_validation.py "$OUT/executed_validation_module.py"
cp -- "$VGP_ROOT_CONFIG" "$OUT/root_config.json"

run_timed() {
    local label=$1
    shift
    echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] START $label: $*"
    command time -v -o "$TELEMETRY/${label}.time.txt" -- "$@"
    echo "[$(date -u +%Y-%m-%dT%H:%M:%SZ)] DONE $label"
}

readonly READ_MANIFEST="$VGP_ROOT/pilot/manifests/vgp_validation_reads_manifest_v1.json"
readonly HIFI_VIEW="$VGP_ROOT/views/accession/P07/SRR25606782/SRR25606782_subreads.fastq.gz"
readonly R1_VIEW="$VGP_ROOT/views/accession/P07/SRR30200290/SRR30200290_1.fastq.gz"
readonly R2_VIEW="$VGP_ROOT/views/accession/P07/SRR30200290/SRR30200290_2.fastq.gz"
readonly H1_BGZF="$VGP_ROOT/derived/freeze1-bgzf/objects/GCA/048/126/635/GCA_048126635.1/GCA_048126635.1.fa.gz"
readonly H1_PROVENANCE="${H1_BGZF%/*}/provenance.json"
readonly CORE="$VGP_ROOT/pilot/outputs/vgp10-auth-20260718-v2/P07/core"
readonly ASSEMBLY_BCF="$CORE/variants/normalized.bcf"
readonly ASSEMBLY_CALLABLE="$CORE/consensus/masks/callable.bed"
readonly ASSEMBLY_JOIN_QC="$CORE/consensus/join_qc.json"
readonly ASSEMBLY_PSMC="$CORE/psmc/replicate-000/unscaled.psmc"

for input in "$READ_MANIFEST" "$HIFI_VIEW" "$R1_VIEW" "$R2_VIEW" "$H1_BGZF" \
    "$H1_PROVENANCE" "$ASSEMBLY_BCF" "$ASSEMBLY_CALLABLE" "$ASSEMBLY_JOIN_QC" "$ASSEMBLY_PSMC"; do
    [[ -r $input ]] || { echo "missing input: $input" >&2; exit 2; }
done

run_timed stage_hifi cp --reflink=auto -- "$HIFI_VIEW" "$WORK/hifi.fastq.gz"
run_timed stage_r1 cp --reflink=auto -- "$R1_VIEW" "$WORK/r1.fastq.gz"
run_timed stage_r2 cp --reflink=auto -- "$R2_VIEW" "$WORK/r2.fastq.gz"
run_timed stage_h1 bash -c 'gzip -cd -- "$1" > "$2"' bash "$H1_BGZF" "$WORK/h1.fa"
samtools faidx "$WORK/h1.fa"
cp -- "$ASSEMBLY_CALLABLE" "$WORK/assembly.callable.bed"
cp -- "$ASSEMBLY_BCF" "$WORK/assembly.bcf"
cp -- "$ASSEMBLY_BCF.csi" "$WORK/assembly.bcf.csi"

python3 - "$READ_MANIFEST" "$H1_PROVENANCE" "$WORK" "$OUT/input_manifest.json" \
    "$HIFI_VIEW" "$R1_VIEW" "$R2_VIEW" "$ASSEMBLY_BCF" "$ASSEMBLY_CALLABLE" "$ASSEMBLY_JOIN_QC" "$ASSEMBLY_PSMC" <<'PY'
import hashlib,json,os,sys
from pathlib import Path
def digest(path):
    value=hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda:handle.read(8*1024*1024),b""):
            value.update(block)
    return value.hexdigest()
manifest_path,provenance_path,work,out,*other=map(Path,sys.argv[1:])
manifest=json.loads(manifest_path.read_text())
if manifest.get("canonical_root") != "/moosefs/erikg/vgp":
    raise SystemExit("read manifest canonical root drift")
expected={row["filename"]:row for row in manifest["objects"] if row["selection_id"] == "P07"}
staged={"SRR25606782_subreads.fastq.gz":work/"hifi.fastq.gz",
        "SRR30200290_1.fastq.gz":work/"r1.fastq.gz",
        "SRR30200290_2.fastq.gz":work/"r2.fastq.gz"}
records=[]
for name,path in staged.items():
    row=expected[name]
    if row.get("status") not in {"verified", "reused"} or not row.get("local_sha256"):
        raise SystemExit(f"P07 object not verified/reused in canonical manifest: {name}")
    observed_digest=digest(path)
    if observed_digest != row["local_sha256"] or path.stat().st_size != row["expected_bytes"]:
        raise SystemExit(f"staged P07 digest/size mismatch: {name}")
    records.append({"role":name,"source":row["accession_view_path"],"staged_sha256":observed_digest,
                    "bytes":path.stat().st_size,"upstream_md5":row["upstream_checksum"]})
provenance=json.loads(provenance_path.read_text())
h1=work/"h1.fa"
h1_digest=digest(h1)
if h1_digest != provenance["derived_decompressed_sha256"]:
    raise SystemExit("decompressed H1 digest mismatch")
for path in other:
    records.append({"role":path.name,"source":str(path),"bytes":path.stat().st_size,
                    "sha256":digest(path)})
records.append({"role":"H1_reference_decompressed","source":provenance["derived_path"],
                "bytes":h1.stat().st_size,"sha256":h1_digest,
                "accession_version":"GCA_048126635.1"})
result={"schema_version":"vgp-read-validation-input-manifest-v1","task_id":"validate-vgp-pilot-reads",
        "selection_id":"P07","canonical_vgp_root":"/moosefs/erikg/vgp",
        "exact_individual":"fSpiSpi1","biosample":"SAMN36735485","inputs":records}
out.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
PY

# Reconstruct the exact final callable denominator used by the inherited π:
# only non-SNP alleles fully contained in the coarse assembly mask generate
# ±10-bp flanks. This must reproduce consensus_callable_bp from join_qc.json.
cut -f1,2 "$WORK/h1.fa.fai" > "$WORK/h1.genome"
bcftools query -i 'TYPE!="snp"' -f '%CHROM\t%POS0\t%END\n' "$WORK/assembly.bcf" | \
    bedtools intersect -a stdin -b "$WORK/assembly.callable.bed" -f 1.0 -u | \
    bedtools slop -b 10 -g "$WORK/h1.genome" | \
    bedtools sort -g "$WORK/h1.genome" | bedtools merge > "$WORK/assembly.indel_flanks.bed"
bedtools subtract -a "$WORK/assembly.callable.bed" -b "$WORK/assembly.indel_flanks.bed" \
    > "$WORK/assembly.pi.callable.bed"
python3 - "$WORK/assembly.pi.callable.bed" "$ASSEMBLY_JOIN_QC" <<'PY'
import json,sys
observed=sum(int(row.split()[2])-int(row.split()[1]) for row in open(sys.argv[1]) if row.strip())
expected=json.load(open(sys.argv[2]))["consensus"]["consensus_callable_bp"]
if observed != expected:
    raise SystemExit(f"final assembly callable denominator mismatch: {observed} != {expected}")
print(f"reproduced inherited final callable denominator: {observed}")
PY

run_timed map_illumina bash -c \
    'minimap2 -t 20 -I 8G -ax sr --secondary=no -R "@RG\\tID:SRR30200290\\tSM:P07" "$1" "$2" "$3" | samtools sort -@ 12 -m 4G -o "$4" -' \
    bash "$WORK/h1.fa" "$WORK/r1.fastq.gz" "$WORK/r2.fastq.gz" "$WORK/illumina.bam"
samtools index -@ 16 "$WORK/illumina.bam"
samtools quickcheck -v "$WORK/illumina.bam"
samtools flagstat -@ 16 "$WORK/illumina.bam" > "$OUT/illumina.flagstat.txt"
samtools stats -@ 16 "$WORK/illumina.bam" > "$OUT/illumina.stats.txt"
samtools coverage "$WORK/illumina.bam" > "$OUT/illumina.coverage.tsv"

run_timed map_hifi bash -c \
    'minimap2 -t 24 -I 8G -ax map-hifi --secondary=no -R "@RG\\tID:SRR25606782\\tSM:P07" "$1" "$2" | samtools sort -@ 8 -m 4G -o "$3" -' \
    bash "$WORK/h1.fa" "$WORK/hifi.fastq.gz" "$WORK/hifi.bam"
samtools index -@ 16 "$WORK/hifi.bam"
samtools quickcheck -v "$WORK/hifi.bam"
samtools flagstat -@ 16 "$WORK/hifi.bam" > "$OUT/hifi.flagstat.txt"
samtools stats -@ 16 "$WORK/hifi.bam" > "$OUT/hifi.stats.txt"
samtools coverage "$WORK/hifi.bam" > "$OUT/hifi.coverage.tsv"

# Jellyfish 2.3.0 rejects gzip containers. Decompress the two concatenable
# FASTQ streams explicitly so the count covers both mates without truncation.
run_timed kmer_count bash -c \
    'gzip -cd -- "$1" "$2" | jellyfish count -C -m 21 -s 1000000000 -t 32 --quality-start=33 --min-quality=20 -o "$3" /dev/fd/0' \
    bash "$WORK/r1.fastq.gz" "$WORK/r2.fastq.gz" "$OUT/kmer/illumina.k21.jf"
jellyfish histo -h 300 -t 32 "$OUT/kmer/illumina.k21.jf" | \
    awk 'BEGIN{OFS="\t"; print "depth","count"} {print $1,$2}' > "$OUT/kmer/illumina.k21.histo.tsv"
python3 -m analysis.vgp_read_validation kmer-model \
    --histogram "$OUT/kmer/illumina.k21.histo.tsv" --k 21 --minimum-depth 5 \
    --output "$OUT/kmer/heterozygosity.json"
run_timed kmer_qv_query bash -c \
    'jellyfish query "$1" -s "$2" | python3 -m analysis.vgp_read_validation kmer-qv --k 21 --trusted-minimum 5 --output "$3"' \
    bash "$OUT/kmer/illumina.k21.jf" "$WORK/h1.fa" "$OUT/kmer/H1.qv.json"

run_timed depth_masks bash -c \
    'samtools depth -aa -q 20 -Q 20 -b "$1" "$2" | python3 -m analysis.vgp_read_validation depth-masks --output-dir "$3" --mask dp5_100:5:100 --mask dp10_60:10:60 --mask dp10_80:10:80 --mask dp15_80:15:80 --mask dp20_80:20:80 --mask dp10_100:10:100' \
    bash "$WORK/assembly.pi.callable.bed" "$WORK/illumina.bam" "$OUT/masks"

run_timed call_variants bash -c \
    'bcftools mpileup --threads 8 -Ou -f "$1" -R "$2" -q 20 -Q 20 -a FORMAT/DP,FORMAT/AD "$3" | bcftools call --threads 8 -m -Ob -o "$4"' \
    bash "$WORK/h1.fa" "$WORK/assembly.pi.callable.bed" "$WORK/illumina.bam" "$WORK/read.raw.bcf"
bcftools norm --threads 16 -f "$WORK/h1.fa" -m -any -Ob -o "$WORK/read.norm.bcf" "$WORK/read.raw.bcf"
bcftools index --threads 16 "$WORK/read.norm.bcf"
bcftools stats "$WORK/read.norm.bcf" > "$OUT/read.normalized.bcftools.stats.txt"

for mask_id in dp5_100 dp10_60 dp10_80 dp15_80 dp20_80 dp10_100; do
    mask="$OUT/masks/$mask_id.bed"
    mask_dir="$OUT/masks/$mask_id"
    mkdir -p "$mask_dir"
    callable_bp=$(python3 - "$OUT/masks/depth_mask_summary.json" "$mask_id" <<'PY'
import json,sys
print(json.load(open(sys.argv[1]))["masks"][sys.argv[2]]["callable_bp"])
PY
)
    printf 'chrom\tposition_1based\tref\talt\n' > "$mask_dir/assembly.snps.tsv"
    bcftools query --regions-overlap 0 -R "$mask" -i 'TYPE="snp" && N_ALT=1' \
        -f '%CHROM\t%POS\t%REF\t%ALT\n' "$WORK/assembly.bcf" >> "$mask_dir/assembly.snps.tsv"

    # Materialize the mask-specific variant VCF once, then query that compact
    # product.  Re-scanning the 2-GiB all-sites BCF separately for the TSV and
    # consensus would duplicate ~8 minutes of I/O per sensitivity mask.
    bcftools view --regions-overlap 0 -R "$mask" -i 'QUAL>=30 && TYPE="snp" && N_ALT=1' \
        -Oz -o "$WORK/$mask_id.vcf.gz" "$WORK/read.norm.bcf"
    bcftools index --threads 8 "$WORK/$mask_id.vcf.gz"
    printf 'chrom\tposition_1based\tref\talt\tquality\tgenotype\tdepth\tallelic_depths\n' > "$mask_dir/read.snps.tsv"
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL[\t%GT\t%DP\t%AD]\n' \
        "$WORK/$mask_id.vcf.gz" >> "$mask_dir/read.snps.tsv"
    python3 -m analysis.vgp_read_validation mask-report \
        --assembly-variants "$mask_dir/assembly.snps.tsv" \
        --read-variants "$mask_dir/read.snps.tsv" --callable-bp "$callable_bp" \
        --mask-id "$mask_id" --output "$mask_dir/comparison.json"

    bedtools complement -i "$mask" -g "$WORK/h1.genome" > "$WORK/$mask_id.noncallable.bed"
    bcftools consensus -s P07 -H I -m "$WORK/$mask_id.noncallable.bed" \
        -f "$WORK/h1.fa" "$WORK/$mask_id.vcf.gz" > "$WORK/$mask_id.consensus.fa"
    python3 -m analysis.vgp_read_validation fasta-to-fastq \
        < "$WORK/$mask_id.consensus.fa" > "$WORK/$mask_id.consensus.fastq"
    # PSMC 0.6.5's fq2psmcfa requires an explicit FASTQ path; it does not
    # accept the sequence stream on stdin.
    fq2psmcfa -q 20 "$WORK/$mask_id.consensus.fastq" > "$WORK/$mask_id.psmcfa"
    run_timed "psmc_$mask_id" psmc -N25 -t15 -r5 -p '4+25*2+4+6' \
        -o "$mask_dir/read.unscaled.psmc" "$WORK/$mask_id.psmcfa"
    python3 -m analysis.vgp_read_validation psmc-compare \
        --assembly-psmc "$ASSEMBLY_PSMC" --read-psmc "$mask_dir/read.unscaled.psmc" \
        --mask-id "$mask_id" --output "$mask_dir/psmc_comparison.json"
done

primary="$OUT/masks/dp10_80"
awk 'BEGIN{OFS="\t"} NR>1 {print $1,$2-1,$2}' "$primary/assembly.snps.tsv" > "$WORK/assembly.snps.bed"
run_timed pileup_illumina samtools mpileup -aa -q 20 -Q 20 -l "$WORK/assembly.snps.bed" \
    -f "$WORK/h1.fa" -o "$WORK/illumina.assembly-sites.pileup" "$WORK/illumina.bam"
python3 -m analysis.vgp_read_validation assembly-evidence \
    --assembly-sites "$primary/assembly.snps.tsv" --pileup "$WORK/illumina.assembly-sites.pileup" \
    --minimum-depth 10 --maximum-depth 80 --output "$OUT/sites/illumina.assembly_evidence.tsv" \
    --summary "$OUT/sites/illumina.assembly_evidence.json"
run_timed pileup_hifi samtools mpileup -aa -q 20 -Q 20 -l "$WORK/assembly.snps.bed" \
    -f "$WORK/h1.fa" -o "$WORK/hifi.assembly-sites.pileup" "$WORK/hifi.bam"
python3 -m analysis.vgp_read_validation assembly-evidence \
    --assembly-sites "$primary/assembly.snps.tsv" --pileup "$WORK/hifi.assembly-sites.pileup" \
    --minimum-depth 10 --maximum-depth 120 --output "$OUT/sites/hifi.assembly_evidence.tsv" \
    --summary "$OUT/sites/hifi.assembly_evidence.json"

python3 - "$OUT" "$STARTED_UTC" "$SLURM_JOB_ID" <<'PY'
import hashlib,json,os,platform,shutil,sys,time
from pathlib import Path
out,started,job=Path(sys.argv[1]),sys.argv[2],sys.argv[3]
profile=Path("/moosefs/erikg/vgp/derived/read-validation/environment/profile").resolve()
tools={}
for name in ("python3","minimap2","samtools","bcftools","bedtools","jellyfish","psmc","fq2psmcfa"):
    resolved=shutil.which(name)
    if resolved is None:
        raise SystemExit(f"missing executable during capture: {name}")
    path=Path(resolved).resolve()
    tools[name]={"path":str(path),"sha256":hashlib.sha256(path.read_bytes()).hexdigest()}
result={"schema_version":"vgp-read-validation-execution-v1","task_id":"validate-vgp-pilot-reads",
        "selection_id":"P07","canonical_vgp_root":"/moosefs/erikg/vgp","slurm_job_id":job,
        "started_utc":started,"completed_utc":time.strftime("%Y-%m-%dT%H:%M:%SZ",time.gmtime()),
        "node":platform.node(),"slurm_cpus":os.environ.get("SLURM_CPUS_PER_TASK"),
        "slurm_mem_per_node":os.environ.get("SLURM_MEM_PER_NODE"),"profile":str(profile),"tools":tools,
        "repository_commit":os.popen("git rev-parse HEAD").read().strip(),
        "validation_module_sha256":hashlib.sha256((out/"executed_validation_module.py").read_bytes()).hexdigest(),
        "method_covariance":"assembly/read pi and PSMC share H1 coordinates, assembly structural mask, and overlapping response data"}
(out/"execution.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
PY

cp -- "$TELEMETRY/worker.stdout" "$TELEMETRY/worker.stdout.snapshot"
cp -- "$TELEMETRY/worker.stderr" "$TELEMETRY/worker.stderr.snapshot"
find "$OUT" -type f ! -name output_manifest.tsv ! -name worker.stdout ! -name worker.stderr -print0 | sort -z | \
    xargs -0 sha256sum | sed "s#  $OUT/#\t#" | awk 'BEGIN{OFS="\t";print "sha256","relative_path"}{print $1,$2}' \
    > "$OUT/output_manifest.tsv"

readonly CANONICAL_PARTIAL="$VGP_ROOT/staging/outputs/vgp-read-validation-P07-${SLURM_JOB_ID}.partial"
readonly CANONICAL_FINAL="$VGP_ROOT/derived/read-validation/runs/P07/slurm-${SLURM_JOB_ID}"
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
