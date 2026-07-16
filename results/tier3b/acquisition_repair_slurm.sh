#!/usr/bin/env bash
#SBATCH --job-name=tier3b-acq-repair
#SBATCH --partition=workers
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --export=NONE
#SBATCH --output=/moosefs/erikg/tier3scratch/tier3b-acquisition/repair-logs/acquisition-repair-%j.out
#SBATCH --error=/moosefs/erikg/tier3scratch/tier3b-acquisition/repair-logs/acquisition-repair-%j.err
set -euo pipefail

# Heavy, fail-closed repair of the two approved Ag1000G tuples.  This job only
# creates a versioned candidate directory.  It never changes the frozen
# acquisition manifest or any of the original 1-Mb staged objects.
REPO=${1:?absolute repository root required}
ENVIRONMENT_RECORD=${2:?absolute Guix environment record required}
ROOT=${3:-/moosefs/erikg/tier3scratch/tier3b-acquisition}
REGION=3R:10000000-30999999
REGION_TAG=3R_10000000_30999999
FINAL=$ROOT/repair-$REGION_TAG-v1
PARTIAL=$ROOT/.repair-$REGION_TAG-v1.partial-${SLURM_JOB_ID:?Slurm job ID required}
GUIX_RUN=$REPO/analysis/slurm/guix_job.sh
COMM=/gnu/store/a5i8avx826brw5grn3n4qv40g514505c-coreutils-9.1/bin/comm
export TIER3_SCRATCH_ROOT=$ROOT/environment

run_guix() { "$GUIX_RUN" "$ENVIRONMENT_RECORD" "$@"; }

AO_SAMPLES=$ROOT/ag1000g_phase3/staged/selected.samples.txt
GM_SAMPLES=$ROOT/ag1000g_phase3_gm_coluzzii/staged/selected.samples.txt
FASTA=$ROOT/ag1000g_phase3/staged/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa
FAI=$FASTA.fai
GFF=$ROOT/ag1000g_phase3/source/Anopheles-gambiae-PEST_BASEFEATURES_AgamP4.12.gff3.gz
SITEFILTER=$ROOT/ag1000g_phase3/source/3R_sitefilters.gamb_colu.dt_20200416.vcf.gz

[[ ! -e "$FINAL" ]] || { echo "refusing to replace existing candidate $FINAL" >&2; exit 1; }
rm -rf "$PARTIAL"
mkdir -p "$PARTIAL/ao" "$PARTIAL/gm" "$PARTIAL/cache/ao" "$PARTIAL/cache/gm"

# Identity gates are checked before any network-heavy work.  These are the
# checksum-locked inputs from the approved 1-Mb acquisition rows.
run_guix python3 - "$AO_SAMPLES" "$GM_SAMPLES" "$FASTA" "$FAI" "$GFF" "$SITEFILTER" <<'PY'
import hashlib, sys
from pathlib import Path
expected = (
    "b52524911b3219bbea9ebd08d129dc269b0b645584a90fcecaf29d62cec627fe",
    "d0d72e84a52453933ad1316eb6783eec46fdca0b84dab619dfb3f0b9e07b90ae",
    "19680ed68a6347f59891ecf0ddc9b54f441bd9b71780db95c02c0dddcd809fe7",
    "75b8cef03d29b8a492430b61dfd62c8cf5447246a3b8dae7a9e0129c7c79d4a3",
    "9329b95dc4daff2ee084674bb693bc8e3647dbf84831e758a1c2ce8179632e3b",
    "25f6a1e77750d9defac12f9d74b225cc7af949a2569a16eefd0a8ae797f9900f",
)
paths = tuple(map(Path, sys.argv[1:]))
for path, wanted in zip(paths, expected):
    observed = hashlib.sha256(path.read_bytes()).hexdigest()
    assert observed == wanted, (path, observed, wanted)
for path in paths[:2]:
    samples = [line for line in path.read_text().splitlines() if line]
    assert len(samples) == len(set(samples)) == 20
print("PASS locked cohort/reference/annotation/site-filter identities")
PY

merge_population() {
  local samples=$1 output=$2 cache=$3
  run_guix bash -c '
set -euo pipefail
samples=$1; output=$2; cache=$3; region=$4
args=()
while IFS= read -r sample; do
  [[ -n "$sample" ]] && args+=("https://vo_agam_output.cog.sanger.ac.uk/${sample}.vcf.gz")
done < "$samples"
[[ ${#args[@]} -eq 20 ]]
cd "$cache"
bcftools merge --regions "$region" --merge none --output-type z --output "$output" "${args[@]}"
bcftools index --tbi "$output"
' bash "$samples" "$output" "$cache" "$REGION"
}

AO_VCF=$PARTIAL/ao/ao_luanda_coluzzii.$REGION_TAG.all_sites.vcf.gz
GM_VCF=$PARTIAL/gm/gm_walikunda_coluzzii.$REGION_TAG.all_sites.vcf.gz
merge_population "$AO_SAMPLES" "$AO_VCF" "$PARTIAL/cache/ao"
merge_population "$GM_SAMPLES" "$GM_VCF" "$PARTIAL/cache/gm"

PASSBED=$PARTIAL/sitefilter_pass.$REGION_TAG.bed
run_guix bcftools query -r "$REGION" -f '%CHROM\t%POS0\t%END\n' -i 'FILTER="PASS"' \
  "$SITEFILTER" > "$PASSBED"

make_callable() {
  local vcf=$1 cohort_bed=$2 callable=$3
  run_guix bcftools query -f '%CHROM\t%POS0\t%END\n' \
    -i 'N_PASS(FMT/GT!="mis" && FMT/DP>=5 && FMT/GQ>=20)>=18' "$vcf" > "$cohort_bed"
  run_guix bash -c '"$4" -12 "$1" "$2" | bedtools merge -i - > "$3"' \
    bash "$PASSBED" "$cohort_bed" "$callable" "$COMM"
}

AO_COHORT=$PARTIAL/ao/cohort_dp5_gq20_ge18.$REGION_TAG.bed
GM_COHORT=$PARTIAL/gm/cohort_dp5_gq20_ge18.$REGION_TAG.bed
AO_CALL=$PARTIAL/ao/ao_luanda_coluzzii.$REGION_TAG.callable.bed
GM_CALL=$PARTIAL/gm/gm_walikunda_coluzzii.$REGION_TAG.callable.bed
make_callable "$AO_VCF" "$AO_COHORT" "$AO_CALL"
make_callable "$GM_VCF" "$GM_COHORT" "$GM_CALL"

for tuple in ao gm; do
  if [[ $tuple == ao ]]; then
    vcf=$AO_VCF; samples=$AO_SAMPLES; callable=$AO_CALL
  else
    vcf=$GM_VCF; samples=$GM_SAMPLES; callable=$GM_CALL
  fi
  run_guix bgzip --test "$vcf"
  run_guix bcftools index --stats "$vcf" > "$PARTIAL/$tuple/index.stats.tsv"
  run_guix bcftools norm --check-ref exit --fasta-ref "$FASTA" --regions "$REGION" \
    --output-type u "$vcf" >/dev/null 2>"$PARTIAL/$tuple/norm.stderr"
  run_guix bcftools query -l "$vcf" > "$PARTIAL/$tuple/observed.samples.txt"
  run_guix python3 - "$samples" "$PARTIAL/$tuple/observed.samples.txt" "$callable" <<'PY'
import sys
from pathlib import Path
expected = [x for x in Path(sys.argv[1]).read_text().splitlines() if x]
observed = [x for x in Path(sys.argv[2]).read_text().splitlines() if x]
assert observed == expected and len(observed) == len(set(observed)) == 20
total = 0
previous = None
blocks = set()
for number, line in enumerate(Path(sys.argv[3]).read_text().splitlines(), 1):
    chrom, start, end = line.split("\t")[:3]
    start, end = int(start), int(end)
    assert chrom == "3R" and 9_999_999 <= start < end <= 30_999_999
    if previous is not None:
        assert start > previous, (number, start, previous)
    previous = end
    total += end - start
    blocks.update(range(start // 1_000_000, (end - 1) // 1_000_000 + 1))
assert total > 0 and len(blocks) >= 20, (total, len(blocks))
print("PASS exact sample order and callable coordinate blocks", total, len(blocks))
PY
done

# Summarize record, nonreference, callability, index, checksum, and coordinate
# gates without trusting filenames.  Frozen 4D power is checked separately by
# the post-acquisition validation job before the manifest is changed.
run_guix python3 - "$PARTIAL" "$FINAL" "$REGION" "$FASTA" "$FAI" "$GFF" "$SITEFILTER" \
  "$AO_SAMPLES" "$GM_SAMPLES" "$AO_VCF" "$GM_VCF" "$AO_CALL" "$GM_CALL" <<'PY'
import hashlib, json, subprocess, sys
from pathlib import Path
out = Path(sys.argv[1]); published = Path(sys.argv[2]); region = sys.argv[3]
fasta, fai, gff, sitefilter, ao_samples, gm_samples, ao_vcf, gm_vcf, ao_call, gm_call = map(Path, sys.argv[4:])
def sha(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()
def lines(command):
    return subprocess.check_output(command, text=True).splitlines()
def count_output_lines(command):
    process = subprocess.Popen(command, stdout=subprocess.PIPE, text=True)
    assert process.stdout is not None
    count = sum(1 for _ in process.stdout)
    return_code = process.wait()
    assert return_code == 0, (command, return_code)
    return count
def callable_bases(path):
    return sum(int(fields[2])-int(fields[1]) for fields in (line.split("\t") for line in path.read_text().splitlines()))
summary = {"schema_version":"tier3b-acquisition-repair-v1", "region":region,
           "slurm_job_id":__import__("os").environ["SLURM_JOB_ID"],
           "locked_inputs":{str(path):sha(path) for path in (fasta,fai,gff,sitefilter,ao_samples,gm_samples)},
           "tuples":{}}
for name, vcf, call, samples in (("ao",ao_vcf,ao_call,ao_samples),("gm",gm_vcf,gm_call,gm_samples)):
    stats = lines(["bcftools","index","--stats",str(vcf)])
    fields = next(line.split("\t") for line in stats if line.startswith("3R\t"))
    # Exact nonreference-genotype records are counted in a streaming query.
    nonref = count_output_lines(["bcftools","query","-r",region,"-f","%POS\\n","-i","N_PASS(FMT/GT!=\"mis\" && FMT/GT!=\"RR\")>0",str(vcf)])
    summary["tuples"][name] = {
        "sample_count":len([x for x in samples.read_text().splitlines() if x]),
        "record_count":int(fields[2]),
        "nonreference_genotype_record_count":nonref, "callable_sites":callable_bases(call),
        "callset_path":str(published / vcf.relative_to(out)), "callset_sha256":sha(vcf),
        "callset_index_path":str(published / Path(str(vcf)+".tbi").relative_to(out)), "callset_index_sha256":sha(Path(str(vcf)+".tbi")),
        "callable_mask_path":str(published / call.relative_to(out)), "callable_mask_sha256":sha(call),
    }
    assert summary["tuples"][name]["record_count"] > 0
    assert nonref > 0 and summary["tuples"][name]["callable_sites"] > 0
(out/"acquisition_base_qc.json").write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n")
PY

# Cache files are reproducible implementation details, not tuple inputs.
rm -rf "$PARTIAL/cache"
mv "$PARTIAL" "$FINAL"
echo "PASS staged repair candidate: $FINAL"
