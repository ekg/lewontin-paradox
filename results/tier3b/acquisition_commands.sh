#!/usr/bin/env bash
set -euo pipefail

# Reproduce the two approved Ag1000G biological tuples and their frozen-power
# repair. Run from repository root on a Slurm login node.
# Large objects remain under MooseFS and are never added to Git.
ROOT=${TIER3B_ACQUISITION_ROOT:-/moosefs/erikg/tier3scratch/tier3b-acquisition}
ENV_RECORD=$PWD/analysis/pilot_results/guix_environment.json
GUIX_RUN=$PWD/analysis/slurm/guix_job.sh
CURL=/gnu/store/b727ryyfiz1cfdywjp8s1wmxd6lzsz8p-curl-7.85.0/bin/curl
COMM=/gnu/store/a5i8avx826brw5grn3n4qv40g514505c-coreutils-9.1/bin/comm
export TIER3_SCRATCH_ROOT=$ROOT/environment

run_guix() { "$GUIX_RUN" "$ENV_RECORD" "$@"; }

run_guix python3 analysis/run_tier3.py audit-environment "$ENV_RECORD"
run_guix bcftools --version
run_guix samtools --version
run_guix bedtools --version
run_guix python3 --version

mkdir -p "$ROOT/ag1000g_phase3/source" "$ROOT/ag1000g_phase3/staged" \
  "$ROOT/ag1000g_phase3_gm_coluzzii/staged" "$ROOT/environment"

run_guix bash -c '
set -euo pipefail
curl_bin=$1; shift
while (( $# )); do
  url=$1; output=$2; shift 2
  "$curl_bin" --fail --location --retry 5 --retry-delay 2 --continue-at - --output "$output" "$url"
done
' bash "$CURL" \
  https://storage.googleapis.com/vo_anoph_temp_us_central1/vo_agam_release/reference/genome/agamp4/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa.gz "$ROOT/ag1000g_phase3/source/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa.gz" \
  https://storage.googleapis.com/vo_anoph_temp_us_central1/vo_agam_release/reference/genome/agamp4/Anopheles-gambiae-PEST_BASEFEATURES_AgamP4.12.gff3.gz "$ROOT/ag1000g_phase3/source/Anopheles-gambiae-PEST_BASEFEATURES_AgamP4.12.gff3.gz" \
  https://storage.googleapis.com/vo_anoph_temp_us_central1/vo_agam_release/v3/site_filters/dt_20200416/vcf/gamb_colu/3R_sitefilters.vcf.gz "$ROOT/ag1000g_phase3/source/3R_sitefilters.gamb_colu.dt_20200416.vcf.gz" \
  https://storage.googleapis.com/vo_anoph_temp_us_central1/vo_agam_release/v3/site_filters/dt_20200416/vcf/gamb_colu/3R_sitefilters.vcf.gz.tbi "$ROOT/ag1000g_phase3/source/3R_sitefilters.gamb_colu.dt_20200416.vcf.gz.tbi" \
  https://raw.githubusercontent.com/malariagen/ag1000g-phase3-data-paper/master/content/tables/per_sample_location_table.csv "$ROOT/ag1000g_phase3/source/per_sample_location_table.csv" \
  https://raw.githubusercontent.com/malariagen/ag1000g-phase3-data-paper/master/content/tables/per_sample_location_table_README.txt "$ROOT/ag1000g_phase3/source/per_sample_location_table_README.txt"

run_guix python3 - "$ROOT/ag1000g_phase3/source/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa.gz" "$ROOT/ag1000g_phase3/staged/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa" <<'PY'
import gzip, shutil, sys
with gzip.open(sys.argv[1], "rb") as src, open(sys.argv[2], "wb") as dst:
    shutil.copyfileobj(src, dst, 16 * 1024 * 1024)
PY
run_guix samtools faidx "$ROOT/ag1000g_phase3/staged/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa"

select_population() {
  local dataset=$1 sample_set=$2 site=$3 year=$4 population=$5 meta=$6 selected=$7
  run_guix python3 - "$ROOT/ag1000g_phase3/source/per_sample_location_table.csv" "$meta" "$selected" "$dataset" "$sample_set" "$site" "$year" "$population" <<'PY'
import csv, sys
from pathlib import Path
from analysis.tier3b_popvcf_collect import deterministic_sample_selection
src, meta, selected = map(Path, sys.argv[1:4])
dataset, sample_set, site, year, population = sys.argv[4:]
rows=[]
with src.open(newline="", encoding="utf-8") as handle:
    for row in csv.DictReader(handle):
        if row["sample_set"]==sample_set and row["site"]==site and row["year"]==year and row["is_coluzzii"]=="True":
            rows.append({"sample_id":row["sample_id"], "population_id":population,
              "qc_pass":"true", "unrelated":"true", "cross_or_progeny":"false",
              "duplicate_or_related":"false", "laboratory_control":"false", "contaminated":"false",
              "partner_sample_id":row["partner_sample_id"], "sex_call":row["sex_call"],
              "aim_fraction_colu":row["aim_fraction_colu"]})
with meta.open("w", newline="", encoding="utf-8") as handle:
    writer=csv.DictWriter(handle, fieldnames=list(rows[0]), dialect="excel-tab")
    writer.writeheader(); writer.writerows(rows)
selection=deterministic_sample_selection(dataset, rows, population)
selected.write_text("".join(sample+"\n" for sample in selection["selected_samples"]), encoding="utf-8")
PY
}

select_population ag1000g_phase3_AO_coluzzii AG1000G-AO Luanda 2009 AO_Luanda_2009_coluzzii \
  "$ROOT/ag1000g_phase3/staged/sample_metadata.tsv" "$ROOT/ag1000g_phase3/staged/selected.samples.txt"
select_population ag1000g_phase3_GM_coluzzii AG1000G-GM-C 'Wali Kunda' 2012 GM_WaliKunda_2012_coluzzii \
  "$ROOT/ag1000g_phase3_gm_coluzzii/staged/sample_metadata.tsv" "$ROOT/ag1000g_phase3_gm_coluzzii/staged/selected.samples.txt"

merge_population() {
  local samples=$1 output=$2 cache=$3
  mkdir -p "$cache"
  run_guix bash -c '
set -euo pipefail
args=()
while IFS= read -r sample; do args+=("https://vo_agam_output.cog.sanger.ac.uk/${sample}.vcf.gz"); done < "$1"
cd "$3"
bcftools merge --regions 3R:10000000-10999999 --merge none --output-type z --output "$2" "${args[@]}"
bcftools index --tbi "$2"
' bash "$samples" "$output" "$cache"
}

AO_VCF=$ROOT/ag1000g_phase3/staged/ao_luanda_coluzzii.3R_10000000_10999999.all_sites.vcf.gz
GM_VCF=$ROOT/ag1000g_phase3_gm_coluzzii/staged/gm_walikunda_coluzzii.3R_10000000_10999999.all_sites.vcf.gz
merge_population "$ROOT/ag1000g_phase3/staged/selected.samples.txt" "$AO_VCF" "$ROOT/ag1000g_phase3/source/remote_index_cache"
merge_population "$ROOT/ag1000g_phase3_gm_coluzzii/staged/selected.samples.txt" "$GM_VCF" "$ROOT/ag1000g_phase3_gm_coluzzii/source/remote_index_cache"

PASSBED=$ROOT/ag1000g_phase3/staged/sitefilter_pass.3R_10000000_10999999.bed
run_guix bcftools query -r 3R:10000000-10999999 -f '%CHROM\t%POS0\t%END\n' -i 'FILTER="PASS"' \
  "$ROOT/ag1000g_phase3/source/3R_sitefilters.gamb_colu.dt_20200416.vcf.gz" > "$PASSBED"

make_callable() {
  local vcf=$1 cohort_bed=$2 callable=$3
  run_guix bcftools query -f '%CHROM\t%POS0\t%END\n' \
    -i 'N_PASS(FMT/GT!="mis" && FMT/DP>=5 && FMT/GQ>=20)>=18' "$vcf" > "$cohort_bed"
  run_guix bash -c '"$4" -12 "$1" "$2" | bedtools merge -i - > "$3"' \
    bash "$PASSBED" "$cohort_bed" "$callable" "$COMM"
}

AO_CALL=$ROOT/ag1000g_phase3/staged/ao_luanda_coluzzii.3R_10000000_10999999.callable.bed
GM_CALL=$ROOT/ag1000g_phase3_gm_coluzzii/staged/gm_walikunda_coluzzii.3R_10000000_10999999.callable.bed
make_callable "$AO_VCF" "$ROOT/ag1000g_phase3/staged/cohort_dp5_gq20_ge18.3R_10000000_10999999.bed" "$AO_CALL"
make_callable "$GM_VCF" "$ROOT/ag1000g_phase3_gm_coluzzii/staged/cohort_dp5_gq20_ge18.3R_10000000_10999999.bed" "$GM_CALL"

FASTA=$ROOT/ag1000g_phase3/staged/Anopheles-gambiae-PEST_CHROMOSOMES_AgamP4.fa
for vcf in "$AO_VCF" "$GM_VCF"; do
  run_guix bgzip --test "$vcf"
  run_guix bcftools index --stats "$vcf"
  run_guix bcftools norm --check-ref exit --fasta-ref "$FASTA" --regions 3R:10000000-10999999 --output-type u "$vcf" >/dev/null
done

smoke() {
  local dataset=$1 full=$2 samples=$3 callable=$4 prefix=$5
  run_guix bcftools view -r 3R:10000000-10099999 -Oz -o "$prefix.vcf.gz" "$full"
  run_guix bcftools index --tbi "$prefix.vcf.gz"
  run_guix python3 - "$callable" "$prefix.callable.bed" <<'PY'
from pathlib import Path
import sys
lo, hi = 9999999, 10099999
with Path(sys.argv[1]).open() as src, Path(sys.argv[2]).open("w") as dst:
    for line in src:
        chrom,start,end=line.rstrip().split("\t")[:3]; start,end=int(start),int(end)
        if chrom=="3R" and max(start,lo)<min(end,hi):
            dst.write(f"{chrom}\t{max(start,lo)}\t{min(end,hi)}\n")
PY
  run_guix python3 analysis/tier3b_popvcf_compute.py --dataset-id "$dataset" --vcf "$prefix.vcf.gz" \
    --fasta "$FASTA" --selected-samples "$samples" --design wild_diploid \
    --denominator-kind cohort_callable_mask --callable-bed "$prefix.callable.bed" --output "$prefix.json"
}

smoke ag1000g_phase3_ao_smoke "$AO_VCF" "$ROOT/ag1000g_phase3/staged/selected.samples.txt" "$AO_CALL" \
  "$ROOT/ag1000g_phase3/staged/ao_luanda_coluzzii.3R_10000000_10099999.smoke"
smoke ag1000g_phase3_gm_smoke "$GM_VCF" "$ROOT/ag1000g_phase3_gm_coluzzii/staged/selected.samples.txt" "$GM_CALL" \
  "$ROOT/ag1000g_phase3_gm_coluzzii/staged/gm_walikunda_coluzzii.3R_10000000_10099999.smoke"

# Reproduce the higher-ranked DGRP attempt. It is preserved but deliberately
# excluded from acquisition_manifest.tsv because the untouched native r5.57
# GFF sequence-region bounds fail the exact FASTA dictionary gate below.
mkdir -p "$ROOT/dgrp_dgn1.1/source" "$ROOT/dgrp_dgn1.1/staged"
run_guix bash -c '
set -euo pipefail
curl_bin=$1; shift
while (( $# )); do
  url=$1; output=$2; shift 2
  "$curl_bin" --fail --location --retry 5 --retry-delay 2 --continue-at - --output "$output" "$url"
done
' bash "$CURL" \
  https://pooldata.genetics.wisc.edu/dgrp_sequences.tar.bz2 "$ROOT/dgrp_dgn1.1/source/dgrp_sequences.tar.bz2" \
  https://s3ftp.flybase.org/genomes/Drosophila_melanogaster/dmel_r5.57_FB2014_03/fasta/dmel-all-chromosome-r5.57.fasta.gz "$ROOT/dgrp_dgn1.1/source/dmel-all-chromosome-r5.57.fasta.gz" \
  https://s3ftp.flybase.org/genomes/Drosophila_melanogaster/dmel_r5.57_FB2014_03/gff/dmel-all-r5.57.gff.gz "$ROOT/dgrp_dgn1.1/source/dmel-all-r5.57.gff.gz" \
  https://johnpool.net/TableS1_individuals.xls "$ROOT/dgrp_dgn1.1/source/TableS1_individuals.xls" \
  https://johnpool.net/TableS2_populations.xls "$ROOT/dgrp_dgn1.1/source/TableS2_populations.xls" \
  https://johnpool.net/masking.zip "$ROOT/dgrp_dgn1.1/source/masking.zip"

run_guix python3 - "$ROOT/dgrp_dgn1.1/source/dgrp_sequences.tar.bz2" <<'PY'
import hashlib, sys
md5=hashlib.md5(); sha=hashlib.sha256()
with open(sys.argv[1], "rb") as handle:
    for block in iter(lambda: handle.read(16*1024*1024), b""):
        md5.update(block); sha.update(block)
assert md5.hexdigest()=="c697730b0720e944ab2be32e391322b0"
assert sha.hexdigest()=="65c63073f14ea02db0421983b5bdd9e8c6434cd8aa1cc5d893d42e6413f02528"
PY

run_guix python3 - "$ROOT/dgrp_dgn1.1/source/dgrp_sequences.tar.bz2" "$ROOT/dgrp_dgn1.1/source/dgrp_Chr2L.tar" <<'PY'
from pathlib import Path
import sys, tarfile
src,out=map(Path,sys.argv[1:])
with tarfile.open(src,"r|bz2") as archive:
    for member in archive:
        if member.name=="dgrp_Chr2L.tar":
            assert member.isfile() and member.size==4717533184
            with archive.extractfile(member) as reader, out.open("wb") as writer:
                for block in iter(lambda: reader.read(16*1024*1024), b""):
                    writer.write(block)
            break
    else: raise SystemExit("DGRP archive lacks dgrp_Chr2L.tar")
PY

run_guix python3 - "$ROOT/dgrp_dgn1.1/source/dmel-all-chromosome-r5.57.fasta.gz" "$ROOT/dgrp_dgn1.1/staged/dmel-all-chromosome-r5.57.fasta" <<'PY'
import gzip, shutil, sys
with gzip.open(sys.argv[1],"rb") as src, open(sys.argv[2],"wb") as dst:
    shutil.copyfileobj(src,dst,16*1024*1024)
PY
run_guix samtools faidx "$ROOT/dgrp_dgn1.1/staged/dmel-all-chromosome-r5.57.fasta"

run_guix python3 - "$ROOT/dgrp_dgn1.1/source/TableS1_individuals.xls" "$ROOT/dgrp_dgn1.1/staged/sample_metadata.tsv" "$ROOT/dgrp_dgn1.1/staged/selected.samples.txt" <<'PY'
import csv,sys
from pathlib import Path
import pandas as pd
from analysis.tier3b_popvcf_collect import deterministic_sample_selection
src,meta,selected=map(Path,sys.argv[1:])
frame=pd.read_excel(src,header=6)
rows=[]
for _,record in frame[frame["Data Group"].astype(str).str.upper()=="DGRP"].iterrows():
    sample=str(record["Stock ID"]).strip()
    represented=str(record["Focal Genome Represented"])
    coverage=float(record["Mb Genomic Coverage (homozygous calls on focal arms)"])
    depth=float(record["Mean Depth"]); inversion=str(record["In(2L)t"]).strip()
    qc=coverage>=100 and depth>=15 and "2L" in represented.split(",") and inversion=="ST"
    rows.append({"sample_id":sample,"population_id":"RAL_Raleigh_NC_USA","qc_pass":str(qc).lower(),
      "unrelated":"true","cross_or_progeny":"false","duplicate_or_related":"false",
      "laboratory_control":"false","contaminated":"false","sra_accession":str(record["SRA Accession"]).strip(),
      "genome_type":str(record["Genome Type"]).strip(),"focal_genome_represented":represented,
      "genomic_coverage_mb":f"{coverage:.6f}","mean_depth":f"{depth:.1f}","In_2L_t":inversion})
with meta.open("w",newline="",encoding="utf-8") as handle:
    writer=csv.DictWriter(handle,fieldnames=list(rows[0]),dialect="excel-tab")
    writer.writeheader(); writer.writerows(rows)
selection=deterministic_sample_selection("dgrp_dgn1.1_ral_2L_standard",rows,"RAL_Raleigh_NC_USA")
assert selection["population_eligible_units"]==67
selected.write_text("".join(x+"\n" for x in selection["selected_samples"]),encoding="utf-8")
PY

DGRP_PLAIN=$ROOT/dgrp_dgn1.1/staged/dgrp_ral.2L_5000000_5999999.all_sites.vcf
DGRP_BED=$ROOT/dgrp_dgn1.1/staged/dgrp_ral.2L_5000000_5999999.callable.bed
run_guix python3 - "$ROOT/dgrp_dgn1.1/source/dgrp_Chr2L.tar" "$ROOT/dgrp_dgn1.1/staged/selected.samples.txt" \
  "$ROOT/dgrp_dgn1.1/staged/dmel-all-chromosome-r5.57.fasta" "$ROOT/dgrp_dgn1.1/source/masking.zip" \
  "$DGRP_PLAIN" "$DGRP_BED" <<'PY'
from pathlib import Path
import sys,tarfile,zipfile
from analysis.tier3_common import read_fasta,fasta_dictionary
inner,selected_path,fasta_path,masking_zip,vcf_path,bed_path=map(Path,sys.argv[1:])
selected=[x for x in selected_path.read_text().splitlines() if x]
assert len(selected)==len(set(selected))==20
fasta=read_fasta(fasta_path); ref=fasta["2L"]; dictionaries=fasta_dictionary(fasta)
sequences={}
with tarfile.open(inner,"r:") as archive:
    for sample in selected:
        seq=archive.extractfile(archive.getmember(f"{sample}_Chr2L.seq")).read().decode("ascii").strip().upper()
        assert len(seq)==len(ref) and not set(seq)-set("ACGTN")
        sequences[sample]=seq
masks={sample:[] for sample in selected}
with zipfile.ZipFile(masking_zip) as archive:
    for filename,ibd in (("ibd_filter_tracts.txt",True),("admixture_filter_tracts.txt",False)):
        for raw in archive.read(filename).decode().splitlines():
            fields=raw.split()
            if len(fields)<4 or fields[0] not in masks or fields[1]!="Chr2L": continue
            start,stop=map(int,fields[2:4])
            lo,hi=(start-100001,stop+100002) if ibd else (start-1,stop+2)
            masks[fields[0]].append((max(0,lo),min(len(ref),hi)))
def masked(sample,pos): return any(lo<=pos<hi for lo,hi in masks[sample])
lo,hi=4_999_999,5_999_999; records=variants=callable=0; interval=None
with vcf_path.open("w",encoding="ascii") as out,bed_path.open("w",encoding="ascii") as bed:
    out.write("##fileformat=VCFv4.2\n##source=DGN1.1_DGRP_SEQ_consensus_plus_provider_IBD_admixture_masks\n")
    out.write("##reference=FlyBase_FB2014_03_Dmel_release_5.57\n")
    out.write('##FILTER=<ID=PASS,Description="DGN consensus site; line filters represented as missing GT">\n')
    out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Haploidized inbred-line consensus genotype">\n')
    for contig,length in dictionaries.items(): out.write(f"##contig=<ID={contig},length={length}>\n")
    out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"+"\t".join(selected)+"\n")
    for pos in range(lo,hi):
        reference=ref[pos]
        calls=[None if sequences[sample][pos] not in "ACGT" or masked(sample,pos) else sequences[sample][pos] for sample in selected]
        is_callable=reference in "ACGT" and sum(base is not None for base in calls)>=18
        if is_callable:
            callable+=1
            if interval is None: interval=pos
        elif interval is not None:
            bed.write(f"2L\t{interval}\t{pos}\n"); interval=None
        if reference not in "ACGT": continue
        alts=sorted({base for base in calls if base is not None and base!=reference})
        allele={reference:"0",**{base:str(i+1) for i,base in enumerate(alts)}}
        genotypes=["." if base is None else allele[base] for base in calls]
        out.write(f"2L\t{pos+1}\t.\t{reference}\t{','.join(alts) if alts else '.'}\t.\tPASS\t.\tGT\t"+"\t".join(genotypes)+"\n")
        records+=1; variants+=bool(alts)
    if interval is not None: bed.write(f"2L\t{interval}\t{hi}\n")
assert (records,variants,callable)==(1000000,25652,831671)
PY
run_guix bgzip --force "$DGRP_PLAIN"
run_guix tabix --preset vcf "$DGRP_PLAIN.gz"
run_guix bgzip --test "$DGRP_PLAIN.gz"
run_guix bcftools index --stats "$DGRP_PLAIN.gz"
run_guix bcftools norm --check-ref exit --fasta-ref "$ROOT/dgrp_dgn1.1/staged/dmel-all-chromosome-r5.57.fasta" \
  --regions 2L:5000000-5999999 --output-type u "$DGRP_PLAIN.gz" >/dev/null

run_guix bcftools view -r 2L:5000000-5099999 -Oz -o "$ROOT/dgrp_dgn1.1/staged/dgrp_ral.2L_5000000_5099999.smoke.vcf.gz" "$DGRP_PLAIN.gz"
run_guix bcftools index --tbi "$ROOT/dgrp_dgn1.1/staged/dgrp_ral.2L_5000000_5099999.smoke.vcf.gz"
run_guix python3 - "$DGRP_BED" "$ROOT/dgrp_dgn1.1/staged/dgrp_ral.2L_5000000_5099999.smoke.callable.bed" <<'PY'
from pathlib import Path
import sys
lo,hi=4_999_999,5_099_999
with Path(sys.argv[1]).open() as src,Path(sys.argv[2]).open("w") as dst:
    for line in src:
        chrom,start,end=line.rstrip().split("\t")[:3]; start,end=int(start),int(end)
        if chrom=="2L" and max(start,lo)<min(end,hi):
            dst.write(f"{chrom}\t{max(start,lo)}\t{min(end,hi)}\n")
PY
run_guix python3 analysis/tier3b_popvcf_compute.py --dataset-id dgrp_dgn1.1_ral_smoke \
  --vcf "$ROOT/dgrp_dgn1.1/staged/dgrp_ral.2L_5000000_5099999.smoke.vcf.gz" \
  --fasta "$ROOT/dgrp_dgn1.1/staged/dmel-all-chromosome-r5.57.fasta" \
  --selected-samples "$ROOT/dgrp_dgn1.1/staged/selected.samples.txt" --design inbred_lines_haploidized \
  --denominator-kind cohort_callable_mask \
  --callable-bed "$ROOT/dgrp_dgn1.1/staged/dgrp_ral.2L_5000000_5099999.smoke.callable.bed" \
  --output "$ROOT/dgrp_dgn1.1/staged/dgrp_ral.2L_5000000_5099999.smoke.json"

run_guix python3 - "$ROOT/dgrp_dgn1.1/source/dmel-all-r5.57.gff.gz" "$ROOT/dgrp_dgn1.1/staged/dmel-all-chromosome-r5.57.fasta" <<'PY'
import gzip,sys
from analysis.tier3_common import fasta_dictionary,read_fasta
reference=fasta_dictionary(read_fasta(sys.argv[2])); declared={}
with gzip.open(sys.argv[1],"rt",encoding="utf-8") as handle:
    for line in handle:
        if line.startswith("##sequence-region "):
            _,_,contig,start,end=line.split()
            declared[contig]=int(end)-int(start)+1
mismatches={key:(reference.get(key),value) for key,value in declared.items() if reference.get(key)!=value}
assert len(declared)==15 and len(mismatches)==15
assert {observed-reference_length for key,(reference_length,observed) in mismatches.items()}=={1,2}
print("DGRP disqualified: native GFF/FASTA sequence-region mismatches",mismatches)
PY

# Reproduce the fail-before-repair evidence and the versioned 21-Mb repair on
# Slurm. The repair job refuses to replace an existing candidate directory;
# use a fresh TIER3B_ACQUISITION_ROOT for a clean end-to-end reproduction.
mkdir -p "$ROOT/repair-logs"
CURRENT_PREFLIGHT_JOB=$(sbatch --parsable --wait \
  results/tier3b/acquisition_repair_validation_slurm.sh "$PWD" "$ENV_RECORD" "$ROOT" current "$ROOT" \
  "$ROOT/repair-logs/current-preflight.json")
REPAIR_ACQUISITION_JOB=$(sbatch --parsable --wait \
  results/tier3b/acquisition_repair_slurm.sh "$PWD" "$ENV_RECORD" "$ROOT")
REPAIR_ROOT=$ROOT/repair-3R_10000000_30999999-v1
REPAIR_POWER_JOB=$(sbatch --parsable --wait \
  results/tier3b/acquisition_repair_validation_slurm.sh "$PWD" "$ENV_RECORD" "$ROOT" repair "$REPAIR_ROOT" \
  "$REPAIR_ROOT/acquisition_power_qc.json")
sacct -j "$CURRENT_PREFLIGHT_JOB,$REPAIR_ACQUISITION_JOB,$REPAIR_POWER_JOB" \
  --format=JobIDRaw,JobName%28,State,Elapsed,TotalCPU,MaxRSS,ReqMem,AllocCPUS,ExitCode,NodeList \
  -n -P > "$REPAIR_ROOT/slurm_sacct.tsv"

# Executed repair provenance (2026-07-16): current preflight 1761600,
# acquisition 1761599, repaired power validation 1761601. All completed on
# octopus07 with exit 0:0; exact sacct output and checksums are recorded in QC.
run_guix python3 - results/tier3b/acquisition_manifest.tsv "$REPAIR_ROOT/acquisition_power_qc.json" <<'PY'
import csv, hashlib, json, math, sys
from pathlib import Path
with Path(sys.argv[1]).open(newline="", encoding="utf-8") as handle:
    rows=list(csv.DictReader(handle, dialect="excel-tab"))
assert len(rows)>=2 and all(r["status"]=="approved" and r["biological"]=="true" for r in rows)
for row in rows:
    for path_key, hash_key in (("reference_path","reference_sha256"),("reference_fai_path","reference_fai_sha256"),
                               ("native_annotation_path","native_annotation_sha256"),("callset_path","callset_sha256"),
                               ("callset_index_path","callset_index_sha256"),("callable_mask_path","callable_mask_sha256"),
                               ("callable_source_path","callable_source_sha256"),("sample_list_path","sample_list_sha256"),
                               ("sample_metadata_path","sample_metadata_sha256")):
        digest=hashlib.sha256(Path(row[path_key]).read_bytes()).hexdigest()
        assert digest==row[hash_key], (row["tuple_id"], path_key)
    assert row["region"]=="3R:10000000-30999999" and int(row["sample_count"])==20
    assert int(row["record_count"])==20817962
    assert int(row["nonreference_genotype_record_count"])>0 and int(row["callable_sites"])>0
power=json.loads(Path(sys.argv[2]).read_text())
assert power["status"]=="PASS" and len(power["tuples"])==2
for result in power["tuples"]:
    assert min(result["exact_callable_fourfold"].values())>=10000
    assert result["eligible_ratio_blocks"]>=20 and result["status"]=="PASS"
print("validated approved tuples", len(rows))
PY
