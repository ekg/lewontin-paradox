# Tier 3A biological diploid QC

## Outcome

All 3 acquired biological tuples produced coding-region estimates through SweepGA 1:1 and IMPG.
Every row is a coding/CDS annotation-panel estimate, not a genome-wide estimate. The panel was selected by stable hash of exact H1-native execution-span identity before mapping or variant inspection.

## Separation of responsibilities

- SweepGA supplied complete whole-H1-versus-H2 mappings and enforced the 1:1 query:target overlap cap. Fixed-point cap validation remains in the acquisition handoff.
- IMPG indexed those complete bounded PAFs, formed native graph partitions (`-w 2000 -d 0`), and queried only partitions intersecting selected native-annotation spans.
- `bcftools norm` normalized and split alleles, panel trimming removed partition context, and exact duplicate records were removed before estimates.

## Tuple audits

### Menidia menidia (`menidia_menidia_fMenMen1`)

Targeted genes/bases: 22,548/364,299,923; deterministic panel: 157/2,038,234; SweepGA-1:1 mapped panel genes/bases: 157/2,038,234.
Across all native targets, SweepGA 1:1 fully covered 20,953 genes and excluded or only partially covered 1,595; its target-gene intersection retained 345,129,609 bases and lost 19,170,314.
IMPG-audited callable coding/CDS bases: 2,038,234/302,933. Representative SNV: H1 CM109610.1:6022075 T versus H2 CM108917.1:5614205 G (strand +); both allele checks passed.
IMPG selected 1,165 native partitions and emitted 1,159 regional VCFs. Normalization yielded 42,509 records before exact panel trim/dedup and 37,353 after; 37,351 overlap callable coding genes.
Cap coverage sensitivity (query,target): 1:1=(0.8593,0.8310), 5:5=(0.8606,0.8311), 10:10=(0.8607,0.8311). Caps 5/10 are coverage sensitivity only and were not passed to IMPG as graph policy.

### Spinachia spinachia (`spinachia_spinachia_SK-2024b`)

Targeted genes/bases: 20,772/275,623,286; deterministic panel: 150/2,018,006; SweepGA-1:1 mapped panel genes/bases: 150/2,018,006.
Across all native targets, SweepGA 1:1 fully covered 20,491 genes and excluded or only partially covered 281; its target-gene intersection retained 272,591,781 bases and lost 3,031,505.
IMPG-audited callable coding/CDS bases: 2,017,806/270,674. Representative SNV: H1 CM106587.1:2892985 C versus H2 CM106669.1:2760454 T (strand +); both allele checks passed.
IMPG selected 1,126 native partitions and emitted 1,122 regional VCFs. Normalization yielded 1,696 records before exact panel trim/dedup and 1,430 after; 1,429 overlap callable coding genes.
Cap coverage sensitivity (query,target): 1:1=(0.9438,0.9394), 5:5=(0.9445,0.9394), 10:10=(0.9445,0.9394). Caps 5/10 are coverage sensitivity only and were not passed to IMPG as graph policy.

### Tautogolabrus adspersus (`tautogolabrus_adspersus_fTauAds1`)

Targeted genes/bases: 21,614/441,517,395; deterministic panel: 123/2,116,227; SweepGA-1:1 mapped panel genes/bases: 123/2,116,227.
Across all native targets, SweepGA 1:1 fully covered 20,167 genes and excluded or only partially covered 1,447; its target-gene intersection retained 424,072,764 bases and lost 17,444,631.
IMPG-audited callable coding/CDS bases: 2,116,227/222,887. Representative SNV: H1 CM036534.1:34981587 C versus H2 JAJGRG010000537.1:207731 A (strand +); both allele checks passed.
IMPG selected 1,164 native partitions and emitted 1,160 regional VCFs. Normalization yielded 5,280 records before exact panel trim/dedup and 4,791 after; 4,785 overlap callable coding genes.
Cap coverage sensitivity (query,target): 1:1=(0.9498,0.9102), 5:5=(0.9502,0.9102), 10:10=(0.9502,0.9102). Caps 5/10 are coverage sensitivity only and were not passed to IMPG as graph policy.

## Boundary, annotation, and uncertainty policy

Coordinates remain zero-based half-open through native feature, SweepGA coverage, and IMPG-partition intersections. Only exact H1 A/C/G/T bases enter denominators. CDS feature row, gene, transcript, protein, locus, strand, and phase identities are retained in each tuple's `variant_feature_audit.tsv`. Fourfold sites require a valid transcript-order phase chain and exclude frame-discordant overlaps.
Uncertainty is a deterministic genomic block bootstrap (50-kb blocks, 1,000 replicates). Direct CIGAR traversal is used only for one representative H1/H2 allele audit per tuple; it does not create the primary call set.

## Reproducibility

The run manifest records Slurm telemetry, exact biological paths, pinned Guix channel commit/profile, and primary artifacts. `diploid_rerun_commands.sh` contains exact commands. Scheduler stdout/stderr and `sacct` telemetry are under `results/tier3a/logs/`.

Slurm measured successful elapsed times of 4:10 (Spinachia), 6:10 (Tautogolabrus), and a cumulative 17:27 partition/query recovery path (Menidia), all on octopus07 with eight allocated CPUs and 64 GiB requested memory. This cluster's accounting plugin returned empty `MaxRSS` and `TotalCPU`; the manifest records that absence explicitly rather than inventing a utilization value.
