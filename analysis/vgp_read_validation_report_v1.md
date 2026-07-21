# VGP pilot raw-read validation v1

Date completed: 2026-07-21 UTC  
WG task: `validate-vgp-pilot-reads`  
Canonical shared VGP root: `/moosefs/erikg/vgp`

## Decision

P07 is validated as an exact-individual paired-method sensitivity result and classified `concrete_haplotype_reconstruction_failure`. The downstream action is to preserve the core artifact and provenance, but do not use its pi or PSMC as validated quantitative evidence. This preserves the original artifact and provenance even when validation changes its evidentiary use. P04 remains a valid completed core result with raw validation pending. P09 supplies only a low-coverage compatibility control and no absent estimate is converted to zero.

| Pair | Raw evidence outcome | Assembly pi | Read/pi comparison | PSMC effect | Decision |
|---|---|---:|---:|---|---|
| P04 | Exact CLR run still pending | 0.004604184796 | Not estimable | Not estimable | Retain core; validation pending |
| P07 | Exact BioSample Illumina + HiFi complete | 0.002147219831 | primary common-mask ratio 0.436332 | theta ratio 0.940281; lambda r -0.144006 | concrete_haplotype_reconstruction_failure |
| P09 | One of seven HiFi cells; 0.810x diploid-equivalent | Incomplete at freeze | Not estimable | Not estimable | Compatibility only |

## P07 common-mask pi concordance

The inherited final denominator was independently reconstructed as 267,379,237 bp after the exact non-SNP flank subtraction; 574,122 SNPs reproduce pi = 0.0021472198306856562. Every comparison below restricts both callsets to the same assembly coordinates and the same read-depth mask. Thus differences are paired sensitivity estimates, not independent replication.
The primary depth mask retains 95.677% of the inherited callable bases but only 38.872% of its assembly differences. The excluded 11,557,905 bp contain 350,950 assembly differences, an assembly-difference density 34.8067 times that inside the primary common mask. Accordingly, primary common-mask assembly/read pi are 0.406281/0.177273 times the inherited core pi. This is a strong mappability/collapse sensitivity signal, not an assertion that every excluded difference is false.

| Mask (inclusive DP) | Common bp | Assembly pi | Read pi | Read/assembly | Shared | Assembly-only | Read-only | Jaccard |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| dp5_100 (5-100) | 257,515,131 | 0.001000947 | 0.00038135235 | 0.380992 | 90,676 | 167,083 | 7,528 | 0.341803 |
| dp10_60 (10-60) | 252,149,011 | 0.00087207163 | 0.00038017996 | 0.43595 | 89,291 | 130,601 | 6,571 | 0.394285 |
| dp10_80 (10-80) | 255,821,332 | 0.00087237447 | 0.00038064457 | 0.436332 | 90,340 | 132,832 | 7,037 | 0.392426 |
| dp15_80 (15-80) | 253,529,700 | 0.00078761581 | 0.00037900096 | 0.4812 | 89,418 | 110,266 | 6,670 | 0.433323 |
| dp20_80 (20-80) | 246,528,943 | 0.00072380548 | 0.0003757936 | 0.519191 | 86,484 | 91,955 | 6,160 | 0.468497 |
| dp10_100 (10-100) | 255,940,536 | 0.00087396472 | 0.00038122918 | 0.436207 | 90,396 | 133,287 | 7,176 | 0.391564 |

Changing only the upper cutoff from DP80 to DP100 changes the read/assembly ratio from 0.436332 to 0.436207, whereas raising the lower cutoff from DP5 to DP20 changes it from 0.380992 to 0.519191. Low-depth/mappability exclusion, rather than high-depth collapse filtering, is therefore the dominant callable-mask sensitivity in this read set.
For the predeclared primary DP10-80 mask, the lower concordant bracket is 0.00035313709 and the union upper bracket is 0.00089988195. Assembly-only calls give a candidate assembly false-positive upper bound of 59.520%; read-only calls give a candidate assembly false-negative upper bound of 7.227%. Neither is a proven error rate: the former includes read-caller misses and the latter includes assembly-caller misses. The read heterozygote allele-balance median is 0.482759; 87.469% lie from 0.30 through 0.70.

## Sequence QV, k-mer heterozygosity, and structural diagnostics

Illumina primary mapping is 97.540% with 38.808x whole-H1 mean depth and 99.966% breadth. HiFi primary mapping is 99.910% with 52.115x mean depth and 99.996% breadth.
The Illumina 21-mer containment estimate gives H1 QV 35.9217 (binomial-only 95% interval 35.91594-35.92746), from 407,531,427 assembly k-mer occurrences and 2,183,232 below the trusted read threshold. Correlated k-mers, read errors, coverage bias, and assembly/read dependence remain systematic and are not hidden inside that narrow interval.
The transparent four-component negative-binomial spectrum model estimates k-mer heterozygosity 0.0026652671, 1.24126 times the inherited core pi, with heterozygous/homozygous peaks at 11.34x/22.68x and fit R-squared 0.997875. This is explicitly a model-based spectrum estimate, not a substituted published number or a population estimate; its disagreement with stringent mapped-read calling is retained rather than averaged away.
Across the inherited assembly callable input, the positive-depth mode is 34x; zero depth is 2.591%, below-half-mode is 5.820%, above 1.5x is 6.417%, and above 2x is 0.352%. Low-depth excess diagnoses duplication/mappability/dropout and high-depth excess diagnoses collapse/repeats; neither alone proves an assembly error.

## Direct read support and error bounds

At the DP10-80 assembly SNPs, Illumina classifies 89,321 as supported heterozygotes and 111,859 as depth-qualified homozygous-reference contradictions. The concrete false-positive lower-bound fraction is 50.122% (binomial Wilson 95% interval 49.915%-50.330%); counting every ambiguous, out-of-mask, or unobserved site produces the deliberately conservative candidate upper bound 59.977%.
Independent HiFi pileups classify 91,598 supported and 118,652 contradicted sites, with concrete lower/upper-candidate fractions 53.166%/58.956% and a lower-bound Wilson 95% interval 52.959%-53.373%. HiFi and Illumina share the individual and reference but differ in chemistry and mapping behavior; their agreement is informative without being fully independent.
The primary read callset has 129 homozygous-alt SNP discrepancies, corresponding to mapping-consensus QV 62.97347 with scope limited to accessible homozygous-alt SNPs. Structural and inaccessible errors are excluded from that QV.

## PSMC sensitivity

| Mask | Theta read/assembly | Lambda Pearson r | log-lambda RMSE | time-grid r |
|---|---:|---:|---:|---:|
| dp5_100 | 0.940646 | -0.144007 | 3.88586 | 0.998103 |
| dp10_60 | 0.923544 | -0.143923 | 3.80672 | 0.998223 |
| dp10_80 | 0.940281 | -0.144006 | 3.78808 | 0.998148 |
| dp15_80 | 0.936788 | -0.143914 | 3.69645 | 0.998127 |
| dp20_80 | 0.921357 | -0.143896 | 3.60665 | 0.998044 |
| dp10_100 | 0.942468 | -0.143992 | 3.80473 | 0.998138 |

Only the final completed 64-interval optimization round is compared. Both trajectories use the same H1 coordinates, PSMC parameterization, inherited structural mask, and overlapping biological response; the correlations therefore measure callable/caller sensitivity and must not be treated as an independent demographic replication. Scaling is intentionally not introduced here.

## P09 low-coverage control

The exact P09 cell mapped 495,129/495,557 reads (99.914%); query-union bases were 99.472% of metadata bases, weighted alignment identity was 97.146%, and H1 breadth was 76.309%. Coverage is only 0.810x diploid-equivalent from one of seven cells. QV, k-mer heterozygosity, callable read pi, and PSMC validation are therefore not estimable. These metrics establish compatibility, not callability.

## Raw evidence, failures, reproducibility, and scope

The cumulative acquisition ledger contains 4 verified/reused objects totaling 31,058,137,613 bytes. P07 contributes one HiFi and paired Illumina objects; P09 contributes the retained cell. Every canonical CAS object was rehashed offline. P04's 42,344,746,693-byte exact CLR run remains planned.
A full-size first P07 R2 transfer failed ENA MD5 and gzip integrity and remains quarantined at `/moosefs/erikg/vgp/quarantine/vgp-validation-reads-v1/2026-07-20T183020Z/P07_SRR30200290_SRR30200290_2.fastq.gz.partial.content-mismatch-ena_md5_match` (SHA-256 `defaed9e929d8acf9d58006a3be51c26dd7c1937079f7d47e9d818420475c965`). A clean independent retry produced the verified canonical SHA-256 `c542f6efd9fc1d8f557c89629743fbe4a39584f24002c84185e479057cb443ac`; the failed payload was never promoted.
GNU Guix is frozen at channel commit `44bbfc24e4bcc48d0e3343cd3d83452721af8c36` with profile `/gnu/store/n3vizxfw5ilggrinaj2mmbmng5ja4d6d-profile` and closure digest `ac0cb3601e56ef62b9ef99419de3659b2a2ba59b2aead29bc5f1928b50c83da2`. Slurm jobs, failed/cancelled precursors, requested resources, elapsed times, and MaxRSS availability are retained in `analysis/vgp_read_validation_sacct_v1.tsv`; executed workers, tool digests, stdout/stderr snapshots, per-stage GNU time telemetry, input manifests, and output manifests remain below `/moosefs/erikg/vgp/derived/read-validation/runs/`.

The result is per individual, not a population mean. Binomial intervals exclude correlated-k-mer and mapping systematics. Read and assembly methods intentionally share the individual, H1 coordinate system, structural mask, and in part the molecular data, so covariance is preserved in every machine comparison. The conservative machine decision rule is: severe sequence error if k-mer QV <20; concrete reconstruction failure only if a majority of primary assembly SNPs are depth-qualified homozygous-reference contradictions in both Illumina and HiFi; otherwise pi ratio outside 0.8-1.25 or PSMC lambda r <0.5 is material method discordance, not a deletion trigger. Validation never deletes the original artifact; only a concrete identity, reconstruction, or severe sequence-error failure changes its downstream quantitative use.
