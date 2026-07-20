# Independent review of the real VGP ten-pair pilot

Review cutoff: 2026-07-20 16:06:19 UTC

Authorization: `vgp10-auth-20260718-v2`

Canonical shared VGP root: `/moosefs/erikg/vgp`

WG task: `review-vgp-real-pilot`

## Decision and numbers first

**CONDITIONAL_GO for a bounded scale-out of the core H1/H2 mapping, exact
variant reconstruction, masks, consensus, and callable-diversity estimator.**
Do not scale the current PSMC bootstrap array until its sampling population is
repaired and revalidated. Annotation and raw-read validation remain separate
confidence decisions and do not veto the two valid core results.

| Quantity | Independently measured result |
|---|---:|
| Completed pairs | 2/10: P04 and P07 |
| Final callable sequence | 1,038,160,202 bp total; 267,379,237–770,780,965 per pair |
| Final heterozygous SNPs | 4,122,940 total; 574,122–3,548,818 per pair |
| Callable phased-haplotype π | 0.0021472198306856562–0.004604184795871289 |
| Finite PSMC primary trajectories | 2/2, each with 64 intervals |
| Numerically finite bootstrap outputs | 400/400, each with 64 intervals |
| Bootstrap uncertainty sets admitted | **0/2**: both resampling distributions are shifted from the primary input population |
| Sensitivity trajectories | 9 scenarios and 576 rows per pair; 1,152 rows total |
| Frozen terminal primary failures | 4: P01, P02, P03, P05 |
| No biological result at cutoff | 4: P06 failed its corrected mapping; P08/P09/P10 were still mapping |

The core statistic is heterozygous SNPs after the ordered callable mask and
the predeclared ±10-bp indel-flank mask, divided by final callable H1 bases.
It is phased-haplotype divergence for one diploid individual, not a population
mean.

## Pair-level recomputation

The reviewer did not import the production implementation. It streamed each
canonical VCF twice, independently summed and indexed the callable BED,
rebuilt the merged callable portion of every ±10-bp indel window, selected the
final SNP numerator, counted the consensus and PSMCFA symbols, and parsed the
final native `RD`/`TR`/`RS` iteration of all 402 PSMC files.

| Measure | P04, *Falco naumanni* | P07, *Spinachia spinachia* |
|---|---:|---:|
| H1 universe | 1,215,719,661 | 407,561,107 |
| Pre-indel callable bp | 781,674,063 | 272,818,693 |
| Pre-indel callable fraction | 0.6429722970 | 0.6693933457 |
| Indel-mask loss | 10,893,098 | 5,439,456 |
| Final callable bp | 770,780,965 | 267,379,237 |
| Final callable fraction | 0.6340120915 | 0.6560469888 |
| Normalized variants | 4,794,096 | 1,488,885 |
| Callable variants before indel flanks | 4,431,550 | 1,316,146 |
| Callable SNPs before indel flanks | 3,804,357 | 944,917 |
| Callable indels | 627,193 | 371,229 |
| Callable SNPs masked by indel flanks | 255,539 | 370,795 |
| Final heterozygous SNPs | 3,548,818 | 574,122 |
| π | 0.004604184795871289 | 0.0021472198306856562 |
| Mask accounting discrepancy | 0 bp | 0 bp |
| Consensus callable symbols | 770,780,965 | 267,379,237 |
| Consensus heterozygous IUPAC symbols | 3,548,818 | 574,122 |
| PSMCFA `K / T / N` | 2,137,417 / 4,907,207 / 5,112,715 | 110,846 / 2,362,279 / 1,602,504 |

The independently counted consensus IUPAC bases equal the final SNP numerator
for both pairs, and callable consensus symbols equal the final denominators.
This closes the variant → mask → consensus inputs without relying on a stage
JSON summary.

## PSMC behavior and bootstrap finding

Both unscaled primary files are finite and internally parseable. P04 has
`theta_0 = 0.524450` per 100-bp bin, `time_2N0` 0–8.056210, and lambda
0.492194–2.874536. Its maximum is early (interval 4), its minimum is interval
40, and its oldest lambda is 0.508550. P07 has `theta_0 = 0.032921`,
`time_2N0` 0–29.843129, and lambda 0.649260–159.317994. Its maximum is in the
deep tail (interval 52), and its oldest lambda remains 51.975015. P07's extreme
old tail is a biological/model outlier to carry as a sensitivity flag; it is
not a reason to erase the finite primary.

Numerical bootstrap completion is real but statistical bootstrap validity is
not:

| Diagnostic | P04 | P07 |
|---|---:|---:|
| Finite files with 64 intervals | 200/200 | 200/200 |
| Primary `theta_0` | 0.524450 | 0.032921 |
| Bootstrap `theta_0` range | 1.866327–2.525863 | 0.052696–0.097470 |
| Bootstrap `theta_0` median | 2.1799945 | 0.075282 |
| Primary / bootstrap median | 0.240574 | 0.437302 |
| Primary inside observed bootstrap range | no | no |

This is not a retrospective effect-size cutoff. A bootstrap must resample the
same input population as its target. The primary PSMCFA retains masked bins,
whereas `bootstrap_units.5mb.bed` is frozen from callable intervals and the
bootstrap encoder samples those intervals as records. Thus baseline masked
bins outside callable intervals cannot be sampled. A reconstructed frozen P04
replicate 1 contained `K=2,180,946`, `T=4,981,528`, and `N=984,769` (12.1%
masked) versus primary `K=2,137,417`, `T=4,907,207`, and `N=5,112,715` (42.1%
masked). Running PSMC with `-b` on that shifted input does not produce an
uncertainty distribution centered on the primary estimand. Finiteness passes
the operational 190/200 gate but cannot cure the sampling mismatch.

Consequently, both primaries may be retained as descriptive unscaled outputs;
the 400 bootstrap curves, bootstrap confidence bands, and any uncertainty
claim based on them are not admitted. The nine mutation-rate × generation-time
grids remain labeled generic sensitivity transformations, never species
calibrations.

## Exact closed-world disposition

| Pair | Review disposition at cutoff |
|---|---|
| P01 | **Hard current-primary invalidity.** IMPG job 1782614, exit `2:0`: frozen graph sequence `CM070280.1` was absent from both staged FASTAs near the end of the 1,150,008-region census. No shard omission/substitution is allowed. |
| P02 | **Hard current-primary invalidity.** IMPG job 1782558, exit `2:0`: `JBDMAU010000296.1`, `JBDMAU010000433.1`, and `JBDMAU010000398.1` were absent from both staged FASTAs after most of the 1,375,084-region census. |
| P03 | **Hard current-primary invalidity.** Mapping job 1787016 failed after 37:25 even though cwd, all temp variables, indexes, and 772 descriptors were proved inside private scratch with 598,585,008,128 bytes free at launch. |
| P04 | **Valid core completion.** Exact π and consensus inputs above. Annotation was not applicable under the frozen exact-binding rules. |
| P05 | **Hard current-primary invalidity.** IMPG job 1782284, exit `2:0`: `JAYKKY010000001.1` and `JAYKKY010000134.1` were absent, and native extraction rejected mixed variant types. |
| P06 | **No result; not a confidence covariate.** Post-upstream corrected mapping job 1787022 failed after 38:41, exit `1:0`, on octopus09. Its stdout/stderr are empty; this review does not project a prior failure cause onto the new job. |
| P07 | **Valid core completion.** Exact π and consensus inputs above; exact-native annotation partitions separately pass. Extreme PSMC deep-tail behavior is an outlier flag, not core invalidity. |
| P08 | **No result at cutoff.** Corrected mapping 1787028 was running for 35:19. |
| P09 | **No result at cutoff.** Corrected mapping 1787065 was running for 15:31. |
| P10 | **No result at cutoff.** Corrected mapping 1787071 was running for 7:46. |

P06/P08/P09/P10 are not zeros, biological outliers, or low-confidence
estimates: no biological estimate existed at the review cutoff.

## Corrected resource model from `sacct`

`sacct` reports elapsed allocation capacity and requested memory on this
cluster. It reports no `MaxRSS` or usable `TotalCPU` rows, so actual peak RSS,
CPU consumption, I/O, and scratch high-water are not imputed. GiB-hours below
are requested-memory × elapsed allocation hours, not measured consumption.

### Outcome-producing P04 template

| Stage | Allocations | Elapsed per allocation | Allocated core-hours | Requested GiB-hours |
|---|---:|---:|---:|---:|
| Mapping (32 CPU, 160 GiB) | 1 | 895 s | 7.956 | 39.778 |
| IMPG (48 CPU, 160 GiB) | 1 | 23,349 s | 311.320 | 1,037.733 |
| Variants (24 CPU, 160 GiB) | 1 | 175 s | 1.167 | 7.778 |
| Consensus (24 CPU, 160 GiB) | 1 | 1,820 s | 12.133 | 80.889 |
| PSMC primary + bootstraps (4 CPU, 16 GiB) | 201 | median 1,799 s; range 1,778–2,518 | 402.901 | 1,611.604 |
| Finalize (4 CPU, 16 GiB) | 1 | 3 s | 0.003 | 0.013 |
| **Total allocation capacity** | **206** | **108.015 summed allocation-hours** | **735.480** | **2,777.796** |

At the authorized PSMC throttle of 21, 201 tasks with the observed median imply
approximately 5.0 hours of array makespan. That expenditure must be held until
bootstrap validity is repaired. For P04-like 2.372-Gbp pairs, retain the
observed mapping allocation and the 48-CPU IMPG allocation plus the measured
600-GiB free-scratch gate. Scratch is a separately observed contract because
`sacct` cannot measure it.

P07's final checkpoint-reusing core allocation was 2.459 hours at 32 CPU and
128 GiB: 78.693 allocated core-hours and 314.773 requested GiB-hours. Its two
successful IMPG rescue allocations add 8.489 core-hours and 33.956 requested
GiB-hours. Because the final canary job reused checkpoints, it must not be
misrepresented as a fresh end-to-end small-pair stage model.

Across all attempts in the frozen telemetry tables, the nine-pair pilot used
3,235.853 allocated core-hours and 16,312.480 requested GiB-hours; the canary
history used 368.616 core-hours and 1,480.142 requested GiB-hours. Those costs
include failed/cancelled work and are reported as operational history, not as
the success template.

The corrected scale-out model is therefore stratified, not a two-point linear
regression:

1. Use the measured P04 stage allocations only for the represented medium
   size stratum and P07 only as a checkpoint-reusing small-pair observation.
2. Preflight graph/FASTA dictionaries and mixed-site behavior before expensive
   IMPG census work; never omit a failed region.
3. Retain private node-local scratch containment and measured-free gates.
4. Run one bounded wave and recalibrate after each newly completed size
   stratum; two completions cannot identify a reliable genome-size slope.
5. Budget 201 PSMC allocations per admitted pair, but do not submit the 200
   bootstraps until the primary-population mismatch is repaired.

## Transparent decision boundary

The frozen core criteria were exact pair/digest provenance, exact bidirectional
1:1 mapping and reconstruction, at least 100,000,000 callable bp, at least 60%
pre-indel callability, zero mask-accounting discrepancy, a heterozygosity-
preserving consensus, a finite unscaled primary, and at least 190/200 finite
bootstrap executions. P04 and P07 pass those declared content/execution gates.

Independent review adds no favorable post-hoc biological threshold. It does
enforce the definition of a bootstrap—resampling the target population—and
therefore separates a numerical completion from valid uncertainty. This makes
the scale decision conditional rather than erasing the valid core estimates.

The scale-out conditions are:

- preserve every frozen core gate and canonical-root-only path;
- repair PSMC bootstraps by sampling blocks from the primary PSMCFA population,
  preserving callable and masked bins, then require centering diagnostics;
- begin only with resource strata represented by a completed sentinel;
- preserve fail-closed dictionary and scratch preflights;
- keep annotation and raw-read decisions as separate confidence columns.

## Evidence handoff: measured versus unvalidated

Measured directly in this review:

- both canonical callable BEDs and normalized VCFs;
- both consensus FASTAs and PSMCFA symbol populations;
- both primary PSMC files and every one of 400 bootstrap PSMC files;
- all sensitivity scenario rows and their generic-source label;
- frozen `sacct` allocation states/capacity plus the 16:06:19 UTC live delta;
- the exact four frozen terminal failures and P06's later scheduler failure.

Still unvalidated:

- raw-read concordance and chemistry for P04/P07;
- assembly QV, BUSCO completeness/duplication, collapse/copy-number, k-mer,
  and standalone repeat evidence;
- PSMC bootstrap uncertainty until resampling is repaired;
- species-specific mutation rates and generation times;
- actual RSS, CPU consumption, I/O, and scratch high-water from `sacct`;
- any population-level diversity/demographic inference; and
- any biological output for the eight non-completed pairs.

Confidence decisions are therefore: core **PASS** for P04/P07; annotation
**PASS** for exact-native P07 and **NOT_APPLICABLE** for P04; raw reads
**UNVALIDATED** for both; assembly validation covariates **UNVALIDATED / tier
C**; PSMC primaries **DESCRIPTIVE ONLY**; bootstrap uncertainty
**STATISTICAL_FAIL pending repair**.

Machine evidence is in `analysis/vgp_real_pilot_independent_review_v1.json`.
The frozen live scheduler delta is
`analysis/vgp_real_pilot_review_live_sacct_v1.tsv`. Both record the canonical
shared VGP root; no data were written below a Lewontin-paradox-named VGP data
directory.
