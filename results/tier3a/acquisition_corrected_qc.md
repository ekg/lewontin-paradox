# Tier 3A origin/main SweepGA correction QC

All three tuples passed the atomic publication gate.

Pinned unmodified origin/main commit: `018e4ce49d2c125820e0ac50dc5feaa02d423683`.
Pinned byte-reproducible GNU Guix binary SHA-256: `fa7f0edb9b7e275c288db254046020e136d4267dd5ee043379227ef80da0573b`.
Production used native `--num-mappings 1:1`; no mapping was replaced by post-hoc filtering.
Higher-cap sensitivities were not run in this correction.

## Tuple gates

### Menidia menidia (`menidia_menidia_fMenMen1`)

Mapping records: 3,353; query/target coverage: 0.9305/0.9512; observed threshold multiplicity: 1:1.
IMPG biological minimum numerator/denominator across reported annotation classes: 412/17133.
Representative allele: H1 CM109610.1:6022075 T to H2 CM108917.1:5614205 G; direct REF/ALT checks passed.
SweepGA mapping ran in Slurm array 1761751_0 for 00:41:41 (the allocation later exited `FAILED` at the superseded initial audit). Corrected downstream gate 1761781_0 completed on octopus07 in 00:18:02 with 8 CPUs and 64G requested; accounting reported TotalCPU `00:00:00` and MaxRSS `unavailable`. Full initial/retry sacct is recorded at `results/tier3a/logs/origin-remap-1761751-1761781.sacct`.

### Spinachia spinachia (`spinachia_spinachia_SK-2024b`)

Mapping records: 753; query/target coverage: 0.9827/0.9870; observed threshold multiplicity: 1:1.
IMPG biological minimum numerator/denominator across reported annotation classes: 10/15067.
Representative allele: H1 CM106587.1:2892985 C to H2 CM106669.1:2760454 T; direct REF/ALT checks passed.
SweepGA mapping ran in Slurm array 1761751_1 for 00:16:55 (the allocation later exited `FAILED` at the superseded initial audit). Corrected downstream gate 1761778_1 completed on octopus08 in 00:04:40 with 8 CPUs and 64G requested; accounting reported TotalCPU `00:00:00` and MaxRSS `unavailable`. Full initial/retry sacct is recorded at `results/tier3a/logs/origin-remap-1761751-1761781.sacct`.

### Tautogolabrus adspersus (`tautogolabrus_adspersus_fTauAds1`)

Mapping records: 1,972; query/target coverage: 0.9491/0.9885; observed threshold multiplicity: 1:1.
IMPG biological minimum numerator/denominator across reported annotation classes: 27/15337.
Representative allele: H1 CM036534.1:34981587 C to H2 JAJGRG010000537.1:207731 A; direct REF/ALT checks passed.
SweepGA mapping ran in Slurm array 1761751_2 for 00:16:45 (the allocation later exited `FAILED` at the superseded initial audit). Corrected downstream gate 1761778_2 completed on octopus07 in 00:06:26 with 8 CPUs and 64G requested; accounting reported TotalCPU `00:00:00` and MaxRSS `unavailable`. Full initial/retry sacct is recorded at `results/tier3a/logs/origin-remap-1761751-1761781.sacct`.
