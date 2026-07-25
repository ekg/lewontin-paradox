# Three-pair VGP bounded-range reliability pilot

## Outcome

Three actual same-individual H1/H2 biological results completed. Every result was derived from bounded H1-coordinate IMPG queries; no all-genome graph, exhaustive all-partition query, global IMPG lace, or global partition-assignment ledger was created. The optional aggregate VCF/BCF is a `bcftools concat` of already normalized, nonoverlapping range products.

| Pair | Failure class | Heterozygous SNPs | Callable bp | Pi | Ranges | PSMC | Annotation |
|---|---|---:|---:|---:|---:|---|---|
| P07 (Spinachia spinachia) | previously_successful_exact_annotation_pair | 316,631 | 270,531,638 | 0.00117040284952 | 109 | 200/200 finite; centered=true | exact_native |
| P03 (Colius striatus) | prior_fastga_execution_failure | 1,632,584 | 875,683,638 | 0.00186435366513 | 407 | 200/200 finite; centered=true | missing_nonblocking |
| P02 (Pseudorca crassidens) | prior_impg_graph_fasta_identifier_failure | 4,195,014 | 1,707,746,195 | 0.00245646221452 | 930 | 200/200 finite; centered=true | missing_nonblocking |

## Architecture correction and canary

The prohibited global chains were canceled: 1797004, 1797005, 1797006, 1797007, 1797008, 1797029, 1797030, 1797031, 1797032, 1797033. Their cancellation is technical provenance and is not a species exclusion.

The fresh P07 canary queried `CM106587.1:0-5000000` (5,000,000 bp) using 2,500 complete native partitions. Its 3,871 normalized IMPG keys exactly matched the like-for-like clean P07 subset (SHA-256 `b7f64cd018a693ade59de7daec49b66d9f57e6ea7baa888c6e9e1c7099da1405`). Callable accounting was 3,116,900 bp, peak local graph state was 4,187,923 bytes, and local graph temporaries were deleted.

The frozen range rule is complete consecutive native H1 partition boundaries, targeting 5 Mb with a hard 20 Mb ceiling. Unaligned H1 ranges remain explicit non-query ranges and enter consensus/PSMC as non-callable sequence; they are never silently omitted.

P02 job 1799308 exposed a graph-allele REF/ALT orientation error in the first range. The corrected job 1804024 completed all 925 query ranges, then correctly rejected the H1-monomorphic `ALT=.` record at CM078529.1:2167395 during reduction. Its zero-copy-failure bundle retained all 7,412 result files plus the exact plan and graph index. Resume job 1804979 removed only post-reconstruction monomorphic records, strict-H1-revalidated all 925 VCF/BCF pairs, issued zero IMPG queries, and completed the biological result. These failures and the resume are technical provenance, not an exclusion.

## Independent validation

### P07

Graph IDs unresolved=0; silently omitted=0. Range duplicate keys=0; unowned callable bp=0; multiply owned callable bp=0.

| Stratum | Range | Variants | Callable bp | Consensus non-N bp | PSMC N/K/T |
|---|---|---:|---:|---:|---|
| early | r000000 (CM106587.1:0-5000000) | 3,871 | 3,096,868 | 3,096,868 | 21,652/890/27,458 |
| middle | r000053 (CM106596.1:15000000-18717046) | 7,769 | 2,797,401 | 2,797,401 | 11,232/1,545/24,394 |
| late | r000107 (JBLUUC010000036.1:0-50000) | 3 | 18,526 | 18,526 | 323/0/177 |

### P03

Graph IDs unresolved=0; silently omitted=0. Range duplicate keys=0; unowned callable bp=0; multiply owned callable bp=0.

| Stratum | Range | Variants | Callable bp | Consensus non-N bp | PSMC N/K/T |
|---|---|---:|---:|---:|---|
| early | r000000 (CM054345.1:0-4999475) | 8,658 | 4,009,965 | 4,009,965 | 11,562/2,093/36,340 |
| middle | r000165 (CM054355.1:9999990-14999595) | 9,094 | 4,956,871 | 4,956,871 | 1,952/5,048/42,997 |
| late | r000406 (JARBXP010000074.1:0-78732) | 21 | 39,726 | 39,726 | 408/0/380 |

### P02

Graph IDs unresolved=0; silently omitted=0. Range duplicate keys=0; unowned callable bp=0; multiply owned callable bp=0.

| Stratum | Range | Variants | Callable bp | Consensus non-N bp | PSMC N/K/T |
|---|---|---:|---:|---:|---|
| early | r000000 (CM078529.1:0-4999147) | 34,527 | 3,019,445 | 3,019,445 | 23,112/4,805/22,075 |
| middle | r000462 (CM078547.1:39992839-44991423) | 5,075 | 3,600,652 | 3,600,652 | 15,764/2,322/31,901 |
| late | r000929 (CM078552.1:0-16382) | 102 | 8,777 | 8,777 | 80/60/24 |

## Exact native annotation partitions

P07’s GFF is exact-native and dictionary-identical to H1. Annotation ranges were queried/intersected directly in H1 coordinates; P02/P03 lack cataloged exact-native annotations, so only their annotation products are absent.

| Partition | Heterozygous variants | Callable bp | Estimate | Estimator |
|---|---:|---:|---:|---|
| CDS | 7,603 | 26,152,171 | 0.000290721561893 | normalized heterozygous allele records / callable CDS bp |
| fourfold | 2,212 | 4,225,276 | 0.000523516096937 | heterozygous SNPs / callable fourfold bp |
| fourfold_W | 987 | 1,452,971 | 0.00067929779741 | heterozygous SNPs / callable fourfold AT bp |
| fourfold_S | 1,225 | 2,772,305 | 0.000441870573404 | heterozygous SNPs / callable fourfold GC bp |
| WS | 868 | 1,452,971 | 0.00059739664453 | AT-to-GC SNPs / callable fourfold AT bp |
| SW | 920 | 2,772,305 | 0.000331853818393 | GC-to-AT SNPs / callable fourfold GC bp |

## Controlled backend comparison

P03’s byte-identical staged-FASTA control has 17,317,077 bp of common target coverage, target-coverage Jaccard 0.970695, and exact raw-variant Jaccard 0.252922. Both reconstructions have zero coordinate/REF-ALT failures. Backend-specific gap placement is retained as a limitation rather than a biological exclusion.

## Resource telemetry

| Job | Name | State | Elapsed | CPUs | MaxRSS | Node |
|---|---|---|---|---:|---:|---|
| 1797004 | vgp10-P02-impg | CANCELLED by 1001 | 02:39:32 | 32 |  | octopus02 |
| 1797005 | vgp10-P02-variants | CANCELLED by 1001 | 00:00:00 | 0 |  | None assigned |
| 1797006 | vgp10-P02-consensus | CANCELLED by 1001 | 00:00:00 | 0 |  | None assigned |
| 1797007 | vgp10-P02-psmc | CANCELLED by 1001 | 00:00:00 | 0 |  | None assigned |
| 1797008 | vgp10-P02-finalize | CANCELLED by 1001 | 00:00:00 | 0 |  | None assigned |
| 1797029 | vgp10-P03-impg | CANCELLED by 1001 | 02:36:49 | 16 |  | octopus02 |
| 1797030 | vgp10-P03-variants | CANCELLED by 1001 | 00:00:00 | 0 |  | None assigned |
| 1797031 | vgp10-P03-consensus | CANCELLED by 1001 | 00:00:00 | 0 |  | None assigned |
| 1797032 | vgp10-P03-psmc | CANCELLED by 1001 | 00:00:00 | 0 |  | None assigned |
| 1797033 | vgp10-P03-finalize | CANCELLED by 1001 | 00:00:00 | 0 |  | None assigned |
| 1797698 | vgp-P07-range-canary | FAILED | 00:00:09 | 8 |  | octopus02 |
| 1797712 | vgp-P07-range-canary | FAILED | 00:00:08 | 8 |  | octopus02 |
| 1797720 | vgp-P07-range-canary | FAILED | 00:03:55 | 8 |  | octopus02 |
| 1797758 | vgp-P07-range-canary | FAILED | 00:03:49 | 8 |  | octopus02 |
| 1797776 | vgp-P07-range-canary | CANCELLED by 1001 | 00:00:37 | 8 |  | octopus02 |
| 1797782 | vgp-P07-range-canary | COMPLETED | 00:03:13 | 8 |  | octopus02 |
| 1797792 | vgp-bounded-P07 | CANCELLED by 1001 | 00:09:09 | 16 |  | octopus02 |
| 1797793 | vgp-bounded-P03 | CANCELLED by 1001 | 00:00:00 | 0 |  | None assigned |
| 1797794 | vgp-bounded-P02 | CANCELLED by 1001 | 00:00:00 | 0 |  | None assigned |
| 1797847 | vgp-bounded-P07 | CANCELLED by 1001 | 00:02:40 | 16 |  | octopus02 |
| 1797848 | vgp-bounded-P03 | CANCELLED by 1001 | 00:00:00 | 0 |  | None assigned |
| 1797849 | vgp-bounded-P02 | CANCELLED by 1001 | 00:00:00 | 0 |  | None assigned |
| 1797867 | vgp-bounded-P07 | FAILED | 04:29:38 | 16 |  | octopus02 |
| 1797868 | vgp-bounded-P03 | CANCELLED by 1001 | 00:00:00 | 0 |  | None assigned |
| 1797869 | vgp-bounded-P02 | CANCELLED by 1001 | 00:00:00 | 0 |  | None assigned |
| 1799306 | vgp-bounded-P07 | COMPLETED | 04:31:36 | 16 |  | octopus02 |
| 1799307 | vgp-bounded-P03 | COMPLETED | 14:20:08 | 16 |  | octopus02 |
| 1799308 | vgp-bounded-P02 | FAILED | 00:18:41 | 32 |  | octopus02 |
| 1804024 | vgp-bounded-P02 | FAILED | 05:25:55 | 32 |  | octopus02 |
| 1804979 | vgp-bounded-P02-resume | COMPLETED | 17:34:02 | 32 |  | octopus02 |

## Remaining limitations

- FastGA remains unreliable for P03 at whole-assembly scale; the pinned WFMASH fallback is infrastructure provenance, not a biological exclusion.
- SweepGA/FastGA and WFMASH differ in raw gap placement on the controlled overlap despite high shared target coverage.
- Pinned IMPG lace with one thread did not progress; each bounded range uses the tested two-thread minimum.
- P02 and P03 lack cataloged exact native annotations; their core range, diversity, consensus, and PSMC outputs continue.
- Assembly-derived same-individual haplotype diversity is not a substitute for population sampling.
- Cluster Slurm accounting did not populate MaxRSS for these jobs; elapsed time, allocated CPUs, state, and node are retained, and direct monitoring observed the serial P02 consensus process at 13,777,560 KiB RSS.

All limitations above are technical or interpretive limitations, not biological exclusions.

## Evidence identities

- Selection freeze SHA-256: `c98ae68e5f1f7435f0bb816c1c04718d1c378a2fbbdda38aab0dbd67912a0400`
- Execution audit SHA-256: `f57b8dc2fcfd90eb63a5b1c03fb5e2193f18b87b7b008d4d29749f88ea7b8c15`
- Bounded transition record SHA-256: `7452ab49ba00ef6fc98403c465c0b809e274e0ec382c18a730596bb532527e50`
- Scheduler telemetry SHA-256: `21a99976dbed34f7e6cbfebeb923fb22b7ac42a94ed77f47e7d28ff99ed7cb24`
