# Real VGP biological canary: P07

Date completed: 2026-07-19 UTC  
WG task: `run-vgp-real-canary`  
Authorization: `vgp10-auth-20260718-v2`  
Canonical shared VGP root: `/moosefs/erikg/vgp`

## Result

**PASS.** A real same-individual H1/H2 biological canary for *Spinachia spinachia* (`P07`, individual `fSpiSpi1`, BioSample `SAMN36735485`) completed on Slurm through the pinned GNU Guix environment. Core Slurm job `1781798` completed with exit `0:0` in `02:27:33` on `octopus02`, using 32 CPUs and 128 GiB. Exact-native annotation job `1781559` completed with exit `0:0` in `00:03:33` on `octopus11`, using 8 CPUs and 128 GiB.

The verified result contains whole-assembly SweepGA mapping, exact bidirectional 1:1 multiplicity, IMPG partitioning and 203,698 region-focused VCF extractions, normalized variants, ordered reason-coded masks, a diploid IUPAC consensus, callable diversity, a PSMCFA, an unscaled PSMC trajectory, 200/200 finite boundary-aware block bootstraps, nine explicitly labeled scaling-sensitivity scenarios, and exact-native annotation partitions. An independent implementation reconstructed the mask and a 5-Mbp variant subset, matched the VCF and BCF allele-for-allele, and passed the promoted result.

The canonical promoted core is:

`/moosefs/erikg/vgp/pilot/outputs/vgp10-auth-20260718-v2/P07/core`

The exact annotation result is:

`/moosefs/erikg/vgp/pilot/outputs/vgp10-auth-20260718-v2/P07/annotation/exact_partitions.json`

## Frozen inputs and runtime

| Item | Value |
|---|---|
| H1 | `GCA_048126635.1`; CAS SHA-256 `2b758e606304f7cb5e795d7939979b08c21bf4f3eac7ea3cf1c6ab0a463733c7` |
| H2 | `GCA_048127205.1`; CAS SHA-256 `6e7e9d88b88a3d030d80009191a06c525fbcc50044db05da4d81a8e9ad97ed40` |
| H1 universe | 407,561,107 bp in 38 sequences |
| Guix closure | SHA-256 `8fcdb32021f1cd8eac839509cff47ab6bdd63b656b30e243fdf78d3c4ba24f9d` |
| SweepGA | SHA-256 `fa7f0edb9b7e275c288db254046020e136d4267dd5ee043379227ef80da0573b` |
| IMPG | SHA-256 `c587dc2326cd24f887b1fcb3938404229ad0f0a27ef0773e90c287b1ade160d4` |
| bcftools 1.14 | SHA-256 `79637872d29b03be83293a56d297519ceb74d861ff61888966d5e17157f57dd4` |
| PSMC 0.6.5 | SHA-256 `2c64b0ce7f68a6251e4795f3b88caa7b53769e12a867b6f59918312f66735f57` |

All new canonical biological data were written beneath `/moosefs/erikg/vgp`. The legacy Lewontin-paradox-named VGP root was used only as migration input by the upstream verified hard-link migration; this run did not create canonical outputs there or redownload accepted objects.

## SweepGA mapping and exact multiplicity

SweepGA ran across the complete H1 and H2 assemblies with `--num-mappings 1:1`. Its native output contained 40,763 records and nevertheless reached overlap depth 20 on both query and target axes. This is a tool-behavior finding, not a biological exclusion: `--num-mappings 1:1` did not alone make the native records geometrically disjoint.

A deterministic length/ANI-ranked bidirectional disjoint filter only removed native records; it did not add or alter alignments. The retained PAF contains 4,297 records and independently verifies maximum query depth 1 and maximum target depth 1. It covers 276,045,662 H1 bp. The retained PAF SHA-256 is `0bcf45366bca9217cd012242cedf79921a6cafd6b1fb1f5247db7e54512d00fa`.

Exact CIGAR reconstruction accounted for 274,174,454 match bp, 1,350,745 mismatch bp, 520,463 H1 deletion bp, and 521,353 H2 insertion bp. Every one of the 4,297 PAF records reconstructed its oriented H2 sequence before normalization; reconstruction failures were zero.

## IMPG partitioning and regional extraction

IMPG indexed 813,035,042 sequence bp across the two assemblies. The durable native partition BED contains 409,805 assembly-specific rows. The H1-focused extraction retained 203,698 uniquely named region rows representing 142,198 unique native partition identifiers.

All 203,698 real regional VCFs were observed nonempty before lacing, totaling 351,783,581 bytes. The regional census was recorded before the transient shards were removed. The final genotype-aware rescue projected nonreference alleles carried by 142,286 H2 region samples, laced 16 coordinate groups with IMPG, reconciled exact boundary overlaps with pinned bcftools, sorted the result, and promoted it atomically. The durable laced VCF is 140,335,516 bytes with SHA-256 `60f4fe1d94d3843e45b064a2ecce6aae1d57dbc1a0b2489bcf8fb7aff70464a6`.

The normalized IMPG comparison callset contains 572,461 records. It is retained as a scientifically useful comparison object. The primary consensus callset is instead derived from the exact retained SweepGA CIGARs because the graph representation and normalization boundary behavior of the laced IMPG comparison could not exactly reconstruct every retained PAF segment. This distinction is explicit; the IMPG extraction was completed and preserved, while exact haplotype reconstruction determined the primary representation.

## Variants, masks, consensus, and diversity

The exact PAF converter produced 1,488,893 raw variant blocks, omitted 105 zero-effect CIGAR blocks, and normalized to 1,488,885 alleles:

| Quantity | Count |
|---|---:|
| Normalized alleles | 1,488,885 |
| Normalized SNPs | 954,193 |
| Normalized indels | 534,692 |
| Pre-indel-mask callable alleles | 1,316,146 |
| Pre-indel-mask callable SNPs | 944,917 |
| Pre-indel-mask callable indels | 371,229 |
| Variants excluded by the declared callable mask | 172,739 |

The independent BCF query returned the same 1,488,885 normalized alleles in the same order as the consensus VCF. The primary normalized BCF SHA-256 is `28458586b7a95883ddeae203bce70d790c48fb092472265106b21916a53a8358`; the normalized VCF.gz SHA-256 is `908491859c2eeae724492e6a9eb9b84c60339def39139489ae6c0f67d6037553`.

The declared reason-coded mask reconciles exactly:

| Primary reason | Excluded bp |
|---|---:|
| `not_1to1` | 131,491,345 |
| `repeat_or_low_complexity_primary` | 2,706,506 |
| `h2_gap_or_N` | 520,463 |
| `h1_gap_or_N` | 24,100 |
| All other declared reasons | 0 |
| Pre-indel callable | 272,818,693 |
| H1 universe | 407,561,107 |

The accounting discrepancy is zero and the pre-indel callable fraction is 0.6693933457. The predeclared ±10-bp indel policy masks a further 5,439,456 H1 bp. It takes precedence over SNP IUPAC encoding: 370,795 otherwise-callable SNPs fall inside those indel masks. The final consensus therefore contains 267,379,237 callable bp and 574,122 encoded heterozygous SNPs.

Genome-wide callable diversity is estimated as final callable heterozygous SNPs per final callable H1 base:

`pi = 574,122 / 267,379,237 = 0.0021472198306856562`

For transparency, the pre-indel callable normalized-allele rate is 1,316,146 / 272,818,693 = 0.0048242515405643414; this is not labeled nucleotide diversity because it includes indels and SNPs later excluded by the indel-flank policy.

The diploid consensus FASTA is 414,354,384 bytes with SHA-256 `d6dc1b0c1b500942464f962ec643219a41b885e1b142adc61457d5f6344fc98e`. The PSMCFA is 4,144,150 bytes with SHA-256 `2e68a0077ed3701e069eecf1c11f6449fc2a7317d807b4417c17ffa14287197a`; it contains 2,362,279 callable homozygous bins (`T`), 110,846 heterozygous bins (`K`), and 1,602,504 masked bins (`N`).

## PSMC and boundary-aware bootstraps

The primary PSMC used `-N25 -t15 -r5 -p 4+25*2+4+6` and remains unscaled. The final native PSMC iteration has fitted `theta_0 = 0.032921`, 64 intervals, unscaled time coordinates from 0 to 29.843129 in units of 2N0, and lambda values from 0.649260 to 159.317994. The primary PSMC SHA-256 is `bd1f85d3e84e6ed3065d3329fc5556289341a0e1b8571d0922a302d143e3dd89`.

The workflow ran 200 deterministic, boundary-aware 5-Mbp block bootstraps. Blocks never cross contigs or callable-mask discontinuities. Every bootstrap completed, parsed into 64 unique intervals, had positive fitted theta, and contained only finite time/lambda values: 200/200 finite, exceeding both the task minimum of 100 and the workflow gate of 190/200.

The unscaled trajectory is primary. Absolute outputs are a separate 3 × 3 sensitivity grid:

- mutation rates: 5e-9, 1e-8, and 2e-8 per nucleotide per generation;
- generation times: 1, 2, and 4 years;
- PSMC bin size: 100 bp, explicitly included in every row;
- source label: `predeclared_generic_sensitivity_grid_not_species_calibration`.

The nine scenario identifiers are fully explicit (`SENS_MU…_G…Y`). They produce 576 labeled rows. They are generic sensitivity transforms, not empirical *S. spinachia* calibrations, and must not be interpreted as independently estimated mutation rates or generation times.

## Exact-native annotation partitions

Annotation `GCA_048126635.1-GB_2025_08_04` is bound to the exact H1 accession. Its 38 sequence regions equal the H1 dictionary, its CAS GFF SHA-256 is `8f640543accd8081d1b7048eda32c9f1eef33b02f321b7b0f8adcf3b01dd6838`, and it contains 20,293 canonical transcripts. A total of 542 frame-discordant overlap positions were excluded under the frozen policy.

| Partition | Callable bp | Heterozygous variants | Estimate |
|---|---:|---:|---:|
| CDS | 26,185,407 | 15,849 | 0.0006052608 |
| Fourfold | 4,230,801 | 3,586 | 0.0008475936 |
| Fourfold W | 1,455,109 | 1,697 | 0.0011662357 |
| Fourfold S | 2,775,692 | 1,889 | 0.0006805510 |
| W→S | 1,455,109 | 1,346 | 0.0009250166 |
| S→W | 2,775,692 | 1,318 | 0.0004748365 |

These are exact matched partitions for this one same-individual pair. They are not by themselves a population-level gBGC test. GC3 was not emitted because no separate frozen GC3 overlap policy was authorized.

## Independent reconstruction

The independent auditor does not import the production pipeline helper. It:

1. rebuilt first-reason-wins masks from raw flag BEDs with an independent event sweep and required exact equality to every emitted BED;
2. recounted all normalized and callable variants;
3. queried the BCF independently with pinned bcftools and required exact ordered allele equality to the VCF;
4. reconstructed the first 5,000,000 bp of `CM106587.1`, finding 3,116,900 callable bp and 7,702 callable normalized variants;
5. independently checked query and target PAF depth;
6. verified PSMCFA symbols, 200 finite bootstraps, scenario labels/sources/bin size, exact annotation, stage digests, promotion, and scheduler telemetry.

The independent subset is `analysis/vgp_real_canary_variant_subset_v1.tsv`, SHA-256 `a60e6008e7cd183e644ba3ffa040728cf81903f15559882e7f7c32a3a16d9243`. The final audit result is `PASS`.

## Checkpoints, promotion, resources, and telemetry

The final checkpoint contains 207 completion sentinels: preflight, mapping, IMPG, variants, consensus, 201 PSMC stages, and finalization. Promotion rehashed 468 payload files totaling 1,267,400,116 bytes, copied them to a same-filesystem partial target, rehashed the copy, and atomically renamed it. The canonical promotion manifest SHA-256 is `c7fc4a45213b7113b5b63d8ccb47812e968dee13c5419c54a108d968ed7d9bc1`; the source checkpoint is preserved.

All biological attempts requested 128 GiB. No attempt failed from OOM or time exhaustion, so the authorized 256/512-GiB escalation packets were not scientifically indicated. The final PSMC run was observed at approximately 377 MiB RSS for the primary and approximately 64 MiB per bootstrap process. Exact annotation peaked near 13.6 GiB RSS. These are labeled live `/proc` observations.

Authoritative login-node `sacct` telemetry records 35 P07 allocations: 4 completed, 19 failed, and 12 deliberately cancelled after a recoverable checkpoint or before invalid promotion. The two final result jobs completed with exit `0:0`. `MaxRSS` is absent because this cluster reports `JobAcctGatherType=(null)`; missing accounting fields were not fabricated. The telemetry table preserves allocation state, elapsed time, requested memory, allocated CPUs, CPU-time capacity, node, exit code, and timestamps.

The retry history reflects observed executable/data errors, not biological exclusions:

- launch and scratch binding were repaired after Slurm spool-path and missing `SLURM_TMPDIR` errors;
- native SweepGA multiplicity was measured and an exact removal-only postfilter was added;
- authenticated IMPG CLI, boundary, temporary-path, duplicate-ID, shard-lifecycle, and lace behaviors were discovered on real data;
- the complete 203,698-shard IMPG result was rescued through a verified 1,903,441,920-byte archive rather than recomputed or redownloaded;
- a site-only IMPG rescue completed but was scientifically superseded by genotype-aware H2 projection;
- normalization/reconstruction boundary defects and consensus accounting were corrected before promotion;
- completed consensus and PSMC stages were reused through verified checkpoints;
- the native PSMC RD/TR/RS format and 100-bp scenario scaling factor were corrected before finalization.

Superseded or invalid checkpoints were moved into labeled recoverable `rejected/` directories under the canonical root. They were not silently deleted and are not used by the promoted result.

## Confidence and interpretation

The core result is valid under the authorization’s hard biological gates. Confidence remains tier C because optional covariates are unavailable: final assembly QV, BUSCO completeness/duplication, k-mer/copy-number audit, raw-read validation/chemistry, a standalone repeat report, and independent Ne evidence. Their absence is a confidence label, not a reason to erase a completed biological result.

This π is the callable divergence between the two phased haplotypes of one diploid individual. It is not a population-sample mean. PSMC is an unscaled single-diploid trajectory; the nine absolute scenarios are sensitivity transforms only. The canary establishes executable feasibility and yields a real biological estimate suitable for the remaining authorized pilot, but comparative inference must use the planned multi-species design and propagate pair-level uncertainty.

## Repository artifacts

| Artifact | SHA-256 |
|---|---|
| `analysis/vgp_real_canary_execution_v1.json` | `1a7a951c3c8bcbfac4931f9b0de99457774639e565bfd2d77f850e9334d77dbd` |
| `analysis/vgp_real_canary_variant_subset_v1.tsv` | `a60e6008e7cd183e644ba3ffa040728cf81903f15559882e7f7c32a3a16d9243` |
| `analysis/vgp_real_canary_sacct_v1.tsv` | `1b4780d46eda66cb43afd0e1363f9818d1349ede695902959ab5e9990deb17ae` |
| `analysis/vgp_real_canary_promotion_v1.json` | `4d9afaa4371d626c5d9337bd8d46842cd2a4e75437efd3f3c4e16c66bec93157` |

The JSON execution manifest is the machine-readable handoff. It records the canonical VGP root and links the promoted core, telemetry, exact annotation, independent subset, environment capture, mapping, IMPG, variants, masks, consensus, diversity, PSMC, and verification result.
