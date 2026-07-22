# Clean VGP P07 canary

Status: **COMPLETE — real clean Slurm result**

Run ID: `vgp-clean-canary-20260722-v1`  
Selection: P07, *Spinachia spinachia*, individual `fSpiSpi1`  
Durable result: `/moosefs/erikg/vgp/pilot/clean-canary/vgp-clean-canary-20260722-v1/P07`  
Slurm job: `1791510`, `COMPLETED`, exit `0:0`, elapsed `04:23:58`, node `octopus02`

## Result

Assembly H1/H2 divergence is the primary heterozygosity result. The clean run
found 574,122 heterozygous SNPs over a final consensus-callable denominator of
267,379,237 bp:

```text
pi = 574122 / 267379237 = 0.0021472198306856562
```

This is exactly the value predeclared from the earlier P07 result. The absolute
pi difference is `0.0` against tolerance `1e-12`; the callable-bp difference is
`0` against tolerance `0`. There is therefore no primary-result difference to
explain. The primary PSMC result is also byte-identical to the earlier P07
primary (`bd1f85d3e84e6ed3065d3329fc5556289341a0e1b8571d0922a302d143e3dd89`).

Reads were not used. They remain optional symmetric sensitivity evidence and
did not veto, select, or modify the assembly result.

## Selection and immutable inputs

The selection was committed before execution in
`analysis/vgp_clean_canary_selection_v1.json` at `2026-07-22T18:20:00Z`.
P07 was selected because it is the smallest authorized P04/P07 pair by summed
compressed input size and has an exact native H1 accession/dictionary annotation.
P04 is larger and its available catalog annotation binds to
`GCF_017639655.2`, not its H1 GCA dictionary.

The exact pair is H1 `GCA_048126635.1` and H2 `GCA_048127205.1`. Immutable CAS
sources were converted once to indexed BGZF by the pinned build command:

```sh
bash analysis/build_vgp_clean_canary_bgzf.sh
```

The resulting source manifest has SHA-256
`21da7542d6ca19b6f8a0c7e35cca1434a956a70fbcc2a78c41e2e325a5325886`.
The staged BGZF inputs are:

| Side | BGZF SHA-256 | Source CAS SHA-256 | Compressed bytes |
| --- | --- | --- | ---: |
| H1 | `2f23c4d60de028f3107658ea472212ae3c73e0ca52b01785f1ce1132ed3564f1` | `2b758e606304f7cb5e795d7939979b08c21bf4f3eac7ea3cf1c6ab0a463733c7` | 129,453,729 |
| H2 | `6e1e22f734c007285beb187dced433baba15dc7d9c2e7d0b34b46f1538d50092` | `6e7e9d88b88a3d030d80009191a06c525fbcc50044db05da4d81a8e9ad97ed40` | 128,431,164 |

The exact annotation is `GCA_048126635.1-GB_2025_08_04`, bound to assembly
`GCA_048126635.1` with equal native sequence dictionaries. Its GFF object is
`8f640543accd8081d1b7048eda32c9f1eef33b02f321b7b0f8adcf3b01dd6838`.

## Clean Slurm execution and scratch proof

The submitted command was:

```sh
sbatch \
  --output=/moosefs/erikg/vgp/pilot/clean-canary/logs/vgp-clean-P07-%j.out \
  --error=/moosefs/erikg/vgp/pilot/clean-canary/logs/vgp-clean-P07-%j.err \
  analysis/slurm/run_vgp_clean_canary.sh
```

The durable target was required not to exist. The runner created a unique
private node-local root, staged and expanded the BGZF inputs there, exported
`TMPDIR`, `TMP`, and `TEMP` to that root, changed directory into it, and ran all
mapping, IMPG, variant, consensus, annotation, and PSMC work without a resume
or checkpoint path. `job.json` records all three temporary variables as
`/scratch/vgp-clean-P07-1791510-KMUN0E` and
`prior_intermediates_reused=false`.

SweepGA was invoked with `--num-mappings 1:1`. FastGA ran below its own private
`/scratch/vgp-map-P07-1791510-sXWE6K` tree. The two-second live guard collected
889 FastGA snapshots. Every observed cwd, temp file, pair descriptor, index,
alignment descriptor, and managed open path resolved below that node-local
tree; `fastga_scratch_contract.json` reports `contract_valid=true`. The final
scratch contract records FastGA PID 3857 and cwd
`/mnt/sdb1/scratch/vgp-map-P07-1791510-sXWE6K/inputs`, the resolved filesystem
target of `/scratch` on `octopus02`.

Native mapping had 40,763 records and overlap depth up to 20. The deterministic
log-length/ANI-ranked bidirectional disjoint greedy filter removed 36,466 and
retained 4,297 records. Both query and target maximum overlap depth are exactly
one. The retained PAF SHA-256 is
`0bcf45366bca9217cd012242cedf79921a6cafd6b1fb1f5247db7e54512d00fa`.

Fresh IMPG work built a new graph index, split 203,698 disjoint and exhaustive
queries across 16 two-thread workers, produced 203,698 nonempty regional VCFs,
and hierarchically laced them with final exact-boundary reconciliation. No prior
IMPG index, partition, query, shard, or laced VCF was read. All 79 sequence IDs
observed in retained PAF and graph-focus records resolve in the staged H1/H2
sequence dictionary to a sequence SHA-256; unresolved IDs are zero.

The root Slurm allocation completed successfully. The full `sacct` capture also
contains short `.0`-style steps created by read-only live `srun --overlap`
inspection. Steps `.4` and `.28` exited 1 because those diagnostic one-liners
had a tail/glob or quoting miss; they were not pipeline stages and did not alter
the root allocation, which remained running and ultimately exited `0:0`.
Cluster accounting did not populate MaxRSS/TotalCPU. The committed stage
telemetry instead records 208 successful stage sentinels, start/end times,
scratch reservations and high-water bytes, retry zero, and a minimum of
849,175,187,456 scratch bytes available at stage start.

## Variants, masks, and consensus

Exact PAF CIGAR reconstruction checked all 4,297 retained alignments and
1,488,893 raw variant records with zero H2 reconstruction failures. Normalized
VCF and BCF were regenerated, indexed, and compared independently in the three
strata below. The reason-coded mask accounts for the entire 407,561,107-bp H1
universe with discrepancy zero. It retains 272,818,693 bp before the consensus
indel-flank rule and 267,379,237 final callable bp after that rule.

The diploid consensus uses IUPAC for callable heterozygous SNPs and `N` for
masked bases. It includes 574,122 callable heterozygous SNPs, 371,229
heterozygous indels, and masks 5,439,456 H1 bp around indels. No non-callable
base is encoded as homozygous reference.

Important output digests are:

| Output | SHA-256 |
| --- | --- |
| IMPG index | `801ed48ccffc486671d8194270b9ee099374e1d401a312abcb4fe4a60cb55990` |
| Laced VCF | `759dd18f570da8a397167c5396f8128c0a18c5a4994f66490ea0691e35e48131` |
| Normalized VCF.gz | `29117558e697a27f7301d32d7c45ae7387560193f545595d197c7154508c7a63` |
| Normalized BCF | `fbc4c1141698cca2af419cc4e7f2774e3cf14ff97bd4f876d18412843a5b91b9` |
| Callable BED | `fbf1c7f9dc214f90a6829f4266351cd40c0afa4e2cda1e685f8331f41125699a` |
| Mask reconciliation | `793a5af347fe4881eb1617f98821523b3c510b975dde336f1ce3e57dd4435394` |
| Diploid consensus | `d6dc1b0c1b500942464f962ec643219a41b885e1b142adc61457d5f6344fc98e` |
| Primary PSMCFA | `2e68a0077ed3701e069eecf1c11f6449fc2a7317d807b4417c17ffa14287197a` |

`analysis/vgp_clean_canary_output_digests_v1.txt` records SHA-256 for all 1,103
files in the 1.1-GiB promoted result tree, not only this concise set.

## PSMC and scaling

The unchanged primary PSMCFA contains exactly 4,075,629 bins:
110,846 `K`, 2,362,279 `T`, and 1,602,504 `N`. The 109 contig-bounded 5-Mb
units preserve those counts exactly and never cross a contig. The deterministic
draw manifest covers exactly replicates 1 through 200 and has SHA-256
`6de6ad6568418e6a5384afac874d344785157729aae9815722f79531e12f6c67`.

All 200 bootstrap PSMC fits are finite and have 64 intervals. Primary theta is
`0.032921` per 100-bp bin. The predeclared nearest-index equal-tail diagnostic
is `[0.028562, 0.039085]`; primary theta lies inside it, so the 200-fit set is
centered and passes. Bootstrap theta summary is minimum `0.027292`, median
`0.033505`, and maximum `0.04099`.

The earlier repaired P07 diagnostic had central bounds
`[0.028314, 0.038858]` and median `0.0333045`. The PSMCFA and deterministic
outer draw manifest are byte-identical between runs, as is the primary fit.
The small bootstrap-distribution difference comes from the pinned PSMC
bootstrap mode (`-b`), which resamples the already frozen block records without
a seed option. Both independently generated 200-fit sets are finite and center
the identical primary. This difference does not affect assembly pi or its
callable denominator.

The primary unscaled trajectory is preserved separately. Nine predeclared
generic mutation-rate/generation-time sensitivity scenarios produce 576 scaled
trajectory rows; they are sensitivity grids, not a claimed species calibration.
The scaled-trajectory SHA-256 is
`358ed66d6b062ecf7d25f1d0bd69e11527e63804af280285140f41d4da2bdecb`.

## Exact annotation partitions

| Partition | Callable bp | Heterozygous variants | Estimate |
| --- | ---: | ---: | ---: |
| CDS | 26,185,407 | 15,849 | 0.000605260785138837 |
| fourfold | 4,230,801 | 3,586 | 0.0008475936353423382 |
| fourfold W | 1,455,109 | 1,697 | 0.0011662356565728066 |
| fourfold S | 2,775,692 | 1,889 | 0.0006805510121439986 |
| W→S | 1,455,109 | 1,346 | 0.0009250166138756615 |
| S→W | 2,775,692 | 1,318 | 0.0004748365452651087 |

The exact partition JSON SHA-256 is
`93500ecd55ea05d1b27f87c6d2b18fa01985cd1ab57756d90bf5daaf5b2facb7`.
It reports 38 sequence regions, 20,293 canonical transcripts, native sequence
dictionary equality, and 542 frame-discordant overlap positions excluded.

## Independent early/middle/late reconstruction

The post-run audit selected positions by cumulative H1 dictionary length, then
queried pinned bcftools independently against normalized VCF and BCF, intersected
the callable mask, sliced the diploid consensus and primary PSMCFA, and parsed
the exact source GFF. All three VCF/BCF comparisons are equal.

After atomic promotion, the same audit was run a second time from a newly
materialized copy of the immutable BGZF inputs, with only a read-only link to
the durable output tree. It independently reproduced pi, callable bp, all three
strata, the 200-fit centering interval, and the 79-ID graph digest ledger. Its
machine-readable record is `analysis/vgp_clean_canary_independent_reaudit_v1.json`
(SHA-256 `c4e24396532789d22d8fd6a21c36de6469fa8d50292a6ec7d8e825e034210ff9`).

| Stratum | Region | Variants | Callable bp | Consensus non-N bp | GFF overlaps | PSMC K/T/N bins |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| early | `CM106587.1:1-5000000` | 9,136 | 3,116,900 | 3,077,548 | 19,721 | 914 / 27,336 / 21,750 |
| middle | `CM106594.1:18790795-19449660` | 77,434 | 604,569 | 292,641 | 1 | 502 / 94 / 5,994 |
| late | `CM106606.1:11237258-14709901` | 10,755 | 2,456,435 | 2,431,661 | 6,068 | 2,359 / 20,708 / 11,661 |

## Pinned validation and evidence files

The repository-wide pinned GNU Guix validation command was:

```sh
analysis/slurm/guix_job.sh \
  "$PWD/analysis/pilot_results/guix_environment.json" \
  python3 -m pytest -q analysis
```

The wrapper verified the committed channel and manifest hashes, isolated the
captured profile, and reported `492 passed in 32.99s`. Focused clean-canary and
authorization tests also passed before submission; shell syntax, Python compile,
JSON parsing, and `git diff --check` passed.

Committed evidence:

- `analysis/vgp_clean_canary_execution_v1.json`: complete machine-readable result
- `analysis/vgp_clean_canary_independent_reaudit_v1.json`: cold post-promotion re-audit
- `analysis/vgp_clean_canary_job_v1.json`: private scratch and job identity
- `analysis/vgp_clean_canary_sacct_v1.tsv`: full scheduler allocation/step capture
- `analysis/vgp_clean_canary_stage_telemetry_v1.tsv`: 208 stage telemetry rows
- `analysis/vgp_clean_canary_output_digests_v1.txt`: closed result-tree digest ledger
- `analysis/vgp_clean_canary_graph_sequence_digests.tsv`: staged sequence dictionary digests
- `analysis/vgp_clean_canary_slurm_1791510.out` and `.err`: exact scheduler logs

The durable execution JSON SHA-256 is
`2454b91837e03fbebf12e6d75b129d4dfd3f918b77889e6f4d09b45a303aec78`;
the closed output-digest ledger SHA-256 is
`325811849e09b0c134a864143b8c5074fe7835f353de2827bfae696bfea7065e`.
