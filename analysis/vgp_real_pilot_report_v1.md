# Real VGP ten-pair pilot: P04 completion and closed-world accounting

Date completed: 2026-07-20 UTC
WG task: `run-vgp-real-pilot`
Authorization: `vgp10-auth-20260718-v2`
Canonical shared VGP root: `/moosefs/erikg/vgp`

## Result

**PASS for the task completion gate.** All nine authorized non-canary primaries crossed the real Slurm submission boundary, and the previously completed P07 canary is included in the ten-pair closed world. P04 (*Falco naumanni*, individual `bFalNau1`, BioSample `SAMN16870685`) is a newly completed biological pair. Together, P04 and P07 provide two positive, nonzero biological estimates. Every other pair retains either its exact hard-input failure or a live corrected infrastructure retry; no alternate was activated.

P04 used H1 `GCA_017639655.1` and H2 `GCA_017639645.1`. Its canonical run is:

`/moosefs/erikg/vgp/pilot/runs/vgp10-auth-20260718-v2-pilot-v1/P04`

All new canonical data and manifests are rooted at `/moosefs/erikg/vgp`. The legacy Lewontin-paradox-named VGP directory was not used as a new canonical destination.

## P04 biological result

P04 completed the frozen whole-assembly workflow with SweepGA `--num-mappings 1:1`, exact removal-only bidirectional multiplicity enforcement, IMPG partition/query/lacing, exact PAF-derived VCF and BCF, ordered masks, a diploid IUPAC consensus, callable diversity, an unscaled PSMC trajectory, 200 boundary-aware bootstraps, and nine separate sensitivity-scaling scenarios. P04 had no eligible exact annotation binding in the frozen authorization, so annotation is explicitly not applicable rather than missing.

| Quantity | Verified value |
|---|---:|
| H1 universe | 1,215,719,661 bp |
| Pre-indel callable | 781,674,063 bp (0.6429722970) |
| Final callable | 770,780,965 bp (0.6340120915) |
| Heterozygous SNPs after masks | 3,548,818 |
| Callable diversity π | 0.004604184795871289 |
| Accounting discrepancy | 0 bp |
| PSMC primary intervals | 64 |
| Finite bootstraps | 200/200 |
| Scaling scenarios | 9 labels, 576 rows |

The diversity estimator is exactly `3,548,818 / 770,780,965`: final callable heterozygous SNPs per final callable H1 base. It is phased-haplotype divergence from one diploid individual, not a population-sample mean.

SweepGA emitted 141,476 native records. The deterministic length/ANI-ranked removal-only filter retained 6,974 records, with independently checked maximum query and target overlap depth both equal to one. Exact CIGAR reconstruction covered 780,758,516 match, 4,219,233 mismatch, 1,196,495 H1-deletion, and 1,204,593 H2-insertion bases; all 6,974 PAF records reconstructed and no record failed.

IMPG produced and audited 608,347 nonempty exact regional VCFs totaling 1,148,891,894 bytes across 407,971 unique native partition identifiers. The durable continuation reused all 608,347 queries and all 24 lace shards; it did not recompute mapping, query extraction, or lacing. Job `1783681` completed site projection and boundary reconciliation in `06:29:09`, producing a 484,058,926-byte laced VCF. The original post-lace decoder failure and every earlier P04 attempt remain in the execution ledger.

Core P04 scheduler outcomes were mapping `1782075` (`00:14:55`, 32 CPUs/160 GiB), IMPG continuation `1783681` (`06:29:09`, 48 CPUs/160 GiB), variants `1783682` (`00:02:55`, 24 CPUs/160 GiB), consensus `1783683` (`00:30:20`, 24 CPUs/160 GiB), PSMC array `1783684`, and finalize `1783685`. The array retained one primary plus 200 independent bootstrap outputs. Per-task resources stayed at 4 CPUs/16 GiB/8 hours; only the array throttle changed from 20 to 21 to dispatch the already-authorized final replicate without canceling or requeueing an allocation.

## Independent reconstruction

The independent audit rebuilt the first-reason-wins mask from raw reason BEDs, required exact emitted-BED equality, queried the BCF independently with pinned bcftools, matched ordered VCF/BCF alleles, rechecked exact PAF multiplicity, and recomputed early/middle/late 5-Mbp dictionary strata. It also verified PSMCFA symbols, all finite bootstraps, the unscaled primary, scenario labels and source, sentinel digests, and canonical roots. Its subset is committed as `analysis/vgp_real_pilot_P04_variant_subset_v1.tsv`; the detailed machine result is `analysis/vgp_real_pilot_P04_execution_v1.json`.

## Closed-world ten-pair disposition

| Pair | Current retained outcome |
|---|---|
| P01 | Hard primary IMPG failure: exact census encountered frozen graph sequence `CM070280.1` absent from both staged FASTAs; completed mapping retained, no region omitted. |
| P02 | Hard primary IMPG failure: exact census encountered frozen graph sequences absent from both staged FASTAs; completed mapping retained, no successful-shard substitution. |
| P03 | Hard post-correction execution failure: unchanged primary retry `1787016` reproduced child loss at 37:25 after live proof that cwd, all temp variables, indexes, and 772 pair descriptors remained in private `/scratch`. |
| P04 | **Complete:** π 0.004604184795871289, 770,780,965 final callable bp, 200/200 finite bootstraps, unscaled trajectory and nine scenarios. |
| P05 | Hard primary IMPG failure: frozen graph sequence IDs absent and a mixed-site variant type rejected; no omission or alternate. |
| P06 | Original FastGA pair-file reasons retained and reclassified retryable infrastructure; unchanged primary retry `1787022`. |
| P07 | **Complete canary:** π 0.0021472198306856562, 267,379,237 callable bp, 200/200 finite bootstraps, exact annotation. |
| P08 | Job `1782082` canceled after explicit live scratch audit because the launch predated the full cwd + `TMPDIR`/`TMP`/`TEMP` contract; no durable alignment existed. Unchanged primary retry `1787028`. |
| P09 | Original FastGA pair-file reason retained and reclassified retryable infrastructure. Corrected job `1787035` stopped before input staging at the 600-GiB measured-free `/scratch` gate on octopus08; unchanged primary retry `1787065` excludes the three nodes with measured headroom/I/O failures. |
| P10 | Original FastGA pair-file reason retained and reclassified retryable infrastructure. Corrected job `1787041` stopped before input staging at the 600-GiB measured-free `/scratch` gate on octopus08; unchanged primary retry `1787071` excludes the three nodes with measured headroom/I/O failures. |

The machine closed world is `analysis/vgp_real_pilot_closed_world_v1.json`. It requires exactly P01–P10, verifies that every non-canary primary has real mapping/IMPG/variant/consensus/PSMC/finalize submissions, rejects a zero-estimate packet, and only labels a pair `failed_primary` when its retained failure remains explicitly terminal.

## Corrected FastGA scratch contract

P08 job `1782082` was inspected before cancellation. FastGA PID 22742 had cwd `/mnt/sdb1/scratch/vgp-map-P08-1782082-UBTHLh/inputs`; every observed `.bps`, `.gix`, `_algn`, and `_uniq` descriptor was under that private tree. `/scratch` resolves to `/mnt/sdb1/scratch` on the node. However, only `TMPDIR` was exported; `TMP` and `TEMP` were absent, and the launcher had not entered the private workdir before invoking SweepGA. Because no mapping sentinel or durable alignment existed, the job was canceled at `1-07:07:34` and retained as a nonbiological infrastructure-contract attempt.

Commit `870715e` added a failing-first regression and a live fail-closed guard. Every mapping now measures free node-local scratch, creates a private `/scratch` directory, stages both FASTAs and generated indexes there, exports all three temp variables, changes cwd before SweepGA, and inspects live FastGA descendants through `/proc`. Any FastGA cwd or managed `_tmp_`, pair, index, or alignment descriptor outside the resolved private tree exits 70 before promotion. Promotion remains atomic and occurs only after successful guarded alignment; cleanup is prefix-guarded.

Real live proofs were captured on two nodes:

- P03 job `1787016`, octopus08, FastGA PID 36733: cwd `/mnt/sdb1/scratch/vgp-map-P03-1787016-74PCXU/inputs`; all three temp variables `/scratch/vgp-map-P03-1787016-74PCXU`; 772 observed `.bps`/`.gix`/`_pair` paths all inside the resolved private tree.
- P06 job `1787022`, octopus09, FastGA PID 42653: cwd `/mnt/sdb1/scratch/vgp-map-P06-1787022-Vx2NGZ/inputs`; all three temp variables `/scratch/vgp-map-P06-1787022-Vx2NGZ`.
- Corrected P08 job `1787028`, octopus07, FastGA PID 28178: cwd `/scratch/vgp-map-P08-1787028-QYOYyn/inputs`; all three temp variables `/scratch/vgp-map-P08-1787028-QYOYyn`; all 1,090 observed `.bps`, `.gix`, `_algn`, `_uniq`, and `_pair` descriptors were inside the private tree.
- Corrected P09 job `1787065`, octopus10, FastGA PID 39733: cwd `/mnt/sdb1/scratch/vgp-map-P09-1787065-V6YK76/inputs`; all three temp variables `/scratch/vgp-map-P09-1787065-V6YK76`; all 2,052 observed `.bps`, `.gix`, and `_pair` descriptors were inside the resolved private tree.

The complete live-path record is `analysis/vgp_real_pilot_fastga_scratch_v1.json`. P03/P06/P08/P09/P10 retain their measured per-genome memory, CPU, scratch, and wall-time envelopes; P08 excludes octopus11 based on measured pair-file I/O history. P09/P10 retries `1787065`/`1787071` additionally exclude octopus08 because jobs `1787035`/`1787041` proved that node did not meet their frozen 600-GiB measured-free gate. Those gate exits occurred before staging or biological computation and remain retryable infrastructure events. No arbitrary global memory or scratch cap was introduced.

## Telemetry and reproducibility

`analysis/vgp_real_pilot_sacct_v1.tsv` is a login-node `sacct` census of every submitted job plus expanded P04 array allocations. Missing `MaxRSS` values are retained as missing because this cluster has no job accounting gatherer; they are not fabricated. The canonical submission and execution-failure ledgers remain under `/moosefs/erikg/vgp/pilot/manifests/` and record the canonical root in every row/object.

At the final live snapshot (2026-07-20 15:57 UTC), corrected mappings P06 `1787022`, P08 `1787028`, and P09 `1787065` were independently running on octopus09, octopus07, and octopus10; P10 `1787071` was pending for resources with its independent descendant chain retained. Completing the task packet does not cancel these continuing primaries.

The pinned GNU Guix analysis suite passed 403 tests before final packet generation. All pilot shell entrypoints parse, JSON inputs validate, the independent audit passes, and the repository diff check is clean. The P04 result satisfies the task’s required new biological completion, positive π/callability, at least 100 bootstraps, trajectories, explicit scaling scenarios, QC, independent stratified recomputation, and scheduler telemetry.

## Interpretation boundary

P04 and P07 are real same-individual phased-haplotype estimates. They establish two usable biological pilot points and operational feasibility, but they are not population means. PSMC primaries remain unscaled. The nine absolute transformations are a predeclared generic mutation-rate × generation-time sensitivity grid, not independently estimated species calibrations. P04 has no authorized exact annotation partition; only P07’s accession-matched annotation may be used for annotation-stratified pilot interpretation.
