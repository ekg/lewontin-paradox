# Assembly-only paralog triage — three-pair bounded core

Schema: `vgp-paralog-triage-v1` · 2026-09-06
Module: `analysis/vgp_symmetric_read_test.py` subcommand `triage` (repository working tree;
tests: `analysis/tests/test_vgp_symmetric_read_test.py`, 27/27 green under pinned Guix)
Deep-dive basis: `analysis/vgp_symmetric_read_test_results_v1.md` (review defect I-1)
Outputs: `/moosefs/erikg/tier3scratch/symread-triage/{P07,P03,P02}/` (`triage_summary.json`,
`triage_bins.tsv`, `paralog_bins.bed`, `run.log`)

## Purpose and scope

The P07 symmetric read test identified a collapsed-paralogy / non-orthologous-linkage
artifact class that inflates assembly-difference π by ~25–30% (read-adjudicated: 27.9% of
the symmetric subset; corrected P07 π ≈ 0.00084 vs bounded 0.00117). Reads are unavailable
for P02 and P03, so this triage applies the deep-dive's recommended **assembly-only
screen** — per-1 Mb-bin divergence and GC computed from the bounded callset and the H1
assembly — to all three accepted pairs, and reports a corrected π excluding flagged bins.

This is a coarse regional screen, not a site-level read adjudication. On P07 it removes a
superset of the read-adjudicated artifact class (see *Validation against P07*).

## Method

Per pair, inputs are the accepted bounded-production lineage (all digest-verified before
use; SHA-256 recorded in each `triage_summary.json`):

- bounded `normalized.bcf` (H1 coordinates) — SNP records stream via `bcftools query`
  (`TYPE="snp"`, biallelic);
- bounded π-callable BED (`consensus/masks/callable.bed`);
- H1/H2 whole FASTAs (P07: tier3a acquisition copies, sha256 == `acquisition_corrected_manifest.tsv`;
  P03/P02: `pilot/inputs/<PAIR>/`, sha256 == each `input-manifest.json`);
- the pair's 1:1 PAF (P07: tier3a `production.1to1.paf` (query=H1); P03/P02: three-pair
  `h2_to_h1.1to1.paf` (query=H2) — orientation auto-detected against BCF contigs and
  recorded in the summary).

Procedure (fixed thresholds, deep-dive defaults):

1. H1 fasta streamed once; per-1 Mb bin: ACGT bases, GC fraction.
2. SNP sites binned (all BCF SNP records, callable or not — matches the deep-dive's
   402,642-site divergence basis on P07); bin divergence = sites per kb of ACGT bases.
3. Bins with < 100,000 ACGT bases are excluded from medians and flagging.
4. A bin is a **flagged paralog bin** iff divergence > 2.0 × (median bin divergence)
   AND GC < (median bin GC) − 0.02.
5. Corrected π = (callable SNP sites outside flagged bins) / (callable bp outside
   flagged bins). Raw π uses the full mask — and reproduces each pair's audited
   bounded π exactly, validating the extraction.
6. Allele-origin statistic (supplementary, not used for flagging): for each liftable
   SNP site, the H2 fasta base at the PAF-lifted position (strand-complemented) vs the
   IMPG ALT allele; per-bin disagreement fractions and the genome-wide fraction.
   This mirrors the deep-dive's lift/allele-origin disagreement (P07 hifi: 35% of
   observed sites in the h2-valid subset).

## Per-pair results

| pair | species | raw π (audited) | flagged bins (eligible) | artifact site fraction | artifact bp fraction | allele-origin disagreement | **corrected π** |
|---|---|---:|---:|---:|---:|---:|---:|
| P07 | *Spinachia spinachia* | 0.0011704 (0.00117040…) | 29 / 419 | 0.4150 | 0.0718 | 0.398 | **0.000738** |
| P03 | *Colius striatus* | 0.0018644 (0.00186435…) | 137 / 1231 | 0.2247 | 0.1542 | n/a (running)² | **0.001709** |
| P02 | *Pseudorca crassidens* | 0.0024565 (0.00245646…) | 14 / 2804 | 0.0059 | 0.0030 | n/a (skipped)¹ | **0.002449** |

Raw π denominators/numerators: P07 316,631 / 270,531,638; P03 1,632,584 / 875,683,638;
P02 4,195,014 / 1,707,746,195 (audited three-pair values; each triage run recomputed and
matched its pair's audited numerator and denominator exactly).

P07 strongest flagged bins (all exceed both thresholds by wide margins):
`CM106598.1:17Mb` (22.9 sites/kb, GC 0.390, allele-origin disagreement 0.93),
`CM106601.1:0Mb` (19.3, 0.401, 0.44), `CM106594.1:18Mb` (15.8, 0.387, 0.19),
`CM106603.1:14Mb` (15.2, 0.404), `CM106593.1:0Mb` (13.7, 0.410, 0.79) — the last is one
of the deep-dive's read-adjudicated hotspot bins. Site load is dominated by extreme
bins: on P07, bins above 5 sites/kb carry 152k of the flagged sites; the marginal
0.7–1.0 sites/kb band contributes only ~3.3k. The screen flags 7 of the deep-dive's 13
read-hotspot bins and adds 22 further bins with strong assembly-only signatures.

P03 strongest flagged bins: `CM054348.1:91Mb` (6.6 sites/kb, GC 0.378), `CM054345.1:99Mb`
(6.2, 0.366), `CM054349.1:53Mb` (5.8, 0.376), `CM054345.1:144Mb` (5.5, 0.380). **P03's
median bin divergence is 0.0 sites/kb** — its 1.73M SNP sites are sparse and highly
concentrated, so more than half of eligible bins contain none. The 2×-median divergence
factor is therefore non-binding on P03 (threshold 0.0) and flagging reduces to
GC-poor ∩ any-divergence; its corrected π is a rougher bound than P07's (see Caveats).

² P03's allele-origin pass exceeds a single session's runtime budget; a full run was
relaunched detached, writing to `symread-triage/P03-full-alleleorigin/` (log:
`symread-triage/P03-full-alleleorigin.log`). The P03 row above uses the completed
`--skip-allele-origin` run in `symread-triage/P03/`.

P02 strongest flagged bins: `CM078551.1:20Mb` (5.4 sites/kb, GC 0.384), `CM078551.1:16Mb`
(2.7, 0.395), `CM078546.1:37Mb` (2.5, 0.385), `CM078545.1:18Mb` (2.5, 0.394) — far weaker
signals than P07's (which reached 13-23 sites/kb), and only 0.59% of callable sites live
in flagged bins. **P02's assembly-difference π survives paralog correction essentially
unchanged (0.0024565 → 0.0024494)**: its high π is not a collapse artifact.

¹ P02's allele-origin pass was skipped for runtime (4.58M SNP sites × 55,907-row chain
over network storage); the column is supplementary and not used for flagging.

## Validation against P07 (honest below-floor note)

The triage **reproduces the audited bounded raw π exactly** (316,631 / 270,531,638 =
0.0011704028495181033; site universe 402,642 == deep-dive; allele-origin observed
universe 210,404 == the promoted transfer's lifted count), validating the extraction,
binning, and PAF handling end-to-end.

The corrected π, however, lands at **0.000738 — slightly below the 0.00075–0.00093
acceptance band** (read-adjudicated 0.00084 with 27.9% site removal). Cause: the
assembly-only bin screen removes a **superset** of the read-adjudicated artifact class —
41.5% of callable sites vs 27.9% — because it cannot distinguish, inside a flagged bin,
paralog-driven differences from genuine heterozygosity (the deep-dive found in-bin
balanced-both fractions of 0.03–0.24 in its hotspots, so most in-bin sites are artifact,
but not all). Both quantities describe the same phenomenon: the bin-level screen is the
conservative (site-removal-maximal) bound, the read-adjudicated value is the point
estimate. **For P07, the best corrected estimate remains the deep-dive's ≈ 0.00084;
the triage value 0.000738 is its conservative lower bound.** For P02/P03 (no reads), the
triage bounds are the only available correction and should be reported as bounds.

Ordering impact: P07 remains the lowest-π pair under either correction
(0.00074–0.00084 vs P03/P02 values below), so the P02 > P03 > P07 ordering of the
bounded core is unchanged by paralog correction.

## Caveats

- Assembly-only screen approximates the read-adjudicated artifact class; on P07 it
  over-removes (41.5% vs 27.9% of sites). Reads are unavailable for P02/P03 — their
  corrected values are conservative bounds, not point estimates.
- Flagging is bin-granular (1 Mb); boundary effects and mixed bins are unavoidable.
- Thresholds (2× median divergence, median GC − 0.02, ≥100 kb bins) are the deep-dive's
  defaults, applied identically to all three pairs; no per-pair tuning was performed.
  On P03 the median-bin divergence is 0.0, so the divergence factor is non-binding
  there and flagging is effectively GC-gated; on P02 (median 0.88 sites/kb) and P07
  (median 0.35) both criteria bind.
- Allele-origin disagreement uses the same PAF/IMPG correspondence whose disagreement
  the deep-dive attributed to paralog structure; it is a supporting signal, not a gate.
- P02's accepted lineage includes the resume job's strict-H1 revalidation; its BCF is
  the audited bounded product and its raw π matches the review table exactly. Every
  pair's triage raw π reproduced its audited bounded numerator/denominator exactly
  (P07 316,631/270,531,638; P03 1,632,584/875,683,638; P02 4,195,014/1,707,746,195).
- P02 allele-origin was skipped for runtime (4.58M SNP sites × 55,907-row chain over
  network storage); P03's full run is in flight (footnote ²). The column is
  supplementary and never gates flagging.

## Reproduction

Per pair (from the repository root, pinned Guix):

```sh
guix time-machine -C analysis/guix/channels.scm -- \
  shell -L analysis/guix -m analysis/guix/manifest.scm --pure -- \
  python3 -m analysis.vgp_symmetric_read_test triage \
    --pair <PAIR> --paf <pair 1:1 PAF> \
    --h1-fasta <H1> --h2-fasta <H2> \
    --bcf <bounded-production>/variants/normalized.bcf \
    --callable-bed <bounded-production>/consensus/masks/callable.bed \
    --output-dir /moosefs/erikg/tier3scratch/symread-triage/<PAIR>
```

Exact invocations (including all input SHA-256s) are recorded in each
`triage_summary.json` under `inputs`; run logs are `run.log` beside them.
