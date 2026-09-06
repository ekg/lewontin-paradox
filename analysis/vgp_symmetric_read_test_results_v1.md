# Symmetric read-vs-assembly validation: deep-dive results (P07)

Schema: `vgp-symmetric-read-test-results-v1` · 2026-09-06
Promoted run: `/moosefs/erikg/vgp/derived/read-validation/runs/P07-symmetric/slurm-2826073` (sha256 in `output_manifest.tsv`; the run's `.complete` marker was never written — a known artifact of the promotion path; `execution.json` records completion, node `octopus10`, job `2826073`, started `2026-09-05T05:23:40Z`, module sha256 `220507…ade5`, repository commit `8860af6`).
Re-derivation scratch: `/moosefs/erikg/tier3scratch/symread-deepdive` (command log: `COMMANDS.md`).

This document supersedes the *interpretation* of the h2-frame columns in the promoted
`symmetric_report.md` (see Finding A). Site counts in that report remain correct.

## Verdict (review defect I-1)

**The bounded P07 assembly callset is substantially vindicated as a heterozygosity
measure, with a quantified, localized artifact class.** On the 172,770-site symmetric
universe (250,487,725 bp), BAQ-free, parser-corrected read evidence shows:

- **54.6%** of assembly-diff sites are balanced heterozygotes in **both** frames
  (hifi; illumina 49.2%), with mean allele fraction 0.485–0.493 — no directional bias
  in either frame. These are real, same-individual hets.
- **20.6%** (hifi) are `H1-only` in **both** frames: reads carry the H1 allele at the
  H1 position and the H2 allele at the H2 position — each assembly is locally
  consistent with reads, but the pair are not two haplotypes of one locus. This is the
  **collapsed-paralogy / non-orthologous-linkage artifact class**.
- The artifact class concentrates in **13 high-divergence, GC-poor 1 Mb bins**
  (11 contigs): median divergence 6.75 sites/kb vs genome 0.50 (13.6×), GC 0.415 vs
  0.448, read depth normal, in-bin paralog-class fraction 0.27–0.95 vs genome-wide
  0.206, in-bin balanced-both 0.03–0.24 vs 0.546. Spearman ρ(divergence, H1-only) = 0.76.
- π accounting on the identical denominator closes the assembly-vs-reads gap:

| quantity (250,487,725 bp) | sites | π |
|---|---:|---:|
| assembly-diff (symmetric subset) | 172,770 | 6.90e-4 |
| hifi read-corroborated (balanced) | 101,982 | 4.07e-4 |
| + read-hets missed by assembly (B-direction) | +17,858 | +0.71e-4 |
| hifi read π (report, reproduced) | 114,320 | 4.56e-4 |
| **artifact class (hifi H1-only)** | **48,262** | **1.93e-4 (27.9%)** |
| assembly π minus artifact | 124,508 | **4.97e-4** |

Removing the artifact class brings assembly π to 4.97e-4 vs read π 4.56e-4 — the two
independent measurements converge within ~9% (illumina strict-calling undercounts; its
frame-pair π's agree with hifi's within 0.5%). **The prior "read/assembly discordance"
label for P07 resolves as: a real 25–30% inflation of assembly-diff π by
paralogy/collapse artifacts, not a wholesale assembly failure.**

## Finding A — two defects in the promoted h2-frame metrics (interpretation only)

1. **Pileup parser bug**: `parse_pileup_bases` counts `.`/`,` (matches the *pileup*
   reference = H2) as the site's `ref` (= the H1 allele) unconditionally. In the h1
   frame pileup-ref == site-ref, so h1-frame evidence is correct; in the h2 frame every
   H2-allele read was counted as H1-allele support. The report's h2 "contradiction"
   columns (94.1% hifi / 68.4% illumina "hom-ref") are therefore **mislabeled**: the
   counts are recoverable, the labels are not. Corrected (BAQ-free, h2-valid subset):
   85.2% balanced (hifi), 84.3% (illumina), H1-only 0.03%, H2-only 8.5%.
2. **Lift/allele-origin disagreement**: at 35.0% (hifi) / 24.0% (illumina) of observed
   sites, the H2 fasta base at the PAF-lifted position equals the *H1* allele — the
   IMPG-derived ALT allele is not carried by H2 at the position the SweepGA 1:1 PAF
   lift selects. Round-trip fidelity is exact (0 mismatches), so the lift is
   self-consistent; the disagreement is between the two H1↔H2 correspondence
   mechanisms — the same paralog structure the read evidence exposes. h2-frame
   adjudication is valid only on the ~65–76% of sites where H2-base == ALT (all h2
   numbers above use that restriction).
3. *(Minor)* **BAQ amplifies reference bias**: with BAQ (as promoted) hifi h1-frame
   balanced was 91,748 and H1-only 64,413; BAQ-free: 101,982 / 48,262. BAQ suppressed
   ~10,234 genuinely balanced sites toward H1-only. The promoted h1-frame numbers are
   therefore conservative about the artifact class, not optimistic.

## Detailed tables (BAQ-free, corrected parser, n = 172,770)

Per-frame classes:

| frame/platform | balanced | H1-only | H2-only | skewed | out-of-depth | not-observed |
|---|---:|---:|---:|---:|---:|---:|
| h1 / hifi | 101,982 (59.0%) | 48,262 (27.9%) | 26 | 17,923 | 4,569 | 8 |
| h1 / illumina | 87,560 (50.7%) | 30,103 (17.4%) | 67 | 4,038 | 23,583 | 27,419 |
| h2 / hifi (valid n=111,457) | 95,002 (85.2%) | 32 | 9,445 (8.5%) | 5,335 | 1,643 | — |
| h2 / illumina (valid n=102,379) | 86,277 (84.3%) | 68 | 5,034 (4.9%) | 2,224 | 8,775 | — |

Cross-frame (hifi, top combinations): `balanced|balanced` 94,238 · `H1-only|H1-only`
35,516 · `skewed|H1-only` 8,777 · `H1-only|skewed` 6,054 · `H1-only|H2-only` 4,071.
Allele fraction at balanced sites: h1-frame frac(H2-allele) mean 0.485 / median 0.491;
h2-frame frac(H1-allele) mean 0.490 / median 0.492 (illumina 0.492/0.493).

Read-het placement (hifi): 114,363 read-hets — 96,505 (84.4%) **at** assembly-diff
sites, 17,858 (15.6%) at assembly-agreeing positions inside the mask (B-direction),
0 outside. Illumina 91.8% / 8.2%. Same-individual phasing expectation is met: read
hets land overwhelmingly on assembly differences.

Flagged-bin anatomy (13 hifi bins; genome-wide baseline in header):

| bin | sites | H1-only | paralog-class | balanced-both |
|---|---:|---:|---:|---:|
| CM106588.1:18Mb | 4,223 | 0.703 | 0.538 | 0.035 |
| CM106589.1:18Mb | 375 | 0.952 | 0.952 | 0.048 |
| CM106590.1:4Mb | 1,218 | 0.654 | 0.567 | 0.216 |
| CM106590.1:11Mb | 3,693 | 0.457 | 0.386 | 0.027 |
| CM106592.1:2Mb | 1,536 | 0.658 | 0.331 | 0.057 |
| CM106593.1:0Mb | 5,087 | 0.667 | 0.369 | 0.027 |
| CM106598.1:17Mb | 788 | 0.882 | 0.692 | 0.080 |
| CM106601.1:16Mb | 119 | 0.908 | 0.908 | 0.084 |
| CM106602.1:0Mb | 81 | 0.951 | 0.951 | 0.049 |
| CM106603.1:11Mb | 2,342 | 0.609 | 0.266 | 0.094 |
| CM106603.1:13Mb | 381 | 0.751 | 0.751 | 0.244 |
| CM106605.1:11Mb | 865 | 0.852 | 0.852 | 0.147 |
| CM106605.1:12Mb | 960 | 0.797 | 0.751 | 0.135 |
| *genome-wide* | *172,770* | *0.279* | *0.206* | *0.546* |

Divergence (all 402,642 assembly sites): flagged median 6.75 sites/kb (max 89.5) vs
baseline 0.50. GC: 0.415 vs 0.448. Depth (hifi): 53.0× vs 52.0×; illumina 29.7× vs
35.4× — depth does not explain the pattern.

## Implications

**(a) Bounded three-pair core, P07 π = 0.00117 (270.5 Mb mask).** The symmetric-subset
artifact fraction is 27.9% (hifi). Extrapolated (the subset covers the 65% of the mask
that lifts through the 1:1 PAF; artifact bins are lift-enriched), the corrected
whole-genome estimate is ≈ 0.00084–0.0009. P07 remains the lowest-π pair, so the
P02 > P03 > P07 ordering is unchanged. The review's P07 caveat "read/assembly
discordance (A-1/A-2)" should be replaced by: "π inflated ~25–30% by localized
collapse/paralogy artifacts; corrected estimate ≈ 0.0009".

**(b) Tier 3A coding panels.** Same assembly-diff method; coding regions are
repeat/paralog-poor, so inflation is probably smaller but is unmeasured. Cheapest
check: intersect each panel's callable intervals with divergence-hotspot bins
(divergence ≥ ~5 sites/kb, GC-poor) — a **no-reads** artifact screen computable from
the assemblies alone; escalate to the symmetric read test only where panels overlap
flagged bins.

**(c) Scale-out QC gating.** Recommended gate per new pair, in escalating cost:
1. *Assembly-only triage* (minutes, no reads): per-1Mb divergence + GC scan; flag bins
   ≥ 2× genome median divergence and GC-poor. Report flagged-bin fraction as the
   collapse-risk metric.
2. *hifi symmetric gate* (~1 node-hour + 2 mappings ≈ 40 min): hifi reads → H1 and H2;
   BAQ-free parallel per-contig pileups (minutes); report balanced/paralog/H2-only mix
   and allele-fraction symmetry. Illumina adds depth QC but no class separation — skip
   for gating. Cost driver is the mappings; reuse existing BAMs where available.
3. *Full symmetric* (this study): both platforms, both frames — reserve for pairs that
   fail gate 2 or for calibration of the artifact-fraction prior.
Also fix before scale-out: the `.`/`,` parser bug (h2 frame), and record which lift
rows carry each ALT allele (or lift through the IMPG graph path) so the
allele-origin/lift disagreement (24–35%) is measurable in-pipeline.

## Reproduction

Commands and validation cross-checks: `/moosefs/erikg/tier3scratch/symread-deepdive/COMMANDS.md`.
Key cross-checks: universe = 172,770 (metrics `in_frame_callable`); promoted h1-frame
classes reproduced exactly (91,748/64,413 hifi; 86,110/31,329 illumina — 3-count delta
from duplicate positions); B-direction 17,858 vs 17,853 (5-count dup delta); π_read
4.564e-4 matches the report table; h1 fasta REF sanity 5,000/5,000; roundtrip 0.
