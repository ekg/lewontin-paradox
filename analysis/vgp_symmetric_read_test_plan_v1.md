# P07 symmetric read-vs-assembly test — run plan v1

Task context: review defect **I-1** (no read comparison symmetric across H1,
H2, and graph representations; prior P07 read evidence was H1-only and
targeted the *rejected* clean/global-lace callset). This run is the direct
execution of the symmetric redesign.

## Design principle

Reads are mapped to **both haplotypes independently**. The assembly π
denominator and SNP sites are lifted H1↔H2 through the corrected 1:1 PAF so
that both frames measure the *same genomic intervals*. No single reference is
treated as a truth oracle; every statistic exists in both frames and the two
frames must agree with each other before any conclusion is drawn against the
assembly callset.

## Exact inputs (all digests verified at submit time)

| Input | Path | Digest / check |
|---|---|---|
| HiFi reads | `/moosefs/erikg/vgp/views/accession/P07/SRR25606782/SRR25606782_subreads.fastq.gz` | CAS-verified view (21.4 Gb) |
| Illumina reads | `/moosefs/erikg/vgp/views/accession/P07/SRR30200290/SRR30200290_{1,2}.fastq.gz` | CAS-verified view (2×17.9 Gb) |
| H1 assembly | `/moosefs/erikg/tier3data/tier3a-acquisition-20260716/spinachia_spinachia_SK-2024b/h1.fna` | sha256 `438faaebe341…1aaf6ba` (re-verified in-job) |
| H2 assembly | `…/h2.fna` | sha256 `0bbd50ea5954…5fb9418` (re-verified in-job) |
| Assembly callset | `…/P07/bounded-production/variants/normalized.bcf(+csi)` | accepted bounded three-pair lineage |
| Assembly callable | `…/bounded-production/consensus/masks/callable.bed` + `join_qc.json` | π-denominator reproduced exactly: 270,531,638 bp |
| 1:1 PAF | `/moosefs/erikg/tier3data/tier3a-origin-remap-20260716/spinachia_spinachia_SK-2024b/mapping/production.1to1.paf` | SweepGA origin/main `--num-mappings 1:1` (753 rows) |
| Tool profile | `/moosefs/erikg/vgp/derived/read-validation/environment/profile` | minimap2 2.24-r1122, samtools 1.14, bcftools 1.14, bedtools 2.30.0 |

## Method decisions

- **Mapping**: minimap2 `-ax sr` (Illumina, both mates) and `-ax map-hifi`,
  `--secondary=no`, RG-tagged, identical parameters for both reference frames.
  Sorted+indexed BAMs are **retained** this time (old run discarded them).
- **Masks**: identical depth-mask suite (`dp5_100, dp10_60, dp10_80,
  dp15_80, dp20_80, dp10_100`), computed from Illumina depth on each frame's
  π-callable bed (inherited convention; `dp10_80` primary, Q≥20/mapQ≥20).
  Both platforms' calls are restricted to the same per-frame `dp10_80` mask so
  π denominators are identical within a frame.
- **Depth caps**: Illumina 10–80 (inherited), HiFi 10–120 (inherited).
- **Calls**: `bcftools mpileup -q 20 -Q 20` → `call -m` → `norm -m -any`,
  filtered QUAL≥30, TYPE=snp, N_ALT=1 per frame and platform.
- **Transfer semantics** (`vgp_symmetric_read_test.py transfer`): CIGAR-exact
  lifting (bisect over per-row prefix sums), strand-aware allele
  complementation, per-position unique ownership on **both** sides. SweepGA's
  1:1 cap permits residual block overlap, so positions in overlap zones
  (~5–6% of bp) are excluded rather than arbitrarily assigned; sites inside
  alignment indel gaps are excluded (structurally unprojectable).
- **Symmetric subsets**: the transfer emits both the lifted H2-frame
  bed/sites and the exact H1-frame sub-bed/sub-sites that lifted, so both
  frames measure identical genomic intervals (this corrects the prior
  design's one-frame privilege).
- **Metrics** per frame×platform: π_read on the dp10_80 mask; genotype
  concordance at assembly SNVs (pileup classification); both-direction
  contradiction rates (A: assembly-het/read-hom-ref; B: read-het QUAL≥30 at
  callable positions where the assembly callset is hom-ref); per-contig and
  per-1Mb-bin localization with bin flagging (≥50 sites, rate ≥2×
  genome-wide). The assembly site universe is **intersected with the frame's
  π-callable bed** (`--frame-callable-bed`; strata counts emitted in the
  metrics JSON and report), so concordance denominators and direction-B
  membership are comparable across frames and against the inherited H1-only
  baselines (Illumina 0.501 / HiFi 0.532).
- **Provenance**: reads verified against
  `analysis/vgp_validation_reads_manifest_v1.json` (`local_sha256` +
  `observed_bytes`); PAF, assembly BCF+CSI, callable.bed, and join_qc.json
  pinned to sha256 constants and re-verified after staging; every consumed
  input recorded in the run's `input_manifest.json` (mismatch = hard fail).
- **BAM-reuse resume mode**: `SYMREAD_REUSE_BAMS=<promoted run>/bams`
  digest-verifies the four sorted BAMs against the source run's
  `output_manifest.tsv`, skips mapping, and reruns transfer+masks+calls+
  pileups+metrics into a **new** run directory with its own manifests.

## Pre-flight results (re-executed 2026-09-04 after the review fix pass)

- π-denominator reconstruction: **270,531,638 bp == join_qc.json** ✓
- Assembly SNP sites extracted: 402,642 (TYPE=snp, N_ALT=1)
- PAF transfer: **256,392,804 bp lifted (94.77% of callable bp)**,
  210,404 sites lifted (52.3%).
- Interval-lift blocker fix, real-data effect: the previous pre-flight
  (256,642,590 bp, 94.87%) over-lifted ~250 kb of intervals that straddled
  seam/interior-overlap zones; the fixed position-ownership lift rejects
  them. Site counts are unchanged (point queries were never affected).
- `bp_delta_source_minus_dest` = **-1,082,714 bp** (destination exceeds
  source inside lifted intervals: net insertion of H2 relative to H1 across
  lifted blocks; see interpretation guide).
- Site rejections now reconcile **exactly** (per-reason counters; the prior
  pre-flight's ~16.6% "unexplained" gap was overlap-zone sites — the old
  ad-hoc estimate counted query-side bp only):

| rejection reason | sites | % of total | % of rejected |
|---|---:|---:|---:|
| overlap (query- or target-side multi-owner zone) | 126,373 | 31.39% | 65.7% |
| indel_gap (inside alignment CIGAR gap) | 28,689 | 7.13% | 14.9% |
| unaligned (no covering row) | 19,082 | 4.74% | 9.9% |
| owner_mismatch (lifted position owned by other row) | 18,094 | 4.49% | 9.4% |
| roundtrip_mismatch | 0 | 0.00% | 0.0% |
| **sum / rejected / lifted** | **192,238 / 210,404** | 47.74% / 52.26% | 100% |

  Every emitted site round-trips by construction; every exclusion carries a
  structural reason. The symmetric H1-frame subset mirrors exactly the lifted
  set, so both frames measure identical intervals and identical site
  universes.
- Unit tests: **16/16 pass** under the pinned Guix channel (`guix time-machine
  -C analysis/guix/channels.scm -- shell -m analysis/guix/manifest.scm
  --pure`), covering transfer correctness (strand flip + complementation,
  reciprocal rejection of inconsistent inverse rows, overlap-zone exclusion,
  **interior-overlap interval rejection, seam-stitch rejection**,
  rejection-reason accounting, bp_delta), metric math (π, both-direction
  contradictions, **frame-callable site-universe strata**), and bin flagging.

## Expected outputs (canonical promotion target)

`/moosefs/erikg/vgp/derived/read-validation/runs/P07-symmetric/slurm-<jobid>/`:
`transfer/` (lifted + symmetric beds/sites + stats), `bams/` (4 retained
BAM+BAI), `h1/` & `h2/` (masks, per-platform calls, assembly-evidence,
metrics.json + .bins.tsv + .contigs.tsv), `symmetric_report.md`,
`execution.json`, `output_manifest.tsv`, telemetry.

## Interpretation guide (review I-1 decision logic)

- read π(H1-frame) ≈ read π(H2-frame) ≈ assembly π (0.0011704) within
  binomial-scale error, and contradiction rates symmetric-low in both frames
  → bounded callset **vindicated**; scale-out GO.
- H1-frame vs H2-frame read π disagree materially → projection/representation
  artifact; localize via flagged bins before any biological claim.
- Both frames contradict the assembly in the *same* regions → genuine
  assembly defect; named contigs/bins become exclusion candidates.
- `bp_delta_source_minus_dest` (−1,082,714 bp on this chain) is the net
  target-minus-query length difference accumulated inside lifted intervals
  from internal CIGAR indels — not an error. It is expected to be nonzero
  whenever lifted blocks span indels; a materially *different* bp_delta
  between the forward and reverse lift of the same bed would instead indicate
  an asymmetry bug (the transfer enforces round-trip identity, so this should
  never appear; investigate if it does).
