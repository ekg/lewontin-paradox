# Tier 3c GC3 control audit and promotion decision

Date: 2026-07-14

Decision: **promote the frozen Tier 3c composition dataset**

Scope: audit and selectively rerun only *D. melanogaster* and *H. sapiens*;
do not recompute the other 133 completed exact-reference analyses.

## Executive finding

The production estimator is valid.  A separate implementation that imports
none of the production estimator code exactly reproduced the GC numerator,
callable-third denominator, retained-gene count, and floating-point ratio for
both controls.  Equal-gene weighting does not recover the old anchors: it
makes GC3 higher for both species.  The failure therefore came from
unsupported, non-comparable original anchor values and bands, not from CDS
phase, strand, genetic code, canonical-transcript selection, pseudogene
filtering, weighting, memory, or timeout behavior.

The original bands remain recorded as failures.  They were not widened or
rewritten.  The post-hoc replacement gate is disclosed in
`analysis/tier3c_control_audit.json`: exact cross-implementation reproduction
from checksum-frozen native inputs, continued provenance/reconstruction hard
gates, and definition-aware consistency with published gene-level GC3.

## Independently reproduced statistics

| control | exact assembly | production and independent pooled GC3 | G/C thirds / callable thirds | retained genes | independent equal-gene mean (median) |
|---|---|---:|---:|---:|---:|
| *D. melanogaster* | GCF_000001215.4 | 0.6312610989373014 | 4,667,727 / 7,394,289 | 13,872 | 0.6404801588873662 (0.6514671924228488) |
| *H. sapiens* | GCF_000001405.40 | 0.5845556551887756 | 7,300,598 / 12,489,141 | 22,162 | 0.6060410199246469 (0.6163579188107490) |

“Pooled” means `sum(G-or-C third bases) / sum(callable third bases)` across
one retained canonical CDS per gene, so longer CDS contribute more positions.
“Equal-gene” means the arithmetic mean of each retained gene's own canonical
CDS GC3, so every gene contributes one value regardless of CDS length.  The
primary Tier 3c statistic remains the predeclared pooled statistic; the
equal-gene result is diagnostic and is not substituted into the data table.

## Independent method audit

`analysis/tier3c_control_audit.py` is a clean implementation: it does not
import `analysis/tier3_common.py` or `analysis/tier3c_ncbi_gc.py`.  It builds
its own read-only FASTA index, parses GFF3 attributes and sequence-region
declarations, reconstructs CDS, selects transcripts, and computes both
weightings.

The audit independently enforced the following rules.

- CDS segments are ordered in biological 5′→3′ order.  Minus-strand segments
  are reverse-complemented.  Phase is trimmed only at the biological start,
  and every later phase is checked against cumulative frame.
- Nuclear genetic code 1 is required.  One terminal stop is excluded and an
  internal stop, ambiguous base, inconsistent phase, overlap, duplicate
  segment, cross-contig transcript, or frameshift invalidates the candidate.
- A unique provider canonical tag is honored.  Otherwise the longest valid
  translated CDS is selected, with bytewise stable-ID tie breaking.
- `pseudo=true`, `exception`, and `transl_except` transcripts/CDS are excluded
  before canonical selection.  Exclusion counts exactly match production:
  Drosophila has 819 translation/pseudogene exceptions and 101 genes without
  valid CDS; human has 6,061 exceptions, 2,022 genes without valid CDS, and
  two invalid candidate transcripts.

The selected gene→transcript mappings are digest-locked in the audit JSON.
Every retained CDS—not only a sample—was reconstructed from the exact FASTA.

## Exact provenance hard gates

Both primary results remain native to the exact reference accession including
version.  No liftover, congener annotation, or de novo prediction is used.

| control | provider/release | FASTA SHA-256 | GFF SHA-256 | contig-map SHA-256 | dictionary audit |
|---|---|---|---|---|---|
| *D. melanogaster* | NCBI RefSeq, 2014-08-08 | `bee09429e48a8e2ce18a30c4f28641693d9a137d5d18f70cab4efb28204d513a` | `662dc13c8c1b1559e50096ae5196cb542b98f00c4ad23079d4245294e9c1b909` | `4152ecde715fbb58c6b7a30bcbc83aa17f7b2091f9de4896fff7490ecee2c1e4` | 1,869 FASTA contigs = 1,869 mapped GFF regions; 11 CDS contigs |
| *H. sapiens* | NCBI RefSeq, 2022-02-03 | `4b16731cc54ebd43f5ea278d01833fc072f0b3c4d8c4bbf2af9d05e06df05fa3` | `f80e5495c58831d23328b2fe593c32d91ba7d569166cd7bd96ce555ed089901b` | `94ec334d3a9dcf3120b2372c1a12300eba204866b2d3d84822b7f4046f7eccbd` | 704 FASTA contigs = 704 mapped GFF regions; 421 CDS contigs |

The committed control QC JSON retains provider, release, accession/version,
FASTA/GFF checksums, full mapping, genetic code, native status, dictionary
validation, and CDS audit.  Projected annotation remains sensitivity-only by
policy; missing native annotation remains explicit GC3 unavailability.

## Anchor trace and diagnosis

Git history traces “dm6 GC3 ≈ 0.55; hg38 GC3 ≈ 0.52” to commit `0846f05`
in `analysis/TIER3_EXECUTION.md`.  That assertion says “published values” but
contains no citation, gene/transcript set, provider release, terminal-stop
policy, or pooled-versus-gene weighting definition.  Commit `6d1a05c` copied
those assertions into the original bands `[0.50, 0.60]` and `[0.47, 0.57]`.
There is therefore no recoverable evidence that the anchors measure the same
statistic as Tier 3c.

The original failures are preserved:

- *D. melanogaster*: 0.6312610989373014 is outside `[0.50, 0.60]`.
- *H. sapiens*: 0.5845556551887756 is outside `[0.47, 0.57]`.

Published gene-level results support the audited scale and contradict treating
0.60/0.57 as safe upper limits:

- Campos et al., [Codon Usage Bias and Effective Population Sizes on the X
  Chromosome versus the Autosomes in *D. melanogaster*](https://academic.oup.com/mbe/article/30/4/811/1066476),
  report mean per-gene GC3 of 0.688 on X and 0.641 on autosomes (95% CIs
  0.683–0.692 and 0.639–0.643).  The audit's equal-gene mean is 0.64048.
- Payne and Alvarez-Ponce, [Codon Usage Differences among Genes Expressed in
  Different Tissues of *D. melanogaster*](https://pmc.ncbi.nlm.nih.gov/articles/PMC6456009/),
  report an all-gene mean of 0.65 and median of 0.66 across 13,088 genes.
- Mitra et al., [Mapping the inter- and intra-genic codon-usage landscape in
  *Homo sapiens*](https://academic.oup.com/nargab/article/8/1/lqag024/8503857),
  use canonical CDS on GRCh38.p14 and report mean/median GC3 of 0.59 with IQR
  0.44–0.73 across about 20,000 genes.
- Gaiti et al., [GC Content of Early Metazoan Genes and Its Impact on Gene
  Expression Levels in Mammalian Cell Lines](https://pmc.ncbi.nlm.nih.gov/articles/PMC5952964/),
  report a mammalian coding-sequence GC3 mean of 0.59 and describe the human
  mean as above 0.55.

These references use different annotation releases and primarily gene-weighted
summaries, so they are scientific context rather than frozen numeric truth.
That is precisely why the replacement gate uses exact independent reproduction
instead of inventing a wider numeric band after seeing the answers.

Diagnosis: estimator wrong—**no**; different weighting explains failure—**no**;
original bands too narrow or misanchored—**yes**.

## Selective rerun and preservation

No method or input changed, so a full recomputation would add scheduler load
without adding scientific evidence.  Only both controls were rerun with the
unchanged production implementation after the audit.  The rerun proof in
`analysis/tier3c_rerun_identity.json` records byte-identical scientific JSON:
Drosophila SHA-256 `3d361b5ccceb4124785db46ca4a9a74d7794e84a66830edf16f596e4bb78deb0`
(job 1755843) and human SHA-256
`6255e08ff903f8d27f76f75b5e7a761579d995e1516c5aab7d3e43d71562a5a0`
(job 1755841).  Operational sidecars changed because they record the new
invocations.

The preserved tables have these SHA-256 values and row counts (header included):

- `analysis/tier3c_data.tsv`: 135 results, 136 lines,
  `443df4891ebc3dbe35b7c0ee92f11731b8ae6eb4fa85fe5f7c5b847155145aee`.
- `analysis/tier3c_manifest.tsv`: 135 exact-reference rows, 136 lines,
  `d8a600e1632f49caad12005d83e1dc2ae599e25283b0a398e168b317180d1961`.
- `analysis/tier3c_failure_ledger.tsv`: 38 reproducible missingness rows, 39
  lines, `62580d9530f6859d3cd40015822c8944dcde3ed2f1d4ac1c70249274039d5750`.

The 38 failures remain 36 species without an exact current same-species NCBI
assembly plus two invalid nuclear host records.  None is converted to a result
by projected annotation or relaxed provenance.

## Scheduler and environment audit

Retry scripts now freeze two non-constrained lanes:

- standard: 2 CPUs, 32 GiB, 2 hours, array throttle no greater than eight;
- observed/predicted outlier: 2 CPUs, 64 GiB, 4 hours, single-task throttle.

The independent audit ran as Slurm arrays 1755822 and 1755840; the final
unchanged production control rerun ran as array 1755841.  All use the
GC-rooted Guix profile
`/gnu/store/z9v2f6faha9cwjz0sm5iphhlzisgi077-profile`, channel commit
`44bbfc24e4bcc48d0e3343cd3d83452721af8c36`, and no node constraint.

`analysis/tier3c_scheduler_telemetry.tsv` joins `sacct` requested memory,
elapsed time, state, exit code, CPUs, time limit, and node to process-level
`getrusage` MaxRSS.  The cluster has `JobAcctGatherType=(null)`, so `sacct`
exposes a `MaxRSS` column but leaves it blank; the report records this limitation
explicitly and never substitutes requested memory for measured RSS.  Every row
must be `COMPLETED` with exit `0:0`, and no `OUT_OF_MEMORY` or `TIMEOUT` state is
accepted.

The telemetry ledger contains 139 records: 135 frozen primary jobs plus two
independent audit jobs and two selective checksum reruns.  All 139 are
`COMPLETED` with exit `0:0`; none is OOM or timed out.  The four retry records
show 32 GiB and 02:00:00 requested.  Their measured process MaxRSS/elapsed are
100,108 KiB/00:00:54 and 572,928 KiB/00:03:00 for the independent audits, then
438,916 KiB/00:00:50 and 9,411,856 KiB/00:10:03 for the production checksum
reruns.  The ledger SHA-256 is
`a1b8b7bcdc059ca542dd8bffaffc8fd40dd17192c3afe95f0eba2fe8afa71cf3`.

Scheduler headroom played no role in promotion.  Promotion follows exact
estimator reproduction, native provenance, definition-aware controls, and
unchanged scientific bytes.
