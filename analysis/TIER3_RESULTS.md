# Tier 3 synthesis results

Status: **executed, biologically recovered, and fail-closed**
(`tier3-decisions-v1`, recovery synthesis 2026-07-16).

The machine-readable result is `tier3_results.tsv`; `fig_tier3.pdf` and
`fig_tier3.png` are regenerated from that committed table. These results do
not establish a causal gBGC effect or saturation. They replace anticipated
claims in the original plan with the estimates, null/negative results, and
structured missingness that passed the frozen gates.

## Inputs and exact join

The synthesis joined all 135 validated Tier 3c exact-species composition rows
to the Buffalo source pinned at commit
`b8f91d5c34675733db8cae8dcab625dcbb55c30a` and SHA-256
`df559451dad94b53ba8675e09811708107a57eeb6ffe8f72b944bcbbf3a1f2eb`.
The join is one-to-one on the full scientific name. Duplicate names,
covariate disagreement, missing exact species, and congener substitution are
fatal. Buffalo supplies the census-size predictor `pred_log10_N` and
class/order/family/genus labels. This is a census-size proxy, not a measured
effective population size.

All 135 exact FASTAs have whole-genome GC. Ninety have GC3 from an annotation
recorded as native to the exact accession; 45 lack an eligible native
annotation and therefore have no primary GC3 value. A projected, lifted,
predicted, mismatched, or congener annotation can never fill those cells.
Whole-genome GC remains a separate, weaker observable and is not substituted
for GC3.

Every GC3 point row carries the annotation provider and release, exact
assembly accession+version, FASTA and GFF SHA-256, contig-map and
sequence-region/dictionary SHA-256, nuclear genetic code,
native-versus-projected status, dictionary-validation result, and the
canonical-CDS reconstruction audit. The synthesis rechecks
accession/hash/dictionary/audit agreement against the committed per-dataset QC
records before fitting.

| Buffalo class | exact FASTAs | eligible native GC3 | missing GC3 |
|---|---:|---:|---:|
| Insecta | 48 | 27 | 21 |
| Mammalia | 32 | 25 | 7 |
| Aves | 12 | 11 | 1 |
| Teleostei | 9 | 6 | 3 |
| Bivalvia | 6 | 6 | 0 |
| other 14 classes | 28 | 15 | 13 |
| **total** | **135** | **90** | **45** |

The upstream Tier 3c point table reports callable-base or callable-CDS-third
denominators but not genomic block intervals. We report those denominators and
do not invent narrow binomial intervals by treating correlated bases as
independent. That is measurement precision. Cross-species intervals instead
come from 10,000 deterministic species-bootstrap replicates, stratified by
taxonomic class for class-fixed models. They are conditional on the observed
species, Buffalo predictor, and labels; they are not measurement,
sample-selection, phylogenetic, or causal uncertainty.

## Composition fits

The effect below is the change in the named GC fraction per one-unit increase
in Buffalo `pred_log10_N`. Primary p-values are adjusted together by
Benjamini–Hochberg; sensitivities form a separate family. Null and negative
effects remain in the table.

| analysis | tier | n / available | effect | 95% species-bootstrap interval | p | BH q | interpretation |
|---|---|---:|---:|---:|---:|---:|---|
| GC3, across classes with class-fixed intercepts | native exact-assembly GC3 | 90 / 135 | +0.00892 | −0.00123, +0.02149 | 0.163 | 0.217 | positive point estimate, interval includes zero |
| GC3, within Aves | native exact-assembly GC3 | 11 / 12 | +0.01459 | +0.00781, +0.02447 | 0.000434 | 0.00217 | positive association in this retained class |
| GC3, within Insecta | native exact-assembly GC3 | 27 / 48 | +0.02232 | −0.00301, +0.05728 | 0.173 | 0.217 | positive but uncertain |
| GC3, within Mammalia | native exact-assembly GC3 | 25 / 32 | **−0.00479** | −0.01179, +0.00010 | 0.103 | 0.217 | negative result, retained; contradicts a uniform positive within-clade pattern |
| GC3, within *Drosophila* | native exact-assembly GC3 | 13 / 17 | +0.00088 | −0.00634, +0.01317 | 0.839 | 0.839 | null result; genus is only a machinery proxy |
| whole-genome GC, class-fixed sensitivity | exact FASTA GC | 135 / 135 | +0.00537 | +0.00194, +0.00976 | 0.00781 | 0.0156 | positive sensitivity, but not GC3 and not a gBGC-specific endpoint |

The predeclared concavity diagnostic fit GC3 to centered
`pred_log10_N` plus its square and class-fixed intercepts (n=90). The
quadratic coefficient was **+0.00226** (95% species-bootstrap interval
+0.000005 to +0.00510; normal-approximation p=0.0602, sensitivity-family
BH q=0.0602). Saturation predicts negative concavity; the observed point
estimate and interval are in the opposite, convex direction. A quadratic is
only a shape diagnostic, not a mechanistic saturation curve, but these data do
not support the predicted GC3 saturation shape.

The primary across-class slope was stable in the sense that its
leave-one-species-out range did not conceal a single-species reversal; the
exact ranges for every primary model are stored in the `sensitivity` column.
Leave-one-class-out ranges are also stored. They do not rescue contradictory
within-class results: the negative mammal and null *Drosophila* estimates are
reported beside the positive avian estimate rather than averaged away.

Taxonomic-class fixed intercepts are a coarse control for deep-clade
differences. They are **not** phylogenetic generalized least squares. No
frozen tree with branch lengths/covariance model was supplied, so a PGLS row
is explicitly unavailable (n=90 eligible GC3 points; effect and interval not
estimable). Buffalo has sparse recombination measurements but no predeclared
recombination-machinery classes, so no post-hoc class was invented and that
model is also explicitly unavailable.

## Recovered diversity observables

Fresh Tier 3a and Tier 3b recovery replaces the former header-only state. The
corrected tables contain 12 direct assembly rows (three identities × four
annotation/statistic rows), three conservatively derived assembly
`pi_S_over_pi_W` rows, and eight population rows (two populations × four
statistics). Every row has a biological identity, explicit eligible `n`,
positive numerator/variant count and callable denominator, a finite estimate,
exclusions, software provenance, and a stated uncertainty method and interval.

The assembly modality is an alignment-conditioned H1/H2 comparison over a
deterministically selected H1-native coding panel (`n=2` haplotypes), not a
whole-genome estimate and not population diversity. Coding-gene diversity is
0.0138517 in *Menidia menidia* (28,233/2,038,234 callable bases; 50-kb
block-bootstrap 95% interval 0.0123545–0.0153404), 0.000397957 in *Spinachia
spinachia* (803/2,017,806; 0.000329512–0.000481235), and 0.00167799 in
*Tautogolabrus adspersus* (3,551/2,116,227; 0.00123835–0.00212181). Each
interval uses 1,000 deterministic genomic bootstrap replicates. The large
range is reported as heterogeneity among three selected assembly pairs; it is
not a population-size effect.

Assembly reference-conditioned `pi_S_over_pi_W` is defined as `(S
variants/S callable fourfold sites)/(W variants/W callable fourfold sites)`.
Its point estimates are 0.672544, 0.316941, and 0.760295. Because paired ratio
bootstrap draws are not present in the published upstream table, the synthesis
reports conservative ratios of marginal bounds: 0.482732–0.932999,
0.0685349–1.70202, and 0.278040–2.12518. The latter two include one.

The population modality uses two *Anopheles coluzzii* cohorts of 20 wild
diploid individuals (40 nominal chromosomes) over the same 21-Mb region with
population-specific callable masks. Population π—the sum of unbiased pairwise
differences divided by callable sites—is 0.0123500 in AO Luanda
(184,914.842/14,972,821; 10,000-replicate chromosome-stratified 1-Mb
block-bootstrap 95% interval 0.0118634–0.0128173) and 0.0158698 in GM Wali
Kunda (237,893.579/14,990,338; 0.0153949–0.0163385). Their intervals do not
overlap. Population `pi_S_over_pi_W` is 0.742656 (0.712917–0.770875) and
0.723555 (0.691945–0.757168), with the same genomic bootstrap. Component
`pi_S` and `pi_W` intervals use a 20-individual delete-one jackknife.

Assembly, population, and exact-single-assembly composition estimates are
never pooled, renamed, or treated as interchangeable. Deposited-call
whole-genome individual heterozygosity remains unavailable. The
reference-conditioned ratio is never relabelled polarized SFS-B; SFS-B remains
deferred because v1 freezes no outgroup, ancestral-error model, demographic
model, or power threshold. Full identity, numerator, denominator, eligibility,
exclusion, uncertainty, and software fields are in
`results/tier3/recovery_evidence_ledger.tsv`.

## Claim boundary

The composition dataset establishes exact-assembly GC and GC3 measurements
for the stated denominators. The recovered diversity tables establish only
their named, modality-specific biological estimates; they do not establish a
cross-species diversity slope, polarized gBGC strength, phylogenetically
corrected effects, or conserved recombination-class effects.
The mixed class-specific GC3 directions and convex rather than concave shape
do not support a general composition-saturation claim. Mutation bias,
selection, life history, deep phylogeny, and assembly/annotation structure
remain competing explanations. Therefore this run makes **no causal claim
that gBGC explains the high-population diversity plateau**.

## Reproduction

Run fitting and the complete analysis tests in the pinned pure Guix
environment:

```sh
guix time-machine -C analysis/guix/channels.scm -- \
  shell -L analysis/guix -m analysis/guix/manifest.scm --pure -- \
  python3 analysis/tier3_fit.py

guix time-machine -C analysis/guix/channels.scm -- \
  shell -L analysis/guix -m analysis/guix/manifest.scm --pure -- \
  python3 -m pytest -q analysis/tests
```

The first command fetches only the immutable pinned Buffalo TSV and verifies
its SHA-256. To regenerate figures without network access or upstream data:

```sh
guix time-machine -C analysis/guix/channels.scm -- \
  shell -L analysis/guix -m analysis/guix/manifest.scm --pure -- \
  python3 analysis/tier3_fit.py --figure-from-results analysis/tier3_results.tsv
```

The headless PNG/PDF renderer uses only declared Python/NumPy dependencies.
No Conda, micromamba, projected annotation, congener fallback, sparse-VCF
denominator, or unpinned plotting environment is involved.
