// GC-biased gene conversion, the recombination requirement, and Lewontin's paradox
// Draft manuscript.  Math is inlined where feasible; the full machine-checked
// argument is in the Lean 4 library at lean/src/LewontinParadox/.

#set page(paper: "us-letter", margin: (x: 1.1in, y: 1.05in))
#set text(font: "New Computer Modern", size: 11pt)
#set par(justify: true, leading: 0.72em, first-line-indent: 1.2em, spacing: 0.8em)
#set heading(numbering: "1.")
#show heading: set text(weight: "bold", size: 12pt)
#show heading: it => v(0.6em) + it + v(0.3em)

#show link: it => text(fill: rgb("#1a5276"), it)

#align(center)[
  #text(size: 16pt, weight: "bold")[
    GC-biased gene conversion, the recombination requirement,\
    and Lewontin's paradox of variation
  ]
  #v(0.5em)
  #text(size: 11pt)[Erik Garrison]
  #v(0.2em)
  #text(size: 9pt, fill: luma(90))[
    Draft. Lean 4 formalization: #link("https://github.com/ekg/lewontin-paradox")[github.com/ekg/lewontin-paradox]
  ]
]
#v(0.6em)

#block(inset: (x: 0.6em, y: 0.4em), stroke: (left: 1.5pt + luma(40)))[
  *Abstract.* Lewontin's paradox — that neutral nucleotide diversity $pi$
  varies by only ~2 orders of magnitude across metazoans whose census
  population sizes span ~7 — is unresolved; current linked-selection models
  do not fully flatten the diversity–size relationship
  (#link("https://doi.org/10.7554/eLife.67509")[Buffalo 2021]); across 172
  metazoan taxa the $pi$–$N_c$ slope is only $~0.11$ (vs. 1 under
  neutrality). We argue that the standard mutation–selection–drift balance at W/S sites
  under GC-biased gene conversion (gBGC) is a candidate contribution to the flattening:
  in large populations it saturates to an $N_e$-independent floor
  $2 mu \/ (u delta)$, where $u$ is the per-site conversion-coverage rate
  and $delta$ the per-conversion GC bias. The bias is always present
  ($delta > 0$, inherent to mismatch repair) and is maintained by selection
  for the sequence homology that meiotic recombination requires. The
  resulting floor is the form previously (and incorrectly) reached via
  "gene-conversion homogenization" arguments; here it is derived by the
  correct mutation–selection–drift route. The argument — including the
  marginal-coalescent refutation of the unbiased case and the
  self-limiting saturation — is formalized and machine-checked in a Lean 4
  library. A new exact-species Tier-3 analysis does not establish the predicted
  saturation or a causal gBGC effect: the across-class GC3 interval includes
  zero and class-specific directions conflict. Biological recovery nevertheless
  yielded three alignment-conditioned H1/H2 coding-panel assembly pairs and two
  population cohorts with explicit callable denominators and genomic
  uncertainty. These assembly, population, and composition modalities are
  reported separately and are insufficient for a cross-species diversity trend.
]

= Introduction

Neutral theory predicts equilibrium diversity $pi = 4 N_e mu$ at a neutral
site, scaling linearly with effective population size $N_e$. Yet across
metazoans, while census sizes $N$ span roughly $10^2$ to $10^(10)$,
within-species $pi$ ranges only ~$10^(-4)$ to $10^(-2)$ — the
"paradox of variation" first noted by #link("https://doi.org/10.1086/419073")[Lewontin (1974)]. The leading
candidate resolution, selection at linked sites (hitchhiking and background
selection), has recently been shown quantitatively insufficient to flatten
the relationship (#link("https://doi.org/10.7554/eLife.67509")[Buffalo 2021]); weak genetic draft
from sparse sweeps gives a sub-linear but not flat scaling
(#link("https://doi.org/10.1101/2023.07.19.549703")[Achaz & Schertzer 2023]).

A recurring but under-developed alternative appeals to non-crossover gene
conversion. The intuition: gene conversion is a *copying* mechanism (it
copies a donor sequence onto a recipient), and copies between
identical-by-descent (IBD) sequences are molecularly silent, so the
observable conversion rate scales with standing heterozygosity — a
"self-regulating homogenization" that, it is claimed, drives an
equilibrium $pi = mu \/ (c dot L)$ decoupled from $N_e$ (with $c$ the
conversion-initiation rate, $L$ the tract length). We show here that this
specific argument is mathematically unsound, but that the destination it
seeks — an $N_e$-independent diversity floor in large populations — is
nonetheless *correct*, reached instead by the standard mutation–selection–
drift balance at W/S sites under gBGC, with the recombination requirement
supplying the evolutionary reason the relevant parameters are maintained.

= Unbiased gene conversion is a martingale: the marginal coalescent

Consider a biallelic `Aa` heterozygote and a gene-conversion tract spanning
the focal site. With no DSB asymmetry, the DSB falls on each homolog with
probability $1\/2$; the donor is the uncut homolog. A *single gamete* from
this meiosis carries `A` with probability

$ "P"("A" | "Aa", "conversion") = 1/2 dot 1/4 + 1/2 dot 3/4 = 1/2, $

exactly Mendelian. The per-gamete variance is therefore $1\/4$ with or
without conversion, the binomial sampling variance $2 N p(1-p)$ is
unchanged, and conversion contributes *zero* drift variance at a single
site. This is the forward-time face of the coalescent result
(#link("https://doi.org/10.1093/genetics/155.1.451")[Wiuf & Hein 2000]): in the coalescent with gene
conversion, the *marginal* genealogy at any single site is the Kingman
coalescent, rate $1\/(2 N_e)$, independent of the conversion rate $c$.
Gene conversion reshuffles the *correlations between sites* (linkage
disequilibrium, the ancestral recombination graph); it does not change the
*marginal* diversity. Hence

$ "E"[pi] = 2 mu dot "E"[T_"MRCA"] = 2 mu dot 2 N_e = 4 N_e mu, quad c-"independent". $

The asserted diversity drain $c dot L dot H^2$ of the "homogenization"
argument is therefore not a valid population-genetic loss term: a
conversion that copies a homozygote across a heterozygous site changes the
*gamete*, not the *population* allele frequency. (Lean:
#link("https://github.com/ekg/lewontin-paradox/blob/main/lean/src/LewontinParadox/Transmission.lean")[`Transmission.lean`],
#link("https://github.com/ekg/lewontin-paradox/blob/main/lean/src/LewontinParadox/Coalescent.lean")[`Coalescent.lean`].)

= GC-biased gene conversion is always on and selection-like

The unbiased case $delta = 0$ (where $delta = b - 1/2$, $b =
"P"("transmit GC" | "W\/S, conversion")$) is a theoretical limit, not the
real regime. The GC bias is inherent to mismatch repair of heteroduplex
mismatches; empirically $b approx 0.68$ ($delta approx 0.18$, often higher
in hotspots). Any $delta > 0$ breaks the martingale and makes gBGC a
directional force favoring GC, with a *discontinuity at zero*: the scaled
strength is zero iff $delta = 0$, and positive otherwise.

Because gBGC only manifests when a conversion tract spans the site, the
effective per-site per-generation selection coefficient is

$ s_"eff" = u dot delta, quad gamma = 4 N_e u delta, $

where $u = c dot L$ is the per-site conversion-coverage rate (the
probability a site lies in a conversion tract per generation). This is
*linear in* $N_e$ — selection-like — so gBGC compresses the large-$N_e$
end of the diversity–size relationship (the right direction for the
paradox), in contrast to classic background selection, which for strong
deleterious selection is $N_e$-*independent* and thus does not compress
the *scaling* (#link("https://doi.org/10.1371/journal.pbio.1002112")[Corbett-Detig et al. 2015]). (Lean:
#link("https://github.com/ekg/lewontin-paradox/blob/main/lean/src/LewontinParadox/GBGC.lean")[`GBGC.lean`].)

A naive estimate with the *genome-wide* $u approx 1.8 times 10^(-3)$ and
$delta approx 0.18$ gives enormous $gamma$ in every species (e.g. $gamma
approx 13$ in humans, $approx 1.3 times 10^3$ in _Drosophila_), which
would crush $pi$ far below observation. That this does not occur is the
self-limitation of gBGC: its rate is proportional to the W/S heterozygosity
$2 f (1-f)$ and vanishes at fixation — gBGC "wins" by depleting AT, at
which point it has no substrate and stops. (Lean:
#link("https://github.com/ekg/lewontin-paradox/blob/main/lean/src/LewontinParadox/Saturation.lean")[`Saturation.lean`].)

= The recombination requirement: why $delta > 0$ and $u$ are maintained
<sec-recomb>

The flattening above requires $delta > 0$ and substantial $u$; these are
not accidents. Meiotic recombination (DSB repair) is essential for
chromosome segregation and fertility. Recent high-resolution work
(#link("https://doi.org/10.1016/j.molcel.2021.08.003")[Ahuja et al. 2021])
shows that both non-crossover (NCO) and crossover (CO) products form by
synthesis-dependent strand annealing, with multiple rounds of strand
invasion and template switching; CO-specific double-Holliday-junction
formation then proceeds via *extensive branch migration* (a median of
$>= 0.66$ kb, in $87%$ of COs). Branch migration through parental duplex
requires near-identity — mismatches impede it and trigger heteroduplex
rejection — so conversion only occurs between near-identical sequences,
and conversion *makes* sequences more identical: a positive feedback
(homology enables conversion, conversion maintains homology). Because
recombination is essential, selection maintains the homology that
conversion maintains, keeping the genome in a "recombino-genic" state
and, with it, the GC-biased repair machinery ($delta > 0$) and high
conversion rates ($u$). This is the evolutionary cause of the parameters
that make the saturation bite, not a separate diversity-reducing force.
(Lean:
#link("https://github.com/ekg/lewontin-paradox/blob/main/lean/src/LewontinParadox/RecombinationHomology.lean")[`RecombinationHomology.lean`].)

= Mutation–selection–drift balance at W/S sites: the saturation

With $delta > 0$ always, gBGC acts as genic selection *against* the AT
allele at every W/S site, with $s = u delta$. In a finite population this
is the textbook mutation–selection–drift balance, with two regimes:

$ "weak gBGC" (2 N_e u delta << 1): quad pi_"W\/S" approx 4 N_e mu quad ("neutral, " N_e "-scaled"); \
 "strong gBGC" (2 N_e u delta >> 1): quad pi_"W\/S" approx 2 mu \/ (u delta) quad (N_e "-independent, saturated"). $

The crossover is at $N_e approx 1 \/ (2 u delta)$. In the strong regime
the AT allele is held at the mutation–selection balance frequency $q^* =
mu \/ s = mu \/ (u delta)$, so the heterozygosity $2 q^* (1 - q^*) approx
2 mu \/ (u delta)$. This floor is (up to a factor of 2 and the bias
$delta$) the form $mu \/ (c dot L)$ sought by the homogenization argument
— but it is reached by mutation–selection–drift balance at W/S sites,
*not* by the flawed "IBD homogenization" drain. The saturation is
$N_e$-independent and is always on (since $delta > 0$); for $N_e$ past the
crossover it pins W/S diversity at the floor, below the neutral
expectation. (Lean:
#link("https://github.com/ekg/lewontin-paradox/blob/main/lean/src/LewontinParadox/MutationSelectionDrift.lean")[`MutationSelectionDrift.lean`],
#link("https://github.com/ekg/lewontin-paradox/blob/main/lean/src/LewontinParadox/Composition.lean")[`Composition.lean`].)

= Magnitude and the effective conversion rate

For the floor to address the paradox, it must lie at the observed
diversity magnitude in large-$N_e$ species. With an effective per-site
coverage $u_"eff" approx 10^(-5)$ (a typical, non-hotspot rate — note the
genome-wide average $approx 1.8 times 10^(-3)$ is dominated by hotspots),
$mu approx 10^(-8)$ and $delta approx 0.18$,

$ 2 mu \/ (u_"eff" delta) approx 2 times 10^(-8) \/ (10^(-5) times 0.18) approx 0.01, $

matching observed $pi$ in large-$N_e$ species. The floor is sensitive to
$u_"eff"$, which is precisely the quantity left uncertain by the
"silent-conversion" measurement problem: most conversion tracts fall in
IBD sequence and leave no molecular footprint, so published genome-wide
rates are lower bounds on the true rate and overstate the average effect.
Quantifying $u_"eff"$ across the genome (from, e.g., UK Biobank
local-heterozygosity versus detected-tract-rate regressions) is the key
empirical input. (Lean:
#link("https://github.com/ekg/lewontin-paradox/blob/main/lean/src/LewontinParadox/ObservableFraction.lean")[`ObservableFraction.lean`],
#link("https://github.com/ekg/lewontin-paradox/blob/main/lean/src/LewontinParadox/Bounds.lean")[`Bounds.lean`].)

= A preliminary test against cross-species data

Using Buffalo's (2021) combined dataset of 172 metazoan taxa (census
sizes $N_c$ estimated from body mass and range; pairwise diversity $pi$
from #link("https://doi.org/10.1371/journal.pbio.1001388")[Leffler et al. 2012],
#link("https://doi.org/10.1371/journal.pbio.1002112")[Corbett-Detig et al. 2015],
#link("https://doi.org/10.1038/nature13685")[Romiguier et al. 2014]), we extend
Buffalo's own regression framework by adding a saturating model and a
recombination-residual sign test that discriminates gBGC from background
selection (BGS).

#figure(
  image("analysis/fig_extend_buffalo.png", width: 100%),
  caption: [Extending Buffalo (2021). _Left:_ $log_(10) pi$ vs $log_(10) N_c$
  for 172 metazoans, with the power law (dashed) and a saturating power law
  (solid). The relationship is a rising-then-plateau: the low-$N_c$ half has
  slope $+0.118$ ($p=0.002$), the high-$N_c$ half slope $+0.018$ ($p=0.58$,
  flat). _Right:_ $pi$ residual (from the power law) against genetic map
  length, both residualized for $N_c$. BGS predicts a positive residual
  (recombination protects diversity); gBGC predicts a negative one (biased
  repair erodes it). The partial correlation is negative ($r=-0.22$), in the
  gBGC direction, though underpowered ($n=39$, $p=0.17$).],
) <fig_buffalo>

*Model comparison (AIC).* On the 172 taxa:

- Power law: $log_(10) pi = -3.24 + 0.110 dot log_(10) N_c$, $R^2 = 0.264$
  (reproducing Buffalo's OLS, $b approx 0.13$).
- Quadratic: weakly concave (negative coefficient), but $Delta "AIC" = +1.0$ —
  the curvature is not justified by the data.
- Pure saturation (Michaelis--Menten, $pi = V_max N_c/(K+N_c)$, i.e. our
  model at $X = 1$): *rejected*, $Delta "AIC" = +179$. It forces a neutral
  slope of $1$ at low $N_c$, but the low-$N_c$ end is far shallower than
  neutral, so it overshoots.
- Saturating power law $pi = A N_c^gamma / (1 + (N_c/K)^gamma)$ (low-$N_c$
  slope $gamma$, high-$N_c$ plateau): $gamma = 0.16$, plateau
  $pi approx 0.026$ above $N_c approx 10^(12)$, $R^2 = 0.272$ — but
  $Delta "AIC" approx 0$, statistically *tied* with the flat power law.

*The structure.* A descriptive split suggests a rising-then-flat pattern
(@fig_buffalo), but the saturating model is statistically tied with a power
law. A high-$N_c$ plateau would be compatible with gBGC saturation if $B =
4 N_e b$ grows with $N_e$; it is not specific evidence for that mechanism.
Likewise the shallow low-$N_c$ rise is compatible with $N_e/N_c$ reduction
(sweepstakes reproduction and demography), but this regression does not
assign causal shares to those processes. The potential division of mechanisms
along $N_c$ is therefore a hypothesis for discriminating observables, not a
result of the Tier-1 fit.

*Under-identification.* The between-species scatter is large ($R^2 = 0.27$),
so the curvature does not beat a straight power law ($Delta "AIC" approx 0$),
and pure saturation is rejected ($Delta "AIC" = +179$). Total-$pi$
cross-species data therefore *cannot* statistically establish gBGC
saturation — it is under-identified. This is precisely why Buffalo (and
prior work) could not see gBGC in total $pi$, and it is why stratification
is essential.

*Recombination sign test.* BGS and gBGC make *opposite* predictions about
recombination: BGS says high recombination $arrow$ *higher* $pi$ (it
protects diversity); gBGC says high recombination $arrow$ *lower* $pi$
(biased repair erodes it). Regressing the power-law residual on genetic map
length, controlling for $N_c$, gives a *negative* partial correlation
($r = -0.22$, $n = 39$, $p = 0.17$); genome size gives $r = -0.23$
($p = 0.17$). The sign is in the gBGC direction, not the BGS direction —
suggestive, though underpowered (Buffalo's recombination column is sparse:
the recombination-rate column has only $n = 18$ joint values).

*Honest read.* Extending Buffalo delivers a qualitative decomposition — gBGC
saturation for the high-$N_c$ plateau, $N_e/N_c$ reduction for the
low-$N_c$ shallow rise — and a recombination-residual sign that points to
gBGC over BGS, both from Buffalo's own data. But total $pi$ is
under-identified: the essential, discriminating test is *W/S-stratified*
$pi$ (the saturation should be sharp at W/S sites and absent at non-W/S
sites, a qualitative prediction not defeated by scatter) together with
*GC-content saturation* (gBGC's hallmark: equilibrium GC plateaus with
$N_e$). The exact-assembly composition compilation and its gated diversity
missingness are reported below. (Code and fits: `analysis/extend_buffalo.py`; data
from #link("https://github.com/vsbuffalo/paradox_variation")[`vsbuffalo/paradox_variation`].)

= Tier 3: exact-species composition and gated diversity

The discriminating plan asks whether GC3 rises and becomes concave with
population size, and whether reference-conditioned $pi_S/pi_W$ or polarized
SFS strength changes in parallel. We joined 135 checksum-validated,
exact-species assemblies to the Buffalo covariate and taxonomy. All 135 have
whole-genome GC; 90 have GC3 from annotation native to the exact accession.
The remaining 45 GC3 cells stay unavailable. A congener, projected annotation,
or whole-genome GC value never fills them.

#figure(
  image("analysis/fig_tier3.png", width: 100%),
  caption: [Tier-3 exact-species composition. Points are exact assemblies,
  colored by Buffalo taxonomic class; lines show the reported class-fixed
  slope through the point centroid. _Left:_ native exact-assembly GC3
  ($n=90$; 45 of 135 composition rows lack eligible native annotation).
  _Right:_ whole-genome GC for all exact FASTAs ($n=135$), retained as a
  weaker sensitivity and never substituted for GC3. Diversity recovery adds
  three alignment-conditioned H1/H2 assembly pairs and two population cohorts;
  they are not plotted against the Buffalo predictor because their modalities
  and sampling frames differ. Polarized SFS-$B$ remains deferred.],
) <fig_tier3>

*Across classes.* In an OLS model with taxonomic-class fixed intercepts, the
GC3 effect per unit of Buffalo `pred_log10_N` was $+0.00892$ ($n=90$;
10,000-species-bootstrap 95% interval $[-0.00123, 0.02149]$; $p=0.163$,
BH $q=0.217$). The interval includes zero. Class intercepts are a coarse
control for deep taxonomy, not phylogenetic generalized least squares: no
frozen branch-length tree or covariance model was supplied. The predictor is
a census-size proxy, not measured $N_e$. The whole-genome GC sensitivity was
positive ($+0.00537$; $n=135$; interval $[0.00194, 0.00976]$; sensitivity
$q=0.0156$), but noncoding composition has many explanations and is not a
gBGC-specific endpoint.

*Within classes and genus.* The GC3 association was positive in Aves
($+0.01459$; $n=11/12$; interval $[0.00781, 0.02447]$; BH $q=0.00217$),
uncertain in Insecta ($+0.02232$; $n=27/48$; interval
$[-0.00301, 0.05728]$; $q=0.217$), negative in Mammalia ($-0.00479$;
$n=25/32$; interval $[-0.01179, 0.00010]$; $q=0.217$), and null within
_Drosophila_ ($+0.00088$; $n=13/17$; interval
$[-0.00634, 0.01317]$; $q=0.839$). Thus these data do not support a uniform
positive within-clade pattern. Leave-one-species and leave-one-class ranges
are retained in the result table; they do not hide the negative mammal or null
_Drosophila_ results.

*Saturation shape.* A centered quadratic class-fixed diagnostic gave a
positive, not negative, quadratic coefficient ($+0.00226$; $n=90$;
bootstrap interval $[0.000005, 0.00510]$; sensitivity $q=0.0602$).
Saturation predicts concavity. A quadratic is not a mechanistic saturation
curve, but its observed convex direction does not support the predicted GC3
saturation shape.

*Recovered assembly modality.* The corrected Tier 3a origin/main SweepGA run
used native `--num-mappings 1:1` mappings followed by fresh IMPG queries. It
produced three biological H1/H2 comparisons over deterministic H1-native coding
panels ($n=2$ eligible haplotypes per pair), not genome-wide deposited-individual
heterozygosity. Coding-gene diversity, defined as H1/H2 alternative alleles per
callable coding base, was $0.013852$ for _Menidia menidia_ (28,233/2,038,234;
1,000-replicate 50-kb block-bootstrap 95% interval
$[0.012354,0.015340]$), $0.00039796$ for _Spinachia spinachia_
(803/2,017,806; $[0.00032951,0.00048123]$), and $0.00167799$ for
_Tautogolabrus adspersus_ (3,551/2,116,227;
$[0.00123835,0.00212181]$). The roughly 35-fold range is heterogeneity, not an
$N_e$ trend: these are three selected assembly pairs without a common
population sampling frame.

The corresponding reference-conditioned $pi_S/pi_W$ point estimates were
$0.6725$, $0.3169$, and $0.7603$. Writing variant counts as $v$ and callable
fourfold-site counts as $d$, they are computed as
$(v_S / d_S) / (v_W / d_W)$. Because paired ratio-bootstrap draws were not published, conservative
intervals formed by dividing the marginal 95% bounds were respectively
$[0.4827,0.9330]$, $[0.0685,1.7020]$, and $[0.2780,2.1252]$; the latter two
include one. This ratio is descriptive and is not polarized SFS-$B$.

*Recovered population modality.* Tier 3b analyzed two _Anopheles coluzzii_
cohorts, each with 20 wild diploid individuals (40 nominal chromosomes), over
the same 21-Mb region but with population-specific callable masks. Population
$pi$ is the sum of unbiased pairwise differences divided by callable sites.
It was $0.0123500$ in AO Luanda (numerator $184,914.842$ over 14,972,821
callable sites; 10,000-replicate chromosome-stratified 1-Mb block-bootstrap 95%
interval $[0.0118634,0.0128173]$) and $0.0158698$ in GM Wali Kunda
($237,893.579/14,990,338$; $[0.0153949,0.0163385]$). The intervals do not
overlap, documenting within-species population heterogeneity. Their
reference-conditioned $pi_S/pi_W$ estimates were $0.74266$
($[0.71292,0.77088]$) and $0.72356$ ($[0.69194,0.75717]$), using the same
10,000-replicate genomic bootstrap. Component $pi_S$ and $pi_W$ uncertainty
uses a 20-individual delete-one jackknife. These population estimates are not
pooled with the assembly comparisons or the single-assembly composition data.
Polarized SFS-$B$ remains deferred because the frozen plan has no outgroup,
ancestral-error model, demographic model, or power threshold.

*Claim boundary.* Measurement denominators identify the precision target
within an assembly; species-bootstrap intervals describe a different,
conditional source of uncertainty. Neither identifies causality. Sparse
Buffalo recombination fields do not define predeclared machinery classes, so
no post-hoc recombination-class fit was invented. Composition alone cannot
separate gBGC from mutation bias, selection, life history, phylogeny, or
assembly/annotation structure. Taken together, the mixed class-specific
directions, convex shape, modality-limited diversity recovery, and deferred SFS
mean that this run *does not establish gBGC saturation or a causal gBGC
explanation of the high-population diversity plateau*. It leaves the formal
mechanism as a testable hypothesis rather than an empirical conclusion.

= Limits and open questions

The saturation flattens *W/S-site* diversity; genome-wide flattening is
partial:

- *W/S fraction.* Only sites that are (or recently were) W/S
  polymorphisms feel gBGC directly. A partial W/S fraction with an
  $N_e$-independent floor and a neutral remainder leaves a mix
  $pi = (1 - X) dot 4 N_e mu + X dot (2 mu \/ (u delta))$, which is flat
  only if $X = 1$ (Lean:
  #link("https://github.com/ekg/lewontin-paradox/blob/main/lean/src/LewontinParadox/GenomeMix.lean")[`GenomeMix.lean`]).
  The linked neutral reduction from gBGC at W/S sites broadens the effect.

- *Effective $u$.* The floor magnitude hinges on $u_"eff"$, uncertain by
  orders of magnitude and concentrated in hotspots.

- *Saturation vs. over-compression.* With the genome-wide $u$, $gamma$ is
  huge; the self-limitation and hotspot concentration are what keep the
  bite finite. The distribution of local $u$ sets the genome-wide average.

- *Comparables.* Whether gBGC's $N_e$-scaled saturation can be
  distinguished from, or combined with, weak genetic draft
  (#link("https://doi.org/10.1101/2023.07.19.549703")[Achaz & Schertzer 2023]) and reassessed linked
  selection (#link("https://doi.org/10.7554/eLife.67509")[Buffalo 2021]) is an open empirical
  question.

= The formalization

The argument is machine-checked in a 16-module Lean 4 library (depending
on Mathlib) at #link("https://github.com/ekg/lewontin-paradox/tree/main/lean/src/LewontinParadox")[`lean/src/LewontinParadox`].
It contains, with no `sorry`, the marginal-coalescent invariance and the
master refutation of the unbiased "homogenization" claim, the gBGC
strength and its $N_e$-linearity, the composition equilibrium, the
self-limiting saturation, the mutation–selection–drift crossover, the
linkage-based distinction between repeat and single-copy diversity, the
genome-mix flatness theorem, and the recombination-homology feedback. The
single imported hypothesis is the Kingman mean $E[T_"MRCA"] = 2 N_e$;
everything else is proven. The build is reproducible:

```sh
cd lean && lake build
```

= References

#text(size: 9pt)[
  - Achaz G, Schertzer E (2023). Weak genetic draft and Lewontin's paradox. _bioRxiv_ 2023.07.19.549703. \
  - Corbett-Detig R et al. (2015). Natural selection constrains neutral diversity across a wide range of species. _PLOS Biol_ 13: e1002112. \
  - Ahuja JS et al. (2021). Repeated strand invasion and extensive branch migration are hallmarks of meiotic recombination. _Mol Cell_ 81: 4258. \
  - Buffalo V (2021). Quantifying the relationship between genetic diversity and population size suggests natural selection cannot explain Lewontin's Paradox. _eLife_ 10: e67509. \
  - Barton HJ, Zeng K (2021). The effective population size modulates the strength of GC biased gene conversion in two passerines. _bioRxiv_ 2021.04.20.440602. \
  - Tisdale et al. (2024). Previously unmeasured genetic diversity explains part of Lewontin's paradox in a k-mer-based meta-analysis of 112 plant species. _Evol Lett_. \
  - Romiguier J et al. (2010). Contrasting GC-content dynamics across 33 mammalian genomes: relationship with life-history traits and chromosome sizes. _Genome Res_ 20: 1001. \
  - Weber CC et al. (2014). Evidence for GC-biased gene conversion as a driver of between-lineage differences in avian base composition. _Genome Biol_ 15: 549. \
  - Galtier N et al. (2018). Codon usage bias in animals: disentangling the effects of natural selection, effective population size, and GC-biased gene conversion. _Mol Biol Evol_ 35: 1092. \
  - Galtier N, Rousselle M (2020). How much does Ne vary among species? _Genetics_ 216: 559. \
  - Lewontin RC (1974). _The Genetic Basis of Evolutionary Change._ Columbia Univ. Press. \
  - Wiuf C, Hein J (2000). The coalescent with gene conversion. _Genetics_ 155: 451–462. \
  - Cole F, Jasin M, Keeney S (2012). Preaching about the converted. _DNA Repair_ 11: 617. \
]
