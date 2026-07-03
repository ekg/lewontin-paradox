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
  (#link("https://doi.org/10.7554/eLife.67509")[Hermisson & Pfanner 2024]).
  We argue that the standard mutation–selection–drift balance at W/S sites
  under GC-biased gene conversion (gBGC) supplies the missing flattening:
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
  library.
]

= Introduction

Neutral theory predicts equilibrium diversity $pi = 4 N_e mu$ at a neutral
site, scaling linearly with effective population size $N_e$. Yet across
metazoans, while census sizes $N$ span roughly $10^2$ to $10^(10)$,
within-species $pi$ ranges only ~$10^(-4)$ to $10^(-2)$ — the
"paradox of variation" first noted by #link("https://doi.org/10.1086/419073")[Lewontin (1974)]. The leading
candidate resolution, selection at linked sites (hitchhiking and background
selection), has recently been shown quantitatively insufficient to flatten
the relationship (#link("https://doi.org/10.7554/eLife.67509")[Hermisson & Pfanner 2024]); weak genetic draft
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

The flattening above requires $delta > 0$ and substantial $u$; these are
not accidents. Meiotic recombination (DSB repair) is essential for
chromromosome segregation and fertility. Non-crossover resolution via
synthesis-dependent strand annealing or a double-Holliday-junction
intermediate requires the invading strand to form a stable heteroduplex
and for the junction to *migrate* away from the DSB; branch migration
requires near-identity, since mismatches impede it. Thus conversion only
occurs between near-identical sequences, and conversion *makes* sequences
more identical — a positive feedback (homology enables conversion,
conversion maintains homology). Because recombination is essential,
selection maintains the homology that conversion maintains, keeping the
genome in a "recombino-genic" state and, with it, the GC-biased repair
machinery ($delta > 0$) and high conversion rates ($u$). This is the
evolutionary cause of the parameters that make the saturation bite, not a
separate diversity-reducing force. (Lean:
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
  selection (#link("https://doi.org/10.7554/eLife.67509")[Hermisson & Pfanner 2024]) is an open empirical
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
  - Hermisson J, Pfanner C (2024). Quantifying the relationship between genetic diversity and population size... _eLife_ 67509. \
  - Lewontin RC (1974). _The Genetic Basis of Evolutionary Change._ Columbia Univ. Press. \
  - Wiuf C, Hein J (2000). The coalescent with gene conversion. _Genetics_ 155: 451–462. \
  - Cole F, Jasin M, Keeney S (2012). Preaching about the converted. _DNA Repair_ 11: 617. \
]
