# Literature Review and Mathematical Analysis

**Subject.** Evaluating the hypothesis in `main.md` — that *non-crossover
gene conversion, acting as a "homogenization" mechanism, drives an
equilibrium diversity `π = μ/(c·L)` decoupled from `Nₑ` and thereby
explains Lewontin's Paradox of Variation* — against the population-genetic
literature, the coalescent-with-gene-conversion theory, and a Lean 4
formalization (in `lean/`).

---

## 1. What the literature actually establishes

### 1.1 Lewontin's paradox is real and unresolved by linked selection alone

- **Lewontin (1974)** framed it: census `N` varies over ~7–10 orders of
  magnitude across metazoans, while neutral `π` varies over only ~2.
- **Leffler et al. (2012, PLoS Biol 10:e1001388)** ("Revisiting an old
  riddle") is the canonical modern survey; `π` ranges ~2 orders across
  taxa despite vast `N` differences. (Note: `main.md` cites this as
  *Am. J. Hum. Genet.*; it is in fact *PLoS Biology*, PMID 22984349.)
- **Ellegren & Galtier (2016, Nat. Rev. Genet.)** review the determinants.
- **Hermisson & Pfanner (2024, eLife)** show quantitatively that current
  *linked-selection* models (hitchhiking + background selection), even
  parameterized to Drosophila, **cannot** flatten the diversity–`N`
  relationship enough — there is a genuine gap.
- **Achaz & Schertzer (2023, bioRxiv 2023.07.19.549703)** derive "weak
  genetic draft" from sparse selective sweeps, giving `π ∝ N^A` with
  `A < 1`, a mathematically rigorous sub-linear scaling. This is a
  *selection/draft* mechanism, not gene conversion.

So: there is a real explanatory gap, and conversion-homogenization is a
genuinely under-examined candidate. The motivation is legitimate.

### 1.2 Gene conversion: mechanism, frequency, and what it actually does

- **Jeffreys & May (2004, Nat. Genet. 36:151)**; **Cole, Keeney & Jasin
  (2010, Mol. Cell 39:700)** and the review **Cole, Jasin & Keeney
  (2012, DNA Repair 11:617)** ("Preaching about the converted",
  PMC3625938) establish that NCO gene conversions are ~10× more frequent
  than crossovers at hotspots (≈90% of DSB repair products), with short
  tracts (median ~100 bp in mammals; ~300 bp in genome-wide human
  estimates).
- **Williams et al. (2015, eLife)**; **Pålsson et al. (2023, Nat. Genet.
  55:1831)**; **Browning & Browning (2024, PLoS Genet. 20:e1011348)**;
  **Masaki & Browning (2025, PLoS Genet. 21:e1011951)** give the
  genome-wide NCO rate ~6×10⁻⁶ /bp/meiosis, mean tract ~100–300 bp,
  and 82.9% single-conversion tracts.

**The crucial mechanistic fact** (Cole et al. 2012, §"Transmission
distortion and allelic shuffling"): a conversion spanning a heterozygous
site produces a **3:1 / 1:3 tetrad** instead of Mendelian 2:2, and this
produces *transmission distortion* **only when the DSB is asymmetrically
placed** (the "hotspot paradox" / meiotic drive, **Coop & Myers 2007**;
**Jeffreys & Neumann 2002**). With no DSB bias, distortion averages to
**exactly zero**. Equally, **gBGC** (the 68% GC-transmission bias,
**Odenthal-Hesse et al. 2014**) is an *allele-content* bias — formally
equivalent to genic selection — not a uniform "copy IBD → IBD"
homogenization.

### 1.3 The coalescent with gene conversion (the decisive theory)

- **Wiuf & Hein (2000, Genetics 155:451)** (*The Coalescent With Gene
  Conversion*) and **Wiuf & Hein (1997, Genetics 147:1459)** establish
  the formal object: gene conversion is a *recombination-like* event
  (two breakpoints, one per tract end) on the ancestral recombination
  graph. The central structural theorem: **the marginal genealogy at any
  single site is the standard Kingman coalescent, independent of the
  gene-conversion rate `c`.** A conversion straddling the focal site
  merely moves its ancestral material to a different genomic background;
  the *coalescence rate* at the focal site is unchanged.

This is the same fact that makes crossover recombination leave
single-site `E[π] = 4Nₑμ` invariant: recombination/conversion reshuffle
the *correlations between sites* (linkage disequilibrium, haplotype
structure), not the *marginal* diversity.

**Consequence.** At a single neutral site,
`E[π] = 2·μ·E[T_MRCA] = 2·μ·(2Nₑ) = 4·Nₑ·μ`, **independent of `c`**.
Diversity is *not* decoupled from `Nₑ` by gene conversion.

---

## 2. The mathematics, studied

### 2.1 The document's derivation (§5.2), reproduced exactly

The document posits a per-site balance with `H = 2p(1−p) ≈ θ`:

```
Rate_input = μ · H            (mutation "creates heterozygosity")
Rate_loss  = c · L · H²       (conversion "erases heterozygosity")
```

and solves `μ·H = c·L·H²` for the nonzero root `H* = μ/(c·L)`, concluding
`π = μ/(c·L)`, independent of `Nₑ` (formalized verbatim in
`lean/.../DocumentModel.lean`; the algebra is correct — see
`balance_roots`).

### 2.2 Where it breaks: the "loss" term is not a population-genetic drain

The `Rate_loss = c·L·H²` term confuses two different levels:

1. **Per-meiosis (individual) level.** A conversion that copies a
   homozygous donor across a heterozygous acceptor site does make that
   *one gamete* homozygous at that site. This is real, and it is the
   basis of the 3:1/1:3 tetrad distortion.
2. **Population (allele-frequency) level.** That gamete-level change does
   **not** remove the allele from the population — the other homolog, the
   other gametes, and other individuals still carry it. Averaged over
   which homolog is cut (and there is no bias by default), the expected
   change in allele frequency is **zero**. The allele-frequency process
   is a *martingale* with respect to conversion.

Formally (`lean/.../Transmission.lean`), the transmission distortion is

```
δ = P(transmit A) − 1/2  =  −1/4   (DSB on the A-homolog)
                          =  +1/4   (DSB on the a-homolog)
```

so the mean distortion is `E[δ] = (1/2)(−1/4) + (1/2)(+1/4) = 0`
(`expectedDistortion_unbiased`, proven by `norm_num`). Equivalently,
`E[δ] = 0  ⇔  bias = 1/2` (`expectedDistortion_eq_zero_iff`). A *strictly
positive* drain `c·L·H² > 0` (`lossRate_positive`) therefore **cannot** be
the mean diversity-loss rate under unbiased conversion
(`loss_term_contradicts_zero_pressure`). The only directional pressures
are

- **DSB asymmetry** (meiotic drive / hotspot paradox): locus-specific,
  self-limiting (it erodes its own hotspot), not a genome-wide,
  diversity-scaling force;
- **gBGC** (`gBGCdistortion2 b = 2b − 1`, zero iff `b = 1/2`):
  selection-like, sign- and locus-type-dependent, not uniform
  homogenization of IBD tracts.

### 2.3 Does it reduce heterozygosity via a variance / drift effect? (No.)

A natural intuition: even if the *mean* distortion is zero, conversion
adds *variance* to segregation (the `3:1`/`1:3` tetrad vs Mendelian
`2:2`), and `H = 2p(1−p)` is concave, so `E[H]` should drop a little — a
drift-like reduction. This intuition is **incorrect** at the level that
feeds the next generation. A *single* gamete from an `Aa` parent is
`Bernoulli(1/2)` whether or not a conversion spans the site, because
`(1/2)·P(transmit A | DSB on A) + (1/2)·P(transmit A | DSB on a) =
(1/2)·(1/4) + (1/2)·(3/4) = 1/2`. The tetrad-level `3:1` distortion only
appears if all four products of one meiosis are kept together; in
Wright–Fisher the gamete pool abstracts this away, so the intra-meiosis
variance is washed out. The coalescent-with-gene-conversion confirms it
exactly: the marginal genealogy at a single site is Kingman, so
`E[π] = 4·Nₑ·μ` with **no** `c`-dependence — not even a small correction.
Unbiased allelic gene conversion reduces single-site heterozygosity by
**zero**.

A related intuition — "the donor (template) is more likely to be the
common allele, so the rare allele gets overwritten" — is also
mechanistically wrong for *meiotic* conversion: the template is the
**sister homolog** in the same heterozygote, i.e. the *other* allele
(one `A`, one `a`), so it is `A` or `a` with probability `1/2`,
**independent of the population frequency** `p`. That picture *is*
correct for **non-allelic / ectopic** gene conversion (concerted
evolution of multigene families: rDNA arrays, duplicated loci), where
the donor can be a more-numerous paralog — but that is a different
mechanism, not the one behind single-copy neutral diversity.

(Earlier versions of this writeup hedged with a "Jensen second-order"
caveat; that hedge was an error and is withdrawn. The Lean docstring in
`Transmission.lean` now states the correct, exact result.)

### 2.4 The coalescent seals it (the master theorem)

From the marginal-coalescent invariance (Wiuf & Hein 2000),
`E[T_MRCA] = 2Nₑ` regardless of `c`, so

```
E[π] = 2·μ·E[T_MRCA] = 4·Nₑ·μ   (independent of c)
```

(`E_pairwise_pi_eq`, `E_pi_independent_of_c` in `Coalescent.lean`; the
single analytic input `E[T_MRCA] = 2Nₑ` is the imported Kingman
hypothesis; everything else is proven). The document's equilibrium and
the coalescent prediction coincide **only** on the measure-zero
relationship `c·L·Nₑ = 1/4`, which has no population-genetic meaning:

```
μ/(c·L) = 4·Nₑ·μ   ⇔   c·L·Nₑ = 1/4      (for μ ≠ 0, c·L ≠ 0)
```

(`doc_eq_coalescent_iff`, fully proven in `Refutation.lean`). The
assembled conclusion is `master_refutation`:

```
E[δ | unbiased] = 0   ∧   E[π] = 4·Nₑ·μ   ∧
( μ/(c·L) = 4·Nₑ·μ  ⇔  c·L·Nₑ = 1/4 )
```

### 2.5 What is *correct* in the document

The **"silent conversion" / observable-fraction observation is real and
important** — but it is a *measurement* caveat, not a diversity-reducing
force. A tract of length `L` is detectable only if it overlaps a
heterozygous site:

```
P(observable) = 1 − (1 − H)^L  ≤  L·H      (union bound, proven: obsFrac_le_Lh)
```

For `H ≈ 0.001`, `L ≈ 300`: `P(observable) ≈ 0.26`, so ~74% of
conversions leave no molecular footprint (`ObservableFraction.lean`).
This means **estimated conversion rates are underestimates of the true
rate** — a genuine and under-appreciated measurement point. But silent
conversions, by definition, change no allele frequencies, so they
manufacture no homogenization *pressure*: detection probability ≠
diversity loss. The document's correct observation and its (incorrect)
mechanistic conclusion must be separated.

---

## 3. Verdict

The hypothesis as formulated in `main.md` §5.2 is **mathematically
unsound** (see §4 for the one setting where the `μ/(conversion rate)`
shape *is* correct).

- its load-bearing "homogenization pressure" is **zero-mean** under
  unbiased conversion (a martingale, not a drain), proven
  (`Transmission.lean`);
- the marginal coalescent gives `E[π] = 4Nₑμ` **independent of `c`**,
  proven (`Coalescent.lean`);
- the claimed `π = μ/(c·L)` agrees with the neutral result only on the
  spurious special relation `c·L·Nₑ = 1/4`, proven (`Refutation.lean`).

Gene conversion **does not** decouple diversity from `Nₑ`, and the
diversity–`N` compression of Lewontin's paradox is **not** explained by
it in this form.

**What survives, and is worth pursuing:**

1. The **silent-conversion / observable-fraction** point is a legitimate
   measurement caveat: published genome-wide NCO rates (~6×10⁻⁶ /bp) are
   *lower bounds* on the true rate. This deserves direct quantification
   (UK Biobank/TOPMed local-heterozygosity vs. detected-tract-rate
   regression), as the README's "In Progress" already plans.
2. **gBGC** is a real diversity-reducing force, but it acts like
   selection (sign- and locus-type-dependent, `Nₑ`-dependent via
   `4Nₑb`), and its genome-wide contribution to Lewontin's paradox is a
   separate, quantitative question — not "homogenization of IBD tracts."
3. The actually-promising `N`-compressing mechanisms in the literature
   are **weak genetic draft** (Achaz & Schertzer 2023) and the linked-
   selection models reassessed by Hermisson & Pfanner (2024) — both
   selection-based, both genuinely sub-linear in `N`. A rigorous
   comparison of those against the *data* (not against a conversion-
   homogenization strawman) is the productive direction.

The Lean formalization in `lean/` makes the assumptions and the refutation
fully explicit and machine-checked; see `lean/README.md`.

---

## 4. Where `μ/(conversion rate)` *is* correct: concerted evolution of
   repeated gene families

The document's `μ/(c·L)` equilibrium is not invented from nothing — it
is the right *shape* for a real phenomenon, just applied to the wrong
*quantity*. Repeated gene families (rDNA arrays, histone/ubiquitin genes,
 paralogs) *do* homogenize within an individual ("concerted evolution";
Dover 1982, "molecular drive"; Ohta & Dover 1984; Walsh 1985). Under
the *same* unbiased-conversion principle, the difference from the
single-copy case is purely coalescent-structural.

**The load-bearing distinction: are the two copies being compared in
the same individual or in different individuals?**

- *Single-copy, between individuals (Lewontin's paradox):* comparing
  allele `A₁` (individual 1) vs `A₂` (individual 2) at one locus, they
  coalesce at the **species drift rate** `1/(2Nₑ)`; allelic conversion
  (homolog-to-homolog within a meiosis) does not create a faster
  between-individual channel. So `T_MRCA ~ 2Nₑ`, `π = 4·Nₑ·μ`, exact,
  `c`-independent (Wiuf & Hein 2000).

- *Repeat array, within one individual:* comparing copy `i` vs copy `j`
  of a tandem family *in the same individual* — they are at different
  paralogous positions, and an **ectopic** (non-allelic) conversion can
  make copy `i`'s ancestry jump onto copy `j`'s lineage, coalescing them
  at rate `g` (per-pair ectopic-conversion rate), much faster than the
  organism-level coalescent. So `T_MRCA(copy i, copy j | same ind.) ~
  1/g`, and
  ```
  π_within-individual  ≈  2μ / g        (conversion-limited, Nₑ-independent)
  ```
  This *is* the homogenization: the `n` copies in one individual's array
  become nearly identical on a timescale `1/g`, independent of `Nₑ`.

**Why "more copies" matters — and why it is not the donor count.**
Under *uniform* random donor choice a tandem array is still a
martingane (`P(Y→X) = P(X→Y) = f(1−f)·n/(n−1)`), so "the common allele is
more often the donor" is not the mechanism. The mechanism is that the
`n` copies are **co-located and linked** (they co-segregate as an array,
not resegregated by outcrossing each generation), forming a closed
finite "copy population" that (i) is connected by a fast ectopic-
conversion coalescence channel and (ii) drifts to fixation within the
array — copy-level neutral drift, fixation probability = the variant's
initial within-array frequency. The direction is **arbitrary**
(zero-mean); a common variant is likelier to be the winner *by drift*,
not by a conversion-direction bias.

**Across species.** Within-individual homogenization (→ `μ/g`) is only
the first half of concerted evolution. For one variant to take over the
whole species requires organism-level spread — under pure unbiased
conversion that is drift (slow, `Nₑ`-dependent); in practice driven by
**unequal crossing over** (copy-number amplification of the winning
block, a genuine replicator asymmetry), **gBGC**, and/or **selection**.

**Geometric requirement.** A "repeat" present as one copy on each
homolog (two alleles at one locus) converts *allelically* — identical to
the single-copy case, no homogenization. Ectopic homogenization requires
the donor to be a copy at a *different* locus; the two copies compared
must be co-located paralogs in the same genome.

**The document's confusion, precisely:**

| Quantity | Coalescence channel | Result |
|---|---|---|
| within-individual, between-copy diversity of a converting repeat family | ectopic conversion (rate `g`) | `~ μ/g` (homogenized, `Nₑ`-independent) |
| between-individual, single-copy species diversity (Lewontin's paradox) | species drift (rate `1/(2Nₑ)`) | `4·Nₑ·μ` (`c`-independent) |

The `μ/(c·L)` of `main.md` §5.2 is the within-individual repeat-family
result applied to between-individual single-copy diversity. The
homogenization of repeats is real and explained by unbiased ectopic
conversion + copy-level drift (+ unequal crossover/bias/selection for
the across-species step); it does **not** compress between-individual,
single-copy diversity, and so does not explain Lewontin's paradox.
