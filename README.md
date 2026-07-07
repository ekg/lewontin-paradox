# Lewontin's Paradox, gBGC, and the Recombination Requirement

A research project investigating whether **GC-biased gene conversion (gBGC),
maintained by the recombination requirement, flattens the diversity–population-size
relationship** and thereby resolves **Lewontin's paradox of variation**.

> Across metazoans, census population sizes span ~7 orders of magnitude, yet
> neutral nucleotide diversity π varies by only ~2. The leading candidate
> (linked selection) does not fully flatten the relationship
> ([Buffalo 2021](https://doi.org/10.7554/eLife.67509)). We argue
> that the standard **mutation–selection–drift balance at W/S sites under gBGC**
> saturates to an Nₑ-independent floor `2μ/(uδ)` in large populations — the form
> previously (and incorrectly) sought via "gene-conversion homogenization,"
> here derived by the correct route, with the recombination requirement
> explaining why the bias (δ>0) and conversion rate (u) are maintained.

## Contents

| Path | What |
|---|---|
| `manuscript.typ` / `manuscript.pdf` | the draft manuscript (Typst) with the math |
| `review.md` | full literature review + mathematical analysis |
| `lean/` | a 16-module **Lean 4** formalization (depends on Mathlib), machine-checked, no `sorry` |
| `main.md` / `main.pdf` | the original research document that motivated the project |

## The argument, in one figure

```
unbiased conversion  →  martingale, marginal Kingman  →  π = 4Nₑμ  (c-independent)
        add δ > 0 (gBGC, always on)  →  selection-like, γ = 4Nₑ·u·δ
   + finite-N drift  →  mutation–selection–drift balance at W/S sites
        strong gBGC (large Nₑ)  →  π_W/S saturates to  2μ/(u·δ)   (Nₑ-independent)
   recombination requires homology  →  selection maintains δ>0 and u  (always on)
```

The destination `μ/(c·L)` (≈ `2μ/(uδ)`, Nₑ-independent) is real — but reached
by **mutation–selection–drift balance at W/S sites**, not by the flawed
"IBD-homogenization" drain. The floor magnitude matches observed π in large-Nₑ
species for a plausible effective conversion rate `u_eff ~ 10⁻⁵`.

## The Lean formalization

16 modules, `lake build` clean, one labeled imported hypothesis
(`E[T_MRCA] = 2Nₑ`, the Kingman mean); everything else is proven. Modules:

- `Heterozygosity` — `H(p)=2p(1−p)`, `θ=4Nₑμ`
- `ObservableFraction` — the *correct* "silent conversion" detection fact (`obsFrac ≤ L·h`)
- `Transmission` — unbiased conversion is a martingale; gBGC directional
- `DocumentModel` — the original §5.2, formalized verbatim (`H*=μ/(c·L)`)
- `Coalescent` — marginal Kingman: `E[π]=4Nₑμ`, `c`-independent
- `Refutation` — `master_refutation`, `doc_eq_coalescent_iff`
- `Selection` — generic `γ=4Nₑs`, `R=exp(−γ)`, BGS Nₑ-independent
- `GBGC` — `γ=4Nₑuδ`, linear in Nₑ, `=0 ⟺ δ=0`
- `Composition` — `f*_AT=μ/(uδ)`: the real home of `μ/(c·L)` (GC composition)
- `Saturation` — gBGC self-limiting: rate `∝f(1−f)`, vanishes at fixation
- `Bounds` — `4Nₑuδ ≥ ln F` to compress diversity by F
- `Repeats` — linkage: `π_within≈2μ/g` vs single-copy `4Nₑμ`
- `GenomeMix` — partial repeat fraction: flat ⟺ X=1 or μ=0
- `Drift` — finite-N: per-gamete var `1/4` both; tetrad var `0` vs `1`
- `MutationSelectionDrift` — the saturation: weak `4Nₑμ` → strong `2μ/(uδ)`, crossover `1/(2uδ)`
- `RecombinationHomology` — the feedback: conversion needs homology; selection maintains it

### Build

```sh
cd lean && lake build      # fetches+builds the Mathlib slice once (~15 min first time)
```

## Status

Draft, with a first empirical pass against Buffalo's (2021) 172-taxon
metazoan dataset (`analysis/extend_buffalo.py`). Findings:

- **Tier 1 (n=172):** the curve is a *rising-then-plateau* — low-Nc half
  slope +0.118 (p=0.002), high-Nc half +0.018 (p=0.58, flat). AIC model
  comparison: power law (slope 0.110, R²=0.264) ties a saturating power law
  (ΔAIC≈0); pure Michaelis-Menten (X=1) is rejected (ΔAIC=+179). Total π is
  *under-identified* for gBGC.
- **Tier 2 (n=18–39):** recombination-residual sign test — BGS predicts +,
  gBGC predicts −. Map-length partial r=−0.22 (ctrl Nc), genome-size
  r=−0.23: the *gBGC direction*, not BGS (underpowered).
- **Tier 3 (the discriminating test, literature):** gBGC's two cross-species
  signatures. *Within clade*, both hold — composition: GC3 rises with Ne
  (mammals, Romiguier 2010: ρ=−0.44 vs body mass; birds, Weber 2014);
  strength: B=4Neb ∝ Ne (passerines, Barton & Zeng 2021, ~2× with 2× Ne).
  *Across clades*, the strength signature vanishes (Galtier 2018, 30 metazoans:
  B ⊥ Ne, because r·l·b0 varies). The §4 homology-maintenance mechanism
  *predicts exactly this* — within-clade machinery conserved → B∝Ne →
  saturation at high Ne; across-clade drift of the recombination landscape →
  attenuation, with selection re-maintaining δ>0 and high u (contra Galtier's
  "inverse b0"). Testable prediction: the cross-species B–Ne correlation
  should reappear within recombination-machinery classes.

**Bottom line:** gBGC saturation is a real, Ne-scaled, within-clade force that
flattens diversity at high Ne — a big and overlooked part of Lewontin's
Paradox, not the whole story. The low-Nc shallow rise is Ne/Nc reduction
(Buffalo's territory); the across-clade attenuation of the strength
signature is the drift of the recombination machinery that selection
repeatedly corrects.

Open empirical questions: the effective per-site conversion rate `u_eff`
(the "silent-conversion" measurement problem), the W/S fraction under strong
gBGC, the linked neutral reduction from gBGC, and the within-recombination-class
B–Ne correlation test.

## Author

Erik Garrison — erik.garrison@gmail.com
