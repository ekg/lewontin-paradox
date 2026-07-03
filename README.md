# Lewontin's Paradox, gBGC, and the Recombination Requirement

A research project investigating whether **GC-biased gene conversion (gBGC),
maintained by the recombination requirement, flattens the diversity–population-size
relationship** and thereby resolves **Lewontin's paradox of variation**.

> Across metazoans, census population sizes span ~7 orders of magnitude, yet
> neutral nucleotide diversity π varies by only ~2. The leading candidate
> (linked selection) does not fully flatten the relationship
> ([Hermisson & Pfanner 2024](https://doi.org/10.7554/eLife.67509)). We argue
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

Draft. Open empirical questions: the effective per-site conversion rate `u_eff`
(the "silent-conversion" measurement problem), the W/S fraction under strong
gBGC, and the linked neutral reduction from gBGC.

## Author

Erik Garrison — erik.garrison@gmail.com
