# Lean 4 formalization

A small, machine-checked companion to `../review.md` and `../main.md`.
It depends on **Mathlib** (only for the real numbers and the `ring` /
`linarith` / `field_simp` / `positivity` tactics); the population-genetic
content is all ours.

## What is proven (no `sorry`, no axioms besides one clearly-labeled
imported hypothesis)

| Module | Statement (selected) | What it says |
|---|---|---|
| `Heterozygosity` | `het p = 2·p·(1−p)`, `het_nonneg`, `het_le_half`, `het_eq_zero_iff`, `theta Nₑ μ = 4·Nₑ·μ` | the elementary heterozygosity / θ facts |
| `ObservableFraction` | `obsFrac h L = 1 − (1−h)^L`, `obsFrac_le_Lh` (`≤ L·h`), `obsFrac_nonneg`, `obs_plus_silent` | the **correct** part of the document: silent-conversion *detection* probability (a measurement caveat, not a diversity force) |
| `Transmission` | `distortion4 (1/2) = 0`, `expectedDistortion_eq_zero_iff` (`= 0 ↔ bias = 1/2`), `gBGC_distortion2_eq_zero_iff` | **unbiased gene conversion has zero mean transmission distortion** — a martingale, not a homogenization drain; gBGC is the only directional bias |
| `DocumentModel` | `balance_roots`, `lossRate_positive` | the document's §5.2 derivation, reproduced *verbatim*: the algebra is correct |
| `Coalescent` | `E_pairwise_pi_eq` (`= 4·Nₑ·μ`), `E_pi_independent_of_c` | from the marginal-coalescent invariance (Wiuf & Hein 2000), `E[π]` is independent of `c` |
| `Refutation` | `doc_eq_coalescent_iff` (`μ/(c·L) = 4·Nₑ·μ ↔ c·L·Nₑ = 1/4`), `loss_term_contradicts_zero_pressure`, `master_refutation` | the document's equilibrium matches the neutral result only on a measure-zero relation; the asserted loss term contradicts the zero mean pressure |

**The one imported hypothesis** (`Coalescent.E_Tmrca_eq`): the Kingman-
coalescent mean pairwise time `E[T_MRCA] = 2Nₑ`. Proving it would require
formalizing the full coalescent process; the *consequences* (including
`c`-independence) are fully proven.

## Build

```bash
cd lean
lake build          # fetches+builds the Mathlib slice once (~15 min first time)
```

The default target is the `LewontinParadox` library; `lake build` builds it
and prints only linter notes (some intentionally meaningful, e.g. the `c`
in `E_pi_independent_of_c` is unused *because* c-independence is the
theorem).

## Layout

```
src/LewontinParadox/
  Heterozygosity.lean     elementary H = 2p(1−p), θ = 4Nₑμ
  ObservableFraction.lean the legitimate "silent conversion" detection fact
  Transmission.lean       transmission-distortion model; the refutation of a
                          directional "homogenization pressure"
  DocumentModel.lean      main.md §5.2, formalized verbatim
  Coalescent.lean         marginal coalescent → E[π] = 4Nₑμ, c-independent
  Refutation.lean         assembles the master theorem
  LewontinParadox.lean    root import
```

Read the module docstrings (the `/-- ... -/` blocks) — they contain the
biological narrative and the citations (Wiuf & Hein 2000; Cole et al.
2012; Achaz & Schertzer 2023; Hermisson & Pfanner 2024).
