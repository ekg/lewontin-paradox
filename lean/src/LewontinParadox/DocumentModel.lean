/-
! # The document's model, formalized exactly as written

This module formalizes the derivation in `main.md` §5.2 verbatim, so that
the assumptions are explicit and the conclusion can be checked against the
coalescent (see `LewontinParadox.Coalescent`,
`LewontinParadox.Refutation`).

The document asserts a per-site balance
```
    Rate_input = μ · H         (mutation creates heterozygosity)
    Rate_loss  = c · L · H²    (conversion erases heterozygosity)
```
and solves `μ · H = c · L · H²` for the equilibrium `H* = μ/(c·L)`,
claiming diversity is decoupled from `Nₑ`.

The algebra here is correct; the problem is the **loss term** (see
`LewontinParadox.Transmission`): a conversion that copies a homozygote
across a heterozygous site changes the *gamete*, not the *population*,
and under symmetric DSB the mean allele-frequency change is zero. So
`Rate_loss = c·L·H²` is not a population-genetic diversity drain.
-/

import Mathlib.Data.Real.Basic
import Mathlib.Tactic.Ring
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.NormNum.Basic
import Mathlib.Tactic.FieldSimp
import LewontinParadox.Heterozygosity

namespace LewontinParadox.DocumentModel

open Real LewontinParadox.Heterozygosity

/-- The "diversity-loss" rate asserted in `main.md` §5.2:
`Rate_loss = c · L · H²`. -/
def lossRate (c L h : ℝ) : ℝ := c * L * h ^ 2

/-- The "diversity-input" rate from mutation asserted in §5.2:
`Rate_input = μ · H`. -/
def inputRate (μ h : ℝ) : ℝ := μ * h

/-- The balance equation `μ · H = c · L · H²` from §5.2. -/
def balance (μ c L h : ℝ) : Prop := inputRate μ h = lossRate c L h

/-- The document's headline equilibrium `π = μ / (c · L)` (§5.2). -/
noncomputable def equilibriumPi (μ c L : ℝ) : ℝ := μ / (c * L)

/-- The asserted balance equation has roots `H = 0` and `H = μ/(c·L)`
(for `c·L ≠ 0`). The algebra of the document is correct. -/
theorem balance_roots (μ c L h : ℝ) (hCL : c * L ≠ 0) :
    balance μ c L h ↔ h = 0 ∨ h = equilibriumPi μ c L := by
  dsimp [balance, inputRate, lossRate, equilibriumPi]
  constructor
  · intro heq
    have key : μ * h - c * L * h ^ 2 = h * (μ - c * L * h) := by ring
    have z : μ * h - c * L * h ^ 2 = 0 := by linarith
    have factored : h * (μ - c * L * h) = 0 := by linarith
    rcases mul_eq_zero.mp factored with h0 | h1
    · left; linarith
    · right
      have h2 : (c * L) * h = μ := by linarith
      have rhs : (c * L) * (μ / (c * L)) = μ := mul_div_cancel₀ μ hCL
      have eq2 : (c * L) * h = (c * L) * (μ / (c * L)) := by linarith
      exact mul_left_cancel₀ hCL eq2
  · rintro (rfl | hr)
    · ring
    · rw [hr]; field_simp; try ring

/-- For strictly positive `c, L, h`, the asserted loss term is strictly
positive — i.e. the document posits a strictly positive diversity drain
whenever there is standing diversity. -/
theorem lossRate_positive {c L h : ℝ} (hc : 0 < c) (hL : 0 < L) (hh : 0 < h) :
    0 < lossRate c L h := by
  dsimp [lossRate]
  exact mul_pos (mul_pos hc hL) (pow_pos hh 2)

end LewontinParadox.DocumentModel
