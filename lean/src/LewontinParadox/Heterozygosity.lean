/-
! # Heterozygosity

Per-site heterozygosity at a biallelic locus with derived-allele
frequency `p` is `H(p) = 2·p·(1−p)`. Under Hardy–Weinberg this is the
expected fraction of diploid individuals that are heterozygous; under
the neutral infinite-sites model it equals the population-scaled
mutation rate `θ = 4·Nₑ·μ`.
-/

import Mathlib.Data.Real.Basic
import Mathlib.Tactic.Ring
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.NormNum.Basic
import Mathlib.Tactic.Positivity

namespace LewontinParadox.Heterozygosity

open Real

/-- Per-site heterozygosity `H(p) = 2·p·(1−p)`. -/
def het (p : ℝ) : ℝ := 2 * p * (1 - p)

/-- `H(p)` is symmetric under allele relabelling: `H(p) = H(1−p)`. -/
theorem het_symm (p : ℝ) : het p = het (1 - p) := by
  dsimp [het]; ring

/-- `H(p) ≥ 0` on the allele-frequency interval `[0,1]`. -/
theorem het_nonneg {p : ℝ} (h0 : 0 ≤ p) (h1 : p ≤ 1) : 0 ≤ het p := by
  dsimp [het]
  have h2p : 0 ≤ 2 * p := mul_nonneg (by norm_num) h0
  have h1mp : 0 ≤ 1 - p := by linarith
  exact mul_nonneg h2p h1mp

/-- `H(p) ≤ 1/2` on `[0,1]`; the maximum is attained at `p = 1/2`. -/
theorem het_le_half {p : ℝ} (h0 : 0 ≤ p) (h1 : p ≤ 1) : het p ≤ 1/2 := by
  have key : het p = 1/2 - 2 * (p - 1/2) ^ 2 := by dsimp [het]; ring
  rw [key]
  have h2 : 0 ≤ 2 * (p - 1/2) ^ 2 :=
    mul_nonneg (by norm_num) (sq_nonneg _)
  linarith

/-- `H(p) = 0  ⇔  p = 0 ∨ p = 1` (no heterozygosity at fixation). -/
theorem het_eq_zero_iff (p : ℝ) : het p = 0 ↔ p = 0 ∨ p = 1 := by
  constructor
  · intro h
    have h2 : (2:ℝ) ≠ 0 := by norm_num
    have hp : p * (1 - p) = 0 := mul_right_cancel₀ h2 (by dsimp [het] at h; linarith)
    rcases mul_eq_zero.mp hp with h1 | h2
    · left; linarith
    · right; linarith
  · rintro (rfl | rfl)
    · dsimp [het]; norm_num
    · dsimp [het]; norm_num

/-- Heterozygosity is the population-scaled mutation rate `θ = 4·Nₑ·μ`. -/
def theta (Nₑ μ : ℝ) : ℝ := 4 * Nₑ * μ

end LewontinParadox.Heterozygosity
