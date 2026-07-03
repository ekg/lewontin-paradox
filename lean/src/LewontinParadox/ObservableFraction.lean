/-
! # Observable fraction of gene-conversion events

This module formalizes the part of `main.md` that is **correct**: a
gene-conversion tract is only *detectable* if it overlaps a heterozygous
site. If a tract of length `L` falls entirely in identical-by-descent
(IBD) sequence, it leaves no molecular footprint ("silent" / "invisible"
conversion).

This is a real *measurement* fact, but — crucially — it is a statement
about **detection probability**, not about a **diversity-reducing force**.
Silent conversions are unobservable, but they also do not change allele
frequencies, so they do not homogenize the population. See
`LewontinParadox.Transmission` and `LewontinParadox.Refutation`.
-/

import Mathlib.Data.Real.Basic
import Mathlib.Tactic.Ring
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.NormNum.Basic
import Mathlib.Tactic.Positivity

namespace LewontinParadox.ObservableFraction

open Real

/-- Probability that a conversion tract of length `L` (in bp) contains at
least one heterozygous site, given per-site heterozygosity `h`:
`obsFrac h L = 1 − (1 − h)^L`. -/
def obsFrac (h : ℝ) (L : ℕ) : ℝ := 1 - (1 - h) ^ L

/-- The observable conversion rate is the true rate scaled by the
observable fraction. -/
def observableRate (trueRate h : ℝ) (L : ℕ) : ℝ := trueRate * obsFrac h L

/-- The silent (invisible) fraction is the complement: `(1−h)^L`. -/
def silentFrac (h : ℝ) (L : ℕ) : ℝ := (1 - h) ^ L

theorem obs_plus_silent (h : ℝ) (L : ℕ) :
    obsFrac h L + silentFrac h L = 1 := by
  dsimp [obsFrac, silentFrac]; ring

/-- Helper: `(1−h)^L ≤ 1` for `h ∈ [0,1]`. Renamed to avoid clashing with
Mathlib's own `pow_le_one`. -/
theorem base_pow_le_one {h : ℝ} (h0 : 0 ≤ h) (h1 : h ≤ 1) (L : ℕ) :
    (1 - h) ^ L ≤ 1 := by
  have bn : 0 ≤ 1 - h := by linarith
  have bl1 : 1 - h ≤ 1 := by linarith
  induction L with
  | zero => simp
  | succ n ih =>
    rw [pow_succ]
    have h0n : 0 ≤ (1 - h) ^ n := pow_nonneg bn n
    nlinarith [bn, bl1, h0n, ih]

/-- On the allele-frequency interval `h ∈ [0,1]`, the observable fraction
lies in `[0,1]`. -/
theorem obsFrac_nonneg {h : ℝ} (h0 : 0 ≤ h) (h1 : h ≤ 1) (L : ℕ) :
    0 ≤ obsFrac h L := by
  dsimp [obsFrac]
  exact sub_nonneg.mpr (base_pow_le_one h0 h1 L)

theorem obsFrac_le_one {h : ℝ} (h0 : 0 ≤ h) (h1 : h ≤ 1) (L : ℕ) :
    obsFrac h L ≤ 1 := by
  dsimp [obsFrac]
  have bn : 0 ≤ 1 - h := by linarith
  have := pow_nonneg bn L
  linarith

/-- Edge cases. -/
theorem obsFrac_zero_h (L : ℕ) : obsFrac 0 L = 0 := by
  dsimp [obsFrac]; rw [sub_zero, one_pow]; norm_num

theorem obsFrac_zero_L (h : ℝ) : obsFrac h 0 = 0 := by
  dsimp [obsFrac]; rw [pow_zero]; norm_num

theorem obsFrac_one_L (h : ℝ) : obsFrac h 1 = h := by
  dsimp [obsFrac]; rw [pow_one]; ring

/-- Recurrence in tract length: `obsFrac h (n+1) = h + (1−h)·obsFrac h n`. -/
theorem obsFrac_succ (h : ℝ) (n : ℕ) :
    obsFrac h (n + 1) = h + (1 - h) * obsFrac h n := by
  dsimp [obsFrac]; ring

/-- The observable fraction is monotone in tract length: longer tracts
hit heterozygous sites more often. -/
theorem obsFrac_monotone_L {h : ℝ} (h0 : 0 ≤ h) (h1 : h ≤ 1) (L : ℕ) :
    obsFrac h L ≤ obsFrac h (L + 1) := by
  have bn : 0 ≤ 1 - h := by linarith
  have bl1 : 1 - h ≤ 1 := by linarith
  have pn : 0 ≤ (1 - h) ^ L := pow_nonneg bn L
  have key : (1 - h) ^ (L + 1) ≤ (1 - h) ^ L := by
    rw [pow_succ]
    nlinarith [bn, bl1, pn]
  dsimp [obsFrac]
  linarith

/-- **Union bound (first-order).** `obsFrac h L ≤ L · h`. This is the
qualitative fact the document relies on ("a tract of length `L` hits a
heterozygous site with probability `≈ L·h`"): it is an *upper bound on
detectability*, valid for all `h ∈ [0,1]`, not a rate of diversity loss. -/
theorem obsFrac_le_Lh {h : ℝ} (h0 : 0 ≤ h) (h1 : h ≤ 1) (L : ℕ) :
    obsFrac h L ≤ (L : ℝ) * h := by
  induction L with
  | zero => dsimp [obsFrac]; rw [pow_zero]; norm_num
  | succ n ih =>
    rw [obsFrac_succ h n]
    have bn : 0 ≤ 1 - h := by linarith
    have bl1 : 1 - h ≤ 1 := by linarith
    have ihn : 0 ≤ (n : ℝ) * h := by positivity
    have s1 : (1 - h) * obsFrac h n ≤ (1 - h) * ((n : ℝ) * h) := by nlinarith [bn, ihn]
    have s2 : (1 - h) * ((n : ℝ) * h) ≤ 1 * ((n : ℝ) * h) := by nlinarith [bl1, ihn]
    have s3 : h + 1 * ((n : ℝ) * h) = ((n + 1 : ℕ) : ℝ) * h := by
      rw [Nat.cast_add_one]; ring
    linarith

end LewontinParadox.ObservableFraction
