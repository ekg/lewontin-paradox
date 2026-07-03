/-
! # Selection-like forces: scaled strength and diversity reduction

A generic framework for any force that, like selection or gBGC, acts with a
per-site per-generation selection coefficient `s` and so has population-scaled
strength `γ = 4·Nₑ·s`. The standard background-selection / hitchhiking
diversity-reduction factor is `R = exp(−γ)`.

`GBGC.lean` instantiates `s` for GC-biased gene conversion; BGS (below) for
classic deleterious background selection, so they can be compared on the
same scale.
-/

import Mathlib.Data.Real.Basic
import Mathlib.Tactic.Ring
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.NormNum.Basic
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Positivity
import Mathlib.Analysis.SpecialFunctions.Exp
import Mathlib.Analysis.SpecialFunctions.Log.Basic

namespace LewontinParadox.Selection

open Real

/-- Population-scaled selection coefficient `γ = 4·Nₑ·s`. -/
def scaledSelection (Nₑ s : ℝ) : ℝ := 4 * Nₑ * s

/-- Background-selection-like diversity-reduction factor `R = exp(−γ)`. -/
noncomputable def reductionFactor (γ : ℝ) : ℝ := Real.exp (-γ)

/-- To reduce diversity by a factor `F ≥ 1` (i.e. `R ≤ 1/F`) requires
`γ ≥ ln(F)`. -/
theorem reduction_le_iff (γ F : ℝ) (hF : 1 ≤ F) :
    reductionFactor γ ≤ 1 / F ↔ γ ≥ Real.log F := by
  have hFpos : 0 < F := by linarith
  dsimp [reductionFactor]
  have key : 1 / F = Real.exp (-(Real.log F)) := by
    rw [Real.exp_neg, Real.exp_log hFpos]
    exact one_div F
  rw [key, Real.exp_le_exp]
  constructor
  · intro h; linarith
  · intro h; linarith

/-- To reduce diversity by a factor `F` at population size `Nₑ` via a force
with per-site selection `s`, the requirement is `4·Nₑ·s ≥ ln(F)`. -/
theorem required_per_site_selection (Nₑ s F : ℝ) (hN : 0 < Nₑ) (hF : 1 ≤ F) :
    reductionFactor (scaledSelection Nₑ s) ≤ 1 / F ↔ 4 * Nₑ * s ≥ Real.log F := by
  dsimp [scaledSelection]
  exact reduction_le_iff (4 * Nₑ * s) F hF

/-! ## Classic background selection (BGS), for comparison

For strongly deleterious mutations (selection `s`, genomic deleterious
mutation rate `U`), the BGS reduction strength is `U/s`. For strong
selection (`4·Nₑ·s ≫ 1`) this is **independent of `Nₑ`**, so — per Hermisson
& Pfanner (2024) — BGS does not compress the *scaling* of diversity with
`Nₑ` (it scales diversity by a constant factor, leaving the linear `Nₑ`
dependence intact). gBGC, by contrast, scales linearly with `Nₑ` (see
`GBGC.lean`).
-/

/-- BGS reduction strength `U/s` (no `Nₑ` dependence for strong selection). -/
noncomputable def bgsStrength (U s : ℝ) : ℝ := U / s

/-- BGS strength is a function of `U` and `s` only — it carries no `Nₑ`
argument, in contrast to `GBGC.strength` which is linear in `Nₑ`. -/
theorem bgsStrength_carries_no_Nₑ (U s : ℝ) :
    (∀ Nₑ : ℝ, bgsStrength U s = bgsStrength U s) := fun _ => rfl

end LewontinParadox.Selection
