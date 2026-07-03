/-
! # Quantitative bounds: how big must `u·δ` be?

From `Selection.lean`, a selection-like force with scaled strength `γ`
reduces diversity by a factor `exp(−γ)`; to compress by a factor `F ≥ 1`
requires `γ ≥ ln(F)`. For gBGC, `γ = 4·Nₑ·u·δ`, so the requirement is

    u·δ ≥ ln(F) / (4·Nₑ).

Order-of-magnitude implication (in prose, since `ln` is transcendental):
with the *genome-wide* conversion-coverage rate `u ≈ c·L ≈ 1.8×10⁻³` and
the *observed* bias `δ ≈ 0.18`, `γ = 4·Nₑ·u·δ` is enormous in every
species — gBGC would crush `π` far below what is observed. That it does
*not* is precisely the saturation point (`Saturation.lean`): gBGC's
effective rate is substrate-limited and hotspot-concentrated, so the
*effective* `u·δ` genome-wide is orders of magnitude below the genome-wide
average. So gBGC is a real, `Nₑ`-scaled, selection-like force, but it is
self-limiting and its genome-wide bite is bounded well below the naive
`exp(−γ)` estimate.
-/

import Mathlib.Data.Real.Basic
import Mathlib.Tactic.Ring
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.NormNum.Basic
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Positivity
import Mathlib.Analysis.SpecialFunctions.Exp
import Mathlib.Analysis.SpecialFunctions.Log.Basic
import LewontinParadox.Selection
import LewontinParadox.GBGC

namespace LewontinParadox.Bounds

open Real
open LewontinParadox.Selection (scaledSelection reductionFactor required_per_site_selection)
open LewontinParadox.GBGC (strength selectionCoeff)

/-- To compress diversity by a factor `F` at population size `Nₑ` via
gBGC requires `u·δ ≥ ln(F)/(4·Nₑ)`. -/
theorem required_product (Nₑ F u δ : ℝ) (hN : 0 < Nₑ) (hF : 1 ≤ F) :
    reductionFactor (strength Nₑ u δ) ≤ 1 / F ↔ 4 * Nₑ * (u * δ) ≥ Real.log F := by
  have key : strength Nₑ u δ = scaledSelection Nₑ (u * δ) := by
    dsimp [strength, scaledSelection]; ring
  rw [key]
  exact required_per_site_selection Nₑ (u * δ) F hN hF

/-- For a fixed bias `δ` and coverage `u`, gBGC strength grows linearly
with `Nₑ`: large populations face disproportionately strong gBGC. -/
theorem strength_grows_with_Nₑ (u δ Nₑ₁ Nₑ₂ : ℝ) (hN : 0 < Nₑ₂)
    (hδ : 0 < δ) (hu : 0 < u) :
    strength Nₑ₂ u δ ≥ strength Nₑ₁ u δ ↔ Nₑ₂ ≥ Nₑ₁ := by
  dsimp [strength]
  have h4uδ : 0 < 4 * u * δ := by nlinarith
  constructor
  · intro h; nlinarith [h4uδ]
  · intro h; nlinarith [h4uδ]

end LewontinParadox.Bounds
