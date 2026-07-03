/-
! # gBGC is self-limiting (saturation)

gBGC only acts at W/S heterozygous sites; its rate is proportional to the
AT/GC heterozygosity `2·f·(1−f)`. As gBGC depletes AT (`f → 0`), the
substrate shrinks and the rate → 0. This is the irony: gBGC "wins" by
driving everything to GC, at which point it has no substrate left and
stops. So gBGC cannot drive `π → 0` on its own; combined with mutation
(which regenerates AT), it reaches the composition equilibrium
`f* = μ/(u·δ)` (`Composition.lean`), not zero.

This self-limitation is why the naive `exp(−γ)` over-suppression
(`Selection.lean`) does not actually occur: the effective rate is
substrate-limited.
-/

import Mathlib.Data.Real.Basic
import Mathlib.Tactic.Ring
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.NormNum.Basic
import Mathlib.Tactic.Positivity
import LewontinParadox.Heterozygosity
import LewontinParadox.GBGC

namespace LewontinParadox.Saturation

open Real
open LewontinParadox.Heterozygosity (het)
open LewontinParadox.GBGC (selectionCoeff)

/-- gBGC substrate: AT/GC heterozygosity `2·f·(1−f)` (= `Heterozygosity.het`). -/
def substrate (f : ℝ) : ℝ := het f

/-- Effective gBGC removal rate: `u·δ·(substrate)`. Vanishes when there is
no AT/GC heterozygosity. -/
def effectiveRate (u δ f : ℝ) : ℝ := u * δ * substrate f

/-- gBGC's effective rate is zero at fixation: when everything is GC
(`f = 0`) or everything is AT (`f = 1`), there is no substrate. This is
the self-limitation. -/
theorem effectiveRate_zero_at_fixation (u δ f : ℝ) (hf : f = 0 ∨ f = 1) :
    effectiveRate u δ f = 0 := by
  rcases hf with rfl | rfl
  · dsimp [effectiveRate, substrate, het]; norm_num
  · dsimp [effectiveRate, substrate, het]; norm_num

/-- The substrate (= heterozygosity) is maximized at `f = 1/2`. -/
theorem substrate_le_half {f : ℝ} (h0 : 0 ≤ f) (h1 : f ≤ 1) :
    substrate f ≤ 1/2 :=
  LewontinParadox.Heterozygosity.het_le_half h0 h1

/-- Negative feedback: for `0 ≤ f ≤ 1`, the substrate is non-negative, so
the gBGC rate is non-negative; and it → 0 as the AT allele is depleted. -/
theorem substrate_nonneg {f : ℝ} (h0 : 0 ≤ f) (h1 : f ≤ 1) :
    0 ≤ substrate f :=
  LewontinParadox.Heterozygosity.het_nonneg h0 h1

/-- The gBGC removal term `(u·δ)·f·(1−f)` in the composition dynamics
(`Composition.dAT`) is exactly the saturated rate: it vanishes at `f = 0`
and `f = 1`, so gBGC cannot push `f` to either fixation alone — mutation
must balance it, giving `f* = μ/(u·δ) > 0`. -/
theorem removal_term_vanishes_at_fixation (u δ f : ℝ) (hf : f = 0 ∨ f = 1) :
    selectionCoeff u δ * f * (1 - f) = 0 := by
  rcases hf with rfl | rfl
  · dsimp [selectionCoeff]; ring
  · dsimp [selectionCoeff]; ring

end LewontinParadox.Saturation
