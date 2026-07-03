/-
! # Refutation: the document's equilibrium is inconsistent with the
! # standard neutral model

This module assembles the pieces and states the central conclusion
formally. The document (`main.md` §5.2, §7) claims that gene-conversion
homogenization drives diversity to a `c`-dependent equilibrium
`π = μ/(c·L)` that is **decoupled from `Nₑ`**, and that this explains
Lewontin's paradox.

Two independently-established facts make this untenable:

1. **No directional pressure.** Under symmetric DSB (the default,
   unbiased case), gene conversion exerts zero mean transmission
   distortion (`LewontinParadox.Transmission.expectedDistortion_unbiased`).
   It cannot drive alleles to fixation, so it cannot be a directional
   "homogenization pressure". The only directional effects are gBGC and
   DSB-asymmetry (meiotic drive), both bias-specific and locus-specific,
   not a genome-wide, diversity-scaling drain.

2. **Marginal coalescent.** Gene conversion does not change the marginal
   genealogy at a single site, so `E[π] = 4·Nₑ·μ` independent of `c`
   (`LewontinParadox.Coalescent`).

Therefore the document's `π = μ/(c·L)` and the coalescent's
`π = 4·Nₑ·μ` coincide only on the special (measure-zero) parameter
relationship `c·L·Nₑ = 1/4`, which has no population-genetic basis. The
"silent conversion" observation remains a valid *measurement* caveat
(`LewontinParadox.ObservableFraction`), but it does not manufacture a
diversity-reducing force.
-/

import Mathlib.Data.Real.Basic
import Mathlib.Tactic.Ring
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.NormNum.Basic
import Mathlib.Tactic.FieldSimp
import LewontinParadox.Heterozygosity
import LewontinParadox.ObservableFraction
import LewontinParadox.Transmission
import LewontinParadox.DocumentModel
import LewontinParadox.Coalescent

namespace LewontinParadox.Refutation

open Real
open LewontinParadox.Heterozygosity (theta)
open LewontinParadox.Transmission (distortion4 gBGCdistortion2)
open LewontinParadox.DocumentModel (lossRate equilibriumPi)
open LewontinParadox.Coalescent (E_pairwise_pi E_pairwise_pi_eq)

/-! ## (1) Unbiased conversion is not a directional pressure -/

/-- Unbiased gene conversion exerts zero mean transmission distortion. -/
theorem unbiased_conversion_no_pressure :
    distortion4 (1/2) = 0 :=
  LewontinParadox.Transmission.expectedDistortion_unbiased

/-- Unbiased gBGC (`b = 1/2`) also exerts zero mean pressure. -/
theorem unbiased_gBGC_no_pressure :
    gBGCdistortion2 (1/2) = 0 :=
  LewontinParadox.Transmission.gBGC_distortion2_unbiased

/-! ## (2) Diversity tracks `Nₑ`, not `c` -/

/-- Expected pairwise diversity equals `θ = 4·Nₑ·μ`, the standard neutral
result, and is independent of the gene-conversion rate. -/
theorem diversity_tracks_Nₑ (Nₑ μ : ℝ) :
    E_pairwise_pi Nₑ μ = theta Nₑ μ :=
  LewontinParadox.Coalescent.E_pi_independent_of_c Nₑ μ 0

/-! ## (3) The document's equilibrium does not generically equal the
   coalescent prediction -/

/-- The document's `μ/(c·L)` equals the coalescent's `4·Nₑ·μ` **only** on
the special relationship `c·L·Nₑ = 1/4` (for `μ ≠ 0`, `c·L ≠ 0`). -/
theorem doc_eq_coalescent_iff (Nₑ μ c L : ℝ) (hμ : μ ≠ 0) (hCL : c * L ≠ 0) :
    equilibriumPi μ c L = 4 * Nₑ * μ ↔ c * L * Nₑ = 1/4 := by
  dsimp [equilibriumPi]
  constructor
  · intro h
    -- h : μ / (c * L) = 4 * Nₑ * μ. Rewrite as μ = 4·Nₑ·μ·(cL).
    rw [div_eq_iff hCL] at h
    have assoc : 4 * Nₑ * μ * (c * L) = (4 * Nₑ * (c * L)) * μ := by ring
    have h2 : 1 * μ = (4 * Nₑ * (c * L)) * μ := by linarith [h, assoc]
    have hB : 4 * Nₑ * (c * L) = 1 := (mul_right_cancel₀ hμ h2).symm
    have hC : 4 * (c * L * Nₑ) = 1 := by
      rw [show 4 * (c * L * Nₑ) = 4 * Nₑ * (c * L) from by ring]; exact hB
    linarith
  · intro h
    -- h : c * L * Nₑ = 1/4
    have hC : 4 * (c * L * Nₑ) = 1 := by linarith
    have hB : 4 * Nₑ * (c * L) = 1 := by
      rw [show 4 * Nₑ * (c * L) = 4 * (c * L * Nₑ) from by ring]; exact hC
    rw [div_eq_iff hCL]
    -- goal : μ = 4 * Nₑ * μ * (c * L)
    have key : 4 * Nₑ * μ * (c * L) = μ := by
      have assoc : (4 * Nₑ * (c * L)) * μ = 4 * Nₑ * μ * (c * L) := by ring
      rw [← assoc, hB]; ring
    linarith [key]

/-! ## (4) The asserted loss term is a strictly positive drain, contradicting
   the zero mean pressure -/

/-- For strictly positive `c, L, H`, the document's loss term
`c·L·H²` is strictly positive, whereas the true mean allele-frequency
pressure from unbiased conversion is zero. Hence the document's loss term
is not a valid population-genetic diversity drain. -/
theorem loss_term_contradicts_zero_pressure
    {c L h : ℝ} (hc : 0 < c) (hL : 0 < L) (hh : 0 < h) :
    lossRate c L h ≠ 0 :=
  ne_of_gt (LewontinParadox.DocumentModel.lossRate_positive hc hL hh)

/-- **Master statement.** Unbiased gene conversion exerts no directional
allele-frequency pressure (`expectedDistortion (1/2) = 0`); the marginal
coalescent gives `E[π] = 4·Nₑ·μ` independent of `c`; and the document's
`μ/(c·L)` agrees with this only on the special relationship
`c·L·Nₑ = 1/4`. Thus gene-conversion "homogenization" does not decouple
diversity from population size and does not, as formulated, explain
Lewontin's paradox. -/
theorem master_refutation
    (Nₑ μ c L : ℝ) (hμ : μ ≠ 0) (hCL : c * L ≠ 0) :
    distortion4 (1/2) = 0 ∧
    E_pairwise_pi Nₑ μ = 4 * Nₑ * μ ∧
    (equilibriumPi μ c L = 4 * Nₑ * μ ↔ c * L * Nₑ = 1/4) :=
  ⟨unbiased_conversion_no_pressure,
   E_pairwise_pi_eq Nₑ μ,
   doc_eq_coalescent_iff Nₑ μ c L hμ hCL⟩

end LewontinParadox.Refutation
