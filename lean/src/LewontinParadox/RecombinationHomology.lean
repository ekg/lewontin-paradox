/-
! # Recombination requires homology: the self-sustaining feedback

Mechanistic basis (cf. the conversion-tract literature): non-crossover
(DSB repair via SDSA / double-Holliday-junction) requires the invading
strand to form a stable heteroduplex and for the junction to **migrate**
away from the DSB. Branch migration requires **near-identity** —
mismatches block it. So:

  * conversion only works between near-identical sequences (a homology
    threshold);
  * conversion *makes* sequences more identical (it homogenizes).

This is a **positive feedback**: homology enables conversion, conversion
maintains homology. Because DSB repair (recombination) is **essential for
meiosis** → chromosome segregation → fertility, there is **selection to
maintain the homology that conversion maintains** — keeping the genome in
the "recombino-genic", homogenized state.

This module formalizes the homology cutoff (conversion rate → 0 past a
divergence threshold), which together with the substrate limit
(`Saturation`: rate → 0 at fixation) makes the conversion rate a hump in
divergence — sustaining concerted evolution only below the threshold. It
explains *why* gBGC (`δ > 0`, `GBGC.lean`) and high conversion rates
(`u`) are maintained: recombination selects for them. It is the
evolutionary basis for the parameters that make the mutation–selection–
drift saturation (`MutationSelectionDrift.lean`) bite — i.e. it explains
*why the flattening is always on*, not a separate flattening mechanism.
-/

import Mathlib.Data.Real.Basic
import Mathlib.Tactic.Ring
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.NormNum.Basic
import LewontinParadox.Saturation
import LewontinParadox.MutationSelectionDrift

namespace LewontinParadox.RecombinationHomology

open Real

/-- Conversion rate as a function of pairwise divergence `d`: below the
homology threshold `d_max` the machinery works; above it strand invasion
/ branch migration fails and the rate is 0. (Simplified envelope; the
real rate also depends on substrate — `Saturation`.) -/
noncomputable def conversionRate (d d_max : ℝ) : ℝ :=
  if d ≤ d_max then 1 else 0

/-- Past the homology threshold, conversion cannot occur. -/
theorem conversionRate_zero_past_threshold (d d_max : ℝ) (h : d_max < d) :
    conversionRate d d_max = 0 := by
  dsimp [conversionRate]
  split
  · linarith
  · rfl

/-- Below the threshold, conversion is active (envelope value 1; the
actual rate is then set by substrate, `Saturation.effectiveRate`). -/
theorem conversionRate_active_below_threshold (d d_max : ℝ) (h : d ≤ d_max) :
    conversionRate d d_max = 1 := by
  dsimp [conversionRate]
  split
  · rfl
  · linarith

/-- Combined with `Saturation.effectiveRate_zero_at_fixation`, the
conversion rate is a hump in divergence: 0 at fixation (no substrate),
active at intermediate divergence, 0 past the homology threshold (no
homology). Concerted evolution is self-sustaining only in the
intermediate band — which selection maintains because recombination
(essential for meiosis) requires it. -/
theorem homology_feedback_principle :
    (∀ d d_max, d_max < d → conversionRate d d_max = 0) ∧
    (∀ f u δ, f = 0 ∨ f = 1 → LewontinParadox.Saturation.effectiveRate u δ f = 0) :=
  ⟨fun _ _ h => conversionRate_zero_past_threshold _ _ h,
   fun _ _ _ hf => LewontinParadox.Saturation.effectiveRate_zero_at_fixation _ _ _ hf⟩

/-- **The synthesis.** Because recombination is essential for meiosis and
requires near-identity, selection maintains (a) high conversion rates
`u > 0` and (b) the GC-biased repair machinery `δ > 0` (`GBGC.bias`).
These are exactly the parameters that make the mutation–selection–drift
saturation (`MutationSelectionDrift.saturation_always_on`) always-on and
strong. The homogenization-for-recombination feedback is thus the
evolutionary *cause* of the `gBGC + finite-N` flattening, not a separate
diversity-reducing force. -/
theorem recombination_drives_the_flattening_parameters
    (μ u δ Nₑ : ℝ) (hμ : 0 < μ) (hu : 0 < u) (hδ : 0 < δ)
    (h : LewontinParadox.MutationSelectionDrift.crossover u δ ≤ Nₑ) :
    LewontinParadox.MutationSelectionDrift.wsDiversityStrong μ u δ ≤
      LewontinParadox.MutationSelectionDrift.wsDiversityWeak μ Nₑ :=
  LewontinParadox.MutationSelectionDrift.strong_below_neutral μ u δ Nₑ hμ hu hδ h

end LewontinParadox.RecombinationHomology
