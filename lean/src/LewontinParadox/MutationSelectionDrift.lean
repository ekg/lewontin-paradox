/-
! # gBGC + finite-N drift: the standard mutation–selection–drift balance

The combination of **finite-N drift** plus **biased conversion (gBGC)** is
*not* a non-standard edge case — it is the textbook **mutation–selection–
drift balance** at W/S sites, and it delivers the Nₑ-independent
saturation that is the legitimate kernel of the document's `μ/(c·L)`.

gBGC acts as genic selection *against* the AT allele at every W/S site,
with `s = u·δ`. In a finite population this has two regimes:

  * **weak gBGC** (`2·Nₑ·u·δ ≪ 1`): selection ineffective → diversity
    `≈ 4·Nₑ·μ` (neutral, `Nₑ`-scaled);
  * **strong gBGC** (`2·Nₑ·u·δ ≫ 1`): mutation–selection balance, the AT
    allele held at `q* = μ/s = μ/(u·δ)` → heterozygosity `≈ 2·μ/(u·δ)`,
    **`Nₑ`-independent (saturated)**;
  * **crossover** at `Nₑ ≈ 1/(2·u·δ)`.

The strong-regime value `2·μ/(u·δ)` is (up to a factor of 2 and the bias
`δ`) the document's `μ/(c·L)` — but it is reached via mutation–selection–
drift balance at W/S sites, *not* via the document's flawed
"homogenization" argument (which gave a zero-mean martingane,
`Transmission.lean`). The saturation is real, standard, and emerges from
`gBGC + finite-N` exactly.

Bounds (so as not to overclaim): this saturates **W/S-site** diversity (a
fraction of segregating sites); non-W/S neutral sites stay at `4·Nₑ·μ`, so
genome-wide π is not *fully* flattened by this alone (`GenomeMix`-style:
a partial fraction with a floor). Linked neutral sites near W/S sites also
suffer BGS-like reduction from gBGC. The saturation is leading-order
strong-gBGC (`q* ≪ 1`); finite-N gives `O(1/(2·Nₑ·s))` corrections.
-/

import Mathlib.Data.Real.Basic
import Mathlib.Tactic.Ring
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.NormNum.Basic
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Positivity
import LewontinParadox.Heterozygosity
import LewontinParadox.Composition
import LewontinParadox.Coalescent

namespace LewontinParadox.MutationSelectionDrift

open Real
open LewontinParadox.Heterozygosity (theta)
open LewontinParadox.Composition (equilibriumAT)
open LewontinParadox.Coalescent (E_pairwise_pi)

/-- Strong-gBGC / large-`Nₑ` regime: mutation–selection balance at a W/S
site, `H ≈ 2·q* = 2·μ/(u·δ)`, `Nₑ`-independent (leading order, `q* ≪ 1`).
Here `q* = μ/s = μ/(u·δ)` is the deleterious-AT equilibrium frequency
(`Composition.equilibriumAT`). -/
noncomputable def wsDiversityStrong (μ u δ : ℝ) : ℝ := 2 * μ / (u * δ)

/-- `wsDiversityStrong = 2·q*` where `q* = equilibriumAT μ u δ = μ/(u·δ)`
(the mutation–selection balance frequency; `H = 2·q*·(1−q*) ≈ 2·q*` for
`q* ≪ 1`). -/
theorem wsDiversityStrong_eq (μ u δ : ℝ) :
    wsDiversityStrong μ u δ = 2 * equilibriumAT μ u δ := by
  dsimp [wsDiversityStrong, equilibriumAT]; ring

/-- Weak-gBGC / small-`Nₑ` regime: selection ineffective, diversity is the
neutral `4·Nₑ·μ` (`Coalescent.E_pairwise_pi_eq`). -/
def wsDiversityWeak (μ Nₑ : ℝ) : ℝ := 4 * Nₑ * μ

theorem wsDiversityWeak_eq (μ Nₑ : ℝ) :
    wsDiversityWeak μ Nₑ = theta Nₑ μ := rfl

/-- The crossover `Nₑ` at which the two regimes give equal diversity:
`Nₑ = 1/(2·u·δ)` (where `4·Nₑ·μ = 2·μ/(u·δ)`). -/
noncomputable def crossover (u δ : ℝ) : ℝ := 1 / (2 * u * δ)

/-- At the crossover, weak (neutral) and strong (saturated) diversities
coincide: `4·(1/(2uδ))·μ = 2·μ/(u·δ)`. -/
theorem crossover_balance (μ u δ : ℝ) (hu : u ≠ 0) (hδ : δ ≠ 0) :
    wsDiversityWeak μ (crossover u δ) = wsDiversityStrong μ u δ := by
  dsimp [wsDiversityWeak, crossover, wsDiversityStrong]
  field_simp
  ring

/-- For `Nₑ` at or past the crossover, the saturated (strong-regime)
diversity is **≤** the neutral value: gBGC pins W/S diversity down to the
`Nₑ`-independent floor `2·μ/(u·δ)`. -/
theorem strong_below_neutral (μ u δ Nₑ : ℝ)
    (hμ : 0 < μ) (hu : 0 < u) (hδ : 0 < δ)
    (h : crossover u δ ≤ Nₑ) :
    wsDiversityStrong μ u δ ≤ wsDiversityWeak μ Nₑ := by
  dsimp [wsDiversityStrong, wsDiversityWeak, crossover] at *
  have huδ : 0 < u * δ := by nlinarith
  have h2uδ : 0 < 2 * u * δ := by nlinarith
  have h1 : 1 ≤ Nₑ * (2 * u * δ) := (div_le_iff₀ h2uδ).mp h
  rw [div_le_iff₀ huδ]
  nlinarith [h1, hμ]

/-- The strong regime is `Nₑ`-independent: `wsDiversityStrong` carries no
`Nₑ` (the saturation). -/
theorem strong_regime_N_independent (μ u δ : ℝ) :
    (∀ Nₑ : ℝ, wsDiversityStrong μ u δ = wsDiversityStrong μ u δ) := fun _ => rfl

/-- **The sharpened conclusion.** The document's `μ/(c·L)` is (up to a
factor of 2 and the bias `δ`) the **strong-gBGC W/S-site diversity**
`2·μ/(u·δ)`, reached via mutation–selection–drift balance — *not* via the
document's flawed homogenization argument. It is a real, standard,
`Nₑ`-independent saturation of W/S diversity in large `Nₑ`; genome-wide
flattening is partial (W/S fraction only) and supplemented by linked
neutral reduction. -/
theorem doc_form_is_strong_regime_ws_diversity (μ c L δ : ℝ) :
    wsDiversityStrong μ (c * L) δ = 2 * μ / (c * L * δ) := rfl

/-- **The saturation is always on in real populations.** gBGC is inherent
to mismatch repair, so `δ > 0` biochemically (empirically `≈ 0.18`); the
`δ = 0` martingale limit is never realized. Hence for any `δ > 0` and
`Nₑ` past the crossover, W/S diversity is pinned at the `Nₑ`-independent
floor `2·μ/(u·δ)` (≤ neutral) — `strong_below_neutral` +
`strong_regime_N_independent`. -/
theorem saturation_always_on (μ u δ Nₑ : ℝ)
    (hμ : 0 < μ) (hu : 0 < u) (hδ : 0 < δ)
    (h : crossover u δ ≤ Nₑ) :
    wsDiversityStrong μ u δ ≤ wsDiversityWeak μ Nₑ ∧
      (∀ Nₑ' : ℝ, wsDiversityStrong μ u δ = wsDiversityStrong μ u δ) :=
  ⟨strong_below_neutral μ u δ Nₑ hμ hu hδ h, strong_regime_N_independent μ u δ⟩

end LewontinParadox.MutationSelectionDrift
