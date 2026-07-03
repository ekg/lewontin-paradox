/-
! # Finite populations, drift, and the tetrad subtlety

The coalescent (`Coalescent.lean`) is *already* the finite-`N` theory: the
Kingman rate `1/(2·Nₑ)` **is** genetic drift. So `E[π] = 4·Nₑ·μ` is not an
infinite-population result — it is the mutation–drift balance in a finite
population. "What if `N` is finite?" is not an escape hatch.

The reason conversion does not add to drift at a single site is the
per-gamete fact (`Transmission.lean`): a single gamete from an `Aa`
meiosis is `Bernoulli(1/2)` whether or not a conversion spans the site
(`(1/2)·(1/4) + (1/2)·(3/4) = 1/2`). So the per-gamete variance is `1/4`
in both cases, the binomial sampling variance `2N·p(1−p)` is unchanged,
and drift stays `1/(2·Nₑ)`, `c`-independent.

There is ONE finite-`N` setting where conversion *does* add variance:
**tetrad preservation** — if all four products of one meiosis survive
together in the next generation, the `3:1`/`1:3` tetrad distortion is not
washed out by sampling. Mendelian tetrads are always exactly `2A:2a`
(zero variance in `#A`); conversion tetrads are `1A:3a` or `3A:1a`
(variance `1`). So conversion adds drift variance *under tetrad
preservation*. But (a) this is non-standard (standard WF / the coalescent
sample one gamete per meiosis), (b) it is a second-order *variance*
(slightly reduced `Nₑ`), not the first-order directional homogenization
of `main.md`, and (c) it does not yield `μ/(c·L)`. It is the same
"linkage" spirit as repeats (`Repeats.lean`): correlated co-segregation of
copies lets conversion's variance survive.

Statements are scaled by `2` to avoid noncomputable division by `2`.
-/

import Mathlib.Data.Real.Basic
import Mathlib.Tactic.Ring
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.NormNum.Basic
import LewontinParadox.Transmission
import LewontinParadox.Coalescent

namespace LewontinParadox.Drift

open Real

/-- Per-gamete probability of transmitting `A` from an `Aa` heterozygote:
`1/2` whether or not a conversion spans the site (the martingale fact of
`Transmission.lean`). -/
noncomputable def perGameteAProb : ℝ := 1/2

/-- Per-gamete (Bernoulli) variance `p·(1−p)`. -/
noncomputable def perGameteVar : ℝ := perGameteAProb * (1 - perGameteAProb)

/-- `4 · perGameteVar = 1`, i.e. per-gamete variance is `1/4` — the same
for Mendelian and conversion meioses (since `P(A | Aa) = 1/2` in both).
Hence standard-WF binomial drift `2N·p(1−p)` is unchanged by conversion. -/
theorem perGameteVar_eq : 4 * perGameteVar = 1 := by
  dsimp [perGameteVar, perGameteAProb]; norm_num

/-- Mendelian tetrad `#A`: always exactly `2`. -/
def mendelianTetradA : ℝ := 2

/-- Conversion tetrad `#A` as a function of DSB location: `1` (DSB on the
A-homolog) or `3` (DSB on the a-homolog). -/
def conversionTetradA (dsbOnA : Bool) : ℝ := if dsbOnA then 1 else 3

/-- The two conversion-tetrad outcomes sum to `4 = 2·2`, so the mean is
`2` — equal to the Mendelian tetrad (no directional effect). -/
theorem conversionTetrad_mean_sum : (1 : ℝ) + 3 = 2 * mendelianTetradA := by
  dsimp [mendelianTetradA]; norm_num

/-- Mendelian tetrad `#A = 2`. -/
theorem mendelianTetrad_mean : mendelianTetradA = 2 := rfl

/-- Mendelian tetrad variance (in `#A`): `0` (always `2`). -/
def mendelianTetradVar : ℝ := 0

/-- Conversion tetrad variance (in `#A`): `E[#A²] − E[#A]² = 1`. -/
noncomputable def conversionTetradVar : ℝ :=
  (1/2) * 1^2 + (1/2) * 3^2 - 2^2

/-- `2 · conversionTetradVar = 2`, i.e. the conversion-tetrad variance is
`1` (vs `0` for Mendelian). -/
theorem conversionTetradVar_eq : 2 * conversionTetradVar = 2 := by
  dsimp [conversionTetradVar]; norm_num

/-- The conversion-tetrad variance is exactly `1`. -/
theorem conversionTetradVar_eq_one : conversionTetradVar = 1 := by
  have h : 2 * conversionTetradVar = 2 := conversionTetradVar_eq
  have h2 : (2 : ℝ) ≠ 0 := by norm_num
  exact mul_left_cancel₀ h2 (by linarith [h] : 2 * conversionTetradVar = 2 * 1)

/-- **Conversion adds tetrad variance** (`1 > 0` vs `0`): under tetrad
preservation, conversion increases drift. (Under standard WF — one gamete
per meiosis — this is washed out: per-gamete variance is `1/4` in both.) -/
theorem conversion_adds_tetrad_variance :
    0 < conversionTetradVar ∧ mendelianTetradVar = 0 := by
  refine ⟨?_, ?_⟩
  · rw [conversionTetradVar_eq_one]; norm_num
  · rfl

/-- The standard-WF drift at a single site is set by the per-gamete
variance (`1/4`, equal for Mendelian and conversion), hence the coalescent
marginal is Kingman (`1/(2·Nₑ)`), `c`-independent — i.e.
`Coalescent.E_pi_independent_of_c`. -/
theorem standard_WF_drift_unchanged_by_conversion :
    4 * perGameteVar = 1 := perGameteVar_eq

/-- Even under tetrad preservation, the added effect is a *variance*
(slightly reduced `Nₑ`), not a directional force: both tetrads have the
same mean (`2`), so conversion changes `Var(Δp)`, not `E[Δp]`. It cannot
produce the directional homogenization `Rate_loss = c·L·H²` of `main.md`. -/
theorem no_directional_effect_even_with_tetrads :
    (1 : ℝ) + 3 = 2 * mendelianTetradA := conversionTetrad_mean_sum

end LewontinParadox.Drift
