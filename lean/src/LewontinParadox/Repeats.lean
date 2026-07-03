/-
! # Why repeats homogenize but single-copy does not (unbiased conversion)

The entire difference is **linkage**: are the two copies being compared on
the *same* chromosome (co-segregating paralogs in a tandem array) or on
*different* homologs (the two alleles at a single locus)?

For two gene copies the pairwise coalescence rate is

    λ = 1/(2·Nₑ)  +  g_link

where `g_link` is the rate at which a conversion event **directly joins
these two specific lineages**.

* **Single-copy alleles** sit on different homologs; going back a
  generation they are in different parents and only meet via the organism
  coalescent. Allelic conversion happens *within one meiosis* (pairing one
  individual's two homologs) and so never directly joins two alleles living
  in different individuals — it only reroutes a lineage between homolog
  backgrounds. Hence `g_link = 0`: the marginal is Kingman (Wiuf & Hein
  2000), `π = 4·Nₑ·μ`. (Forward-time: outcrossing resegregates the homologs
  each generation, so `P(A | Aa) = ½` and no within-pair copying
  accumulates — see `Transmission.lean`.)

* **Linked paralogs** `i, j` sit on the *same* chromosome and co-segregate
  as a block; an ectopic conversion between positions `i` and `j` makes
  copy `i`'s ancestry derive from copy `j`'s, directly joining the two
  lineages at rate `g`. Hence `g_link = g`, `E[T_MRCA] ≈ 1/g`, and
  `π_within ≈ 2μ/g` (conversion-limited, `Nₑ`-independent). This is the
  homogenization — and it requires no GC bias.

Caveat: this compresses *within-individual, between-copy* diversity only.
The *between-individual* diversity at a fixed repeat position (copy `i` in
individual 1 vs copy `i` in individual 2) is still the single-copy marginal
(`g_link = 0`, Kingman, `4·Nₑ·μ`): those two copies are on different
homologs in different individuals. So even for repeats, species-level `π`
is not compressed; only the within-array diversity is.
-/

import Mathlib.Data.Real.Basic
import Mathlib.Tactic.Ring
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.NormNum.Basic
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Positivity
import LewontinParadox.Heterozygosity
import LewontinParadox.Coalescent

namespace LewontinParadox.Repeats

open Real
open LewontinParadox.Heterozygosity (theta)
open LewontinParadox.Coalescent (E_Tmrca)

/-- Pairwise coalescence rate: organism rate `1/(2·Nₑ)` plus the
conversion rate `g_link` that directly joins the two lineages. -/
noncomputable def coalRate (Nₑ g_link : ℝ) : ℝ := 1 / (2 * Nₑ) + g_link

/-- Expected pairwise `T_MRCA = 1/λ`. -/
noncomputable def eTmrca (Nₑ g_link : ℝ) : ℝ := 1 / coalRate Nₑ g_link

/-- Expected pairwise diversity `2·μ·E[T_MRCA]`. -/
noncomputable def ePi (μ Nₑ g_link : ℝ) : ℝ := 2 * μ * eTmrca Nₑ g_link

/-- **Single-copy** (`g_link = 0`): `π = 4·Nₑ·μ`, the standard neutral
result, with no conversion term. -/
theorem singleCopy_pi (μ Nₑ : ℝ) (hN : 0 < Nₑ) :
    ePi μ Nₑ 0 = theta Nₑ μ := by
  dsimp [ePi, eTmrca, coalRate, theta]
  field_simp
  ring

/-- **Linked paralogs** (`g_link = g`): `π = 2μ / (1/(2Nₑ) + g)`. For
`g ≫ 1/(2Nₑ)` this is `≈ 2μ/g` (conversion-limited, `Nₑ`-independent). -/
theorem paralog_pi (μ Nₑ g : ℝ) (hN : 0 < Nₑ) (hg : 0 ≤ g) :
    ePi μ Nₑ g = 2 * μ / (1 / (2 * Nₑ) + g) := by
  dsimp [ePi, eTmrca, coalRate]
  field_simp

/-- The linking term strictly lowers within-individual between-copy
diversity: repeats homogenize relative to single-copy. -/
theorem paralog_pi_le_singleCopy (μ Nₑ g : ℝ)
    (hN : 0 < Nₑ) (hg : 0 < g) (hμ : 0 < μ) :
    ePi μ Nₑ g ≤ ePi μ Nₑ 0 := by
  dsimp [ePi, eTmrca, coalRate]
  have hR : 0 < 1 / (2 * Nₑ) := by positivity
  have hL : 0 < 1 / (2 * Nₑ) + g := by positivity
  have hL0 : 0 < 1 / (2 * Nₑ) + 0 := by positivity
  have hAB : 1 / (2 * Nₑ) ≤ 1 / (2 * Nₑ) + g := by linarith
  field_simp
  nlinarith [hR, hL, hL0, hAB, hμ, hg, hN]

/-- In the conversion-limited regime `g ≫ 1/(2·Nₑ)`, within-individual
between-copy diversity approaches `2μ/g`, independent of `Nₑ`. -/
theorem paralog_pi_approx (μ Nₑ g : ℝ) (hN : 0 < Nₑ) (hg : 0 < g) :
    ePi μ Nₑ g = 2 * μ / g - 2 * μ / (g * (1 + 2 * Nₑ * g)) := by
  -- 2μ / (1/(2N) + g) = 2μ/g · 1/(1 + 1/(2Ng)) = 2μ/g · (2Ng)/(1+2Ng)
  -- = 2μ/g · (1 - 1/(1+2Ng)) = 2μ/g - 2μ/(g(1+2Ng))
  dsimp [ePi, eTmrca, coalRate]
  field_simp
  ring

end LewontinParadox.Repeats
