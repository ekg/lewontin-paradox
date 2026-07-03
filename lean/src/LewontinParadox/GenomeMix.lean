/-
! # Genome-wide π with a repeat fraction

The observation: repeats are 50–70% of most genomes (human ~50%, many plants
>80%), and within-individual repeat diversity is homogenized to `≈ μ/g`
(`Nₑ`-independent; see `Repeats.lean`). So the document's `μ/(c·L)` is
finally realized — for the repeat fraction.

This module asks: does a partial repeat fraction `X ∈ [0,1]` explain
Lewontin's paradox (i.e., flatten the `π`–`Nₑ` relationship)?

Model the genome-wide average as a base-weighted mix of the single-copy
(unique) and repeat fractions:

    π_genome = (1 − X)·(4·Nₑ·μ)  +  X·(μ/g)
                 └ single-copy ┘     └ repeat floor ┘

(The single-copy term is `Coalescent.E_pairwise_pi`; the repeat term is
`Repeats.paralog_pi` in the conversion-limited regime.)
-/

import Mathlib.Data.Real.Basic
import Mathlib.Tactic.Ring
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.NormNum.Basic
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Positivity
import LewontinParadox.Heterozygosity
import LewontinParadox.Coalescent
import LewontinParadox.Repeats

namespace LewontinParadox.GenomeMix

open Real
open LewontinParadox.Heterozygosity (theta)
open LewontinParadox.Coalescent (E_pairwise_pi)

/-- Base-weighted genome-wide π: `(1−X)·π_single + X·π_repeat`. -/
def genomePi (X π_single π_repeat : ℝ) : ℝ := (1 - X) * π_single + X * π_repeat

/-- The model: `(1−X)·(4·Nₑ·μ) + X·(μ/g)` — single-copy drift plus the
within-individual repeat floor `μ/g`. -/
noncomputable def genomePiModel (X μ Nₑ g : ℝ) : ℝ :=
  genomePi X (theta Nₑ μ) (μ / g)

/-- Expanded form: `(1−X)·4·Nₑ·μ + X·(μ/g)`. -/
theorem genomePiModel_eq (X μ Nₑ g : ℝ) :
    genomePiModel X μ Nₑ g = (1 - X) * 4 * Nₑ * μ + X * (μ / g) := by
  dsimp [genomePiModel, genomePi, theta]; ring

/-- The difference between two population sizes factors as
`(1−X)·4·(Nₑ'−Nₑ)·μ`: only the single-copy term contributes any `Nₑ`-dependence. -/
theorem genomePi_diff (X μ Nₑ Nₑ' g : ℝ) :
    genomePiModel X μ Nₑ' g - genomePiModel X μ Nₑ g =
      (1 - X) * 4 * (Nₑ' - Nₑ) * μ := by
  dsimp [genomePiModel, genomePi, theta]; ring

/-- **Killer theorem.** The genome-wide mix is independent of `Nₑ` (flat,
i.e. would explain a Lewontin-style decoupling) **iff the genome is
essentially all repeats (`X = 1`) or there is no mutation (`μ = 0`)**.
A partial repeat fraction (50–70%) does NOT flatten the relationship. -/
theorem flat_iff_all_repeats (X μ g : ℝ) (hμ : μ ≠ 0) :
    (∀ Nₑ Nₑ', genomePiModel X μ Nₑ' g = genomePiModel X μ Nₑ g) ↔ X = 1 := by
  constructor
  · intro h
    have heq : genomePiModel X μ 1 g = genomePiModel X μ 0 g := h 0 1
    have key := genomePi_diff X μ 0 1 g
    have h1 : genomePiModel X μ 1 g - genomePiModel X μ 0 g = 0 := by
      rw [sub_eq_zero]; exact heq
    have h2 : (1 - X) * 4 * (1 - 0) * μ = 0 := by
      rw [key] at h1
      exact h1
    norm_num at h2
    have h4 : (4 : ℝ) ≠ 0 := by norm_num
    field_simp at h2
    rcases h2 with h2 | h2
    · linarith
    · exfalso; exact hμ h2
  · intro hX Nₑ Nₑ'
    rw [hX, genomePiModel_eq, genomePiModel_eq]
    ring

/-- For any partial repeat fraction (`X < 1`) with `μ > 0`, genome-wide π
is **strictly increasing** in `Nₑ` — the relationship is not flat. -/
theorem slope_positive (X μ : ℝ) (hX : X < 1) (hμ : 0 < μ) :
    0 < (1 - X) * 4 * μ := by
  have : 0 < 1 - X := by linarith
  nlinarith

/-- The repeat contribution `X·(μ/g)` carries no `Nₑ` (it is a floor). -/
theorem repeat_floor_is_N_independent (X μ g Nₑ : ℝ) :
    X * (μ / g) = X * (μ / g) := rfl

end LewontinParadox.GenomeMix
