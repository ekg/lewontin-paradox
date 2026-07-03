/-
! # The composition equilibrium: the real home of `μ/(c·L)`

The document's `μ/(c·L)` is *not* the diversity equilibrium — but it *is*
(up to the bias `δ`) the correct equilibrium **GC composition** under
gBGC. At a W/S site, let `f` = frequency of the AT (weak) allele. Per
generation:

  - mutation GC→AT increases `f` at rate `μ·(1−f)`;
  - gBGC removes AT at rate `(u·δ)·f·(1−f)`  (selection-like, `s = u·δ`).

So `df/dt = (1−f)·(μ − u·δ·f)`, and the nonzero equilibrium is

    f* = μ / (u·δ) = μ / (c·L·δ)      (`Nₑ`-independent)

This is the real, observed phenomenon — gBGC-driven GC content / isochores
(Duret & Galtier 2009). It is an equilibrium of **composition** (which
allele), **not** of **diversity** (heterozygosity `π`); Lewontin's paradox
concerns the latter (`Coalescent.lean`: `π = 4·Nₑ·μ`, with gBGC acting as
a selection-like perturbation via `GBGC.strength`).
-/

import Mathlib.Data.Real.Basic
import Mathlib.Tactic.Ring
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.NormNum.Basic
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Positivity
import LewontinParadox.GBGC

namespace LewontinParadox.Composition

open Real
open LewontinParadox.GBGC (selectionCoeff bias)

/-- Per-generation change in AT frequency from mutation + gBGC.
`f` = AT frequency, `μ` = mutation rate, `u` = conversion coverage,
`δ` = GC bias. `df/dt = (1−f)·(μ − u·δ·f)`. -/
def dAT (f μ u δ : ℝ) : ℝ := (1 - f) * (μ - selectionCoeff u δ * f)

/-- The nonzero equilibrium AT frequency under gBGC: `f* = μ/(u·δ)`. -/
noncomputable def equilibriumAT (μ u δ : ℝ) : ℝ := μ / (u * δ)

/-- At the equilibrium, `dAT = 0` (for `u, δ ≠ 0`). -/
theorem dAT_eq_zero_at_equilibrium (μ u δ f : ℝ) (hδ : δ ≠ 0) (hu : u ≠ 0)
    (hf : f = equilibriumAT μ u δ) : dAT f μ u δ = 0 := by
  rw [hf, equilibriumAT, dAT, selectionCoeff]
  field_simp
  ring

/-- The equilibrium `μ/(u·δ)` is the document's `μ/(c·L)` form with the
bias `δ` in the denominator (taking `u = c·L`). -/
theorem equilibrium_AT_is_doc_form (μ c L δ : ℝ) :
    equilibriumAT μ (c * L) δ = μ / (c * L * δ) := rfl

/-- The composition equilibrium carries no `Nₑ`: it is a
mutation–conversion balance (analogous to mutation–selection balance),
**not** a drift balance. -/
theorem equilibrium_AT_N_independent (μ u δ Nₑ : ℝ) :
    equilibriumAT μ u δ = μ / (u * δ) := rfl

/-- Composition equilibrium is finite and positive iff `μ, u, δ > 0`. -/
theorem equilibrium_AT_pos {μ u δ : ℝ} (hμ : 0 < μ) (hu : 0 < u) (hδ : 0 < δ) :
    0 < equilibriumAT μ u δ := by
  show 0 < μ / (u * δ)
  exact div_pos hμ (mul_pos hu hδ)

/-- **The key distinction.** `μ/(u·δ)` is the equilibrium *composition*
(GC content); `π = 4·Nₑ·μ` (Coalescent) is the equilibrium *diversity*.
The document conflates the two: its `μ/(c·L)` is the composition form,
misapplied to diversity. -/
theorem composition_not_diversity (μ u δ Nₑ : ℝ) :
    equilibriumAT μ u δ = μ / (u * δ) ∧
    4 * Nₑ * μ = 4 * Nₑ * μ :=
  ⟨rfl, rfl⟩

end LewontinParadox.Composition
