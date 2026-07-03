/-
! # GC-biased gene conversion (gBGC) as a selection-like force

The GC bias of gene conversion — transmit GC with probability `b > 1/2` at
W/S heterozygous sites — acts during *every* conversion spanning a W/S
site, whether allelic (homolog-to-homolog) or ectopic
(paralog-to-paralog): the biochemistry is the same. So the bias measured
from substitution patterns / sperm typing is the same bias that acts at
single-copy loci.

gBGC only manifests when a conversion spans the site. With `u` the
per-site conversion-coverage rate (prob a site is in a tract per
generation) and `δ = b − 1/2` the per-conversion GC bias, the effective
per-site per-generation selection coefficient is `s_eff = u·δ`, with
scaled strength `γ = 4·Nₑ·u·δ`.

This is **linear in `Nₑ`** (selection-like), so it does *not* produce the
document's flat, `Nₑ`-independent `μ/(c·L)`; it produces BGS-like
sub-linear compression of *diversity*, and is bounded by self-limitation
(`Saturation.lean`). The document's `μ/(c·L)` is instead the
**composition** equilibrium (`Composition.lean`), not the diversity
equilibrium.
-/

import Mathlib.Data.Real.Basic
import Mathlib.Tactic.Ring
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.NormNum.Basic
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Positivity
import LewontinParadox.Heterozygosity
import LewontinParadox.Selection

namespace LewontinParadox.GBGC

open Real
open LewontinParadox.Selection (scaledSelection reductionFactor)

/-- Per-conversion GC bias `δ = b − 1/2` (`b = P(transmit GC | W/S, conv)`). -/
noncomputable def bias (b : ℝ) : ℝ := b - 1/2

/-- Effective per-site per-generation gBGC selection `s_eff = u·δ`
(`u` = conversion-coverage rate, `δ` = per-conversion GC bias). -/
def selectionCoeff (u δ : ℝ) : ℝ := u * δ

/-- Scaled gBGC strength `γ = 4·Nₑ·u·δ`. -/
def strength (Nₑ u δ : ℝ) : ℝ := 4 * Nₑ * u * δ

/-- gBGC strength equals the generic scaled-selection at `s = u·δ`. -/
theorem strength_eq_scaledSelection (Nₑ u δ : ℝ) :
    strength Nₑ u δ = scaledSelection Nₑ (selectionCoeff u δ) := by
  dsimp [strength, scaledSelection, selectionCoeff]; ring

/-- gBGC exerts a force iff there is a nonzero bias (`δ = 0`), for
`Nₑ, u > 0`. (Unbiased conversion — `δ = 0` — is the martingane case of
`Transmission.lean`.) -/
theorem strength_zero_iff_no_bias {Nₑ u δ : ℝ} (hN : 0 < Nₑ) (hu : 0 < u) :
    strength Nₑ u δ = 0 ↔ δ = 0 := by
  dsimp [strength]
  constructor
  · intro h
    have hkey : (4 * Nₑ * u) * δ = 0 := h
    rcases mul_eq_zero.mp hkey with hz | hδ
    · exfalso; nlinarith [hN, hu, hz]
    · exact hδ
  · intro h
    rw [h]; ring

/-- gBGC strength is **linear in `Nₑ`**: `γ = (4·u·δ)·Nₑ`. Hence it is
`Nₑ`-dependent and cannot produce the document's flat, `Nₑ`-independent
`μ/(c·L)` equilibrium for *diversity*. -/
theorem strength_linear_in_Nₑ (Nₑ u δ : ℝ) :
    strength Nₑ u δ = (4 * u * δ) * Nₑ := by dsimp [strength]; ring

/-- With fixed `u, δ`, doubling `Nₑ` doubles gBGC strength (selection-like,
unlike the document's claimed `Nₑ`-independence). -/
theorem strength_doubles_with_Nₑ (u δ Nₑ : ℝ) :
    strength (2 * Nₑ) u δ = 2 * strength Nₑ u δ := by
  dsimp [strength]; ring

end LewontinParadox.GBGC
