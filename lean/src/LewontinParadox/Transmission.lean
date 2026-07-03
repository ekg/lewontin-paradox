/-
! # Transmission distortion: does conversion push alleles to fixation?

This is the central question that decides whether gene conversion is a
"homogenization pressure" (the document's §5.2) or merely a recombination
mechanism. A gene-conversion tract spanning a heterozygous site changes
the segregation ratio of the two alleles **in that one meiosis**. Whether
that change has a *sign* (a directional / mean effect on allele
frequency) is determined entirely by **bias**:

* **DSB asymmetry** (the "hotspot paradox" / meiotic drive): if the DSB
  preferentially falls on one homolog, the other homolog's allele is
  over-transmitted. This is locus-specific and self-limiting (it drives
  hotspots to extinction); it is *not* a genome-wide force scaling with
  standing diversity.
* **gBGC** (allele-content bias): conversion preferentially transmits the
  GC allele. This is formally equivalent to genic selection — directional,
  sign-dependent, type-specific — not a uniform homogenization of IBD
  tracts.
* **No bias** (the default, symmetric case): the mean transmission
  distortion is exactly zero. The allele-frequency process is a
  martingale with respect to conversion; there is **no** directional
  homogenization pressure, and hence no `c`-driven, `Nₑ`-independent
  equilibrium.

(Precise variance point: a *single* gamete from an `Aa` parent is
`Bernoulli(1/2)` whether or not a conversion spans the site, since
`(1/2)·(1/4) + (1/2)·(3/4) = 1/2`. The tetrad-level `3:1` distortion
exists, but in Wright–Fisher sampling gametes are pooled, so that
intra-meiosis variance is washed out — and the coalescent confirms there
is *no* effect on `E[π]`, not even a small one. The donor template is the
sister homolog, not a random population allele, so it is `A` or `a` with
probability `1/2` regardless of population frequency `p`; there is no
frequency-weighting toward the common allele. A genuine homogenization
force of the kind intuited here exists only for *non-allelic / ectopic*
gene conversion — concerted evolution of multigene families — which is
not the mechanism behind single-copy neutral diversity.)

For ease of exact reasoning, the transmission distortions are scaled by a
positive constant (`4` for DSB asymmetry, `2` for gBGC) so that all
coefficients are integers; the *sign* and *zero-set* — the only things
that matter for "is there a directional pressure?" — are preserved.
-/

import Mathlib.Data.Real.Basic
import Mathlib.Tactic.Ring
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.NormNum.Basic
import LewontinParadox.Heterozygosity

namespace LewontinParadox.Transmission

open Real

/-! ## DSB-asymmetry model

For a heterozygote `Aa`, a conversion event spanning the focal site in one
meiosis produces a `3:1` (or `1:3`) tetrad instead of the Mendelian `2:2`,
depending on which homolog received the DSB. The transmission distortion
`δ = P(transmit A) − 1/2` is `−1/4` (DSB on A) or `+1/4` (DSB on a); the
`4×`-scaled value is `−1` or `+1`.
-/

/-- `4×` the transmission distortion `δ`, as a function of which homolog
was cut. `bias` is `Pr(DSB on A-homolog)`; `bias = 1/2` is the unbiased
case. -/
noncomputable def distortion4 (bias : ℝ) : ℝ := 1 - 2 * bias

/-- Closed form already: `4·E[δ] = 1 − 2·bias`. -/
theorem distortion4_eq (bias : ℝ) : distortion4 bias = 1 - 2 * bias := rfl

/-- **Key fact (unbiased case).** With no DSB asymmetry (`bias = 1/2`),
the mean transmission distortion is exactly zero: gene conversion does not
change the expected allele frequency of a segregating allele. -/
theorem expectedDistortion_unbiased : distortion4 (1/2) = 0 := by
  dsimp [distortion4]; norm_num

/-- Mean transmission distortion vanishes iff there is no DSB bias. -/
theorem expectedDistortion_eq_zero_iff (bias : ℝ) :
    distortion4 bias = 0 ↔ bias = 1/2 := by
  dsimp [distortion4]
  constructor
  · intro h; linarith
  · intro h; linarith

/-- In a population at allele frequency `p`, a conversion event at a
heterozygous locus perturbs the gametic allele frequency by the
transmission distortion `δ`. With no DSB bias the **expected** change is
zero, so the allele-frequency process is a martingale w.r.t. conversion:
no directional homogenization. -/
theorem expected_delta_p_zero_unbiased (p : ℝ) :
    (distortion4 (1/2) : ℝ) = 0 := expectedDistortion_unbiased

/-! ## GC-biased gene conversion (gBGC)

At a heterozygous W/S (AT/GC) site, conversion is biased toward
transmitting the GC allele with probability `b ∈ [1/2, 1]`. This is a
directional, allele-content-dependent bias — formally equivalent to genic
selection — **not** a uniform homogenization of IBD tracts. The `2×`-scaled
distortion is `2·b − 1`.
-/

/-- `2×` the gBGC distortion `δ_GC = b − 1/2`: `gBGCdistortion2 b = 2·b − 1`.
Zero iff `b = 1/2`. -/
def gBGCdistortion2 (b : ℝ) : ℝ := 2 * b - 1

/-- gBGC exerts a nonzero mean effect on allele frequency iff `b ≠ 1/2`. -/
theorem gBGC_distortion2_eq_zero_iff (b : ℝ) :
    gBGCdistortion2 b = 0 ↔ b = 1/2 := by
  dsimp [gBGCdistortion2]
  constructor
  · intro h; linarith
  · intro h; linarith

/-- In the unbiased case (`b = 1/2`) gBGC also exerts no mean pressure. -/
theorem gBGC_distortion2_unbiased : gBGCdistortion2 (1/2) = 0 := by
  dsimp [gBGCdistortion2]; norm_num

end LewontinParadox.Transmission
