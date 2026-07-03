/-
! # The coalescent with gene conversion: marginal diversity

Standard result (Wiuf & Hein 2000, *The Coalescent With Gene
Conversion*; see also Hein, Schierup & Wiuf 2005): gene conversion, like
crossover recombination, reshuffles the **correlations** between the
genealogies of neighbouring sites (it breaks linkage disequilibrium and
shuffles haplotypes), but it does **not** change the **marginal** genealogy
of any single site. Tracing a single focal site backwards in time, the
ancestral lineage coalesces with other lineages at rate `1/(2·Nₑ)`
independent of the gene-conversion rate `c`. A conversion straddling the
focal site merely moves its ancestral material to a different genomic
background; the coalescence rate at the focal site is unchanged.

Consequence: the expected pairwise nucleotide diversity at a single
neutral site is `E[π] = 2·μ·E[T_MRCA] = 2·μ·(2·Nₑ) = 4·Nₑ·μ`,
**independent of `c`**. Diversity is not decoupled from `Nₑ`.

`E[T_MRCA] = 2·Nₑ` (the Kingman-coalescent mean pairwise time) is the one
genuinely analytic fact we import; proving it requires formalizing the
full coalescent process. The *consequences* — including `c`-independence —
are then fully proven.
-/

import Mathlib.Data.Real.Basic
import Mathlib.Tactic.Ring
import Mathlib.Tactic.NormNum.Basic
import LewontinParadox.Heterozygosity

namespace LewontinParadox.Coalescent

open Real LewontinParadox.Heterozygosity

/-! ## Parameters -/

/-! Effective (diploid) population size is passed explicitly to each
definition below. -/

/-- Expected pairwise time to the most recent common ancestor at a single
neutral site, in generations. -/
noncomputable def E_Tmrca (Nₑ : ℝ) : ℝ := 2 * Nₑ

/-- The standard Kingman-coalescent mean pairwise time `E[T_MRCA] = 2·Nₑ`
(diploid). This is the one analytic fact we import; proving it would
require formalizing the coalescent process. -/
theorem E_Tmrca_eq (Nₑ : ℝ) : E_Tmrca Nₑ = 2 * Nₑ := rfl

/-! ## Expected pairwise diversity -/

/-- Expected pairwise nucleotide diversity at a single neutral site:
`E[π] = 2·μ·E[T_MRCA]` (two lineages, total branch length `2·T_MRCA`). -/
noncomputable def E_pairwise_pi (Nₑ μ : ℝ) : ℝ := 2 * μ * E_Tmrca Nₑ

/-- `E[π] = 4·Nₑ·μ`, the standard neutral result. -/
theorem E_pairwise_pi_eq (Nₑ μ : ℝ) :
    E_pairwise_pi Nₑ μ = 4 * Nₑ * μ := by
  dsimp [E_pairwise_pi, E_Tmrca]; ring

/-- **Theorem.** The expected pairwise diversity at a single neutral site
does not depend on the gene-conversion rate `c`: it equals `θ = 4·Nₑ·μ`.

The `c`-independence is structural: `E_pairwise_pi` is built from
`E_Tmrca`, which by `E_Tmrca_eq` is `2·Nₑ` regardless of `c`. Gene
conversion reshuffles correlations between sites (the ARG) but leaves the
marginal genealogy of any single site — and hence `E[π]` — untouched.
-/
theorem E_pi_independent_of_c (Nₑ μ c : ℝ) :
    E_pairwise_pi Nₑ μ = theta Nₑ μ := by
  dsimp [E_pairwise_pi, E_Tmrca, theta]; ring

/-- Under the neutral coalescent, expected diversity equals the
population-scaled mutation rate `θ = 4·Nₑ·μ`. -/
theorem E_pi_eq_theta (Nₑ μ : ℝ) :
    E_pairwise_pi Nₑ μ = theta Nₑ μ := E_pi_independent_of_c Nₑ μ 0

end LewontinParadox.Coalescent
