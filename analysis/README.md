# analysis/ — empirical work for the gBGC-saturation test of Lewontin's Paradox

The formal theory is in `../lean/`; the manuscript is `../manuscript.typ`.
This directory holds the **empirical** Tier-3 work and its supporting plans.

## Plans & surveys

| file | what |
|---|---|
| `TIER3_PLAN.md` | the science: the discriminating cross-species test (composition + W/S-stratified π), the 3a/3b/3c tier design, and the within-vs-across-clade framing |
| `TIER3B_VCF_SURVEY.md` | data availability: ENA WGS run counts + named population-VCF resources (DGRP, Ag1000G, CaeNDR, etc.) per high-Nₑ species |
| `TIER3_EXECUTION.md` | the ops: moving Tier 3 to a cluster (environment, storage budget, SLURM job templates, phased execution, reproducibility). *Cluster to be chosen.* |
| `TIER3_RESULTS.md` | final fail-closed synthesis: exact joins, fits, uncertainty, negative/null results, structured missingness, and claim boundaries |

## Scripts (existing)

| file | what | status |
|---|---|---|
| `extend_buffalo.py` | Tier 1: power-law / saturating / MM fits + recombination sign test on Buffalo's 172 taxa | done, committed |
| `make_fig.py` | Tier 1 figure (`fig_extend_buffalo.{pdf,png}`) | done |
| `proxy_pinps.py` | proxied Tier-3 test: πN/πS vs Nc from Buffalo's bundled Romiguier 2014 table | done, committed (supportive, confounded) |

## Tier 3 scripts

- `tier3_manifest.py` — join Buffalo's 173 core species to VGP / pop-VCF / NCBI accessions
- `tier3c_ncbi_gc.py` — Phase 1: GC3 + GC% from NCBI assemblies (array job)
- `tier3a_vgp_*.py` — Phase 3: GC3 + individual π + W/S-stratified π from VGP phased assemblies
- `tier3b_popvcf_*.py` — Phase 2: population π + W/S-stratified π + (SFS-B) from population VCFs
- `tier3_fit.py` — Phase 4: exact-species merge, observable-separated fits,
  small result table, and headless figure (`fig_tier3.{pdf,png}`); done

## Figures

- `fig_extend_buffalo.{pdf,png}` — Tier 1 (rising-then-plateau + recombination sign)
- `fig_proxy_pinps.{pdf,png}` — proxied πN/πS test (confounded)
- `fig_tier3.{pdf,png}` — final exact-species composition points and the
  explicit diversity/SFS missingness boundary

## Data

Raw data is **not committed** (gitignored). Buffalo's `combined_data.tsv` lives
at `/tmp/paradox_variation/data/` on this machine (not in repo). Cluster-staged
assemblies/VCFs live in a scratch dir with an accession `manifest.tsv`
(reproducible, re-fetchable). Only small result TSVs and figures are committed.

## Final Tier 3 synthesis

The recovered biological inputs are `results/tier3a/diploid_diversity.tsv`
and `results/tier3b/population_diversity.tsv`; the legacy header-only
`tier3a_data.tsv` and `tier3b_data.tsv` are no longer synthesis inputs. Tier 3a
has three corrected origin/main SweepGA + IMPG H1/H2 coding-panel assembly
pairs. Tier 3b has two 20-individual *Anopheles coluzzii* populations. The
final synthesis keeps these alignment-conditioned assembly estimates,
population diversity, and exact-single-assembly composition distinct; it does
not pool them into generic π.

`tier3_results.tsv` now carries all 23 recovered diversity observations with
identity, modality, numerator, callable denominator, explicit eligible `n`,
exclusions, software provenance, and uncertainty, alongside the 270
composition point rows and model/claim rows. Null/negative composition effects
and explicit unavailable PGLS, recombination-class, deposited-individual, and
SFS-B results remain visible. The recovery summary, complete evidence ledger,
independent headline audit, and checksummed artifact index are under
`results/tier3/`.

Run the fits and all analysis tests through the pinned environment:

```sh
guix time-machine -C analysis/guix/channels.scm -- \
  shell -L analysis/guix -m analysis/guix/manifest.scm --pure -- \
  python3 analysis/tier3_fit.py

guix time-machine -C analysis/guix/channels.scm -- \
  shell -L analysis/guix -m analysis/guix/manifest.scm --pure -- \
  python3 -m pytest -q analysis/tests
```

Offline figure regeneration needs only the committed small result table:

```sh
guix time-machine -C analysis/guix/channels.scm -- \
  shell -L analysis/guix -m analysis/guix/manifest.scm --pure -- \
  python3 analysis/tier3_fit.py \
    --figure-from-results analysis/tier3_results.tsv
```
