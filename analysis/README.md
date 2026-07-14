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

The population-VCF execution audit is recorded in `TIER3B_RUN.md`.  No surveyed
population tuple currently passes the frozen exact-reference, native-annotation,
sample, and invariant-denominator gates, so `tier3b_data.tsv` is intentionally
header-only.  Treat it as structured missingness and consult
`tier3b_failure_ledger.tsv` and `tier3b_qc_provenance.json`; do not interpret
the absence of rows as zero diversity.

Tier 3a is likewise header-only after its frozen gates. The final synthesis
keeps population π, deposited-call individual heterozygosity, and
alignment-conditioned individual heterozygosity distinct; it does not pool
them into generic π. `tier3_results.tsv` has 270 observable point rows plus
model and claim rows, including null/negative effects and explicit unavailable
PGLS, recombination-class, callable-policy, and SFS-B results.

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
