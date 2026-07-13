# analysis/ — empirical work for the gBGC-saturation test of Lewontin's Paradox

The formal theory is in `../lean/`; the manuscript is `../manuscript.typ`.
This directory holds the **empirical** Tier-3 work and its supporting plans.

## Plans & surveys

| file | what |
|---|---|
| `TIER3_PLAN.md` | the science: the discriminating cross-species test (composition + W/S-stratified π), the 3a/3b/3c tier design, and the within-vs-across-clade framing |
| `TIER3B_VCF_SURVEY.md` | data availability: ENA WGS run counts + named population-VCF resources (DGRP, Ag1000G, CaeNDR, etc.) per high-Nₑ species |
| `TIER3_EXECUTION.md` | the ops: moving Tier 3 to a cluster (environment, storage budget, SLURM job templates, phased execution, reproducibility). *Cluster to be chosen.* |

## Scripts (existing)

| file | what | status |
|---|---|---|
| `extend_buffalo.py` | Tier 1: power-law / saturating / MM fits + recombination sign test on Buffalo's 172 taxa | done, committed |
| `make_fig.py` | Tier 1 figure (`fig_extend_buffalo.{pdf,png}`) | done |
| `proxy_pinps.py` | proxied Tier-3 test: πN/πS vs Nc from Buffalo's bundled Romiguier 2014 table | done, committed (supportive, confounded) |

## Scripts (to be written, on the cluster — see `TIER3_EXECUTION.md`)

- `tier3_manifest.py` — join Buffalo's 173 core species to VGP / pop-VCF / NCBI accessions
- `tier3c_ncbi_gc.py` — Phase 1: GC3 + GC% from NCBI assemblies (array job)
- `tier3a_vgp_*.py` — Phase 3: GC3 + individual π + W/S-stratified π from VGP phased assemblies
- `tier3b_popvcf_*.py` — Phase 2: population π + W/S-stratified π + (SFS-B) from population VCFs
- `tier3_fit.py` — Phase 4: merge tiers, fit, figure (`fig_tier3.{pdf,png}`)

## Figures

- `fig_extend_buffalo.{pdf,png}` — Tier 1 (rising-then-plateau + recombination sign)
- `fig_proxy_pinps.{pdf,png}` — proxied πN/πS test (confounded)

## Data

Raw data is **not committed** (gitignored). Buffalo's `combined_data.tsv` lives
at `/tmp/paradox_variation/data/` on this machine (not in repo). Cluster-staged
assemblies/VCFs live in a scratch dir with an accession `manifest.tsv`
(reproducible, re-fetchable). Only small result TSVs and figures are committed.

## Current Tier 3b execution

The population-VCF execution audit is recorded in `TIER3B_RUN.md`.  No surveyed
population tuple currently passes the frozen exact-reference, native-annotation,
sample, and invariant-denominator gates, so `tier3b_data.tsv` is intentionally
header-only.  Treat it as structured missingness and consult
`tier3b_failure_ledger.tsv` and `tier3b_qc_provenance.json`; do not interpret
the absence of rows as zero diversity.
