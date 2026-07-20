# analysis/ — empirical work for the gBGC-saturation test of Lewontin's Paradox

The formal theory is in `../lean/`; the manuscript is `../manuscript.typ`.
This directory holds the **empirical** Tier-3 work and its supporting plans.

## VGP integration and authorization lineage

The comprehensive VGP branch/commit audit is recorded in
`vgp_comprehensive_integration_audit.md` and its machine-readable companion
`vgp_comprehensive_integration_manifest.json`. Run
`python3 analysis/validate_vgp_comprehensive_integration.py` to verify the
exact 10-primary/6-alternate rosters, required SweepGA/IMPG/PSMC/Slurm/mirror
tooling, integration order, and historical artifact digests.

Prior `NO_GO`, `CONDITIONAL_GO`, zero-job, zero-core-pass, and zero-estimate
outputs are timestamped in `vgp_historical_run_registry.json`. They are
immutable historical evidence, not current execution authorization and not an
inheritable authorization policy. A new run must bind a separate versioned
current authorization packet to its immutable inputs and scientific preflight
gates.

The current packet is `vgp_pilot_authorization_v2.json`, with live read-only
preflight in `vgp_pilot_authorization_preflight_v2.json` and the exact P07
128/256/512-GiB canary packets under
`slurm/vgp_10_pilot/authorized/v2.0.0/`. It authorizes all ten frozen primary
pairs for core diversity and unscaled PSMC. Missing optional QC lowers
confidence; it does not reduce the authorized pair count. See
`vgp_pilot_authorization_handoff_v2.md` for the exact execution boundary and
handoff command.

All active VGP CAS, mirror, staging, raw/validation, derived-output, and Slurm
paths resolve under the canonical shared root `/moosefs/erikg/vgp`. The former
project-named root is migration input only; verified objects were hard-linked
into the canonical CAS without redownload as recorded in
`vgp_data_root_migration_v1.json`.

### Completed real P07 canary

The authorized P07 *Spinachia spinachia* canary completed on real Slurm jobs
through the pinned Guix environment. See `vgp_real_canary_report_v1.md` for
the scientific and operational report and `vgp_real_canary_execution_v1.json`
for the machine-readable independent audit. The small independent variant
subset, authoritative allocation telemetry, and atomic-promotion receipt are
`vgp_real_canary_variant_subset_v1.tsv`, `vgp_real_canary_sacct_v1.tsv`, and
`vgp_real_canary_promotion_v1.json`. The promoted biological payload remains
under `/moosefs/erikg/vgp/pilot/outputs/vgp10-auth-20260718-v2/P07/`; it is not
duplicated in Git.

### Real ten-pair pilot and P04 completion

The remaining nine authorized primaries were submitted as independent,
resumable Slurm chains. P04 *Falco naumanni* is the first new completed pair:
770,780,965 final callable bp, 3,548,818 heterozygous SNPs,
π = 0.004604184795871289, an unscaled PSMC trajectory, and 200 finite
bootstraps. See `vgp_real_pilot_report_v1.md` and the machine-readable
`vgp_real_pilot_closed_world_v1.json`. The detailed independent P04 audit,
early/middle/late variant subset, complete Slurm census, and live FastGA
scratch audit are `vgp_real_pilot_P04_execution_v1.json`,
`vgp_real_pilot_P04_variant_subset_v1.tsv`, `vgp_real_pilot_sacct_v1.tsv`,
and `vgp_real_pilot_fastga_scratch_v1.json`.

All pilot products remain under `/moosefs/erikg/vgp`. Hard primary input
failures and retryable infrastructure attempts remain explicit in the
canonical failure ledger; no alternate was activated. Every current SweepGA/
FastGA retry stages both FASTAs and indexes in a private `/scratch` tree,
exports `TMPDIR`, `TMP`, and `TEMP`, enters the tree before alignment, and
fails closed if live `/proc` cwd or managed pair/index/intermediate paths
escape node-local scratch.

### Independent validation reads

The independent raw-read validation subset is frozen in
`vgp_validation_read_plan_v1.json`: P07 is the small later-generation fish,
P09 is the very large/repeat-sensitive shark, and P04 is the early CLR/
TrioCanu bird. Every object repeats the exact BioSample, individual, run,
experiment, platform/library, and both assembly accession.versions. ENA MD5
and byte counts are verified before a local SHA-256 CAS promotion; the
post-promotion object is rehashed and exposed read-only under the configured
accession view (or `raw_reads` view when a future contract defines one).
`vgp_validation_reads_manifest_v1.json` accounts for every
planned, transferred, verified, reused, missing, pending, and quarantined
object and byte. It is also atomically copied to the configured canonical
`pilot_manifests` directory so downstream projects do not depend on this
checkout.

Acquisition is intentionally separate from the assembly-derived pilot run and
is safe to resume. On this cluster, use the host Python so HTTPS uses the site
trust configuration (the pure Guix Python correctly refuses the site's
self-signed TLS chain):

```sh
python3 analysis/acquire_vgp_validation_reads.py
```

Interrupted payloads remain only in the configured canonical staging tree;
content mismatches move atomically to the configured quarantine tree and are
never exposed through a read view.

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
- `vgp_10_pilot.py` — strict ten-pair H1/H2 SweepGA→IMPG→mask/consensus→PSMC
  workflow accounting; see `vgp_10_pilot_workflow_handoff.md`, the pinned
  `guix/vgp_10_pilot/` environment, and resumable `slurm/vgp_10_pilot/` entry
  points. Its deterministic fixtures prove software behavior only.

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
