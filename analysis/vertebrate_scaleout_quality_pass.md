# Vertebrate scale-out planning quality pass

Date: 2026-07-17 (UTC)
WG task: `.quality-pass-vertebrate-scaleout-plan`

## Scope and disposition

This quality pass reviewed the definitions and graph relationships for all three
tasks in the vertebrate scale-out planning batch:

1. `plan-all-vertebrate-tier3`
2. `plan-vertebrate-ne-strategy`
3. `synthesize-vertebrate-scaleout-plan`

The definitions were tightened before planning dispatch. This pass authorizes
planning only. It does not authorize bulk downloads, Slurm submission, a full
vertebrate run, or expensive demographic inference.

At review time, `sync-origin-scaleout-guidance` was
`failed-pending-eval` and had not registered its inventory artifact. Both
planning tasks therefore now have that synchronization task as a direct
dependency, as well as the quality-pass dependency. This prevents a planner
from running without the synchronized guidance even if the quality task is
released early by graph-state recovery.

## Edits reviewed and applied

### `plan-all-vertebrate-tier3`

- Named `analysis/vertebrate_scaleout_origin_inventory.md` and its identified
  guidance file as controlling inputs; missing guidance must be reported as a
  blocker rather than replaced with assumptions.
- Added a row-level eligibility inventory with source URLs, retrieval dates,
  accessions, versions, checksums, reference relationships, provenance,
  uncertainty, eligibility tier, and rejection reason.
- Made exact H1/H2 reference identity, original H1-native annotation
  compatibility, and callable/queryable denominators explicit execution gates.
  Annotation lift-over or substitution to another reference is prohibited.
- Preserved validated SweepGA and IMPG responsibilities and required pinned GNU
  Guix channels, manifests, derivation hashes, build commands, and runtime
  capture. Ambient package or tool substitutions cannot enter production
  lineage without a separate validated plan change.
- Required per-stage and per-species extrapolation of core-hours, memory, wall
  time, scratch, persistent storage, temporary amplification, inode counts,
  MooseFS reads/writes, metadata operations, and peak bandwidth. Estimates must
  include formulas, calibration rows, low/base/high bounds, catalog totals,
  concurrency constraints, and retention policy.
- Expanded the Slurm/I/O design to require staged read-only inputs, local
  scratch, bounded arrays, atomic output promotion, checksums, idempotent
  sentinels, dependencies, retries, `sacct` telemetry, throttling, and
  failed-element recovery.
- Prespecified outputs, sampling units, exclusions, missingness, uncertainty,
  phylogenetic covariance, clade balance, multiple testing, sensitivity
  analyses, and descriptive-versus-causal claim boundaries.
- Required inventory/preflight, stratified pilot, reviewed expansion waves, and
  full-catalog phases with quantitative go/no-go conditions and separate
  authorization before expensive execution.
- Added `sync-origin-scaleout-guidance` as a direct dependency.

Manual result: the definition contains the required eligibility inventory,
resource extrapolation, Slurm/I/O design, provenance, phylogenetic/statistical
design, and pilot/full phases.

### `plan-vertebrate-ne-strategy`

- Separated method contracts rather than treating sequentially Markovian
  coalescent methods as interchangeable:
  - PSMC requires one callable diploid genome from one biological individual,
    diploid consensus/heterozygous sites, and a matched callability mask; its
    primary output is a pairwise coalescent-rate/history in scaled units before
    external calibration.
  - MSMC2 requires multiple accurately phased haplotypes and matched masks; it
    estimates within/cross-coalescence-rate histories and related scaled
    summaries, subject to phasing and time-resolution limits. Legacy MSMC is
    not silently substituted for MSMC2.
  - SMC++ requires multiple same-population diploid samples represented as
    matched genotype data and masks, with sample/population assignments; its
    SFS/linkage model estimates population-size histories or selected split
    models. Phasing or polarization is required only when the selected mode
    actually requires it.
- Required a data-availability audit that distinguishes a phased assembly, raw
  reads, a callable diploid genotype, multiple phased genomes, and a
  multi-sample population VCF. Deposited H1/H2 assemblies cannot be assumed to
  be method-valid diploid genotype input without same-individual, reference,
  heterozygosity, mask, and input-generation evidence.
- Made the circularity prohibition operational: an Ne computed from the same
  heterozygosity, SFS, or pi response—including `pi/(4*mu)`—cannot be treated as
  an independent predictor of that response. PSMC/MSMC2/SMC++ histories remain
  eligible only for clearly labeled descriptive, temporal-shape, validation,
  decomposition, or sensitivity uses that acknowledge shared data.
- Separated coalescent-scaled quantities from absolute Ne and calendar time.
  Absolute conversion requires an explicit convention/ploidy, mutation rate,
  generation time, relevant recombination assumptions, calibration source,
  taxonomic transfer rule, and propagated uncertainty or scenarios.
- Defined independent alternatives and covariates, including independent
  temporal/LD/pedigree estimates and ecological/life-history measurements.
- Required literature/ecological provenance down to estimand definition,
  population/geography, measurement date, units, uncertainty, sample size,
  source, curator, taxonomic match distance, transfer rationale, and quality
  flag. Contemporary, variance, inbreeding, linkage, long-term coalescent, and
  census quantities cannot be collapsed into one undifferentiated Ne field.
- Required a tiered data model, phylogenetically aware primary analysis,
  method-specific availability and resource estimates, pinned GNU Guix
  environments, truth-case pilots, and quantitative go/no-go review.
- Added `sync-origin-scaleout-guidance` as a direct dependency.

Manual result: the definition compares viable methods, audits data
availability, controls circularity and absolute scaling, records uncertainty,
and centers independent covariates in the primary analysis.

### `synthesize-vertebrate-scaleout-plan`

- Confirmed both planners are direct dependencies and required both completed
  artifact sets to be reconciled with the synchronized guidance.
- Required an explicit crosswalk covering guidance, eligibility,
  exact-reference/native-annotation/callability gates, GNU Guix, provenance,
  Slurm/I/O, resource budgets, statistical design, and demographic tiers.
- Separated assembly/composition, population-data, coalescent-history, and
  independent-predictor analysis strata and constrained the claims supported by
  each.
- Required nine staged gates: guidance/schema freeze, eligibility preflight,
  dry-run validation, bounded pilot, independent review, expansion waves,
  full-catalog decision, demographic subset, and final synthesis/release
  review. (The final two execution streams may be ordered or joined explicitly
  in the proposed graph, but neither is silently authorized.)
- Each transition now requires quantitative evidence, stop/rollback criteria,
  a reviewer, an approved resource/download envelope, and a resulting WG state.
- Required low/base/high totals for compute, memory, time, scratch, storage,
  inode count, MooseFS I/O and metadata operations, downloads, bandwidth, and
  concurrency, plus traceability to pilot telemetry.
- Required explicit user decision tasks before any over-budget pilot,
  expansion wave, full-catalog submission, large download, or expensive
  demographic inference. Generated scripts remain inert; plan completion is
  not execution approval, and no executable task may be ready across an
  unapproved boundary.
- Required an executive decision table whose states are bounded proceed,
  explicit authorization required, or deferred.

Manual result: the synthesis joins both plans and contains a structural,
explicit execution-authorization boundary with staged go/no-go decisions.

## Validation record

- All three requested task IDs were reviewed, and every edit is listed above
  and logged on `.quality-pass-vertebrate-scaleout-plan`.
- The scale-out planner covers eligibility, resources, Slurm/I/O, provenance,
  phylogenetic/statistical design, and inventory/pilot/expansion/full phases.
- The Ne planner distinguishes PSMC/MSMC2/SMC++ inputs and estimands; compares
  other viable methods and data coverage; separates scaled and absolute
  quantities; addresses circularity; and defines independent, provenance-rich
  covariates.
- The synthesis has direct dependencies on both planners and an explicit,
  non-ready authorization boundary before expensive execution.
- Automated graph assertions are recorded in the WG task log after execution.
- No compute job or bulk download was launched by this quality pass.
