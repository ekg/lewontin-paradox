# Tier 3 Slurm execution

These scripts target the frozen `linux` Slurm cluster and `workers` partition.
They use a two-part Guix topology because the compute-node daemon socket is not
usable:

1. `prepare_guix.sh ABSOLUTE_STATE_DIR` runs on `octopus01`. It realizes
   `analysis/guix/manifest.scm` through the committed time-machine channels,
   registers the profile derivation as a persistent GC root, and records the
   resolved channels, derivations, complete closure, executable store paths,
   and versions in `environment.json`.
2. `guix_job.sh` runs on compute nodes. It checks the channel/manifest hashes
   and exact shared-store profile from that record, constructs a clean
   environment, and directly executes profile programs. It never asks the
   unavailable compute-node daemon to build.

`submit_tier3.sh ABSOLUTE_WORKFLOW ABSOLUTE_STATE_DIR` performs realization,
validates every dataset before scheduling expensive work, enforces the frozen
500-GiB free-space gate, submits the compute smoke, and makes each non-empty
array depend on the smoke succeeding. Individual scripts are also accepted by
`sbatch --test-only`.

## Workflow contract

A workflow is JSON with `schema_version: "1.0"`,
`decision_version: "tier3-decisions-v1"`, absolute `scratch_root` and
`output_root`, and a `datasets` list. Each dataset declares:

- `dataset_id`, `tier`, `mode`, `reference_accession`, and optional `pilot`;
- `preflight` paths and provenance. Diversity modes require the applicable
  exact-reference callable mask; direct WFMASH requires H1, H2, extended-CIGAR
  PAF, phase identity QC, collapse QC, and native H1 annotation provenance;
- an ordered `stages` list. A stage contains `name`, an argv list, and relative
  non-empty `outputs`. No command is interpreted by a shell.

The exact stage sequences are:

| mode | atomic stages |
|---|---|
| `composition` | `annotation_4d`, `qc` |
| `population_vcf` | `mapping`, `normalized_bcf`, `annotation_4d`, `qc` |
| `deposited_vcf` | `mapping`, `normalized_bcf`, `annotation_4d`, `qc` |
| `direct_wfmash` | `alignment`, `mapping`, `normalized_bcf`, `annotation_4d`, `qc` |

Arguments may use `{stage_dir}`, `{dataset_dir}`, `{repo_root}`,
`{workflow_dir}`, `{input:FIELD}`, or `{stage:NAME}`. Environment expansion is
limited to `${SLURM_*}`, `${TIER3_*}`, `${SCRATCH}`, and `${TMPDIR}`. A stage
is reusable only when its fingerprint and every declared output checksum still
match. Abandoned `.partial-*` directories are removed on resubmission; a
successful stage is atomically renamed and never inferred from file presence
alone.

Native annotation provenance must include provider, release, exact assembly
accession+version, FASTA/GFF SHA-256, explicit contig mapping, genetic code,
and native/projected status. The runner verifies sequence-region dictionaries
and reconstructs every retained CDS against the exact FASTA before any primary
annotation/4D stage. Projected annotation is never accepted as primary.

## Optional IMPG

IMPG remains execution-disabled in the current manifest. If a later frozen
policy enables it, its stage must retain the existing implementation contract:
padded queries with disjoint ownership cores, `vcf:poa`, merge distance zero,
`lace` with at least two threads, global FASTA normalization, exact allele
deduplication, normalized `POS-1` ownership exactly once, BCF/CSI, and an
explicit phase-orientation audit. It is concordance QC only and never supplies
the callable denominator.

A Guix pack is not used while the shared store works. If a future compute smoke
cannot execute the recorded profile, the only permitted compatibility path is
a time-machine-built Guix pack whose manifest checksum, runtime store path,
pack path, and pack SHA-256 are all added to the environment record before use.

## Production Tier 3c freeze and run

NCBI discovery is not itself a frozen analysis input.  Production composition
therefore uses `analysis/tier3c_batch.py` in two scheduler phases:

1. `discover` joins the checksum-verified Buffalo cohort to the deterministic
   same-species NCBI ranking policy.  It never substitutes a congener.
2. `tier3c_stage_array.sh` acquires each selected accession, verifies NCBI's
   provider MD5s, excludes only assembly-report-labelled non-nuclear units,
   BGZF-compresses the exact derived FASTA/GFF, and records SHA-256, byte size,
   native-annotation provenance, Guix closure, RAM, storage, and walltime.
3. `freeze` verifies every staged checksum and creates the only batch manifest
   accepted by `tier3c_run_batch_array.sh`.  Missing stage records remain
   explicit pre-analysis failures rather than silently shrinking the cohort.
4. `collect` writes `analysis/tier3c_manifest.tsv`, `analysis/tier3c_data.tsv`,
   `analysis/tier3c_failure_ledger.tsv`, `analysis/tier3c_qc_summary.json`, and
   per-species `analysis/tier3c_qc/*.json`.  Raw assemblies and scheduler logs
   remain in the external scratch/output roots.

The analysis result JSON is independent of scheduler walltime and is left
untouched by an identical rerun.  Per-invocation telemetry is written to a
separate `.job.json` sidecar, so byte-identity checks compare scientific result
bytes without discarding operational provenance.  Primary GC3 is emitted only
when the NCBI GFF is native to the identical accession/version and all retained
CDS reconstruct against that filtered exact FASTA; missing GFF yields explicit
GC3 unavailability while nuclear whole-genome GC remains eligible.

Production analysis uses two non-constrained resource lanes.  Standard tasks
run through `tier3c_run_batch_array.sh` at 2 CPUs, 32 GiB, 2 hours, and a Slurm
array throttle of at most eight tasks.  Observed or predicted memory/time
outliers use `tier3c_run_outlier_retry.sh` at 2 CPUs, 64 GiB, 4 hours, and a
single-task throttle; its final argument is a newline-delimited list of frozen
batch indices.  Requested capacity is scheduler headroom, not a scientific QC
criterion.

The post-run control audit uses `tier3c_control_audit_array.sh` at the standard
profile.  It runs a separate FASTA/GFF3 parser that does not import the
production estimator, reconstructs both exact native control annotations, and
reports pooled-CDS-base and equal-gene-weight GC3 separately.  Its output is
combined and reviewed before any promotion decision; a larger memory request
cannot turn a control mismatch into a scientific pass.
