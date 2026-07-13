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
