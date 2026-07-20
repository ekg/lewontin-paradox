# VGP ten-pair Slurm entry points

These scripts never submit by default. `submit.sh --dry-run` prints the exact
dependency-linked commands; `--submit` is an explicit operational action for
the separately authorized run task. Every pair requires an acquired
`input-manifest.json`, a measured `resources.json`, a captured pinned Guix
environment, and stage success sentinels.

The pair path is `pilot/runs/<run_id>/<selection_id>`. `pair_stage.sh` owns the
preflight, whole-assembly SweepGA mapping, exact IMPG index/partition/query/lace
sequence, and bcftools normalization stages. `psmc_array.sh` reserves task zero
for the unscaled primary run and tasks 1–200 for job-local inputs generated
by resampling contig-bounded blocks of the primary PSMCFA itself. The frozen
unit/sample manifests therefore retain the primary population of masked `N`
and callable `K/T` bins instead of sampling callable BED fragments. `psmc_finalize.sh`
requires the full array, checks the 190/200 finite-output gate, and writes
unscaled and scenario-scaled tables separately. Consensus/mask/bootstrap materialization is intentionally a
pair-level join: it must be produced only after REF, H2 reconstruction, and
reason-mask tests pass.

No script embeds a global memory, byte, or wall ceiling. The submitter reads
job-specific limits derived from actual FASTA byte counts, sequence lengths,
contig counts, native partition counts, and prior telemetry. Node jobs use a
pre-realized `/gnu/store` profile and do not access the network.
