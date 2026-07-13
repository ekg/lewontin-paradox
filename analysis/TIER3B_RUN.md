# Tier 3b population-VCF run status

Run ID: `tier3b-fail-closed-20260713`

## Outcome

No Tier 3b population job was submitted.  This is the required fail-closed
outcome, not a zero-diversity result.  At execution time there was no frozen
analysis manifest and no checksum-locked population tuple under the declared
Tier 3 scratch root.  The committed pilot registry still labels both DGRP and
Ag1000G ineligible pending exact-reference and invariant-denominator inputs.
The frozen policy requires those pilots to reproduce exactly before any
expansion, so *D. simulans*, *D. pseudoobscura*, *Aedes*, and *Daphnia* were not
scheduled.

The resulting [tier3b_data.tsv](tier3b_data.tsv) contains its schema header and
zero rows.  It must be interpreted as structured missingness.  It does not
contain an SFS-B column because the polarization gate is deferred in
`tier3-decisions-v1`.  [tier3b_failure_ledger.tsv](tier3b_failure_ledger.tsv)
records the six candidate-specific gate failures, and
[tier3b_qc_provenance.json](tier3b_qc_provenance.json) records the full
candidate audit, pinned environment, compute smoke, intended sampling policy,
and absence of resource telemetry for jobs that were never launched.

In particular, an available SNP-only VCF was not treated as evidence for
invariant callability, a generic accessibility mask was not substituted for an
exact selected-cohort callable mask, and a native annotation was not inferred
from species or assembly names.  No raw-read calling, Conda, micromamba,
outgroup download, or projected annotation was used.

## Reproduction

The artifacts were generated through the previously validated, GC-rooted Guix
profile and compute-node store-path record:

```sh
TIER3_SCRATCH_ROOT=/moosefs/erikg/tier3scratch/tier3-run-popvcf \
analysis/slurm/guix_job.sh "$PWD/analysis/pilot_results/guix_environment.json" \
  python3 analysis/tier3b_finalize_run.py \
  --pilot-registry analysis/pilot_results/pilot_registry.json \
  --environment analysis/pilot_results/guix_environment.json \
  --compute-smoke analysis/pilot_results/compute_smoke.json \
  --data-output analysis/tier3b_data.tsv \
  --failure-ledger analysis/tier3b_failure_ledger.tsv \
  --qc-output analysis/tier3b_qc_provenance.json
```

The finalizer verifies that both pilot entries remain explicitly pending, that
the Guix and compute-smoke policy versions match, and that login and compute
nodes used the same profile store path.  It refuses to emit an empty table if
either pilot is promoted, preventing a stale missingness ledger from hiding a
newly runnable dataset.  Repeated finalization is byte-identical.

## Promotion gate

Before any population is added to `tier3b_data.tsv`, a frozen manifest must
supply an immutable deposited VCF/BCF release, exact FASTA and native GFF,
assembly accession including version, file checksums and sizes, explicit
contig mapping, genetic code, provider CDS reconstruction evidence, one
qualifying population, exact sample and exclusion checksums, and an all-sites,
gVCF, or exact selected-cohort callable denominator.  DGRP and Ag1000G must
then independently reproduce their approved numerators and denominators before
the expansion candidates can be submitted.

For *Daphnia pulex*, supplying files alone is insufficient: the currently
surveyed cyclical-parthenogen/temporal design is excluded from primary
population pi by the frozen policy and requires a new, separately approved
population design decision.
