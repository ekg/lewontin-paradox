# Tier 3a VGP execution audit

Status: **no eligible VGP individual tuple**, as of the pinned environment
record dated 2026-07-13.

The frozen pilot registry does not identify an individual, assembly pair, or
checksum-locked deposited-call/callable-mask/native-annotation tuple.  Its sole
Tier 3a pilot remains explicitly `ineligible_pending`.  Consequently no VGP
array job was submitted and `tier3a_data.tsv` is intentionally header-only.
This is structured missingness, not an estimate of zero heterozygosity.

`tier3a_finalize_run.py` checks the registry against the audited Guix profile
and real compute-node smoke before producing three small committed artifacts:

- `tier3a_data.tsv` defines the complete result schema, including the exact H1
  reference tuple, unique callable H1 bases, total/4D-W/4D-S numerators and
  denominators, bootstrap uncertainty, annotation provenance, and telemetry;
- `tier3a_failure_ledger.tsv` records every considered candidate and its
  reason-coded preflight exclusions;
- `tier3a_qc_provenance.json` records the native-annotation, phasing/collapse,
  deposited-callability, direct-WFMASH, dual-modality, IMPG, Guix, and compute
  smoke gates.

The local inventory contains a primary/alternate *Myotis daubentonii* pair
(`GCF_963259705.1` / `GCA_963242275.1`), but it was not promoted.  It is absent
from the frozen eligible manifest, has no exact Buffalo species covariate, no
frozen deposited-call plus callable-mask tuple, no staged native-H1 GFF/CDS
audit, and no frozen phase-identity or collapse/duplication audit.  Existing
legacy alignments therefore cannot substitute for a pinned WFMASH run or an
explicit denominator.

The QC record identifies the audited VGP Phase 1 freeze and Buffalo core table
by SHA-256.  Those off-repository inventories were used only to establish
missingness; the deterministic finalizer does not silently reread mutable
external paths, and it labels the stored inventory hashes as not reverified on
later regeneration.

The compute smoke did establish that the shared GC-rooted profile uses WFMASH
commit `e040aa10e87cab44ed5a4db005e784be62b0bd21`, emits extended CIGAR, handles
unique mappings, builds BCF/CSI, and reproduces the 1/100 denominator truth.
That validates execution readiness, not eligibility of a biological row.
IMPG remains source-only and execution-disabled; its real-candidate boundary,
deduplication, indexed-query, and phase-orientation gates were not claimed.

Regenerate the committed artifacts from the repository root with the pinned
pure Guix job wrapper:

```sh
analysis/slurm/guix_job.sh analysis/pilot_results/guix_environment.json \
  python3 analysis/tier3a_finalize_run.py \
  --pilot-registry analysis/pilot_results/pilot_registry.json \
  --environment analysis/pilot_results/guix_environment.json \
  --compute-smoke analysis/pilot_results/compute_smoke.json \
  --data-output analysis/tier3a_data.tsv \
  --failure-ledger analysis/tier3a_failure_ledger.tsv \
  --qc-output analysis/tier3a_qc_provenance.json
```

If the pilot registry is later promoted, the finalizer rejects the empty-table
path.  The frozen tuple must then run through the Tier 3a Slurm workflow and
pass deposited-versus-direct common-callable concordance before expansion.
