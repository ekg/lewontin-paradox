# Tier 3 diversity recovery task quality review

Review task: `.quality-pass-tier3-diversity-recovery`
Review date: 2026-07-16

This is a coordination-only review of WG task definitions. It contains no
biological acquisition, computation, or project-result analysis.

## Reviewed task IDs and edits

1. `acquire-tier3b-population-tuples`
   - Made the public-source order explicit: DGRP, Ag1000G, then ranked
     authoritative high-Ne alternatives.
   - Required continued pursuit after unavailable or disqualified resources,
     with attempted locations and rejection reasons recorded.
   - Strengthened exact-reference, original native-annotation, callability,
     checksum, Guix, nonempty-input, and exact acquisition-artifact contracts.
   - Restricted tracked writes to `results/tier3b/acquisition_*` and made prior
     results/workflows read-only.

2. `acquire-tier3a-vgp-tuples`
   - Added an explicit ranked sequence from authoritative VGP releases to other
     public phased diploid alternatives.
   - Required continued pursuit after the prior local failure or any unavailable
     candidate and retained a mandatory minimum of one biological tuple.
   - Strengthened exact H1/H2 reference, original H1-native annotation,
     haplotype, callability, checksum, Guix, nonempty-input, and exact
     acquisition-artifact contracts.
   - Restricted tracked writes to `results/tier3a/acquisition_*` and made prior
     results/workflows read-only.

3. `run-tier3b-biological-recovery`
   - Made acquisition artifacts read-only and assigned run outputs to
     `results/tier3b/population_*` and `results/tier3b/run_logs/`.
   - Required at least two biological tuple results and at least one non-header
     row per accepted tuple; `n=0`, all-failed, zero-denominator, missing-
     uncertainty, and synthetic outputs are explicit failures.
   - Required per-estimate numerator, callable denominator, eligible sample
     size, uncertainty method/unit/replicates/interval, and a programmatic
     nonempty-output assertion.
   - Made all environments Guix-only and enumerated exact run deliverables.

4. `run-tier3a-biological-recovery`
   - Made acquisition artifacts read-only and assigned run outputs to
     `results/tier3a/diploid_*` and `results/tier3a/run_logs/`.
   - Required at least one non-synthetic biological result with a non-header
     data row; `n=0`, all-failed, zero-denominator, missing-uncertainty, and
     synthetic-only outputs are explicit failures.
   - Required per-estimate numerator, callable/aligned denominator, eligible
     sample size, uncertainty method/unit/replicates/interval, and a
     programmatic nonempty-output assertion.
   - Made all environments Guix-only and enumerated exact run deliverables.

5. `synthesize-tier3-diversity-recovery`
   - Made the two-parent join explicit: synthesis waits for both biological run
     tasks to be Done and for both QC gates to pass.
   - Named both result tables and both QC reports as required read-only inputs.
   - Added pre-synthesis assertions for biological rows, positive callable
     denominators, finite estimates, eligible `n`, and uncertainty.
   - Assigned new outputs to `results/tier3/recovery_*` and named four exact
     synthesis artifacts, including an index of every changed manuscript,
     figure, table, command, environment, and supporting result path.

## Validation outcome

- Both acquisition tasks require pursuit of ranked public alternatives and
  prohibit stopping at the first unavailable resource.
- Exact-reference, original native-annotation, haplotype/sample, callability,
  checksum, and coordinate provenance survive acquisition and downstream use.
- GNU Guix is mandatory for every acquisition, analysis, and synthesis
  environment.
- Both biological run tasks require nonempty biological result rows, positive
  callable denominators, eligible sample sizes, explicit uncertainty, and
  programmatic failure on header-only or `n=0` outputs.
- Acquisition and run writers are sequential and have disjoint filename
  prefixes; synthesis treats both upstream result trees as read-only and owns
  only the `results/tier3/recovery_*` namespace.
- `synthesize-tier3-diversity-recovery` directly depends on both
  `run-tier3a-biological-recovery` and `run-tier3b-biological-recovery`.
- The synthesis inputs, QC preconditions, manuscript-update requirement, and
  four exact machine-readable/document deliverables are unambiguous.

All quality-pass checks passed after the task-definition repairs above.
