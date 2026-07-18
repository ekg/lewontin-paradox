# VGP pilot lineage integration ledger

Date audited: 2026-07-18 UTC
Task: `integrate-vgp-pilot`

## Integration decision

The bounded VGP pilot lineage was already integrated on live `main` and
`origin/main` when this audit began. Both refs resolved to
`9480b64f848d168ce3e44dbf27924589b9a0170d`. The integrated history contains
one patch-equivalent replay of every required task-owned commit in dependency
order. Replaying or cherry-picking the task branches again would duplicate
their changes, so no lineage commit was replayed.

This document is the only new project artifact from the integration audit. The
original task branches and commit objects were treated as immutable evidence.
No reset, rebase, clean, force push, amend, broad staging, biological download,
Slurm submission, or Ne inference was performed.

## Live-ref and user-change audit

- Audit-start `HEAD`, local `main`, and `origin/main`:
  `9480b64f848d168ce3e44dbf27924589b9a0170d`.
- The root `main` worktree was clean except for the pre-existing untracked
  `.wg-worktrees/` directory. It was preserved and never staged.
- The isolated `wg/agent-212/integrate-vgp-pilot` worktree was clean before
  this ledger was created and had no commits ahead of `main`.
- A live `git ls-remote --heads origin` query found no surviving remote refs
  for the eight required task branches. WG logs record that each branch was
  pushed before integration, and the local WG worktree refs and exact commit
  objects remain available. Their local branch tips are enumerated below.
- The frozen external manifest is the only biological-scale source object
  found beneath the authorized external root. The immutable-object, staging,
  and pilot-output trees each contained zero files at audit time.
- Git contains only the four pre-existing, tiny synthetic truth fixtures
  (`truth.fa`, `truth.gff3`, `truth.vcf`, and `truth.normalized.bcf`) matching
  biological file suffixes; no VGP payload is tracked.

## Controlling scale-out ancestry

The live first-parent ancestry from the synchronized Tier 3 recovery handoff
through the VGP synthesis is linear. The plan tasks that were logically
parallel were serialized only to make a safe Git history; their WG dependency
edges remain satisfied.

```text
dec7266d3427987d04f65bb1300032c5f233cb95  integrated Tier 3 recovery handoff
  -> 14fe57d1d0593a2859672279be60aa858e822410  frozen-source guidance root
  -> 234b9ca97d63964d130bbe0f9c0ca1bb6eb5d51f  rerun-sync-origin
  -> d586f9e8d33bb9583a2aa94ddeca3d7a7efd53fc  scale-out quality audit replay
  -> ee383af1c9c65cc563afc695d9268397d3ba0284  synchronized quality handoff
  -> e0bfd66e6452a007b1d90d7c613b149b8f868cb6  plan-vertebrate-ne-strategy
  -> 11ee6657637ad48065f189c27e7c5e2e2b7e3d40  plan-all-vertebrate-tier3
  -> 0337dabaf2d3e4288c37dc747cac26880d467661  synthesize-vertebrate-scaleout-plan
  -> 9e3322e715fe739388850f7eadd4c0d7f8298234  integrate-vertebrate-scaleout
  -> 04183df47bd7587612002ced63c332c2c6bcdb11  quality-vgp-pilot
  -> cfa6f34209ad7cd4c2ab9364f6dfdad570f3f3e0  prepare-vgp-data
  -> 4fe9d3f3c0c652f1aadd73251b4347a62b4454e9  freeze-vgp-manifest
  -> b084df1cdc1eef328bec67a2f35741bb2404d750  inventory-pilot-independent
  -> 003a757323bbe592fe5d044ae4691db58b9b06ce  gate-vgp-pilot
  -> 41d8aa1418e42e457822ebbd83fb5a7312e75309  acquire-vgp-pilot
  -> 2b12df1d76d4882b21c0af8ee7a3463dc30421d2  run-vgp-tier3
  -> cbac04e234dad93b2acd07e5f251c37da404a133  review-vgp-pilot
  -> 9480b64f848d168ce3e44dbf27924589b9a0170d  synthesize-vgp-pilot
```

WG declares both `inventory-pilot-independent` and `gate-vgp-pilot` after
`freeze-vgp-manifest`. The synthesis depends on both the inventory and the
review, while acquisition and run depend only on the gate path. Serializing the
inventory replay before the gate replay is therefore dependency-correct and
preserves both disjoint results. No required task-owned patch occurs more than
once on `main`.

## Scale-out task commit equivalence

Stable patch IDs were computed from the complete task-owned and integrated
commit patches. Matching IDs prove that the integrated replay contains the
same textual changes even though its commit and parent IDs differ.

| Task | Task-owned commit | Integrated commit | Stable patch ID | Task-owned artifacts |
| --- | --- | --- | --- | --- |
| `rerun-sync-origin` | `00f72fc86150404bba0712b2672c320077c7b962` | `234b9ca97d63964d130bbe0f9c0ca1bb6eb5d51f` | `c509d70f948de4b6064b222f06a09377f9cc79d5` | `analysis/vertebrate_scaleout_origin_inventory.md` |
| `.quality-pass-vertebrate-scaleout-plan` | `70c17e73e0874151cf13ed935c7daef6e96d46e1` | `d586f9e8d33bb9583a2aa94ddeca3d7a7efd53fc` | `9374def4a6746ff92cbd19a1953efe7e2d45655d` | `analysis/vertebrate_scaleout_quality_pass.md` |
| `.quality-pass-vertebrate-scaleout-plan` synchronized handoff | task follow-up on live main | `ee383af1c9c65cc563afc695d9268397d3ba0284` | not a task-branch replay | updates `analysis/vertebrate_scaleout_quality_pass.md` with the synchronized inventory handoff |
| `plan-vertebrate-ne-strategy` | `f728e73cb434d2e68723913ed5c4d4989b744458` | `e0bfd66e6452a007b1d90d7c613b149b8f868cb6` | `5f5189f9ca92b028bfe75976d62d988b08452dbb` | `analysis/vertebrate_ne_strategy.md`; `analysis/vertebrate_ne_method_matrix.tsv`; `analysis/vertebrate_ne_source_schema.tsv`; `analysis/vertebrate_ne_execution_graph.tsv` |
| `plan-all-vertebrate-tier3` | `1f4e86a870fb79e9a922e622ea67e26f8ca4a714` | `11ee6657637ad48065f189c27e7c5e2e2b7e3d40` | `e1c27b61be50b2a73d79d490cf4aebba1ad976e7` | `analysis/vertebrate_scaleout_plan.md`; candidate schema; resource budget; WG graph |
| `synthesize-vertebrate-scaleout-plan` | `7f10182045ec4c45012b275a4d17eb21caa94436` | `0337dabaf2d3e4288c37dc747cac26880d467661` | `ededb2aee1723df14cac8edf56b83c46c962de13` | `analysis/vertebrate_scaleout_execution_plan.md`; execution graph; decisions table |
| `integrate-vertebrate-scaleout` | `9e3322e715fe739388850f7eadd4c0d7f8298234` | same commit | same object | `analysis/vertebrate_scaleout_integration_handoff.md` |
| `quality-vgp-pilot` | `53525c3be2d70504aa9d0e2adde6ac8d2f343a9d` | `04183df47bd7587612002ced63c332c2c6bcdb11` | `854160301205f32a95d33bd958b33f5c1a461e95` | `analysis/vgp_pilot_execution_quality.md` |

The controlling artifacts establish the external-root boundary, at-most-six
pilot boundary, exact-reference/native-annotation rules, resource caps,
fail-closed gate, and separation of independent Ne evidence from circular
genomic quantities.

## Required VGP task ledger

Each row lists the task-owned branch tip, its single patch-equivalent replay on
`main`, the replay's direct parent, the stable patch ID, and every path changed
by the task-owned commit. Paths registered in WG are a subset of these exact
commit-owned paths where a task also changed a shared validator.

### 1. `prepare-vgp-data`

- Task branch: `wg/agent-176/prepare-vgp-data`
- Task-owned commit: `362f8e0cf0dcc63d084a3bca2259a446ba593028`
- Integrated commit: `cfa6f34209ad7cd4c2ab9364f6dfdad570f3f3e0`
- Integrated parent: `04183df47bd7587612002ced63c332c2c6bcdb11`
- Stable patch ID on both commits:
  `92a2d9d51790dadcbe31a0965017aded97b8c49f`
- Exact task-owned paths:
  `analysis/prepare_vgp_data_root.py`,
  `analysis/tests/test_prepare_vgp_data_root.py`,
  `analysis/vgp_data_root_config.json`,
  `analysis/vgp_data_root_contract.md`, and
  `analysis/vgp_data_root_validation.json`.
- Contract result: exact root
  `/moosefs/erikg/lewontin-paradox-data/vgp/phase1-freeze-1.0`, mode `2770`,
  owner/group `erikg:erikg`, non-world-writable, same-filesystem atomic
  promotion demonstrated, and acquisition readiness `false` because the quota
  interface is unavailable.

### 2. `freeze-vgp-manifest`

- Task branch: `wg/agent-179/freeze-vgp-manifest`
- Task-owned commit: `c0925628c23d902b6582790c4da2f685ffd61958`
- Integrated commit: `4fe9d3f3c0c652f1aadd73251b4347a62b4454e9`
- Integrated parent: `cfa6f34209ad7cd4c2ab9364f6dfdad570f3f3e0`
- Stable patch ID on both commits:
  `8daf2ba114b524f63be9c9ffe6056d3eb13be73f`
- Exact task-owned paths:
  `analysis/freeze_vgp_manifest.py`,
  `analysis/tests/test_freeze_vgp_manifest.py`,
  `analysis/validate_tier3_guix.sh`,
  `analysis/vgp_phase1_freeze_provenance.json`,
  `analysis/vgp_pilot_manifest.tsv`,
  `analysis/vgp_pilot_rejections.tsv`, and
  `analysis/vgp_pilot_size_budget.tsv`.
- Registered external artifact:
  `/moosefs/erikg/lewontin-paradox-data/vgp/phase1-freeze-1.0/manifests/VGPPhase1-freeze-1.0.commit-dc1b2af5a7741b97d66fb10cb2bce97f41765cdf.tsv`.
- External-object audit: source commit
  `dc1b2af5a7741b97d66fb10cb2bce97f41765cdf`, 327,466 bytes, 717 lines,
  SHA-256
  `9c58420484a8b76a2d6175b7c26bf709e68bdc726a67fc7541b8c2b5a2fc13a4`.
  The live file's size, line count, and digest match the committed provenance.
- Resolver result: 74 candidate rows, 74 rejection rows, and zero selected rows.

### 3. `inventory-pilot-independent`

- Task branch: `wg/agent-184/inventory-pilot-independent`
- Task-owned commit: `e73f4c14941fe2799da6f903182ebe42c8772126`
- Integrated commit: `b084df1cdc1eef328bec67a2f35741bb2404d750`
- Integrated parent: `4fe9d3f3c0c652f1aadd73251b4347a62b4454e9`
- Stable patch ID on both commits:
  `11417ab0425713bf56a50d27aeb1b509fd7df126`
- Exact task-owned paths:
  `analysis/build_vgp_pilot_ne_inventory.py`,
  `analysis/validate_vgp_pilot_ne_inventory.py`,
  `analysis/tests/test_vgp_pilot_ne_inventory.py`,
  `analysis/validate_tier3_guix.sh`,
  `analysis/vgp_pilot_ne_inventory.md`,
  `analysis/vgp_pilot_ne_sources.tsv`, and
  `analysis/vgp_pilot_population_data_availability.tsv`.
- Schema result: both TSVs retain their full schema headers and zero data rows,
  matching the zero selected-species denominator. The validator rejects
  circular primary rows and requires manual curation before selected species
  can be populated. The controlling source schema and method classification
  are also present in `analysis/vertebrate_ne_source_schema.tsv` and
  `analysis/vertebrate_ne_method_matrix.tsv`.

### 4. `gate-vgp-pilot`

- Task branch: `wg/agent-183/gate-vgp-pilot`
- Task-owned commit: `56ab031aed77af520bdcbcd4ab42f93c8dbfe58e`
- Integrated commit: `003a757323bbe592fe5d044ae4691db58b9b06ce`
- Integrated parent: `b084df1cdc1eef328bec67a2f35741bb2404d750`
- Stable patch ID on both commits:
  `ba96bc7ea78b633f6cf174f035f04e437ee39cc8`
- Exact task-owned paths:
  `analysis/gate_vgp_pilot.py`,
  `analysis/tests/test_gate_vgp_pilot.py`,
  `analysis/vgp_pilot_gate.json`, and
  `analysis/vgp_pilot_gate_review.md`.
- Preserved authorization state: `NO_GO`, decision SHA-256
  `ae67031eaa4781984cb9da31ee7a9cd18ee4f2667f98e2a201918fc75ae57284`,
  manifest SHA-256
  `f27b81f369af18caf97a1ccee1b14d8a8b050d665832e6e99f188041275f49ce`,
  root-contract SHA-256
  `912e47c66eeacaf8db3ff535fb857b296eb890d56e2ba0f38af3fdeb4ec12989`,
  and cap-vector SHA-256
  `d51827af5600976f52549ce1a16a7f38bef03f671814d351ff2270ebe7b1b59e`.

### 5. `acquire-vgp-pilot`

- Task branch: `wg/agent-190/acquire-vgp-pilot`
- Task-owned commit: `41925371aa13e8bd6d5c8ad30f964ade816ee83b`
- Integrated commit: `41d8aa1418e42e457822ebbd83fb5a7312e75309`
- Integrated parent: `003a757323bbe592fe5d044ae4691db58b9b06ce`
- Stable patch ID on both commits:
  `7537b2f75de1a969922ebd533bb9c27fbb5e7605`
- Exact task-owned paths:
  `analysis/acquire_vgp_pilot.py`,
  `analysis/gate_vgp_pilot.py`,
  `analysis/tests/test_acquire_vgp_pilot.py`,
  `analysis/tests/test_gate_vgp_pilot.py`,
  `analysis/vgp_pilot_acquisition_manifest.tsv`, and
  `analysis/vgp_pilot_acquisition_report.md`.
- Preserved acquisition state: `refused_preflight`, gate `NO_GO`, refused
  before the first biological byte, zero bytes transferred, zero immutable
  objects promoted, zero quarantine objects, and zero accession views.

### 6. `run-vgp-tier3`

- Task branch: `wg/agent-193/run-vgp-tier3`
- Task-owned commit: `234f31131ddc3ad77498161b9b1e91034af29176`
- Integrated commit: `2b12df1d76d4882b21c0af8ee7a3463dc30421d2`
- Integrated parent: `41d8aa1418e42e457822ebbd83fb5a7312e75309`
- Stable patch ID on both commits:
  `d9cbe8854df52511bf3d3b02c5f0522ee0ee0211`
- Exact task-owned paths:
  `analysis/run_vgp_pilot.py`,
  `analysis/tests/test_run_vgp_pilot.py`,
  `analysis/vgp_pilot_results.tsv`,
  `analysis/vgp_pilot_run_manifest.tsv`, and
  `analysis/vgp_pilot_slurm_telemetry.tsv`.
- Preserved run state: one refusal-summary telemetry row,
  `status=refused_preflight`, `final_state=NOT_SUBMITTED`, blank Slurm job and
  array IDs, no `sbatch` command, `failure_code=GATE_NO_GO`, and zero CPU,
  wall, scratch, read, write, and metadata-operation metrics.

### 7. `review-vgp-pilot`

- Task branch: `wg/agent-197/review-vgp-pilot`
- Task-owned commit: `d224f409791d994a47c144d3c267ad797b778a86`
- Integrated commit: `cbac04e234dad93b2acd07e5f251c37da404a133`
- Integrated parent: `2b12df1d76d4882b21c0af8ee7a3463dc30421d2`
- Stable patch ID on both commits:
  `44850dd521fc0e0de8079862902a44cd58ad420a`
- Exact task-owned paths:
  `analysis/review_vgp_pilot.py`,
  `analysis/tests/test_review_vgp_pilot.py`,
  `analysis/vgp_pilot_qc.tsv`,
  `analysis/vgp_pilot_resource_calibration.tsv`, and
  `analysis/vgp_pilot_review.md`.
- Preserved review: overall `FAIL`, promoted gate `NO_GO`, promoted Slurm
  state `NOT_SUBMITTED`, QC counts `PASS=13`, `FAIL=6`, and resource counts
  `PASS=9`, `INCONCLUSIVE=3`. Zero use is treated as refusal evidence, not as
  completed-pilot performance telemetry.

### 8. `synthesize-vgp-pilot`

- Task branch: `wg/agent-201/synthesize-vgp-pilot`
- Task-owned commit: `255ce708e2d70a9605565b557e1b43ecdd234e9d`
- Integrated commit: `9480b64f848d168ce3e44dbf27924589b9a0170d`
- Integrated parent: `cbac04e234dad93b2acd07e5f251c37da404a133`
- Stable patch ID on both commits:
  `088ed6d94f6d94ed5a0955678316b901fae9540f`
- Exact task-owned paths:
  `analysis/synthesize_vgp_pilot.py`,
  `analysis/tests/test_synthesize_vgp_pilot.py`,
  `analysis/vgp_pilot_next_decision.tsv`,
  `analysis/vgp_pilot_paper_table.tsv`, and
  `analysis/vgp_pilot_synthesis.md`.
- Preserved conclusion: `stop_repair`; zero validated executable species;
  no promoted cross-species diversity, composition, population-genomic, or
  demographic result; all non-zero future scenarios require a new explicit
  authorization.

## Refusal evidence that must not be rewritten

The following values are scientific results of this bounded attempt, not
temporary placeholders:

- Frozen catalog: 714 unique species and 717 raw lines.
- Candidate manifest: 74 data rows; selected rows: 0.
- Rejections: 74 data rows, one for every candidate.
- Composition-ready rows: 0; diversity-ready rows: 0.
- Gate: `NO_GO` with the immutable decision digest listed above.
- Acquisition: refused at preflight before the first biological byte; zero
  transferred bytes and zero promoted objects.
- External immutable-object, staging, and pilot-output file counts: 0, 0, 0.
- Run: `NOT_SUBMITTED`; no Slurm job ID, array ID, dependency, or `sbatch`
  command; all execution resource metrics are zero.
- Independent Ne/ecological source rows: 0; population-data availability
  rows: 0. No Ne was inferred.
- Review: `FAIL`; synthesis: `stop_repair`.

These records must not be changed to imply that a biological pilot ran. The
sentinel SweepGA/IMPG smoke evidence is toolchain evidence only and is not a
pilot-species result.

## Known repair blockers

The current signed gate records five independent blockers:

1. `SOURCE_COUNT_DISCREPANCY_UNRESOLVED`: observed completed-plus-annotated
   count 118 versus expected 120; triple-eligible count 38 versus expected 40;
   triple-eligible fish 14 versus expected 13; completed-plus-RefSeq fish 49
   versus expected 46. The 714 unique-species and 223 completed counts match.
2. `NO_SELECTED_ROWS`: no row met every bounded-pilot selection requirement.
3. `ZERO_COMPOSITION_ELIGIBLE_ROWS`: no exact H1/native-annotation row had
   all finite checksum, size, callable-base, and queryable-gene denominators.
4. `ZERO_DIVERSITY_ELIGIBLE_ROWS`: no paired same-individual diversity row
   met every exact identity, checksum, size, and denominator requirement.
5. `QUOTA_UNAVAILABLE`: filesystem and inode headroom were visible, but no
   exact user-visible quota interface was available; acquisition remains
   fail-closed despite free-space evidence.

The row audit further records 74 missing/non-native annotation cases, 74
unresolved callable-base values, 74 unresolved queryable-gene counts, 74
unresolved queryable-gene-base values, missing H1 SHA-256 evidence for all 74
rows, ten declared download-size mismatches, three invalid H2 accessions, four
pair mismatches, and other row-level H2 checksum/size/URL deficiencies.

Repair must resolve these facts from authoritative, versioned metadata and a
reported quota interface. It must not backfill an ineligible row, reinterpret
the zero selection as execution, or use the zero execution metrics as an
empirical biological or cluster calibration.

## Validator and artifact coverage

The integrated main includes all required executable generators and
validators:

- external-root generator and test:
  `analysis/prepare_vgp_data_root.py`,
  `analysis/tests/test_prepare_vgp_data_root.py`;
- frozen-manifest resolver and test:
  `analysis/freeze_vgp_manifest.py`,
  `analysis/tests/test_freeze_vgp_manifest.py`;
- fail-closed gate and mutation/refusal tests:
  `analysis/gate_vgp_pilot.py`,
  `analysis/tests/test_gate_vgp_pilot.py`;
- acquisition generator and test:
  `analysis/acquire_vgp_pilot.py`,
  `analysis/tests/test_acquire_vgp_pilot.py`;
- run generator and test:
  `analysis/run_vgp_pilot.py`,
  `analysis/tests/test_run_vgp_pilot.py`;
- independent Ne inventory generator, explicit schema validator, and test:
  `analysis/build_vgp_pilot_ne_inventory.py`,
  `analysis/validate_vgp_pilot_ne_inventory.py`,
  `analysis/tests/test_vgp_pilot_ne_inventory.py`;
- review generator and test:
  `analysis/review_vgp_pilot.py`,
  `analysis/tests/test_review_vgp_pilot.py`;
- synthesis generator and test:
  `analysis/synthesize_vgp_pilot.py`,
  `analysis/tests/test_synthesize_vgp_pilot.py`;
- pinned complete analysis entry point:
  `analysis/validate_tier3_guix.sh`, using
  `analysis/guix/channels.scm` and `analysis/guix/manifest.scm`.

The versioned external-root contract and validation, frozen provenance,
candidate/rejection/budget tables, gate/acquisition/run/review/synthesis
artifacts, refusal telemetry, paper table, decision table, and both Ne
inventory schemas are all present on the audited main lineage.

## Validation executed during integration

All validation was performed without downloading biological data, submitting
jobs, or inferring Ne.

1. `./analysis/validate_tier3_guix.sh`
   - Result: PASS.
   - Reproducible profile:
     `/gnu/store/z9v2f6faha9cwjz0sm5iphhlzisgi077-profile`.
   - Result inside the pinned pure shell: 27 tests passed.
   - Ne schema validator result:
     `VGP_PILOT_NE_INVENTORY_OK selected_count=0 source_rows=0 availability_rows=0`.
2. Pinned targeted command using the same channels and manifest:
   `python3 -m pytest -q` over
   `test_prepare_vgp_data_root.py`, `test_freeze_vgp_manifest.py`,
   `test_vgp_pilot_ne_inventory.py`, `test_gate_vgp_pilot.py`,
   `test_acquire_vgp_pilot.py`, `test_run_vgp_pilot.py`,
   `test_review_vgp_pilot.py`, and `test_synthesize_vgp_pilot.py`.
   - Result: PASS, 24 tests passed.
3. Stable patch-ID comparison for all eight required task commits and the
   controlling plan/quality commits.
   - Result: every task-owned patch has exactly one matching integrated patch.
4. External frozen-manifest `stat`, SHA-256, and line-count verification.
   - Result: exact match to committed provenance.
5. Read-only inspection of the authorized external object, staging, and pilot
   output subtrees.
   - Result: zero files in each subtree.

The required repair lineage can therefore proceed from the preserved
`stop_repair`/`NO_GO` state without reconstructing, duplicating, or falsifying
the bounded pilot history.
