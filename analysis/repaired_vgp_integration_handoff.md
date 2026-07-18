# Repaired VGP evidence integration handoff

Date audited: 2026-07-18 UTC
WG task: `integrate-repaired-vgp`

## Integration disposition

The complete repaired-VGP pilot lineage was already present on live `main` and
`origin/main` when this audit began. After an explicit `git fetch origin main`,
all four audited refs resolved to the same commit:

```text
HEAD        f3d9bd23a088e9be8cb68a08cb71e90db1223248
main        f3d9bd23a088e9be8cb68a08cb71e90db1223248
origin/main f3d9bd23a088e9be8cb68a08cb71e90db1223248
FETCH_HEAD  f3d9bd23a088e9be8cb68a08cb71e90db1223248
```

Each of the seven required task-owned commits has exactly one stable-patch-ID
match in the linear `main` history. Replaying a task branch would therefore
duplicate already integrated evidence and tooling. No task commit was replayed;
this handoff is the only new project artifact produced by the integration task.

No reset, rebase, clean, force push, amend, stash, broad staging, biological
download, Slurm submission, or demographic inference was performed.

## Worktree and user-change audit

The assigned worktree
`/moosefs/erikg/lewontin-paradox/.wg-worktrees/agent-249` was clean, had no
commit ahead of `main`, and was already recorded by WG as merged to `main`
before this document was created.

The live root worktree was on `main...origin/main` with no tracked or staged
change. It contained pre-existing untracked state:

- the WG-managed `.wg-worktrees/` directory; and
- seven `analysis/tests/*.pyc` files.

Those user/unrelated paths were not staged, edited, removed, or otherwise
altered. Each source task worktree (`agent-218`, `221`, `222`, `227`, `230`,
`233`, and `236`) remained at its exact task-owned commit with only the
WG-owned untracked `.wg-cleanup-pending` marker. Their former remote task
branches had been removed after integration, but the local refs and commit
objects remain intact. The corresponding WG logs record that every source
commit was pushed before task completion.

## Dependency-correct integration order

The controlling prerequisite chain was already integrated before the repaired
lineage:

```text
04183df47bd7587612002ced63c332c2c6bcdb11  quality-vgp-pilot
  -> cfa6f34209ad7cd4c2ab9364f6dfdad570f3f3e0  prepare-vgp-data
  -> 4fe9d3f3c0c652f1aadd73251b4347a62b4454e9  freeze-vgp-manifest
  -> b084df1cdc1eef328bec67a2f35741bb2404d750  inventory-pilot-independent
  -> 003a757323bbe592fe5d044ae4691db58b9b06ce  gate-vgp-pilot
  -> 41d8aa1418e42e457822ebbd83fb5a7312e75309  acquire-vgp-pilot
  -> 2b12df1d76d4882b21c0af8ee7a3463dc30421d2  run-vgp-tier3
  -> cbac04e234dad93b2acd07e5f251c37da404a133  review-vgp-pilot
  -> 9480b64f848d168ce3e44dbf27924589b9a0170d  synthesize-vgp-pilot
  -> 138d340111efc71c09b50b7de746f02ed91663b6  integrate-vgp-pilot
  -> a91c8b6a49825d1604fa8d0e48cc276c4942d099  quality-vgp-pilot-2
```

The exact branch commits, stable patch IDs, and path inventories for the
original bounded-pilot lineage are preserved in
`analysis/vgp_pilot_lineage_integration.md`. The repair quality prerequisite
was task-owned commit
`31b3db3dd80ca76ff484aa8cfb9e85b3d3e8fcff`, replayed exactly once as
`a91c8b6a49825d1604fa8d0e48cc276c4942d099`; both have stable patch ID
`2a3fe2649e225873efe28c6b78d6c9b5f5b8ecfe`. Its three artifacts are
`analysis/assert_vgp_pilot_repair_wg.py`,
`analysis/tests/test_assert_vgp_pilot_repair_wg.py`, and
`analysis/vgp_pilot_repair_quality.md`.

The repaired lineage then appears on `main` in this dependency-correct order:

```text
a91c8b6  quality-vgp-pilot-2
  -> c23972c  repair-vgp-candidate
       |\
       | -> 9a38270  audit-vgp-demography
       -> b5135ee  regate-vgp-pilot
            -> 5d92d12  acquire-repaired-vgp
                 -> b3beb98  run-repaired-vgp
  9a38270 + b3beb98
            -> 9d90c31  review-repaired-vgp
                 -> f3d9bd2  synthesize-repaired-vgp
```

WG defines `audit-vgp-demography` and `regate-vgp-pilot` as parallel children
of `repair-vgp-candidate`. Serializing the audit commit before the regate
commit is dependency-correct. `review-repaired-vgp` is the join: it depends on
both the run path and the demography audit. The final synthesis depends on that
review.

## Exact task commit and artifact ledger

Stable patch IDs were calculated with `git patch-id --stable` over each whole
commit. Matching source and integrated IDs prove that each task patch occurs
once on `main`, even though WG integration changed its parent and commit ID. A
complete `git rev-list main` scan recomputed the stable patch ID of every
commit and returned count `1` for each of the seven IDs below, at exactly the
listed integrated commit.

| task | pushed task-owned commit | integrated `main` commit | stable patch ID | changed paths |
| --- | --- | --- | --- | ---: |
| `repair-vgp-candidate` | `71b5815b5fe8007c87edeb81f5a47595d84f6328` | `c23972c5b6dceea0f9c197f5016d43b7306f60ac` | `7ac9164a88b4ea9c6b708570e084b8eee2d931a6` | 70 |
| `audit-vgp-demography` | `23f7da933d323416c7ac19d85b52483a1fe25146` | `9a3827074de031245a54feeca276cc268cb8d515` | `2092b1708c61eeba8a86d94182fde7895449a9c0` | 17 |
| `regate-vgp-pilot` | `d694f5a64d355869cd4b7600fbce13753f9615bb` | `b5135eefeea5962827c76e04606655f32f109d22` | `22b3e43fcbdf43e7d9df5e8bbdc6050d31da5f47` | 8 |
| `acquire-repaired-vgp` | `eb40db96d6b00cb84c2efcc19c0e40a12ed90368` | `5d92d12022ad1ad876614741b7a61bea5c0569af` | `51cc6f3304dd23b7da8a807dfd977ef366696ade` | 10 |
| `run-repaired-vgp` | `f2837b751e19f849c11211bf685017bb81a111ab` | `b3beb9881f8f0ae96be95ed2fec0bfeedadc4c1d` | `40858bad97c60389c50b7d57b7c93c2c44dcfb3e` | 11 |
| `review-repaired-vgp` | `61afc23a7f21f543399858f445f0d1db72277060` | `9d90c316a825746d112f5f6bf20d11ad0560e55e` | `bf5e0bd085434d16de163b312d7c4150c3575651` | 5 |
| `synthesize-repaired-vgp` | `4a0c4103a32aef2baeeaec7c2d04cfa6d118051c` | `f3d9bd23a088e9be8cb68a08cb71e90db1223248` | `2d3b6b0eb3a3ff1dc2800d946467206d70a0b0b9` | 5 |

The complete task-owned artifact inventory follows. Shared tooling and tests
are listed because they are part of the task patches, not just the final report
files.

### `repair-vgp-candidate` — 70 paths

Tooling, tests, and primary evidence (13 paths):

- `.gitattributes`
- `analysis/resolve_vgp_candidates.py`
- `analysis/tests/test_acquire_vgp_pilot.py`
- `analysis/tests/test_gate_vgp_pilot.py`
- `analysis/tests/test_resolve_vgp_candidates.py`
- `analysis/tests/test_review_vgp_pilot.py`
- `analysis/tests/test_run_vgp_pilot.py`
- `analysis/tests/test_synthesize_vgp_pilot.py`
- `analysis/validate_tier3_guix.sh`
- `analysis/vgp_pilot_manifest.tsv`
- `analysis/vgp_pilot_rejections.tsv`
- `analysis/vgp_pilot_resolution_report.md`
- `analysis/vgp_pilot_size_budget.tsv`

Immutable metadata-cache evidence (57 paths):

- `analysis/vgp_resolution_cache/checkpoint.json`
- `analysis/vgp_resolution_cache/index.json`
- all 26 keyed JSON records in `analysis/vgp_resolution_cache/entries/`
- all 26 corresponding immutable response blobs in
  `analysis/vgp_resolution_cache/responses/`
- the three preserved prior-refusal TSVs in
  `analysis/vgp_resolution_cache/prior_refusal/`

The cache patterns above are exhaustive for this commit: 26 entries, 26
responses, 3 prior-refusal ledgers, one checkpoint, and one index.

### `audit-vgp-demography` — 17 paths

- `analysis/build_vgp_demography_audit.py`
- `analysis/refresh_vgp_demography_metadata.py`
- `analysis/tests/test_vgp_demography_audit.py`
- `analysis/validate_vgp_demography_audit.py`
- `analysis/vgp_demography_audit.md`
- `analysis/vgp_demography_input_audit.tsv`
- `analysis/vgp_independent_ne_sources.tsv`
- `analysis/vgp_demography_cache/checkpoint.json`
- `analysis/vgp_demography_cache/index.json`
- all 4 keyed JSON records in `analysis/vgp_demography_cache/entries/`
- all 4 corresponding immutable metadata responses in
  `analysis/vgp_demography_cache/responses/`

The directory patterns are exhaustive for this commit: 4 entries, 4
responses, one checkpoint, and one index, in addition to the seven named
tooling/report paths.

### `regate-vgp-pilot` — 8 paths

- `analysis/gate_vgp_pilot.py`
- `analysis/review_vgp_pilot.py`
- `analysis/tests/test_acquire_vgp_pilot.py`
- `analysis/tests/test_gate_vgp_pilot.py`
- `analysis/tests/test_review_vgp_pilot.py`
- `analysis/tests/test_run_vgp_pilot.py`
- `analysis/vgp_pilot_gate.json`
- `analysis/vgp_pilot_gate_review.md`

### `acquire-repaired-vgp` — 10 paths

- `analysis/acquire_vgp_pilot.py`
- `analysis/gate_vgp_pilot.py`
- `analysis/tests/test_acquire_vgp_pilot.py`
- `analysis/tests/test_gate_vgp_pilot.py`
- `analysis/vgp_pilot_acquisition_manifest.tsv`
- `analysis/vgp_pilot_acquisition_refusal.json`
- `analysis/vgp_pilot_acquisition_report.md`
- `analysis/vgp_pilot_gate.json`
- `analysis/vgp_pilot_gate_review.md`
- `analysis/vgp_pilot_immutable_object_inventory.tsv`

### `run-repaired-vgp` — 11 paths

- `analysis/run_vgp_pilot.py`
- `analysis/slurm/run_repaired_vgp.sh`
- `analysis/tests/test_review_vgp_pilot.py`
- `analysis/tests/test_run_vgp_pilot.py`
- `analysis/vgp_pilot_compute.py`
- `analysis/vgp_pilot_exclusions.tsv`
- `analysis/vgp_pilot_refusals.tsv`
- `analysis/vgp_pilot_results.tsv`
- `analysis/vgp_pilot_run_manifest.tsv`
- `analysis/vgp_pilot_run_report.md`
- `analysis/vgp_pilot_slurm_telemetry.tsv`

### `review-repaired-vgp` — 5 paths

- `analysis/repaired_vgp_pilot_qc.tsv`
- `analysis/repaired_vgp_pilot_review.md`
- `analysis/repaired_vgp_resource_calibration.tsv`
- `analysis/review_repaired_vgp.py`
- `analysis/tests/test_review_repaired_vgp.py`

### `synthesize-repaired-vgp` — 5 paths

- `analysis/repaired_vgp_next_decision.tsv`
- `analysis/repaired_vgp_paper_table.tsv`
- `analysis/repaired_vgp_pilot_synthesis.md`
- `analysis/synthesize_repaired_vgp.py`
- `analysis/tests/test_synthesize_repaired_vgp.py`

## Exact repaired metadata evidence

The repaired manifest has exactly one header plus six candidate rows. All six
have an exact-version RefSeq H1 and an official native annotation whose
reference accession equals that H1. Every row resolves only to
`tier3c_composition`; no row is promoted as a Tier3A diploid-diversity or
demography input.

| candidate | TaxId | exact H1 and native-annotation reference | proposed compressed bytes |
| --- | ---: | --- | ---: |
| *Camelus dromedarius* | 9838 | `GCF_036321535.1` | 697,371,781 |
| *Colius striatus* | 57412 | `GCF_028858725.1` | 382,723,741 |
| *Candoia aspera* | 51853 | `GCF_035149785.1` | 493,154,703 |
| *Dendropsophus ebraccatus* | 150705 | `GCF_027789765.1` | 664,270,825 |
| *Lepisosteus oculatus* | 7918 | `GCF_040954835.1` | 355,464,988 |
| *Heterodontus francisci* | 7792 | `GCF_036365525.1` | 1,557,354,377 |

The manifest SHA-256 is
`b7cbe6cf22287d58dae9e270a81dee2970aaa3d8c73a2270323fa6210b276988`.
The catalog units remain distinct and exact: 717 physical lines, 1 header,
716 data rows, and 714 unique species. The two-row data-row excess is explained
by *Lophostoma evotis* and *Micronycteris microtis*, each with multiplicity 2.

## Exact refusal and zero-use evidence

The repaired gate is `NO_GO`, decision SHA-256
`9f39b13be5e0b1999c4cd98498399aee8700a7487d57fa81dfb7c59c29ff867d`.
It has exactly three blockers:

1. `CAP_MOOSEFS_READ_GB_EXCEEDED`: proposed finite worst-case read volume
   406.7274 GB exceeds the strict 180 GB limit.
2. `CAP_SCRATCH_GIB_EXCEEDED`: proposed finite worst-case scratch
   205.5066 GiB exceeds the strict 139.69838619232178 GiB limit.
3. `QUOTA_UNAVAILABLE`: no enforceable user-visible quota report was
   available; filesystem free space was correctly not substituted for quota.

The task-wide ceilings were also preserved, with every stricter integrated cap
winning: at most 6 species, 120 GiB compressed inputs, 750 GiB scratch, 1,500
core-hours, 2 concurrent species, and 256 GiB per job. In this gate the
effective scratch, core-hour, and per-job-memory caps are stricter:
139.69838619232178 GiB, 280 core-hours, and 96 GiB respectively.

Acquisition contains one refusal summary plus one row per blocker. Every row
records zero expected, observed, and cumulative transferred biological bytes;
the immutable-object inventory has a header and zero objects. The acquisition
refusal digest is
`240b90e70b940ea67efae8ba2cf0276752efdf5efb924fabc25f6641e1ab844b`.

The run ledger contains one refusal summary plus the same three blockers. Its
summary records `selected_species=0`, `promoted_objects=0`,
`compressed_input_bytes=0`, `slurm_jobs_submitted=0`, and
`final_state=NOT_SUBMITTED`. The sole telemetry row has a blank `sbatch`
command and blank job/array/dependency fields, with zero CPU, wall, scratch,
I/O, metadata-operation, and network use. This is refusal evidence, not
performance calibration and not a biological result.

## Annotation-pilot refusal is not PSMC infeasibility

The refused pilot and a future assembly-derived PSMC program answer different
questions and have different input contracts.

The current `NO_GO` applies to a six-species **annotation/composition pilot**.
Its candidates passed exact H1/native-annotation metadata checks, but the
proposed batch failed two resource caps and lacked enforceable quota evidence.
Therefore no biological payload was acquired and no composition analysis ran.
Nothing in that control outcome shows that the species lack usable diversity,
that their demographic histories are biologically uninformative, or that PSMC
is intrinsically infeasible.

The metadata-only demography audit separately found all six current rows
ineligible **on the evidence presently assembled** for PSMC, MSMC2, and SMC++.
For PSMC, a haploid H1 assembly and a linked H2 accession do not establish a
heterozygosity-retaining diploid callable consensus, a compatible mask, or
sufficient read/individual/reference/coverage provenance. MSMC2 additionally
needs validated mutually comparable phased haplotypes and their relationships;
SMC++ needs a multi-sample exact-reference population genotype set, population
definitions, masks, and QC. Those are method-specific missing-input findings,
not biological limitations and not permanent impossibility findings.

Consequently, an assembly-derived PSMC scale-out may be scientifically
feasible only after a new bounded design identifies and validates diploid
heterozygosity-bearing inputs, masks, exact identity/provenance, coverage, and
scaling scenarios. It also requires new explicit authorization for any
metadata refresh, biological acquisition, consensus/mask construction, and
inference. The annotation-pilot `NO_GO` must neither be reused as a PSMC
feasibility verdict nor bypassed to authorize that new work.

No PSMC/MSMC2/SMC++ inference was performed here. Coalescent-scaled output must
remain separate from absolute Ne/time until explicit mutation-rate and
generation-time scenarios are justified. Only six dromedary LD-Ne population
records qualify as independent numeric Ne in the bounded audit; the six
same-response `pi/(4mu)` policy rows remain circular and excluded.

## Tooling presence on the integrated tree

The final tree contains and tests every required boundary:

- repaired resolver and immutable/resumable cache:
  `analysis/resolve_vgp_candidates.py` and `analysis/vgp_resolution_cache/`;
- independent gate and digest binding: `analysis/gate_vgp_pilot.py` and
  `analysis/vgp_pilot_gate.json`;
- fail-closed zero-byte acquisition: `analysis/acquire_vgp_pilot.py`;
- fail-closed `NOT_SUBMITTED` execution:
  `analysis/run_vgp_pilot.py`, `analysis/vgp_pilot_compute.py`, and
  `analysis/slurm/run_repaired_vgp.sh`;
- method-specific demography audit and cache:
  `analysis/build_vgp_demography_audit.py`,
  `analysis/refresh_vgp_demography_metadata.py`, and
  `analysis/vgp_demography_cache/`;
- independent review: `analysis/review_repaired_vgp.py` and
  `analysis/repaired_vgp_pilot_review.md`; and
- final synthesis: `analysis/synthesize_repaired_vgp.py` and
  `analysis/repaired_vgp_pilot_synthesis.md`.

## Pinned GNU Guix validation

The full repository analysis suite was executed through the already-realized
pinned profile using the committed environment record:

```sh
analysis/slurm/guix_job.sh \
  "$PWD/analysis/pilot_results/guix_environment.json" \
  python3 -m pytest -q analysis/tests
```

Observed result: **241 passed in 6.32 seconds**.

The wrapper verified the committed channel and manifest digests before
execution. The recorded profile is
`/gnu/store/z9v2f6faha9cwjz0sm5iphhlzisgi077-profile`, derived from
`/gnu/store/x7vw3qvf5yxsff7x5cjhxs713m90ni6n-profile.drv`, Guix channel commit
`44bbfc24e4bcc48d0e3343cd3d83452721af8c36`, channels SHA-256
`45c055cd1d9010a72eacbb720037a22bccb2d8d6891dbd11b5d66365f29b3a17`,
and manifest SHA-256
`2fb05e87aa2ac45ce51d4dcf93b232cb98627f525adace98357629ee3f15720a`.

No test command used a cache-refresh flag, acquisition entrypoint, submission
entrypoint, biological payload, or demographic inference.

## Final authorization boundary

This integration preserves evidence; it does not grant execution authority.
No data was downloaded and no job was submitted. Full-catalog work, raw
population bulk acquisition, genotype construction, assembly-derived PSMC,
MSMC2, SMC++, or any other demographic inference still requires a new explicit
and numerically bounded authorization. The existing `NO_GO`, zero-byte,
zero-object, empty-result, and `NOT_SUBMITTED` records remain the authoritative
outcome of the repaired annotation pilot.
