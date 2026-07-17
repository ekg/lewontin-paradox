# VGP pilot execution quality review

Date: 2026-07-17 UTC

WG task: `quality-vgp-pilot`

## Scope

This review covered the live WG task definitions for:

- `prepare-vgp-data`
- `freeze-vgp-manifest`
- `gate-vgp-pilot`
- `acquire-vgp-pilot`
- `run-vgp-tier3`
- `inventory-pilot-independent`
- `review-vgp-pilot`
- `synthesize-vgp-pilot`

Controlling integrated artifacts reviewed before editing:

- `analysis/vertebrate_scaleout_integration_handoff.md`
- `analysis/vertebrate_scaleout_execution_plan.md`
- `analysis/vertebrate_scaleout_execution_graph.tsv`
- `analysis/vertebrate_scaleout_decisions.tsv`
- `analysis/vertebrate_scaleout_resource_budget.tsv`

This task did not create the external data root, download data, or submit jobs.

## Findings and edits

### 1. Resource-cap drift existed in the live graph

Before this review, the live pilot tasks had the correct overall dependency
shape, but `gate-vgp-pilot`, `acquire-vgp-pilot`, and `run-vgp-tier3` were
using partially hard-coded pilot ceilings that did not clearly derive from the
integrated decision and budget artifacts. The task set also mixed the older
eight-slot planning language with the now-authorized live ceiling of at most
six species.

Edits made:

- `prepare-vgp-data`
  - clarified that the exact external root is
    `/moosefs/erikg/lewontin-paradox-data/vgp/phase1-freeze-1.0`
  - made quota/inode unknown state a downstream blocker rather than something
    that can be guessed
  - made the transfer contract explicit that compute arrays never download
- `freeze-vgp-manifest`
  - reconciled the live task with the six-species ceiling
  - required finite size/checksum evidence for selected rows
  - required per-row and aggregate budget inputs for the later GO gate
- `gate-vgp-pilot`
  - removed reliance on ad hoc caps
  - required every numeric cap to be derived from the integrated decisions,
    execution plan, resource-budget rows, selected manifest, and live quota
    evidence, with the more restrictive value winning per dimension
  - required NO_GO when any required size/checksum/quota/inode/telemetry input
    is unavailable
- `acquire-vgp-pilot`
  - made `analysis/vgp_pilot_gate.json` the authoritative cap vector consumed
    before first byte
  - made acquisition the only task allowed to download biological payloads
  - required mechanical enforcement of gate-derived byte/I-O/object caps
- `run-vgp-tier3`
  - required Slurm requests to come from genome size plus empirical
    high-bound telemetry and stage formulas in
    `analysis/vertebrate_scaleout_resource_budget.tsv`
  - explicitly prohibited arbitrary low memory fallback
  - required network isolation for compute arrays and fail-closed behavior when
    safe memory/scratch requests cannot be derived
- `inventory-pilot-independent`
  - anchored the task to the selected pilot manifest of at most six species
  - reinforced that it may not create a ready demographic execution task
- `review-vgp-pilot`
  - added explicit failure treatment for missing telemetry, network use inside
    arrays, cap overruns, quota-headroom loss, and unresolved reference or
    denominator defects
  - reinforced that review outputs may not create a ready expansion/full-catalog
    task
- `synthesize-vgp-pilot`
  - reconciled synthesis with the live `<=6`-species pilot
  - required any next-step recommendation to stay behind explicit future
    authorization tasks analogous to `A50`/`A60`/`A70`/`A71`

### 2. Dependency graph shape was correct and preserved

The live graph already matched the required fork/join structure and did not
need dependency edits:

- `integrate-vertebrate-scaleout -> quality-vgp-pilot -> prepare-vgp-data ->
  freeze-vgp-manifest -> gate-vgp-pilot -> acquire-vgp-pilot ->
  run-vgp-tier3 -> review-vgp-pilot`
- `inventory-pilot-independent` branches after `freeze-vgp-manifest`
- `synthesize-vgp-pilot` joins `review-vgp-pilot` and
  `inventory-pilot-independent`

No phantom dependency or direct executable bypass around the gate was present
after review.

### 3. Authorization boundaries remain intact

The edited live tasks preserve these boundaries:

- full-catalog execution remains unauthorized
- demographic acquisition and demographic inference remain unauthorized
- bulk biological download occurs only in `acquire-vgp-pilot`
- compute arrays are network-isolated and may not download
- missing quota, inode, size, or checksum metadata now fail closed rather than
  allowing silent continuation

Inference from the user-provided live subgraph:

The integrated plan models explicit `A30`/`A31` pilot decision nodes, but the
requested live execution chain for this review does not include those tasks.
The edited live graph therefore keeps the user-authorized 2026-07-17 bounded
pilot behind the machine-readable `analysis/vgp_pilot_gate.json` boundary
rather than introducing new graph nodes outside the requested chain.

## Automated assertions

I ran automated WG state and description assertions directly against
`.wg/graph.jsonl` using `python3`.

Assertions covered:

- presence of all required tasks
- exact `after` edges for the reviewed subgraph
- exact fork/join shape for `inventory-pilot-independent` and
  `synthesize-vgp-pilot`
- absence of bypass dependencies around `gate-vgp-pilot`
- exact external-root path presence
- six-species ceiling language across the reviewed live tasks
- gate-derived cap-source language
- acquisition-only biological download authorization
- compute-array no-download / network-isolation language
- future authorization boundaries for expansion, full catalog, and demographic
  work

Observed result:

```text
ASSERTIONS_OK
checked_tasks=quality-vgp-pilot,prepare-vgp-data,freeze-vgp-manifest,gate-vgp-pilot,acquire-vgp-pilot,run-vgp-tier3,inventory-pilot-independent,review-vgp-pilot,synthesize-vgp-pilot
phrase_checks=23
```

## Validation summary

- Reviewed and logged all eight named task definitions.
- Reconciled the live graph to an at-most-six-species pilot.
- Replaced ambiguous pilot ceilings with artifact-sourced cap derivation and
  explicit more-restrictive-wins language.
- Preserved the exact dependency chain and fork/join structure.
- Preserved the external-root boundary and bulk-data-outside-Git requirement.
- Preserved the full-catalog and demographic authorization boundaries.
- Confirmed this review task performed no filesystem setup, biological
  acquisition, or Slurm execution.
