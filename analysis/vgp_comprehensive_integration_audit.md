# Comprehensive VGP integration audit

## Result

The comprehensive VGP work is present as a linear, dependency-correct sequence
on `main`. Seventeen source commits from thirteen WG task branches were
promoted as thirteen commits. Stable patch IDs match for every complete
source-branch delta and its promoted commit, including the four branches with
two source commits. No source delta was dropped, and all source-touched blobs
match their promoted versions.

The exact machine-readable mapping, parent chain, source branches, source
commits, dependencies, and stable patch IDs are frozen in
`analysis/vgp_comprehensive_integration_manifest.json`.

## Exact integration order

| Order | Task branch | Source commit(s) | Promoted commit | Required predecessors |
|---:|---|---|---|---|
| 1 | `wg/agent-252/quality-vgp-psmc` | `94d99c766d97d62ec575158ce53b197e74e36d5f` | `089368bba833dc875e5147caf52026cf569145a5` | `integrate-repaired-vgp` |
| 2 | `wg/agent-255/design-vgp-comprehensive` | `37db03165121227dc0543929b0865a8a7d76b1e7` | `3f4b093f080d1ff69d56b4e8bf02b4282ea94b29` | quality |
| 3 | `wg/agent-259/design-gbgc-evidence` | `d94d52317b0b65eaf1a164aa81ba355bd87178d3` | `b6203cceec17a49f54669875ac31d4ee7e6f40d5` | comprehensive design |
| 4 | `wg/agent-258/mirror-vgp-freeze1` | `66a3aaeee3f99fd9ef9cb20f31152fea535e9058` | `e8d6529cfa016c3af505a2cc04af1fb7b5db501f` | comprehensive design |
| 5 | `wg/agent-263/implement-vgp-10-pilot` | `feede1d4dd7fb4d389e09d3695e95a5d93e52e7a` | `0cb473331a8736f3043a4ba8b84b5887fd79e98c` | comprehensive design |
| 6 | `wg/agent-268/pilot-vgp-phylo-gbgc` | `bd95cb7ad5e0ab539a958f7bcb091cd37d439e90`, `5ceb16a71aac288a853aa8041f1d783935b634f6` | `e60e35aa2da5a72bce37bc46bb89ca0aafdc7d13` | gBGC design and mirror |
| 7 | `wg/agent-262/acquire-vgp-10-pilot` | `188613371517615b8eac62b736643b02c0b21a3f` | `b20abfcee603bf243854ecccf8cf8ead7e100867` | comprehensive design |
| 8 | `wg/agent-277/run-vgp-10-pilot` | `0eb2cc80a1a4bd2aebf492993f01d6a8cc95fd36`, `e6cab17bd66d4ff5970c3a76d055851cda28d1a8` | `debfd267a764e7bf5cd6c4e270ed73d9e85a9870` | implementation and acquisition |
| 9 | `wg/agent-275/pilot-pedigree-gbgc` | `25d30504958653eb6a78c471b36f8f621c31908c`, `9ab54a7210b6fb9a1bfd8b68235c2b15acb61ca6` | `877ff60e6c1af5474aa3b8e4a6303791f7fba5aa` | gBGC design and acquisition |
| 10 | `wg/agent-280/review-vgp-10-pilot` | `604ddca09e9e3952cf1298e6cc1b5146fa7f8f55` | `3371eb7f073d073e9fe0fb422dfe321f4d32e565` | ten-pair run |
| 11 | `wg/agent-286/scale-vgp-core` | `f50c0bacfd9320c044d0f1c81e5b7cc8eae17324` | `7fac34c65bdce4c941b5a087b9ede8874065f117` | review and mirror |
| 12 | `wg/agent-283/repair-stale-vgp` | `8aadd57374d5704ca2cd0596f8627a56c9be5edf` | `d71ff7f3cc5c36c22e368e6d319a15ea693bc9c5` | pedigree pilot |
| 13 | `wg/agent-289/synthesize-vgp-program` | `2e3de8fc330ab4b6d3ce39a3574c0f2a5aace256`, `7cf33828bc9e84d01e81871dd2fbf728d4796163` | `f6b337dad2d471003e7714875a0009741f848a7c` | scale-out, phylogenetic pilot, pedigree pilot, and stale-fixture repair |

Siblings were serialized without adding false scientific dependencies. The
stale-fixture repair was deliberately placed after scale-out and before final
synthesis: its pedigree dependency was satisfied, its legacy refusal files did
not overlap the scale-out outputs, and the synthesis was therefore applied to
the repaired tree.

## Reconciliation decisions

- The four two-commit branches were promoted as one commit each. Patch IDs
  were computed over the complete ranges, not only their formatting follow-ups.
- Later siblings were applied over already-promoted, disjoint siblings. File
  inventories and final blob IDs were compared so a later generated artifact
  could not silently replace another task's output.
- The repaired legacy review artifacts were retained byte-for-byte at their
  promoted digests. The repair changes only the stale, host-dependent
  recomputation contract and its legacy evidence packet.
- The old `NO_GO`, zero-job, zero-core-pass, and zero-estimate packets remain
  immutable evidence, but they are explicitly non-binding historical runs in
  `analysis/vgp_historical_run_registry.json`. They are not a current
  authorization packet and their decisions cannot be inherited by a new run.
  A separate versioned packet must establish current execution authorization.

## Integrated execution surface

- The frozen rosters contain exactly ten primaries (`P01`–`P10`) and six
  ranked alternates (`A01`–`A06`).
- The core workflow includes whole-assembly SweepGA `--num-mappings 1:1`,
  IMPG index/partition/query/lace, reason-coded masks, consensus construction,
  PSMC plus 200 boundary-aware bootstraps, and separate scaling scenarios.
- Pinned Guix channel, manifest, package, realization, and executable identity
  records are present.
- Slurm dry-run/submission entry points, measured per-stage resources,
  sentinels, telemetry, bounded arrays, and atomic promotion are present.
- The Freeze 1 mirror has a frozen 43,371-file inventory, resumable worker,
  checksum quarantine, capacity gate, and atomic promotion tooling.
- Acquisition, design, mirror, workflow, review, stale-recomputation,
  scale-out, synthesis, and integration validators and tests are present.

## Operational boundary

This integration performed repository validation only. It submitted no Slurm
jobs and transferred no biological or mirror payload. Historical outcomes do
not grant or deny current execution; they document what happened at their UTC
timestamps.
