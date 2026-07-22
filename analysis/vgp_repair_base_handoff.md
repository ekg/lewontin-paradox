# Corrective VGP implementation base

Status: **CURRENT IMPLEMENTATION BASE**

WG task: `integrate-vgp-repair-base`

Audited main head: `0c20b121b25a8ff139e8c9704b6fa9f4f0de743f`

Biological jobs submitted by this integration: **0**

## Outcome

The useful mirror, BGZF, SweepGA/IMPG, exact-annotation, FastGA scratch,
PSMC-bootstrap-repair, and raw-read sensitivity implementations are retained on
main for a clean corrective canary. Prior execution packets have not been
deleted or rewritten. Their exact bytes and explicit `HISTORICAL` or
`SUPERSEDED` status are bound by
`analysis/vgp_repair_base_artifact_status.tsv`.

This base does **not** adopt any of the following inferences:

- P07 is assembly-invalid. The prior raw-read comparison is a historical,
  paired method/callability sensitivity result. It does not by itself establish
  assembly invalidity or authorize biological exclusion.
- A mapping, scratch, IMPG, scheduler, resource, or other technical failure is
  a biological exclusion. Such outcomes remain retryable or terminal technical
  provenance for that attempt only.
- Two completed pair packets constitute VGP scale-out. The prior two-result
  accounting and synthesis reports are preserved but superseded as scale-out
  conclusions.

## Dependency-correct selected order

1. `f83edbdc8c84c0993ac7b6a3e0239d734ea5e09f` — canonical shared-root and
   immutable authorization contracts.
2. `8e62848e11d7974f674a2bb24043a4fcb116ba2e` — completed resumable mirror,
   exception isolation, provenance, and tests.
3. `47c9d91953dd6fcbd491c72973166251ee2dea83` — BGZF conversion, resource
   classes, pinned Guix environment, atomic triplet promotion, and tests.
4. `3aea32edb52a2bdf9b2eef7f1329dba743d5ff26` — canary SweepGA/IMPG workflow
   prerequisite; its prior outputs are historical.
5. `a017809f44814f660093bec55b986b91daf9bcce` — generalized SweepGA/IMPG,
   annotation, audit, telemetry, and private FastGA scratch workflow; prior
   pilot outputs are historical.
6. `ee451c91bdcabebdfcfd17d77bcce91f428f6361` — repaired PSMC bootstrap
   constructor sampling frozen blocks of the primary PSMCFA N/K/T population.
7. `2044765ca50fc6526c7656da161bcc3f113d1276` — final requested/resolved
   `/scratch` contract, live `/proc` evidence, alias-safe cleanup, and realized
   Guix-profile fixes. Its two-result scale-out packet is superseded.
8. `8a605e310d3ff55fc4960a0a667335d558b6b7a7` — immutable raw-read acquisition
   and checksum tooling.
9. `d3ca516541f37698bd8afbb8bad64d6de54ba160` — pinned raw-read mapping,
   callable-mask, pileup, k-mer, paired-covariance, and PSMC sensitivity tools.
   Its prior pair disposition and interpretive report are superseded.

The independent review integration
`f7e48cba9294ebf7a233dd4798e443f347bdf240` is retained as historical evidence
that identified the bootstrap defect. The synthesis integration
`0c20b121b25a8ff139e8c9704b6fa9f4f0de743f` is not selected into this
corrective base: its files remain on main with explicit sidecar status, but its
conclusions have no gating authority.

For reproducibility, the machine manifest also records the observed live-main
order before this integration (`f83edbd`, `8a605e3`, `3aea32e`, `8e62848`,
`47c9d91`, `a017809`, `f7e48cb`, `ee451c9`, `2044765`, `d3ca516`,
`0c20b12`). The numbered order above is the dependency-correct selection order
for the corrective base, not a claim that the historical squashes arrived in
that order.

## Source-patch audit

`analysis/vgp_repair_base_manifest.json` enumerates every audited source-branch
commit in source order, including selected implementation patches,
implementation patches with historical/superseded outputs, already-integrated
dependencies, merge-only topology, historical-output-only patches, and skipped
conclusion-only patches. The manifest is generated once from the retained
source refs by `analysis/build_vgp_repair_base.py`; validation thereafter uses
the frozen enumeration and does not depend on branch names surviving cleanup.

## Contracts available to the corrective canary

- FastGA creates a private lexical `/scratch` work tree, stages FASTAs and
  indexes there, exports `TMPDIR`, `TMP`, and `TEMP`, changes cwd, audits live
  descendant cwd/open paths through `/proc`, distinguishes requested and
  resolved scratch roots, and performs prefix- and alias-guarded cleanup.
- PSMC bootstrap inputs are assembled from frozen, boundary-aware primary
  PSMCFA N/K/T bins. A replicate resamples those units rather than reconstructing
  a shifted callable-only population.
- The BGZF view has size-derived CPU, memory, scratch, walltime, and concurrency
  classes plus BGZF/FAI/GZI identity and random-access validation under the
  pinned Guix channel.
- SweepGA/IMPG has pinned companions, sidecar construction, exact multiplicity,
  resumable partition/query/lace stages, lossless compressed-stream decoding,
  and exact annotation dictionary gates.
- Raw-read validation preserves same-individual covariance and exposes depth
  masks, pileup evidence, k-mer QV/heterozygosity, paired callset bounds, and
  PSMC sensitivity without converting method discordance into a biological
  exclusion.

Run `python3 analysis/build_vgp_repair_base.py check` for the fail-closed local
contract check. The pinned full-suite invocation remains the project Guix
time-machine command documented in the task validation log; it executes tests
only and does not submit Slurm or biological work.
