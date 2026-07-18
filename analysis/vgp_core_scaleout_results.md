# VGP Freeze 1 core scale-out: closed-world gated result

Generated: `2026-07-18T13:15:23Z`. Release: VGP Phase 1 Freeze 1 commit
`dc1b2af5a7741b97d66fb10cb2bce97f41765cdf`, catalog SHA-256 `9c58420484a8b76a2d6175b7c26bf709e68bdc726a67fc7541b8c2b5a2fc13a4`.

## Outcome

**No biological scale-out was authorized or run.** The independent pilot review
is `CONDITIONAL_GO` for bounded repair and ten-slot re-pilot only, explicitly
sets full scale-out to false, and reports zero core pilot passes. The Freeze 1
mirror is a complete metadata inventory but all 47,870 objects remain `planned`;
zero payload objects are verified or reusable. This run therefore submitted
zero Slurm jobs, completed zero pairs, computed zero callable-diversity or PSMC
estimates, and promoted zero pair packets. Technical non-execution is not a
low-diversity result.

## Closed-world accounting

The manifest contains all **716 catalog rows** plus all **569
catalog-linked haplotype entries** (264 labeled other-high-quality and 305
labeled alternate). Catalog-row dispositions are: {'excluded': 187, 'failed': 529}.
Linked-entry dispositions are: {'excluded': 3, 'failed': 566}. Three linked
entries repeat the H1 accession and are excluded as self-links, not pairs. The
remaining 566 distinct pair candidates are failed operationally at the upstream
authorization gate; they are not declared scientifically ineligible. Each pair
has exactly one manifest row and one QC row (569 = 569).

The catalog label “other high-quality” is discovery metadata, not proof that a
pair passes exact-individual provenance, final-sequence QV, completeness,
duplication/collapse, mutual comparability, whole-assembly 1:1, or measured
callability. Those gates remain unmeasured. Hi-C/trio catalog signals are retained
as confidence context; their absence is explicitly not a core veto. Annotation
absence is never a core veto.

## Biological outputs and interpretation

No callable denominator, diversity estimate, mask, consensus, unscaled PSMC,
or bootstrap exists. PSMC requires 200 predeclared boundary-aware replicates
(minimum 190 finite; blocks never cross contigs or mask discontinuities) after a
core pass. The scenario table retains unscaled PSMC as primary and refuses to
invent mutation or generation values. There are no approved confidence tiers,
so there are no scaled trajectories to report.

CDS, fourfold, nonsynonymous, synonymous, WS, SW, and GC3 partitions were not
run. They remain optional post-core outputs requiring an exact native annotation
accession/version and equal sequence dictionary, or a separately validated,
manifest-bound liftover. No valid non-annotated core result was deleted or
downgraded—none yet exists.

Assembly-derived PSMC from an H1/H2 pair is descriptive demographic context.
It reuses the same individual, haplotypes, variants, callable mask, and consensus
as same-pair diversity and **is not statistically independent evidence** for
that diversity.

## Independent validation and sensitivity

Independent biological recomputation, raw-read/k-mer checking, and literature
triangulation were not performed: no pair passed the core gates, the mirror has
no verified sequence payload, and the authorization forbids full processing.
The stratified biological sample is therefore empty rather than post hoc. The
sensitivity table reports QV, completeness, duplication/collapse, long-range
phasing, callability, genome size, generation, bootstrap, and scaling effects as
not estimable. It is invalid to infer any sensitivity from universal technical
non-execution.

## Resumable wave contract

`GATE-0` is the only current state and has a finite zero-job limit. The low/base/
high rows are non-authorized planning templates copied from the reviewed pilot
resource envelope, whose mapping/PSMC prediction error was not estimable. A
future GO requires a new immutable pair audit and wave manifest; the templates
cannot authorize themselves. Each future wave is capped at 25 pairs, uses finite
array concurrency and per-job resources, allows two transient retries and one
resource re-estimation retry, permits zero scientific-threshold relaxation,
checkpoints stage digests, and promotes only a complete output packet by atomic
rename after digest/completeness verification. A job stops at a hard scientific
gate or 1.5× its then-reviewed high resource estimate. These are operational
limits, never a global byte/memory eligibility ceiling. A new assembly release
requires a new manifest and cannot drift into Freeze 1.

## Digest verification

Inputs:

- `catalog`: `9c58420484a8b76a2d6175b7c26bf709e68bdc726a67fc7541b8c2b5a2fc13a4`
- `design`: `17217bef8f6c2c5efaef1111eebef6f97444f5d39b11b5403ccb76c2e48c185c`
- `guix_channels`: `fd7d00be59ace5d82eb093ab2bd4efab45c30a166f3c3e1572d6e2d4d628bed4`
- `guix_manifest`: `afd14da446d29d4e66b64837194dfae524342dc5b8de3c2778ec85c1cca09a3d`
- `mirror_manifest`: `935d18268c8b3b2ca5eb68332282a39db29fef563164adc8b6a493131b4675a0`
- `mirror_summary`: `801ab3ffce1eca0f43d08ac406ded2604839f783221b4ff91cdc856699596a4c`
- `review_decision`: `561a327d2a7d83c004a011d9bac1ec4b45cf44067858a29a3036bb772a57fc7b`
- `review_resources`: `efc3d69119744116ab73cc75a438d2a2ccd9dcd10eeb9b06c95b02cd9b8b456e`
- `review_scaleout`: `4bed1367c2f521a72b3f7ec2e625b141b7e35e6ce2103a8f29b33f0a824805fc`

Outputs (computed after atomic write):

- `independent_validation`: `b280b493e9c23bc6c3ebb6f9d39d32f0cade2ad3898f28e4bf013467d7eebefb`
- `manifest`: `5bf2f6102cbd66e881d5fcaedfb22bfdff518e3f49c2dfbb01bbd493e0b461f6`
- `paper_pairs`: `d76e1c242a26e08328f280b557677b4eece0ef64d0a7a3890a077b34f3243de1`
- `paper_summary`: `24f4f644c33d08a3b417deee0071b2898026d1eee51d7f3c7d5edc6dfc86da6c`
- `qc`: `14f37e1ed85176ab49ddac69ab19f8d2b5a3a22837957482f2b0a21f8218b321`
- `scaling`: `959e79777df0cc334b2a8a95d896d1717da0b66898c1e6fe6d8b441c745641df`
- `sensitivity`: `0cd004efa435383da1714e499cc4cf0607b9bfc9b44fe38e421ee545739f7cca`
- `telemetry`: `9f423f2a296210b36059e66f2db392e822849dc97c3603053050960f6b8d66a3`
- `waves`: `3003ec7239376d79189f4ba4e3801835dcaf3818bd35bc3e9b65a2ed67754806`
