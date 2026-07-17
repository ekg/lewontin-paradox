# VGP pilot review

- Overall decision: `FAIL`
- QC counts: `PASS=13`, `FAIL=6`, `INCONCLUSIVE=0`
- Resource calibration counts: `PASS=9`, `FAIL=0`, `INCONCLUSIVE=3`
- Promoted gate decision: `NO_GO` with decision SHA-256 `ae67031eaa4781984cb9da31ee7a9cd18ee4f2667f98e2a201918fc75ae57284`
- Promoted Slurm terminal state: `NOT_SUBMITTED`

## Headline findings

- `source_catalog_counts`: The frozen catalog still disagreed with planning headline counts, so the gate correctly stayed NO_GO. (analysis/vgp_phase1_freeze_provenance.json)
- `selected_manifest_rows`: No rows were selected for the bounded pilot. (analysis/vgp_pilot_manifest.tsv)
- `composition_ready_rows`: No manifest row independently satisfied exact H1/native-annotation/denominator requirements. (analysis/vgp_pilot_manifest.tsv)
- `diversity_ready_rows`: No manifest row independently satisfied paired same-individual diversity requirements. (analysis/vgp_pilot_manifest.tsv)
- `quota_interface`: Filesystem free space existed, but no exact user-visible quota interface was available, so the gate failed closed. (analysis/vgp_data_root_validation.json)
- `scientific_validity`: The promoted manifest never crossed the scientific validity threshold for pilot execution. (analysis/vgp_pilot_manifest.tsv; analysis/vgp_pilot_gate.json)

## Independent recomputation

- Fresh gate recomputation matched the promoted gate on stable fields, including blocker set, cap-vector digest, manifest/root digests, row-audit summary, and quota evidence.
- Fresh refusal-path reruns of `analysis/run_vgp_pilot.py` reproduced the promoted run manifest, refusal telemetry, and result rows after normalizing timestamps, run IDs, and absolute worktree prefixes.
- The promoted refusal path therefore appears immutable and internally consistent even though it did not authorize any biological execution.

## Scientific and execution evidence

- `analysis/sweepga_origin_main_build.json` still records the accepted native `--num-mappings 1:1` SweepGA build with SHA-256 `fa7f0edb9b7e275c288db254046020e136d4267dd5ee043379227ef80da0573b`.
- `analysis/sweepga_impg_observed.json` still records native exact-assembly annotation linkage, callable denominator `14507`, queryable gene count `3`, and 1:1 mapping depth.
- The manifest never crossed the scientific gate: `selected_rows=0`, `composition_ready=0`, `diversity_ready=0`.
- Dominant unresolved row defects remained `ANNOTATION_NOT_NATIVE=74`, `CALLABLE_BASES_UNRESOLVED=74`, `QUERYABLE_GENE_COUNT_UNRESOLVED=74`, and `QUERYABLE_GENE_BASES_UNRESOLVED=74`.
- `analysis/vgp_data_root_validation.json` still reports filesystem headroom but no user-visible quota interface, so quota evidence remained fail-closed rather than advisory.

## Resource calibration

- Aggregate observed CPU, wall, scratch, read, write, and metadata usage stayed at zero because the gate refused preflight before any `sbatch` command was issued.
- Zero observed usage stayed within the zero executable cap that resulted from the empty selected manifest, even though the broader stratified pilot envelopes in `analysis/vertebrate_scaleout_resource_budget.tsv` remained non-zero.
- Per-job request and peak metrics such as CPUs, RSS, and bandwidth are `INCONCLUSIVE`, not because telemetry is missing from executed jobs, but because no jobs existed to observe after the refusal boundary fired.

## Array/network and recommendation

- No Slurm array job ID, dependency, or submission command was recorded in `analysis/vgp_pilot_slurm_telemetry.tsv`; no compute-array network activity could therefore have occurred.
- No authorized species existed, so no coding/CDS/fourfold pilot outputs were scientifically promoted beyond the refusal artifacts.
- Recommendation: keep the pilot at `FAIL` / `NO_GO`. Do not authorize another wave, full-catalog consideration, or demographic follow-up until the source-count discrepancy, quota interface, and row-level reference/annotation/denominator defects are resolved under a new explicit authorization step.
