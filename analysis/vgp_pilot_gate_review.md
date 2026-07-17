# VGP pilot gate review

- Decision: `NO_GO`
- Decision SHA-256: `ae67031eaa4781984cb9da31ee7a9cd18ee4f2667f98e2a201918fc75ae57284`
- Manifest digest: `f27b81f369af18caf97a1ccee1b14d8a8b050d665832e6e99f188041275f49ce`
- Root contract digest: `912e47c66eeacaf8db3ff535fb857b296eb890d56e2ba0f38af3fdeb4ec12989`
- Cap vector digest: `d51827af5600976f52549ce1a16a7f38bef03f671814d351ff2270ebe7b1b59e`

## Reproduced counts

- Raw frozen catalog lines: `717`
- Raw frozen candidate seeds: `74`
- Selected rows in manifest: `0`
- Independently composition-ready rows: `0`
- Independently diversity-ready rows: `0`

## Winning caps

- `aggregate_core_hours`: `0.0` core-h (winner: `resource_budget_selected_stage_sum_high`)
- `aggregate_wall_hours`: `0.0` h (winner: `resource_budget_selected_stage_sum_high`)
- `cpus_per_element`: `8.0` count (winner: `execution_plan_small_cap`)
- `file_inodes`: `0.0` count (winner: `resource_budget_selected_stage_sum_high`)
- `metadata_operations`: `0.0` count (winner: `resource_budget_selected_stage_sum_high`)
- `moosefs_read_gb`: `0.0` GB (winner: `resource_budget_selected_stage_sum_high`)
- `moosefs_write_gb`: `0.0` GB (winner: `resource_budget_selected_stage_sum_high`)
- `pause_fraction`: `0.8` fraction (winner: `decision_D018`)
- `peak_bandwidth_mib_s`: `0.0` MiB/s (winner: `resource_budget_selected_stage_max_high`)
- `peak_local_scratch_gb`: `0.0` GB (winner: `resource_budget_selected_stage_max_high`)
- `peak_memory_gib_per_job`: `0.0` GiB (winner: `resource_budget_selected_stage_max_high`)
- `persistent_input_gb`: `0.0` GB (winner: `selected_manifest_exact`)
- `persistent_output_gb`: `0.0` GB (winner: `resource_budget_selected_stage_sum_high`)
- `quota_headroom_fraction`: `0.25` fraction (winner: `decision_D018`)
- `selected_pair_species_count`: `0.0` count (winner: `selected_pair_count`)
- `selected_species_count`: `0.0` count (winner: `selected_manifest_count`)
- `stop_fraction`: `0.95` fraction (winner: `decision_D018`)
- `tier3a_concurrency`: `2.0` count (winner: `execution_plan_small_cap`)
- `tier3c_concurrency`: `8.0` count (winner: `execution_plan_small_cap`)
- `transfer_concurrency`: `2.0` count (winner: `execution_plan_small_cap`)

## Blockers

- `SOURCE_COUNT_DISCREPANCY_UNRESOLVED`: the frozen raw VGP catalog still disagrees with the planning headline counts and no explicit signed discrepancy resolution is bundled here (analysis/vgp_phase1_freeze_provenance.json)
- `NO_SELECTED_ROWS`: the frozen pilot manifest selects zero rows, so no bounded pilot can be authorized (analysis/vgp_pilot_manifest.tsv)
- `ZERO_COMPOSITION_ELIGIBLE_ROWS`: no manifest row independently satisfies exact H1/native-annotation/denominator requirements (analysis/vgp_pilot_manifest.tsv)
- `ZERO_DIVERSITY_ELIGIBLE_ROWS`: no manifest row independently satisfies the paired same-individual diversity requirements (analysis/vgp_pilot_manifest.tsv)
- `QUOTA_UNAVAILABLE`: No user-visible quota command was available from the current environment. (analysis/vgp_data_root_validation.json)

## Dominant row-level failures

- `ANNOTATION_MD5_MISSING`: `74` rows
- `ANNOTATION_NOT_NATIVE`: `74` rows
- `ANNOTATION_SHA256_MISSING`: `74` rows
- `ANNOTATION_SIZE_MISSING`: `74` rows
- `ANNOTATION_URL_MISSING`: `74` rows
- `CALLABLE_BASES_UNRESOLVED`: `74` rows
- `DECLARED_DOWNLOAD_SIZE_MISMATCH`: `10` rows
- `H1_SHA256_MISSING`: `74` rows
- `H2_ACCESSION_INVALID`: `3` rows
- `H2_MD5_MISSING`: `3` rows
- `H2_SHA256_MISSING`: `25` rows
- `H2_SIZE_MISSING`: `3` rows
- `H2_URL_MISSING`: `3` rows
- `PAIR_NOT_SAME_INDIVIDUAL`: `4` rows
- `QUERYABLE_GENE_BASES_UNRESOLVED`: `74` rows
- `QUERYABLE_GENE_COUNT_UNRESOLVED`: `74` rows
