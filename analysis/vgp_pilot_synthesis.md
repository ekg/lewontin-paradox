# VGP pilot synthesis handoff

Date: 2026-07-17 UTC

## Executive readout

- Recommended next decision: `stop_repair`.
- Frozen pilot manifest rows reviewed: `74`; pilot-selected rows: `0`; validated executable species: `0`.
- Actual bounded pilot execution never crossed the authorization boundary: `final_state = NOT_SUBMITTED`, `sbatch_command = none`, and every executable cost metric remained zero.
- Review QC remained fail-closed (`PASS=13`, `FAIL=6`) and resource calibration stayed refusal-only (`PASS=9`, `INCONCLUSIVE=3`).
- Independent Ne/ecological inventory remained intentionally empty because the selected-species denominator was zero; no literature Ne, raw-read, VCF, or callable-genotype claims were populated.

## Scope and provenance

- This synthesis joins `analysis/vgp_pilot_manifest.tsv`, `analysis/vgp_pilot_rejections.tsv`, `analysis/vgp_pilot_results.tsv`, `analysis/vgp_pilot_qc.tsv`, `analysis/vgp_pilot_resource_calibration.tsv`, `analysis/vgp_pilot_ne_sources.tsv`, and `analysis/vgp_pilot_population_data_availability.tsv`.
- The join key is versioned taxon identity `candidate_id|h1_accession_version`; every one of the 74 manifest rows appears in `analysis/vgp_pilot_paper_table.tsv` and every joined Ne/availability count is zero because both inventory files are header-only.
- No additional download, Slurm submission, or demographic inference was launched while producing this handoff.

## Which strata were successfully analyzed

- No vertebrate class or stratum crossed the executable pilot gate. The `<=6`-species bounded pilot analyzed zero species biologically and emitted zero coding/CDS/fourfold estimates.
- What did survive review was only precondition evidence: the accepted SweepGA build remained byte-identical, and the IMPG smoke artifact still showed native 1:1 mapping depth with `callable_bp = 14507`, `queryable_gene_count = 3`, and `queryable_gene_bp = 14507` on the exact-native sentinel assembly.
- Those sentinel mapping/queryability/callability checks are not promoted to any pilot-species diversity, composition, population-genomic, or demographic result.

### Frozen candidate strata reviewed

- `Actinopteri`: 38 candidate rows, 9 Tier3A-capable, 29 Tier3C-only.
- `Amphibia`: 4 candidate rows, 4 Tier3A-capable, 0 Tier3C-only.
- `Aves`: 6 candidate rows, 6 Tier3A-capable, 0 Tier3C-only.
- `Chondrichthyes`: 10 candidate rows, 5 Tier3A-capable, 5 Tier3C-only.
- `Cladistia`: 1 candidate rows, 0 Tier3A-capable, 1 Tier3C-only.
- `Lepidosauria`: 4 candidate rows, 4 Tier3A-capable, 0 Tier3C-only.
- `Mammalia`: 8 candidate rows, 8 Tier3A-capable, 0 Tier3C-only.
- `UNRESOLVED`: 3 candidate rows, 2 Tier3A-capable, 1 Tier3C-only.

## Exact eligibility and failure reasons

- Eligibility stayed zero for all four modalities: every row has `assembly_composition_eligible = no`, `assembly_diversity_eligible = no`, `population_genomic_eligible = no`, and `demographic_eligible = no`.
- Rejection ledger size equals manifest size: `74` rejected rows, `8` mammals, `6` birds, `4` lepidosaurs, `4` amphibians, `38` actinopterygians, `10` chondrichthyans, and `1` cladistian row, plus `3` rows whose class label was unresolved in the frozen source.
- Blocking requirement counts were exact: `B03` alone for `36` Tier3C-only rows and `B02;B03` for `38` Tier3A-capable rows.
- Annotation evidence failed for every row: `annotation_file_status = missing` for all 74 rows and `annotation_native_status = unresolved_or_missing` for every rejection row.
- Same-individual pairing was present for `31` rows but absent for `43` rows; even the paired rows remained blocked because native exact-H1 annotation and denominator evidence were missing.

### Rejection reason ledger

- `missing_exact_h1_annotation_or_file_size_checksum_evidence`: `36` rows.
- `missing_exact_h1_annotation_or_file_size_checksum_evidence;ok`: `31` rows.
- `missing_exact_h1_annotation_or_file_size_checksum_evidence;reject_pair_mismatch`: `4` rows.
- `missing_exact_h1_annotation_or_file_size_checksum_evidence;reject_pair_unresolved`: `2` rows.
- `missing_exact_h1_annotation_or_file_size_checksum_evidence;reject_same_individual_unresolved`: `1` rows.

### Gate blockers copied into the refusal results

- `SOURCE_COUNT_DISCREPANCY_UNRESOLVED`: the frozen raw VGP catalog still disagreed with earlier planning headline counts.
- `NO_SELECTED_ROWS`: the frozen pilot manifest selected zero rows, so no bounded pilot could be authorized.
- `ZERO_COMPOSITION_ELIGIBLE_ROWS`: no row independently satisfied exact-H1/native-annotation/denominator requirements.
- `ZERO_DIVERSITY_ELIGIBLE_ROWS`: no row independently satisfied paired same-individual diversity requirements.
- `QUOTA_UNAVAILABLE`: the environment exposed free space but no user-visible quota command, so the storage gate failed closed.

## Diversity, composition, and uncertainty

- Supported paper claim: the bounded pilot produced **no** promoted cross-species diversity or composition estimate. `validated_species_count = 0` is the exact result, not a missing-value placeholder.
- Supported paper claim: the only numeric mapping/queryability/callability evidence is the sentinel IMPG smoke artifact described above; it remains separate from biological pilot outputs.
- Unsupported paper claim: any coding diversity ratio, CDS/fourfold composition estimate, or class-level summary derived from this bounded pilot. None was authorized or computed.
- Unsupported paper claim: any population-genomic or demographic inference from VGP H1/H2 pairs. No callable diploid genotype set, no raw-read cohort, no VCF, and no phased demographic input was reviewed into existence here.

## Mapping, queryability, and callability behavior

- Global precondition behavior remained stable: exact-native sentinel mapping was 1:1 and queryable/callable denominators stayed positive in the handoff smoke artifact.
- Candidate-row behavior remained unresolved: every row in `analysis/vgp_pilot_paper_table.tsv` carries unresolved or missing `callable_bases`, `queryable_gene_count`, and `queryable_gene_bases` fields, so no row can support composition or diversity claims.
- This separation matters for the paper: sentinel assembly-engineering evidence supports toolchain integrity only; it does not support taxon-level biological inference.

## Independent Ne and ecological metadata inventory

- Usable independent Ne/life-history rows now: `0`. Population-data availability rows now: `0`.
- Supported paper claim: no independent Ne estimate, life-history covariate, or ecological covariate entered the reviewed pilot because the selected-species denominator was zero and the inventory was intentionally header-only.
- Supported paper claim: no circular predictor slipped in. The independent inventory explicitly rejected `pi/(4mu)`-style algebraic back-calculations, shared-sample genomic histories, and any claim that VGP H1/H2 pairs constitute callable diploid genotypes.
- Unsupported paper claim: that later PSMC, MSMC2, SMC++, or population-VCF work is already feasible for any reviewed pilot species. Every such method remains unsupported here because there are zero curated population-data rows.

## Supported and unsupported paper claims

### Supported

- The bounded VGP pilot remained a fail-closed refusal on July 17, 2026 UTC.
- Zero species were selected, zero species were executed, and zero diversity/composition outputs were promoted.
- All 74 reviewed vertebrate candidates failed exact eligibility for documented annotation/denominator and, where relevant, pairing reasons.
- Current executable cost under the frozen gate is zero across download, compute, memory, scratch, storage, inode, and I/O dimensions.

### Unsupported

- Any biological estimate for the bounded pilot beyond `validated_species_count = 0`.
- Any demographic or population-genomic claim from assembly pairs alone.
- Any use of algebraically derived Ne or overlapping genomic histories as independent predictors.
- Any interpretation that this synthesis authorizes another wave, a full catalog, bulk acquisition, or PSMC/MSMC2/SMC++ work.

## Resource calibration and cost consequences

- Current executable cost is exactly zero: download `0` GB, compute `0` core-h, wall `0` h, persistent input `0` GB, persistent output `0` GB, inodes `0`, read `0` GB, write `0` GB, metadata ops `0`.
- Contingent post-repair <=6-species wave proxy: `4.0332` / `8.5788` / `17.6184` core-h. Metrics that lacked executed-job telemetry remain explicitly labeled as unchanged planning envelopes in `analysis/vgp_pilot_next_decision.tsv`.
- Contingent post-repair reported full catalog proxy: `30.888` / `67.992` / `160.016` core-h. Storage and I/O for the overlap-deduplicated full catalog remain planning envelopes because the refusal pilot produced zero transfer telemetry.
- The key change from this bounded pilot is therefore negative but important: executable cost is now known to be zero under the present gate, while any non-zero future budget still depends on a repair step and a new explicit authorization packet.

## Terminal-state and retention verification

- Job/download terminal state from committed artifacts: `analysis/vgp_pilot_slurm_telemetry.tsv` records a single `run_summary` row with `final_state = NOT_SUBMITTED`, blank job identifiers, and `failure_code = GATE_NO_GO`.
- No pilot retrieval or staging manifest was populated, and no biological output tree crossed the boundary; the refusal artifacts themselves are the only retained pilot outputs.
- Retention and quarantine policy remains the reviewed one in `analysis/vertebrate_scaleout_decisions.tsv`: provenance, final refusal artifacts, failure ledgers, and telemetry are retained; unresolved/partial execution objects do not exist here and therefore there is nothing to quarantine or resume.

## Recommendation

- Choose `stop/repair`, not `another bounded expansion wave`, `full eligible-catalog consideration`, or `deferred`.
- Repair means resolving the source-count discrepancy, exposing a user-visible quota interface, and producing row-level native-annotation and denominator evidence before any new selection occurs.
- If a human later wants more work, it must be represented as a new explicit authorization task with numeric species scope and budgets analogous to `A50`/`A60`/`A70`/`A71`. This synthesis does not create or imply a ready executable node across that boundary.
