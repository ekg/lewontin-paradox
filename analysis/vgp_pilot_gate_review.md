# Independently regenerated VGP pilot gate review

- Decision: `NO_GO`
- Decision SHA-256: `9f39b13be5e0b1999c4cd98498399aee8700a7487d57fa81dfb7c59c29ff867d`
- Authorization tuple SHA-256: `410c42c2894d69d88b70361ebe20a3c806e203524595be746562be97343f9240`
- Downstream rule: only the literal decision `GO`, with every bound digest reverified, authorizes acquisition or compute.

## Catalog units and duplicate identities

- 717 physical lines = `1` header + `716` data rows.
- `714` unique species; the data-row excess over unique species is `2`.
- `Lophostoma evotis` multiplicity: `2`.
- `Micronycteris microtis` multiplicity: `2`.

## Independent row closure

- Manifest rows audited: `6`.
- Exact H1/native-annotation metadata-ready rows: `6`.
- Exact paired Tier3A-ready rows: `0`.
- Seed/rejection rows independently closed: `74` / `74`; all match: `true`.
- Independently expected selected rows: `none`.

## Strictest finite cap vector

- `species`: observed `6.0`; limit `6.0` count; pass `true`; winner `regate_task_bound`.
- `compressed_inputs_gib`: observed `3.8653049999999998`; limit `120.0` GiB; pass `true`; winner `regate_task_bound`.
- `scratch_gib`: observed `205.5066`; limit `139.69838619232178` GiB; pass `false`; winner `integrated_execution_plan_small_cap`.
- `core_hours`: observed `11.1353`; limit `280.0` core-hours; pass `true`; winner `integrated_execution_plan_small_cap`.
- `concurrent_species`: observed `2.0`; limit `2.0` count; pass `true`; winner `regate_task_bound`.
- `memory_per_job_gib`: observed `32.0`; limit `96.0` GiB; pass `true`; winner `integrated_execution_plan_small_cap`.
- `file_inodes`: observed `22596.0`; limit `60000.0` count; pass `true`; winner `integrated_execution_plan_small_cap`.
- `moosefs_read_gb`: observed `406.72740000000005`; limit `180.0` GB; pass `false`; winner `integrated_execution_plan_small_cap`.
- `moosefs_write_gb`: observed `3.1634`; limit `200.0` GB; pass `true`; winner `integrated_execution_plan_small_cap`.
- `metadata_operations`: observed `271151.0`; limit `500000.0` count; pass `true`; winner `integrated_execution_plan_small_cap`.
- `peak_bandwidth_mib_s`: observed `120.0`; limit `120.0` MiB/s; pass `true`; winner `integrated_execution_plan_small_cap`.
- `aggregate_wall_hours`: observed `8.7311`; limit `17.5` hours; pass `true`; winner `integrated_execution_plan_small_cap`.
- `persistent_input_gb`: observed `4.150340415`; limit `160.0` GB; pass `true`; winner `integrated_execution_plan_small_cap`.
- `persistent_output_gb`: observed `0.0`; limit `32.0` GB; pass `true`; winner `integrated_execution_plan_small_cap`.
- `cpus_per_element`: observed `8.0`; limit `8.0` count; pass `true`; winner `integrated_execution_plan_small_cap`.

## External-root and storage contract

- Root: `/moosefs/erikg/lewontin-paradox-data/vgp/phase1-freeze-1.0`; live identity digest is bound in the decision.
- Filesystem free bytes/inodes: `101224219672576` / `1004580301`.
- Enforceable allocation/quota status: `unknown`.
- Filesystem free space and inodes are reported independently from enforceable per-user/allocation limits. Unknown quota never counts as adequate and cannot be overridden by free space.
- Required worst-case byte and inode capacity includes at least 25% headroom; every stricter integrated limit wins.

## Retrieval, checksums, pairs, and denominators

- Pre-download-ready exact-version rows: `6`.
- Each exact official payload is staged, finite-size checked, checked against every available official checksum, locally SHA-256 hashed, reverified, and atomically promoted read-only. A missing remote SHA-256/MD5 is not itself a pre-download blocker.
- Pair evidence digest: `73e641104bad3f5a5f7c269a303e72dbd492141b9be3a05d6b60c0d143f1f7b3`; Tier3A requires exact versioned H2 plus affirmative same-individual and phasing evidence.
- Denominator contract: `vgp_post_alignment_denominators_v1`; pre-download prerequisite `false`.
- Callable/queryable denominators are measured after alignment. Missing or sub-threshold measurements exclude the affected downstream result.

## Bound decision components

- `catalog_provenance_digest`: `dc3fdcaea438635595a4d1c50b2d3146846101a154e6bbc1406a5da077cae58b`
- `data_root_storage_contract_digest`: `1f761b8a8df248f50ea8c61d98189c15113687439ecf81ca519696c91005f7fb`
- `environment_digest`: `fe8e45bb716ffc9560cf2cda37d4974c46905279635ca5e8c590666d1f9b354e`
- `cap_vector_digest`: `b6d15c85feee41edac606dba1d51e5d288016de0b99b0922292d0244437bac4c`
- `retrieval_checksum_obligations_digest`: `3e8ef805ded1ac58f1e1446d96a45922a40b5ddc34fff6424fd183be1cfe5be9`
- `pair_evidence_digest`: `73e641104bad3f5a5f7c269a303e72dbd492141b9be3a05d6b60c0d143f1f7b3`
- `measurement_contract_digest`: `6fa68a0a657877953d5d92dc99f850492fbe5d823672d38d81cd806653f84018`
- `row_dispositions_digest`: `dc6357145a61b4692d8db30b0e3a19dcde683281b9a3abfa7dd59bf887eb90fe`
- `manifest_digest`: `b7cbe6cf22287d58dae9e270a81dee2970aaa3d8c73a2270323fa6210b276988`
- `root_contract_digest`: `1f761b8a8df248f50ea8c61d98189c15113687439ecf81ca519696c91005f7fb`
- `input_bundle_digest`: `5595cb4cf38504fba4b6c862ee46283e8d48c7cafaeacbc813f5804390f135be`
- `authorization_tuple_digest`: `410c42c2894d69d88b70361ebe20a3c806e203524595be746562be97343f9240`

## Blockers

- `CAP_MOOSEFS_READ_GB_EXCEEDED`: proposed finite worst-case moosefs_read_gb=406.72740000000005 exceeds strictest limit 180.0 GB (strict cap vector)
- `CAP_SCRATCH_GIB_EXCEEDED`: proposed finite worst-case scratch_gib=205.5066 exceeds strictest limit 139.69838619232178 GiB (strict cap vector)
- `QUOTA_UNAVAILABLE`: No user-visible quota command was available from the current environment. (analysis/vgp_data_root_validation.json)

## Authorization exclusions

- Full-catalog acquisition is unauthorized.
- Raw population bulk download is unauthorized.
- Demographic inference is unauthorized.
- This gate launched zero downloads and zero jobs.
