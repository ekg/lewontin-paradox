# Repaired VGP pilot run report

- Run ID: `vgp-pilot-run-20260718T090728Z-14273`
- Gate decision: `NO_GO`
- Final state: `NOT_SUBMITTED`
- Failure: `GATE_NO_GO` — gate decision is NO_GO; repaired VGP compute is not authorized

The repaired authorization boundary closed before any Slurm or provider command. Exactly zero sbatch commands, jobs, compute starts, core-seconds, scratch bytes, I/O bytes, network bytes, provider requests, full-catalog downloads, population-bulk downloads, and demographic inferences occurred. No callable/queryable denominator or biological result was imputed.

## Bound and observed digests

- `gate_file_sha256`: `fc284eb617aadc6af4d348beaf47377df55fc4798c5e1460583dc6b9ef7115d6`
- `decision_sha256`: `9f39b13be5e0b1999c4cd98498399aee8700a7487d57fa81dfb7c59c29ff867d`
- `authorization_tuple_digest`: `410c42c2894d69d88b70361ebe20a3c806e203524595be746562be97343f9240`
- `manifest_sha256`: `b7cbe6cf22287d58dae9e270a81dee2970aaa3d8c73a2270323fa6210b276988`
- `acquisition_manifest_sha256`: `55c8270abe62d6508898f89e4f66b24f338ede8fbc5d8aa9b1593104ce581d88`
- `immutable_inventory_sha256`: `f8d835008801cb32b0a4cc75990cfad4642d5e0f62bcb27df1c3814fd399a248`
- `root_config_sha256`: `912e47c66eeacaf8db3ff535fb857b296eb890d56e2ba0f38af3fdeb4ec12989`
- `input_bundle_digest`: `5595cb4cf38504fba4b6c862ee46283e8d48c7cafaeacbc813f5804390f135be`
- `root_contract_digest`: `1f761b8a8df248f50ea8c61d98189c15113687439ecf81ca519696c91005f7fb`
- `cap_vector_digest`: `b6d15c85feee41edac606dba1d51e5d288016de0b99b0922292d0244437bac4c`
- `retrieval_digest`: `3e8ef805ded1ac58f1e1446d96a45922a40b5ddc34fff6424fd183be1cfe5be9`
- `pair_evidence_digest`: `73e641104bad3f5a5f7c269a303e72dbd492141b9be3a05d6b60c0d143f1f7b3`
- `measurement_contract_digest`: `6fa68a0a657877953d5d92dc99f850492fbe5d823672d38d81cd806653f84018`
- `environment_digest`: `fe8e45bb716ffc9560cf2cda37d4974c46905279635ca5e8c590666d1f9b354e`
- `sweepga_build_sha256`: `ead0fc6bf770f573a3ca8f2532f04c57a4653a9248e31963b564971221f78e36`
- `impg_handoff_sha256`: `ff088828a3dc7336beec4e727fa7598840ec0f8b1234758f64f4480319364b61`
- `worker_sha256`: `5123e25783c511cacdbd43425365a0f02131c29e397cb091040fbdf31dead90e`

## Gate blockers

- `CAP_MOOSEFS_READ_GB_EXCEEDED` — proposed finite worst-case moosefs_read_gb=406.72740000000005 exceeds strictest limit 180.0 GB
- `CAP_SCRATCH_GIB_EXCEEDED` — proposed finite worst-case scratch_gib=205.5066 exceeds strictest limit 139.69838619232178 GiB
- `QUOTA_UNAVAILABLE` — No user-visible quota command was available from the current environment.

## Promoted evidence tables

- `run_manifest`: `/moosefs/erikg/lewontin-paradox/.wg-worktrees/agent-230/analysis/vgp_pilot_run_manifest.tsv` (`82e8316b91acf7d7645e2381232b181d76ff7b17c919eca5696e2a3417a02f2f`)
- `telemetry`: `/moosefs/erikg/lewontin-paradox/.wg-worktrees/agent-230/analysis/vgp_pilot_slurm_telemetry.tsv` (`a1b94966d31c615628069f455a1412210c7b58e94fc6d2576e60fedcb530ef5e`)
- `results`: `/moosefs/erikg/lewontin-paradox/.wg-worktrees/agent-230/analysis/vgp_pilot_results.tsv` (`4c4b27babd0603986ecb801964abc47c207e27264ab2636cb9dabc9aada011ef`)
- `exclusions`: `/moosefs/erikg/lewontin-paradox/.wg-worktrees/agent-230/analysis/vgp_pilot_exclusions.tsv` (`a503f143b00b35c86419478883e7a3879f387d0ab6cfa59dbd4e7041fbac5839`)
- `refusals`: `/moosefs/erikg/lewontin-paradox/.wg-worktrees/agent-230/analysis/vgp_pilot_refusals.tsv` (`569e96eb35f3c6c05d048aa6a48b011f9a38544c5335fa18019a0213c77b3407`)

The Slurm worker contract remains dormant. If a future independently regenerated exact GO passes every repaired binding and complete acquisition check, it requires pinned GNU Guix, node-local `$SLURM_TMPDIR`, no compute-node network, SweepGA whole-haplotype `--num-mappings 1:1`, IMPG native partition/query VCF/BCF contracts, atomic promotion, bounded retries, sentinels, dependency records, `sacct`, and scratch/I/O telemetry. VGP H1/H2 is never population or demographic input.
