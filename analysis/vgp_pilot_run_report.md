# Repaired VGP pilot run report

- Run ID: `vgp-pilot-run-20260718T131009Z-32686`
- Gate decision: `NO_GO`
- Final state: `NOT_SUBMITTED`
- Failure: `GATE_NO_GO` — gate decision is NO_GO; repaired VGP compute is not authorized

The repaired authorization boundary closed before any Slurm or provider command. Exactly zero sbatch commands, jobs, compute starts, core-seconds, scratch bytes, I/O bytes, network bytes, provider requests, full-catalog downloads, population-bulk downloads, and demographic inferences occurred. No callable/queryable denominator or biological result was imputed.

## Bound and observed digests

- `gate_file_sha256`: `0c7b37078ce6997023923fa0fc0572e4cac2362316a9f223531ef257502730a6`
- `decision_sha256`: `6ba08bae6be0903b04cbf2b7a73c5201554ade82d2ea10a22280eff63289eb16`
- `authorization_tuple_digest`: `8c4aa1927cf48ae6fcee3587e271f4d729586511e26e23d17e237fbb637b3ee2`
- `manifest_sha256`: `b7cbe6cf22287d58dae9e270a81dee2970aaa3d8c73a2270323fa6210b276988`
- `acquisition_manifest_sha256`: `f44054e2a72e1aa2047ce50367c60b751ecf852a17fdc7075aef6d413c474c28`
- `immutable_inventory_sha256`: `f8d835008801cb32b0a4cc75990cfad4642d5e0f62bcb27df1c3814fd399a248`
- `root_config_sha256`: `912e47c66eeacaf8db3ff535fb857b296eb890d56e2ba0f38af3fdeb4ec12989`
- `input_bundle_digest`: `5595cb4cf38504fba4b6c862ee46283e8d48c7cafaeacbc813f5804390f135be`
- `root_contract_digest`: `05846dd869c05986e7a155be78fd9de03d4042a50cc488b34522621b6c600c2a`
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

- `run_manifest`: `/moosefs/erikg/lewontin-paradox/.wg-worktrees/agent-283/analysis/vgp_pilot_run_manifest.tsv` (`fb84735c8d4f9dda7b7f4a77cff4714802d3d067820a46ffb87f972e12027e75`)
- `telemetry`: `/moosefs/erikg/lewontin-paradox/.wg-worktrees/agent-283/analysis/vgp_pilot_slurm_telemetry.tsv` (`f7191a075de40d9c51f94092951657a0f2f0bd43841a5a31bca635768727cb76`)
- `results`: `/moosefs/erikg/lewontin-paradox/.wg-worktrees/agent-283/analysis/vgp_pilot_results.tsv` (`48de42e2c36df4f8fe27856117c17bca945b98244e5c11e2c5498da2c12f07f7`)
- `exclusions`: `/moosefs/erikg/lewontin-paradox/.wg-worktrees/agent-283/analysis/vgp_pilot_exclusions.tsv` (`7ad02cceff7b300504a947b7c1649d34849893c3ffbb380b8cd42dc435a62cde`)
- `refusals`: `/moosefs/erikg/lewontin-paradox/.wg-worktrees/agent-283/analysis/vgp_pilot_refusals.tsv` (`d098e9b516ca712c20f36485e26bd3d685ff64125a8ffb1b830594e01290d67a`)

The Slurm worker contract remains dormant. If a future independently regenerated exact GO passes every repaired binding and complete acquisition check, it requires pinned GNU Guix, node-local `$SLURM_TMPDIR`, no compute-node network, SweepGA whole-haplotype `--num-mappings 1:1`, IMPG native partition/query VCF/BCF contracts, atomic promotion, bounded retries, sentinels, dependency records, `sacct`, and scratch/I/O telemetry. VGP H1/H2 is never population or demographic input.
