# Repaired VGP pilot acquisition report

- Run ID: `vgp-pilot-acquire-20260718T085124Z-22120`
- Generated at: `2026-07-18T08:51:24Z`
- Gate path: `/moosefs/erikg/lewontin-paradox/.wg-worktrees/agent-227/analysis/vgp_pilot_gate.json`
- Gate decision: `NO_GO`
- Acquisition status: `refused_preflight`
- Refused before first biological byte: `true`
- Provider requests attempted: `0`
- Biological payload bytes transferred: `0`
- Verified immutable objects promoted: `0`
- Quarantine objects written: `0`
- Slurm environment detected: `false`
- Output manifest: `/moosefs/erikg/lewontin-paradox/.wg-worktrees/agent-227/analysis/vgp_pilot_acquisition_manifest.tsv`
- Immutable-object inventory: `/moosefs/erikg/lewontin-paradox/.wg-worktrees/agent-227/analysis/vgp_pilot_immutable_object_inventory.tsv`
- Refusal evidence: `/moosefs/erikg/lewontin-paradox/.wg-worktrees/agent-227/analysis/vgp_pilot_acquisition_refusal.json`

## Authorization boundary

- Manifest path: `/moosefs/erikg/lewontin-paradox/.wg-worktrees/agent-227/analysis/vgp_pilot_manifest.tsv`
- Manifest digest: `b7cbe6cf22287d58dae9e270a81dee2970aaa3d8c73a2270323fa6210b276988`
- Root contract path: `/moosefs/erikg/lewontin-paradox/.wg-worktrees/agent-227/analysis/vgp_data_root_config.json`
- Root/storage digest: `1f761b8a8df248f50ea8c61d98189c15113687439ecf81ca519696c91005f7fb`
- Environment digest: `fe8e45bb716ffc9560cf2cda37d4974c46905279635ca5e8c590666d1f9b354e`
- Cap-vector digest: `b6d15c85feee41edac606dba1d51e5d288016de0b99b0922292d0244437bac4c`
- Retrieval/checksum digest: `3e8ef805ded1ac58f1e1446d96a45922a40b5ddc34fff6424fd183be1cfe5be9`
- Pair-evidence digest: `73e641104bad3f5a5f7c269a303e72dbd492141b9be3a05d6b60c0d143f1f7b3`
- Measurement-contract digest: `6fa68a0a657877953d5d92dc99f850492fbe5d823672d38d81cd806653f84018`
- Authorization-tuple digest: `410c42c2894d69d88b70361ebe20a3c806e203524595be746562be97343f9240`

## Refusal reason

- `GATE_NO_GO`: gate decision is NO_GO; acquire is not authorized

## Gate blockers

- `CAP_MOOSEFS_READ_GB_EXCEEDED`: proposed finite worst-case moosefs_read_gb=406.72740000000005 exceeds strictest limit 180.0 GB (strict cap vector)
- `CAP_SCRATCH_GIB_EXCEEDED`: proposed finite worst-case scratch_gib=205.5066 exceeds strictest limit 139.69838619232178 GiB (strict cap vector)
- `QUOTA_UNAVAILABLE`: No user-visible quota command was available from the current environment. (analysis/vgp_data_root_validation.json)

## Validation notes

- The gate's nested authorization tuple was verified before the downloader boundary; a literal GO would additionally trigger live recomputation of every bound input and digest.
- Literal GO, exact selected rows/assets, strictest cap values, 25% storage headroom, and explicit scope exclusions are mechanical preconditions.
- Every authorized payload uses resumable `.part` staging, finite-size and cumulative-byte checks, any official checksum, local SHA-256, immediate SHA-256 reverification, and same-filesystem atomic read-only promotion.
- Mismatches are quarantined and excluded from the immutable inventory.
- Exact-reference/native-annotation/Tier3A-pair checks are gate-bound to the pinned GNU Guix environment; on refusal no staged biological report exists to re-run.
- Exact-reference/native-annotation linkage validation under pinned GNU Guix was not re-run because zero assets were authorized or acquired.
- No full catalog or raw population bulk acquisition, Slurm submission, or demographic inference is performed by this entrypoint.
