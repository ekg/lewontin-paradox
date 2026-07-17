# VGP pilot acquisition report

- Run ID: `vgp-pilot-acquire-20260717T192115Z`
- Generated at: `2026-07-17T19:21:15Z`
- Gate path: `/moosefs/erikg/lewontin-paradox/.wg-worktrees/agent-190/analysis/vgp_pilot_gate.json`
- Gate decision: `NO_GO`
- Acquisition status: `refused_preflight`
- Refused before first biological byte: `true`
- Slurm environment detected: `false`
- Biological payload bytes transferred: `0`
- Verified immutable objects promoted: `0`
- Quarantine objects written: `0`
- Accessions promoted into views: `0`
- Output manifest: `/moosefs/erikg/lewontin-paradox/.wg-worktrees/agent-190/analysis/vgp_pilot_acquisition_manifest.tsv`

## Authorization Boundary

- Manifest path: `/moosefs/erikg/lewontin-paradox/.wg-worktrees/agent-190/analysis/vgp_pilot_manifest.tsv`
- Recorded manifest digest: `f27b81f369af18caf97a1ccee1b14d8a8b050d665832e6e99f188041275f49ce`
- Root contract path: `/moosefs/erikg/lewontin-paradox/.wg-worktrees/agent-190/analysis/vgp_data_root_config.json`
- Recorded root contract digest: `912e47c66eeacaf8db3ff535fb857b296eb890d56e2ba0f38af3fdeb4ec12989`
- Recorded cap vector digest: `d51827af5600976f52549ce1a16a7f38bef03f671814d351ff2270ebe7b1b59e`

## Refusal Reason

- `GATE_NO_GO`: gate decision is NO_GO; acquire is not authorized for this manifest/root/cap vector

## Gate Blockers

- `SOURCE_COUNT_DISCREPANCY_UNRESOLVED`: the frozen raw VGP catalog still disagrees with the planning headline counts and no explicit signed discrepancy resolution is bundled here (analysis/vgp_phase1_freeze_provenance.json)
- `NO_SELECTED_ROWS`: the frozen pilot manifest selects zero rows, so no bounded pilot can be authorized (analysis/vgp_pilot_manifest.tsv)
- `ZERO_COMPOSITION_ELIGIBLE_ROWS`: no manifest row independently satisfies exact H1/native-annotation/denominator requirements (analysis/vgp_pilot_manifest.tsv)
- `ZERO_DIVERSITY_ELIGIBLE_ROWS`: no manifest row independently satisfies the paired same-individual diversity requirements (analysis/vgp_pilot_manifest.tsv)
- `QUOTA_UNAVAILABLE`: No user-visible quota command was available from the current environment. (analysis/vgp_data_root_validation.json)

## Validation Notes

- Live gate authorization was executed against the current manifest, root contract, and recomputed cap vector.
- The current gate is not exactly `GO`, so no provider request, staging write, partial file, checksum verification, quarantine action, or promotion step was attempted.
- Exact-reference/native-annotation linkage validation under pinned GNU Guix was not re-run because zero assets were authorized or acquired.
- No Slurm command was submitted by this task.
- No biological bulk file was added to Git by this task.
