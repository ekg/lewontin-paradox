# Vertebrate scale-out integration handoff

Date: 2026-07-17 UTC

WG task: `integrate-vertebrate-scaleout`

## Executive disposition

The required vertebrate planning lineage was already present on live `main`
before this task wrote any repository content. At audit time:

- `main = origin/main = 0337dabaf2d3e4288c37dc747cac26880d467661`
- the live root worktree was on `main...origin/main` with no tracked/staged
  changes; the only pre-existing user/unrelated delta was the untracked
  `.wg-worktrees/` directory in the root worktree
- the assigned integration worktree
  `wg/agent-167/integrate-vertebrate-scaleout` was clean

The five required upstream tasks were all completed and their required
artifacts were present exactly once on `main`, either as direct squash-merge
commits or as replayed-equivalent documentation commits. No missing planning
artifact required a content merge from another branch.

This task therefore performs a documentation-only integration on top of the
already integrated planning tree by adding this handoff file. No biological
asset download, Slurm submission, pilot execution, or demographic inference was
performed.

## Starting-state audit

### Live main worktree

- Path: `/moosefs/erikg/lewontin-paradox`
- Branch state: `main...origin/main`
- Pre-existing user/unrelated delta: `?? .wg-worktrees/`
- Tracked/staged delta: none

### Relevant task worktrees

These were inspected for residual user/task-owned state and were not modified:

| worktree | observed status |
|---|---|
| `agent-135` (`rerun-sync-origin`) | branch upstream shown as `[gone]`; only `?? .wg-cleanup-pending` |
| `agent-153` (`plan-vertebrate-ne-strategy`) | branch upstream shown as `[gone]`; only `?? .wg-cleanup-pending` |
| `agent-155` (`plan-all-vertebrate-tier3`) | branch upstream shown as `[gone]`; only `?? .wg-cleanup-pending` |
| `agent-160` (`synthesize-vertebrate-scaleout-plan`) | branch upstream shown as `[gone]`; only `?? .wg-cleanup-pending` |
| `agent-118` (`.quality-pass-vertebrate-scaleout-plan`) | worktree already cleaned up; path absent |

No tracked user changes were present in those worktrees. This task did not
reset, rebase, clean, amend, or otherwise alter them.

## Dependency-correct integration ledger

### Required upstream tasks and exact task-owned commits

| task | task-owned branch commit(s) | integrated `main` commit(s) | artifact(s) on final tree | proof |
|---|---|---|---|---|
| `rerun-sync-origin` | `00f72fc86150404bba0712b2672c320077c7b962` | `234b9ca97d63964d130bbe0f9c0ca1bb6eb5d51f` | `analysis/vertebrate_scaleout_origin_inventory.md` | `git diff --quiet 00f72fc 234b9ca -- analysis/vertebrate_scaleout_origin_inventory.md` exited `0` |
| `.quality-pass-vertebrate-scaleout-plan` | `70c17e73e0874151cf13ed935c7daef6e96d46e1` | `d586f9e8d33bb9583a2aa94ddeca3d7a7efd53fc`, then `ee383af1c9c65cc563afc695d9268397d3ba0284` | `analysis/vertebrate_scaleout_quality_pass.md` | `git diff --quiet 70c17e7 d586f9e -- analysis/vertebrate_scaleout_quality_pass.md` exited `0`; `git diff --quiet ee383af 0337dab -- analysis/vertebrate_scaleout_quality_pass.md` exited `0` |
| `plan-all-vertebrate-tier3` | `1f4e86a870fb79e9a922e622ea67e26f8ca4a714` | `11ee6657637ad48065f189c27e7c5e2e2b7e3d40` | `analysis/vertebrate_scaleout_candidate_schema.tsv`, `analysis/vertebrate_scaleout_plan.md`, `analysis/vertebrate_scaleout_resource_budget.tsv`, `analysis/vertebrate_scaleout_wg_graph.md` | `git diff --quiet 1f4e86a 11ee665 -- ...` exited `0` on all four paths |
| `plan-vertebrate-ne-strategy` | `f728e73cb434d2e68723913ed5c4d4989b744458` | `e0bfd66e6452a007b1d90d7c613b149b8f868cb6` | `analysis/vertebrate_ne_execution_graph.tsv`, `analysis/vertebrate_ne_method_matrix.tsv`, `analysis/vertebrate_ne_source_schema.tsv`, `analysis/vertebrate_ne_strategy.md` | `git diff --quiet f728e73 e0bfd66 -- ...` exited `0` on all four paths |
| `synthesize-vertebrate-scaleout-plan` | `7f10182045ec4c45012b275a4d17eb21caa94436` | `0337dabaf2d3e4288c37dc747cac26880d467661` | `analysis/vertebrate_scaleout_decisions.tsv`, `analysis/vertebrate_scaleout_execution_graph.tsv`, `analysis/vertebrate_scaleout_execution_plan.md` | `git diff --quiet 7f10182 0337dab -- ...` exited `0` on all three paths |

### Pre-handoff ancestry already on `main`

The integration-critical ancestry path present on `main` before this handoff
commit was:

1. `14fe57d1d0593a2859672279be60aa858e822410`  
   upstream VGP guidance, `results/tier3/vgp_freeze_analysis.md`
2. `234b9ca97d63964d130bbe0f9c0ca1bb6eb5d51f`  
   squash merge of `rerun-sync-origin`
3. `d586f9e8d33bb9583a2aa94ddeca3d7a7efd53fc`  
   replay of the quality-pass audit artifact
4. `ee383af1c9c65cc563afc695d9268397d3ba0284`  
   synchronized planning handoff update to the quality-pass artifact
5. `e0bfd66e6452a007b1d90d7c613b149b8f868cb6`  
   squash merge of `plan-vertebrate-ne-strategy`
6. `11ee6657637ad48065f189c27e7c5e2e2b7e3d40`  
   squash merge of `plan-all-vertebrate-tier3`
7. `0337dabaf2d3e4288c37dc747cac26880d467661`  
   squash merge of `synthesize-vertebrate-scaleout-plan`

That ancestry proves the synchronized guidance commit, the rerun inventory, the
quality audit, both planner outputs, and the synthesis outputs were all already
on live `main` before this task added the handoff artifact.

## Final-tree artifact inventory

The following required planning and guidance artifacts were confirmed on the
integrated tree with `git ls-tree -r --long HEAD`:

| path | blob | size (bytes) |
|---|---|---:|
| `analysis/vertebrate_scaleout_origin_inventory.md` | `cc2dfddeb2f9ab85953ac4ddfd690a9e0a39ccd4` | 21840 |
| `analysis/vertebrate_scaleout_quality_pass.md` | `00876114ad1df14a3de015dd20c5de48d0188af9` | 9713 |
| `analysis/vertebrate_scaleout_candidate_schema.tsv` | `bb5f0b801e02d1955a0df374ae8debbd74c33731` | 23917 |
| `analysis/vertebrate_scaleout_plan.md` | `e249fc8e7014947bb37dc3cc52373e716e3ad4d3` | 49139 |
| `analysis/vertebrate_scaleout_resource_budget.tsv` | `1e57a0c9b6d2154f0686531b9adcd44366445b4f` | 12362 |
| `analysis/vertebrate_scaleout_wg_graph.md` | `1c151211b125df89fc20d20bc66c65a8290c2f5a` | 7444 |
| `analysis/vertebrate_ne_execution_graph.tsv` | `2dba1895d33acdfe1407cf9cf328eb560179cdb9` | 7563 |
| `analysis/vertebrate_ne_method_matrix.tsv` | `10e1a61adca96b4fbc0671fa08e0b5c0d77ec311` | 23811 |
| `analysis/vertebrate_ne_source_schema.tsv` | `4de200d4d66d8df006f6ae8ed70889c1cc58d448` | 15662 |
| `analysis/vertebrate_ne_strategy.md` | `126c0823a1a78030e26de127879ce26bd8c55a89` | 51936 |
| `analysis/vertebrate_scaleout_decisions.tsv` | `6311b81c80a509eec5b34a2f6be5b78203c68239` | 18255 |
| `analysis/vertebrate_scaleout_execution_graph.tsv` | `b78d40ae43daa69f6aca0594dc4b64242f9e42b0` | 34819 |
| `analysis/vertebrate_scaleout_execution_plan.md` | `704b04b6c21ce5e49969b594e2f74f105016e3b6` | 48357 |
| `results/tier3/vgp_freeze_analysis.md` | `d37faec861ad100685d3a69fa0b68c8bd7f3a7bc` | 4481 |

No required artifact was missing, and no task-owned artifact needed a second
copy under an alternate path.

## Validation evidence

### Git ancestry and equivalence checks

All five dependency checks named in the ledger above were executed directly on
the repository and exited `0`, proving that the task-owned branch commits were
already represented on `main` either by exact squash-equivalent content or by a
replayed/finalized documentation path.

Additional ancestry observations:

- `main = origin/main = 0337dabaf2d3e4288c37dc747cac26880d467661` before this
  task's commit
- the path `dec7266..0337dab` contains the replayed quality-pass commits and
  both planner squash merges in dependency order
- no tracked worktree or index change from another user/agent was overwritten

### Pinned GNU Guix foundational smoke

Command:

```bash
./analysis/validate_tier3_guix.sh
```

Observed result:

- reproducible profile path:
  `/gnu/store/z9v2f6faha9cwjz0sm5iphhlzisgi077-profile`
- scientific imports reported expected versions:
  `Bio 1.80`, `jsonschema 4.5.1`, `numpy 1.23.2`, `pandas 1.4.4`,
  `pyfaidx 0.7.2.1`, `pysam 0.20.0`, `scipy 1.10.1`
- foundational smoke tests: `19 passed in 0.37s`
- required toolchain surfaced in the pure profile: `bcftools 1.14`,
  `samtools 1.14`, `bgzip 1.16`, `bedtools 2.30.0`, `VCFtools 0.1.16`
- the pure profile did **not** leak an executable `impg`

### Full proportionate analysis suite

Successful command:

```bash
ROOT=$(pwd)
guix time-machine -C analysis/guix/channels.scm -- \
  shell -L analysis/guix -m analysis/guix/manifest.scm --pure -- \
  bash --noprofile --norc -c 'cd "$1" && python3 -m pytest -q analysis/tests' _ "$ROOT"
```

Observed result:

- `153 passed in 2.86s`

Notes:

- an initial `bash -lc` wrapper attempt failed before tests ran because
  `/etc/profile` inside the pure shell dropped `id` and `pytest` from the
  runtime command path
- no project test failed; the corrected `--noprofile --norc` command above is
  the authoritative full-suite validation result

## Authorization boundary

This task did **not** authorize or perform any execution beyond repository
audit, documentation, and validation:

- no VGP biological asset was downloaded
- no Slurm command, `sbatch`, `srun`, or scheduler submission was run
- no pilot, expansion wave, or full-catalog execution was launched
- no demographic inference was run

The execution boundary defined by the integrated planning artifacts remains in
force. Expensive actions are still structurally blocked behind explicit human
authorization nodes in:

- `analysis/vertebrate_scaleout_execution_graph.tsv`
- `analysis/vertebrate_scaleout_decisions.tsv`
- `analysis/vertebrate_ne_execution_graph.tsv`

## Integration conclusion

The dependency-correct vertebrate planning lineage was already integrated on
live `main` before this task wrote the handoff. The only repository mutation
needed here is this handoff artifact itself, which documents:

- exact upstream task-owned commits
- exact `main` replay/squash commits
- artifact presence on the integrated tree
- worktree-safety and unrelated-change preservation
- pinned Guix smoke and full-suite results
- the still-active authorization boundary preventing premature execution
