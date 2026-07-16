# Tier 3 recovery integration: final handoff

Status: **PASS — the corrected Tier 3A lineage, powered Tier 3B lineage, recovery synthesis, manuscript, figures, evidence ledger, numerical audit, and artifact index are integrated on `main`, validated from the integrated tree, and pushed without rewriting history.**

Validation date: 2026-07-16 UTC  
Integration task: `integrate-tier3-recovery-main`  
Repository remote: `git@github.com:ekg/lewontin-paradox.git`

## Starting repository state and safety boundary

The integration audit began after `git fetch origin --prune` with:

| item | starting value |
|---|---|
| main HEAD | `043f460fb2eda6f53bd28b4afbf4e7bce8fd760d` |
| origin/main | `043f460fb2eda6f53bd28b4afbf4e7bce8fd760d` |
| task-worktree HEAD | `043f460fb2eda6f53bd28b4afbf4e7bce8fd760d` |
| task-worktree status | clean |
| main-worktree status | tracked files clean; untracked `.wg-worktrees/` administrative directory present |

The untracked main-worktree directory was inspected rather than removed. It contains the WG-managed linked worktrees and `.merge-lock`, not an uncommitted repository result. No reset, clean, stash, forced checkout, destructive restore, force push, or deletion of a dirty worktree was performed.

## Actual ancestry and integration decision

The four named scientific task tips are independent commit objects: none of `43bb47c`, `dd23c1b`, `03d2259`, or `ced8329` is an ancestor of any of the other three. In particular, the synthesis task tip does **not** contain the other named task tips as Git ancestors.

WG had already integrated each validated delta onto `main` as a single-parent task commit before this final release task began. Patch-range identity and exact tip/integration tree identity were checked, rather than inferred from subjects:

| dependency order | validated task lineage | validated tip | main integration commit | evidence |
|---:|---|---|---|---|
| 1 | repaired/powered Tier 3B prerequisite | `205dae93d716e5819b1626c66e7374e7c9827970` | `ba693e5e4b2f4453f80bf80c0b0740e7047b2235` | task-range and main-commit patch IDs equal |
| 2 | `run-tier3b-biological-recovery` | `03d2259d0060a44d837125ad2559b07d1dda7df6` plus prerequisites | `55712f4469b45eae1358cdb253e3c291e9630456` | exact trees equal: `84da08c2edd85e5f8bbb3fd7f526d3b634dc3f7b` |
| 3 | origin/main SweepGA build prerequisite | `4cb38abf1dddb05828d1b46a65fc72eca87fd185` | included as equivalent commit `f9e11ee0f7fd57632b05c05a61fc42a76bad9bf4` in the remap range | stable patch ID `f91970dbd49288b167cffa05e13d46d784805b5f` for both commits |
| 4 | `remap-tier3a-origin-sweepga` | `43bb47c15787770328100d98b131a62a8469bfc9` and its five-commit range from `55712f4` | `1f7f61253562fc1d815bff9bae3bd56477309241` | exact trees equal: `8318c545567231a86079adc5d53f3cc3c2574967` |
| 5 | `rerun-tier3a-origin-sweepga` | `dd23c1b9a51eb07b38b34ea91c614f2f5765a4d1` | `bc2e846846722d1274abfce6d89b8303c3b6c44e` | exact trees equal: `cda4ba6dc234e5a222cf527dc4fccb565986d8b9` |
| 6 | `synthesize-tier3-diversity-recovery` | `ced8329a887ce1c96ef1aaff8ccffbfe815dcb6a` | `043f460fb2eda6f53bd28b4afbf4e7bce8fd760d` | exact trees equal: `1b055bfbabd076899181d37a559df7f5083ae3cc` |

The minimal safe integration decision was therefore to retain the existing main lineage and add only this final handoff. Re-merging or cherry-picking any named scientific range would have duplicated already integrated patches. No conflict prerequisite was needed.

The resulting first-parent scientific integration chain is:

```text
ba693e5  repair-tier-3b
  -> 55712f4  run-tier3b-biological-recovery
  -> 1f7f612  remap-tier3a-origin-sweepga
  -> bc2e846  rerun-tier3a-origin-sweepga
  -> 043f460  synthesize-tier3-diversity-recovery
```

## Corrected Tier 3A provenance gate

All three accepted H1/H2 tuples trace to:

- fetched SweepGA origin/main source commit `018e4ce49d2c125820e0ac50dc5feaa02d423683`;
- reproducibly built SweepGA binary SHA-256 `fa7f0edb9b7e275c288db254046020e136d4267dd5ee043379227ef80da0573b`;
- pinned Guix closure `/gnu/store/z9v2f6faha9cwjz0sm5iphhlzisgi077-profile`;
- native commands containing the exact option `--num-mappings 1:1`;
- observed native query and target multiplicities of one;
- fresh IMPG partition/query outputs and H1-native annotation panels.

`results/tier3a/diploid_lineage_audit.json` reports `status: passed` and empty arrays for all three fail-closed fields:

```text
old_result_checksum_intersection = []
superseded_checksum_intersection = []
superseded_path_intersection = []
```

The machine audit additionally hashed every selected Tier 3A artifact and verified that no selected checksum intersects `results/tier3a/diploid_superseded_results.tsv`. Superseded paths and hashes remain only in explicitly labeled negative-control/supersession ledgers; they are absent from the active result lineage and synthesis inputs.

## Headline biological evidence

The independent audit recomputed 23/23 published headline values from committed numerators and denominators at `1e-12` relative/absolute tolerance.

### Tier 3A: alignment-conditioned H1/H2 coding panels

| assembly pair | coding-gene diversity (95% CI; variants/callable) | CDS diversity (95% CI; variants/callable) | reference-conditioned pi_S/pi_W (conservative interval) |
|---|---|---|---|
| *Menidia menidia* | 0.0138517 [0.0123545, 0.0153404]; 28,233/2,038,234 | 0.00674737 [0.00589925, 0.00776835]; 2,044/302,933 | 0.672544 [0.482732, 0.932999] |
| *Spinachia spinachia* | 0.000397957 [0.000329512, 0.000481235]; 803/2,017,806 | 0.000232752 [0.000148300, 0.000338608]; 63/270,674 | 0.316941 [0.0685349, 1.70202] |
| *Tautogolabrus adspersus* | 0.00167799 [0.00123835, 0.00212181]; 3,551/2,116,227 | 0.000515957 [0.000344761, 0.000713813]; 115/222,887 | 0.760295 [0.278040, 2.12518] |

These are selected, alignment-conditioned assembly comparisons, not population diversity or genome-wide deposited-individual heterozygosity. Uncertainty for direct estimates is a 1,000-replicate 50-kb genomic block bootstrap. The ratio bounds conservatively divide marginal component intervals.

### Tier 3B: powered population estimates

| population | population pi (95% CI; pairwise numerator/callable) | reference-conditioned pi_S/pi_W (95% block-bootstrap CI) |
|---|---|---|
| `AO_Luanda_2009_coluzzii` | 0.01235003 [0.01186339, 0.01281725]; 184,914.842/14,972,821 | 0.7426563 [0.7129171, 0.7708752] |
| `GM_WaliKunda_2012_coluzzii` | 0.01586979 [0.01539494, 0.01633848]; 237,893.579/14,990,338 | 0.7235550 [0.6919446, 0.7571679] |

Each population contains 20 wild diploid individuals. Population pi and ratio intervals use 10,000 chromosome-stratified 1-Mb block-bootstrap replicates; component pi_S and pi_W intervals use a 20-unit delete-one-individual jackknife.

### Claim boundary and remaining scientific limitations

- Tier 3A contains three selected coding-panel assembly pairs, not a population-size series.
- Tier 3B contains two populations from one species, so it measures within-species population heterogeneity rather than an across-species slope.
- Reference-conditioned pi_S/pi_W is descriptive and is not polarized SFS-B.
- Tier 3C composition is an exact-single-assembly observable and is not pooled with either diversity modality.
- The class-fixed GC3 slope remains `+0.00892005` per Buffalo proxy unit, with a 10,000-species-bootstrap 95% interval spanning zero (`-0.00123` to `+0.02149`; `n=90`; BH `q=0.217`). Class-specific signs remain heterogeneous and the quadratic coefficient has the opposite sign from predicted concavity.
- These recovered estimates do not establish that gBGC causally resolves Lewontin's paradox.

## Release validation from the integrated tree

### Pinned-Guix profile and full suite

```sh
bash analysis/validate_tier3_guix.sh

guix time-machine -C analysis/guix/channels.scm -- \
  shell -L analysis/guix -m analysis/guix/manifest.scm --pure -- \
  python3 -m pytest -q analysis/tests
```

Observed results:

- two independent profile realizations resolved to the identical store closure `/gnu/store/z9v2f6faha9cwjz0sm5iphhlzisgi077-profile`;
- reproducible-profile validator: **19 passed**;
- full analysis suite: **153 passed** in 2.93 seconds;
- the pure profile exposed the pinned scientific tools and did not expose an unapproved `impg` executable.

These match the passed synthesis-task log (19 profile tests and 153 full-suite tests).

### Focused provenance and numerical gates

```sh
guix time-machine -C analysis/guix/channels.scm -- \
  shell -L analysis/guix -m analysis/guix/manifest.scm --pure -- \
  python3 analysis/tier3_recovery_audit.py

guix time-machine -C analysis/guix/channels.scm -- \
  shell -L analysis/guix -m analysis/guix/manifest.scm --pure -- \
  python3 -m pytest -q \
    analysis/tests/test_tier3_recovery_audit.py \
    analysis/tests/test_tier3a_origin_remap.py \
    analysis/tests/test_tier3a_biological.py \
    analysis/tests/test_tier3b_recovery.py
```

Observed results:

- recovery audit: **23 recovered observations; 23 independent headline checks PASS**;
- focused suite: **22 passed**;
- artifact-index verification: **23/23 indexed files exist, have the recorded byte size, and match the recorded SHA-256**;
- rerunning the recovery audit reproduced all generated ledgers byte-for-byte and left the tracked tree clean.

### Deterministic manuscript build

Typst 0.13.1 was used. A fixed creation timestamp is mandatory because Typst otherwise intentionally embeds the current wall-clock time and a derived PDF identifier.

```sh
SOURCE_DATE_EPOCH=1784225429
typst compile --root . --creation-timestamp "$SOURCE_DATE_EPOCH" \
  manuscript.typ /tmp/manuscript-a.pdf
typst compile --root . --creation-timestamp "$SOURCE_DATE_EPOCH" \
  manuscript.typ /tmp/manuscript-b.pdf
cmp /tmp/manuscript-a.pdf /tmp/manuscript-b.pdf
```

The two fixed-timestamp builds were byte-identical. The complete rendered-object payload (bytes before PDF object 101, where document metadata begins) was also byte-identical to committed `manuscript.pdf`. The accepted PDF was originally built at the same UTC second without `--creation-timestamp`; consequently its XMP timestamp is serialized as `+00:00`, while the reproducible-build option serializes it as `Z`. This causes metadata length, xref offsets, and trailer bytes to differ without changing rendered content. The accepted, indexed PDF was deliberately preserved rather than rewritten during repository integration.

## Data and evidence locations

| purpose | committed location |
|---|---|
| corrected mapping manifest | `results/tier3a/acquisition_corrected_manifest.tsv` |
| corrected mapping commands/QC | `results/tier3a/acquisition_corrected_commands.tsv`, `results/tier3a/acquisition_corrected_qc.md` |
| mapping supersession negative control | `results/tier3a/acquisition_sweepga_supersession.tsv` |
| corrected Tier 3A estimates | `results/tier3a/diploid_diversity.tsv` |
| Tier 3A run provenance and fail-closed audit | `results/tier3a/diploid_run_manifest.tsv`, `results/tier3a/diploid_lineage_audit.json` |
| powered Tier 3B estimates | `results/tier3b/population_diversity.tsv` |
| Tier 3B run provenance and QC | `results/tier3b/population_run_manifest.tsv`, `results/tier3b/population_qc.md` |
| integrated result table | `analysis/tier3_results.tsv` |
| updated figures | `analysis/fig_tier3.png`, `analysis/fig_tier3.pdf` |
| updated manuscript | `manuscript.typ`, `manuscript.pdf` |
| recovery summary and full evidence ledger | `results/tier3/recovery_evidence_summary.md`, `results/tier3/recovery_evidence_ledger.tsv` |
| independent numerical audit | `results/tier3/headline_audit.tsv` |
| checksummed artifact index | `results/tier3/artifact_index.tsv` |
| rerun commands and scheduler evidence | `results/tier3a/diploid_rerun_commands.sh`, `results/tier3a/logs/`, `results/tier3b/population_rerun_commands.sh`, `results/tier3b/run_logs/` |

Large immutable acquisition inputs and run products remain at the absolute paths recorded in the acquisition and run manifests. The committed repository contains their checksums, commands, QC, scheduler telemetry, and small publication artifacts; it does not duplicate all large staged biological inputs.

## All-worktree audit and cleanup status

Every registered WG worktree was audited with `git status --porcelain`, branch/head inspection, commit-range inspection, and task-range versus main-integration patch-ID comparison. Non-ancestor counts below are expected because WG squashed task deltas into single-parent main commits; `patch-equal` proves that the branch delta is represented by the named main commit.

| agent | task | branch tip | commits not literal main ancestors | integration/content disposition | worktree status |
|---:|---|---|---:|---|---|
| 103 | remap-tier3a-origin-sweepga | `43bb47c15787` | 5 | patch-equal to `1f7f61253562`; exact tree | only `.wg-cleanup-pending` |
| 106 | rerun-tier3a-origin-sweepga | `dd23c1b9a51e` | 1 | patch-equal to `bc2e84684672`; exact tree | only `.wg-cleanup-pending` |
| 109 | synthesize-tier3-diversity-recovery | `ced8329a887c` | 1 | patch-equal to `043f460fb2ed`; exact tree | only `.wg-cleanup-pending` |
| 112 | integrate-tier3-recovery-main | starting `043f460fb2ed` | 0 | active final-integration worktree | clean before this handoff |
| 16 | validate-impg-bcf | `6b3de795566e` | 1 | unique commit object is an empty checkpoint commit; no file delta | only `.wg-cleanup-pending` |
| 20 | tier3-freeze-decisions | `1e2e0902e20e` | 1 | patch-equal to `5fc0feaea9fb` | only `.wg-cleanup-pending` |
| 23 | tier3-foundations-manifest | `413d45f6571b` | 1 | patch-equal to `4ea9b6527802` | only `.wg-cleanup-pending` |
| 26 | tier3-implement-vgp | `f65b44ff8e6a` | 1 | patch-equal to `3d37de2def21` | only `.wg-cleanup-pending` |
| 27 | tier3-implement-popvcf | `3cf07c49a566` | 1 | patch-equal to `13755ecf0488` | only `.wg-cleanup-pending` |
| 28 | tier3-implement-composition | `338974bd9f0f` | 1 | patch-equal to `6d1a05ce5365` | only `.wg-cleanup-pending` |
| 35 | tier3-integrate-cluster | `9a4a05ce5666` | 1 | patch-equal to `a52b34dab836` | only `.wg-cleanup-pending` |
| 38 | tier3-run-composition | `223758cd4134` | 23 | final task range patch-equal to `e0d079406e91`; retry history retained | only `.wg-cleanup-pending` |
| 39 | tier3-run-popvcf | `b0b731cc2562` | 1 | patch-equal to `c9b29f3eac6c` | only `.wg-cleanup-pending` |
| 40 | fix-pinned-wfmash | `9f7e3edbce38` | 1 | patch-equal to `810198fee229` | only `.wg-cleanup-pending` |
| 43 | tier3-run-vgp | `6397903515f3` | 1 | patch-equal to `5f914f879370` | only `.wg-cleanup-pending` |
| 51 | tier3-synthesize-results | `cf0809927705` | 1 | patch-equal to `333d4255781f` | only `.wg-cleanup-pending` |
| 59 | acquire-tier3a-vgp-tuples | `e2dd9434d903` | 1 | patch-equal to `03a4c5bed391` | only `.wg-cleanup-pending` |
| 60 | acquire-tier3b-population-tuples | `bdd1098b1933` | 2 | patch-equal to `a30b267c77d2` | only `.wg-cleanup-pending` |
| 80 | verify-sweepga-impg-handoff | `3dcf88544dca` | 1 | patch-equal to `40e6f4f0912e` | only `.wg-cleanup-pending` |
| 83 | run-tier3b-biological-recovery | `03d2259d0060` | 6 | patch-equal to `55712f4469b4`; exact tree | only `.wg-cleanup-pending` |
| 84 | repair-tier-3b | `205dae93d716` | 3 | patch-equal to `ba693e5e4b2f` | only `.wg-cleanup-pending` |
| 90 | superseded run-tier3a-biological-recovery | `7791b4528477` | 2 | intentionally retained only as rejected lineage; active selectors and audits exclude it | only `.wg-cleanup-pending` |
| 93 | build-sweepga-origin-main | `4cb38abf1ddd` | 1 | stable patch-equal to prerequisite `f9e11ee0f7fd`, included in `1f7f61253562` | only `.wg-cleanup-pending` |
| 99 | resolve-sweepga-origin-interface | `55712f4469b4` | 0 | points directly at integrated main history | clean |

No worktree contained an uncommitted scientific artifact. No worktree was removed: 22 completed worktrees carry WG's untracked zero-byte cleanup sentinel and the task explicitly forbids discarding any dirty worktree; agent 112 was active; agent 99 was clean but retained for WG to archive through its normal lifecycle. This conservative outcome preserves all branch objects and rejected-lineage evidence while leaving no unidentified user work.

## Safe release procedure

The final release uses only a fast-forward of the clean main branch followed by a normal push after re-fetching and verifying remote ancestry:

```sh
git fetch origin --prune
git merge-base --is-ancestor origin/main wg/agent-112/integrate-tier3-recovery-main
git -C /moosefs/erikg/lewontin-paradox merge --ff-only \
  wg/agent-112/integrate-tier3-recovery-main
git -C /moosefs/erikg/lewontin-paradox push origin main
git -C /moosefs/erikg/lewontin-paradox fetch origin
test "$(git -C /moosefs/erikg/lewontin-paradox rev-parse main)" = \
     "$(git -C /moosefs/erikg/lewontin-paradox rev-parse origin/main)"
```

The exact final pushed main hash and post-push equality check are recorded in the WG task log, avoiding a self-referential commit hash inside the commit that contains this document.
