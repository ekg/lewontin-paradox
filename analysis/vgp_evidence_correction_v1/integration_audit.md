# Reviewed repair-lineage integration audit (correct-vgp-evidence)

Audit date: 2026-09-03/04 UTC · task: `correct-vgp-evidence`

## Question

Is the independently reviewed repair lineage fully integrated onto `main` in
dependency order, with no task-owned content missing, and is `main` equal to
`origin/main`?

## Verdict

**Yes.** Every task lineage that produced or reviewed the corrected evidence is
present on `main` in dependency-correct order, content-complete, and
`main == origin/main == 599100e57ef1f7f42e291d3c656ca429db32dc57` at audit
time (verified after `git fetch origin main`). No task commit needed replaying;
replaying would have duplicated already integrated content. This audit, the
correction publication, and its tests are the only new artifacts.

## Integration order on main (verified)

```text
ee451c9  feat: repair-vgp-psmc (agent-368)              2026-07-20 18:35:38 +0000
  -> 2044765  feat: scale-vgp-real (agent-371)          2026-07-20 20:58:38 +0000
  -> d3ca516  feat: validate-vgp-pilot-reads (agent-365) 2026-07-21 18:12:26 +0000
  -> 0c20b12  feat: synthesize-vgp-real (agent-376)     2026-07-21 18:33:38 +0000
  -> 04659d1  feat: integrate-vgp-repair-base (agent-385) 2026-07-22 17:33:38 +0000
  -> ed1c82f  feat: catalog-vgp-gff (agent-384)         2026-07-22 18:06:01 +0000
  -> 65f422e  feat: run-vgp-clean-canary (agent-390)    2026-07-22 22:54:57 +0000
  -> bee7757  feat: run-vgp-three-pair (agent-394)      2026-07-25 15:15:52 +0000
  -> 59cf8e4  docs: land independent review of bounded VGP three-pair pilot
              (review-vgp-three-pair, agent-397)        2026-09-03 21:19:13 +0000
```

This order is dependency-correct: the PSMC repair precedes the scale-out whose
results it repairs; read validation precedes the synthesis that classifies it;
the repair-base integration and annotation catalog precede the clean canary
that consumes both; the clean canary precedes the three-pair run that reuses
its closed evidence; the independent review is the join that gates this
correction. The review landed as a true merge (`59cf8e4` parents
`bee7757 fdc0e37`), so both the bounded run and the review content are first
parents of the current evidence state.

## Content-completeness verification

Integration used WG's squash-landing convention: each task branch's full
content landed as one `feat:/docs: ... (task, agent-N)` commit on `main`.
Completeness was verified per branch, in increasing strictness:

1. **Tree-identity against the landing commit** (`git diff <landing> <branch>
   --stat` empty ⇒ the landing commit's tree contains the branch tip exactly):

   | task branch | landing commit | result |
   |---|---|---|
   | `wg/agent-368/repair-vgp-psmc` | `ee451c9` | tree-identical |
   | `wg/agent-371/scale-vgp-real` | `2044765` | tree-identical |
   | `wg/agent-376/synthesize-vgp-real` | `0c20b12` | tree-identical |
   | `wg/agent-390/run-vgp-clean-canary` | `65f422e` | tree-identical |
   | `wg/agent-393/run-vgp-three-pair` | `bee7757` | tree-identical |
   | `wg/agent-397/review-vgp-three-pair` | `59cf8e4` | branch tip is the merge itself |

2. **Per-file blob containment** for the two branches whose landing commits
   interleaved with later main content
   (`wg/agent-384/catalog-vgp-gff` -> `ed1c82f`,
   `wg/agent-365/validate-vgp-pilot-reads` -> `d3ca516`): every file version in
   each branch tip appears somewhere in `main`'s history of that file (verified
   by blob-hash search over `git log main -- <file>`); no branch-only file
   content exists.

Replay was therefore refused as duplicative, following the precedent of
`analysis/repaired_vgp_integration_handoff.md` ("Replaying a task branch would
therefore duplicate already integrated evidence and tooling").

## Branch hygiene

- No reset, rebase, amend, stash, force push, clean, or broad staging was used.
- Older non-merged worktree branches unrelated to this lineage (tier3-era
  branches) were left untouched; their content was superseded by later
  recovery tasks and is out of scope here.
- The shared root worktree (`/moosefs/erikg/lewontin-paradox`) was on
  `main...origin/main`, clean, and was not modified by this task.

## Post-landing verification contract

After this correction lands on `main` and is pushed:

```sh
git fetch origin main
git rev-parse main origin/main   # must be equal
```

The result is recorded in the `correct-vgp-evidence` WG task log. No VGP
scale-out was launched by this integration; the review's broad-scale NO-GO
stands until its GO conditions are met.
