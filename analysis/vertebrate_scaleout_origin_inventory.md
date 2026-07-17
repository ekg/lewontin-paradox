# Vertebrate scale-out origin inventory

Date: 2026-07-17 UTC

WG task: `rerun-sync-origin`

Scope: read-only origin reconciliation audit plus planning inventory; no data
download, Slurm submission, or scale-out computation was performed.

## Executive disposition

The live project root was already synchronized with `origin/main` before this
task fetched. Both refs were
`14fe57d1d0593a2859672279be60aa858e822410`, the upstream VGP guidance commit,
with ahead/behind `0 0`. The previously reported one-commit-per-side
divergence therefore no longer existed. In particular, the known local
quality-pass commit
`70c17e73e0874151cf13ed935c7daef6e96d46e1` was not the live `main` tip and
was not an ancestor of live `main`.

The task's merge authorization was explicitly conditional on the live state
still having the previously understood divergence. That condition failed, so
no merge was attempted. Creating a new merge from the active quality-pass task
branch would have been a materially different operation. The fetch changed
only Git's fetch metadata; no tracked file, index entry, local commit, branch
tip, or worktree content in the live root was changed by the reconciliation
audit.

The final live relation remains `main == origin/main == 14fe57d1...` with
ahead/behind `0 0`. The upstream guidance commit is an ancestor of `main`; the
quality-pass commit is not. Both commit objects remain intact, and the latter
is protected by local and remote task-branch refs
`wg/agent-118/quality-pass-vertebrate-scaleout-plan` and
`origin/wg/agent-118/quality-pass-vertebrate-scaleout-plan`.

## Git safety and provenance record

### Pre-fetch state

All values below were read in the live root
`/moosefs/erikg/lewontin-paradox` before `git fetch origin main`.

| field | observed value |
|---|---|
| live branch | `main`, tracking `origin/main` |
| pre-fetch `HEAD` | `14fe57d1d0593a2859672279be60aa858e822410` |
| pre-fetch `main` | `14fe57d1d0593a2859672279be60aa858e822410` |
| pre-fetch `origin/main` | `14fe57d1d0593a2859672279be60aa858e822410` |
| origin fetch URL | `git@github.com:ekg/lewontin-paradox.git` |
| origin push URL | `git@github.com:ekg/lewontin-paradox.git` |
| merge base of live refs | `14fe57d1d0593a2859672279be60aa858e822410` |
| `main...origin/main` ahead/behind | `0 0` |
| tracked worktree delta | none |
| staged/index delta | none |
| porcelain status exception | only `? .wg-worktrees/`, the directory containing WG-managed worktrees; it was not staged, modified, or removed |

`git cat-file -t` returned `commit` for both known SHAs. The commits have the
same parent, the previously reported common base
`dec7266d3427987d04f65bb1300032c5f233cb95`:

| side from the historical divergence | commit | parent | commit/path inventory |
|---|---|---|---|
| upstream guidance | `14fe57d1d0593a2859672279be60aa858e822410` | `dec7266d3427987d04f65bb1300032c5f233cb95` | one commit; adds `results/tier3/vgp_freeze_analysis.md` |
| local quality pass | `70c17e73e0874151cf13ed935c7daef6e96d46e1` | `dec7266d3427987d04f65bb1300032c5f233cb95` | one commit; adds `analysis/vertebrate_scaleout_quality_pass.md` |

The path sets are disjoint. A read-only historical
`git merge-tree dec7266d... 70c17e7... 14fe57d...` reported only the guidance
file as “added in remote” and emitted no conflict record or conflict marker.
That proves the two historical trees are mechanically mergeable, but it does
not override the failed live-state condition or authorize importing an active
task branch into `main`.

### Fetch and post-fetch proof

The only origin operation was:

```text
git fetch origin main
```

It completed successfully with `main -> FETCH_HEAD` and did not advance
`origin/main`. Post-fetch observations were:

| field | observed value |
|---|---|
| post-fetch live `HEAD` | `14fe57d1d0593a2859672279be60aa858e822410` |
| post-fetch `main` | `14fe57d1d0593a2859672279be60aa858e822410` |
| post-fetch `origin/main` | `14fe57d1d0593a2859672279be60aa858e822410` |
| post-fetch merge base | `14fe57d1d0593a2859672279be60aa858e822410` |
| post-fetch ahead/behind | `0 0` |
| `14fe57d...` ancestor of `main` | yes; `git merge-base --is-ancestor` exit 0 |
| `70c17e7...` ancestor of `main` | no; `git merge-base --is-ancestor` exit 1 |
| `origin/main` ancestor of `main` | yes; exit 0 |
| `main` ancestor of `origin/main` | yes; exit 0 |

No `merge`, `pull`, `rebase`, `reset`, `clean`, `checkout`, `restore`,
`commit --amend`, force update, or push to `main` was run. The precise
non-destructive blocker to the requested merge is: live `main` and
`origin/main` are already equal at the upstream commit, so they do not present
the authorized one-commit-per-side divergence; the other known commit is an
active task-branch commit, not a live-main-side commit.

### Exact audit commands

Commands were run from the live root unless noted otherwise:

```text
git status --porcelain=v2 --branch
git diff --name-status
git diff --cached --name-status
git remote get-url origin
git remote get-url --push origin
git rev-parse HEAD
git rev-parse refs/heads/main
git rev-parse refs/remotes/origin/main
git merge-base refs/heads/main refs/remotes/origin/main
git rev-list --left-right --count refs/heads/main...refs/remotes/origin/main
git show --no-patch --format=fuller 14fe57d1d0593a2859672279be60aa858e822410
git show --no-patch --format=fuller 70c17e73e0874151cf13ed935c7daef6e96d46e1
git cat-file -t 14fe57d1d0593a2859672279be60aa858e822410
git cat-file -t 70c17e73e0874151cf13ed935c7daef6e96d46e1
git merge-base --is-ancestor 14fe57d1d0593a2859672279be60aa858e822410 refs/heads/main
git merge-base --is-ancestor 70c17e73e0874151cf13ed935c7daef6e96d46e1 refs/heads/main
git merge-base --is-ancestor 70c17e73e0874151cf13ed935c7daef6e96d46e1 refs/remotes/origin/main
git log --left-right --cherry-mark --oneline refs/heads/main...refs/remotes/origin/main
git diff-tree --no-commit-id --name-status -r 14fe57d1d0593a2859672279be60aa858e822410
git diff-tree --no-commit-id --name-status -r 70c17e73e0874151cf13ed935c7daef6e96d46e1
git fetch origin main
git rev-parse 14fe57d1d0593a2859672279be60aa858e822410^
git rev-parse 70c17e73e0874151cf13ed935c7daef6e96d46e1^
git diff --name-only dec7266d3427987d04f65bb1300032c5f233cb95..14fe57d1d0593a2859672279be60aa858e822410
git diff --name-only dec7266d3427987d04f65bb1300032c5f233cb95..70c17e73e0874151cf13ed935c7daef6e96d46e1
git merge-tree dec7266d3427987d04f65bb1300032c5f233cb95 70c17e73e0874151cf13ed935c7daef6e96d46e1 14fe57d1d0593a2859672279be60aa858e822410
git merge-base --is-ancestor refs/remotes/origin/main refs/heads/main
git merge-base --is-ancestor refs/heads/main refs/remotes/origin/main
```

## Upstream guidance inventory

The controlling upstream artifact is
`results/tier3/vgp_freeze_analysis.md` (98 lines, 4,481 bytes, SHA-256
`0ad6fa03ceeec9d07c39c5456ddc4702c54c586115e578af3599d07f29e5316d`).
It was introduced by `14fe57d...` and is tracked at the final live and task
branch trees.

### Reported candidate counts

These are reported counts, not locally reproducible counts, because the raw
TSV is absent.

| reported population/filter | count | interpretation |
|---|---:|---|
| unique species in freeze | 714 | species-level denominator |
| completed assemblies (`status=4`) | 223 | assembly-complete subset |
| annotation status `Completed NCBI` or `Ready` | 248 | annotation-status subset, not necessarily assembly-complete |
| non-empty paired-haplotype field | 271 | pair-label subset, not necessarily completed/annotated |
| completed assembly + annotation | 120 | all taxa |
| completed + annotation + paired | 40 | guidance's Tier 3A-ready headline |
| triple-eligible fish | 13 | Actinopterygii/Sarcopterygii as grouped by the guidance |
| triple-eligible amphibians/reptiles/mammals/birds/other | 4 / 3 / 9 / 6 / 5 | sums with fish to 40 |
| completed + RefSeq-annotated fish | 46 | Tier 3C composition candidates regardless of pair status |
| RefSeq-annotated fish with a pair | 13 | strongest shared Tier 3A/Tier 3C fish subset |

### Catalog eligibility versus execution eligibility

The guidance's catalog rule uses TSV columns 10, 13, 16, 17, 21, and 26:
scientific name; completed status `4`; H1 GCA; RefSeq GCF; non-empty paired
assembly ID; and annotation status. A row belongs to the reported 40 only when
completed, annotated, and paired conditions intersect. A row belongs to the
reported 46-fish Tier 3C pool when it has a completed assembly and RefSeq
annotation; a pair is not required for composition.

That catalog screen is necessary but insufficient for execution. Before a row
can enter a frozen Tier 3 workflow, planning must resolve and validate:

1. exact versioned H1 and H2 assembly accessions, reciprocal relationship,
   same biological individual/isolate, haplotype roles, release metadata, and
   FASTA provider hashes plus repository SHA-256;
2. for direct diploid assembly diversity, an H1-reference alignment/callability
   mask and an auditable denominator, rather than interpreting a pair label as
   callable diploid genotype evidence;
3. a native annotation on the identical H1 accession/version for
   annotation-derived results, with exact GFF/FASTA hashes, provider/release,
   genetic code, sequence-region/contig mapping, and reconstructed CDS audit;
4. explicit native-versus-projected status—projected or accession-mismatched
   annotation cannot be a primary GC3/4D input;
5. checksummed commands, pinned toolchain/profile, resource telemetry,
   failure/missingness status, and atomic output/QC records.

This is consistent with the local validators. `analysis/tier3_manifest.py`
requires the annotation accession to equal the reference accession and native
status for annotation-derived or diversity-eligible results. The direct
WFMASH mode requires an H1-reference alignable denominator.
`analysis/tier3a_vgp_compute.py` fails closed on incomplete provenance,
non-native annotation, reference-accession mismatch, checksum mismatch,
contig-dictionary mismatch, absent callability, and zero callable bases.

### The 13 fish rows and local resolution status

| species | guidance H1 GCA | guidance H2 label | locally synchronized row assets |
|---|---|---|---|
| *Acipenser ruthenus* | `GCA_902713425.2` | `fAciRut3.pat` | none outside guidance |
| *Amia calva* | `GCA_036373705.1` | `fAmiCal2.hap2` | none outside guidance |
| *Lepisosteus oculatus* | `GCA_040954835.1` | `fLepOcu1.hap1` | no scale-out tuple; species occurs separately in the existing Tier 3C manifest |
| *Lampris incognitus* | `GCA_029633865.1` | `fLamInc1.hap1` | none outside guidance |
| *Syngnathus acus* | `GCA_901709675.2` | `fSynAcu2.pri` | none outside guidance |
| *Enoplosus armatus* | `GCA_043641665.1` | `fEnoArm2.hap2` | none outside guidance |
| *Pempheris klunzingeri* | `GCA_042242105.1` | `fPemKlu1.hap2` | none outside guidance |
| *Cyclopterus lumpus* | `GCA_009769545.1` | `fCycLum2.pri` | none outside guidance |
| *Heterodontus francisci* | `GCA_036365525.1` | `sHetFra1.hap1` | none outside guidance |
| *Hemiscyllium ocellatum* | `GCA_020745735.1` | `sHemOce1.mat` | none outside guidance |
| *Heptranchias perlo* | `GCA_035084215.1` | `sHepPer1.hap2` | none outside guidance |
| *Pristiophorus japonicus* | `GCA_044704955.1` | `sPriJap1.hap2` | none outside guidance |
| *Narcine bancroftii* | `GCA_036971445.1` | `sNarBan1.hap1` | none outside guidance |

Before this inventory was added, repository searches found each of those 13
H1 accessions and H2 labels only in the guidance file. The table omits the
candidates' H2 GCA accession
versions and RefSeq GCF annotation accession versions, so it cannot itself be
used as an acquisition or execution manifest.

## Referenced assets and caveats

| asset/reference | local verification | consequence |
|---|---|---|
| `results/tier3/VGPPhase1-freeze-1.0.tsv` | **missing** and absent from `git ls-files` | the reported 717-line/716-row source and all aggregate counts cannot be recomputed locally |
| VGP GitHub blob/raw URLs | recorded in guidance; not fetched in this task | the moving `main` URL is not a frozen, checksummed input |
| source date | guidance says `2025-01-XX` | placeholder is not adequate retrieval provenance |
| `results/tier3/vgp_freeze_analysis.md` | present, tracked, hashed above | usable as guidance and reported-count provenance only |
| PRJNA489243 | named as Tier 3A acquisition project | a project identifier does not replace row-level versioned accessions, pair/individual proof, or checksums |
| 13 fish H1 GCAs | present only as text in guidance | H1 binaries were not staged or downloaded |
| 13 H2 assembly labels | present only as text in guidance | labels must be resolved to versioned accessions and reciprocal/same-individual evidence |
| RefSeq GCF annotations for the 13 | asserted but GCF values absent from the table | exact reference-coupled FASTA/GFF selection remains a preflight requirement |
| 46 Tier 3C fish | only an aggregate count; no row list | cannot produce a frozen 46-row composition manifest from the committed guidance |
| current three Tier 3A prototypes | present as three rows in `results/tier3a/acquisition_corrected_manifest.tsv` | useful calibration only; they do not supply the 13-row candidate tuples |
| current Tier 3C production manifest | 135 rows in `analysis/tier3c_manifest.tsv`; includes *L. oculatus* | existing production evidence is not a synchronized VGP-freeze candidate inventory |

There is a direct contradiction in the upstream text: line 5 says the raw
manifest “is committed alongside this document,” while the file is absent and
the introducing commit message says it is not committed. The provided command
also fetches the “latest” branch copy rather than a content-addressed freeze.
Downstream planning must not silently execute that command. It should first
define an approved acquisition step that records immutable source revision,
retrieval timestamp, byte size, provider checksum, repository SHA-256, row
count, schema/header hash, and derived-query code/hash. This task deliberately
did not repair the upstream guidance, download the TSV, or infer missing rows.

The three current prototypes need careful interpretation. The local corrected
manifest includes *Spinachia spinachia*, *Menidia menidia*, and
*Tautogolabrus adspersus* with extensive exact-input provenance. The guidance
classifies the first two as paired but not RefSeq-annotated and the third as
having only an alternate assembly rather than a reciprocal pair. Therefore
prototype execution success must not be mistaken for satisfaction of the
freeze's catalog-level “paired + RefSeq annotation” rule.

## GNU Guix implications and smoke gate

The production baseline is pinned by:

- `analysis/guix/channels.scm`, Guix commit
  `44bbfc24e4bcc48d0e3343cd3d83452721af8c36`, SHA-256
  `45c055cd1d9010a72eacbb720037a22bccb2d8d6891dbd11b5d66365f29b3a17`;
- `analysis/guix/manifest.scm`, SHA-256
  `2fb05e87aa2ac45ce51d4dcf93b232cb98627f525adace98357629ee3f15720a`;
- the already-realized recorded profile
  `/gnu/store/z9v2f6faha9cwjz0sm5iphhlzisgi077-profile`, verified present;
- `analysis/pilot_results/guix_environment.json`, which records the same
  channel/manifest hashes and profile store path.

Scale-out planning must reuse the validated package lineage and record exact
profile/store paths for each stage. Ambient Python, alignment, annotation, or
VCF tools are not substitutes. The current general manifest deliberately has
no executable IMPG; if assembly scale-out requires the validated SweepGA/IMPG
path, planners must select the already-audited supplemental manifest and exact
commits rather than silently adding tools to the baseline.

A proportionate no-download smoke used `analysis/slurm/guix_job.sh`, which
first verified the recorded profile and committed channel/manifest hashes,
then executed only tool/import checks and the two foundational test modules.
The exact successful command was:

```sh
analysis/slurm/guix_job.sh "$PWD/analysis/pilot_results/guix_environment.json" bash -c '
set -euo pipefail
python3 --version
pytest --version
wfmash --version
bcftools --version
samtools --version
python3 - <<"PY"
from importlib.metadata import version
import Bio, numpy, pandas, pyfaidx, pysam, scipy
print("scientific_versions", Bio.__version__, version("jsonschema"), numpy.__version__, pandas.__version__, pyfaidx.__version__, pysam.__version__, scipy.__version__)
PY
python3 -m pytest -q analysis/tests/test_common.py analysis/tests/test_manifest.py
if command -v impg >/dev/null; then
  echo "ERROR: unapproved impg executable present" >&2
  exit 1
fi
printf "impg=absent_as_required\n"
'
```

Observed: exit 0; Python 3.10.7; pytest 7.1.3; bcftools 1.14; samtools 1.14;
Biopython 1.80, jsonschema 4.5.1, NumPy 1.23.2, pandas 1.4.4, pyfaidx
0.7.2.1, pysam 0.20.0, SciPy 1.10.1; `19 passed in 0.36s`; unapproved
`impg` absent as required. `wfmash --version` exited successfully but emitted
no version text. An initial wrapper attempt exited 127 before tests because it
called `sed`, which is intentionally absent from the pure manifest; the
corrected command removed that ambient-tool assumption. No Guix realization,
package download, compute-node smoke, or biological analysis was launched.

## Cluster phases for downstream planning

The committed cluster framework provides patterns, not authorization to run.
The downstream execution plan should preserve these ordered gates:

1. **Guidance/schema freeze.** Acquire the raw catalog only after explicit
   approval; content-lock it, validate the 717-line claim and columns, emit a
   row-level eligibility/rejection ledger, and resolve all accessions and
   reference relationships. Stop if counts cannot be reproduced.
2. **Asset preflight and resource model.** Resolve provider metadata without
   staging bulk sequence, estimate low/base/high bytes, inodes, MooseFS I/O,
   memory, walltime, core-hours, and temporary amplification per modality,
   establish quotas/retention, and define download and execution approvals.
3. **Pinned environment preparation.** On `octopus01`, realize and GC-root the
   exact time-machine profile; record channel/manifest hashes, derivations,
   closure, tool paths, versions, and shared-store identity. Compute nodes use
   `guix_job.sh` and must not build against an unavailable daemon.
4. **Approved staging.** Use bounded arrays to stage exact H1/H2 FASTA and
   native H1 GFF assets, verify provider checksums and SHA-256, create contig
   maps and provenance, and write explicit failures rather than shrinking the
   cohort. For composition, follow the existing `discover -> stage-one ->
   freeze` contract before `run-one`.
5. **Smoke and stratified pilot.** Run a small truth-case/pilot across genome
   sizes and clades. Require exact-reference, phase/pair, native-annotation,
   callability, checksum, storage, runtime, and scientific QC to pass. Review
   measured resource bounds before any expansion.
6. **Reviewed expansion waves.** Submit bounded, throttled arrays only after an
   explicit go decision, with `afterok` dependencies, local scratch, immutable
   inputs, atomic promotion, idempotent fingerprints, retry ledgers, and
   `sacct` telemetry. Separate Tier 3A paired-assembly work from Tier 3C
   annotation/composition work so ineligible rows do not cross modalities.
7. **Collection and independent QC.** Collect every eligible, failed, and
   unavailable row; retain checksums and telemetry; run independent controls;
   compare observed to budget; and require another review before a full
   catalog or demographic extension.

The existing generic submission path performs Guix preparation, workflow
preflight, a 500-GiB free-space gate, compute smoke, and dependency-gated Tier
3 arrays. The Tier 3C production path further separates discovery, staging,
freeze, bounded standard/outlier run lanes, collection, and independent
control audit. These controls should be adapted only through an approved,
versioned plan; none was invoked here.

## Downstream planning handoff

For `plan-all-vertebrate-tier3`:

- treat 714/223/248/271/120/40 and the 46/13 fish counts as reported guidance
  until the missing raw TSV is content-locked and independently recomputed;
- produce one row per catalog species with source revision, retrieval date,
  accessions/versions, pair and individual proof, annotation coupling,
  checksum/provenance fields, modality eligibility, uncertainty, and rejection
  reason;
- distinguish the 40 all-taxon triple candidates, the 13 paired fish, and the
  46 composition-capable fish rather than using one “ready” flag;
- calibrate costs from the three current Tier 3A tuples and existing Tier 3C
  telemetry, but do not generalize without genome-size/clade-stratified bounds;
- preserve the approval boundary before catalog download, staging, pilot,
  expansion waves, and full execution.

For `plan-vertebrate-ne-strategy`:

- do not treat deposited H1/H2 assemblies as a callable diploid genotype or a
  population sample by default;
- separately inventory raw reads, same-individual relationship, heterozygous
  genotype generation, callability masks, phase evidence, population sample
  count, and mutation/generation-time calibration;
- keep PSMC, MSMC2, and SMC++ eligibility method-specific and keep independent
  ecological/Ne covariates separate from estimates derived from the response
  heterozygosity itself;
- make missing demographic inputs explicit rather than lowering assembly or
  annotation gates.

For synthesis and execution governance, planning artifacts should crosswalk
every guidance claim to a frozen row/field, local validation rule, resource
estimate, stage owner, go/no-go threshold, and authorization task. Completion
of a plan is not authorization to download data or submit jobs.
