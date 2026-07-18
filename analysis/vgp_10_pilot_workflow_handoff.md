# Ten-pair VGP pilot workflow handoff

Implementation freeze: 2026-07-18 UTC

Task: `implement-vgp-10-pilot`

Authorization state: implementation and deterministic software validation only;
no biological acquisition and no Slurm submission were performed.

## Delivered contract

The executable accounting and validation library is
`analysis/vgp_10_pilot.py`. The machine-readable stage contract is
`analysis/vgp_10_pilot_workflow.json`, and every final pair index must validate
against `analysis/vgp_10_pilot_output_schema.json`. The Slurm entry points are
under `analysis/slurm/vgp_10_pilot/`; their default submit mode is dry-run.

One unit always means one exact reciprocal, same-individual pair from
`analysis/vgp_10_pair_manifest.tsv`. H1 is the reference coordinate system and
H2 is the query. The implementation refuses a floating accession, a different
BioSample/individual, a reversed orientation, an unmeasured core QC field, or
an input digest/dictionary mismatch. Raw-read, k-mer, long-range switch, and
published-estimate checks are tracked as selective confidence evidence. Their
absence can distinguish confidence tier B from A, but is not silently promoted
to a universal core veto.

## Exact stage boundary

The only valid production order is:

```text
immutable whole H2 FASTA (query) + immutable whole H1 FASTA (target)
  -> SweepGA --num-mappings 1:1 --scaffold-jump 0
  -> independent H2-query and H1-target overlap-depth audit
  -> exact unchanged PAF
  -> IMPG index
  -> IMPG partition -> native partitions.bed
  -> exact native rows intersecting eligible query regions -> focus.native.bed
  -> IMPG query --sequence-files H1.fa H2.fa -b focus.native.bed -o vcf:poa
  -> IMPG lace -> H1 source-coordinate VCF
  -> bcftools norm -f H1 -m -any
  -> bcftools view -R verified eligible H1 1:1 BED
  -> bcftools norm -f H1 -d exact
  -> VCF.gz/TBI and BCF/CSI
  -> H1 REF audit and manifest-bound aligned-H2 reconstruction
  -> ordered reason-coded H1 mask and exact denominator
  -> diploid H1-coordinate consensus and masked PSMCFA
  -> unscaled PSMC + 200 frozen 5-Mb block bootstraps
  -> separate 1/10-Mb block sensitivities and separately sourced scaling scenarios
```

SweepGA receives both whole assemblies. It owns no VCF, variant, partition,
annotation, mask, or consensus operation. Both retained coordinate axes are
interval-swept independently after SweepGA; query or target depth greater than
one is fatal even if SweepGA exited zero. `--scaffold-jump 0` prevents SweepGA
scaffolding from trespassing on graph partition ownership.

The mapping join derives the H1 1:1 BED and its `not_1to1` complement directly
from that audited PAF. It also scans H1 for non-ACGT runs and projects H2
non-ACGT runs plus CIGAR deletion spans onto H1 on both `+` and `-` mappings.
CIGAR consumption must reconcile exactly to each PAF coordinate span. These
four core reason inputs cannot be replaced by an acquisition-supplied BED.

IMPG consumes the exact PAF used by the audit. The implemented command order is
`index`, `partition`, `query`, `lace`. The focus file is copied only from IMPG's
four-column native partition rows that intersect eligibility intervals. A
three-column user BED, a manufactured fixed-window scheme, or an empty native
selection fails. Both sequence files are mandatory at query time. IMPG emits
regional and laced VCF; it is not credited with BCF serialization.

bcftools owns decomposition, left normalization against exact H1, 1:1-region
trimming, exact deduplication, VCF/BCF serialization, and TBI/CSI creation.
Every REF substring is checked against H1. Normalized alternate alleles are
then applied in non-overlapping order to either whole manifest-matched contigs
or explicit strand-aware concordance regions. Every normalized record must be
assigned to exactly one such region, and the result must equal H2 sequence.

## Mask and consensus invariants

The H1 universe is reconciled using the following first-reason-wins order:

1. `not_eligible_contig`
2. `organellar`
3. `sex_linked_primary_exclusion`
4. `unplaced_or_unlocalized`
5. `h1_gap_or_N`
6. `h2_gap_or_N`
7. `not_1to1`
8. `mapping_breakpoint`
9. `low_base_accuracy`
10. `repeat_or_low_complexity_primary`
11. `duplication_or_collapse`
12. `phase_uncertain`
13. `other_predeclared`

Overlapping diagnostic evidence may be retained separately, but the accounting
complement is disjoint. Callable bp plus the 13 reason totals must equal the
declared H1 universe exactly; the output stores the integer denominator,
fraction, reason order, per-reason BEDs, and zero discrepancy. The primary
autosomal mask excludes sex-linked sequence; an explicitly labeled all-nuclear
sensitivity may use a different pre-registered flag set without changing the
primary object.

The consensus has H1 length and coordinates. Callable invariant bases carry
H1 sequence, validated heterozygous SNPs carry diploid IUPAC symbols, and every
non-callable base is `N`. Primary indels mask the H1 REF span plus 10 bp on both
sides; 0- and 50-bp sensitivities are frozen separately. A variant outside the
callable contract is fatal. The implementation reconciles non-`N` consensus
bases to callable bp minus the union of indel masks, so a masked base cannot be
quietly encoded as homozygous H1.

PSMCFA uses the native `fq2psmcfa` alphabet: `T` for a fully callable
homozygous 100-bp bin, `K` for a fully callable bin containing an IUPAC
heterozygote, and `N` if a bin is incomplete or contains any masked base. This
strict rule does not turn a partly callable bin into invariant sequence.

## Bootstrap and scenario freeze

The primary plan contains exactly 200 attempts. Five-megabase units are cut
within each continuous callable interval, never across a contig or mask
discontinuity. The unit BED, deterministic manifest-derived seed, and sampled
unit indices are durable. A replicate PSMCFA is generated atomically in
job-local scratch; each sampled unit is a separate FASTA record so neither a
100-bp encoding bin nor PSMC sequence crosses its boundary. One- and ten-
megabase unit BEDs are frozen as sensitivities. Attempts below 100 are refused.

The primary PSMC object is unscaled. Mutation rate, generation time, both
sources, and scenario ID occur only in a separate scaled table. Scenario
conversion never mutates or replaces the unscaled interval/time/lambda rows.

## Annotation branch

Annotation is optional for genome-wide diversity and PSMC. `annotation_stage.sh`
always writes an audit. A missing annotation writes a non-blocking
`not_available` disposition. An available native annotation must name the
exact H1 assembly accession and version and contain a sequence dictionary
identical in name, order, length, and sequence MD5. A different source assembly
is refused unless its liftover record is manifest-bound, targets the exact H1,
contains fixed chain/validation/manifest SHA-256 values, and passed validation.

Only after that gate may callable/core variants intersect `CDS`, `fourfold`,
`nonsynonymous`, `synonymous`, `WS`, `SW`, and `GC3` feature BEDs. Every feature
summary retains its own callable H1 denominator; zero denominators are null,
not zero diversity. WS means A/T REF to G/C ALT and SW the reverse for canonical
SNVs. Annotation failure cannot rewrite a valid core result.

## GNU Guix identity

All runtime packages are selected through
`analysis/guix/vgp_10_pilot/channels.scm` at Guix commit
`44bbfc24e4bcc48d0e3343cd3d83452721af8c36`. The production manifest SHA-256
is `afd14da446d29d4e66b64837194dfae524342dc5b8de3c2778ec85c1cca09a3d`.
`environment-lock.json` binds that manifest, the SweepGA build manifest, the
IMPG handoff manifest/package, source commits, submodules, companion hashes,
and the PSMC recipe.

PSMC 0.6.5 is pinned to commit
`b37b1cfa05b89c67c2ad1b63c699a27600d5516e`; its archive SHA-256 is
`4f6dff8a4d5a4cf98928100c87487d2fec5657da186da0f36795f46e70c48040`.
Guix built `psmc`, `fq2psmcfa`, and `splitfa` and repeated the identical
derivation with `guix build --check`. The captured profile is
`/gnu/store/3c2mxm30rbzvnw7qsi235mrkk3m38fym-profile`, its derivation is
`/gnu/store/g05w4h08h891dcaf7v4gy3irf2zwbwa7-profile.drv`, and the 151-item
transitive closure-list SHA-256 is
`8fcdb32021f1cd8eac839509cff47ab6bdd63b656b30e243fdf78d3c4ba24f9d`.

`analysis/guix/vgp_10_pilot/realization.json` records every executable path and
digest. Critical accepted hashes are:

| executable | SHA-256 |
| --- | --- |
| SweepGA | `fa7f0edb9b7e275c288db254046020e136d4267dd5ee043379227ef80da0573b` |
| IMPG | `c587dc2326cd24f887b1fcb3938404229ad0f0a27ef0773e90c287b1ade160d4` |
| gfaffix | `4bc1c5e236a8fe6aa1dbcff6e6cf515e8a70c808549990a515c3c5212776a627` |
| PSMC | `2c64b0ce7f68a6251e4795f3b88caa7b53769e12a867b6f59918312f66735f57` |
| fq2psmcfa | `9fc53b98c9111c390d2f652df6fbdeb7405a5ddc5f9272219c6ddcec566f6442` |
| splitfa | `8e015b166824ec660ed38bcff784c654fafe7f5db6c2fece37f8120eeb9b329a` |

The capture also binds `ALNtoPAF`, `FAtoGDB`, `FastGA`, `GIXmake`, `GIXrm`,
`ONEview`, `PAFtoALN`, and `wfmash`; the Slurm job PATH is constructed only
from captured executable directories plus the captured Guix profile. Ambient
modules, Conda, virtualenvs, containers, rustup, and network installers are not
accepted.

## Slurm operation

Required exported values are `VGP_DATA_ROOT`, `VGP_RUN_ID`, and
`VGP_ENVIRONMENT_CAPTURE`. Outside Slurm, `VGP_SCRATCH_ROOT` is also required;
on nodes, `$SLURM_TMPDIR` is used. Each pair needs:

- `pilot/inputs/Pxx/input-manifest.json` with immutable pair/QC/evidence data;
- whole `h1.fa` and `h2.fa` plus their exact dictionary/hash records;
- the H1 universe, eligible query regions, and manifest-bound reason BEDs for
  eligibility, organelles, sex linkage, placement, breakpoints, base accuracy,
  repeats/complexity, collapse, and phase (1:1 and both N/gap masks are derived
  from the audited PAF/FASTA at mapping time);
- explicit concordance regions if whole contig reconstruction is not valid;
- `resources.json` with per-job CPUs, Slurm time, memory, and approved scratch
  high estimate; and
- optional, strictly bound annotation and feature BED records.

Preview only:

```bash
export VGP_DATA_ROOT=/authorized/absolute/root
export VGP_RUN_ID=vgp10-run-id
export VGP_ENVIRONMENT_CAPTURE="$PWD/analysis/guix/vgp_10_pilot/realization.json"
bash analysis/slurm/vgp_10_pilot/submit.sh --dry-run
```

`--submit` is reserved for the downstream authorized run task. The submitter
reads pair-specific resource limits; scripts contain no arbitrary global byte,
memory, or wall ceiling. Preflight, mapping, IMPG, variants, and consensus are
dependency-linked. PSMC primary task zero and bootstrap tasks 1–200 start only
after the consensus sentinel. A pair-level finalize job requires all 201 array
tasks, preserves the unscaled table, enforces at least 190/200 finite bootstrap
trajectories, and writes scenario-scaled trajectories separately. Annotation
is a separate optional branch.

Every stage writes to a job-specific `.partial` directory, emits failure or
success telemetry, hashes all files into `.complete.json`, and promotes with a
single rename. A final directory without a sentinel is fatal. Rerunning a
complete stage is a no-op; stale partials are retained with a timestamp rather
than confused with success. Telemetry records job/array ID, elapsed epoch,
maximum child RSS, scratch path, retry count, and disposition. The run task must
augment this with `sacct` CPU, wall, read/write, scratch high-water, metadata
operations, and estimator error before program adjudication.

## Validation and refusal coverage

`analysis/tests/test_vgp_10_pilot.py` uses only tiny deterministic fixtures for
software truth. It covers whole-input/native 1:1 commands; query and target
multiplicity independently; exact IMPG order and native partition selection;
normalization/dedup/index commands; explicit tool-role refusals; ordered mask
accounting; denominators; REF/H2 SNP and indel reconstruction; IUPAC/indel/N
consensus behavior; strict PSMCFA; 200 boundary-aware bootstraps; unscaled versus
scenario output; annotation accession/dictionary/liftover gates; confidence;
resource estimation; atomic resume; Guix locks; Slurm syntax/order; and output
schema validity.

The full pinned analysis suite is run locally inside the pre-existing frozen
analysis profile. This implementation validation submits no Slurm job and uses
no biological result as a fixture or conclusion.

## Downstream run checklist

1. Copy the frozen design manifests into the authorized data-root run manifest
   and verify their hashes before inspecting results.
2. Complete all ten input manifests or record explicit pre-result failures; do
   not activate an alternate without the versioned same-clade amendment rule.
3. Re-run `capture_environment.sh` on the execution host and compare every
   identity to this capture. A changed store path is acceptable only if all
   locked inputs, derivations, closure contents, and executable bytes reconcile.
4. Generate measured per-pair resource records and confirm scheduler quota,
   storage headroom, scratch, and I/O authorization.
5. Dry-run `submit.sh`, review exact dependencies/limits, then use `--submit`
   only under the `run-vgp-10-pilot` authorization.
6. Stop a pair on any hard gate. Never retry a scientific failure with relaxed
   mapping, mask, concordance, or annotation rules.
7. Preserve unscaled trajectories and bootstrap failures. Add scaling scenarios
   only through a separately versioned source table.
8. Validate every pair output index against the schema and ensure all ten slots
   end as pass or explicit failure before evaluating the 8/10 program gate.

Nothing in this implementation is a biological result. Passing tiny fixtures
proves software invariants only.
