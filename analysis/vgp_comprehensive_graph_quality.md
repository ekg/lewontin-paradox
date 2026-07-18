# Comprehensive VGP program graph quality review

Date reviewed: 2026-07-18 UTC

WG task: `quality-vgp-psmc`

Scope: eleven task definitions from design through synthesis

Execution boundary: metadata and repository review only; no biological data
were downloaded and no Slurm job was submitted

## Verdict

The reviewed graph is fit to begin design. All eleven task definitions,
scientific dependencies, timeouts, and deliverables were rewritten or
reasserted. The graph now treats exact, same-individual, mutually comparable,
high-quality H1/H2 assemblies as a potentially valid source of heterozygosity:
after whole-haplotype 1:1 alignment, variant extraction, and a reason-coded
callable mask, they can support callable-genome diversity and assembly-derived
PSMC. This is a positive input contract, not an automatic eligibility claim
for every linked pair.

The earlier negative finding remains valid only for the six rows and evidence
state it actually audited. Those rows lacked a validated comparable H2,
callable consensus, and mask at audit time; that does not establish that a
future exact H1/H2 program is conceptually invalid. The integrated handoff
already makes the key distinction under “Annotation-pilot refusal is not PSMC
infeasibility” at `analysis/repaired_vgp_integration_handoff.md:282` and records
that the earlier task transferred no data and submitted no job at lines
368–370.

The graph also now makes these boundaries executable:

- annotation is never a gate for genome-wide diversity or PSMC;
- Hi-C is evidence about long-range phasing, not a universal eligibility
  requirement;
- HiFi/base accuracy, QV, completeness, duplication or collapse, exact pair
  provenance, mutually comparable haplotypes, and measured callability govern
  core confidence;
- assembly-derived PSMC is descriptive demographic evidence from the same
  H1/H2 pair and is not independent validation of diversity from that pair;
- direct meiotic events, population-frequency gBGC, historical phylogenetic
  substitution bias, and non-allelic conversion among paralogs are four
  different estimands; and
- H1/H2-only observations are never direct evidence of conversion direction
  or biased transmission.

The machine-enforced audit is
`analysis/assert_vgp_comprehensive_wg.py`, with regression tests in
`analysis/tests/test_assert_vgp_comprehensive_wg.py`.

## Evidence reconciled

### Frozen release identity

Freeze 1 membership is anchored to the immutable catalog already recorded in
`analysis/vgp_phase1_freeze_provenance.json:121-135`:

| Field | Reviewed value |
| --- | --- |
| Repository | `VGP/vgp-phase1` |
| Commit | `dc1b2af5a7741b97d66fb10cb2bce97f41765cdf` |
| Catalog object | `VGPPhase1-freeze-1.0.tsv` |
| SHA-256 | `9c58420484a8b76a2d6175b7c26bf709e68bdc726a67fc7541b8c2b5a2fc13a4` |
| Object bytes | 327,466 |
| Physical lines | 717 |
| Data rows | 716 |
| Unique species | 714 |

The two excess data rows relative to unique species are legitimate catalog
multiplicities, not permission to deduplicate the release silently. Scale-out
must account for all 716 rows and separately dispose every linked H1/H2 pair.

Authoritative public references reviewed for source identity were:

- the [immutable VGP catalog commit](https://github.com/VGP/vgp-phase1/commit/dc1b2af5a7741b97d66fb10cb2bce97f41765cdf);
- the [official UCSC VGP assembly hub](https://hgdownload.soe.ucsc.edu/hubs/VGP/),
  which identifies itself as containing VGP-released primary assemblies and
  links each assembly to accession, BioSample, BioProject, and source data;
- the VGP flagship paper’s data-availability statement, which identifies
  GenomeArk and NCBI/EBI as archives and the UCSC VGP hub as the browser/hub
  resource; and
- the [UCSC download guidance](https://hgdownload.soe.ucsc.edu/goldenPath/help/ftp.html),
  which recommends rsync for large files and documents partial-transfer use.

The reviewed rsync root is
`rsync://hgdownload.soe.ucsc.edu/hubs/VGP/`. The current hub is a moving
transport/source. It does not define Freeze 1 membership; the pinned catalog
does. Mirror execution must first freeze a recursive source inventory and
reconcile its accession/version paths to the catalog’s closed-world set.

### Size claim correction

The approximately 967 GB whole-collection and approximately 520 GB FASTA-only
figures were retained only as historical planning estimates. They are not
treated as an authoritative inventory, checksum, acceptance threshold, quota,
or download ceiling. No remote listing or payload transfer was launched during
this quality task.

Before any mirror payload transfer, `mirror-vgp-freeze1` must create an
immutable metadata inventory with source endpoint, retrieval time, relative
path, object type, size, mtime, and link target. That inventory must provide:

1. exact file and byte totals for every Freeze 1 hub product in scope;
2. exact sequence-only totals, with FASTA and 2bit categories distinguished;
3. exact staging, durable, quarantine, checksum, and inode requirements;
4. storage/quota evidence and explicit operational headroom derived from those
   totals; and
5. a closed-world reconciliation from all 716 catalog rows to source objects.

This satisfies the need to verify size without enshrining a stale estimate or
an arbitrary laptop-scale total-byte ceiling.

### Checksum, resume, storage, and provenance contract

The mirror definition now requires all of the following:

- source-relative partial staging rather than writing into final paths;
- bounded cluster-friendly metadata/network concurrency and retry backoff;
- interruption-safe per-object resume;
- expected byte verification before promotion;
- an upstream digest where the frozen source publishes one;
- local SHA-256 after complete staging when no authoritative upstream digest
  exists;
- local SHA-256 revalidation immediately before and after atomic promotion;
- quarantine of mismatches without deleting the last verified object;
- preservation of useful source timestamps without accepting mtime as an
  integrity check;
- no unconstrained `--delete` against the durable mirror;
- mutually exclusive object and byte states: planned, transferred, verified,
  reused, missing, superseded, or quarantined; and
- exact source URL/rsync root, relative path, accession/version, catalog
  commit/hash, retrieval UTC, source size, checksum provenance, local SHA-256,
  license, and promotion state in the mirror manifest.

The pilot can consume any individually verified mirror object immediately.
It does not wait for the full mirror. Full-catalog scale-out and the
phylogenetic branch do wait for mirror completion. This is what “full official
Freeze 1 assembly mirror in the background” means operationally.

## Scientific correction: what H1/H2 can and cannot support

An H1 and H2 label alone proves very little. A usable core pair must be tied to
the same BioSample/individual or isolate, have exact accession versions and
roles, be mutually comparable assemblies of the two homologous haplotypes,
and carry enough assembly evidence to distinguish biological differences from
base error, collapse, duplication, gaps, and phasing switches.

When those conditions hold, aligning H2 to H1 yields candidate heterozygous
differences. The inference becomes a callable-genome estimate only after the
workflow excludes non-1:1 sequence and every other unsupported denominator.
The result is assembly-derived because its sensitivity is governed by how the
two assemblies and masks were made. Selective raw reads, k-mers, or published
estimates calibrate confidence; they are not mandatory for every core pair.

Hi-C helps establish chromosome-scale phasing and detect switches. It does not
by itself establish base accuracy, and its absence does not by itself prove a
pair unusable. Trio binning or other long-range phase evidence may supply a
different route. Conversely, Hi-C cannot rescue low QV, incomplete sequence,
collapse, duplication, wrong-individual pairing, or poor callability.

PSMC uses the spatial pattern of heterozygous and callable sites in a diploid
consensus. Therefore the graph requires a validated H1-coordinate diploid
consensus and mask, not annotation. Annotation partitions become eligible only
after the core result exists and only for an exact accession/version plus
sequence-dictionary match or a separately validated, manifest-bound liftover.

PSMC and genome-wide diversity from the same pair reuse the same aligned
haplotypes, heterozygous sites, and callable mask. The PSMC trajectory may
describe historical coalescence patterns relevant to interpretation, but it
is not an independent predictor or validation datum for that pair’s diversity.
Cross-species models must preserve that covariance.

## Ten-pair pilot contract

The design task must select exactly ten primary slots and at least five
pre-ranked alternates. The pilot is a ten-pair execution even if a unit fails:
every primary remains in the result ledger. An alternate is activated only by
a predeclared trigger and versioned manifest amendment; it does not erase or
silently replace the failed primary.

### Required stratification

The ten primaries must jointly cover:

| Axis | Required coverage |
| --- | --- |
| Major vertebrate clades | Mammals, birds, fishes, reptiles, and amphibians |
| Fish depth | Ray-finned and/or cartilaginous fish represented explicitly |
| Assembly generation | Earlier CLR/trio or primary/alternate generation and later HiFi/Hi-C hap1/hap2 generation |
| Genome size | Predeclared small, medium, and large strata from exact assembly lengths |
| Expected diversity | Predeclared low, medium, and high evidence strata; not inferred from pilot output |
| Pair confidence | Exact same-individual identity, reciprocal linkage, technologies, QV, completeness, duplication/collapse, and phasing evidence |
| Annotation | Exact-native available and unavailable pairs both represented so optionality is tested |

Each primary and alternate row must freeze the catalog row, TaxId, species,
BioSample, individual/isolate, H1/H2 accession versions, haplotype roles,
reciprocal pair evidence, assembly names/dates/generation/technologies, HiFi
and raw-read provenance, Hi-C/trio/other phase evidence, QV, completeness,
duplication/collapse, contiguity, exact assembly length, expected-diversity
source, annotation accession/status, URL, expected bytes, checksum policy,
license, retrieval state, and alternate trigger.

The quality task does not choose species. That is deliberately left to
`design-vgp-comprehensive`, which must reconcile current metadata and produce
auditable manifests rather than copy the earlier six annotation candidates.

## Exact core workflow contract

The verified repository handoff already establishes tool ownership and the
observed command order in `analysis/sweepga_impg_handoff.md:159-168` and the
biological proof at lines 296–352. The new implementation task must preserve
that contract at full-genome scale.

1. Orient H1 as reference and H2 as query.
2. Run SweepGA over both whole haplotypes with native
   `--num-mappings 1:1`.
3. Verify maximum query and target overlap multiplicity of one; reject any
   retained non-1:1 base.
4. Give the exact PAF to `impg index`.
5. Run `impg partition` and retain its native partition provenance.
6. Select only native partitions intersecting an eligible query region; never
   fabricate arbitrary substitute windows.
7. Run `impg query` with both sequence files and the selected native
   partitions.
8. Run `impg lace` to restore source coordinates and combine regional VCFs.
9. Normalize against H1 before exact target trimming; deduplicate deterministically.
10. Verify every REF against H1 and reconstruct the aligned H2 allele sequence
    from variants as a truth check.

SweepGA owns whole-haplotype mapping. IMPG owns index, partition, query, local
graph/VCF extraction, and lacing. `bcftools` owns normalization, trimming,
deduplication, and indexing. A task must fail if one tool silently assumes
another’s responsibility.

### Callable mask and denominator

The primary callable set is the H1-coordinate intersection of:

- verified whole-haplotype 1:1 alignment;
- eligible nuclear sequence;
- non-gap and non-N sequence on H1 and the aligned H2;
- resolved assembly breakpoints and coordinate continuity; and
- predeclared filters or sensitivities for base error, gross divergence,
  low complexity/repeats, duplication/collapse, sex chromosomes, organelles,
  unlocalized/unplaced sequence, and contig eligibility.

Every excluded base receives a reason code. The callable union plus the
reason-coded complement must reconcile exactly to the declared reference
universe without double counting. Repeat and other debatable filters should be
predeclared as primary or sensitivity masks, not changed after seeing PSMC.

The diploid consensus must encode validated heterozygous SNPs and indels in H1
coordinates and must not turn masked sequence into homozygous reference.
Variant/consensus checks must include REF validation, H2 allele reconstruction,
mask-denominator reconciliation, VCF normalization, and exact sequence-length
effects of indels.

### PSMC and bootstrap output

For every passing pair the workflow must retain:

- the exact masked diploid PSMC input and its digest;
- the unscaled PSMC trajectory;
- at least 100 block-bootstrap attempts;
- bootstrap units that do not cross contig or mask boundaries;
- bootstrap block-length sensitivity fixed before results;
- the number and fraction of finite successful replicates;
- mutation-rate scenarios with sources and uncertainty; and
- generation-time scenarios with sources and uncertainty.

Unscaled output is the primary reproducible object. Absolute time and effective
population size are scenario-derived and must not be reported without their
mutation and generation assumptions.

## Quantitative review gates

The design manifest must fix all quantitative thresholds before the pilot is
unblinded. The review task may apply them or tighten them with an explicit
technical reason; it may not relax them to rescue a result.

Hard zero-tolerance gates are:

- unresolved same-individual or H1/H2 accession identity;
- input or environment digest drift;
- retained query or target alignment multiplicity above one;
- unexplained mask-accounting discrepancy;
- H1 REF mismatch or H2 allele-reconstruction failure;
- annotation accession/sequence-dictionary mismatch for annotation outputs;
- a non-callable site encoded as homozygous reference in PSMC input; or
- conflation of unscaled and scenario-scaled PSMC.

Minimum program-level gates now encoded in `review-vgp-10-pilot` are:

| Gate | Required outcome |
| --- | --- |
| Primary accounting | 10/10 primary slots completed or failed explicitly |
| PSMC resampling | At least 100 attempts and at least 95% finite successful replicates per passing pair |
| Core GO yield | At least 8/10 primaries pass all pre-registered core gates |
| Coverage of design | Every required major clade and both assembly-generation strata remain among passers |
| Systematic bias | No assembly-technology stratum has systematic technical failure |
| Resource prediction | Absolute percentage error no worse than 25% at median and 50% at 95th percentile |
| Provenance/multiplicity/sequence truth | Zero hard-gate violations |

The exact callable-fraction floor and any assembly-specific minimum lengths
remain design-manifest parameters because one universal post hoc value would
be false precision across very different vertebrate assemblies. The review
must report measured callability continuously in addition to pass/fail.

Core, annotation, direct, population, phylogenetic, and non-allelic branches
receive separate `GO`, `CONDITIONAL_GO`, `NO_GO`, or
`NOT_RUN/DESIGN_ONLY` decisions. A missing annotation cannot change a core GO.
An absent independent Ne estimate cannot change a core GO.

## Gene-conversion evidence audit

The rewritten graph uses the following non-interchangeable claim matrix:

| Branch | Observations | Identifiable estimands | Not identifiable from branch alone | Current execution task |
| --- | --- | --- | --- | --- |
| Direct pedigree/gamete | Parent-to-offspring transmitted haplotypes, complete pedigrees, sperm or gametes | Callable meiotic event rate, tract length, crossover association, GC transmission distortion | Cross-vertebrate rate without transfer assumptions; historical substitution bias | `pilot-pedigree-gbgc` |
| Population frequency spectrum | Multi-individual same-population genotypes and polarized/polarization-aware WS/SW states | Model-dependent population gBGC strength/B and WS/SW frequency asymmetry | Direct event count, tract rate, parent-of-origin transmission | Design-only; no execution node in this graph |
| Historical phylogenetic | Ortholog-controlled alignments of close clades plus outgroups, including semi-complete genomes | Branch-specific WS/SW substitution asymmetry, clustered historical signatures, nonstationary gBGC-like strength | Direct events, current population B, pedigree tract rate | `pilot-vgp-phylo-gbgc` |
| Non-allelic/paralog | Copy-resolved paralogs or segmental duplications | Copy-homogenization/conversion candidates and tract-like patterns under paralog nulls | Meiotic allelic conversion or GC-biased transmission | Design-only; no execution node in this graph |

H1/H2-only data can count heterozygous WS/SW states and flag candidate
clusters. Without parent-of-origin/gamete transmission they do not identify
direction or transmission bias. Without multiple individuals they do not
identify an allele-frequency spectrum. Without an outgroup they do not
polarize ancestral and derived states. Without copy-resolved paralogy they do
not identify non-allelic conversion.

The current graph deliberately executes one direct pilot and two phylogenetic
clade pilots. `design-gbgc-evidence` must still design population and
non-allelic branches and their alternates, but synthesis must mark those rows
`NOT_RUN/DESIGN_ONLY` unless separate authorized tasks later produce conforming
manifests. This closes the previous loophole in which absent branches could be
filled rhetorically with H1/H2 summaries.

Semi-complete genomes are acceptable for historical phylogenetic inference
when callable orthology, fragmentation, missingness, alignment error,
ancestral-state uncertainty, and alternate outgroups are modeled. Completeness
affects denominators and confidence; it is not a requirement for direct event
claims because those claims are not made by this branch at all.

## Task-by-task graph audit

All timeout changes are worker hard timeouts, not runtime estimates or resource
ceilings. Longer data tasks received enough orchestration time to inventory,
resume, validate, and report rather than being killed at a laptop-oriented
deadline.

| Task | Scientific dependencies after review | Timeout | Required deliverable group | Key quality correction |
| --- | --- | ---: | --- | --- |
| `design-vgp-comprehensive` | `quality-vgp-psmc` | 1d | research plan, ten-primary manifest, alternate manifest, analysis JSON | Exact ten pairs; H1/H2 core validity; annotation/Hi-C boundaries; four evidence branches |
| `mirror-vgp-freeze1` | `design-vgp-comprehensive` | 14d | frozen source inventory, mirror manifest/summary/handoff | Pinned catalog defines release; current rsync hub is transport; inventory replaces 967/520 GB estimates |
| `implement-vgp-10-pilot` | `design-vgp-comprehensive` | 2d | code/tests, Guix/Slurm entrypoints, schema, handoff | Exact SweepGA→IMPG→mask→consensus→PSMC ownership and tests |
| `acquire-vgp-10-pilot` | `design-vgp-comprehensive` | 5d | acquisition/object/control manifests and handoff | Ten primaries, five alternates, exact pair provenance, selective raw reads, no silent substitution |
| `design-gbgc-evidence` | `design-vgp-comprehensive` | 1d | evidence plan, dataset/estimand/claim manifests | Four disjoint branches and explicit design-only status |
| `run-vgp-10-pilot` | implementation + acquisition | 7d | result/QC/telemetry manifests and report | All ten adjudicated; 1:1, masks, consensus, ≥100 bootstraps; measured resources |
| `review-vgp-10-pilot` | pilot run | 1d | independent review, gate TSV, decision JSON, corrected scale manifests | Pre-registered quantitative core/branch gates and independent recomputation |
| `scale-vgp-core` | independent review + completed mirror | 21d | closed-world scale manifest/QC/telemetry/results | Every eligible high-quality pair across all 716 rows; annotation subset only; no legacy global cap |
| `pilot-vgp-phylo-gbgc` | gBGC design + completed mirror | 7d | clade manifest/QC/results/report | Historical substitution signal only; semi-complete allowed with uncertainty |
| `pilot-pedigree-gbgc` | gBGC design + acquired direct control | 7d | direct dataset/tract/summary/report | Directional transmission is mandatory; H1/H2 is not a replacement |
| `synthesize-vgp-program` | scale-out + phylogenetic pilot + direct pilot | 2d | synthesis, claim ledger, final manifest, paper assets | Same-pair covariance and four-way claim ledger; absent branches stay design-only |

WG’s `.assign-*` lifecycle dependencies remain present and are checked for
referential integrity by the assertion script. They are excluded from the
scientific-edge table.

## Dependency rationale and background execution

The domain graph contains 17 exact scientific edges:

```text
quality-vgp-psmc
  -> design-vgp-comprehensive
       -> mirror-vgp-freeze1 -------------------+----> scale-vgp-core --+
       -> implement-vgp-10-pilot --+            |                      |
       -> acquire-vgp-10-pilot -----+-> run -> review -----------------+-> synthesis
       -> design-gbgc-evidence -----+-> pedigree pilot ----------------+
                                    +-> phylogenetic pilot <- mirror ---+
```

The actual graph has two important parallel paths:

- the full mirror starts after design and runs concurrently with pilot
  implementation and pilot acquisition; and
- the direct and phylogenetic specialized branches can proceed independently
  once their own inputs are ready.

Pilot compute intentionally does not depend on full mirror completion because
acquisition can reuse verified mirror objects or stage the exact pilot objects
into the same content-addressed contract. Scale-out depends on both independent
pilot review and full mirror completion. The phylogenetic pilot depends on the
full mirror because it may need several clade and outgroup assemblies. The
direct pilot depends on acquisition because that task owns its selected public
pedigree/gamete control.

No new population or non-allelic execution task was created. The user-assigned
graph named eleven tasks, and those two evidence types are currently design
and claim-boundary requirements. If the design recommends execution, new
authorized tasks must be added and synthesis must depend on them before their
results can be claimed.

## Resource and GNU Guix policy

The repaired six-species annotation pilot’s values—six species, 120 GiB
compressed input, 750 GiB scratch, 1,500 core-hours, two concurrent species,
256 GiB per job, and stricter run-specific I/O limits—remain valid provenance
for that refused scope. They are forbidden as global scientific eligibility
ceilings in this program.

The new graph instead requires:

- inventory-derived mirror storage and staging needs;
- per-pair estimates from exact object sizes and workflow scaling;
- measured pilot peak RSS, CPU, wall time, scratch, I/O, metadata operations,
  and retry behavior;
- per-job and per-wave scheduler safety limits with allocation headroom;
- resource-model validation at median and 95th-percentile errors; and
- bounded concurrency and stop conditions that protect the cluster without
  excluding large genomes by arbitrary global byte or memory caps.

Every design, build, acquisition, analysis, review, and synthesis definition
requires pinned GNU Guix. The already integrated full-suite result was 241
tests passing through the recorded Guix profile
(`analysis/repaired_vgp_integration_handoff.md:351`). This quality task adds a
read-only standard-library assertion; its regression test is part of the same
pinned full suite.

## WG graph assertions

Run the live graph validator from the repository root through the already
recorded pinned GNU Guix environment:

```sh
analysis/slurm/guix_job.sh \
  "$PWD/analysis/pilot_results/guix_environment.json" \
  python3 analysis/assert_vgp_comprehensive_wg.py
```

Expected output:

```text
VGP_COMPREHENSIVE_WG_ASSERTIONS_OK tasks=11 scientific_edges=17 pilot=10 alternates>=5 mirror=full-freeze1 gbgc=4-disjoint-branches
```

The validator asserts:

- all eleven controlled tasks exist;
- every dependency resolves and every `.assign-*` lifecycle dependency remains;
- exact forward scientific edges and reverse controlled-child edges;
- the implementation+acquisition join before pilot compute;
- the review+full-mirror join before scale-out;
- the scale+phylogenetic+direct join before synthesis;
- every reviewed timeout;
- required source, H1/H2, tool, mask, PSMC, annotation, pilot, quantitative
  gate, and gene-conversion claim text;
- every named deliverable;
- absence of the repaired pilot’s legacy global caps; and
- the no-bypass rules for core and specialized evidence.

The test module mutates a copy of the live graph to prove rejection of an
acquisition bypass, mirror timeout regression, annotation-as-core-gate
regression, missing directional-transmission language, reintroduced legacy
global cap, and missing machine-readable review decision.

## Validation checklist

- [x] Reviewed and edited all eleven definitions.
- [x] Reviewed and reasserted all 17 scientific dependencies and lifecycle
  edge integrity.
- [x] Reviewed and edited all eleven timeouts.
- [x] Named machine-readable and narrative deliverables for all eleven tasks.
- [x] Pinned the official Freeze 1 catalog and official UCSC rsync transport;
  made live inventory reconciliation a pre-transfer hard gate.
- [x] Reclassified 967/520 GB as historical estimates pending exact frozen
  inventory; required exact object/byte counts, storage headroom, checksums,
  resume, quarantine, atomic promotion, and provenance.
- [x] Required exactly ten primaries, at least five alternates, clade,
  assembly-generation, genome-size, and diversity strata, and exact pair
  provenance.
- [x] Locked SweepGA `--num-mappings 1:1`, bidirectional multiplicity, exact
  IMPG partition/query extraction, mask, consensus, PSMC, bootstrap, scenario,
  and annotation-accession contracts.
- [x] Removed legacy laptop-scale global memory/byte ceilings from the new
  program and added an automated regression assertion.
- [x] Separated direct, population, phylogenetic, and non-allelic claims and
  stated all relevant non-identifiability boundaries.
- [x] Preserved annotation as optional for core analysis and Hi-C as phasing
  evidence rather than a universal gate.
- [x] Preserved same-pair PSMC/diversity dependence in review and synthesis.
- [x] Added a read-only WG graph validator and mutation-based regression tests.
- [x] Launched no biological download and no Slurm job.

## Authorization boundary

This quality task authorizes the downstream design task to start when WG
completes its lifecycle. It does not itself authorize mirror transfer, pilot
acquisition, raw-read acquisition, Guix realization, or Slurm execution.
Each downstream task must satisfy its own immutable manifests, preflight gates,
licenses, storage evidence, and branch-specific decision before taking the
action named by that task.
