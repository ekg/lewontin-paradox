# Vertebrate scale-out execution program

Date: 2026-07-17 UTC

WG task: `synthesize-vertebrate-scaleout-plan`

Status: **reviewed specification and proposed graph only; no execution authorization**

## 1. Executive decision

The program is ready to freeze local schemas and review already-held metadata.
It is not authorized to retrieve the missing VGP catalog, download biological
assets, realize a new shared Guix closure, submit a Slurm job, run a pilot,
perform demographic inference, or expand to the full catalog. Those actions
remain behind the human-decision nodes in
`analysis/vertebrate_scaleout_execution_graph.tsv`. Silence, completion of a
planner, existence of an `sbatch` script, or success of a predecessor is a
**no**, not implied approval.

The synthesis preserves two different scientific programs:

1. assembly/composition and assembly-pair diversity for exact, versioned,
   original H1/H2 assets using Guix-built SweepGA and IMPG; and
2. population diversity, independent predictors, and method-valid demographic
   subsets whose inputs cannot be inferred from the assembly catalog.

The programs join only at provenance-aware statistical synthesis. A deposited
H1/H2 pair is not a population sample or a demographic input. A PSMC, MSMC2,
SMC++, SFS, or site-pattern result made from diversity-bearing observations is
not an independent predictor of that diversity. Absolute effective population
size or calendar time is released only as a calibration-scenario distribution
with mutation-rate and generation-time provenance and propagated uncertainty.

### 1.1 Mandatory dependency proof and frozen planning inputs

`wg show synthesize-vertebrate-scaleout-plan` was manually checked on
2026-07-17. Both `plan-all-vertebrate-tier3` and
`plan-vertebrate-ne-strategy` are direct, completed dependencies (the automatic
assignment task is an additional dependency). The following completed
deliverables were read and reconciled; none may be bypassed:

| branch | artifact | SHA-256 at synthesis |
|---|---|---|
| synchronized origin audit | `analysis/vertebrate_scaleout_origin_inventory.md` | `e7ad18e6b8fceb79b9fb58e4ae5adfc345b4ccb6ce507f85a0d0e4069b31313b` |
| synchronized guidance | `results/tier3/vgp_freeze_analysis.md` | `0ad6fa03ceeec9d07c39c5456ddc4702c54c586115e578af3599d07f29e5316d` |
| Tier 3 plan | `analysis/vertebrate_scaleout_plan.md` | `22d307842df84b54d24ceb24e0290f560730a0e08b8fb3dccfd6b9de0597731c` |
| Tier 3 schema | `analysis/vertebrate_scaleout_candidate_schema.tsv` | `575f35f84daaa4bb476f2692f4cae860a5f4f1fa864f05efc2b9e36bd04b9ae0` |
| Tier 3 budget | `analysis/vertebrate_scaleout_resource_budget.tsv` | `d0659b1dfc29fd078de7b5f6f82399fb43d125d96403b88a774515d4539e7692` |
| Tier 3 graph | `analysis/vertebrate_scaleout_wg_graph.md` | `0efdef977e0d7b2fd4a85324fdb7c149729666a86c2d0aecdf87ecc7f53fce97` |
| Ne strategy | `analysis/vertebrate_ne_strategy.md` | `a6967a169ca58a38f53b6a899620e0ee603c4549367eb13c0943ff8b8baec643` |
| Ne methods | `analysis/vertebrate_ne_method_matrix.tsv` | `d158f3f6e7f8b13d7a8267ea092e47da72df75537280f7dcb2d131dd7af3c905` |
| Ne source schema | `analysis/vertebrate_ne_source_schema.tsv` | `60561382d0bec030defc3e74a4382345af6741c3a1e1207e3a50e6dae194af05` |
| Ne graph | `analysis/vertebrate_ne_execution_graph.tsv` | `fd6e461b37d56128b668af5247f879a30da57bdfc3ba1c7ea696048a45079c32` |

The quality-pass artifact is advisory review evidence, not a replacement for
the synchronized guidance or either dependency. Its prior off-main history is
retained in the origin inventory; the current tracked file was not used to
waive any gate. The contradiction and every other material reconciliation are
recorded in `analysis/vertebrate_scaleout_decisions.tsv`.

### 1.2 Executive stage classification

Only the three status strings below are permitted. “Proceed” always names a
hard cap; it never means unlimited work.

| stage/action | status | owner | recorded scope and evidence required |
|---|---|---|---|
| 0: freeze schemas, IDs, overlap policy, and this proposed graph from local files | proceed within recorded bounded scope | Program Architect + Data Provenance Reviewer | ≤8 core-h, ≤4 CPU, ≤16 GiB, ≤8 wall-h, ≤2 GB scratch/persistent, zero network; hashes above and clean schema validation |
| 1a: retrieve the missing immutable VGP TSV | requires explicit user authorization | Human project owner | exact commit/object URL, ≤1 MB, one object, retrieval window and license record in `A00` |
| 1b: metadata/API inventory after 1a | proceed within recorded bounded scope | Inventory Curator | ≤30 core-h, 4 concurrent metadata requests, ≤16 GiB/job, ≤24 wall-h, ≤4 GB scratch, ≤2 GB persistent, no FASTA/GFF/VCF/read download |
| 2a: local lint, fixture validation, and dry-run command rendering using existing profiles | proceed within recorded bounded scope | Reproducibility Reviewer | ≤8 core-h, ≤16 GiB, ≤8 wall-h, ≤5 GB scratch, zero biological download and zero Slurm submission |
| 2b: realize/build a new shared Guix closure or run compute-node smoke | requires explicit user authorization | Human project owner + Reproducibility Reviewer | exact manifests/derivations and ≤30 GB store-growth envelope in `A20` |
| 3: stage and run the eight-slot assembly/composition pilot | requires explicit user authorization | Human project owner | separate exact-asset `A30` and compute `A31` decisions; hard small-pilot cap in Section 7 |
| any pilot continuation or retry that would exceed the small-pilot cap | requires explicit user authorization | Human project owner + HPC Storage Reviewer | overage task `A32` records delta and revised rollback; the default pilot cannot spend it |
| 4: independent review of a completed pilot packet | proceed within recorded bounded scope | Independent Scientific, Circularity, Reproducibility, and HPC Storage Reviewers | ≤16 core-h, ≤32 GiB, ≤24 wall-h, ≤10 GB; read-only review, no automatic continuation |
| 5: each expansion wave | requires explicit user authorization | Human project owner + three reviewers | distinct `A50`, `A5N` decision for exact rows and telemetry-refit envelope; never blanket approval |
| 6: full-catalog submission | requires explicit user authorization | Joint human panel | distinct `A60`; cumulative assembly high bound is a ceiling, not a request |
| 7a: bounded public bibliographic/independent-covariate curation | proceed within recorded bounded scope | Demographic Source Curator | 1–4 CPU, ≤16 GiB, ≤40 wall-h plus manual review, ≤20 GB; licensed/bulk exports excluded |
| 7b: bulk/very-large population or read acquisition | requires explicit user authorization | Human project owner + Data Steward | exact URLs/accessions, terms, bytes, bandwidth, storage and retention in `A70` |
| 7c: PSMC/MSMC2/SMC++/SFS/BPP or expensive LD/temporal/pedigree inference | requires explicit user authorization | Human project owner + Demographic Methods Reviewer | separate `A71`, ≤2,000 core-h pilot ceiling; no assembly authorization carries over |
| demographic catalog-scale or unbounded method expansion | deferred | Human project owner | no defensible totals or frozen eligible inputs exist; reconsider only after reviewed pilot telemetry |
| 8: frozen primary fit and release | requires explicit user authorization | Human project owner + Independent Results Reviewer | `A80` analysis decision then `A81` release decision; no unreviewed or circular joins |

## 2. Scientific questions, estimands, and analysis strata

### 2.1 Preregistered hypotheses

- **H1 (primary Lewontin-paradox comparison):** population-level nucleotide
  diversity is associated with independently measured census abundance,
  range, density, life history, and separately justified mutation-rate records
  after phylogeny, measurement error, sampling, and missingness are modeled.
- **H2 (assembly/composition descriptive):** exact-H1 assembly-pair diversity
  and original-H1-native GC3 covary across eligible species. The estimand is a
  descriptive phylogenetic slope; it is not a population-Ne or causal effect.
- **H3 (demographic context):** native coalescent-scaled curve shapes vary
  across eligible taxa and methods. These are genomic outcomes and can support
  historical context and robustness comparisons, not independent explanation
  of the response diversity.
- **H4 (secondary independent demography):** population/date-matched LD,
  temporal, pedigree, or sibship estimates that use no response samples/sites
  may covary with population diversity. Each estimand remains method-specific;
  there is no generic `Ne` coefficient.

### 2.2 Four non-interchangeable strata

| stratum | estimand and sampling unit | conjunctive eligibility | uncertainty and missingness | phylogenetic/statistical design | circularity | claims permitted / prohibited |
|---|---|---|---|---|---|---|
| 1. Vertebrate-wide assembly/composition | GC3 numerator/queryable codons; CDS/coding/4D composition; where a true pair passes, `diploid_haplotype_assembly_diversity = H1/H2 alternatives / callable H1 bases`; unit is one exact H1 assembly or one same-individual phased H1/H2 pair, summarized once per species | immutable H1; for diversity exact versioned homologous H2 and same individual; original H1-native GFF; sequence-set/contig bijection; reconstructed CDS/genetic code; whole-FASTA SweepGA 1:1; IMPG calls; positive original/queryable/callable denominators | block-bootstrap SE/CI; mapping/collapse/phase sensitivity; structural exclusions and zero denominators are missing, never zero; calibration rows excluded from scientific cohort | species unit; clade-balanced median/IQR; Brownian PGLS of log10 CDS diversity on logit GC3; Pagel-λ/OU/mixed/order-robust and leave-clade-out sensitivity; inferential fit only at ≥20 species, ≥4 classes, ≥8 orders | not an Ne predictor; assembly diversity remains distinct from population pi | may describe exact assemblies and a cross-species association; cannot claim population diversity, contemporary/long-term Ne, demographic history, or causation |
| 2. High-quality population-data subset | population pi or site-class-specific diversity with callable denominator; unit is a biologically delimited population sample, with multiple populations nested in species | exact reference accession/hash/dictionary; versioned VCF/BCF or all-sites evidence; population/sample/date/geography; unrelatedness/ploidy/sex; filter provenance; reference-matched callable mask; adequate sample count and simulation-derived callability | block/jackknife/bootstrap or genotype-likelihood uncertainty; within-species covariance; no mean imputation; failures classified by reference, mask, sample, QC, resource, or access | Bayesian phylogenetic errors-in-variables model with species random effect and tree ensemble; complete-case PGLS baseline; clade, missingness, quality, population/date and callable-site sensitivities | response stratum; predictors must have `response_dataset_overlap=none` for primary fit | may support comparative population-diversity associations for represented populations; cannot generalize pooled geography, treat sparse absent VCF sites as invariant, or mix assembly-pair pi as the same response |
| 3. Coalescent-scaled and scenario-scaled histories | PSMC pairwise history from one callable diploid individual; MSMC2 within/cross coalescence rates from ordered phased haplotypes; SMC++/SFS population histories from valid cohorts; BPP theta/tau for audited loci; unit is method × individual/population/comparison/clade | every method-specific gate in the method matrix; exact reference/masks; callable invariant sequence; phase/switch evidence where required; population assignment; model/optimization/convergence validation; no H1/H2 shortcut | retain native scaled parameters, bootstrap/chains/starts and identifiable windows; absolute scenario draws propagate mu, generation time, recombination and transfer uncertainty; invalid/nonidentifiable rows remain in ledger | descriptive curve/shape clustering and method agreement; separate labeled shared-genomic sensitivity only; no entry into Tier-1 primary design matrix | genomic-shared unless independently measured from disjoint data; PSMC/MSMC2/SMC++ and SFS derived from response-bearing data are not independent | may describe conditional historical shape and scenario ranges; cannot be a causal/independent predictor of the same pi, contemporary Ne by default, or an absolute history without calibration provenance |
| 4. Independent ecological/life-history or measured demographic predictors | census/range/density/body mass/generation/fecundity/longevity/mode/dispersal/direct mutation measurement, or method-specific independent `Ne_LD`, `N_eV`, pedigree/variance/inbreeding Ne or `N_b`; unit is source × estimand × taxon × population × time interval | original DOI/URL/version/hash/license and locator; exact taxon concept; native units/definition; uncertainty/sample effort; geography/date; population match; response overlap audit; no selection conditioned on pi | source intervals/posteriors enter measurement model; interval censoring/latent values; taxonomic/geographic/temporal transfer variance; no zero/mean imputation; availability modeled and reported by clade/traits without selecting on pi | primary Tier 1 phylogenetic errors-in-variables model; each independent demographic estimand enters a separate secondary model; tree ensemble, leave-clade-out, grade/transfer/captive/stale-record sensitivities | primary requires overlap `none`; `unknown` fails closed; disjoint population/sample records only in preregistered secondary tiers | may support conditional comparative association/prediction; cannot prove causation, convert census to Ne silently, or collapse unlike demographic quantities into a generic Ne |

### 2.3 Demographic-covariate tiers

| tier | disposition and claim boundary |
|---|---|
| independent measured demographic | Direct raw/reanalyzed LD, temporal, pedigree/sibship or variance/inbreeding estimate with exact population, dates, samples, loci, method, interval, and no response overlap. Secondary predictor; primary only if preregistered independence and alignment are exact. |
| literature/ecological Ne | Same evidence burden even when copied from a paper; audit underlying samples/sites. Retain the paper's estimand code and measurement window. Generic or untraceable `Ne` is rejected. |
| coalescent-scaled genomic summary | Native PSMC/MSMC2/SMC++/SFS/BPP parameters and prespecified curve features. Descriptive/shared-data sensitivity only. |
| scenario-scaled absolute history | Coalescent result plus foreign-keyed mutation-rate, generation-time, ploidy/heredity and recombination assumptions drawn jointly. Conditional scenario output, never an independent predictor. |
| excluded/circular | `pi/(4*mu)`, algebraic equivalents, theta or curves from the response sites/samples offered as predictors, invalid H1/H2 conversion, hidden borrowed calibration, unknown overlap, or unresolved taxon/population/reference. Retain only in rejection/missingness audit. |

For absolute output, record the exact scaling formula, genomic compartment and
ploidy; per-generation mutation-rate distribution and assay; generation-time
definition and distribution; population/taxon and transfer distance; any
mutation–generation covariance; recombination map/rho treatment; calibration
record IDs; and all scenario draws. If those fields fail, the scaled result is
retained and absolute Ne/calendar time is absent.

## 3. Synchronized-guidance implementation crosswalk

The IDs are inherited from the Tier 3 plan. “Gate” means fail closed; a review
cannot substitute a newer, projected, alternate, or more convenient asset.

| ID | controlling requirement | implementation stage | exact eligibility gate | required artifact | validation decision |
|---|---|---|---|---|---|
| G01 | named freeze TSV and 717 lines | 0/1 | immutable source object | `source_catalog.json`, frozen TSV | 717 lines/716 rows exactly or stop |
| G02 | guidance says committed but file is absent | 0/1 | absence/contradiction explicit | source reconciliation record | user-authorized retrieval only; never assume local |
| G03 | moving `main` and placeholder date invalid | 1 | commit/object, real UTC, bytes/hash | retrieval provenance | reject moving revision or placeholder |
| G04 | reported 714/223/248/271/120/40 | 1 | deterministic parser and closed ledger | `catalog_counts.json` | all reproduce or discrepancy review stops |
| G05 | 13/4/3/9/6/5 taxonomic split | 1 | TaxId/lineage snapshot | taxonomic count table | exact sum/group match or stop |
| G06 | columns 10/13/16/17/21/26 control screen | 0/1 | ordinal plus header assertion | parser/schema hash | schema drift/duplicate ambiguity rejected |
| G07 | 13 fish H1 accessions/H2 labels | 1 | all seeds retained; labels not accessions | row inventory | none omitted; unresolved labels remain blocked |
| G08 | 46 composition fish lack row list | 1 | derive from frozen catalog | exact 46-row ledger | reproduce 46 or stop |
| G09 | exact H1/H2, reciprocal same-individual roles | 1/3/5/6 | versioned accessions and pair proof | `pair_evidence.json` | 100% fields/proof or reject diversity |
| G10 | provider and repository hashes | 1/3/5/6/7 | provider checksum plus SHA-256 and sequence-set digest | asset manifest | all agree; mismatch quarantines lineage |
| G11 | pair label is not callable evidence | 1/3/5/6/7 | positive H1-reference callable/queryable denominator | mask/denominator packet | zero/unknown denominator rejects result |
| G12 | annotation identical to H1 content | 1/3/5/6 | declared original native release, sequence-set identity and contig bijection | annotation linkage packet | any disagreement rejects row |
| G13 | no projected/mismatched primary annotation | all compute | original H1-native only | native/projected audit | projected/lifted/newer substitute rejected |
| G14 | genetic code/regions/CDS reconstruction audited | 1/3/5/6 | valid code, phase, codons and nonzero CDS | CDS audit/target BED | any frame/reference failure rejects composition |
| G15 | three completed fish are calibration only | 0/8 | `record_role=calibration_only` anti-join | cohort ledger | never count in freeze estimands |
| G16 | pinned Guix channel/manifest control | 2 onward | exact hashes, derivations, closure/store path | environment record | no ambient/container/module substitute |
| G17 | baseline lacks IMPG | 2 onward | reviewed supplemental combined profile | IMPG/SweepGA build packet | ambient/unrecorded IMPG is fatal |
| G18 | SweepGA whole H1/H2 native 1:1 | 2/3/5/6 | whole FASTAs, `--num-mappings 1:1`, multiplicity ≤1 | bounded PAF and argv | multiplicity/coverage audit must pass |
| G19 | `--scaffold-jump 0` assigns partitioning to IMPG | 2/3/5/6 | command fingerprint contains exact option | mapping record | absent/different option rejects lineage |
| G20 | IMPG owns index/partition/query | 2/3/5/6 | native partitions intersect original H1 targets | graph index, partitions, focus BED, regional VCF | other partitioner/target transfer rejected |
| G21 | IMPG VCF; bcftools normalized VCF/BCF/index | 2/3/5/6 | REF=H1, queryable TBI/CSI, same normalized records | VCF.gz/TBI and BCF/CSI | 100% REF and record concordance |
| G22 | original H1-native targets/denominators | all summaries | no alternate annotation or target lift | target/denominator packet | lineage anti-join rejects substitutes |
| G23 | metadata preflight before bulk | 1 before 3/5/6/7 | exact size/terms/eligibility before transfer | frozen preflight manifest | no biological staging until authorization |
| G24 | bounded staging and retained failures | 3/5/6/7 | immutable row list and closed disposition | asset/failure ledgers | accepted+rejected+unavailable+failed=input |
| G25 | truth/stratified pilot before expansion | 2/3/4 | fixture truth then eight-slot strata | pilot review packet | all quantitative gates in Section 7 |
| G26 | dependency gates, scratch, atomic promotion, retry, sacct | 3/5/6/7 | run fingerprint and telemetry complete | run/retry/telemetry ledgers | missing telemetry or partial artifact stops |
| G27 | Tier 3A/Tier 3C eligibility separate | 0 onward | independent modality flags | candidate ledger | no generic `ready` selector |
| G28 | population/demography not inferred from H1/H2 | 0/1/7 | deposited method-valid data and masks | method eligibility ledger | assembly status never promotes population/demography |
| G29 | every row collected | every collector | closed-world reconciliation | complete status ledger | exact set equality required |
| G30 | independent QC/resource review before full catalog | 4/5/6 | reviewer packet and zero bypass | signed decision record | `A60` remains non-ready until review |
| G31 | project ID is not a row asset | 1/3/5/6/7 | versioned accession/URL/hash per object | retrieval manifest | PRJNA489243 alone is insufficient |
| G32 | planning is not approval | every boundary | human decision with scope/budget/timestamp | WG decision record | no executable descendant ready without it |

The population/demographic extension adds equally hard gates: exact
sample/reference/mask identity; PSMC novel-call coverage initially ≥15× plus
downsampled heterozygote-FNR audit; novel cohorts initially 10–15× with ≥90%
individual callability at retained sites, with final thresholds frozen by
simulation; measured MSMC2 phase/switch error; SMC++ proof that absent sparse
VCF sites are distinguishable from callable invariant sequence and that
FILTER/QUAL was applied upstream; optimization/chains/bootstraps; and
`response_dataset_overlap=none` for primary predictors.

## 4. Frozen software and provenance contract

### 4.1 Assembly/composition lineage

The controlling Guix channel is commit
`44bbfc24e4bcc48d0e3343cd3d83452721af8c36`; `analysis/guix/channels.scm`
SHA-256 is
`45c055cd1d9010a72eacbb720037a22bccb2d8d6891dbd11b5d66365f29b3a17`.
The baseline manifest SHA-256 is
`2fb05e87aa2ac45ce51d4dcf93b232cb98627f525adace98357629ee3f15720a`
and its recorded profile is
`/gnu/store/z9v2f6faha9cwjz0sm5iphhlzisgi077-profile`.

Production assembly work additionally requires the reviewed combined packet:

- SweepGA manifest SHA-256
  `ea9ae1ba3e51ac3302d93add158532befec8fb3c09d188f524ac29237bab17d1`,
  SweepGA commit `018e4ce49d2c125820e0ac50dc5feaa02d423683`, accepted
  binary SHA-256
  `fa7f0edb9b7e275c288db254046020e136d4267dd5ee043379227ef80da0573b`;
- WFMASH commit `e040aa10e87cab44ed5a4db005e784be62b0bd21`;
- combined runtime manifest SHA-256
  `c0ef5afd6c988341da8446ff3f70af274dd12f5514bc053d3d4e6f0cbdcee521`;
- IMPG commit `101df81eb28a809c8fac97d297acd9fcfbbfa048`, binary
  SHA-256
  `c587dc2326cd24f887b1fcb3938404229ad0f0a27ef0773e90c287b1ade160d4`,
  `gfaffix` `460e0dd798a9da7d12aef4f9181419d71489da95`, and `syng`
  `dd00f52b688c0fb78cb7f25336ef9ac9f6a3e109`;
- bcftools 1.14 from its recorded store path.

Two clean SweepGA builds must be byte-identical and match the accepted hash.
The packet records source/submodule archives, lock files, patches, commands,
build logs, binary hashes/realpaths/linkage, `.drv` paths, resolved channels,
closure paths/hash/size, GC root, license, login and compute smoke, and fixture
results. The baseline's absence of IMPG is intentional; discovering an ambient
IMPG does not satisfy the supplemental lineage.

### 4.2 Demographic environments

PSMC, MSMC2 (never legacy MSMC), SMC++, phasing/helpers, LD/temporal/pedigree,
SFS and BPP each receive a separate reviewed
`analysis/guix/ne/<method>-manifest.scm` and package definition. A tag named
`latest`, ambient Conda, a module, or an opaque container is rejected. Each
method record has the same channel/source/derivation/closure/binary/license
contract plus deterministic synthetic truth and failure fixtures. No method
profile may be treated as covered by the baseline merely because Python or
bcftools is present.

### 4.3 Object-level provenance and supersession

Every object is addressed by `(inventory_release, lineage_id, run_id,
candidate_id, attempt_id, logical_role)` and records accession.version, exact
URL/revision, retrieval UTC, provider checksum, compressed and uncompressed
SHA-256, sequence-set SHA-256 where relevant, byte size, permissions, source
terms, producing argv JSON/hash, environment hash, input/output hashes, and
validation status. References also carry `.fai` and contig-dictionary hashes;
annotations carry release, declared reference, native/projected status,
sequence-region map, genetic code and CDS audit; population objects carry
sample-list, population/date/geography, filters, mask and callable denominator.

Consumers use only the signed manifest, never wildcards or “latest” paths.
`superseded.tsv` marks every old/mismatched/legacy object `consumable=no` with
replacement and reason. A later accession, alternate haplotype, projected GFF,
reference-mismatched VCF/mask, stale sentinel, unlisted file, or `.partial`
object cannot be selected to rescue a row.

## 5. Workflow, Slurm, I/O, and recovery contract

Generated shell or `sbatch` text is inert. This synthesis did not invoke it.
After authorization, the assembly path is exactly:

```text
immutable whole H1 + whole H2
 -> SweepGA/WFMASH --num-mappings 1:1 --scaffold-jump 0
 -> bounded PAF + multiplicity audit
 -> IMPG index -> native partition
 -> intersection with original H1-native targets
 -> IMPG regional query/lace VCF
 -> bcftools normalized VCF.gz/TBI + BCF/CSI
 -> original/queryable/callable denominator summaries
```

The array topology is stage → exact preflight → compute, with `aftercorr` for
corresponding row elements. Collectors use `afterany` so failure cannot vanish;
release/evaluation uses `afterok` only after the collector proves a complete
ledger and zero unexplained failure. Tier 3A, Tier 3C, population calling and
each demographic method have distinct manifests, result directories and
arrays. Their reviewed collectors are explicit joins.

Inputs are checksum-verified and read-only on MooseFS, then copied once to
`$SLURM_TMPDIR/<run>/<candidate>/<attempt>`. Whole genomes are never streamed
repeatedly from MooseFS. Work, indexes, checkpoints and intermediate chains
remain local. Promotion writes `.partial.<jobid>` on the destination
filesystem, validates format/denominators/reference, computes SHA-256, fsyncs
file and directory, and atomically renames into a new run directory. A
cross-filesystem promotion is copy → checksum compare → destination rename.
Accepted objects are never overwritten.

A success sentinel is validated JSON containing schema/row/command/environment
fingerprints, all hashes, Slurm IDs, attempt, times, status and validations.
Resume skips only an exactly matching sentinel. All other work is quarantined
into a new attempt. Every job records five-second application/cgroup CPU, RSS,
I/O and scratch high-water samples plus `sacct` fields `ElapsedRaw`,
`TotalCPU`, `MaxRSS`, `MaxDiskRead`, `MaxDiskWrite`, `AllocCPUS`, `ReqMem`,
state and exit. Blank scheduler fields require the validated cgroup/getrusage
sidecar, not an estimate.

Retries are limited as follows:

- network/node/filesystem transient: two automatic retries with exponential
  backoff;
- OOM/time: one reviewer-approved outlier retry at no more than 2× the row's
  authorized memory/time and never beyond the stage cap;
- checksum/reference/native-annotation/denominator/sample/taxon/circularity or
  scientific failure: zero automatic retries and no alternate asset;
- code/command/environment mismatch: quarantine the entire lineage;
- demographic identifiability/mask/phase failure: one prespecified remediation
  after independent review, then stop that method.

`retry_manifest.tsv` contains only failed frozen row fingerprints. Successful
elements are not rerun unless a recorded supersession invalidates them.
Cleanup of node-local scratch occurs after verified promotion; failed partials
are held at most 14 days for triage. Reproducible PAF/index/regional
intermediates are held through independent QC plus 90 days. Source catalog,
provenance, Guix closure records, final indexed calls, denominators, results,
failure ledgers, sentinels and telemetry are held through publication plus
seven years subject to institutional policy. Compressed source assets remain
through final reproducibility review and may be evicted only by a
checksum-listed, human-approved cleanup with a tested rehydration recipe.

## 6. Resource estimates and assumptions

### 6.1 Evidence base

Assembly low/base/high values below reproduce
`analysis/vertebrate_scaleout_resource_budget.tsv`. The only whole-pair
calibration comprises three fish: 0.622–2.404 allocated core-h at 8 CPUs, with
0.296–0.786 GB output and 1.975–2.611 GB staged input. The 135-row Tier 3C
batch used 18.223 allocated core-h, mean 3.398 GiB and maximum 15.631 GiB RSS.
Tier 3A actual RSS/CPU efficiency, scratch, MooseFS bytes/operations and
bandwidth were not captured. Thus all catalog totals are planning bounds until
the stratified pilot collects those dimensions and refits a robust log-linear
model using genome size, contig count/N50, annotation bytes/CDS rows, mapping
records and callable fraction. Scheduling retains the upper prediction
interval and reports leave-one-out error.

The inherited eight-CPU Tier 3A scheduling equation is retained verbatim as a
provisional bound, not promoted to a fitted biological model:

```text
W_hours = k * [0.033 + 0.083*G + 0.067*max(0, M-1)] * F_C * F_A * F_Q
core_hours = 8 * W_hours
```

Here `G` is H1 Gbp, `M` is bounded PAF records in thousands, and `F_C`,
`F_A`, and `F_Q` are fragmentation, annotation/target-size, and inverse
callability/queryability pressure factors. Missing factors remain flagged; they
are never silently imputed. Base storage assumes
`persistent_input = 1.6 * source_bytes`,
`local_scratch = 3 * source_bytes`, and about 0.5 GB persistent output per
species, where `source_bytes` is uncompressed H1 + H2 + GFF. Pilot measurements
replace these factors before expansion.

The demographic branch has zero frozen vertebrate inputs and no local pilot
telemetry. Its low/base/high values are **provisional envelope partitions** of
the upstream hard pilot cap, using the upstream 2–4× scratch/download staging
and 3–6× read/write amplification ranges. They are not quota requests. Phase 7
metadata preflight must replace every number from provider sizes and synthetic
power/resource runs before `A70` or `A71` can be approved.

### 6.2 Conditional low/base/high totals

All units are decimal GB except demographic scratch/download also show the
binary-cap interpretation. Elapsed time is machine elapsed under the stated
concurrency and excludes human waits. Peak memory is concurrent aggregate,
with the per-job maximum in parentheses. “Program” means the reported
assembly catalog plus the bounded demographic pilot, executed sequentially;
it is not a demographic catalog-scale estimate.

| envelope | core-h | peak memory GiB | machine elapsed | peak local scratch | persistent storage | files/inodes | MooseFS reads | MooseFS writes | metadata operations | bulk download | peak bandwidth and concurrency |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| assembly low: inventory + reported 40 Tier3A/120 Tier3C | 28 | 128 aggregate (32/job) | 0.8 h | 24 GB | 128 GB (120 input + 8 output) | 45,000 | 123 GB | 83 GB | 255,000 | ≤120 GB | 50 MiB/s; Tier3A 4, Tier3C 16, transfers 2 |
| assembly base | 116.2 | 640 aggregate (64/job) | 1.7 h | 111 GB | 425 GB (400 + 25) | 100,000 | 370 GB | 190 GB | 970,000 | ≤400 GB | 120 MiB/s; Tier3A 10, Tier3C 16, transfers 2 |
| assembly high | 1,482 | 960 aggregate (96/job) | 19.4 h | 495 GB | 1,672 GB (1,500 + 172) | 320,000 | 1,670 GB | 1,020 GB | 3,680,000 | ≤1,500 GB | 200 MiB/s ceiling; Tier3A 10, Tier3C 16, transfers 2 |
| demographic pilot low | 300 | 128 aggregate (64/job) | 72 h | 512 GB | 75 GB | 20,000 | 750 GB | 250 GB | 100,000 | 256 GB | 40 MiB/s; compute 2, transfer 1 |
| demographic pilot base | 1,000 | 512 aggregate (128/job) | 240 h | 1,024 GB | 180 GB | 75,000 | 3,000 GB | 1,200 GB | 500,000 | 1,024 GB | 80 MiB/s; compute 4, transfers 2 |
| demographic pilot high hard cap | 2,000 | 1,536 aggregate (256/job) | 504 h | 2,048 GB (2 TiB) | 300 GB | 200,000 | 8,000 GB | 3,000 GB | 2,000,000 | 2,048 GB (2 TiB) | 100 MiB/s; compute 6, transfers 2 |
| conditional program low | 328 | 128 aggregate | 72.8 h | 512 GB | 203 GB | 65,000 | 873 GB | 333 GB | 355,000 | 376 GB | 50 MiB/s; stage-specific caps above |
| conditional program base | 1,116.2 | 640 aggregate | 241.7 h | 1,024 GB | 605 GB | 175,000 | 3,370 GB | 1,390 GB | 1,470,000 | 1,424 GB | 120 MiB/s; stage-specific caps above |
| conditional program high | 3,482 | 1,536 aggregate | 523.4 h | 2,048 GB | 1,972 GB | 520,000 | 9,670 GB | 4,020 GB | 5,680,000 | 3,548 GB | 200 MiB/s; stage-specific caps above |

Catalog species counts are reported, not reproduced; asset overlap is assumed
checksum-deduplicated. Persistent inputs are used as the conservative assembly
download proxy. Demographic provisional reads/writes apply 3–6× amplification
to staged bytes and bootstrap/model reuse. File counts assume attempt-scoped
manifests, masks, indexes, bootstrap summaries and telemetry; the pilot must
measure them. If provider bytes, file counts, or any high scenario exceeds
available quota, the corresponding authorization is no-go rather than an
invitation to raise the cap silently. No full demographic scale-out total is
given because eligibility and method mix are unknown; that work is deferred.

### 6.3 Throttling and pressure controls

The assembly pilot starts with Tier 3A concurrency 2, Tier 3C concurrency 8,
metadata requests 4, remote transfers 2 and 100 MiB/s transfer token bucket.
MooseFS read plus promotion traffic is capped at 120 MiB/s five-minute average
and 200 metadata operations/s. At ≥80% of either cap for ten minutes, new
launches pause and concurrency steps down. At ≥95%, running jobs finish safely
but no resume occurs without HPC Storage Reviewer approval. Reviewed waves may
raise Tier 3A to 10 and Tier 3C to 16 only after telemetry supports it.
Demographic compute starts at two jobs and transfers at one; it cannot exceed
the envelope-specific values above. Quota projection reserves 25% headroom.

## 7. Staged go/no-go transitions

Every “go” advances only to the named proposed WG state. It does not complete
the next human decision. Exact-reference, original-native-annotation,
checksum, denominator, taxon, sample and circularity gates can never be lowered
to make a quota or sample-size target.

| transition | evidence and quantitative go threshold | stop/rollback | accountable reviewer | absolute stage maximum | retention and resulting WG state |
|---|---|---|---|---|---|
| 0 → 1 guidance/schema freeze | both dependency hashes and direct edges verified; G01–G32 mapped; candidate/source/result/status vocabularies frozen; graph/TSVs parse; 100% required fields have explicit missingness | dependency/artifact absent, schema conflict, permissive resolution, or unreviewed ID change | Program Architect + Data Provenance Reviewer | 4 CPU; 8 core-h; 16 GiB; 8 h; 2 GB scratch + 2 GB persistent; zero network/download/I/O beyond local files | retain schemas/decision log 7 years; `A00` becomes askable, not approved |
| 1 → 2 inventory/preflight | after `A00`, exact object ≤1 MB; 717/716 and column assertions; 714/223/248/271/120/40, taxonomic 13/4/3/9/6/5 and 46-fish subset reproduced or a signed discrepancy no-go; 100% rows closed; proposed rows have TaxId, exact assets where needed, terms, size and modality disposition | moving/unlicensed source, count/schema ambiguity, unresolved TaxId, H2 label treated as accession, projected/mismatched annotation, unknown denominator, or silent row loss | Inventory Curator + Scientific + Provenance Reviewers | 2 CPU/element; 30 core-h; 16 GiB; 24 h; 4 GB scratch; 2 GB persistent; catalog ≤1 MB; ≤4 API requests; 20 GB reads/writes, 80k metadata ops; no biological assets | catalog/provenance retained 7 years; accepted metadata manifest freezes; phase-2 local dry run enabled |
| 2 → 3 no-cost/dry-run | existing-profile hashes match; rendered commands contain exact SweepGA/IMPG responsibility flags; 100% schema/unit/fixture and circularity negative tests pass; stale/partial/superseded fixtures fail; no network, Slurm, biological run or new realization | any ambient binary, hash/argv drift, legacy MSMC, fixture/reference/mask/circularity acceptance, or compute profile unavailable | Reproducibility + Methods Reviewers | local: 4 CPU, 8 core-h, 16 GiB, 8 h, 5 GB scratch/persistent; new shared realization separately `A20`, ≤8 CPU, 32 GiB, 8 h, 30 GB store | retain test/env packet; `A30/A31` askable only after `R20`, not automatically ready |
| 3 → 4 bounded assembly/composition pilot | `A30` exact 8 rows/assets and `A31` compute scope; ≥6 pair-eligible and all 8 composition-eligible unless reviewer redesigns before outcomes; 100% exact reference/native annotation/checksum/REF; positive denominators and base target queryability ≥50% or prespecified stratum; first-attempt technical failures ≤1/8, unexplained 0; telemetry 100%; median absolute resource error ≤35%, p95 ≤75%, no row >2× base; RSS/scratch/I/O <80% cap | any hard-gate contamination, zero denominator, mapping multiplicity >1, missing telemetry, >1/8 technical or any unexplained failure, ≥95% sustained I/O, >100% quota; stop and rollback partials; `A32` required for any overage | Scientific + Reproducibility + HPC Storage Reviewers | **small cap:** 8 slots; 8 CPU/element; 280 core-h; 96 GiB/job; 17.5 h; 150 GB scratch; 160 GB input/download + 32 GB output; 60k files; 180 GB read/200 GB write; 500k ops; 120 MiB/s; Tier3A 2/Tier3C 8/transfers 2 | sources through review; intermediate +90 days; collector enters `Q30/R30`; no expansion node ready |
| 4 → 5 independent review | closed ledger 100%; hashes/provenance concordant; call/record concordance ≥99.9%; zero circular predictors or invalid H1/H2 conversions; missingness and selection audit; pilot quality and resource thresholds above; quota headroom ≥25%; no unresolved storage incident | reviewer dissent, severe failure, resource model outside thresholds, callability/taxon ambiguity, or I/O incident | independent Scientific, Statistical/Circularity, Reproducibility and HPC Storage Reviewers | 4 CPU; 16 core-h; 32 GiB; 24 h; 10 GB; no download/Slurm | immutable review retained; only `A50` can authorize wave 1 |
| 5 → next wave | each distinct `A50/A5N` names ≤10 Tier3A and ≤32 Tier3C; 100% hard gates/ledgers; cumulative technical failure ≤5%, unexplained ≤1%; median resource error ≤25%, p95 ≤50%; no OOM/lineage mismatch/storage incident; outputs/checksums/telemetry ≥99% and ledgers 100% | two similar unexplained failures; any hard-gate breach; quota >approved or projected >80%; I/O ≥80% pause/≥95% stop; failure/resource threshold breach rolls back only unpromoted attempt | same independent panel + human project owner | each exact approval is lower; cumulative across all assembly work may never exceed high ceiling: 1,482 core-h, 96 GiB/job/960 aggregate, 19.4 machine-h, 495 GB scratch, 1,672 GB persistent, 1.67 TB read/1.02 TB write, 3.68m ops, 1.5 TB download, 200 MiB/s | wave packet retained; collector/review joins; next `A5N` remains a separate human task |
| 5 → 6 full-catalog decision | all waves reviewed; 100% candidates resolved; zero exact/native/denominator violations; retry backlog 0; independent hashes 100% and records ≥99.9%; model/tree frozen; final base/high prediction and ≥25% headroom | any blocker, dissent, high scenario beyond quota, unresolved retry, or storage/I/O saturation | joint human owner + Scientific + HPC Storage + Reproducibility Reviewers | `A60` may approve no more than the remaining portion of the same cumulative assembly high ceiling; no demographic spend | decision retained permanently; only approved Tier3A/Tier3C branches become ready; their collectors join at `J60` |
| 1/2 → 7 independent-covariate readiness | ≥95% double-entry agreement before resolution and 100% after; exact taxon/tree joins; every accepted record has source, unit, estimand, population/time, uncertainty and overlap disposition; primary overlap is `none`; simulation interval coverage ≥90% in identifiable scenarios | generic Ne, unresolved license/taxon/population, pi-conditioned inclusion, unknown/shared response overlap in primary, borrowed calibration without distribution | Source Curator + Statistical/Circularity Reviewer | bounded curation: 4 CPU, 16 GiB, 40 h plus manual work, 20 GB; no licensed/bulk export | source ledger retained; independent primary table may proceed to freeze; acquisitions still behind `A70` |
| 7 acquisition → inference | `A70` exact bytes/terms/assets and 100% checksum/reference/sample/mask gates; `A71` exact methods/species/profiles; ≥90% expected classifiable; high bound fits quota; scaled outputs first | H1/H2 shortcut, sparse VCF denominator failure, pooled/ambiguous population, missing phase evidence, invalid calibration, >10% hard-input failure, storage >110%, or circularity breach | Demographic Methods + Data Steward + HPC + Circularity Reviewers | demographic pilot ceiling: ≤6 species/method, 2,000 core-h, 256 GiB/job, 2 TiB scratch/download, 300 GB persistent, 8 TB read/3 TB write, 2m ops, 100 MiB/s, 6 compute/2 transfers | raw/intermediate per approved retention; scaled and scenario outputs/failures retained; `R71` review only, no automatic expansion |
| 7 → 8 demographic review | ≥90% eligible rows finish; truth coverage ≥90%; finite intervals in identifiable window; mask/phase/parameter/replicate stability; median resource error ≤35%, p95 ≤75% and no stage >2× base; scaled/absolute separation and calibration draws complete; zero primary leakage | optimizer/mixing failure, phase control changes conclusion, unflagged mask crash, lineage/reference/circularity breach, or one failed prespecified remediation stops method | independent Demographic, Statistical/Circularity and Reproducibility Reviewers | 4 CPU, 16 core-h, 32 GiB, 24 h, 10 GB review only | tier-labeled review packet; catalog-scale demographic expansion remains deferred unless a new graph/program is authorized |
| 8 → release | accepted+rejected+unavailable+failed equals inventories; primary independent fit frozen; tree ensemble, measurement error, missingness, multiple-testing and leave-clade-out diagnostics pass; every panel labels unit/scale/overlap/calibration; independent reproduction matches hashes/results | count/hash mismatch, unresolved tree tip, nonconvergence, circular join, hidden absolute scaling, superseded/partial input, or result changes under same hashes | Independent Results Reviewer + human project owner | analysis ≤16 CPU, 64 GiB, 48 h, 100 GB scratch/persistent, no new download; release/cleanup separately `A81` | provenance/release retained 7 years; cleanup only checksum-led; WG program closes after release review |

## 8. Universal stop conditions

The following stop the affected row; systemic occurrences quarantine the
lineage or wave. They are not “warnings”:

- no original H1-native annotation, projected/lifted annotation, sequence-set
  or annotation/reference/checksum mismatch;
- unresolved H1/H2 relationship, taxonomic ambiguity, population/sample
  ambiguity, or reuse terms;
- SweepGA query or target multiplicity above one, wrong whole-FASTA inputs,
  wrong scaffold-jump, or IMPG responsibility bypass;
- zero or below-prespecified callable/queryable denominator, low callability,
  REF mismatch, sparse-VCF invariant-site ambiguity, invalid demographic
  input, missing phase/switch evidence, or invalid optimizer/convergence;
- response-derived, shared, or unknown-overlap Ne offered as an independent
  predictor; generic Ne; or absolute Ne/time without `mu` and generation-time
  scenario provenance;
- median/p95 resource prediction error above the stage threshold, technical or
  hard-input failure above its threshold, telemetry incompleteness, <25% quota
  headroom, >100% authorized spend, or storage >110% in the demographic gate;
- MooseFS/metadata pressure ≥80% pauses launches and ≥95% stops resumption;
  any metadata-server incident stops the wave;
- any wildcard/legacy/unlisted/partial/superseded artifact selected as input.

Low callability is handled without inventing a universal cross-method number.
Assembly pilot base queryability is ≥50%, with preregistered 25/50/75%
sensitivities. Population/demographic minimum callable spans are frozen by
method-specific power simulations before biological outcomes are inspected.
Positive but inadequate denominators are ineligible, not zero-valued results.

## 9. Statistical synthesis and release outputs

The canonical long table has one row per candidate/population, sampling unit,
target class and estimand. It contains identity/taxonomy/tree, population and
dates, exact H1/H2/reference/annotation, sample/haplotype counts, original/
queryable/callable genes and bases, numerator/denominator, estimate and
uncertainty, mapping/phase/mask QC, missingness/status, predictor independence,
calibration scenario, input/output/environment/command hashes, attempt and
telemetry/review IDs.

The primary population model is a Bayesian phylogenetic errors-in-variables
model on log pi using only Tier-1 independent predictors, with species random
effects and tree covariance. Range, density and their product-derived
abundance do not enter one main model without a prespecified collinearity/
derived-variable choice. PGLS is the transparent baseline. Tree topology and
branch uncertainty, source measurement error, multiple populations, and
within-species covariance are propagated. Prior predictive simulation,
posterior predictive checks, effective sample size/convergence, influential
taxa, leave-one-major-clade-out, alternative covariance, clade balance,
complete-case and missingness/quality/transfer sensitivities are reported.

Assembly PGLS remains a separate response stratum. Independent LD/temporal/
pedigree quantities enter separate secondary models. Genomic curve features
enter only a table/model labeled `shared-genomic sensitivity`; they cannot be
described as independent evidence. The primary coefficient family is small and
prespecified; secondary tests report raw and Benjamini–Hochberg-adjusted
values. No causal vocabulary is allowed without an external identification
design.

Release contains the frozen inventories, accepted/rejected/unavailable/failed
ledgers, source and method schemas, tree ensemble/crosswalk, independent
predictor table, assembly and population responses, native scaled histories,
calibration draws and scenario histories, model code/specifications,
diagnostics, resource/I/O telemetry, Guix derivations/closures, artifact
checksums, claim-permission table and checksum-led retention/cleanup decision.

## 10. Proposed WG graph and structural assertion

The ordered, machine-readable proposal is
`analysis/vertebrate_scaleout_execution_graph.tsv`. Parallel nodes have
file-disjoint scopes; every branch has a named join. Authorization nodes are
human decision tasks whose records must include approver, timestamp, exact row
and artifact fingerprints, maximum downloads/core-hours/wall/memory/scratch/
storage/I/O/files/concurrency/bandwidth, retention, expiry and rollback owner.
No shell task can auto-complete them. No executable node has an alternate path
around its authorization predecessor.

The graph is a proposal and was **not instantiated** during synthesis. No
`sbatch`, `srun`, scheduler submission, demographic command, catalog retrieval,
bulk download, Guix realization, or biological inference was run. Repository
and process review at completion must confirm no new job IDs, download
manifests, staged assets, or demographic outputs were created. Full-scale
execution remains structurally blocked by `A60`; bulk demographic acquisition
and expensive inference remain separately blocked by `A70` and `A71`.
