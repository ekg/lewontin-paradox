# Vertebrate effective-population-size strategy

**Status:** execution-ready plan; no sequence acquisition or demographic inference authorized

**WG task:** `plan-vertebrate-ne-strategy`

**Plan date and source-retrieval date:** 2026-07-17 UTC

**Controlling local inputs:** `analysis/vertebrate_scaleout_origin_inventory.md`, `analysis/vertebrate_scaleout_quality_pass.md`, and the synchronized guidance `results/tier3/vgp_freeze_analysis.md`

## Decision

The practical vertebrate-wide analysis should **not** use genomic effective population size as an independent explanation of nucleotide diversity. Its primary design should model the existing, independently measured population-diversity response using provenance-complete ecological, life-history, census, and—where genuinely separate—contemporary demographic covariates in a phylogenetic measurement-error model. PSMC, MSMC2, SMC++, SFS, and site-pattern estimates use the same kinds of heterozygosity or allele-frequency observations that define the response. They belong in separate genomic-outcome tiers for descriptive history, method validation, temporal-shape comparison, non-causal decomposition, or prespecified shared-data sensitivity analyses.

The current repository contains **zero frozen method-valid vertebrate inputs** for PSMC, MSMC2, or SMC++. The reported VGP freeze counts—714 unique species, 223 completed assemblies, 248 annotated rows, 271 rows with a non-empty pair field, 120 completed plus annotated, and 40 completed plus annotated plus paired—are guidance-reported catalog counts because their underlying raw TSV is missing. The 40 are not demographic-ready. The repository has two accepted 20-individual population VCF tuples, but both are *Anopheles coluzzii* and therefore not vertebrates. Eight vertebrate model organisms named in the older population-data survey are candidates, not frozen inputs. These facts set the present eligible count at zero, not 40 or eight.

The recommended program therefore has three useful scales:

1. A vertebrate-wide primary comparative analysis targeting a realized 100–300 species with independently curated predictors, exact taxonomy/tree joins, and explicit missingness. The reported 714-species catalog is a denominator to audit, not a promised sample size.
2. A smaller absolute-calibration subset anchored by direct mutation-rate and generation-time evidence. Bergeron et al. measured de novo mutation rates from 151 parent–offspring trios across 68 vertebrate species (36 mammals, 18 birds, eight ray-finned fishes, and six reptiles), providing a strong candidate source whose exact overlap still must be joined and licensed rather than assumed ([Bergeron et al. 2023](https://doi.org/10.1038/s41586-023-05752-y)).
3. Explicitly eligible demographic pilots: PSMC 4–6 species, MSMC2 2–4, SMC++ 4–6, plus 2–5 independent LD/temporal/pedigree truth cases. Expansion ranges are conditional on pilot review and a separate human authorization.

Completion of this document is not permission to fetch the absent VGP freeze, stage sequence, create a Slurm array, or run an inference.

## Non-negotiable circularity and leakage policy

Let `pi_dataset_id` identify the response observations, including their biological samples, population, reference, callable mask, site class, genomic intervals, and computation lineage. Every proposed predictor record must declare `response_dataset_overlap` and `independence_tier` using `analysis/vertebrate_ne_source_schema.tsv`.

The following operation is prohibited:

```text
Ne_hat = pi_response / (4 * mu)
fit pi_response ~ Ne_hat
```

The prohibition includes algebraically equivalent transformations, a theta estimated from the same pi, and PSMC/MSMC2/SMC++/SFS outputs driven by the same heterozygous sites, SFS, individuals, masks, or ascertainment frame. A change of label, smoothing, time averaging, use of another program, or conversion with a mutation rate does not create independence. Such a regressor is a deterministic or noisy recoding of the response and cannot be interpreted as independent support for a Lewontin-paradox mechanism.

### Enforcement in the data model

- `response_dataset_overlap=none` is required for the primary design matrix. `unknown` fails closed.
- `same_species_different_population` may be an independent secondary record only when the population mismatch is substantively acceptable and prespecified; it is not silently species-wide.
- `same_population_different_samples` may be secondary after temporal/spatial alignment and non-overlapping sampling proof.
- `shared_samples_different_loci`, `shared_sites_or_samples`, and `derived_from_response` are excluded from the primary matrix. The last is `circular_excluded`; the others are `partially_shared` or `genomic_shared` and permitted only in labeled sensitivity/descriptive output.
- A source cannot be promoted merely because its point estimate was copied from literature. Curators must audit the literature's underlying samples and sites against the pi lineage.
- The frozen analysis table should be produced by a validator that anti-joins any predictor whose overlap is not `none`, rejects a generic `Ne` field, rejects mutation/generation scaling without foreign-keyed calibration records, and writes every rejection to a retained ledger.

### Enforcement in models and reports

- The primary model formula may reference only fields with `independent_primary`.
- A separate secondary model may add independent LD, temporal, pedigree, or literature Ne with measurement error and population/date matching.
- A shared-genomic sensitivity model is fit from a separate table and named accordingly. Its coefficients are not compared as if they were primary causal effects, and its table/figure caption states the shared-data path.
- PSMC/MSMC2/SMC++ curves may be reduced to prespecified shape features—such as standardized timing of troughs, number of robust change points, or curve clusters—in **native scaled coordinates** for descriptive comparison. They cannot be called independent drivers of pi.
- Species inclusion, method tuning, calibration choice, and time-window selection must not be chosen after looking for agreement with pi. All are frozen before the response-model join.

## What each available data form can and cannot support

| available object | what is inferable now | additional proof needed | prohibited shortcut |
|---|---|---|---|
| One deposited phased H1/H2 assembly pair | Assembly/composition properties; pairwise alignment and an assembly-defined queryable denominator when the existing Tier 3A rules pass | Versioned reciprocal H1/H2 accessions; same BioSample/individual/isolate and homolog relationship; identical-reference coordinate projection; evidence retained heterozygosity; collapse/duplication and phase QC; method-valid diploid VCF/consensus; reference-matched callable mask | Calling two FASTAs a callable diploid genotype; treating assembly differences as all biological heterozygotes; feeding H1/H2 labels directly to PSMC; calling two haplotypes a population |
| Raw reads without a validated genotype | Coverage, library, contamination, mapping and prospective power/QC after metadata preflight; no demographic history from the read files themselves | One biological individual/sample per library; ploidy and sex; exact reference; validated mapping/calling; genotype likelihood/error and depth audit; callable mask; population/date assignment for multi-sample designs | Counting read mismatches as heterozygosity; treating nominal WGS run count as an eligible individual count; pooling runs/populations before identity checks |
| H1/H2 plus raw reads from the same individual | Potential PSMC input after exact-reference mapping/calling and mask construction; potential phase validation | Library/individual identity, sufficient coverage, contamination/sex/ploidy checks, genotype error, depth/mappability masks, downsampled heterozygote-FNR calibration | Assuming assembly callability or invariant sites from alignment alone |
| One callable diploid genotype/consensus | PSMC scaled pairwise history; within-individual two-haplotype MSMC2 run only if phase/input contract is met | Exact reference and mask, genotype QC, natural diploidy, no sample mixture, callable span threshold, stable bootstrap/parameter grid | Population SFS, contemporary Ne, SMC++ population history, or cross-population history from one individual |
| Multiple accurately phased diploid genomes | MSMC2 within/cross-coalescence histories; replicated PSMC; possibly other haplotype methods | Same-population definitions, sample-specific masks, common reference, long-range phase and switch-error audit, relatives/admixture QC, balanced pair indices | Treating independent assemblies or phase blocks as homologs; substituting legacy MSMC for MSMC2 |
| Multi-sample population VCF/BCF | SMC++ and SFS only after cohort, invariant/callable denominator, reference, and mask gates; LD or temporal Ne if their distinct sampling contracts pass | Same-population sample assignments and dates, unrelatedness, all-sites/gVCF or exact cohort mask, exact reference, genotype/filter provenance; phasing/polarization only for the selected mode | Treating absent sparse-VCF sites as invariant; treating pooled geography as one population; claiming unfolded SFS without ancestral validation |
| Literature LD/temporal/pedigree estimate | A potentially independent contemporary predictor at its original population and time | Exact estimand vocabulary, original source, sample/locus counts, geography/dates, interval, method assumptions, population and response-overlap audit | Copying it into a species-level generic `Ne` column or mixing N_b, N_eLD, N_eV, inbreeding Ne, and long-term Ne |

For an H1/H2 assembly pair, all five gates are conjunctive: (1) the assemblies are homologs from the same biological individual, (2) one exact versioned coordinate reference is named, (3) heterozygous state was retained and can be reconstructed without treating assembly error or structural difference as genotype, (4) callable and missing sites can be constructed against that reference, and (5) the generated file satisfies the selected method's format and biological assumptions. Failure of any gate leaves the row in assembly/composition status. The synchronized inventory shows that the 13 reported paired fish rows lack H2 GCA versions and exact RefSeq GCF values in the guidance table; none can cross these gates from the committed text alone.

## Method contracts and fit to vertebrate data

The complete comparison is the machine-readable `analysis/vertebrate_ne_method_matrix.tsv`. The distinctions below control implementation.

### PSMC

PSMC analyzes one callable diploid genome from one biological individual. It needs the positions of heterozygous and callable/noncallable sequence, not two generic assemblies and not a population VCF. The official implementation converts a diploid consensus to PSMCFA bins, recommends depth filtering, reports parameters scaled to `2N0`, and documents both absolute scaling and the safer mutation-scaled alternatives `d_k` and `theta_k` ([PSMC software documentation](https://github.com/lh3/psmc); [Li and Durbin 2011](https://doi.org/10.1038/nature10231)). The documentation also identifies loss of heterozygotes at low coverage as a common pitfall. Thus the pilot requires at least 15× median autosomal depth for newly called inputs as an initial threshold, plus a downsampling-derived heterozygote false-negative check; this threshold is a pilot rule, not a universal biological constant.

PSMC's estimand is a piecewise history of pairwise coalescence for the two homologs under its SMC model. Population structure, introgression, linked selection, reference bias, and genotype error can look like population-size change. Recent time resolution is limited because a single diploid genome contains few very recent coalescences. PSMC is therefore a compact ancient/intermediate-history descriptor, not contemporary Ne and not an independent predictor of pi.

### MSMC2, not legacy MSMC

MSMC2 estimates within- and cross-population coalescence-rate histories from ordered haplotypes. Its multihetsep input contains phased haplotypes and distances through callable sequence; sample-specific masks and a common mappability mask must be on the same reference. The official MSMC2 options explicitly define haplotype pair indices for multiple phased diploid genomes and distinct cross-population pairs ([MSMC2 documentation](https://github.com/stschiff/msmc2)). MSMC2 must be pinned and invoked as MSMC2; the legacy MSMC executable is not an acceptable silent substitute. The associated MSMC2 work describes coalescence-rate-based population separation analyses ([Wang et al. 2020](https://doi.org/10.1371/journal.pgen.1008552)).

Cross-individual analyses require accurate long-range phase and measured switch error. Statistical phase, phase blocks from independent assemblies, and assembly labels are not interchangeable. Switch errors erase/introduce haplotype transitions and distort especially recent inference; physically phased human analyses have shown that statistical-phase artifacts can change inferred separation patterns. Sample count improves recent information but increases pair design and computation, while very recent and very deep intervals remain weak. A relative cross-coalescence curve is a summary of within/cross rates under the sampling model, not an identified instant of reproductive isolation.

### SMC++

SMC++ combines linkage along distinguished lineages with conditional SFS information from many unphased genomes. The core population-history mode accepts unphased diploid VCF data; phasing becomes relevant when distinguished alleles are taken across individuals. The official converter takes an indexed biallelic diploid VCF, contig, named populations and samples, optional distinguished pair, and a BED mask. Critically, the converter ignores VCF FILTER/QUAL, assumes the reference is ancestral during conversion, and treats sites absent from the VCF as homozygous ancestral unless they are masked. This makes an exact cohort mask and upstream filtering mandatory, not optional ([SMC++ documentation](https://github.com/popgenmethods/smcpp); [Terhorst, Kamm, and Song 2017](https://doi.org/10.1038/ng.3748)).

Polarization is mode-specific. A folded/uncertain analysis can use the documented default polarization-error treatment; `--unfold` is allowed only with validated ancestral alleles and an outgroup/reference audit. Phasing is likewise not asserted globally. Split mode estimates the selected clean two-population split model after marginal histories; it does not identify arbitrary isolation-with-migration. Long ROH from true inbreeding, selective sweeps, or unmarked missing data are confounded, and regularization/timepoint sensitivity must be reported.

### Contemporary LD, temporal, and pedigree methods

These quantities are valuable precisely because their sampling units and time windows differ from long-term coalescent histories.

- LD estimators use multilocus association in one contemporaneous population sample and estimate `N_eLD` or, for a cohort, often `N_b`. Small samples relative to true Ne yield imprecise, infinite, or one-sided estimates; physical linkage, relatives, population mixture, migration, nonrandom sampling, and overlapping generations matter. NeEstimator v2 implements LD and temporal estimators, but software availability does not make their contracts interchangeable ([Do et al. 2014](https://doi.org/10.1111/1755-0998.12157); [official NeEstimator page](https://www.molecularfisherieslaboratory.com/neestimator-software/)).
- Temporal methods estimate variance effective size from allele-frequency change between dated samples. They require at least two timepoints, consistent loci/assays, a stated sampling plan, and a generation/age-structure model. Bias from low-frequency alleles and finite samples is method-specific ([Jorde and Ryman 2007](https://doi.org/10.1534/genetics.107.075481)); overlapping generations and cohort sampling require explicit correction ([Waples and Yokota 2007](https://doi.org/10.1534/genetics.106.065300)).
- Pedigree/sibship methods infer inbreeding/variance Ne or `N_b` from pedigrees, reproductive success, parentage, or sibship frequencies. They require family/cohort coverage, mating-system and sex information, genotyping-error models, and missing-parent assumptions. A single-sample sibship estimator is described by [Wang (2009)](https://doi.org/10.1111/j.1365-294X.2009.04175.x). Captive, hatchery, or managed pedigree estimates are not wild-species values without a scenario-labeled transfer.

An LD, temporal, or pedigree estimate can enter the secondary predictor tier only when its samples and data are independent of the pi response and its population/time match is recorded. It never loses its specific estimand code.

### SFS and site-pattern/multispecies-coalescent methods

Folded or unfolded SFS inference (dadi/moments class) uses multiple population genotypes and a prespecified parametric model to estimate model-specific theta, relative sizes, divergence times, and migration parameters. It needs the invariant/callable denominator and an ascertainment model. Phasing is not required; ancestral polarization is required only for an unfolded SFS. Model nonidentifiability, local optima, linked selection, population structure, and mutation scaling dominate interpretation. The foundational diffusion and moments approaches are described by [Gutenkunst et al. (2009)](https://doi.org/10.1371/journal.pgen.1000695) and [Jouganous et al. (2017)](https://doi.org/10.1534/genetics.117.200493). These are genomic outcomes and therefore circular when computed from the response data.

Multispecies site-pattern/coalescent methods such as BPP use many independent locus alignments, multiple alleles per species, and a species/population map. They estimate branch/population `theta` and divergence `tau` under a multispecies coalescent; they are neither contemporary Ne nor a generic site-pattern shortcut. The BPP documentation states the distinct theta and tau parameters and their locus/species assumptions ([BPP manual](https://bpp.github.io/bpp-manual/bpp-4-manual/); [Yang 2015](https://doi.org/10.1093/czoolo/61.5.854)). Orthology, paralogy, recombination, introgression, heredity/ploidy scalars, clock priors, and MCMC convergence make this a small-clade validation option, not the vertebrate-wide first pass.

## Scaled quantities, absolute values, and calendar time

Native scaled results are first-class outputs, not incomplete versions to overwrite. Each result row records method, formula/convention, genomic compartment, ploidy, and the exact software output field.

| method/output | retain without borrowed calibration | permitted absolute conversion |
|---|---|---|
| PSMC | `theta0`, `t_k`, `lambda_k`, `d_k = 2 mu T_k`, `theta_k = 4 N_k mu` | For diploid autosomes and PSMC bin size `s`: `N0 = theta0/(4 mu s)`, `T_k = 2 N0 t_k` generations, `N_k = N0 lambda_k`; years `= T_k g` |
| MSMC2 | time boundaries and within/cross coalescence-rate or lambda summaries in the implementation's native scale; RCCR labeled as derived ratio | Apply the documented MSMC2 convention with autosomal diploidy and sampled `mu`; calendar years additionally sample `g`; record rho/mu or fixed-rho settings |
| SMC++ | model JSON, native/coalescent coordinates when plotted without generation scaling, likelihood/regularization settings | Absolute history is conditional on the per-generation `mu` supplied to estimation; years additionally require `g`; rho interpretation is separate and retains recombination assumptions |
| SFS | model theta, relative sizes, relative times/migration in the model's native units | Apply the exact model definition, commonly diploid-autosomal `theta=4Ne mu L` or per-site counterpart, with `mu`, sequence-length convention, and `g` |
| BPP-class | `theta=4Ne mu` for diploid autosomal loci and `tau` in expected substitutions/site, with locus/heredity scalars | Sample locus-rate/mutation and generation-time distributions; apply each locus/ploidy convention, never a single borrowed scalar across compartments |

Every absolute record must include:

1. the algebraic convention and whether Ne is diploid-autosomal, haploid, sex-linked, mitochondrial, or another heredity scalar;
2. mutation rate per site per **generation**, its population/species, assay (direct trio, phylogenetic, model-predicted), original interval/posterior, and taxonomic transfer distance;
3. generation time in years per generation, definition (mean parental age/generation interval versus database generation length), population, sex/management context, and uncertainty;
4. recombination map/rate or rho/mu treatment where the fitted method uses it, including reference/map version and any uniform-rate assumption;
5. calibration source record IDs, population and taxon, any domesticated-to-wild distinction, and a written transfer rationale;
6. uncertainty propagation through joint posterior draws or a prespecified scenario grid.

The default scenario grid should draw `mu` and `g` from their reported distributions. If covariance is known—generation time and mutation rate often are not independent—it must be retained. If only ranges exist, use interval/scenario distributions rather than midpoint point estimates. Direct same-population values receive the highest calibration quality; same-species values are next; congener/family or model-predicted values are scenario-only. No absolute point estimate may hide a borrowed rate. If calibration fails, the result remains scaled.

## Tiered data model and supported claims

| tier | contents | primary permissions | claims supported |
|---|---|---|---|
| 1: independent primary covariates | Census abundance, range/EOO/AOO, density, body mass, fecundity, longevity, generation length/time, reproductive mode, dispersal, mating system, direct mutation rate, and other measurements demonstrably independent of pi | Primary phylogenetic model after exact taxon/source/unit/uncertainty and response-overlap audit | Comparative association and predictive contribution conditional on measured covariates, tree, sampling, and missingness; not direct proof of historical Ne causation |
| 2: independent secondary demographic records | Literature or reanalyzed `Ne_LD`, `N_eV`, pedigree/inbreeding/variance Ne, or `N_b`, each population/date/method specific | Secondary measurement-error model; primary only if exact independence and alignment were preregistered | Whether independently measured contemporary breeding/drift quantities covary with pi across the eligible subset |
| 3: coalescent-scaled genomic histories | PSMC/MSMC2/SMC++ native histories, SFS theta/relative parameters, site-pattern theta/tau | Separate tables; descriptive curves/features; method validation; non-causal decomposition; shared-data sensitivity | Historical shape, clustering, method agreement, and plausible demographic context under method assumptions |
| 4: derived absolute scenario histories | Tier 3 plus explicit mutation/generation/recombination calibration draws | Scenario distributions only; never replaces scaled tier | Conditional ranges of Ne and calendar time and calibration sensitivity |
| 5: excluded/ambiguous | `pi/(4mu)` from response, shared-response predictors offered as independent, generic Ne, invalid H1/H2 conversions, unresolved taxa/populations/references, hidden borrowed calibration | Retained rejection ledger; no model input | Missingness, availability, and bias audit only |

The source table remains long-form: one row per source × estimand × population × measurement interval × calibration scenario. Contemporary, variance, inbreeding, linkage, breeding, long-term coalescent, and census quantities cannot be collapsed into one `Ne` column.

## Independent sources and curation

Candidate independent sources should prioritize original studies and stable, licensed releases. AnAge provides curated longevity and related life-history records but remains a secondary compilation whose original citations and field definitions should be captured ([de Magalhães, Costa, and Church 2009](https://doi.org/10.1111/j.1420-9101.2009.01783.x); [AnAge](https://hagr.ageing-map.org/species/)). IUCN records may contribute range, population trend, census/range definitions, assessment date, and taxon concept when access terms permit, but categories and trends are not numerical abundance and revised taxonomy can change the represented concept ([IUCN Red List API](https://api.iucnredlist.org/)). Bergeron et al.'s direct trio rates should be retained with their trio count, parental ages, population/captive status, coverage, callable sites, uncertainty, and species concept rather than copied as one unqualified rate.

For every literature/ecological record, `analysis/vertebrate_ne_source_schema.tsv` captures the original DOI/URL, version/hash/license, retrieval date, source locator, extraction and verifier, quoted/source-specific estimand definition, method, native units, interval/error, individuals/families/loci, life stage and sex, population/geography, sampling and measurement dates, taxon authority/ID, synonym path, match distance, transfer rationale, and quality. Quantitative truth cases use independent double entry; discrepancies in value, unit, interval interpretation, or population are resolved against the original source.

### Taxonomy, population, reference, and time joins

- **Taxon:** join through a versioned authority ID and taxon concept, never a normalized binomial alone. Exact concept or documented synonym is allowed in the primary tier. Subspecies-to-species aggregation must match the response's concept and becomes sensitivity unless coverage is representative. Congener/higher transfer is scenario-only.
- **Population:** retain a stable population ID, verbatim label, polygon/coordinates and geography method. An exact same population is preferred. A species-level ecological trait can join at species level; a population Ne cannot be generalized across a species without a hierarchical model and transfer variance.
- **Reference:** genomic inputs join by versioned accession plus FASTA hash and exact contig dictionary. VCF/mask/reference mismatch is a hard exclusion, not a covariate.
- **Time:** record measurement ranges and precision. The primary species-trait model may use slowly varying adult body mass/reproductive mode if definitions match; census, density, LD, temporal or pedigree measures require explicit overlap/proximity to response sampling. A preregistered score combines date separation in generations and geographic/population match. Low-score records move to sensitivity rather than being silently averaged.

### Missingness and measurement error

No missing value is zero. Missingness codes distinguish not measured, not reported, inaccessible, failed QC, structural exclusion, and unmatched taxon. The main analysis uses complete or partially observed records under a joint/hierarchical measurement model; it does not mean-impute and then treat imputed values as known. Continuous source intervals/posteriors enter as likelihoods or latent values. Counts may use lognormal/negative-binomial or interval-censored measurement models as appropriate; range and density share errors if used to derive abundance. Categorical reproductive modes retain uncertain/compound states.

The missingness analysis reports availability by clade, body size, conservation status, data modality, and pi magnitude **without using pi to decide inclusion**. Complete-case PGLS is a transparent baseline; multiple/model imputation is a labeled sensitivity with imputed predictors excluded from the main evidentiary fit. Results are repeated excluding domesticated/captive taxa, taxonomic transfers, stale measurements, and records below each quality grade.

## Primary phylogenetically aware analysis

### Response and unit

The sampling unit is a biological population where population pi exists and a species otherwise; these are not mixed without a hierarchy. The response is the frozen Tier 3 population diversity statistic and interval/replicate information, with reference, callable site class and population ID. Multiple populations per species receive a species random effect and within-species covariance. Assembly-pair diversity is a distinct response stratum, not silently pooled with population pi.

### Preregistered model family

For positive pi, fit a Bayesian phylogenetic errors-in-variables model on log scale:

```text
log(pi_j) ~ Normal(alpha
                   + beta_N log(Nc_j or abundance_j)
                   + beta_R log(range_j)
                   + beta_D log(density_j)
                   + beta_M log(body_mass_j)
                   + beta_G log(generation_time_j)
                   + beta_F fecundity_j
                   + beta_L log(longevity_j)
                   + reproductive_mode effects
                   + dispersal effects
                   + clade-balanced prespecified interactions
                   + species_random[j],
                   sigma_observation_j)

species_random ~ MVN(0, sigma_phylo^2 * C(tree, lambda) + sigma_species^2 * I)
```

Do not put range, density, and their product-derived abundance in one main model without a collinearity/derived-variable decision; use prespecified alternative formulations. Direct mutation rate is an independent predictor or offset only in a separately motivated model because mutation rate is part of the process generating pi; it is not used to manufacture Ne. Estimate across a frozen tree ensemble or propagate branch/topology uncertainty. PGLS is a frequentist baseline; phylogenetic generalized models explicitly account for shared ancestry and can include measurement error ([Martins and Hansen 1997](https://doi.org/10.1086/286013)).

Use shrinkage/regularizing priors selected by prior predictive simulation, standardized predictors, variance-inflation and posterior-correlation checks, and a limited prespecified interaction set. Evaluate calibration and effective sample size, posterior predictive residuals by clade, leave-one-clade-out prediction, influential species/populations, and sensitivity to Brownian versus justified alternative covariance. Multiple testing is controlled by defining a small primary coefficient family and treating exploratory traits as shrinkage/exploratory outputs.

### Prespecified analyses

1. **Primary:** independent Tier 1 predictors only.
2. **Secondary independent-Ne:** add each specific LD/temporal/pedigree estimand separately with its error distribution and population/date match; never a pooled generic Ne.
3. **Mutation/generation:** direct same-species trio mutation rates and independently sourced generation times; borrowed/model rates excluded or scenario-labeled.
4. **Missingness/quality:** complete-case grade A/B, hierarchical missing-data model, and grade/transfer exclusions.
5. **Phylogeny:** tree ensemble, leave-one-major-clade-out, and clade-stratified coefficients where power permits.
6. **Response construction:** population pi versus assembly-pair response strata; callable-site and neutral-site definitions; no mixing of denominators.
7. **Shared genomic sensitivity:** standardized PSMC/MSMC2/SMC++ shape features in a separately labeled, explicitly non-independent model; interpretation limited to concordance/non-causal decomposition.

The primary model can support statements such as “pi is associated with range and life history after accounting for shared ancestry and measurement error.” It cannot establish census causation, identify a single universal Ne, or validate a mechanism using an Ne algebraically/genomically derived from the same response.

## Availability targets and pilot truth cases

Counts below are planning envelopes, not claims of current eligibility.

| stream | current frozen vertebrate count | pilot target | reviewed expansion target | truth-case design |
|---|---:|---:|---:|---|
| Independent ecological/life history | 0 completed joined table | 30–50 species spanning mammals, birds, reptiles, amphibians and fishes | 100–300 realized species | 10 records double-extracted from original studies/databases; exact taxon/tree joins; known unit/date edge cases |
| Direct mutation rate | 0 joined; 68-species external candidate study | 12–20 exact-overlap species | all exact/synonym overlap with uncertainty, no taxonomic borrowing in main tier | One high-quality trio-based model species per major represented clade; reproduce source value/interval from supplement |
| Literature LD/temporal/pedigree | 0 inventoried | 10–20 species/records; 2–5 raw reanalysis truth cases | 20–60 records if independence/provenance passes | Simulated known Ne plus a vertebrate population with published genotype data and published estimate; blind agreement interval |
| PSMC | 0 | 4–6 species/individuals across genome size, coverage and clade | 10–20 species; catalog ceiling ≤40 after proof | Human or another exceptionally documented benchmark; one VGP same-individual raw-read case; downsample and mask perturbation |
| MSMC2 | 0 | 2–4 species with physical/trio phase and measured switch error | 4–10 species | Human truth set with physical phase; synthetic switch-error series; no assembly-only truth case |
| SMC++ | 0 | 4–6 model vertebrate populations with frozen VCF/masks | 8–15 species/populations | Public benchmark population, simulated bottleneck/growth, mask-induced false-ROH negative control |
| SFS | 0 | 3–5 populations after SMC++ preflight | 6–12 | Simulated identifiable model; folded/unfolded comparison; optimization replicate recovery |
| Site-pattern/BPP | 0 | 1–2 small clades | 2–4 clades only after convergence review | Simulated loci plus published small-clade dataset with known theta/tau convention |

Pilot selection is stratified before inference by clade, genome size (small/median/large), repeat content, coverage, data modality, and calibration quality. Human can be a software/input truth case but must not set universal parameters. At least one fish and one bird/reptile should test non-mammalian mutation/generation scaling; an amphibian enters only if exact input/calibration gates pass. A failure to find an eligible clade is reported as structured missingness, not repaired by a congener.

## GNU Guix, compute, storage, and I/O

The validated baseline is pinned by `analysis/guix/channels.scm` at Guix commit `44bbfc24e4bcc48d0e3343cd3d83452721af8c36` (file SHA-256 `45c055cd1d9010a72eacbb720037a22bccb2d8d6891dbd11b5d663f29b3a17`) and `analysis/guix/manifest.scm` (SHA-256 `2fb05e87aa2ac45ce51d4dcf93b232cb98627f525adace98357629ee3f15720a`), with recorded profile `/gnu/store/z9v2f6faha9cwjz0sm5iphhlzisgi077-profile`. It supplies the established Python/scientific and bcftools/samtools lineage, but it does **not** establish PSMC, MSMC2, SMC++, NeEstimator, phasing, dadi/moments, or BPP.

Each method receives a separate `analysis/guix/ne/<method>-manifest.scm` and any reviewed package definition. Do not edit the baseline to make an inference appear covered. Freeze source URL, tag and commit, source hash, patches, compiler/runtime, package derivation, closure hash/size, profile store path, tool-reported version, wrapper and commands. For SMC++ the documented `latest` container is not reproducible provenance; a specific source revision or versioned image ingredients must be packaged and verified. For MSMC2, pin the D compiler/GSL and helper-tool commit. For licensed tools, record license approval before packaging.

The cluster execution pattern, if later authorized, uses the existing `analysis/slurm/guix_job.sh`, a GC-rooted profile realized on `octopus01`, verified shared store paths on compute nodes, immutable read-only inputs on MooseFS, node-local scratch, bounded arrays, atomic output promotion, idempotent fingerprints, checksums, `afterok` dependencies, and `sacct` plus explicit application-I/O telemetry. No tool builds or large streaming reads run from compute jobs.

### Quantitative planning envelope

| unit | CPU / memory / wall | scratch / persistent | I/O notes |
|---|---|---|---|
| Metadata/source curation batch | 1–4 CPU / ≤16 GiB / ≤8 h | ≤10 GiB / ≤10 GiB | API/export snapshots only after terms review; manual curation dominates |
| PSMC from accepted callset, 100 bootstraps | 2–8 CPU / 8–16 GiB / 1–6 h | 5–30 GiB / <2 GiB | One sequential reference/consensus read plus bootstrap segments; raw calling is separate |
| Raw-read mapping/calling per individual | 16 CPU / 32–64 GiB / 12–48 h | 50–300 GiB / 20–100 GiB retained callset+mask | Approx. 2–4× compressed-read bytes in scratch and 3–6× read/write amplification; measured pilot replaces bound |
| MSMC2 prep/inference per comparison | 8–32 CPU / 32–128 GiB prep, 16–64 GiB infer / 12–72 h prep, 6–48 h infer | 10–100 GiB / 1–10 GiB | Phasing and VCF scans dominate; bootstrap multiplies inference reads |
| SMC++ per population/model | 8–32 CPU / 64–256 GiB / 12–96 h | 20–100 GiB plus 0.5–2× compressed VCF / 2–20 GiB | Converter/inference scale with analyzed length and composite count; stage once locally |
| LD/temporal/pedigree per population | 1–64 CPU / 8–128 GiB / 0.5–72 h | 5–100 GiB / <10 GiB | Locus pairs or family assignments may be superlinear; pruning/power pilot required |
| SFS per model | 8–32 CPU / 16–128 GiB / 2–72 h × 50–500 starts/bootstraps | 20–200 GiB / 1–20 GiB | Cache SFS; do not reread whole VCF per optimizer start |
| BPP small clade | 8–64 CPU / 32–256 GiB / 1 day–weeks | 50–500 GiB / 10–100 GiB | MCMC chains checkpoint locally then atomically promote; no scale-out estimate until pilot |

The proposed demographic pilot is capped at six species per method, **≤2,000 aggregate core-hours, ≤256 GiB per job, ≤2 TiB scratch, and ≤300 GiB persistent output** before a human changes those numbers. Raw reads may exceed the storage cap; such candidates are excluded or require a revised authorization. Resource accounting includes provider download bytes, compressed and decompressed bytes, local-scratch amplification, MooseFS reads/writes, peak bandwidth, inode count, and metadata operations. Low/base/high estimates come from provider sizes and pilot telemetry, not genome length alone.

## Phased acquisition and inference plan

The inert proposed graph is `analysis/vertebrate_ne_execution_graph.tsv`. It deliberately contains a non-ready `authorize-vertebrate-ne-execution` human-decision task between metadata/preflight and sequence staging, and a second authorization before any expansion wave.

### Phase 0 — freeze contracts and audit metadata

Freeze the estimand vocabulary, response IDs, overlap policy, taxonomy/reference/population joins, source schema, and method matrix. Acquire only approved metadata and bibliographic records. If the VGP raw freeze is later approved for acquisition, content-lock it, reproduce its row counts, resolve versioned H2/GCF fields, and retain every rejected row. No FASTA/VCF/read download occurs.

**Go:** 100% candidate rows have a method-specific disposition; every accepted source has DOI/URL, source locator, taxon ID/concept, unit, estimand and uncertainty status; reported freeze counts reproduce or discrepancies are resolved.

**No-go:** moving/unlicensed source, absent raw catalog, unresolved schema/count contradiction, or generic pair/Ne flags used as eligibility.

### Phase 1 — independent covariates, taxonomy/tree, and environments

Curate the Tier 1/2 source table; double-enter truth records; freeze the tree ensemble and population/time joins. Build separate Guix profiles and run no-data version/unit tests. Construct tiny synthetic PSMCFA, multihetsep, SMC++, SFS, LD/temporal/pedigree, and BPP inputs with known failure cases.

**Go:** ≥95% double-entry agreement before resolution and 100% afterward; all exact-taxon/tree joins pass; all packages reproduce identical hashes/results on login and compute smoke; circular records fail closed.

**No-go:** a tool is available only as unpinned `latest`, a legacy method is substituted, license is unresolved, or a synthetic reference/mask/overlap violation is accepted.

### Phase 2 — candidate preflight and power/resource model

Use provider metadata/HEAD requests only to resolve sample counts, individual relationships, reference versions, masks, phase evidence, source byte sizes, licenses, and calibration sources. Simulate each design over expected Ne, coverage, phase/genotype error, callable span, marker count, and demographic misspecification. Produce an exact pilot manifest and low/base/high I/O/compute budget.

**Go:** every proposed row passes all hard input gates; predicted interval coverage is ≥90% in prespecified identifiable scenarios; no method has >20% expected unclassifiable rows; base resources fit the cap and high scenario fits available quota/retention.

**No-go:** H1/H2 is the only genotype evidence, sparse VCF lacks callable denominator, MSMC2 lacks phase error evidence, population assignments are pooled/ambiguous, or borrowed calibration lacks a distribution.

### Explicit execution authorization

A human must approve the exact species/sample/accession manifest, methods, URLs/access terms, provider and scratch byte caps, core-hour/memory caps, retention/deletion policy, Guix derivations, and rollback owner. The task is proposed, not created or ready. No response by a deadline means **no**. Plan completion and environment success are not implicit approval.

### Phase 3 — staged truth/pilot, only after authorization

Stage only authorized objects, verify provider and repository hashes, construct exact-reference inputs, and run stratified truth cases before biological rows. Run scaled inference first. Generate absolute histories only from joined mutation/generation scenario draws. Capture application and scheduler telemetry.

**Go:** 100% reference/sample/mask/hash hard gates pass; ≥90% of eligible pilot rows finish; simulations/benchmark summaries cover truth in ≥90% of prespecified scenarios; bootstrap/scenario intervals are finite in the identifiable window; replicate-individual curve discrepancy is within a prespecified simulation-derived tolerance; median wall/scratch/I/O is ≤2× base and no high bound/quota is exceeded.

**No-go:** any invalid assembly shortcut; any circular result enters the primary table; mask perturbation produces an unflagged recent crash; phase-error control changes MSMC2 qualitative conclusion; optimizer/mixing failure; >10% hard-input failures; runtime >2× base without reviewed model; storage >110% authorized cap.

### Phase 4 — independent review and optional waves

Review science, exclusions, circularity, calibration, resources, and missingness by method. A second human task authorizes an exact next wave. Arrays are bounded and clade/genome-size stratified; failed elements remain in the denominator.

**Stop the method** if truth recovery, mask/phase robustness, or identifiability fails after one prespecified remediation; do not tune until agreement with pi. **Stop the wave** at any provenance/reference/circularity breach, >10% hard-input failure, or quota breach. A method can be retained for a smaller descriptive subset even when vertebrate-wide expansion is no-go.

### Phase 5 — primary fit and synthesis

The independent primary model does not wait for demographic inference. It proceeds after source/tree/model freeze and data review. Demographic results join only later as Tier 3/4 descriptive panels or the prespecified shared-data sensitivity. Release includes accepted and excluded records, native scaled histories, scenario draws, software/input hashes, telemetry, and a claim-permission table.

## Validation thresholds and failure ledger

Every candidate receives one of `eligible`, `ineligible`, `unavailable`, `failed_preflight`, `failed_inference`, or `excluded_circular`; no failure is dropped from denominators. Hard thresholds at initial pilot freeze are:

- exact reference dictionary and allele audit: 100% match;
- input/provider/repository checksums and sample-list hashes: 100% present and verified;
- PSMC novel-call coverage: median autosomal ≥15×, with downsampling FNR and depth-tail audit;
- novel multi-sample calling: target ≥10–15× and ≥90% individuals callable at retained sites, with the exact threshold simulation-validated;
- genotype/mask callable bases: positive and above a method-specific simulation-derived minimum; no universal genome fraction is invented;
- MSMC2: measured switch-error and phase-block distribution; cross-individual run fails if phase evidence is missing;
- SMC++: all absent sparse-VCF sequence distinguishable from callable invariant sequence; upstream FILTER/QUAL application proven;
- replicate/parameter/bootstrap stability: robust interval and shape over the prespecified grid, judged only within the simulated identifiable time window;
- optimization/MCMC: multiple starts/chains, convergence diagnostics and effective sample size thresholds frozen per tool in Phase 1;
- calibration: `mu` and `g` have source record, taxonomic distance, interval/distribution and scenario draws; otherwise scaled-only;
- circularity validator: zero primary rows with overlap other than `none`;
- curation: 100% quantitative accepted records verified and 100% exclusions retain reason.

Thresholds based on coverage, callable span, switch error, and curve agreement are revised only from blinded truth simulations and recorded before biological-response comparison. Resource thresholds are already quantitative and cannot be loosened inside a running array.

## Claims and decision table

| decision | proceed now | requires explicit authorization | deferred/no-go now |
|---|---|---|---|
| Freeze schemas, curate primary literature metadata, build taxonomy/estimand crosswalk, design simulations | Yes, within approved local/metadata scope | Licensed database exports or missing VGP catalog retrieval may need source-specific approval | None |
| Build method-specific Guix definitions and tiny synthetic tests | After source/license and quota review; no biological inputs | Compute-node realization if it changes shared store/quota | Ambient, unpinned or `latest` software |
| Fit primary independent-covariate model | After predictor/response/tree/model freeze and analysis review | Normal analysis execution review | Any genomic-shared predictor in primary formula |
| Stage FASTA/VCF/raw reads or run PSMC/MSMC2/SMC++/SFS/BPP | No | `authorize-vertebrate-ne-execution` exact-manifest human decision | Using current H1/H2 catalog as input |
| Expand beyond pilot | No | Separate `authorize-vertebrate-ne-wave` after pilot review | Automatic wave/full-catalog dispatch |
| Report absolute Ne/time | Only for already available, source-complete scenario calculations | New inference follows execution authorization | Borrowed-rate point estimates without uncertainty |

The strongest defensible paper result is a phylogenetically aware association between pi and independent ecological/life-history/demographic evidence, accompanied by a transparent data-availability and measurement-error analysis. Genomic histories can enrich the biological narrative and reveal temporal heterogeneity, but because they reuse diversity-bearing observations they are supporting descriptions, not an independent escape from Lewontin's paradox.

## Authoritative sources consulted

All web sources below were retrieved 2026-07-17 UTC. Local controlling artifacts were read from their committed versions.

- Li H, Durbin R. 2011. “Inference of human population history from individual whole-genome sequences.” *Nature* 475:493–496. [doi:10.1038/nature10231](https://doi.org/10.1038/nature10231). Official implementation and scaling/input documentation: [lh3/psmc](https://github.com/lh3/psmc).
- Schiffels S, Wang K. 2020. “MSMC and MSMC2: The Multiple Sequentially Markovian Coalescent.” *Methods in Molecular Biology* 2090:147–166. [doi:10.1007/978-1-0716-0199-0_7](https://doi.org/10.1007/978-1-0716-0199-0_7). Official implementation: [stschiff/msmc2](https://github.com/stschiff/msmc2).
- Wang K, Mathieson I, O'Connell J, Schiffels S. 2020. “Tracking human population structure through time from whole genome sequences.” *PLOS Genetics* 16:e1008552. [doi:10.1371/journal.pgen.1008552](https://doi.org/10.1371/journal.pgen.1008552).
- Terhorst J, Kamm JA, Song YS. 2017. “Robust and scalable inference of population history from hundreds of unphased whole genomes.” *Nature Genetics* 49:303–309. [doi:10.1038/ng.3748](https://doi.org/10.1038/ng.3748). Official input/mask/polarization documentation: [popgenmethods/smcpp](https://github.com/popgenmethods/smcpp).
- Do C et al. 2014. “NeEstimator v2: re-implementation of software for the estimation of contemporary effective population size from genetic data.” *Molecular Ecology Resources* 14:209–214. [doi:10.1111/1755-0998.12157](https://doi.org/10.1111/1755-0998.12157). [Official software page](https://www.molecularfisherieslaboratory.com/neestimator-software/).
- Jorde PE, Ryman N. 2007. “Unbiased estimator for genetic drift and effective population size.” *Genetics* 177:927–935. [doi:10.1534/genetics.107.075481](https://doi.org/10.1534/genetics.107.075481).
- Waples RS, Yokota M. 2007. “Temporal estimates of effective population size in species with overlapping generations.” *Genetics* 175:219–233. [doi:10.1534/genetics.106.065300](https://doi.org/10.1534/genetics.106.065300).
- Wang J. 2009. “A new method for estimating effective population sizes from a single sample of multilocus genotypes.” *Molecular Ecology* 18:2148–2164. [doi:10.1111/j.1365-294X.2009.04175.x](https://doi.org/10.1111/j.1365-294X.2009.04175.x).
- Gutenkunst RN et al. 2009. “Inferring the joint demographic history of multiple populations from multidimensional SNP frequency data.” *PLOS Genetics* 5:e1000695. [doi:10.1371/journal.pgen.1000695](https://doi.org/10.1371/journal.pgen.1000695).
- Jouganous J et al. 2017. “Inferring the joint demographic history of multiple populations: beyond the diffusion approximation.” *Genetics* 206:1549–1567. [doi:10.1534/genetics.117.200493](https://doi.org/10.1534/genetics.117.200493).
- Yang Z. 2015. “The BPP program for species tree estimation and species delimitation.” *Current Zoology* 61:854–865. [doi:10.1093/czoolo/61.5.854](https://doi.org/10.1093/czoolo/61.5.854). Official method documentation: [BPP manual](https://bpp.github.io/bpp-manual/bpp-4-manual/).
- Bergeron LA et al. 2023. “Evolution of the germline mutation rate across vertebrates.” *Nature* 615:285–291. [doi:10.1038/s41586-023-05752-y](https://doi.org/10.1038/s41586-023-05752-y).
- de Magalhães JP, Costa J, Church GM. 2009. “A database of vertebrate longevity records and their relation to other life-history traits.” *Journal of Evolutionary Biology* 22:1770–1774. [doi:10.1111/j.1420-9101.2009.01783.x](https://doi.org/10.1111/j.1420-9101.2009.01783.x). [AnAge database](https://hagr.ageing-map.org/species/).
- Martins EP, Hansen TF. 1997. “Phylogenies and the comparative method: a general approach to incorporating phylogenetic information into the analysis of interspecific data.” *American Naturalist* 149:646–667. [doi:10.1086/286013](https://doi.org/10.1086/286013).
- GNU Guix. Reproducible profiles require both a manifest and a channel specification; the project records both and their realized profile. [GNU Guix cookbook, reproducible profiles](https://guix.gnu.org/cookbook/en/html_node/Reproducible-profiles.html).
