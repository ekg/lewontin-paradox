# Gene-conversion and GC-biased gene-conversion evidence program

**Design freeze:** 2026-07-18 UTC<br>
**WG task:** `design-gbgc-evidence`<br>
**Execution boundary:** repository design and metadata reconciliation only. This task downloaded no biological data, realized no new software environment, submitted no local or scheduler analysis, and launched no direct, population, phylogenetic, or non-allelic computation.

## Decision

Gene conversion is not one interchangeable measurement. This program freezes four separate evidence branches with distinct observations, denominators, estimands, nulls, and claim boundaries:

1. **Direct allelic evidence** observes transmissions in complete meiotic products, pedigrees, sperm, or gametes. It can estimate detectable events per callable meiosis, interval-censored tract lengths, crossover association, and—only at informative W/S mismatches—GC transmission distortion.
2. **Population evidence** compares polarized or polarization-aware W→S (`WS`) and S→W (`SW`) allele-frequency spectra in multiple individuals from the same population. It can estimate a model-dependent population parameter such as `B = 4 N_e b` under an explicitly declared convention.
3. **Historical phylogenetic evidence** assigns substitutions to branches of close-clade/outgroup alignments. It can estimate branch-specific WS/SW substitution asymmetry and a long-term gBGC-like fixation component.
4. **Non-allelic evidence** searches resolved paralogs or segmental duplications for copy homogenization beyond gene-tree and sequence-evolution nulls. It concerns exchange among copies, not allelic meiotic gBGC.

The executable choices are:

- direct pilot `D01`, the public Arabidopsis study `PRJEB4500` / `ERP003793`, with 13 complete four-product tetrads;
- phylogenetic pilot `H01`, a five-species stickleback/pipefish panel anchored by the VGP *Spinachia spinachia* assembly `GCA_048126635.1`; and
- phylogenetic pilot `H02`, a five-species falcon panel anchored by the VGP *Falco naumanni* assembly `GCA_017639655.1`.

`D02` (the 17-member Platinum pedigree) is the direct alternate. `H03` (delphinids) is the first phylogenetic alternate; `H04` (syngnathids) is second but cannot activate until a second distinct outgroup is frozen. Alternates activate only after a pre-result identity, access, callability, topology, or orthology failure. The failed primary stays in the ledger.

The current WG graph authorizes downstream execution only for the direct and phylogenetic pilots (`pilot-pedigree-gbgc` and `pilot-vgp-phylo-gbgc`). Population (`P01`, `P02`) and non-allelic (`N01`, `N02`) work is **`NOT_RUN_DESIGN_ONLY`**. No synthesis may fill those cells from VGP H1/H2 data. It must say **population gBGC not measured** and **non-allelic conversion not measured** unless separately authorized tasks are added.

The collection-level inventory, exact accessions, sample structures, sizes, annotations, terms, sources, callability requirements, and alternates are machine-readable in `analysis/gene_conversion_dataset_manifest.tsv`. Exact estimands and controls are in `analysis/gene_conversion_estimand_manifest.tsv`; permitted and forbidden language is in `analysis/gene_conversion_claim_matrix.tsv`.

## Common vocabulary and the hard identifiability boundary

`W` means a weak base, A or T. `S` means a strong base, G or C. After strand-complement normalization, `WS` means that the ancestral or donor state is W and the derived or resolved state is S; `SW` is the reverse. `WW` and `SS` changes are GC-conservative controls. W/S labels do not by themselves supply direction: direction comes from parental transmission in the direct branch, ancestral-state inference in the population/historical branches, or a rooted copy tree in the non-allelic branch.

### H1/H2 assemblies are state observations, not conversion observations

A same-individual H1/H2 assembly pair can identify callable heterozygous states. Inside the frozen 1:1 mask it may report:

- unpolarized W/S heterozygous state counts;
- local density and reason-coded `candidate_cluster` objects; and
- context, annotation, and recombination-proxy correlations labeled descriptive.

It cannot, without additional observations, identify which allele is ancestral, which homolog donated sequence, which allele was transmitted, whether a difference resulted from mutation or conversion, or whether S was preferentially transmitted. H1/H2 pairs have neither a multi-individual frequency spectrum nor an outgroup-polarized substitution history. Similarity between duplicated H1/H2 regions also cannot establish non-allelic conversion unless the copies themselves are resolved and assigned.

Therefore the following H1/H2-only phrases are forbidden: “conversion event,” “conversion direction,” “conversion tract rate,” “parent-of-origin,” “transmission distortion,” “population B,” “historical substitution bias,” “paralog gene conversion,” and “allelic gBGC.” The allowed terms are `state_count` and `candidate_cluster`.

### Cross-branch quantities must never be numerically substituted

| Evidence object | Identifiable quantity | It is not |
| --- | --- | --- |
| direct meiosis | detectable event rate, tract bounds, CO association, S/W resolution distortion | population B or historical substitution count |
| population SFS | WS/SW frequency asymmetry and model-dependent B | event count, tract rate, or transmission observation |
| phylogenetic alignment | branch substitution asymmetry and long-term gBGC-like model parameter | direct event count, current B, or current transmission distortion |
| paralog panel | copy homogenization/candidate non-allelic tracts under a copy-aware null | allelic meiotic conversion or allelic gBGC |

Historical substitution bias is not a count of direct events. Population B is model-dependent and is not a direct tract rate. Paralog homogenization is not allelic gBGC. A human or any other single-species direct estimate may enter a comparative model only as a named, bounded sensitivity prior with a transfer function; it is not a measured rate for vertebrates.

## Dataset selection and provenance

The direct choice is based on observation structure, not taxonomic proximity to VGP. The Wijnker et al. resource resequenced all four products from 13 Arabidopsis meioses, plus parents and doubled haploids; five tetrads were sequenced at a mean 54× and eight at 14×. Complete tetrads reconstruct reciprocal crossovers and 3:1/1:3 segregation directly, while the paper also documents short-read and rearrangement filters that are essential for gene-conversion calling ([study record](https://www.ebi.ac.uk/ena/browser/view/PRJEB4500), [open article](https://pmc.ncbi.nlm.nih.gov/articles/PMC3865688/)). This is preferable for a small executable direct pilot to an isolated diploid H1/H2 pair.

The Platinum alternate contains a three-generation 17-member pedigree sequenced at approximately 50×, with 11 full siblings and public inheritance/call resources; Illumina lists the sample IDs and accessions `phs001224`, `PRJEB3381`, and `PRJEB3246` ([Platinum Genomes](https://www.illumina.com/platinumgenomes.html)). Its raw-access split, short-read tract resolution, cell-line artifacts, and only 26 potential transmissions make it a benchmark/sensitivity alternate, not an automatic replacement and not a vertebrate rate reference.

The population primary design uses the high-coverage 1000 Genomes resource. The expanded study comprises 3,202 genomes at targeted 30× and includes 602 trios; `PRJEB31736` covers the original 2,504 Phase-3 individuals and `PRJEB36890` supplies expanded samples ([study description](https://pmc.ncbi.nlm.nih.gov/articles/PMC9439720/)). The design selects unrelated individuals within one frozen same-population stratum—YRI first—rather than pooling continental labels. The 1,135-accession Arabidopsis collection `PRJNA273563` is the alternate, but only a local genetic/geographic group may be analyzed; global pooling would turn structure into a false frequency-spectrum signal ([BioProject](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA273563)).

The phylogenetic choices reuse exact VGP assemblies from the upstream frozen design and add close species/outgroups from the archived VGP `PRJNA489243` NCBI report snapshot. `H01` is favored for its 0.36–0.51 Gb genomes and five distinct species. `H02` supplies an independent avian clade with five distinct falcons. `H03` has well-annotated but much larger delphinid genomes. `H04` is inexpensive but remains non-executable until a second outgroup is exact-versioned. Exact base counts and annotation releases are not silently updated; they are frozen in the dataset manifest.

The non-allelic primary uses 47 phased diploid HPRC Release-1 genomes (94 haplotypes) plus `GCA_009914755.4` CHM13v2.0. The HPRC describes the release as freely public and provides phased assemblies and annotations ([release description](https://humanpangenome.org/release-timeline/)); CHM13v2.0 is 3,117,292,070 bp in 25 contigs ([assembly record](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_009914755.4/)). The HGSVC `PRJEB36100` alternate provides 64 haplotypes from 32 genomes, but older assemblies leave many highly identical duplications unresolved ([BioProject](https://www.ncbi.nlm.nih.gov/bioproject/PRJEB36100), [resource article](https://pmc.ncbi.nlm.nih.gov/articles/PMC8026704/)). That limitation becomes missingness, never evidence of no conversion.

All acquisition, if later authorized, must bind accession version, provider URL, expected bytes, provider digest when available, local SHA-256, retrieval UTC, exact reference/annotation, consent or reuse terms, and a content-addressed path. NCBI states its own molecular databases generally impose no use restriction but does not transfer submitter rights; assembly-specific publication/embargo status and VGP/G10K attribution must therefore be rechecked. Public access is not treated as a universal copyright or consent waiver.

## Branch 1: direct allelic events

### Observation and pilot unit

`D01` uses each four-product tetrad as one independent meiosis. Col and Ler are the heterozygous donor-parent haplotypes; Cvi is the receptor contribution used to recover the recombinant male product in each offspring. The workflow reconstructs the recombinant Col/Ler haplotype in all four products and calls reciprocal CO switches and 3:1/1:3 segregation only where the parental state is known.

The primary pilot includes every one of the 13 complete tetrads. It does not discard low-coverage tetrads after seeing events. Coverage stratum is frozen as a covariate and sensitivity. A claim-grade within-cross rate requires at least 50 independent complete meioses and at least 20 confirmed events; a claim-grade GC-resolution distortion requires at least 100 independent informative W/S converted mismatches across at least 50 meioses. The 13-meiosis pilot will therefore usually be descriptive for gBGC even if it can validate event reconstruction.

### Callable opportunities and event definition

For tetrad `m`, the callable universe is a reason-coded interval set where:

- both Col and Ler parents have a high-confidence biallelic informative state;
- all four tetrad products pass depth, genotype likelihood, allele-balance, and mapping gates;
- Cvi contribution is separable;
- sequence is uniquely mappable, copy-normal, non-rearranged, and outside unresolved repeats;
- parental phase is continuous; and
- segregation outside a candidate tract is consistent with one meiosis.

An event is a coherent 3:1 or 1:3 departure from 2:2 segregation, supported by flanking Mendelian markers and not explained by de novo mutation, sample exchange, copy change, alignment error, or a known Col/Ler rearrangement. CO-associated conversion overlaps the uncertainty interval of a reciprocal haplotype switch. NCO-associated conversion has no reciprocal CO switch in the declared interval. Complex and ambiguous patterns are retained in separate reason-coded classes.

Two denominators are reported and never conflated:

1. `events / callable_complete_meioses`, where a meiosis passes a preregistered genome-wide callable-fraction gate; and
2. `events / callable_informative_marker_opportunities`, stratified by marker density and context.

The first is an assay-level detectable event count, not corrected to all molecular events. A detection-corrected rate is secondary and exists only if spike-in simulation recovers known tracts across the empirical marker/callability map with calibrated coverage. It must carry the injected tract model in its name.

### Estimands

- `D_EVT`: CO and NCO event rates with exact or cluster-bootstrap intervals.
- `D_TRACT`: inner and outer tract bounds. Endpoints between markers are interval-censored. Midpoint lengths are never treated as observations. With at least 20 events, fit a Turnbull/nonparametric interval-censored distribution.
- `D_CO`: conversion enrichment near CO breakpoints compared with within-meiosis, chromosome-, marker-density-, GC-, and callability-matched opportunity.
- `D_GCBIAS`: among directly informative W/S mismatches in confirmed events, the probability that resolution favors S. The event, not each tightly linked marker, is the conservative cluster unit.

### Nulls, error controls, and multiplicity

Mendelian null tetrads are simulated through the measured depth, genotype error, parental marker spacing, structural variants, phase switches, and callable mask. Tracts and COs are then injected at known positions to estimate sensitivity and length censoring. Reciprocal W→S and S→W spike-ins with the same context and allele balance quantify asymmetric detection. A de novo singleton is not a conversion; an inherited CNV or rearrangement is not a conversion.

The primary total-rate test is single and preregistered. Holm correction covers CO/NCO, parental-direction, and W/S subgroup contrasts. Candidate loci use BH `q <= 0.05` with a BY sensitivity under dependence. Chromosome blocks or meioses—not markers—are bootstrap/permutation units.

### Direct branch stop/go

Identity, all-four-product completeness, parental genotype, and callability reconciliation are zero-tolerance gates. The pilot is `GO` for technical direct-event feasibility if all 13 tetrads are adjudicated, simulation false-discovery is at most 5% in the passing call class, injected-event sensitivity is reported by tract length/marker density, and every accepted event is auditable. It may report a rate only with its observed callable denominator and uncertainty. If informative W/S events do not reach the threshold, `D_GCBIAS` is `NOT_ESTIMABLE_POWER`, not zero bias.

## Branch 2: population frequency-spectrum gBGC — design only

### Observation and minimum sample

The unit is one unrelated diploid individual sampled from one population. The minimum is 50 after kinship, ancestry, coverage, missingness, and contamination QC; the target is the full passing YRI stratum from the exact 1000 Genomes metadata. Admixed labels are not pooled with YRI. Other populations are independent replication strata with their own demographic and polarization fits.

The primary neutral callable universe contains autosomal, biallelic, copy-two SNVs and monomorphic opportunities callable in at least 95% of retained individuals. It excludes coding/conserved sites, segmental duplications, low mappability, extreme depth, known selected regions in the primary mask, and CpG-prone contexts in the primary mutation-robust analysis. Every polymorphic count has the matching ancestral W or S opportunity denominator; a VCF-only numerator is insufficient.

### Polarization and mutation controls

Ancestral state uses two exact-versioned primate outgroups and a posterior, not a single forced allele. Sites at which outgroups disagree are excluded from the primary, modeled with context-specific mispolarization probability `epsilon`, and reintroduced only in sensitivity. The analysis is repeated folded, with high-posterior sites, with transversions only, and after swapping outgroup topology. A sign that disappears under folded or epsilon sensitivity is not a robust gBGC result.

Mutation bias is not gBGC. A context-dependent 3-mer/5-mer mutation matrix is fit jointly, using rare variants and an external de novo spectrum only as a sensitivity. WS, SW, WW, and SS invariant and polymorphic opportunities all enter. Reference allele is never equated with ancestral allele.

### Demography, structure, linkage, recombination, and selection

Demography is fit from GC-conservative neutral spectra with flexible size-change models; alternative histories are retained in the B sensitivity envelope. Kinship and ancestry PCs/local ancestry are audited before retaining a stratum. Background/linked selection masks and a broader mask sensitivity are reported. Genotype likelihood analysis is a sensitivity to hard-called genotypes.

The exact sex-averaged genetic map and its digest must be frozen before execution. Recombination quintiles are defined before seeing WS/SW effects. LD-aware composite likelihood, map-defined blocks, and at least 5-Mb block jackknife prevent SNP pseudoreplication. Physical distance to hotspots is a labeled proxy, not measured recombination where no map exists.

`P_SFS` reports normalized WS/SW derived-count bins. `P_B` reports `B = 4 N_e b` under the declared diploid convention, with the fitted demography, mutation matrix, epsilon, recombination model, and selection mask. `b` and `N_e` are not separately claimed from the SFS. At least 1,000 polymorphic WS and 1,000 SW sites are required for a continuous B fit.

Parametric null simulations set `B=0` but retain fitted mutation bias, demography, recombination, linkage, missingness, and polarization error; alternative simulations add known B. Posterior predictive checks must reproduce WW/SS and folded spectra before WS/SW inference. Likelihood-ratio calibration is by parametric bootstrap. There is one primary B per population, Holm correction across preregistered populations and interactions, and BH/BY correction for exploratory bins or annotations.

### Identifiability and forbidden population claims

Population data may identify an SFS asymmetry and model-dependent B. It cannot identify parent-of-origin, individual meiotic events, donor/recipient homologs, tract lengths, or a direct event rate. Population B must not be inferred from VGP H1/H2 states. The branch is `NOT_RUN_DESIGN_ONLY` in the current graph.

## Branch 3: historical close-clade evidence

### Panels and minimum phylogenetic structure

Each executable panel has five distinct species: at least three ingroup taxa and two taxa that provide independent ancestral-state information. `H01` uses *Spinachia*, *Pungitius*, and *Gasterosteus* with *Syngnathus acus* and *S. typhle*. `H02` uses five *Falco* species and fits all branches on a preregistered species tree, with the focal *F. naumanni* branch primary. A branch is interpreted only if its ancestral state is identifiable under the two-outgroup posterior or explicit polarization-error mixture.

The primary panel tests exactly one neutral focal-branch contrast. Across H01 and H02 these two tests receive Holm correction. Other branches and partitions are exploratory with BH `q <= 0.05` and BY sensitivity.

### Alignment, partitions, and semi-complete genomes

Whole-genome alignment is followed by reciprocal 1:1 synteny and single-copy orthology. No many-to-many or copy-ambiguous column enters the allelic historical branch. Each aligned base receives one primary exclusion reason: absent taxon, assembly gap/N, non-1:1, paralogous, repeat/low complexity, alignment uncertainty, topology uncertainty, ancestral ambiguity, annotation-phase mismatch, or other preregistered reason.

Semi-complete genomes are allowed because the estimand uses explicit callable opportunities. Missing sequence is missing; it is never reconstructed as ancestral reference and never counted as no substitution. A claimed neutral branch/partition requires at least 10,000 callable aligned bases and at least 20 WS plus 20 SW assigned substitutions. A fitted continuous historical-bias coefficient requires at least 100 directionally informative substitutions. Cluster inference requires 50 focal-branch WS substitutions and 20 Mb callable sequence.

Four non-interchangeable partitions are emitted:

- **neutral primary:** single-copy, noncoding, nonconserved, repeat-masked, non-CpG sequence with ancestral context;
- **fourfold:** codon-phase-audited fourfold-degenerate sites only where homologous CDS and reading frame are valid in every required taxon;
- **coding:** other coding sites, stratified by degeneracy and consequence, explicitly subject to selection; and
- **noncoding:** intron/intergenic/conserved/nonconserved strata, never all labeled neutral.

No native annotation is silently projected by coordinate. Where an exact GCA lacks annotation, fourfold/coding results for columns requiring that taxon are absent unless a downstream task performs and audits orthology-aware annotation transfer. This cannot veto the neutral whole-genome pilot.

### Historical estimands and nulls

For branch `b`, the transparent count summary is

```text
A_b = log[(N_WS,b / O_W,b) / (N_SW,b / O_S,b)]
```

where `O_W,b` and `O_S,b` are callable ancestral W and S opportunities in the same context/partition. `H_SUB` reports counts, opportunities, rates, and `A_b`. `H_GBGC` fits a context-dependent nonstationary substitution model in which a branch fixation-bias component changes WS relative to SW rates. It may be described as a **long-term gBGC-like substitution bias**, not population B, unless a population generative model is explicitly added. `H_CLUSTER` tests excess spatial clustering of same-branch WS substitutions; accepted rows remain candidate historical signatures.

The null is a fitted context-dependent nonstationary model with no WS/SW fixation component beyond mutation opportunity. Parametric simulations use fitted branch lengths, context composition, empirical gaps/masks, annotation errors, and multinucleotide mutation. An inhomogeneous Poisson/renewal null for clusters preserves chromosome, local opportunity, callable gaps, and local mutation rate. Cluster-scale scanning uses maxT permutation.

### Polarization, mutation, recombination, and phylogenetic correction

Ancestral-state posteriors integrate both outgroups. The primary excludes discordant outgroups; sensitivities vary epsilon, retain only high-posterior states, exclude CpG, use transversions, and perturb branch lengths/topology. Context-dependent mutation rates and WW/SS changes constrain mutation asymmetry. Neutral and selected partitions are never pooled to increase significance.

The species tree is primary. Local gene trees and a multispecies-coalescent/ILS sensitivity quantify discordance; known introgression is modeled or masked. Chromosome/synteny blocks are the resampling units. A frozen linkage map, if available, may define recombination strata. Otherwise local GC, distance to chromosome ends, or syntenic map transfer is labeled a recombination **proxy**, and association with that proxy is supporting heterogeneity—not proof of gBGC.

### Historical branch stop/go and forbidden claims

Both pilots must pass exact accession identity, at least 3+2 taxa, orthology reconciliation, callable-universe equality, stable focal-branch topology, two-outgroup polarization, null calibration, and minimum counts before interpretation. A sign reversal under reasonable topology or polarization sensitivity is `UNSTABLE`, not a discovery. Historical results cannot be reported as direct event counts, current transmission distortion, current population B, or measured conversion tracts.

## Branch 4: non-allelic conversion — design only

### Observation and copy-aware minimum

The unit is a paralog orthogroup in one haplotype, not a homologous allele pair. A tested family requires at least three distinguishable copies, unique flanking anchors, long-read/k-mer/read-depth copy validation, and five independent diploid individuals in the same copy-state class. Genome-wide comparison requires at least ten independent orthogroups. Collapsed, missing, or copy-number-ambiguous haplotypes are missing observations.

For every haplotype the workflow must:

1. inventory copy number and unique flanks;
2. distinguish homologous alleles from paralogous loci;
3. reconcile a copy gene tree with the species/haplotype tree;
4. align copies without forcing nonorthologous sequence;
5. run a GENECONV-like maximal-identical-fragment scan; and
6. validate candidates against raw long reads, k-mers, depth, assembly graph, and alternate haplotype.

`N_TRACT` reports interval-bounded candidate copy-homogenization tracts. Its null permutes polymorphic-site labels conditional on copy tree, divergence, context, gaps, and searched length; simulations additionally include recurrent mutation, unequal crossing-over, copy loss, assembly collapse, and error. `N_HOM` measures excess within-species between-copy identity or tree discordance against copy-age-, context-, constraint-, and linkage-matched nulls. Tandem and dispersed duplicates are separate strata.

Donor/recipient direction is reported only with a rooted copy tree and an informative outgroup copy. Otherwise it is unknown. Even a significant GENECONV-like tract establishes a candidate non-allelic homogenization event under the model, not meiotic allelic conversion, biased allelic transmission, or population B. MaxT corrects tract scans within an orthogroup, Holm corrects preregistered orthogroups/classes, and BH/BY applies only to exploratory scans.

Both `N01` and `N02` remain `NOT_RUN_DESIGN_ONLY`. HPRC H1/H2 allelic differences must not be repurposed as this branch without copy resolution and separate authorization.

## gBGC sensitivity is separate from demographic and selection interpretation

gBGC can raise S alleles in frequency and fixation in a way that resembles positive selection. It can also reshape diversity, site-frequency, coding-substitution, GC-content, and branch-length summaries that feed demographic or selection analyses. The program therefore never “controls for gBGC” with a single fitted number.

Every demographic/selection synthesis has two separately labeled paths:

1. the original demographic or selection model using its declared sites and masks; and
2. a gBGC sensitivity that excludes WS/SW changes, stratifies by recombination proxy, or propagates a bounded direct/population/historical bias range appropriate to that evidence level.

If the demographic conclusion changes sign, category, or interval coverage under a reasonable gBGC sensitivity, the conclusion is `GBGC_SENSITIVE`. The gBGC result is not then called selection, and the demographic result is not used to validate gBGC. Same-pair VGP PSMC remains descriptive and cannot infer population-frequency B.

## Frozen workflow handoff

### Direct downstream task

`pilot-pedigree-gbgc` may acquire only `D01` after recording an object manifest and accession/run inventory. It must emit:

- `direct/input_manifest.tsv` and exact run/sample/pedigree mapping;
- callable and exclusion BED/TSV objects for each tetrad;
- parental-marker and phase QC;
- `per_meiosis_callable.tsv`;
- `events.tsv` with 3:1/1:3, CO/NCO, donor/recipient allele, inner/outer bounds, and reason codes;
- `rate_summary.tsv`, `tract_summary.tsv`, `crossover_association.tsv`, and `gc_transmission.tsv`;
- injection/null calibration, seeds, software/environment digests, and a branch decision.

The task must not activate D02, download controlled phs001224 data, or expand sample scope without authorization. A low event count is a valid pilot outcome.

### Phylogenetic downstream task

`pilot-vgp-phylo-gbgc` may acquire H01 and H02 exact assemblies only after license/publication audit and object manifests. It must emit:

- species-tree and local-tree inputs with sources;
- per-taxon sequence/annotation dictionaries;
- orthology/paralogy and 1:1 callability ledgers;
- neutral, fourfold, coding, and noncoding alignment partitions, with absent partitions explicit;
- `branch_substitution.tsv`, `historical_bias.tsv`, `candidate_historical_clusters.tsv`;
- polarization, topology, mutation-context, null, block-bootstrap, and multiple-testing tables; and
- separate H01 and H02 decisions.

H03 activates only if a primary fails a pre-result gate. H04 cannot activate until a second outgroup is exact-versioned. No alternative is selected after viewing the focal WS/SW result.

### Synthesis state table

| Branch | Current state | Synthesis requirement |
| --- | --- | --- |
| direct | authorized downstream pilot, not measured by this design task | use only `pilot-pedigree-gbgc` result; otherwise say not yet measured |
| population | `NOT_RUN_DESIGN_ONLY` | say population B and frequency-spectrum gBGC not measured |
| historical phylogenetic | authorized downstream pilots, not measured by this design task | use only passing H01/H02 outputs and call them historical |
| non-allelic | `NOT_RUN_DESIGN_ONLY` | say paralog/segmental-duplication conversion not measured |

## Reproducibility, amendment, and validation contract

All future random seeds derive from the run-manifest digest. Every simulation records generator, parameters, replicate count, seed, and output digest. Callable plus mutually exclusive primary exclusion reasons must reconcile exactly to each declared universe; diagnostic overlaps live in a non-accounting table. Windows, markers, sites, and copies are not independent replicates; meioses, individuals, chromosome blocks, species-tree branches, or haplotypes are used as appropriate.

An amendment is append-only and contains parent-manifest SHA-256, changed field, old/new value, reason, author, UTC, authorization task, and whether any result was unblinded. Accession or annotation versions never float. A primary may be replaced only for a preregistered pre-result failure.

Repository validation is read-only and runs under the pinned GNU Guix channel in `analysis/guix/channels.scm` and manifest in `analysis/guix/manifest.scm`. It checks the four required artifacts, exact branch/state vocabulary, selected pilot/alternate counts, non-identifiability phrases, accessions, minimum sample/outgroup/copy rules, controls, and TSV integrity. It performs no network request, biological download, environment mutation, or job submission.

## Forbidden claims checklist

The following statements fail this design regardless of statistical significance:

- “VGP H1/H2 clusters are gene-conversion tracts” without parents/gametes or a separately valid historical/copy-aware observation.
- “The favored allele was GC” from an unpolarized heterozygous W/S state.
- “Population B is the direct conversion rate” or “B counts tracts.”
- “Historical WS excess counts meiotic conversion events.”
- “Paralog homogenization is allelic gBGC.”
- “The Arabidopsis or human direct rate is the vertebrate rate.”
- “No detected events means conversion is absent” without branch-specific power and callability.
- “gBGC confirms selection” or “selection confirms gBGC.”
- “Population and non-allelic evidence were measured” in the present execution graph.

Only the wording and evidence mappings in `analysis/gene_conversion_claim_matrix.tsv` may enter the cross-branch synthesis.
