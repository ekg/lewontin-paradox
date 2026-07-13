# Tier 3 frozen decisions

Status: **frozen for implementation** (`tier3-decisions-v1`), 2026-07-13.

This document supersedes the open choices and the Conda/Apptainer examples in
`TIER3_EXECUTION.md`.  “MUST”, “MUST NOT”, “REQUIRED”, and “INELIGIBLE” are
normative.  A later change to a frozen value requires a new decision version,
an explanation in the manifest, and revalidation of affected results.  It
must never be made silently in a run script.

The governing principle is simple: a diversity numerator is not a diversity
estimate until its invariant denominator is identified on the same reference.
No VCF row, no matter how reputable its source, is allowed to manufacture
invariant sites by assuming that every absent VCF record is callable reference.

## 1. Answers to every open decision in `TIER3_EXECUTION.md` section 11

| open decision | frozen answer |
|---|---|
| Cluster and scheduler | Use the institutional UTHSC cluster named `linux`, with SLURM, `octopus01.uthsc.edu` as the controller/head node and the `workers` partition for ordinary Tier 3 jobs. Scheduler-specific scripts target SLURM. |
| Compute-node internet | Direct HTTPS egress is allowed, and was demonstrated from `octopus07` to NCBI in SLURM job `1753211`. Every remote object is nevertheless staged to `/moosefs/erikg/tier3data` before analysis and verified by SHA-256. If a provider blocks compute-node traffic, prestage on `octopus01`; streaming an unpinned remote object is forbidden. |
| Population downsample | Select exactly 20 unrelated diploid individuals (40 autosomal chromosome copies), or 20 independently derived inbred lines treated as one haploid consensus each, before computing the cohort mask. Selection is deterministic as specified in section 5. Rows with fewer than 20 qualifying units are ineligible for the primary population-diversity fit. Polarized SFS-B is deferred, not rescued by choosing a larger sample ad hoc. |
| Typst and Lean on cluster | No. Tier 3 cluster jobs build data products only. Manuscript compilation and Lean validation remain separate repository validation steps and are not added to the Tier 3 Guix closure. |
| Congener substitution | No congener is substituted for a focal species in a primary result. A congener can be entered only as its own, separately labelled sensitivity row and cannot inherit the focal species’ diversity, annotation, or Buffalo predictor. Missing focal data means that observable is unavailable. |

Section 3 of the old execution plan proposed Conda/mamba.  That is also closed:
GNU Guix is the sole environment manager (section 9).  Conda, micromamba,
`pip --user`, environment modules that shadow Guix tools, and hand-built
containers are not approved.

## 2. Observables are distinct

The result table MUST use these names and meanings.  Values from different
rows below MUST NOT share one generic `pi` column.

| result name | unit and estimator | admissible source | interpretation |
|---|---|---|---|
| `population_pi` | Mean unbiased pairwise nucleotide differences per callable nuclear site among a declared population sample. At a callable site with called allele counts \(n_i\) and total \(n\), the contribution is \(1-\sum_i n_i(n_i-1)/(n(n-1))\); callable invariant sites contribute zero. | Exact-reference population VCF/BCF plus all-sites, gVCF, or exact-cohort callable mask. Phase is not required. | Population diversity. It is not an individual’s heterozygosity. |
| `individual_snv_heterozygosity` | Heterozygous single-nucleotide differences divided by callable/alignable A/C/G/T reference bases for one diploid individual. With two haplotypes, an accepted H1/H2 SNP contributes one and an invariant aligned base contributes zero. | Exact-reference individual VCF/BCF plus its callable mask, or the direct WFMASH H2-to-H1 fallback. | One individual, conditional on callability or assembly alignability. It is never renamed population π. |
| `pi_S_over_pi_W` | `pi_S / pi_W`, where `pi_S` and `pi_W` are the preceding source-appropriate estimator restricted to callable nuclear 4D sites whose **forward reference base** is S=G/C or W=A/T. | A diversity-eligible row plus native exact-reference annotation. | A descriptive, reference-conditioned 4D ratio. It is not mutation-direction counts and is not B. The reciprocal may be displayed but is never the stored statistic. |
| `polarized_sfs_B` | A demographic-model parameter inferred from a W→S versus S→W ancestral-state-polarized SFS. | None in `tier3-decisions-v1`. | **Deferred.** It requires a separately frozen outgroup, ancestral-error model, population design, demographic model, and power threshold. No reference-conditioned ratio may be relabelled B. |

Primary nucleotide diversity includes invariant callable sites, permits
multi-allelic SNVs with single-base A/C/G/T alleles through the formula above,
and excludes indels, structural variants, organellar sequence, ambiguous bases,
and sex/heterogametic chromosomes.  Indels and structural variants are retained
as QC outputs but are not in a primary numerator or denominator.

For every estimate, store the diversity sum (the sum of per-site pairwise
contributions), callable site count, nominal and observed chromosome counts,
point estimate, bootstrap interval, source modality, and exclusion counts.
`pi_S_over_pi_W` is unavailable when either class has fewer than 10,000
callable sites or when `pi_W` is zero; it is never coerced to zero or infinity.

## 3. Native annotation and exact coordinate policy

Primary GC3 or 4D analysis has a hard provenance gate.  The annotation MUST be
native to the exact reference FASTA assembly accession **including version**.
The manifest records:

- annotation provider and release;
- assembly accession and version on both FASTA and GFF/GTF records;
- SHA-256 and byte length of FASTA, annotation, and indexes;
- provider sequence-region declarations, the FASTA contig dictionary, and an
  explicit contig-name mapping checksum;
- `native`, `projected`, or `predicted` status and the nuclear genetic code;
- results of exact dictionary validation and sampled CDS reconstruction.

The primary status MUST be `native`.  A lifted/projected annotation or de novo
gene prediction is a separately named sensitivity analysis, with its own
validation and coordinate provenance, and never fills a primary GC3/4D cell.
If native annotation is absent, GC3 and 4D results are unavailable; whole-genome
GC of the exact FASTA may remain eligible.  A congener annotation is never
native.  For H1/H2, H1 carries its native annotation and only H1 coordinates
and bases are projected through the H1-to-H2 alignment; the GFF itself is not
lifted onto H2 for a primary result.

Before analysis, validators MUST establish all of the following:

1. Every annotation `sequence-region`/contig used by a retained feature maps
   one-to-one to an exact FASTA contig with the same length. Undeclared aliases
   are errors, not best-effort name stripping.
2. At least 100 deterministically selected retained transcripts, or all when
   fewer exist, are reconstructed from the exact FASTA and compared with a
   provider transcript/CDS sequence when one is deposited. The selection is
   by ascending SHA-256 of `dataset_id + transcript_id`. Any base discrepancy
   blocks annotation-derived primary results.
3. All retained CDS, not merely the sample, pass strand-aware phase,
   translation, internal-stop, and length checks described below.

### Canonical transcript, phase, and code

- Retain protein-coding, non-pseudogene nuclear features only. Choose the
  provider-designated canonical transcript when exactly one is declared.
  Otherwise choose the valid transcript with the longest translated CDS;
  ties are resolved by bytewise ascending stable transcript ID.
- Order CDS segments in biological 5′→3′ order. Interpret GFF3/GTF phase as
  the number of bases to skip at the biological start of that segment to reach
  the next codon boundary, including on the minus strand. Concatenated CDS
  length after phase handling MUST be divisible by three.
- Exclude the terminal stop codon from GC3 and 4D. Reject a transcript with an
  internal stop, inconsistent phase, missing/duplicate segment, frameshift,
  or ambiguous codon. Do not repair it heuristically.
- The manifest explicitly names the NCBI nuclear translation table. Use table
  1 unless a native provider annotation explicitly and validly declares a
  different nuclear code. Mitochondrial, plastid, and other organellar CDS are
  excluded.
- A third position is 4D only when substituting each of A, C, G, and T under
  the frozen translation table yields the same non-stop amino acid. If
  overlapping retained canonical CDS disagree on frame or 4D status, exclude
  the genomic site; concordant duplicates are counted once.
- W/S class comes from the forward reference FASTA base. Complementing a
  minus-strand CDS preserves W versus S, so no strand-specific relabelling is
  needed.

GC3 is the pooled fraction of G/C third positions across retained canonical
CDS, with numerator, denominator, genes, transcripts, and exclusions reported.
Whole-genome GC is reported separately and is not substituted for GC3.

## 4. Required same-reference tuple and invariant denominator

An **eligible diversity row** is one checksum-locked tuple:

1. reference FASTA and indexes;
2. native GFF3/GTF and its sequence-region/contig map;
3. normalized VCF/BCF and index, **or** raw WFMASH PAF plus H2 FASTA;
4. an explicit all-sites, gVCF, cohort callable mask, or H1-reference
   alignable mask;
5. sample/population/ploidy declaration and exact sample-list checksum;
6. assembly accession+version and a successful exact-reference audit.

All payloads have SHA-256 and size; upstream MD5 may be retained as metadata but
does not replace SHA-256. VCF `##contig` names/lengths and every REF allele are
checked against FASTA. GFF and callable coordinates are validated against the
same dictionary. A row with a mismatch is ineligible; liftover is not automatic.

Allowed denominator kinds are frozen as follows:

| denominator kind | required semantics |
|---|---|
| `all_sites_vcf` | Contains invariant and variant genotype records for the exact selected samples, with sufficient genotype/reference-confidence fields to apply the declared QC. Presence of a record alone is not proof of callability. |
| `gvcf` | Per-sample or joint reference-confidence blocks for the exact selected samples. Expand blocks only after frozen DP/GQ/filter rules; intersect called samples to form the cohort denominator. |
| `cohort_callable_mask` | A deposited or reproducibly generated BED/bitmap for the exact assembly and selected cohort. The generator, filters, input checksums, sample-list checksum, and mask checksum are mandatory. A generic assembly “accessible genome” mask may be an additional exclusion but cannot stand in for genotype callability unless its provider explicitly gives cohort-callable semantics. |
| `h1_reference_alignable_mask` | Generated by the frozen direct WFMASH CIGAR traversal in section 6. It includes accepted invariant `=` bases and accepted `X` bases, not just variant records. It is alignability-conditioned and used only for individual heterozygosity. |

A sparse variant-only VCF/BCF without one of the first three denominator sources
is **ineligible** for π and `pi_S_over_pi_W`. It may be inventoried, but absent
records are never assigned zero. IMPG VCF is also sparse and never supplies a
denominator.

## 5. Population, sample, ploidy, and downsampling rules

### Population datasets

1. Use one predefined biological population per dataset: the largest single
   named locality/population having at least 20 unrelated QC-passing sampling
   units and an exact cohort denominator. Do not pool localities. Ties are
   broken by bytewise ascending provider population ID.
2. Exclude crosses/progeny, duplicate/related samples, laboratory controls,
   contaminated samples, and samples failing provider QC. Record the exact
   exclusion list checksum, without committing protected identifiers.
3. For wild diploids, rank stable public sample IDs by
   `SHA256(dataset_id + "\0" + population_id + "\0" + sample_id)` and take
   the first 20. Use autosomes and 40 nominal chromosomes. This hash ranking,
   not an implementation RNG, is the downsampling rule.
4. For established inbred-line panels such as DGRP, select 20 lines by the
   same hash rule. Treat each line as one haploid consensus: a heterozygous or
   uncertain line call is missing, never arbitrarily phased or randomly
   sampled. At least 18 of 20 line alleles must be called at a site.
5. For diploids, at least 36 of 40 chromosome copies must be called at a site.
   The unbiased formula uses the actual called count. Site callability is
   recomputed for the selected sample set; no post-hoc substitution occurs.
6. Primary cross-species population π excludes selfers, mixed/unknown ploidy,
   haplodiploid mixtures, clonal/cyclical-parthenogen panels, pooled sequencing,
   and fewer-than-20 designs. Such resources may support explicitly labelled
   exploratory values after a new design decision, and may still support
   exact-assembly composition.
7. A haplodiploid resource can become eligible only from 20 unrelated diploid
   females with sex verified and a female-specific cohort mask. Haploid males
   are not mixed with them. Sex chromosomes are excluded in all primary rows.

Population structure is not “corrected” by pooling. If multiple populations
qualify, the frozen selected population is the primary row; other populations
are sensitivity rows. Population π does not require phased genotypes.

### Assembly individuals

H1 and H2 MUST be phased haplotypes from the same diploid individual, with
stable sample/haplotype identity and no evidence of haplotype collapse,
contamination, duplicated haplotigs, or gross H1/H2 coverage imbalance. The
result is `individual_snv_heterozygosity`. Family/trio origin, sex, phase method,
assembly QV, phased-block evidence, and all QC outcomes are recorded. A failed
phase/identity audit makes the assembly diversity row ineligible, although H1
composition may remain eligible.

## 6. Frozen modality hierarchy

### 6.1 Primary: deposited exact-reference variants plus mask

The primary route is a deposited VCF/BCF for the exact FASTA reference plus an
explicit denominator from section 4. Normalize/split representation against
the locked FASTA, retain original and normalized checksums, and compute allele
counts independently from normalized genotypes. Deposited individual calls do
not need phase for heterozygosity; deposited population calls do not need phase
for π. Reference bias and provider filters remain reported limitations.

### 6.2 Assembly fallback: direct pinned WFMASH parsing

When no defensible deposited call+mask tuple exists, a phased H1/H2 individual
may use direct WFMASH parsing:

- H1 is target/reference and carries the native GFF; H2 is query. Run the
  Guix-store WFMASH pinned in section 9 with base alignment enabled and an
  extended `cg:Z` CIGAR. `M`, missing `cg`, approximate `-m`, self mappings,
  secondary/multiple coverage, and unvalidated filtering are rejected.
- Production starts with `-p 90 -s 5k -l 25k -4`; the complete command,
  store path, raw PAF checksum, and WFMASH metrics are manifest fields. A
  threshold change creates a new policy version.
- Traverse both FASTAs and CIGAR, including reverse-strand records. Retain a
  target base only when it has exactly one accepted H2 projection, both bases
  are unambiguous A/C/G/T, and the operation is `=` or `X`. Exclude gaps,
  overlapping/multiple mappings, repeats declared by the H1 accessibility
  mask, 100 bp at each alignment edge, and 10 bp on each side of an indel or
  structural-breakpoint anchor. Preserve reason-coded exclusion masks.
- Count only single-base `X` operations in the heterozygosity numerator;
  accepted `=` and `X` target bases form the denominator. Indels/SVs are QC.
  The estimator is explicitly alignability-conditioned and may be biased low
  in divergent/rearranged regions.

### 6.3 IMPG: optional orthogonal concordance only

IMPG is not a primary caller or denominator. Version `tier3-decisions-v1`
allows it only on predeclared QC loci/individuals after a Guix package passes
the synthetic truth fixture. Frozen behavior is:

- input PAF must pass the direct WFMASH gate and use extended `= X I D` CIGAR;
- query H1 in disjoint 1,000,000-bp ownership cores with 10,000-bp padding on
  both sides, clipped at contig ends; `--merge-distance 0`;
- `query -o vcf:poa`, then `lace -t 2` or more, then FASTA-based
  `bcftools norm`, BCF conversion/indexing, exact duplicate removal, and only
  then keep records whose normalized 0-based anchor (`POS-1`) belongs to the
  core. Every eligible edge variant must be owned exactly once;
- the direct PAF callable BED remains mandatory. IMPG output never defines
  invariant bases;
- audit PanSN sample/ploidy and phase against known H1=REF/H2=ALT sites. The
  validated build emitted `1|0` where path order implied `0|1`, so IMPG phase
  orientation is **not trusted** and no phase-sensitive result is eligible;
- the observed `lace -t 1` deadlock makes one-thread lace forbidden.

The heavier WFMASH→seqwish→vg-deconstruct path and SweepGA are not approved
estimators. They require a future policy version if graph-complex variation
becomes an explicit objective.

### 6.4 Common-callable concordance gates

Compare modalities only after exact normalization and restriction to the same
reference coordinates and the **intersection** of their callable masks.
All gates below must pass; report the un-intersected mask sizes as well.

| metric | frozen pass tolerance |
|---|---|
| Common-mask coverage | Intersection is at least 80% of each modality’s eligible callable bases and at least 10,000,000 total bases (or the entire contig for a synthetic fixture). |
| Exact SNV allele/genotype agreement | Precision ≥0.99, recall ≥0.99, and heterozygous/non-heterozygous genotype concordance ≥0.99. Synthetic truth requires 1.000 for all three. |
| Total heterozygosity or π | Absolute difference ≤ `max(5e-5, 0.05 * mean_of_pair)`. |
| `pi_S` and `pi_W` | Each absolute difference ≤ `max(5e-5, 0.05 * mean_of_pair)`. |
| `pi_S_over_pi_W` | Absolute difference ≤ `max(0.05, 0.10 * mean_of_pair)`. |
| IMPG boundary/dedup QC | 100% of synthetic/core-edge truth variants retained exactly once; zero duplicate normalized alleles. |
| Phase audit | 100% expected orientation before any phase-sensitive use. Current IMPG fails, so phase-sensitive eligibility remains false. |

A concordance failure is recorded, investigated by exclusion stratum, and does
not get averaged away. It blocks promotion of the fallback for the affected
row; it does not rewrite the primary deposited estimate.

## 7. Current diversity inventory gate

The availability survey names promising resources but did not establish exact
assembly tuples or invariant denominators. Therefore no surveyed diversity row
is pre-approved merely by this decision document. Foundations must populate
the schema and turn a row eligible only after all gates pass.

| candidate row | required variant/alignment source | denominator decision | frozen eligibility now |
|---|---|---|---|
| VGP phased individuals | Exact-H1 deposited individual VCF/BCF+mask; otherwise direct H2→H1 WFMASH PAF | Deposited individual callable mask or generated `h1_reference_alignable_mask` | **INELIGIBLE_PENDING_TUPLE** per individual until H1 FASTA, native H1 GFF, H2/calls, mask, accessions, checksums, and phase QC are locked. |
| Ag1000G gambiae/arabiensis/merus/coluzzii | Deposited VCF/BCF, not a streamed implicit “latest” Zarr | Exact selected-population all-sites/gVCF/cohort mask | **INELIGIBLE_PENDING_TUPLE**; survey-level SNP availability does not identify invariant sites or exact reference releases. |
| DGRP Freeze 2.0 | Deposited VCF/BCF, 20 hash-selected inbred lines | Exact-line all-sites/gVCF/cohort mask | **INELIGIBLE_PENDING_DENOMINATOR_AND_REFERENCE_AUDIT**. Variant-only Freeze 2 calls are not assumed callable elsewhere. |
| *D. simulans* 170-line panel | Deposited VCF/BCF, 20 hash-selected lines | Exact-line all-sites/gVCF/cohort mask | **INELIGIBLE_PENDING_DENOMINATOR_AND_REFERENCE_AUDIT**. |
| *D. pseudoobscura* panel | Deposited release VCF/BCF | Exact selected-population all-sites/gVCF/cohort mask | **INELIGIBLE_PENDING_STABLE_RELEASE_TUPLE**; a mutable/preprint pointer is insufficient. |
| *Aedes aegypti* 1206 panel | Deposited release VCF/BCF, one locality, 20 diploids | Exact selected-population all-sites/gVCF/cohort mask | **INELIGIBLE_PENDING_TUPLE_AND_STRUCTURE_QC**. |
| *Daphnia pulex* panels | Deposited calls | A qualifying exact-cohort mask | **INELIGIBLE_FOR_PRIMARY_POPULATION_PI** because the surveyed designs are cyclical-parthenogen/temporal; composition may qualify independently. |
| CaeNDR *C. elegans* | Deposited calls | Would require exact-release selected-isolate mask | **INELIGIBLE_FOR_PRIMARY_POPULATION_PI** because of predominant selfing; composition may qualify independently. |
| *Nasonia vitripennis* | Deposited calls if found | Exact mask for 20 verified diploid females only | **INELIGIBLE_PENDING_SEX_PLOIDY_AND_TUPLE**; mixed haploid/diploid data are forbidden. |
| *Apis mellifera* and raw-read-only Drosophila/other surveys | Raw reads only in the present inventory | None | **INELIGIBLE** in v1. Raw-read joint calling is not an approved shortcut to a denominator. A future calling policy is required. |
| Tier 3c single assemblies | No diversity source | None | **INELIGIBLE_BY_DESIGN** for all diversity observables; eligible only for exact-assembly whole-genome GC and, with native annotation, GC3. |

This table is deliberately conservative. “Pending” means “do not run the
diversity estimator,” not “fill later with an assumption.”

## 8. Bootstrap and reporting policy

- Use non-overlapping 1,000,000-bp blocks in reference coordinates, split at
  contig boundaries. A final partial block is retained. Only blocks containing
  at least one callable site for the statistic participate.
- Perform 10,000 chromosome/contig-stratified block-bootstrap replicates,
  resampling blocks with replacement within each contig while preserving the
  number of blocks per contig. Recompute numerators, denominators, π values,
  and the ratio in every replicate; never bootstrap precomputed site ratios.
- Seed a counter-based generator from the first 64 bits of
  `SHA256("tier3-decisions-v1\0" + dataset_id + "\0" + statistic_name)` in
  big-endian order. Record the full seed digest and RNG implementation/version.
- Report the point estimate and the 2.5th/97.5th percentile interval. If fewer
  than 20 eligible blocks or fewer than 10,000 callable S or W sites exist,
  report the point estimate when otherwise valid but mark the interval or ratio
  unavailable; do not switch to a site bootstrap.
- Bootstrap intervals describe genomic block uncertainty conditional on the
  frozen sample. They do not represent phylogenetic, sample-selection, or
  causal uncertainty. Cross-species fits must retain clade and modality labels.

Polarized SFS-B is deferred in full. No outgroup download, polarization,
demographic fit, or B column belongs in v1 outputs.

## 9. GNU Guix, source, access, and cluster policy

### Immutable pins

| component | frozen source identity |
|---|---|
| GNU Guix channel | `https://git.savannah.gnu.org/git/guix.git`, commit `44bbfc24e4bcc48d0e3343cd3d83452721af8c36`; authenticate from Guix introduction commit `9edb3f66fd807b096b48283debdcddccfea34bad` and fingerprint `BBB0 2DDF 2E5F 5EAA C74C 3F7F 27D5 86A4 552F 384F`. No inherited user channels. |
| WFMASH | `https://github.com/waveygang/wfmash.git`, commit `e040aa10e87cab44ed5a4db005e784be62b0bd21`, recursive submodule `deps/WFA2-lib` at `49c255df126ee536fe92caff7a9f7c183ec3ff29`; recursive Guix/NAR SHA-256 base32 `0xvrgb7aimsdyq7kn3lzf7rrxxljs54dv81kv3apsfpcyp6nxnwj` (hex `92db6ecdf5ec3a7dd5d833a0dd48d192f69ef3719f0e3b0ff64dd7a8ce7a7977`). The tested non-Guix binary SHA-256 `14a6d5c7ac7be8890e904d11121341df118fad0c11193d0e91a9899e18a53d60` is validation evidence, not an approved runtime path. |
| Optional IMPG source | `https://github.com/pangenome/impg.git`, commit `101df81eb28a809c8fac97d297acd9fcfbbfa048`; recursive submodules `vendor/gfaffix`=`460e0dd798a9da7d12aef4f9181419d71489da95`, `vendor/syng`=`68ac19745201a7d2a17d9bb190671ef7d3ac8c29`; recursive Guix/NAR SHA-256 base32 `1g7z0xwaryz6cbkc2mhsrysi9cvg8599k50n3qh6ny4i0rx1z4di` (hex `b1911f7a0691786b201e16949952416fb314b5cf1a56c1e662e6fbac7807ffbc`). Validated development binary SHA-256 `39f47bc41a4023a3227f692ade55fb7bea7f541020c24306a04effe8b9781df0` is evidence only. |

Optional IMPG remains execution-disabled until its Guix package vendors and
hash-locks the complete `Cargo.lock` graph. The frozen git revisions include
POVU `0c9f21a40e7755293f5b21f33a7e2c894dc6cb15`, seqwish
`cf256a3e223dba7f268b6d1306c5dccc01619c25`, SweepGA
`ddd31d39b6a68fc972025b048076032341b66835`, wfmash-rs
`d47b7e3eba7f1d0f2bfdb72629cc667c2d5f8382`, allwave
`6b2798aa9fd405abe680cfb1e644331eff9d58e0`, fastga-rs
`e5037d5ef818f0ed1eef68e5678324a3b6f9111a`, and the remaining lockfile
revisions. A branch reference, rustup build, PATH binary, or un-hashed Cargo
fetch cannot enable IMPG.

Use `guix time-machine -C analysis/guix/channels.scm -- shell -m
analysis/guix/manifest.scm --pure`. Foundations will create a persistent GC
root, record derivations and complete store paths, and invoke tools by those
paths. The user profile’s WFMASH 0.12.5 is specifically forbidden.

### Compute-node truth and store topology

Two SLURM compute jobs on 2026-07-13 froze the actual topology:

- job `1753211` on `octopus07.uthsc.edu` saw Guix 1.4.0 at channel
  `44bbfc24e4bcc48d0e3343cd3d83452721af8c36`, direct NCBI HTTPS egress, and the shared MooseFS mount. It
  also proved that the visible daemon socket refuses compute-node connections.
- `octopus01` realized `/gnu/store/rlb2gljax8lzmhhidbvbzp3al1ad1mww-bcftools-1.14`;
  job `1753216` on `octopus07` executed that exact shared-store binary and
  reported bcftools 1.14/HTSlib 1.16, while a compute-side `guix build`
  predictably failed with “Connection refused.”

Therefore `octopus01` is the only build/realization host. Compute jobs neither
build nor contact the daemon; they execute manifest-recorded shared store paths.
Before arrays, a compute smoke MUST execute WFMASH and the synthetic fixture
from the realized GC root. If shared-store execution later fails, the sole
compatibility fallback is a checksum-recorded `guix pack -f squashfs` built by
the same pinned time-machine environment, with any runtime itself declared by
Guix. A Dockerfile, hand-built Apptainer/Singularity image, Conda container, or
downloaded third-party image is forbidden.

### Scratch, egress, retention, and data access

- Stage raw inputs under `/moosefs/erikg/tier3data`, transient intermediates
  under `/moosefs/erikg/tier3scratch`, and results under
  `/moosefs/erikg/tier3out`. Create directories with user-only write access.
- Job `1753211` observed 101,887,888,916,480 bytes free on the shared
  207,344,084,451,328-byte MooseFS and created a sentinel that remained visible
  after job completion. No user quota-query utility is installed, so the
  project freezes a **500 GiB effective Tier 3 budget**, not a claim of an
  administrator-enforced quota. Every array submission MUST preflight at least
  536,870,912,000 bytes free and stop rather than exceed that budget.
- MooseFS demonstrates post-job persistence but no automatic deletion contract
  was discoverable. The project therefore owns retention: inputs and
  intermediates remain through result acceptance plus 30 days, then are
  deleted from scratch; manifests, logs, checksums, small results, and Guix
  provenance are retained permanently. A weekly budget/age report is required
  while runs are active. Nothing relies on `/tmp` surviving a job.
- Public, redistribution-compatible accessioned data requiring no credentials
  are approved. Restricted/controlled data are ineligible until a separate
  access decision names the DUA, secure location, and approved operators.
  Credentials, tokens, signed URLs, cookies, private sample identifiers, and
  raw data MUST NOT enter git, logs, the manifest, or result TSVs.
- Downloads use accessioned immutable URLs where possible, retry/resume, stage
  to a temporary name, verify SHA-256 and expected size, and atomically rename.
  Any mutable “latest” endpoint must be resolved to a release and checksum
  before use. Egress is never used to stream unstaged data into an estimator.

## 10. Machine-readable contract and implementation gate

`schemas/tier3_manifest.schema.json` is the machine-readable contract. It
freezes the manager, pins, observable names, sample counts, modality choices,
IMPG query/core policy, bootstrap policy, native-annotation gate, and
compute-preflight evidence fields. JSON Schema cannot prove that two files
really share coordinates, so implementation must additionally perform content
validation (dictionary, REF allele, CDS reconstruction, and mask bounds) and
write the audit results the schema requires.

A row is runnable only when:

1. it validates against the schema;
2. `diversity_eligibility` is `eligible` for a diversity run;
3. every referenced local object exists and matches SHA-256/size;
4. the exact-reference/native-annotation audits pass;
5. the sample and denominator construction audits pass;
6. the Guix store paths resolve to the frozen derivations on a compute node;
7. the scratch/egress preflight passes; and
8. no credential or raw payload is tracked by git.

Fail closed. An unavailable observable stays unavailable. That is a valid Tier
3 result; inventing reference sites, annotations, samples, or provenance is not.
