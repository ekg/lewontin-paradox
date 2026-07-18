# Comprehensive VGP callable-diversity and PSMC research program

**Design freeze:** 2026-07-18 UTC

**WG task:** `design-vgp-comprehensive`

**Execution boundary:** design, metadata reconciliation, and repository tests only. This task downloaded no biological payload, realized no new environment, and submitted no local, batch, or Slurm analysis job.

## Decision

The pilot is frozen as exactly ten primary same-individual VGP haplotype pairs and six pre-ranked, clade-matched alternates. The panel spans Mammalia, Aves, Reptilia, Amphibia, Actinopterygii, and Chondrichthyes; three genome-size strata; low, medium, and high pre-result diversity expectations; and both early CLR/trio and later HiFi/Hi-C assembly generations. The machine-readable rows are in `analysis/vgp_10_pair_manifest.tsv` and `analysis/vgp_10_pair_alternates.tsv`.

The core scientific unit is not an annotation, a population sample, or a pair of assembly labels. It is an exact, reciprocal, same-BioSample H1/H2 pair whose two assemblies remain mutually comparable after independent base-accuracy, completeness, duplication/collapse, contiguity, and phasing review. H2 is aligned to H1 across the eligible nuclear sequence. Only verified 1:1 sequence enters a reason-coded callable mask. Normalized differences within that mask estimate callable-genome heterozygosity, and the same variants plus the same mask construct the diploid consensus used by PSMC.

This is why a high-quality pair can contain usable heterozygosity even though each deposited assembly is labeled haploid: the biological differences are between the two homologous assemblies of one diploid individual. It is also why annotation is not a core gate. Gene models can partition a core result, but neither genome-wide difference counting nor PSMC requires gene labels. Conversely, annotation cannot repair a wrong individual, poor base accuracy, collapsed sequence, retained many-to-many alignment, or low callability.

PSMC from this construction is descriptive demographic evidence, not independent validation of the diversity estimate. Both values reuse the same individual, haplotypes, aligned differences, and callable mask. Cross-species modeling must retain that covariance; it may not put same-pair PSMC on one side of a regression and same-pair diversity on the other as though they were independent observations.

## Immutable evidence reconciliation

### Freeze 1 identity

The release membership anchor is exactly:

| Property | Frozen value |
| --- | --- |
| Repository | `VGP/vgp-phase1` |
| Commit | `dc1b2af5a7741b97d66fb10cb2bce97f41765cdf` |
| Catalog | `VGPPhase1-freeze-1.0.tsv` |
| SHA-256 | `9c58420484a8b76a2d6175b7c26bf709e68bdc726a67fc7541b8c2b5a2fc13a4` |
| Bytes | 327,466 |
| Physical lines | 717 |
| Header | 1 line |
| Data rows | 716 |
| Unique scientific names | 714 |
| Legitimate duplicate names | *Lophostoma evotis* and *Micronycteris microtis*, each twice |

The catalog row—not a deduplicated species list, a moving NCBI query, or the current UCSC hub—defines membership. The mirror and later scale-out must account for all 716 rows. The current official NCBI and UCSC records are evidence and transport, not authority to silently change this release.

### Repaired evidence

The prior repaired six-row effort resolved exact current GCA/GCF records, BioSamples, isolates, annotations, checksums, and resource obligations, but it remained a composition/annotation pilot. Its `NO_GO` was driven by unavailable enforceable quota and scope-specific caps. It downloaded zero new biological bytes and submitted zero work. Its absence of a validated H2, whole-genome mask, and consensus at that time is not evidence that an exact H1/H2 pair cannot support callable diversity or PSMC.

This design reuses the repaired identities for camel, mousebird, viper boa, hourglass treefrog, and horn shark. It adds exact records from the already archived VGP project report snapshot at `/moosefs/erikg/tier3data/tier3a-acquisition-20260716/authoritative_vgp_PRJNA489243_dataset_report.json`, SHA-256 `8de3b6a74fe13e92c1e863f153f200fa163d3da80a27935b26f90eb015d78460`. That snapshot has 2,793 assembly reports and records its three source-page URLs. Every selected pair has two exact GCA versions, one shared BioSample and isolate, and reciprocal linked-assembly evidence. Each individual report entry is canonically serialized and SHA-256-bound in the TSV.

The catalog and exact reports sometimes differ. Those discrepancies are retained rather than erased. Most importantly, the catalog calls *Inia geoffrensis* CLR/FALCON while its exact selected accessions report HiFi/Omni-C Hifiasm; therefore A01 is classified by the exact assembly reports and carries the catalog discrepancy. The horn-shark catalog calls the generation HiFi, while the report text says only PacBio Sequel; P09 remains provisional until acquisition resolves exact read chemistry. No pair or accession identity is unresolved.

### Verified SweepGA/IMPG handoff

The implementation contract is bound to `analysis/sweepga_impg_handoff.md` (SHA-256 `27c16b54f6e63aab5eb3790cc923295a1582c948adbb33c0dd58a12db4565363`) and `analysis/sweepga_impg_observed.json` (SHA-256 `ff088828a3dc7336beec4e727fa7598840ec0f8b1234758f64f4480319364b61`). The biological completion proof used an exact native-annotation *Spinachia spinachia* excerpt, observed query and target overlap depth one, selected native IMPG partitions, produced normalized indexed VCF and BCF, validated every REF against H1, and reconstructed aligned H2 from the observed alleles. IMPG emits VCF; `bcftools`, not IMPG, owns normalization, BCF serialization, and indexing.

## Frozen pilot and alternates

| Slot | Species | Clade | Generation | Max haplotype size | Size | Expected diversity | Exact pair |
| --- | --- | --- | --- | ---: | --- | --- | --- |
| P01 | *Camelus dromedarius* | Mammalia | later HiFi/Omni-C, trio | 2.307 Gb | medium | medium | `GCA_036321535.1` / `GCA_036321565.1` |
| P02 | *Pseudorca crassidens* | Mammalia | later HiFi/Omni-C | 2.679 Gb | medium | medium | `GCA_039906515.1` / `GCA_039906525.1` |
| P03 | *Colius striatus* | Aves | later HiFi/Hi-C, trio | 1.212 Gb | medium | medium | `GCA_028858725.2` / `GCA_028858625.2` |
| P04 | *Falco naumanni* | Aves | early CLR/TrioCanu | 1.216 Gb | medium | medium | `GCA_017639655.1` / `GCA_017639645.1` |
| P05 | *Candoia aspera* | Reptilia | later HiFi/Hi-C | 1.531 Gb | medium | medium | `GCA_035149785.1` / `GCA_035125265.1` |
| P06 | *Dendropsophus ebraccatus* | Amphibia | early CLR/TrioCanu | 2.353 Gb | medium | medium | `GCA_027789765.1` / `GCA_027789725.1` |
| P07 | *Spinachia spinachia* | Actinopterygii | later HiFi/Hi-C | 0.408 Gb | small | low | `GCA_048126635.1` / `GCA_048127205.1` |
| P08 | *Menidia menidia* | Actinopterygii | later HiFi/Hi-C | 0.571 Gb | small | high | `GCA_048628825.1` / `GCA_048544195.1` |
| P09 | *Heterodontus francisci* | Chondrichthyes | later Hifiasm/Hi-C, chemistry audit | 6.013 Gb | large | medium | `GCA_036365525.1` / `GCA_036365495.1` |
| P10 | *Hemiscyllium ocellatum* | Chondrichthyes | early CLR/TrioCanu | 4.149 Gb | large | medium | `GCA_020745735.1` / `GCA_020745765.1` |

Expected diversity is frozen before the whole-genome result. The low and high anchors use prior targeted gene-span estimates—0.000397957 for *Spinachia* and 0.0138517 for *Menidia*—only as ranking evidence. They are not silently promoted to genome-wide estimates, are not responses in the pilot, and will not be overwritten after whole-genome analysis. Other rows use broad pre-result life-history categories and explicit missingness rather than invented numerical priors.

Six alternates are ranked one through six: A01 *Inia geoffrensis* (mammal), A02 *Lonchura striata domestica* (bird), A03 *Anolis sagrei* (reptile), A04 *Xenopus petersii* (amphibian), A05 *Syngnathus typhle* (ray-finned fish), and A06 *Hydrolagus colliei* (cartilaginous fish). An alternate may be activated only when a primary fails a pre-result hard gate, only within the same clade, and only through a signed/versioned amendment made before that slot's diversity or PSMC result is unblinded. Activation never deletes or renames the failed primary row. If two primaries in a clade fail and only one alternate exists, the second remains an explicit failure; there is no post-result shopping.

## Hypotheses and estimands

The hypotheses are deliberately narrower than “explain Lewontin's paradox” in one pilot.

- **H1, compressed diversity range:** callable autosomal heterozygosity varies substantially across the pre-registered vertebrate strata but spans a narrower range than plausible long-term census-size differences. The null for the pilot is no monotonic association with the ordinal expected-diversity strata. Confirmatory scale-out waits for independent review.
- **H2, technical robustness:** after the same 1:1 and mask contract, diversity ranks are not driven by assembly generation, genome-size stratum, callable fraction, QV, contiguity, duplication, or collapse. A systematic generation-specific failure is a program stop, not a biological result.
- **H3, descriptive histories:** unscaled PSMC curves differ in shape among passing individuals and may identify candidate bottleneck/expansion regimes. They are descriptive same-pair summaries and cannot independently validate H1.
- **H4, annotation optionality:** core diversity and PSMC acceptance is invariant to annotation availability. Annotation affects only optional partitioned summaries.
- **H5, specialized gBGC claims:** allelic conversion or biased transmission is identifiable only in a directional pedigree/gamete design. H1/H2-only state counts cannot establish event direction or transmission bias.

Primary estimands are:

1. `callable_heterozygosity_snv = heterozygous normalized biallelic SNV sites / callable diploid H1-coordinate bases`;
2. `callable_diversity_all_variant_bp`, a separately labeled sensitivity numerator that includes indel-affected bases under a frozen accounting rule;
3. callable fraction and the complete, mutually exclusive exclusion-reason composition of the eligible H1 universe;
4. unscaled PSMC interval parameters and trajectory from the exact masked consensus;
5. block-bootstrap uncertainty for heterozygosity and PSMC, with blocks confined within contigs and continuous callable runs; and
6. scenario-scaled PSMC time and effective-size trajectories, each tied to an explicit mutation-rate and generation-time source and never replacing the unscaled object.

The unit of replication is one diploid individual/pair. H1-coordinate windows and bootstrap blocks quantify within-pair uncertainty; they are not independent individuals. Species-level scale-out must use hierarchical models with pair and technical strata, not pseudoreplication of windows.

## Confidence tiers

Manifest confidence is provisional until measured QC exists.

| Tier | Requirements | Permitted claim |
| --- | --- | --- |
| A | exact reciprocal pair; exact input digests; QV, completeness, duplication/collapse, contiguity, 1:1, callability, consensus, and reproducibility gates all pass | primary callable diversity and unscaled PSMC; scenario-scaled PSMC with explicit assumptions |
| B | all Tier A core gates pass, but raw-read/k-mer calibration or long-range switch evidence is incomplete | core estimate with widened technical uncertainty; no top-confidence cross-generation contrast |
| C | exact pair is resolved but any quantitative core gate is not yet measured | metadata/design candidate only; no diversity or PSMC result |
| X | wrong/unresolved identity, non-comparable pair, retained multiplicity, sequence-truth failure, or irreconcilable mask | excluded core row with reason; never replaced silently |

Hi-C, Omni-C, Bionano, 10X, or trio data strengthen long-range phase evidence. Hi-C is not a universal gate because trio binning or other evidence can support phasing, and Hi-C cannot establish base accuracy. Its absence alone does not force X; poor QV, collapse, wrong pairing, or inadequate callability does.

## Gates fixed before unblinding

### Zero-tolerance provenance and truth gates

Any one of these immediately fails the pair:

- anything other than 0 unresolved catalog row, TaxId, BioSample, individual/isolate, accession version, or haplotype role fields;
- failure of shared-BioSample plus reciprocal-linked-assembly evidence;
- input, catalog, metadata-snapshot, software, command, environment, or output digest drift;
- any retained query-overlap depth above 1 or target-overlap depth above 1 after SweepGA's native `--num-mappings 1:1`;
- any nonzero base assigned to two retained 1:1 blocks;
- any unexplained base in callable-mask universe reconciliation;
- any normalized VCF REF mismatch against H1;
- any failure to reconstruct the aligned H2 allele sequence from normalized variants;
- any masked or non-callable base encoded as homozygous H1 reference in the PSMC consensus;
- annotation accession/dictionary mismatch for annotation-branch output; or
- conflation of unscaled and scenario-scaled PSMC.

### Numerical assembly and pair gates

- Each haplotype consensus QV must be at least 40. QV is measured or frozen from an exact, method-described Merqury/k-mer report; a blank catalog QV is not imputed.
- BUSCO complete fraction must be at least 0.90 and missing fraction at most 0.05 on each haplotype, using the same lineage and version. The pairwise difference in complete fraction must be at most 0.05.
- BUSCO duplicated fraction must be at most 0.05 on each haplotype. The k-mer spectrum and copy-number audit must show no unresolved collapse/duplication mode; a method name such as `purge_dups` is supporting evidence, not a pass by itself.
- Each haplotype must be at least 250 Mb and contig N50 at least 1 Mb. The H1/H2 total-length ratio must remain within 0.80–1.25. These are coarse preflight gates; passing them does not bypass measured callability.
- The exact read chemistry and phasing evidence must be resolved before compute. P09 therefore cannot leave acquisition while “PacBio Sequel” versus catalog “HiFi” remains unresolved.
- Sex chromosomes, mitochondria, unlocalized/unplaced sequence, and assembly-specific anomalous contigs are excluded from the primary autosomal universe and reported separately. A pre-registered all-nuclear sensitivity may include eligible sex-linked sequence.

### Mapping, mask, and consensus gates

- Primary callability must be at least 0.60 of the declared eligible H1 autosomal universe and at least 100,000,000 bp. The 0.50–0.60 interval is sensitivity-only and receives core `NO_GO`; below 0.50 is X.
- At least 50 non-overlapping 1-Mb windows must each be at least 80% callable for window and bootstrap summaries.
- Every excluded base receives exactly one primary reason from the ordered vocabulary: `not_eligible_contig`, `organellar`, `sex_linked_primary_exclusion`, `unplaced_or_unlocalized`, `h1_gap_or_N`, `h2_gap_or_N`, `not_1to1`, `mapping_breakpoint`, `low_base_accuracy`, `repeat_or_low_complexity_primary`, `duplication_or_collapse`, `phase_uncertain`, or `other_predeclared`.
- Callable plus all primary exclusion intervals must be disjoint and sum exactly to the eligible H1 universe. Overlapping diagnostic flags are retained in a separate non-accounting table.
- Variants are left-normalized and decomposed against exact H1, exact duplicates are removed deterministically, and only sites entirely supported by the callable contract reach the consensus.
- Masked bases become `N`; they never become homozygous reference. Heterozygous SNPs are encoded with IUPAC diploid codes. The indel policy is frozen as masked ±10 bp in the primary PSMC input, with 0 and 50 bp flank sensitivities.

### Bootstrap, reproducibility, and resource gates

- Run exactly 200 primary block-bootstrap attempts per passing pair, exceeding the program minimum of 100. Blocks are 5 Mb within contigs and continuous callable runs; 1 Mb and 10 Mb are frozen sensitivities. At least 190/200 (95%) must produce finite PSMC output.
- Heterozygosity confidence intervals use the same 5-Mb constrained blocks with exactly 10,000 deterministic replicates and a manifest-derived seed. Bootstrap units never cross contig or mask discontinuities.
- Two sentinel pairs (P07 low/small and P08 high/small) are recomputed independently from immutable inputs in separate scratch directories; every normalized VCF, mask, denominator, consensus, and unscaled PSMC digest must be identical. All other pairs rerun one deterministic chromosome shard and must match byte-for-byte.
- Scheduler prediction records CPU-hours, wall time, peak RSS, scratch high-water mark, read/write bytes, metadata operations, and retry count. Across completed primaries, absolute percentage error must be no worse than 25% at the median and 50% at the 95th percentile for wall, CPU, RSS, scratch, and I/O.
- A job is stopped if measured RSS or scratch exceeds 1.5 times its approved high estimate, if it approaches its Slurm limit within 10%, or if shared-storage I/O exceeds the approved per-wave rate. It may be re-estimated once without changing scientific gates.

### Program decisions

- All 10 primary slots must end as pass or explicit failure.
- Core pilot GO requires at least 8/10 primary passes, representation among passers of mammals, birds, reptiles, amphibians, fishes, both early and later generations, and no hard-gate violation.
- No assembly-generation stratum may have a systematic technical failure (at least half its attempted members failing the same technical gate). If it does, stop scale-out and repair the method.
- Core, annotation, PSMC, direct, population, phylogenetic, and non-allelic branches each receive independent `GO`, `CONDITIONAL_GO`, `NO_GO`, or `NOT_RUN_DESIGN_ONLY` decisions. Missing annotation or independent Ne cannot change a core GO.

## Exact SweepGA-to-PSMC contract

1. Bind immutable H1 and H2 FASTA digests, sequence dictionaries, exact accession versions, and the frozen pair row. Orient H1 as reference and H2 as query.
2. Run the accepted SweepGA origin/main binary at commit `018e4ce49d2c125820e0ac50dc5feaa02d423683` using the native long option `--num-mappings 1:1`, not the rejected short `-n` spelling. Preserve full command, binary digest, log, and PAF digest.
3. Independently interval-sweep query and target coordinates. Reject any retained overlap depth above one. Construct a disjoint H1-coordinate 1:1 BED and preserve excluded mappings with reason.
4. Pass the exact PAF without conversion to `impg index`. Run `impg partition` and preserve its native `partitions.bed`. Select only native partitions that intersect eligible query regions; never manufacture substitute windows.
5. Run `impg query` with both exact FASTAs, the index, PAF, and selected native partitions. Run `impg lace` to restore source coordinates and combine regional VCFs.
6. Normalize against exact H1 with `bcftools norm -f H1 -m -any`, trim to the exact eligible 1:1 regions, remove exact duplicates deterministically, serialize BGZF VCF+TBI and BCF+CSI, and bind all digests.
7. Validate REF against H1 and reconstruct the aligned H2 sequence from variants. Any mismatch is fatal.
8. Intersect the verified 1:1 set with eligible contigs, non-gap/non-N sequence in both haplotypes, breakpoints, base-accuracy, low-complexity/repeat policy, duplication/collapse, and phase-confidence filters. Emit the disjoint callable BED, one reason-coded complement BED per primary reason, and the universe reconciliation JSON.
9. Count primary SNVs and sensitivity indels only inside the callable contract. Emit numerator, denominator, fraction, window table, and constrained bootstrap table.
10. Construct the H1-coordinate diploid consensus, encoding validated heterozygous SNPs, applying the frozen indel flank mask, and writing `N` outside callable sequence. Reconcile consensus callable positions to the callable BED exactly.
11. Generate PSMC input from that consensus and mask. Preserve the unscaled trajectory as primary, execute 200 5-Mb bootstraps plus frozen block-length sensitivities, and only then apply separately versioned mutation-rate/generation-time scenarios.

SweepGA owns whole-haplotype mapping and 1:1 limiting. IMPG owns indexing, partitioning, regional graph/VCF extraction, and lacing. `bcftools` owns normalization, trimming, BCF/VCF encoding, deduplication, and indexing. Mask/consensus code owns the denominator and diploid encoding. PSMC owns the unscaled inference. No stage silently assumes another stage's responsibility.

## Evidence branches and exact outputs

| Branch | Required observations | Identifiable output | Forbidden inference | Pilot state |
| --- | --- | --- | --- | --- |
| Core | exact same-individual H1/H2, whole-haplotype 1:1 alignment, mask | `core/pair_summary.tsv`, normalized variants, callable/exclusion BEDs, window and bootstrap tables, QC JSON | population allele frequencies, transmission direction, conversion events | execute for ten adjudicated slots |
| Annotation | exact H1-native annotation or audited liftover with dictionary match | `annotation/feature_callable.tsv`, 4D/CDS/intron/intergenic summaries, mapping audit | vetoing core; calling projected sequence native | optional subset after core |
| PSMC | exact core variants, mask, consensus | consensus FASTA/PSMCFA digest, unscaled `.psmc`, 200 bootstraps, scenario table | independent validation of same-pair diversity; unqualified absolute Ne/time | execute for core passers |
| Direct pedigree/gamete | complete pedigree or gametes with transmitted parental haplotypes | callable event/tract TSV, rate denominator, GC transmission distortion, crossover association | cross-vertebrate rate without transfer assumptions | separate authorized pilot; not H1/H2 |
| Population frequency-spectrum gBGC | multi-individual, same-population genotypes and polarized or uncertainty-aware WS/SW states | frequency bins, polarization QC, model-dependent B and WS/SW asymmetry | direct event count, parent-of-origin, tract rate | `NOT_RUN_DESIGN_ONLY` until new task |
| Historical phylogenetic | ortholog-controlled close-clade alignment plus outgroups | branch-specific WS/SW substitutions, clustered signatures, outgroup sensitivity | current transmission bias, direct events, current population B | separate authorized phylogenetic pilot |
| Non-allelic/paralog | copy-resolved paralogs/segmental duplications and copy-number nulls | paralog candidate tracts, copy-homogenization score, null/sensitivity table | allelic meiotic conversion or biased transmission | `NOT_RUN_DESIGN_ONLY` until new task |

H1/H2-only outputs may report unpolarized heterozygous WS and SW states and candidate clusters. They must be named `state_counts` or `candidate_clusters`, never direct conversion, conversion direction, transmission distortion, or gBGC strength. Without parents/gametes there is no transmitted direction; without multiple individuals there is no frequency spectrum; without an outgroup there is no ancestral polarization; without copy-resolved paralogs there is no non-allelic claim.

## Immutable manifests and schemas

The authoritative design objects are:

- `analysis/vgp_10_pair_manifest.tsv`, exactly ten rows, SHA-256 `bf6c9ff647aed332bfc002bf803e8307203b51432343f2eca6d95a6c80d82997` at design freeze;
- `analysis/vgp_10_pair_alternates.tsv`, six ranked rows, SHA-256 `55127b18f0f17f6673cc0367c60736207a7e2184198cd3855a7e6ea83f39c52e` at design freeze;
- `analysis/vgp_analysis_manifest.json`, program contract and gates;
- `analysis/schemas/vgp_analysis_manifest.schema.json`, JSON Schema; and
- `analysis/assert_vgp_comprehensive_design.py` plus its regression tests.

Every downstream amendment is append-only and must include parent-manifest SHA-256, changed fields, reason code, author, UTC time, pre-result/unblinded state, and replacement trigger. Accession versions never float. An object manifest records source URL, expected bytes, upstream digest/provenance, staged local SHA-256, post-promotion SHA-256, license, retrieval UTC, and immutable content-addressed path. A run manifest binds the design manifest, object manifest, commands, environment profile, binary digests, Slurm IDs, inputs, outputs, random seeds, telemetry, and final disposition.

## Data-root layout

All execution remains under the authorized root `/moosefs/erikg/lewontin-paradox-data/vgp/phase1-freeze-1.0`:

```text
manifests/
  source/                 # pinned catalog and frozen source inventories
  design/                 # copied design TSV/JSON/schema + amendments
  acquisition/            # object, read, checksum, license manifests
  runs/                   # immutable run manifests and decisions
objects/sha256/aa/<hash>   # read-only content-addressed payloads
views/accession/           # non-authoritative links by exact accession
views/version/             # non-authoritative links by manifest version
staging/acquisition/       # resumable source-relative partials
staging/partials/
quarantine/
locks/
logs/
pilot/
  inputs/<selection_id>/
  runs/<run_id>/<selection_id>/
    mapping/ impg/ variants/ masks/ consensus/ psmc/ annotation/ qc/ telemetry/
  outputs/<run_id>/core/
  outputs/<run_id>/psmc/
  outputs/<run_id>/annotation/
scale/
  waves/<wave_id>/
specialized/
  direct/ population/ phylogenetic/ non_allelic/
```

Scratch is job-local or in an explicitly approved scratch root. Final paths are populated only by checksum-verified atomic promotion. A failed job cannot write a “complete” marker. No unconstrained delete or rsync `--delete` is permitted against durable objects.

## GNU Guix and software identity

All commands run through the authenticated single channel in `analysis/guix/channels.scm`, Guix commit `44bbfc24e4bcc48d0e3343cd3d83452721af8c36`, file SHA-256 `45c055cd1d9010a72eacbb720037a22bccb2d8d6891dbd11b5d66365f29b3a17`.

- Metadata, schemas, tests, statistics, samtools/bcftools, and report generation use `analysis/guix/manifest.scm`, SHA-256 `2fb05e87aa2ac45ce51d4dcf93b232cb98627f525adace98357629ee3f15720a`.
- SweepGA reproducible build uses `analysis/guix/sweepga_origin_main_manifest.scm`, SHA-256 `ea9ae1ba3e51ac3302d93add158532befec8fb3c09d188f524ac29237bab17d1`, and the byte-identical binary SHA-256 `fa7f0edb9b7e275c288db254046020e136d4267dd5ee043379227ef80da0573b` recorded in `analysis/sweepga_origin_main_build.json`.
- SweepGA→IMPG execution uses `analysis/guix/sweepga_impg_smoke_manifest.scm`, SHA-256 `c0ef5afd6c988341da8446ff3f70af274dd12f5514bc053d3d4e6f0cbdcee521`, with SweepGA commit `018e4ce49d2c125820e0ac50dc5feaa02d423683` and IMPG commit `101df81eb28a809c8fac97d297acd9fcfbbfa048` / observed version 0.4.1.
- PSMC must be added by the implementation task as a source-pinned Guix package, with repository revision, source hash, build recipe hash, binary hash, and two-build reproducibility proof. Ambient `psmc`, containers, modules, and network-installed packages are forbidden.

The implementation task may extend manifests but may not change the scientific gates without a versioned design amendment made before results.

## Slurm strategy and resource prediction

No jobs are submitted by this design. The implementation will create a dry-run preflight and then use dependency-linked arrays:

1. one lightweight metadata/input-integrity array over ten primaries;
2. one mapping array, with a separate resource class for small (<1 Gb), medium (1–3 Gb), and large (>3 Gb) pairs;
3. IMPG partition/query shards by native partition groups, never arbitrary preprocessing windows;
4. pair-level lace/normalize/mask/consensus joins after all shards pass;
5. one PSMC primary array plus bounded bootstrap arrays in batches of 20; and
6. annotation arrays only for rows with exact annotation approval.

Initial concurrency is two small, one medium, or one large mapping pair at a time, subject to the stricter live allocation and measured storage/I/O limits. This is scheduler safety, not a scientific eligibility cap. Per-pair high estimates derive from exact compressed/uncompressed bytes, total sequence length, contig count, native partition count, and calibrated coefficients. The six-row pilot's historical 120 GiB/750 GiB/1,500 core-hour/256 GiB limits are not global ceilings.

Each `sbatch` is preceded by a no-network preflight that checks manifest digests, available enforceable quota, 25% storage headroom after high estimates, scratch, environment realization, exact binary hashes, scheduler limits, and absence of an existing success marker. Jobs stage read-only inputs to node-local scratch where possible, batch metadata operations, checkpoint native partitions, trap signals, and emit `sacct` plus application telemetry. Retries are capped at one for resource re-estimation and two for transient infrastructure failure; scientific failures are never retried with relaxed gates.

## Stop conditions

Stop the affected pair immediately on any hard gate, checksum mismatch, wrong accession, reciprocal-link loss, QV/completeness/collapse failure, multiplicity, REF/H2 reconstruction error, mask discrepancy, consensus discrepancy, or annotation dictionary mismatch. Quarantine corrupt objects; never remove the last verified copy.

Stop the entire pilot before further mapping if catalog/design/environment digests drift, any pair substitution lacks an authorized amendment, enforceable quota is unknown, projected headroom is below 25%, or the same unexplained pipeline failure occurs in two pairs. Stop scale-out if fewer than 8/10 primaries pass, any required clade or generation disappears from passers, resource prediction misses its program gate, either sentinel is not byte-reproducible, or technical failure is systematic by generation. Specialized branches stop independently and cannot downgrade a valid core result.

## Pilot-to-scale waves

- **Wave 0—implementation proof:** controlled truth plus the already verified *Spinachia* excerpt; no new biological acquisition.
- **Wave 1—small sentinels:** P07 and P08. Validate the expected low/high anchors, deterministic reruns, masks, consensus, telemetry, and PSMC mechanics before larger jobs.
- **Wave 2—medium later generation:** P01, P02, P03, and P05. Test cross-clade behavior with modern phasing.
- **Wave 3—early generation:** P04 and P06. Admit only if QV/completeness/collapse and 60% callability gates pass; compare technical residuals to later pairs.
- **Wave 4—large cartilaginous fish:** P09 and P10, one at a time after resource model recalibration. Resolve P09 chemistry before transfer or compute.
- **Independent review:** adjudicate all ten slots and branch decisions without relaxing gates.
- **Scale wave S1:** every Freeze 1 exact pair that passes closed-world metadata/pair preflight, grouped by size/generation and preserving all 716 catalog-row dispositions.
- **Scale wave S2:** optional annotation partitions for exact-native/dictionary-audited subsets only.
- **Scale wave S3:** phylogenetic clade pilots after mirror completion; direct pedigree/gamete pilot on its separate dataset. Population and non-allelic branches remain design-only until explicitly authorized.

Full scale-out is not a fixed species count. It is the closed-world set of every catalog-linked pair that passes the pre-registered exact-pair and quantitative gates, with every catalog row receiving `eligible`, `ineligible`, `superseded`, `missing`, or `not_a_pair` and a reason. No result-dependent alternate selection or annotation gate may change that universe.

## Authorization boundary

This document authorizes downstream implementation, selective acquisition, and mirror planning to use the frozen identities and gates. It does not authorize a download, raw-read transfer, Guix realization, Slurm submission, demographic inference, direct-event claim, population gBGC estimate, phylogenetic substitution result, or paralog-conversion claim. Each action requires its own immutable manifest, capacity evidence, license review, and task-specific GO decision.
