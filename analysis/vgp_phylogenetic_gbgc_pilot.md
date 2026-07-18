# VGP phylogenetic gBGC pilot: fail-closed execution report

**Disposition:** `NOT_EXECUTED_INPUT_GATE`
**Preflight:** pinned GNU Guix channel `44bbfc24e4bcc48d0e3343cd3d83452721af8c36`, immutable profile `/gnu/store/ff3msd4h7prl407s3d74a4i5cxbvn694-profile`, authorized Slurm job `1781129`
**Frozen VGP catalog:** commit `dc1b2af5a7741b97d66fb10cb2bce97f41765cdf`, SHA-256 `9c58420484a8b76a2d6175b7c26bf709e68bdc726a67fc7541b8c2b5a2fc13a4`

## Outcome

Neither H01 nor H02 produced a biological estimate. The upstream Freeze 1 mirror stopped at
`quota_visibility_unavailable_fail_closed`: zero bulk Slurm jobs launched, zero objects transferred,
and zero objects verified. H01 has released UCSC FASTA/2bit inventory entries for 3/5 exact
assemblies; H02 has entries for 1/5. All eight represented sequence objects remain `planned` with
zero observed bytes. The remaining exact close relatives or unreleased catalog entries have no
sequence object in the frozen released-source inventory. No sequence digest, sequence dictionary,
coordinate-compatible annotation, assembly-specific license audit, versioned/checksummed topology,
reciprocal 1:1 orthology product, alignment, callable base, or branch-polarized substitution was
therefore available.

Running phastBias or another nonstationary model against missing or floating inputs would violate
the design. The Slurm job performed metadata-only preflight generation and fail-closed validation;
it did not download sequence, align genomes, infer ancestors, shuffle biological labels, simulate a
fitted biological null, or estimate WS/SW asymmetry, historical bias, or clusters.

## Exact clades and polarization design

| Panel | Ingroup | Outgroup paths | Frozen identity result |
|---|---|---|---|
| H01 | *Spinachia spinachia* focal; *Pungitius pungitius*; *Gasterosteus aculeatus* | *Syngnathus acus* and *S. typhle* | 5 accessions frozen in the design/catalog; only 3 represented in released source inventory; 0 local verified |
| H02 | *Falco naumanni* focal; *F. tinnunculus*; *F. peregrinus* | *F. cherrug* and *F. punctatus* | 3 VGP catalog assemblies plus 2 exact close relatives frozen in design; only 1 represented in released source inventory; 0 local verified |

The row-level manifest records exact accession/version, catalog membership (including the H01
*S. acus* other-haplotype selection), source URL, metadata assembly span/fragmentation, annotation
release, license terms, inventory paths, checksum state, coordinate evidence, and exclusion reason.
The source inventory SHA-256 is `f2374dc274c0bb90ac233bdaa0c4e6ca979c977697b30507b93b73ca272404ac`; the mirror manifest SHA-256 is
`935d18268c8b3b2ca5eb68332282a39db29fef563164adc8b6a493131b4675a0`; the mirror summary SHA-256 is `801ab3ffce1eca0f43d08ac406ded2604839f783221b4ff91cdc856699596a4c`.

## Pre-registered gates and missingness accounting

The minimum is three callable ingroup species plus two callable outgroups. H_SUB additionally needs
10,000 callable neutral aligned bases and at least 20 WS plus 20 SW substitutions per interpreted
branch/partition. H_GBGC needs at least 100 directionally informative substitutions. H_CLUSTER needs
at least 50 high-posterior branch-polarized WS substitutions and 20 Mb of continuous callable
single-copy alignment. Callable input is zero because no verified sequence was available;
biological substitution counts are not estimable, so all three gates fail before modeling.

Assembly span and contig count are retained as design metadata, but alignment missingness,
fragmentation loss, callable orthology, and excluded bases cannot be calculated without sequence.
They are reported as `NOT_ESTIMABLE_NO_ALIGNMENT`, never as zero missingness. Exclusion is accounted
for by reason: `NO_VERIFIED_SEQUENCE_PAYLOAD`, `MISSING_RELEASED_SEQUENCE_SOURCE`,
`CHECKSUM_UNRESOLVED`, `ANNOTATION_DICTIONARY_UNAVAILABLE`, `NO_ORTHOLOGY_PRODUCT`, and
`NO_CALLABLE_ALIGNMENT`.

## Models and controls not run

All planned controls remain explicit `NOT_RUN_UPSTREAM_GATE` rows: negative controls, branch/taxon
label shuffling, context/opportunity-matched null simulation, alternate outgroups, ancestral-state
posterior and polarization-error sensitivity, alignment filters, 3-mer/5-mer mutation-bias models,
CpG sensitivity, GC equilibrium and multiple hits, ILS/introgression and local-gene-tree sensitivity,
branch length and phylogenetic correction, chromosome block resampling, recombination proxies,
chromosome size, local GC, paralogy/segmental-duplication exclusion, and Holm/BH/BY correction.
Neutral, fourfold, coding, and noncoding partitions are separately ledgered. No coding partition is
treated as neutral. No native annotation or liftover is admitted until its exact sequence dictionary
and coordinate mapping validate.

## Historical results and uncertainty

Genome-wide H_SUB and H_GBGC rows, tract-like H_CLUSTER rows, and fragmentation/missingness
sensitivity rows are present for both panels. Their counts, opportunities, estimates, intervals, and
callable denominators are `NOT_ESTIMABLE_INPUT_GATE`. This is not a null result and supplies no
evidence for or against historical gBGC. Candidate clusters are not conversion tracts even when a
future run can estimate them.

This branch can only estimate long-term, model-dependent substitution signatures. It does not
observe direct conversion events or biased transmission, does not estimate a pedigree tract rate,
and does not substitute for multi-individual population-frequency evidence.
Direct, population, and non-allelic estimands are explicitly `NOT_MEASURED`.
H1/H2 heterozygous WS/SW states were not used as direct evidence and could only enter a future
QC/sensitivity analysis.

## Safe activation rule

Do not rerun this blocked-output generator as a biological analysis. First obtain enforceable quota
evidence and an authorized selective acquisition amendment covering all ten exact sequences and
needed annotations; verify provider digests plus local SHA-256; bind native annotation or a separately
validated liftover to each exact sequence dictionary; build reciprocal 1:1 single-copy syntenic
alignments; account for excluded bases by reason; then submit a separately fingerprinted Slurm
analysis under a pinned environment that includes the declared nonstationary implementation. Any
accession substitution, outgroup change, or annotation remap requires a pre-result amendment.
