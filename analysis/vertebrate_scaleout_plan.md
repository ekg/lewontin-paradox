# All-vertebrate Tier 3 scale-out plan

Date: 2026-07-17 UTC  
WG task: `plan-all-vertebrate-tier3`  
Status: **execution-ready planning contract; no acquisition or execution authorization**

## 1. Decision and scope

This plan scales the validated Tier 3A/Tier 3C design only after a fail-closed
inventory, an explicitly authorized small pilot, and reviewed waves. It does
not authorize downloading the VGP catalog, staging a catalog of assemblies,
submitting Slurm jobs, or running demographic analyses. The only currently
resolved biological rows are the three completed fish tuples, and they are
calibration rows rather than members of the synchronized freeze's reported
40-row triple-eligible cohort.

The controlling local inputs are:

1. `analysis/vertebrate_scaleout_origin_inventory.md`, the synchronization and
   limitation audit; and
2. `results/tier3/vgp_freeze_analysis.md`, the guidance file it identifies,
   SHA-256
   `0ad6fa03ceeec9d07c39c5456ddc4702c54c586115e578af3599d07f29e5316d`.

The plan treats the guidance counts—714 unique species, 223 completed
assemblies, 248 with the named annotation statuses, 271 with a pair field, 120
completed plus annotated, 40 completed plus annotated plus paired, and 46
completed plus RefSeq-annotated fish—as **reported, not reproduced**. The raw
717-line TSV is absent locally. GitHub's official VGP repository page was
inspected on 2026-07-17 and displayed 717 lines and 320 KB, but its `main` URL
is moving and was not downloaded. Therefore current execution eligibility is
zero rows from the synchronized freeze. This is a blocker, not an estimate of
zero biological availability.

The plan has four concrete artifacts:

- this plan;
- `analysis/vertebrate_scaleout_candidate_schema.tsv`, a row contract seeded
  with all 13 named fish and three calibration-only rows;
- `analysis/vertebrate_scaleout_resource_budget.tsv`, observed calibration
  plus low/base/high estimates; and
- `analysis/vertebrate_scaleout_wg_graph.md`, an inert proposed execution
  graph with explicit authorization tasks.

## 2. Controlling-guidance crosswalk

Every actionable statement in the synchronized inventory and its guidance is
mapped below. `I` refers to a named plan element, `O` to an output contract,
and `V` to a validation gate. A blocker ID means the requirement is not yet
satisfied and must not be replaced by an assumption.

| ID | controlling requirement or reported claim | plan element | output | validation gate / disposition |
|---|---|---|---|---|
| G01 | freeze source is `VGPPhase1-freeze-1.0.tsv`, described as 717 lines | I01 content-addressed catalog acquisition | O01 `source_catalog.json` and raw immutable TSV | V01 exact 717 lines, 716 data rows, header hash; B01 until authorized retrieval |
| G02 | guidance says the raw TSV is committed, while synchronization proves it is absent | I01 source reconciliation | O01 includes local-Git and upstream-revision evidence | V01 stops on absence or contradiction; B01 |
| G03 | moving `main` fetch command and placeholder `2025-01-XX` date are inadequate | I01 immutable revision and real UTC timestamp | O01 URL, Git object/revision, bytes, SHA-256, retrieval time | V01 rejects `main` without resolved commit and rejects placeholder dates |
| G04 | reported counts 714/223/248/271/120/40 | I02 deterministic catalog recomputation | O02 `catalog_counts.json` plus complete acceptance/rejection ledger | V02 all counts match or reviewed discrepancy stops the program; B02 |
| G05 | triple-eligible taxonomic counts are fish 13, amphibian 4, reptile 3, mammal 9, bird 6, other 5 | I02 taxonomy-aware recomputation | O02 count table keyed by TaxId and lineage snapshot | V02 sum 40 and each group reproduced; B02 |
| G06 | catalog screen uses columns 10, 13, 16, 17, 21, 26 | I02 versioned parser with header-name and ordinal assertions | O02 parser commit/hash and field-level derivation | V02 reject schema drift, duplicate species ambiguity, or silent column fallback |
| G07 | 13 fish H1 accessions and H2 labels are named | I03 seed inventory and resolver | O03 candidate TSV contains all 13 | V03 no seed omitted; labels are not accessions; B03 |
| G08 | 46 Tier 3C fish are reported but no row list is available | I02 catalog derivation, not hand reconstruction | O02 emits the exact 46-row subset | V02 must reproduce 46; B04 |
| G09 | versioned H1 and H2 accessions, reciprocal relationship, same individual, roles, and release metadata are required | I04 pair identity preflight | O04 `pair_evidence.json` per row | V04 100% required fields and reciprocal/same-individual evidence; B03 |
| G10 | provider hashes and repository SHA-256 are required | I05 two-level checksum capture | O05 asset manifest | V05 provider checksum plus local SHA-256 and canonical sequence-set digest agree |
| G11 | a pair label alone is not callable diploid evidence | I06 assembly-callability contract | O06 callable BED and denominator record | V06 positive, measurable H1-reference denominator; no inference from pair status |
| G12 | native annotation must refer to identical H1 sequence content | I07 exact-reference annotation preflight | O07 annotation/reference linkage and contig map | V07 sequence-set identity, native status, contig bijection, reconstructed CDS pass |
| G13 | projected or accession-mismatched annotation is not a primary input | I07 no-lift policy | O07 explicit `native/projected` status | V07 any projected/lifted/mismatched row is rejected, never repaired by substitution |
| G14 | genetic code, sequence regions, and CDS reconstruction must be audited | I07 annotation audit | O07 genetic-code field, sequence dictionary, CDS audit | V07 nonzero valid CDS, phase consistency, stop on reference disagreement |
| G15 | three completed tuples are calibration only and do not satisfy freeze triple rules | I08 calibration isolation | O03 rows have `record_role=calibration_only` | V08 cohort collector rejects calibration IDs from scale-out estimands; B06 |
| G16 | GNU Guix channel and manifest are controlling | I09 pinned environment | O08 environment/derivation/closure record | V09 exact hashes and store paths; no ambient tools |
| G17 | general manifest does not contain executable IMPG | I09 combined production-lineage capture | O08 separate supplemental profile and binary provenance | V09 an IMPG discovered outside the approved record is fatal; B08 until production capture is reviewed |
| G18 | SweepGA performs whole-H1/H2 mapping with native `--num-mappings 1:1` | I10 mapping stage | O09 bounded PAF, command, native multiplicity audit | V10 whole FASTA inputs; query and target overlap depth at most one |
| G19 | `--scaffold-jump 0` leaves partition ownership to IMPG | I10 responsibility boundary | O09 mapping command | V10 command fingerprint includes `--scaffold-jump 0` |
| G20 | IMPG owns index, native implicit-graph partitioning, and annotation-selected regional VCF queries | I11 IMPG stage | O10 graph index, native partitions, focus BED, regional VCFs | V11 focus derives from IMPG partitions intersected with original H1-native targets |
| G21 | IMPG emits VCF; bcftools creates normalized VCF/BCF and indexes | I11 serialization boundary | O10 laced VCF, normalized VCF/TBI, BCF/CSI | V11 REF equals H1; indexes query; formats contain the same normalized records |
| G22 | original H1-native annotations define targets and denominators | I12 summary stage | O11 target BEDs and denominator table | V12 no alternate annotation or target lift enters the lineage |
| G23 | metadata preflight precedes bulk staging | I13 staged governance | O12 metadata-only candidate freeze | V13 all pilot rows resolve and budget before any asset download |
| G24 | stage exact assets in bounded arrays and retain failures | I14 staging array | O05 asset ledger plus failure ledger | V14 100% input rows end accepted/rejected/unavailable; never shrink silently |
| G25 | truth case and stratified pilot precede expansion | I15 eight-slot pilot | O13 pilot review packet | V15 quantitative thresholds in Section 11 |
| G26 | reviewed waves use dependency gates, immutable input, local scratch, atomic promotion, retries, and `sacct` | I16 resumable Slurm design | O14 run ledger and telemetry | V16 thresholds in Sections 8 and 11 |
| G27 | Tier 3A and Tier 3C eligibility remain separate | I17 modality state machine | O03 four independent eligibility flags | V17 collectors select explicit modality status, not a general `ready` flag |
| G28 | population and demographic validity cannot be inferred from H1/H2 | I18 analysis-tier separation | O03 population/demographic fields | V18 deposited matched data and method-specific inputs required; otherwise ineligible |
| G29 | every eligible, failed, and unavailable row is collected | I19 closed-world collection | O02/O14 row-count reconciliation | V19 input IDs equal accepted + rejected + unavailable + failed IDs exactly |
| G30 | independent QC and resource review precede full catalog | I20 reviewed waves | O13/O14 review decisions | V20 no full-catalog task becomes ready before `authorize-full-vertebrate-execution` |
| G31 | reported source is PRJNA489243 and NCBI/VGP resources require accession-level resolution | I03/I05 primary-repository retrieval | O03/O05 source URLs and accessions | V03/V05 project IDs never substitute for versioned row assets |
| G32 | planning completion is not execution approval | I21 authorization firewall | O15 WG graph | V21 all expensive nodes depend on an explicit human authorization task |

## 3. Blocker register

| ID | blocker | consequence | resolution evidence |
|---|---|---|---|
| B01 | raw freeze TSV is absent locally; moving upstream revision and SHA-256 are unresolved | no catalog parser run and no bulk acquisition | approved metadata-only retrieval, immutable Git commit/object, URL, UTC date, byte size, SHA-256, header and line counts |
| B02 | all headline and taxonomic counts are reported, not locally reproduced | no claim of 40 or 120 execution-eligible rows | independently reproduced counts with discrepancy review |
| B03 | 13 seed fish lack H2 accession versions, pair/individual proof, exact checksums, and exact annotation tuples | all 13 remain `blocked_unresolved` | complete O03/O04/O05/O07 fields and pass V03–V07 |
| B04 | 46-fish composition row list is absent | no 46-row manifest | deterministic catalog query reproduces and emits it |
| B05 | row-level reuse terms were not recorded in guidance | no asset staging for a row with unresolved terms | repository policy URL, retrieval date, any submitter-specific notice, reviewer decision |
| B06 | the three completed species contradict the freeze triple-candidate screen | calibration cannot populate the target cohort | keep `calibration_only`; do not waive freeze rules |
| B07 | Tier 3A run omitted actual RSS/CPU, scratch peak, filesystem bytes, metadata operations, and bandwidth | resource estimates in those dimensions are planning bounds | pilot captures cgroup/`sacct`, `/usr/bin/time -v`, local du/inodes, and MooseFS counters |
| B08 | existing profiles and binaries are pinned, but there is no single newly captured production combined-profile review packet | no production lineage realization for scale-out | create O08 with all derivations, closure, GC root, binary hashes, two-build comparison, and smoke results |
| B09 | cluster `JobAcctGatherType` left `sacct` MaxRSS and TotalCPU blank for Tier 3A | cannot calibrate actual memory or CPU efficiency | administrator-enabled accounting or sidecar cgroup/getrusage capture validated against `sacct` |
| B10 | non-fish names in the reported 40-row cohort are unavailable | pilot slots outside fish cannot yet be named | resolve raw catalog and select deterministic strata |
| B11 | the final individual human assignees for reviewer roles are not stored in this plan | authorization tasks remain non-ready | WG authorization tasks record assignee, decision, timestamp, envelope, and signature/message reference |
| B12 | no frozen dated species phylogeny matching the final TaxIds exists | no confirmatory phylogenetic model | reviewed tree source, topology/version/hash, taxon reconciliation, branch-length policy |

`UNRESOLVED` is a valid inventory value only before an acceptance gate. It is
never coerced to an empty string, a current/latest accession, or a taxonomic
synonym.

## 4. Candidate inventory and eligibility

### 4.1 Expected scale and current state

The inventory procedure expects 716 source rows and 714 unique species. Its
reported candidate strata are 40 Tier 3A catalog candidates and 120 Tier 3C
catalog candidates across taxa; the 40 are expected to overlap the 120. There
are 13 named fish in the 40 and 46 reported fish in the 120. These are expected
scales, not guaranteed eligible totals. The number passing exact-reference,
native-annotation, and denominator gates is unknown and may be much smaller.

The seed TSV is a real row-level starting inventory: all 13 named fish are
present, with their reported H1 accession and pair label, and every unavailable
field is explicit. The three completed fish carry exact checksums and telemetry
links but are `calibration_only`. No seed row is currently accepted for
scale-out execution.

### 4.2 Inventory procedure

The following procedure is executable only after the small catalog-acquisition
authorization in the WG graph. It is metadata-first and does not stage genome
FASTA/GFF data.

1. Resolve the VGP repository's `main` ref to an immutable commit using the
   repository API. Retrieve exactly the TSV object at that commit to
   `sources/<commit>/VGPPhase1-freeze-1.0.tsv.partial`, record HTTP URL/status,
   ETag if supplied, UTC start/end, byte size, and SHA-256, then rename. Never
   use the moving raw URL as the frozen identity.
2. Assert 717 lines, 716 data rows, unique header names, expected ordinal/header
   pairs for columns 10/13/16/17/21/26, UTF-8, and no embedded delimiter
   corruption. Record the header SHA-256 and parser Git commit.
3. Emit one ledger row per source row and a species-level reconciliation table.
   Duplicate species are retained with source row numbers; deterministic
   aggregation produces the 714-species denominator.
4. Resolve the name through NCBI Taxonomy and retain TaxId, current name,
   source name, synonyms/merged IDs, full lineage, and metadata retrieval time.
   NCBI documents TaxId as a stable unique numerical identifier, while merged
   identifiers remain searchable; therefore both submitted and current IDs are
   retained ([NCBI Taxonomy](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/data-processing/taxonomy-processing/taxonomy/), inspected 2026-07-17).
5. Query NCBI Datasets metadata by the **exact versioned** H1 accession. Resolve
   the H2 label to a versioned accession only with reciprocal assembly and
   BioSample/individual evidence. Record release status, dates, assembly units,
   lengths, contiguity, submitter, BioProject, BioSample, and retrieval URL.
6. Resolve annotation releases whose declared assembly is H1 or its NCBI
   GCA/GCF paired assembly with identical component sequences. Older releases
   remain addressable where NCBI exposes them; no newer release is silently
   selected. NCBI states that assembly accession.version identifies a precise
   sequence set and increments on sequence updates
   ([NCBI assembly model](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/data-processing/policies-annotation/data-model/), inspected 2026-07-17).
7. Inventory deposited VCF/BCF, genotype, mask, raw-read, and sample resources
   without presuming compatibility. Record accession, version, reference,
   sample unit/population, URL, provider checksum, expected bytes, and terms.
8. Use NCBI Datasets `--include none` or dehydrated/preview functionality for
   metadata and size planning. NCBI supports accession-specific packages and
   metadata-only packages ([genome download documentation](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/how-tos/genomes/download-genome/), inspected 2026-07-17). Full assets wait for a separate staging authorization.
9. Recompute all reported counts and write every acceptance/rejection reason.
   A mismatch stops the freeze; a reviewer may accept a documented upstream
   change only by creating a new inventory release, never by changing the
   expected result in place.

### 4.3 Exact-reference identity

For each FASTA, store the provider MD5, compressed-byte SHA-256,
uncompressed-byte SHA-256, `.fai` SHA-256, and a `sequence_set_sha256` computed
from a sorted stream of `(sequence accession.version, length,
SHA-256(uppercase sequence))`. Byte identity alone is insufficient when
providers alter wrapping or headers; sequence-set identity alone is
insufficient for retrieval provenance, so both are required.

An annotation passes only when:

- its release metadata declares H1 or the paired RefSeq/GenBank assembly;
- the annotation reference and staged H1 have identical component sequence
  digests after an explicit one-to-one accession map;
- every GFF `sequence-region` maps bijectively to the H1 dictionary with the
  same length;
- the provider labels it native/original for that assembly, not projected;
- CDS reconstruction respects strand, phase, translation table, and complete
  codons; and
- the GFF, linkage report, contig map, and audit all have hashes.

A later H1 version, a different haplotype, a GCF annotation whose sequence set
does not match the named GCA, lift-over, projection, or locally transferred
features fails. It is not repaired within this plan.

### 4.4 Four independent eligibility states

| state | minimum acceptance rule | supported claim | common rejection |
|---|---|---|---|
| assembly/composition | valid TaxId; immutable H1; exact original H1-native annotation; valid genetic code/CDS reconstruction; positive original/queryable CDS and gene denominators | GC3 and annotation composition on that exact H1 | annotation mismatch, projection, zero queryable denominator |
| assembly-diversity | composition rule plus resolved same-individual H1/H2 roles; comparable span/collapse audit; whole-haplotype SweepGA 1:1; IMPG calls; positive mapping/callable denominator and explicit exclusions | diversity between two deposited haplotype assemblies, not population diversity | label-only pair, alternate from different individual, zero/unknown denominator |
| population-genomic | exact reference-linked multi-individual genotype/call set; sample/population metadata; filters and matched callability; adequate sample count defined by the downstream model | population-level diversity for the represented population | one assembly pair, missing mask, unmatched VCF reference, pooled/unknown samples |
| demographic | method-specific callable diploid or multi-sample inputs, masks, phasing where required, population labels, and reviewed mutation/generation calibration | only the prespecified demographic method's estimand | assuming H1/H2 is PSMC/MSMC2/SMC++ input, circular `pi/(4mu)` predictor |

The population and demographic flags cannot become `yes` merely because an
assembly or composition flag is `yes`.

### 4.5 License and reuse evidence

Every row stores both a normalized summary and the exact policy or record URL.
NCBI states that it places no restrictions on molecular-data use but also
notes that submitter or country-of-origin rights may exist and are not
transferred to NCBI
([NCBI policies](https://www.ncbi.nlm.nih.gov/home/about/policies/), inspected
2026-07-17). Consequently `NCBI public` is not silently converted into a
blanket license. Unresolved repository or submitter terms block staging (B05).

## 5. Pinned software lineage and responsibilities

### 5.1 Frozen environment facts

| item | frozen identity |
|---|---|
| Guix channel | `https://git.savannah.gnu.org/git/guix.git` at `44bbfc24e4bcc48d0e3343cd3d83452721af8c36` |
| channel file | `analysis/guix/channels.scm`, SHA-256 `45c055cd1d9010a72eacbb720037a22bccb2d8d6891dbd11b5d66365f29b3a17` |
| general manifest | `analysis/guix/manifest.scm`, SHA-256 `2fb05e87aa2ac45ce51d4dcf93b232cb98627f525adace98357629ee3f15720a`, profile `/gnu/store/z9v2f6faha9cwjz0sm5iphhlzisgi077-profile`, derivation `/gnu/store/x7vw3qvf5yxsff7x5cjhxs713m90ni6n-profile.drv` |
| SweepGA build manifest | `analysis/guix/sweepga_origin_main_manifest.scm`, SHA-256 `ea9ae1ba3e51ac3302d93add158532befec8fb3c09d188f524ac29237bab17d1`, profile `/gnu/store/yfffyhdm3a9bsah4gzw9dzri623af3f6-profile`, derivation `/gnu/store/p0clwjvxf0rx09b11gz7dzlhr8mcrc02-profile.drv` |
| handoff runtime manifest | `analysis/guix/sweepga_impg_smoke_manifest.scm`, SHA-256 `c0ef5afd6c988341da8446ff3f70af274dd12f5514bc053d3d4e6f0cbdcee521`, profile `/gnu/store/8x4hx7d9hnv187yprjrzqyg0kxj2z32k-profile`, derivation `/gnu/store/zfl2zzf4i1q6i3m1apmhrggqfw0ijgpr-profile.drv` |
| SweepGA | commit `018e4ce49d2c125820e0ac50dc5feaa02d423683`, version 0.1.1, accepted binary SHA-256 `fa7f0edb9b7e275c288db254046020e136d4267dd5ee043379227ef80da0573b` |
| WFMASH backend | commit `e040aa10e87cab44ed5a4db005e784be62b0bd21`, store `/gnu/store/w9x6axr2w0hhvzzm10gzlp06jg07806d-wfmash-tier3-0.24.2-12.e040aa1` |
| IMPG | commit `101df81eb28a809c8fac97d297acd9fcfbbfa048`, version 0.4.1, binary SHA-256 `c587dc2326cd24f887b1fcb3938404229ad0f0a27ef0773e90c287b1ade160d4` |
| IMPG submodules | `gfaffix` `460e0dd798a9da7d12aef4f9181419d71489da95`; `syng` runtime checkout `dd00f52b688c0fb78cb7f25336ef9ac9f6a3e109` |
| bcftools | 1.14, store `/gnu/store/rlb2gljax8lzmhhidbvbzp3al1ad1mww-bcftools-1.14` |

The Guix manual describes `guix time-machine` as selecting software from a
specific channel commit for reproducible computation
([GNU Guix manual](https://guix.gnu.org/manual/devel/en/guix.pdf), inspected
2026-07-17). Production uses exactly the committed authenticated channel; it
does not append default channels.

### 5.2 Required realization/build capture

On the authorized Guix host only, the environment task will run the following
shape of commands, with absolute output paths inside a new immutable run
directory:

```bash
guix time-machine -C analysis/guix/channels.scm -- \
  package -L analysis/guix -p RUN/profile \
  -m analysis/guix/sweepga_impg_smoke_manifest.scm --no-grafts
guix gc --derivers "$(readlink -f RUN/profile)" > RUN/profile-deriver.txt
guix time-machine -C analysis/guix/channels.scm -- \
  build -L analysis/guix -m analysis/guix/sweepga_impg_smoke_manifest.scm \
  --no-grafts --derivations > RUN/derivations.txt
guix gc --requisites "$(readlink -f RUN/profile)" | sort -u > RUN/store-paths.txt
guix time-machine -C analysis/guix/channels.scm -- \
  describe --format=channels > RUN/resolved-channels.scm
```

SweepGA is rebuilt twice from a clean archive with
`analysis/sweepga_origin_main_rebuild.sh`; the script hash is
`6138d89a96120aa39272893e64636749c5612ee743c1c81c51f8eac401c5d752`.
It runs `cargo fetch --locked`, `cargo build --release --locked --bin
sweepga`, and deterministic stripping. Production requires the two SweepGA
binaries to be byte-identical and to match the accepted hash. The combined
packet must analogously record IMPG source/submodule archives, Cargo lock,
build logs, binaries, `ldd`, store closure, profile derivation, commands, and
truth smoke. Existing evidence is sufficient to define the target but B08
remains until that packet is reviewed for scale-out.

Each job starts through a hash-verified pure wrapper equivalent to
`analysis/slurm/guix_job.sh` (SHA-256
`e0175953ad57b43c0ff95d352ef3a94c8d0bbf520cf70b752f3ee2c669dbbe12`).
It records `env -0`, profile/store paths, executable realpaths and SHA-256,
versions, command argv as JSON, host/kernel/CPU, Slurm IDs, and UTC times.
Ambient Conda, modules, containers, PATH tools, or a rebuilt “equivalent”
binary are prohibited. A tool or channel change requires a separate plan
change, two-build reproduction, truth smoke, scientific equivalence review,
and a new lineage ID.

### 5.3 Responsibility boundary

```text
whole H1 FASTA + whole H2 FASTA
  -> SweepGA/WFMASH, --num-mappings 1:1, --scaffold-jump 0 -> bounded PAF
  -> IMPG index -> implicit graph index
  -> IMPG partition -> native partitions.bed
  -> original H1-native annotation intersection -> selected native focus BED
  -> IMPG query vcf:poa -> regional VCF
  -> IMPG lace -> H1 source-coordinate VCF
  -> bcftools norm/view/dedup -> VCF.gz+TBI and BCF+CSI
  -> original H1-native coding/CDS/fourfold targets + callable intersection
  -> diversity/composition summary and uncertainty
```

SweepGA is not an annotation partitioner. IMPG is not a whole-haplotype
mapping substitute and does not natively emit BCF. The annotation target
builder does not create mappings. These responsibilities are testable file
contracts, not interchangeable implementation suggestions.

## 6. Workflow stages and file contracts

Each run has `inventory_release`, `lineage_id`, and `run_id`; consumers open
only paths listed in the signed/frozen manifest for those IDs.

| contract | producer -> consumer | required content and validation |
|---|---|---|
| O01 source catalog | authorized source task -> parser | URL, immutable revision, SHA-256, bytes, header hash, line counts, retrieval UTC |
| O02 eligibility ledger | parser/resolver -> reviewers | one source row plus species reconciliation; all fields from candidate schema; reason codes; closed-world counts |
| O04 pair evidence | metadata resolver -> preflight | versioned H1/H2, BioSample/individual, reciprocal links, roles, release dates, evidence URLs/hashes |
| O05 asset manifest | staging -> all compute | logical role, accession.version, source URL, provider checksum, byte and sequence-set SHA-256, size, read-only path, license evidence |
| O07 annotation packet | annotation audit -> composition/mapping | native linkage, GFF hash, contig map, genetic code, target BEDs, CDS audit, denominators |
| O08 environment record | Guix preparation -> all jobs | channels/manifests hashes, all `.drv` paths, store closure, GC root, executable hashes, build/smoke packet |
| O09 bounded mapping | SweepGA -> IMPG | standard PAF with `cg:Z`, whole-input hashes, exact argv, 1:1 recheck and coverage/multiplicity QC |
| O10 variant packet | IMPG/bcftools -> summary | graph index, partitions, focus BED, regional and laced VCF, normalized VCF/TBI and BCF/CSI, REF audit |
| O11 denominator packet | native annotation + mapping/calls -> summary | original/queryable/callable genes and bases for coding, CDS, S/W fourfold classes; all exclusions |
| O12 row preflight | validator -> arrays | immutable row fingerprint, accepted modalities, size/resource prediction, zero unresolved required fields |
| O13 pilot review | pilot collector -> authorization | complete outputs/failures, telemetry, prediction errors, QC metrics, I/O saturation, decision |
| O14 run ledger | each array + collector -> release | status, attempt, dependency, checksums, telemetry, error class, sentinel fingerprint, reviewer disposition |

Large objects have `.partial.<jobid>` names in node-local scratch. Promotion
requires successful format validation, nonzero expected records/denominators,
SHA-256 calculation, `fsync` of file and directory, and atomic rename into a
new run-scoped directory. Cross-filesystem copy promotion first writes a
temporary file on the destination filesystem, verifies its hash, then renames.
Nothing writes over an accepted object.

A success sentinel is JSON—not an empty `touch` file—and contains schema
version, row fingerprint, every input/output hash, command hash, environment
record hash, attempt, Slurm IDs, start/end, exit status, and validation
results. Resume accepts the sentinel only if it validates against the current
frozen manifest. Otherwise the element is quarantined and recomputed into a
new attempt directory.

Legacy paths, filenames discovered by wildcard, and prior run directories are
never inputs. `superseded.tsv` records old hash, reason, replacement hash, and
`consumable=no`. Collectors fail if a listed artifact is missing, if an
unlisted artifact is offered, or if any `.partial` path exists in the selected
lineage.

## 7. Resource model

### 7.1 Calibration

The completed whole-pair SweepGA→IMPG runs provide the following usable
calibration:

| dataset | H1 Gbp | H1 contigs | GFF GB | bounded mappings | wall s | CPUs | allocated core-h | staged input dir GB | remap output dir GB | inodes input/output |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| *Menidia menidia* | 0.571 | 457 | 0.443 | 3,353 | 1,082 | 8 | 2.404 | 2.611 | 0.786 | 82 / 1,211 |
| *Spinachia spinachia* | 0.408 | 38 | 0.552 | 753 | 280 | 8 | 0.622 | 1.975 | 0.335 | 99 / 1,174 |
| *Tautogolabrus adspersus* | 0.724 | 64 | 0.313 | 1,972 | 386 | 8 | 0.858 | 2.550 | 0.296 | 84 / 1,212 |

Core-hours above are allocated CPU × wall time, not measured CPU utilization.
The jobs requested 64 GiB but `sacct` recorded blank MaxRSS and zero TotalCPU.
The 135-row Tier 3C primary batch did record process RSS and wall time: 32,801
total seconds at 2 CPUs (18.223 allocated core-hours), 177-second median,
957-second maximum, 3.398-GiB mean RSS, and 15.631-GiB maximum RSS. This batch
includes non-vertebrates and is a stage calibration, not a taxonomic sample.

Actual Tier 3A peak/resident memory, local scratch peak, temporary
amplification, MooseFS read/write bytes, metadata operations, and aggregate
bandwidth were not captured. The resource TSV marks them `UNAVAILABLE`; the
pilot resolves B07/B09.

### 7.2 Formula and covariates

For Tier 3A at eight CPUs, the provisional base wall-time model is:

```text
W_hours = k * [0.033 + 0.083*G + 0.067*max(0, M-1)] * F_C * F_A * F_Q
core_hours = 8 * W_hours
```

`G` is H1 Gbp; `M` is bounded PAF records in thousands; `F_C` is a
contiguity/fragmentation factor; `F_A` is annotation/query target size; and
`F_Q` is inverse callable/queryable fraction pressure. `M`, `F_A`, and `F_Q`
are measured in pilot preflight or replaced by their stated scenario, never
imputed without a flag. The three points cannot support a stable multivariable
fit; coefficients are scheduling bounds. Low/base/high use increasing `k` and
factors as specified in the resource TSV. The pilot refits a robust log-linear
model with genome bases, contig count/N50, GFF bytes/CDS rows, mapping records,
and callable fraction; it reports leave-one-out error and retains an upper
prediction interval for scheduling.

Storage formulas are:

```text
source_bytes = uncompressed H1 + uncompressed H2 + uncompressed GFF
persistent_input = source_bytes * 1.6 base (compressed originals, indexes, provenance)
local_scratch_per_element = source_bytes * 3 base
persistent_output = measured stage files, budgeted 0.5 GB/species base
temporary_amplification = maximum live temporary bytes / final persistent bytes
```

Every stage records `du --bytes`, file/inode counts before/after, local scratch
high-water samples, and process I/O. MooseFS client/server counters are sampled
at stage boundaries where available. `sacct` requests `ElapsedRaw`, `TotalCPU`,
`MaxRSS`, `MaxDiskRead`, `MaxDiskWrite`, `AllocCPUS`, `ReqMem`, state and exit
code; Slurm documents that these fields depend on the configured accounting
plugin ([Slurm sacct](https://slurm.schedmd.com/sacct.html), inspected
2026-07-17). A sidecar cgroup/getrusage record is mandatory when any field is
blank.

### 7.3 Totals, concurrency, and retention

The base reported-catalog envelope is 116.2 allocated core-hours, 400 GB
persistent input, 25 GB persistent output, 111 GB peak local scratch, about
100,000 inodes, 370 GB MooseFS reads, 190 GB writes, 970,000 metadata
operations, and 120 MiB/s aggregate bandwidth. With 25% quota headroom, request
500 GB persistent and 140 GB scratch only after the pilot confirms the base
case. The high envelope is not an allocation request.

Initial pilot caps are two Tier 3A tasks, four metadata requests, two concurrent
downloads, and eight Tier 3C tasks. Reviewed waves may raise Tier 3A to 10 and
Tier 3C to 16 only if memory and I/O gates pass. Ten Tier 3A tasks at 64 GiB
require 640 GiB aggregate requested memory. A token bucket limits staging to
two remote transfers and 100 MiB/s; local promotion plus reads are capped at
120 MiB/s five-minute average and 200 metadata operations/s. At 80% of either
cap for ten minutes, launch pauses; at 95%, running tasks finish but the wave
does not resume without the HPC storage reviewer.

Retention classes are:

- immutable source catalog, asset checksums/metadata, Guix closure records,
  final BCF/CSI, target/denominator packets, summaries, sentinels, failures,
  and telemetry: retain through publication plus seven years, subject to
  institutional policy;
- compressed H1/H2/GFF inputs: retain through final reproducibility review;
  after review they may be evicted only if accession, URL, provider checksum,
  repository SHA-256, and a tested rehydration recipe remain;
- PAF, graph indexes, regional VCFs, and other reproducible intermediates:
  retain through independent QC and 90 days after release, then delete only
  via an approved checksum-listed cleanup task;
- node-local scratch and `.partial` files: remove after verified promotion or
  after 14-day failure triage. Cleanup never follows unresolved globs and
  always writes a deletion ledger.

## 8. Resumable Slurm and I/O design

The array topology is:

```text
freeze inventory
  -> environment smoke
  -> stage assets[%row] --array ...%2
  -> exact preflight[%row] --dependency=aftercorr:stage
  -> SweepGA map[%row] --dependency=aftercorr:preflight --array ...%2 pilot
  -> IMPG/summary[%row] --dependency=aftercorr:map
  -> collect --dependency=afterany:<all arrays>
  -> independent QC
  -> human wave/full authorization
```

`aftercorr` couples corresponding elements; collection uses `afterany` so a
single failure cannot hide the ledger. Slurm documents `aftercorr` and notes
that array-wide `afterok` depends on every element's success
([Slurm arrays](https://slurm.schedmd.com/job_array.html), inspected
2026-07-17). Array-wide release/QC uses `afterok` only after the collector has
materialized a zero-unexplained-failure manifest.

Inputs are staged once, checksum-verified, permissions changed read-only, and
copied to `$SLURM_TMPDIR/<run>/<candidate>` before compute. Compute does not
stream whole genomes from MooseFS. Outputs remain local until complete, then
are promoted in batches to reduce metadata load. Logs and telemetry are small
and append only within an attempt directory.

Retry classes are explicit:

- transient network, node loss, or filesystem timeout: exponential backoff,
  at most two automated retries;
- OOM/time limit: one reviewed retry at at most 2× memory/time, placed in the
  outlier lane and fed back to the budget;
- checksum, exact-reference, annotation, denominator, sample-identity, or
  scientific QC failure: no automated retry and no alternate asset; reject or
  open a new inventory amendment;
- code/command/environment mismatch: quarantine the whole affected lineage,
  not merely the failed element.

The collector emits `retry_manifest.tsv` containing only failed IDs and their
frozen row fingerprints. Recovery submits that manifest as a new array with
new attempt IDs. Successful elements are never rerun unless an explicit
supersession task invalidates their lineage.

## 9. Scientific output and estimands

### 9.1 Species-result table

The final long table has one row per `(candidate_id, sampling_unit_id,
target_class, estimand)` and the following groups of fields:

| group | required fields |
|---|---|
| identity | inventory release, candidate/run/lineage IDs, source and current scientific name, TaxId, lineage ranks, clade, H1/H2 accessions, BioSample, individual, sex, population/geography where known |
| sampling | data modality, sampling unit (`one phased diploid assembly pair`, `one individual`, or population sample), sample/haplotype/chromosome counts, ascertainment |
| targets | original/queryable/callable gene counts and bases for coding gene, CDS, fourfold reference-S and reference-W; GC3 numerator/denominator; genetic code |
| diversity | variant numerator, callable denominator, coding π, CDS π, reference-conditioned π_S, π_W, π_S/π_W, block-bootstrap SE/CI, bootstrap unit/count |
| mapping/QC | H1/H2 lengths/contigs, mapped bases/fraction, bounded PAF records, observed query/target multiplicity, excluded partial/repeat/ambiguous bases, REF mismatches, phase/frame discordance |
| missingness | eligibility state, reason code, stage, whether structural/resource/QC missing, no zero substitution |
| phylogeny | tree ID/hash, tip ID, reconciliation status, branch length, order/class, clade-balance stratum/weight |
| provenance | input/output hashes, annotation release/linkage, tool lineage, commands, attempt, telemetry and review IDs |

π is the number of observed H1/H2 alternative alleles divided by callable H1
bases for that target class. These assembly-pair values are explicitly named
`diploid_haplotype_assembly_diversity`; they are not labeled population
heterozygosity. π_S and π_W are reference-conditioned fourfold classes and the
ratio retains paired-block uncertainty where possible. GC3 is computed only
from successfully reconstructed original H1-native CDS and reports both genes
and codons in its denominator.

### 9.2 Prespecified analysis

Primary descriptive estimands are:

1. the clade-balanced median and interquartile range of log10 CDS π among
   assembly-diversity-eligible species; and
2. the phylogenetically adjusted association between log10 CDS π and logit
   GC3 among species passing both workflows, reported as a slope with 95%
   interval under a Brownian PGLS covariance.

The first is a distributional summary, not a test of effective population
size. The second is a descriptive cross-species association and does not imply
that GC3 causes diversity or vice versa.

Species is the primary independent unit. If multiple individuals later exist,
the confirmatory analysis first produces one prespecified species estimate
using inverse-variance weighting with a species-level robust variance; it does
not count individuals as independent species. One representative per species
and inverse eligible-count weights by order form the clade-balanced analysis;
unweighted and one-per-family results are sensitivities.

The primary phylogeny is frozen before outcome modeling, reconciled by TaxId,
and hashed (B12). Brownian PGLS is primary; Pagel-λ, OU, phylogenetic mixed
model with clade random intercepts, and order-cluster robust regression are
sensitivity models. Results are reported only when at least 20 species, at
least four classes, and at least eight orders contribute; otherwise analysis
is descriptive tables/plots only. No tree tip is guessed from string
similarity; unresolved tips are excluded with a reason.

Missing data are never mean-imputed. Exact-reference/annotation failures are
structural exclusion; zero callable bases are undefined, not zero diversity;
resource failures remain pending until retry resolution; missing covariates
reduce only models requiring them. The report compares eligibility and
missingness by class/order, genome size, contiguity and annotation size. A
selection model or inverse-probability weighting is sensitivity-only because
the availability mechanism is unlikely to be ignorable.

Secondary analyses include coding π, π_S, π_W, π_S/π_W, GC3, mapping
fraction/multiplicity, callable fraction, and interactions prespecified by
major class. The two primary estimands are labeled primary; all secondary
coefficient tests are corrected together with Benjamini–Hochberg FDR at 0.05,
and raw and adjusted values are both shown. Post hoc clade contrasts are
exploratory. Sensitivities include callable-fraction thresholds 0.25/0.50/0.75,
mapping-coverage thresholds, exclusion of high-fragmentation/high-multiplicity
rows, leave-one-order-out fits, alternative tree/covariance models, ratio
analysis on log π_S − log π_W, and bounds using bootstrap endpoints.

No causal vocabulary is permitted without an external identification design.
Assembly selection, annotation availability, genome architecture, mutation,
recombination, life history, and shared ancestry are plausible confounders.
Population-genomic and demographic analyses use separate eligible subsets and
must not promote `pi/(4*mu)` or another response-derived Ne as an independent
predictor.

## 10. Pilot design

After B01–B05 are resolved, choose eight rows deterministically from the
accepted metadata inventory, before inspecting diversity outcomes:

- one cartilaginous fish and two ray-finned fish spanning the low/high genome
  size and contiguity cells among the 13 named fish;
- one amphibian, one reptile, one bird, one mammal, and one reported “other”
  vertebrate from the 40-row cohort;
- across the eight, cover genome-size tertiles, contig-count tertiles,
  annotation-size tertiles, and expected low/high mapping complexity; and
- at least six must be pair-eligible for Tier 3A and all eight must be
  composition-eligible. If the resolved inventory cannot fill these slots, the
  pilot is redesigned and reviewed; an ineligible row is not substituted from
  another reference.

This design is provisional because the non-fish rows are unavailable (B10).
Selection code, seed, strata and tie-break by TaxId are frozen in O12. The
completed *Spinachia* truth case is rerun as a technical positive control but
is not counted among the eight scale-out species or scientific results.

The pilot records five-second RSS/CPU/I/O/scratch samples, `sacct`, file and
inode deltas, MooseFS client/server bytes and metadata RPCs where available,
stage wall times, input covariates, mapping counts, callable fractions, and
all failures. It is the calibration experiment for B07/B09, not a miniature
publication analysis.

## 11. Stages, quantitative gates, and reviewers

Reviewer names below are WG role names. B11 requires an individual human
assignee to be recorded before the corresponding task becomes ready.

| stage | go criteria | stop/rollback criteria | named reviewer / authorization boundary |
|---|---|---|---|
| 0 guidance/schema freeze | source 717/716 lines; all six column assertions; counts 714/223/248/271/120/40 and taxonomic breakdown reproduced; 13/46 subsets reproduced; 100% rows closed-world; source/license provenance complete | any unexplained count/schema/source mismatch; moving revision; placeholder date | `Inventory Curator`; `authorize-catalog-metadata-retrieval` permits only the small catalog object |
| 1 metadata/preflight | 100% of proposed pilot rows have TaxId, exact versioned H1/H2 where needed, pair evidence, exact annotation candidate, license evidence, size estimate, and explicit modality status; zero unresolved required fields | any silent substitute, projected annotation, unresolved terms, or missing row | `Tier3 Scientific Reviewer` + `Data Provenance Reviewer`; `authorize-pilot-asset-download` sets bytes/accessions |
| 2 environment/dry run | every channel/manifest/derivation/closure hash captured; two-build SweepGA identity; IMPG/submodule identity; controlled fixture and H1-native truth case 100%; zero ambient executables | any binary/store mismatch, missing derivation, truth discrepancy, or compute-node profile failure | `Reproducibility Reviewer`; no biological compute until signed O08 |
| 3 eight-species pilot | 100% exact-reference/native-annotation gates; positive denominators; 100% emitted REF alleles match H1; ≥50% queryable target bases base threshold or reasoned sensitivity stratum; first-attempt technical failure ≤1/8 and zero unexplained failures; median core-hour/wall/storage prediction error ≤35%, 95th percentile ≤75%; RSS/scratch <80% request; metadata and bandwidth <80% caps; complete telemetry 100% | any reference/annotation contamination; >1/8 technical failures; zero denominator; missing telemetry; sustained ≥95% I/O cap | `Tier3 Scientific Reviewer`, `HPC Storage Reviewer`, `Reproducibility Reviewer`; `review-pilot-go-no-go` |
| 4 reviewed wave 1 | up to 8 Tier3A and 24 Tier3C; accepted rows 100% gates; technical failure ≤5%, unexplained ≤1%; closed-world completeness 100%; median resource error ≤25%, p95 ≤50%; no OOM; I/O <80% caps | threshold breach, quota > approved envelope, systemic error in ≥2 rows | same three reviewers; `authorize-expansion-wave-1` records row IDs and budget |
| 5 later waves | each wave ≤10 Tier3A and ≤32 Tier3C; cumulative technical failure ≤5%; unexpected scientific failure classified; outputs/checksums/telemetry ≥99% and 100% ledgers; no MFS incident; prediction error remains within wave-1 bounds | two similar unexplained failures, any lineage mismatch, >80% quota after projected next wave | reviewers issue a separate `authorize-expansion-wave-N`; no blanket approval |
| 6 full catalog decision | all prior waves reviewed; 100% candidates resolved; exact/native/denominator violations zero; final expected resources within approved quota +25% headroom; retry backlog zero; independent QC concordance 100% on hashes and ≥99.9% record-level calls; statistical/phylogeny plan frozen | any unresolved blocker, storage/I/O saturation, resource high bound needed, or reviewer dissent | **`authorize-full-vertebrate-execution`**, jointly approved by human project owner, Scientific Reviewer, HPC Storage Reviewer, and Reproducibility Reviewer |
| 7 final release | accepted+rejected+unavailable+failed equals inventory; 100% provenance; no `.partial`/legacy/superseded input; primary table reproduces; multiple-testing and sensitivity reports complete | count mismatch, untracked artifact, result changes under same hashes | `Independent Results Reviewer` and human project owner |

Quality/completeness, resource prediction error, failure rate, and I/O
saturation therefore have explicit thresholds at every transition. A gate may
lower concurrency or narrow a wave. It may not lower exact-reference,
native-annotation, checksum, denominator, or provenance requirements.

## 12. Primary repositories and retrieval provenance

These official sources define later metadata/acquisition procedures. They were
inspected on 2026-07-17; no full catalog or assembly package was downloaded by
this planning task.

| source | authoritative use | URL | captured uncertainty |
|---|---|---|---|
| VGP Phase 1 repository | freeze TSV location and displayed 717-line/320-KB object | https://github.com/VGP/vgp-phase1/blob/main/VGPPhase1-freeze-1.0.tsv | moving `main`, commit and content SHA unresolved (B01) |
| NCBI Datasets genome docs | accession-specific metadata/packages, include/dehydrated options | https://www.ncbi.nlm.nih.gov/datasets/docs/v2/how-tos/genomes/download-genome/ | exact CLI version must be frozen in Guix |
| NCBI assembly model | assembly accession/version and GCA/GCF semantics | https://www.ncbi.nlm.nih.gov/datasets/docs/v2/data-processing/policies-annotation/data-model/ | GCA/GCF relationship still requires exact sequence-set proof per row |
| NCBI sequence report | component accession and length dictionary | https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/data-reports/genome-sequence/ | report is metadata, not a content checksum |
| NCBI package validation | built-in validation and package `md5sum.txt` | https://www.ncbi.nlm.nih.gov/datasets/docs/v2/how-tos/validation/ | MD5 is provider transport evidence; repository SHA-256 is also required |
| NCBI Genomes FTP | old annotation releases and per-directory checksums | https://www.ncbi.nlm.nih.gov/datasets/docs/v2/data-processing/policies-annotation/genomeftp/ | availability varies by release and must be recorded |
| NCBI Taxonomy | TaxId/current name/lineage | https://www.ncbi.nlm.nih.gov/datasets/docs/v2/data-processing/taxonomy-processing/taxonomy/ | merged/synonym IDs retained, never silently rewritten |
| NCBI policies | molecular-data reuse statement | https://www.ncbi.nlm.nih.gov/home/about/policies/ | submitter rights may remain; row-level review required |
| GNU Guix manual | time-machine and channel reproducibility semantics | https://guix.gnu.org/manual/devel/en/guix.pdf | old channels do not receive later security fixes; production network exposure minimized |
| Slurm official docs | arrays, dependency semantics, accounting fields | https://slurm.schedmd.com/job_array.html and https://slurm.schedmd.com/sacct.html | metric availability depends on cluster plugins (B09) |

## 13. Authorization conclusion

The program may proceed now only to review these planning artifacts. The next
possible external action is the small, separately approved acquisition of the
320-KB-class freeze TSV at an immutable revision. No pilot asset, bulk genome
download, Slurm array, expansion wave, full-catalog execution, or demographic
inference is authorized by this plan. Generated future scripts remain inert
until the exact authorization node in the proposed WG graph is completed with
an explicit scope and resource envelope.
