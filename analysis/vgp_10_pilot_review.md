# Independent review of the ten-genome VGP pilot

Review ID: `vgp10-review-20260718-v1`
Run reviewed: `vgp10-20260718-preflight-v1`
Program decision: **CONDITIONAL_GO — bounded repair and re-pilot only; full biological scale-out is not authorized**

## Executive finding

All ten immutable primary slots are accounted for, retain the exact frozen H1/H2 accession.version and BioSample/individual identity, and pass live content-addressed checksum review. All ten nevertheless failed before mapping because mandatory assembly-QC evidence was absent: exact-final QV for both haplotypes (10/10), H2 BUSCO (10/10), a manifest-bound k-mer/copy-number audit (10/10), and a repeat/low-complexity mask (10/10). P07/P08 additionally lack H1 BUSCO and P09 lacks resolved exact read chemistry. These are technical packet failures. They are not low diversity, biological outliers, mapping failures, or failed PSMC histories.

Zero biological jobs were submitted. Therefore whole-assembly 1:1 mapping, multiplicity, IMPG extraction, normalized variants, REF/alternate reconstruction, callable masks, consensus, diversity, unscaled PSMC, scaling scenarios, and bootstraps are **not reached**, not passing zeros. With 0/10 callable+consensus passers, zero represented clades/generations among passers, systematic failure in both generation/technology strata, and unestimable resource-model APE, a core GO is prohibited. Exact identity/checksum success and the repairable common cause support CONDITIONAL_GO rather than permanent NO_GO.

Annotation state never caused a core failure. P03/P04 annotation mismatches remain branch-local, and absence of annotation or independent Ne evidence cannot veto a future technically valid core result. PSMC, if later produced, is descriptive evidence from the same pair and is not independent validation of same-pair diversity.

## Ten primary slots

| slot | species | exact H1 / H2 | individual | generation; technology and phasing | annotation branch | classification | reasons |
|---|---|---|---|---|---|---|---|
| P01 | *Camelus dromedarius* | `GCA_036321535.1` / `GCA_036321565.1` | `SAMN39296380`; `mCamDro1` | `later_hifi_hic_hap1_hap2`; H1: PacBio Sequel II HiFi; Dovetail OmniC; H2: PacBio Sequel II HiFi; Dovetail OmniC; trio-binning plus Omni-C and YaHS | `available_paired_refseq_pending_dictionary_audit` | technical failure; repairable; no biological classification | `MISSING_EXACT_FINAL_SEQUENCE_QV;MISSING_H2_BUSCO;MISSING_KMER_COPY_NUMBER_AUDIT;MISSING_REPEAT_OR_LOW_COMPLEXITY_MASK` |
| P02 | *Pseudorca crassidens* | `GCA_039906515.1` / `GCA_039906525.1` | `SAMN41253811`; `mPseCra1` | `later_hifi_hic_hap1_hap2`; H1: PacBio Sequel II HiFi; Bionano Genomics DLS; Dovetail Omni-C; H2: PacBio Sequel II HiFi; Bionano Genomics DLS; Dovetail Omni-C; Hi-C phasing plus Omni-C; Bionano DLS; YaHS | `available_paired_refseq_pending_dictionary_audit` | technical failure; repairable; no biological classification | `MISSING_EXACT_FINAL_SEQUENCE_QV;MISSING_H2_BUSCO;MISSING_KMER_COPY_NUMBER_AUDIT;MISSING_REPEAT_OR_LOW_COMPLEXITY_MASK` |
| P03 | *Colius striatus* | `GCA_028858725.2` / `GCA_028858625.2` | `SAMN33339572`; `bColStr4` | `later_hifi_hic_hap1_hap2`; H1: PacBio Sequel II HiFi; Bionano DLS; Arima Hi-C v2; H2: PacBio Sequel II HiFi; Bionano DLS; Arima Hi-C v2; trio phasing plus Arima Hi-C and Bionano; Salsa | `failed_mismatch` | technical failure; repairable; no biological classification | `MISSING_EXACT_FINAL_SEQUENCE_QV;MISSING_H2_BUSCO;MISSING_KMER_COPY_NUMBER_AUDIT;MISSING_REPEAT_OR_LOW_COMPLEXITY_MASK` |
| P04 | *Falco naumanni* | `GCA_017639655.1` / `GCA_017639645.1` | `SAMN16870685`; `bFalNau1` | `early_clr_triocanu_parental`; H1: PacBio Sequel I CLR; Illumina NovaSeq; Arima Genomics Hi-C; Bionano Genomics DLS; H2: PacBio Sequel I CLR; Illumina NovaSeq; Arima Genomics Hi-C; Bionano Genomics DLS; TrioCanu parental phasing plus Arima Hi-C; Bionano; Salsa2 | `failed_mismatch` | technical failure; repairable; no biological classification | `MISSING_EXACT_FINAL_SEQUENCE_QV;MISSING_H2_BUSCO;MISSING_KMER_COPY_NUMBER_AUDIT;MISSING_REPEAT_OR_LOW_COMPLEXITY_MASK` |
| P05 | *Candoia aspera* | `GCA_035149785.1` / `GCA_035125265.1` | `SAMN37159891`; `rCanAsp1` | `later_hifi_hic_hap1_hap2`; H1: PacBio Sequel II HiFi; Arima Hi-C v2; H2: PacBio Sequel II HiFi; Arima Hi-C v2; Hi-C phasing plus Arima Hi-C and YaHS | `available_paired_refseq_pending_dictionary_audit` | technical failure; repairable; no biological classification | `MISSING_EXACT_FINAL_SEQUENCE_QV;MISSING_H2_BUSCO;MISSING_KMER_COPY_NUMBER_AUDIT;MISSING_REPEAT_OR_LOW_COMPLEXITY_MASK` |
| P06 | *Dendropsophus ebraccatus* | `GCA_027789765.1` / `GCA_027789725.1` | `SAMN32145295`; `aDenEbr1` | `early_clr_triocanu_parental`; H1: PacBio Sequel I CLR; 10X Gemonics linked reads; Bionano Genomics DLS; Arima Genomics Hi-C v1; Illumina WGS on parents; H2: PacBio Sequel I CLR; 10X Gemonics linked reads; Bionano Genomics DLS; Arima Genomics Hi-C v1; Ilumina WGS on parents; TrioCanu parental phasing plus 10X; Bionano; Arima Hi-C; Salsa | `available_paired_refseq_pending_dictionary_audit` | technical failure; repairable; no biological classification | `MISSING_EXACT_FINAL_SEQUENCE_QV;MISSING_H2_BUSCO;MISSING_KMER_COPY_NUMBER_AUDIT;MISSING_REPEAT_OR_LOW_COMPLEXITY_MASK` |
| P07 | *Spinachia spinachia* | `GCA_048126635.1` / `GCA_048127205.1` | `SAMN36735485`; `fSpiSpi1` | `later_hifi_hic_hap1_hap2`; H1: HiFi; HiC; H2: HiFi; HiC; Hi-C phasing plus YaHS; Samba | `available_exact_native_pending_dictionary_audit` | technical failure; repairable; no biological classification | `MISSING_EXACT_FINAL_SEQUENCE_QV;MISSING_H1_BUSCO;MISSING_H2_BUSCO;MISSING_KMER_COPY_NUMBER_AUDIT;MISSING_REPEAT_OR_LOW_COMPLEXITY_MASK` |
| P08 | *Menidia menidia* | `GCA_048628825.1` / `GCA_048544195.1` | `SAMN46987722`; `fMenMen1` | `later_hifi_hic_hap1_hap2`; H1: PacBio HiFi; Arima Hi-C; H2: PacBio HiFi; Arima Hi-C; Arima Hi-C phasing plus YaHS | `available_exact_native_pending_dictionary_audit` | technical failure; repairable; no biological classification | `MISSING_EXACT_FINAL_SEQUENCE_QV;MISSING_H1_BUSCO;MISSING_H2_BUSCO;MISSING_KMER_COPY_NUMBER_AUDIT;MISSING_REPEAT_OR_LOW_COMPLEXITY_MASK` |
| P09 | *Heterodontus francisci* | `GCA_036365525.1` / `GCA_036365495.1` | `SAMN39432692`; `sHetFra1` | `later_hifi_hic_hap1_hap2`; H1: PacBio Sequel; H2: PacBio Sequel; Hi-C phasing plus Bionano and YaHS | `available_paired_refseq_pending_dictionary_audit` | technical failure; repairable; no biological classification | `MISSING_EXACT_FINAL_SEQUENCE_QV;MISSING_H2_BUSCO;MISSING_KMER_COPY_NUMBER_AUDIT;UNRESOLVED_EXACT_READ_CHEMISTRY;MISSING_REPEAT_OR_LOW_COMPLEXITY_MASK` |
| P10 | *Hemiscyllium ocellatum* | `GCA_020745735.1` / `GCA_020745765.1` | `SAMN22550098`; `sHemOce1` | `early_clr_triocanu_parental`; H1: PacBio Sequel II CLR; Illumina NovaSeq; Bionano Genomics DLS; Arima Genomics Hi-C v1; H2: PacBio Sequel II CLR; Illumina NovaSeq; Bionano Genomics DLS; Arima Genomics Hi-C v1; TrioCanu parental phasing plus Bionano; Arima Hi-C; Salsa2 | `available_paired_refseq_pending_dictionary_audit` | technical failure; repairable; no biological classification | `MISSING_EXACT_FINAL_SEQUENCE_QV;MISSING_H2_BUSCO;MISSING_KMER_COPY_NUMBER_AUDIT;MISSING_REPEAT_OR_LOW_COMPLEXITY_MASK` |

There are no low-confidence-but-usable core results and no biological outliers because no biological measurement crossed preflight. Long-range phasing evidence is present for every pair and is retained as later confidence evidence; it cannot substitute for QV, completeness, collapse, mapping, or callability measurements.

## Alternates and failure accounting

All six declared alternates remain `standby_not_triggered`; none has an amendment, none replaced a failed primary, and every failed primary remains in the result ledger. All 10 artifact packets and their 90 files pass digest verification. The exact reason totals reconcile across result, QC, pair failure artifacts, independent validation, and run summary: 10 QV, 2 H1 BUSCO, 10 H2 BUSCO, 10 k-mer/copy-number, 10 repeat-mask, and 1 chemistry reasons (43 total reason instances). There are no unknown warning codes, retries, silent drops, scheduler IDs, or Slurm dependency edges.

| alternate | replaces clade | species | exact H1 / H2 | individual | generation | disposition |
|---|---|---|---|---|---|---|
| A01 | Mammalia | *Inia geoffrensis* | `GCA_036417435.1` / `GCA_036417475.1` | `SAMN32797734`; `mIniGeo1` | `later_hifi_omnic_hap1_hap2` | standby; no trigger or versioned amendment; not reviewed as a biological result |
| A02 | Aves | *Lonchura striata domestica* | `GCA_046129695.1` / `GCA_046129705.1` | `SAMN44779081`; `bLonStr1` | `later_hifi_hic_parental` | standby; no trigger or versioned amendment; not reviewed as a biological result |
| A03 | Reptilia | *Anolis sagrei* | `GCA_037176765.1` / `GCA_037176775.1` | `SAMN40144551`; `rAnoSag1` | `later_hifi_hic_parental` | standby; no trigger or versioned amendment; not reviewed as a biological result |
| A04 | Amphibia | *Xenopus petersii* | `GCA_038501925.1` / `GCA_038501915.1` | `SAMN39187339`; `aXenPet1` | `later_hifi_hic_parental` | standby; no trigger or versioned amendment; not reviewed as a biological result |
| A05 | Actinopterygii | *Syngnathus typhle* | `GCA_048301445.1` / `GCA_048301605.1` | `SAMN36735486`; `fSynTyp1` | `later_hifi_hic_hap1_hap2` | standby; no trigger or versioned amendment; not reviewed as a biological result |
| A06 | Chondrichthyes | *Hydrolagus colliei* | `GCA_035084275.1` / `GCA_035084065.1` | `SAMN39156054`; `sHydCol1` | `later_hifi_hic_hap1_hap2` | standby; no trigger or versioned amendment; not reviewed as a biological result |

## Independent stratified recomputation

P07 and P08 are the predeclared full-independent sentinels. In the independent output root, the reviewer rehashed and size-checked all 20 of their H1/H2 core objects, then compared production and independent eligibility for normalized variants, masks and reason totals, callable denominators, diversity, diploid consensus, PSMC input/unscaled result, and bootstraps. Every object matched. Every biological quantity matched only as explicit nonmaterialization after the same immutable preflight failure. A numerical biological recomputation would require bypassing preregistered gates and was correctly not attempted. No Slurm or biological job was submitted.

## Gate application

| gate | requirement | observed | result | consequence |
|---|---|---|---|---|
| H01 pair/individual identity resolved | unresolved=0 | 0 | **PASS** | none |
| H02 exact accession.version identity resolved | unresolved=0 | 0 | **PASS** | none |
| H02A taxid/BioSample/individual/haplotype-role fields resolved | unresolved=0 | 0 | **PASS** | none |
| H03 immutable checksum drift | events=0 | 0 | **PASS** | none |
| H04 retained query multiplicity >1 | count=0 | no retained mapping | **NOT_REACHED** | cannot support GO |
| H05 retained target multiplicity >1 | count=0 | no retained mapping | **NOT_REACHED** | cannot support GO |
| H05A retained non-1:1 bases | bp=0 | no retained mapping | **NOT_REACHED** | cannot support GO |
| H06 unexplained mask accounting | discrepancy_bp=0 | not measured | **NOT_REACHED** | cannot support GO |
| H07 H1 REF reconstruction | failures=0 | not measured | **NOT_REACHED** | cannot support GO |
| H08 H2 alternate reconstruction | failures=0 | not measured | **NOT_REACHED** | cannot support GO |
| H08A non-callable bases encoded homozygous reference | bp=0 | not measured | **NOT_REACHED** | cannot support GO |
| H08B annotation sequence-dictionary mismatch retained | count=0 | no annotation output | **NOT_REACHED** | cannot support annotation GO |
| H09 unscaled/scaled PSMC separation | conflations=0 | no PSMC outputs | **NOT_REACHED** | cannot support GO |
| A01 exact-final-sequence QV each haplotype | QV>=40 each | 0/10 measured | **FAIL** | repair required |
| A02 H1/H2 BUSCO completeness and missingness | complete>=0.90 and missing<=0.05 each | H1 missing 2/10; H2 missing 10/10 | **FAIL** | repair required |
| A02A H1/H2 BUSCO completeness difference | absolute difference<=0.05 | 0/10 pair differences measurable | **FAIL** | repair required |
| A02B H1/H2 BUSCO duplication | duplicated<=0.05 each | 0/10 both-haplotype audits measurable | **FAIL** | repair required |
| A03 copy-number/k-mer collapse audit | passing both haplotypes | 0/10 measured | **FAIL** | repair required |
| A04 repeat/low-complexity mask | exact manifest-bound mask | 0/10 measured | **FAIL** | repair required |
| A05 exact read chemistry | resolved | 9/10 resolved | **FAIL** | P09-specific repair |
| A06 minimum haplotype span | >=250,000,000 bp each | 20/20 pass | **PASS** | does not override missing QC |
| A07 minimum contig N50 | >=1,000,000 bp each | 20/20 pass | **PASS** | does not override missing QC |
| A08 H1/H2 length ratio | 0.80..1.25 | 10/10 pass; 0.941..1.157 | **PASS** | does not override missing QC |
| C01 primary callable fraction | >=0.60 | 0/10 measured | **NOT_REACHED** | cannot support GO |
| C01A sensitivity callable fraction | >=0.50 | 0/10 measured | **NOT_REACHED** | cannot support GO |
| C01B minimum callable bases | >=100,000,000 bp | 0/10 measured | **NOT_REACHED** | cannot support GO |
| C01C well-callable windows | >=50 1-Mb windows at >=80% callable | 0/10 measured | **NOT_REACHED** | cannot support GO |
| C01D ordered disjoint reason assignment | 13 predeclared reasons; exact complement | 0/10 measured | **NOT_REACHED** | cannot support GO |
| C01E universe reconciliation | callable + reasons = universe exactly | 0/10 measured | **NOT_REACHED** | cannot support GO |
| C02 masked and heterozygous encoding | masked=N; heterozygous SNP=IUPAC | 0/10 produced | **NOT_REACHED** | cannot support GO |
| C03 primary indel flank mask | 10 bp | 0/10 produced | **NOT_REACHED** | cannot support GO |
| C04 indel flank sensitivities | 0 bp and 50 bp | 0/10 produced | **NOT_REACHED** | cannot support GO |
| C05 exact-reference normalization and duplicate removal | exact H1; exact duplicates removed | 0/10 produced | **NOT_REACHED** | cannot support GO |
| B00 pre-registered primary bootstrap attempts | exactly 200 | no passing PSMC unit | **NOT_REACHED** | cannot support GO |
| B01 minimum bootstrap attempts | >=100 | no passing PSMC unit | **NOT_REACHED** | cannot support GO |
| B02 finite bootstrap fraction | >=0.95 | no passing PSMC unit | **NOT_REACHED** | cannot support GO |
| B02A finite primary bootstrap successes | >=190/200 | no passing PSMC unit | **NOT_REACHED** | cannot support GO |
| B03 boundary-aware block construction | 5-Mb primary; 1/10-Mb sensitivities; never cross contig/mask boundary | no bootstrap units | **NOT_REACHED** | cannot support GO |
| B04 heterozygosity bootstrap replicates | 10,000 | 0/10 attempted | **NOT_REACHED** | cannot support GO |
| D01 P07/P08 full independent biological recomputation | 2 sentinels; digest mismatches=0 | 20 inputs rehashed; biological outputs ineligible | **NOT_REACHED** | re-pilot must complete biological comparison |
| D02 other-pair deterministic shard comparison | 1 shard each; digest mismatches=0 | 0/8 biological shards eligible | **NOT_REACHED** | re-pilot must complete biological comparison |
| P00 primary slots adjudicated | exactly 10 | 10 | **PASS** | none |
| P01 primary callable+consensus passes | >=8/10 | 0/10 | **FAIL** | core GO prohibited |
| P02 required major clades among passers | 6/6 | 0/6 | **FAIL** | core GO prohibited |
| P03 assembly generations among passers | 2/2 | 0/2 | **FAIL** | core GO prohibited |
| P04 no technology-stratum systematic failure | failure fraction <0.50 | 1.00 in both early CLR and later HiFi strata | **FAIL** | core GO prohibited |
| P05 zero hard-gate violations | 0 | measured identity/checksum violations=0; downstream hard gates not reached | **NOT_REACHED** | core GO prohibited |
| R01 median absolute percentage error | <=0.25 | not estimable | **FAIL_NOT_ESTIMABLE** | core GO prohibited |
| R02 95th percentile absolute percentage error | <=0.50 | not estimable | **FAIL_NOT_ESTIMABLE** | core GO prohibited |
| R03 storage headroom | >=0.25 | 0.25 retained in corrected envelope | **PASS** | planning only; does not rescue APE failure |
| R04 per-job stop multiple | 1.5x reviewed high estimate | 1.5x retained; no job invoked | **PASS** | planning only; not a global ceiling |
| L01 all primary slots reconciled | 10/10 | 10/10 | **PASS** | none |
| L02 all alternates reconciled | 6/6 | 6/6 standby; 0 activated | **PASS** | none |
| L03 failure/warning reason totals | exact | 43 reason instances across 6 codes | **PASS** | none |
| N01 annotation absence cannot veto core | zero annotation-driven core failures | 0 | **PASS** | none |
| N02 same-pair PSMC is not independent validation | no independent claim | no PSMC result and no claim | **PASS** | none |

The bootstrap 100-attempt/95%-finite rules were not violated by a passing PSMC unit—there is no passing unit—but they were also not demonstrated. The program requires positive evidence, so vacuous truth is not used to support GO.

## Resource review and upper-bound scale sensitivity

Observed preflight revalidated 100 objects and 12,567,760,437 logical bytes using 9.403 CPU-seconds, 10.112 elapsed seconds, 40.7 MiB peak RSS, zero scratch, 466,944 filesystem-read bytes, 20,793 logical report bytes, zero filesystem-write bytes, and 300 metadata operations. A transparent 716-pair upper-bound extrapolation for this checksum-only step is 7,160 objects, 899,851,647,289 logical bytes, 0.187 core-hours, and 0.201 serial wall-hours.

The pilot contains no observed mapping, IMPG, consensus, PSMC/bootstrap, scratch, or cluster-I/O telemetry, so a fitted full biological resource model—and hence median/p95 APE—does not exist. The corrected resource manifest retains a low/base/high 716-pair sensitivity by scaling the preregistered 40-pair Tier3A envelope (which was informed by three earlier calibration tuples), with eligibility explicitly unresolved. It is a planning lower bound because PSMC plus 200 bootstraps were absent from that older envelope:

| scenario | minimum durable objects | persistent input/output GB | operational inodes | core-hours lower bound | memory/job GiB | scratch/job; aggregate GB | read/write GB | concurrency | wall-hours lower bound |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| low | 29,356 | 1,145.6 / 107.4 | 537,000 | 322.2 | 32 | 3 / 12 | 895.0 / 1,432.0 | 4 | 10.02 |
| base | 29,356 | 2,864.0 / 358.0 | 1,074,000 | 1,575.2 | 64 | 7.5 / 75 | 2,237.5 / 3,222.0 | 10 | 19.69 |
| high | 29,356 | 14,320.0 / 2,864.0 | 3,580,000 | 24,845.2 | 96 | 37.5 / 375 | 13,425.0 / 17,900.0 | 10 | 310.56 |

Operational headroom is explicit rather than a global eligibility ceiling: 25% for storage/inodes, producing 21,480.0 GB and 4,475,000 inodes at the high upper-bound sensitivity; per-job stopping is 1.5× a reviewed high estimate (144 GiB for a 96-GiB mapping estimate). Initial repaired-pilot concurrency remains 2. Full-scale concurrency 10 is only a sensitivity and is not authorized. A repaired ten-pair run must capture all six APE dimensions and PSMC/bootstrap cost before any full-scale resource authorization.

## Branch decisions

| branch | decision | reason |
|---|---|---|
| `core_diversity_psmc` | **CONDITIONAL_GO** | 0/10 passed mandatory assembly QC, callable, consensus, bootstrap, and resource-model gates; identities and immutable core checksums passed. |
| `exact_annotation_partitions` | **CONDITIONAL_GO** | Run only after a core pass and exact annotation accession/version plus sequence-dictionary or validated-liftover binding; annotation absence remains non-vetoing for core. |
| `direct_conversion` | **CONDITIONAL_GO** | H1/H2 differences are not transmitted conversion events; await the separately authorized pedigree/gamete pilot review. |
| `population_gbgc` | **NOT_RUN/DESIGN_ONLY** | No multi-individual population-frequency execution evidence exists in the current graph. |
| `phylogenetic_substitution_bias` | **CONDITIONAL_GO** | The separately authorized H01/H02 pilot produced a pinned metadata preflight but zero verified sequences, callable alignments, substitutions, or biological estimates. |
| `non_allelic_conversion` | **NOT_RUN/DESIGN_ONLY** | No copy-resolved non-allelic execution evidence exists in the current graph. |

## Bounded repair / re-pilot authorization boundary

CONDITIONAL_GO permits only preparation of a versioned repair packet and a newly reviewed ten-slot re-pilot. It does not authorize this review task to acquire new biology or submit jobs, and it does not authorize `scale-vgp-core` to begin full catalog processing.

Before re-pilot submission, every retained primary must have exact-final QV for both haplotypes; BUSCO completeness/missing/duplication for both; manifest-bound k-mer/copy-number audits; exact repeat/low-complexity masks; P09 chemistry resolution; immutable hashes; reviewed resources; and an amendment only if a same-clade alternate is triggered before results. The re-pilot must then demonstrate at least 8/10 callable+consensus passes, all required clades, both generations, no technology-stratum systematic failure, at least 100 PSMC bootstraps and 95% finite success per passing PSMC unit, zero hard violations, and median/p95 resource APE no worse than 25%/50%. No threshold may be relaxed after outcome inspection.

Scale-out may be reconsidered only after an independent review of that repaired immutable packet. Annotation and specialized branches remain orthogonal to core validity.

## Reproducibility

The review generator, its tests, and live sentinel rehash run through the same pinned GNU Guix time-machine/channel and production manifest. The independent output root contains the object ledger, requested-quantity comparison, and validation JSON. The corrected manifests are immutable derivations of the reviewed source digests and explicitly set `biological_jobs_authorized=false` and `full_scaleout_authorized=false`.
