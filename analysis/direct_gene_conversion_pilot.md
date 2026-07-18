# Direct meiotic gene-conversion pilot: fail-closed execution report

**Disposition:** `NOT_EXECUTED_RAW_EVIDENCE_AND_CALLABILITY_GATE`<br>
**Dataset:** D01, *Arabidopsis thaliana*, PRJEB4500 / ERP003793 / eLife 01426 v1<br>
**Execution:** pinned GNU Guix channel `44bbfc24e4bcc48d0e3343cd3d83452721af8c36`, immutable profile `/gnu/store/3c2mxm30rbzvnw7qsi235mrkk3m38fym-profile`, authorized Slurm job `1781172`<br>
**Manifest:** `direct-gene-conversion-pilot-v1.0.0`

## Outcome

D01 passes the relationship and direction gates: the frozen study contains 13 complete male
Col/Ler four-product tetrads (52 products), exact Col/Ler/Cvi parent libraries, and published
reciprocal haplotype, 2:2, 3:1/1:3, crossover, and converted-marker tables. This is direct
transmission structure, not a trio label, static H1/H2 pair, unpolarized genotype pair, population
frequency, or phylogenetic substitution comparison.

The pilot does **not** produce a validated allelic event rate. The upstream versioned acquisition
deliberately superseded all 951 ENA FASTQ objects (356984558868 bytes),
including the bounded 13-tetrad-plus-parent subset of 697 objects
(265809761864 bytes). It acquired no BAMs or per-product genotype/callability masks.
Consequently no candidate can be re-called from all four products, checked for sample exchange or
de novo mutation, inspected against pileups, or passed through candidate-specific mapping,
rearrangement, copy-number, repeat, segmental-duplication, and paralog filters. No empirical error
surface exists for a defensible Mendelian-null simulation or reciprocal spike-in sensitivity.

The Slurm job therefore performs an independent structural reconciliation of the immutable
published tables and emits a fail-closed ledger. It inventories 58 published directional
candidates—44 crossover-associated and 14 non-crossover—but marks every candidate
`PUBLISHED_DIRECTIONAL_CANDIDATE_NOT_RAW_REVALIDATED` and excludes all 58 from direct-rate, tract-distribution, crossover-association,
and GC-bias estimates. An admitted count of zero is not a biological zero.

## Relationship, sex, direction, and ascertainment audit

- The biological units are 13 male meioses recovered as qrt1 pollen tetrads after a Col/Ler donor
  parent was crossed to a Cvi receptor. Every tetrad ID (29, 34, 35, 38, 39, 40, 51, 53, 58, 62, 68, 69, 76) has products 1–4.
  The read table contains 52 tetrad products, three biological parents, and 10 doubled haploids;
  separate Col/Ler DH and tetrad libraries are technical/library provenance, not extra parents.
- All four products establish reciprocal Col/Ler haplotypes and distinguish 2:2 segregation from
  3:1/1:3 conversion. The published NCO table identifies background and converted parental allele;
  the CO marker table identifies the resolved Col or Ler allele in reciprocal products.
- Exact tract-resolution ascertainment is restricted to five high-depth tetrads
  (29, 38, 40, 58, 62). The eight shallow tetrads have approximate crossover rows and no
  equivalent NCO discovery/callability table. Discarding those eight or pretending equal
  sensitivity would bias a per-meiosis rate.
- The ENA snapshot contains 472 runs and 951 provider-MD5-bound FASTQ objects,
  all Illumina HiSeq 2000 genomic WGS. URLs, byte counts, and MD5 cardinalities reconcile, but the
  manifest is provenance, not local raw evidence.
- The material is non-human, so human-participant consent is not applicable. ENA is public and the
  eLife supplements are CC BY 3.0; ENA/EMBL-EBI terms and study citation apply, and public access
  does not transfer third-party rights.

## Assembly and callable opportunity

The exact reference is TAIR10.1 / GCF_000001735.4 with matching RefSeq annotation. The five nuclear
chromosomes contain 119146348 reference bases. Supplement 2B has
137339 filtered Col/Ler markers, of which 137338 are
simple SNVs and 102299 are W/S parental differences. The arithmetic design grid is
1785407 marker–meiosis cells.

None of those quantities is a callable denominator. Callability requires, for each tetrad, the
intersection of confident biallelic parental state; adequate genotype likelihood/depth in all four
products; separable Cvi contribution; continuous phase; unique diploid/copy-normal mapping; and
repeat, rearrangement, and paralog exclusion. Those per-product states are absent. The nuclear span
cannot be used as callable bases, and the common marker list cannot be multiplied by 13 and used as
callable marker opportunities. Events per meiosis, per base-meiosis, and per callable marker are
therefore `NOT_ESTIMABLE_RAW_EVIDENCE_AND_CALLABILITY_GATE` with no confidence interval.

## Candidate tract reconciliation and validation boundary

Supplement 1F contains 71 exactly localized crossovers in the five high-depth tetrads; 44 have a
positive inner converted-tract bound and a corresponding converted-marker inventory. Supplement 1D
contains 14 unique NCO IDs (18 marker rows because four events contain multiple converted markers).
Supplement 2D contains the published tetrad-NCO PCR primer inventory. The parser checks event IDs,
tetrad/product identity, parental direction, inner/outer bounds, marker containment, and duplicated
NCO-row consistency. Inner and outer bounds remain interval-censored; midpoints are not analyzed as
observed lengths.

This is independent computational **table reconciliation only**. It is neither a blinded manual
review nor raw read-backed validation. Candidate-level pileups, all-four-product Mendelian recalls,
copy/paralog masks, and rearrangement/mapping evidence are absent, so the predeclared manual/raw
subset has size zero and all candidates remain excluded. The tract table preserves this status and
reason codes row by row so published candidates cannot inflate the allelic estimate.

## Error model, simulations, and sensitivity

False-discovery and false-negative bounds are not estimable. Genotyping error, depth, allele
balance, phase switches, mapping ambiguity, de novo mutation, crossover uncertainty, structural
rearrangements, copy state, and marker deserts cannot be measured from the selected objects. A
simulation with invented error parameters would not validate this dataset. Mendelian-null
simulation, tract spike-ins, reciprocal W→S/S→W spike-ins, and length/marker-density recovery are
therefore explicit `NOT_RUN_UPSTREAM_GATE` results, as are sensitivities to phasing, mapping, tract
definition, paralog mask, and informative-site ascertainment. The documented five-high-depth versus
eight-shallow split is an observed ascertainment failure, not a corrected sensitivity analysis.

## Direct estimands and GC transmission

No candidate passes the raw-evidence gate, so validated event counts, event rates, interval-censored
tract distributions, and near-CO enrichment are not estimated. Across published candidate rows the
parser can structurally tally simple W/S markers as S=45, W=37, with 1
ambiguous/non-SNV marker rows, but this is not a GC-transmission estimate: linked markers are not
independent event clusters, raw genotype support is absent, reciprocal detection is uncalibrated,
and the preregistered threshold of 100 informative converted mismatches across 50 meioses cannot be
met by 13 meioses. `D_GCBIAS` is `NOT_ESTIMABLE_POWER_AND_INPUT_GATE`, never zero bias.

Population-frequency gBGC, historical phylogenetic substitution bias, and non-allelic conversion
are `NOT_MEASURED`. No H1/H2 evidence is used. No Arabidopsis number is transferred as a human,
vertebrate, or VGP rate, and this blocked pilot supplies no bounded cross-species sensitivity prior.

## Alternate and safe continuation

D02 remains the pre-approved alternate but is not activated. D01 is public and directionally valid;
the present failure is the deliberate omission of raw payload and a per-tetrad callable product, not
an intrinsic access or transmission failure. Activating D02 now would require a versioned pre-result
manifest amendment plus a fresh access/consent and acquisition audit. It must never be replaced by a
VGP H1/H2 pair.

A valid D01 continuation must first authorize and checksum-bind the bounded raw objects (or an exact
study BAM/pileup release), construct per-tetrad callable and exclusion masks on GCF_000001735.4,
re-call all four products and parents, apply assembly-aware copy/paralog/rearrangement filters,
predeclare a blinded candidate subset, and submit a separately fingerprinted pinned-Guix Slurm
analysis with empirical null and reciprocal spike-ins. The blocked-output generator must refuse to
run if such inputs later appear; biological execution then requires a new manifest version.
