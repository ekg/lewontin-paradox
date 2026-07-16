# Tier 3A phased-diploid acquisition QC

Acquisition date: 2026-07-16 UTC. Status: **three eligible, non-synthetic
biological individuals are staged**. Large inputs and derived objects live at
`/moosefs/erikg/tier3data/tier3a-acquisition-20260716`; Git contains only the
acquisition records and parsers. The production handoff is whole-H1/H2
SweepGA mapping with an explicit query:target hit cap, followed by IMPG-native
graph partitioning, annotation-region query, and VCF extraction. BCF encoding
is a downstream bcftools responsibility. Standalone WFMASH products and the
PAF-derived mask/BCF are retained only as mapping/callability sensitivity
evidence, not as the production variant extractor.

## Outcome and ranked public search

The prior read-only audit ended with no eligible individual. Its unnamed pilot
lacked every exact tuple field. The local `mMyoDau2.1` candidate mixed an
annotated H1 with an alternate that is now superseded. NCBI released
`mMyoDau2.hap1.2`/`hap2.2` on 2026-07-15, but did not provide a native
annotation for the exact new H1 version; the historical H2 also lacks a
current Datasets record. This was rejected and the ranked search continued.

The complete NCBI Datasets inventory for VGP parent BioProject PRJNA489243 was
fetched through all three pages: 2,793 reports, 2,761 current reports, and 115
current reports having both an annotation and a linked assembly. Ranking was:
exact reciprocal H1/H2 accessions; same BioSample/isolate; original exact-H1
annotation; explicit haplotype roles; deposited calls and callability; then
phase/collapse and continuity evidence. `GCA_013358685.2` was rejected because
the advertised annotation names `.1`. The first three qualifying tuples were:

1. *Spinachia spinachia* SK-2024b: `GCA_048126635.1` haplotype 1 and
   `GCA_048127205.1` haplotype 2.
2. *Menidia menidia* fMenMen1: `GCA_048628825.1` haplotype 1 and
   `GCA_048544195.1` haplotype 2.
3. *Tautogolabrus adspersus* fTauAds1: `GCA_020745685.1` principal and
   `GCA_020745675.1` alternate pseudohaplotype.

`acquisition_sources.tsv` records every authoritative endpoint attempted,
including rejected candidates and the GIAB dual-modality search. Rejections
are never counted as biological inputs.

## Exact references, annotations, and contigs

For every accepted tuple, reciprocal NCBI assembly links identify the other
exact accession and both reports name the same isolate and BioSample.
Compressed H1/H2 FASTAs and the original H1 GFF were checked against the MD5
in the NCBI FTP checksum ledger. Decompressed FASTAs, FAI files, GFFs, mapping
files, PAFs, BEDs, BCFs, and indexes have SHA-256 values in
`acquisition_manifest.tsv`. Each manifest row also points to the common
`staged_object_inventory.tsv`, whose rows give the relative path,
absolute path, byte size, and SHA-256 of every other file under the acquisition
root; validation compares that set to a live recursive inventory.

The exact H1-native VGL annotations are:

- `GCA_048126635.1-GB_2025_08_04`;
- `GCA_048628825.1-GB_2025_08_13`;
- `GCA_020745685.1-GB_2025_12_08`.

No projected, lifted, congener, or acquisition-generated annotation was used.
All use vertebrate nuclear genetic code 1. Every GFF `##sequence-region`
dictionary is an exact name-and-length bijection to its full H1 FAI: 38/38,
457/457, and 64/64 records, respectively. The contig maps are therefore
identity maps over exact INSDC accessions, not inferred aliases. Each staged
`h1.annotation_provenance.json` binds the annotation release, H1 hash, GFF
hash, genetic code, and contig-map hash.

Each query inventory also contains `h1_h2_contig_map.tsv`, derived from the
selected direct SweepGA PAF rather than from naming assumptions. It records
every selected H1/H2 contig relationship, strand set, mapping-record count,
aligned union on both haplotypes, and every unmapped H1 or H2 contig. Both
FAI dictionaries and every PAF-declared sequence length are checked while the
map is built.

Each tuple has a nonempty `annotation_query/` inventory compiled directly from
that exact native GFF. `gene_manifest.tsv` records every coding target and its
disposition. `impg_query_manifest.tsv` preserves deterministic source-line,
gene, transcript, CDS, protein, locus, strand, coordinate, and phase identity
for every CDS row, including excluded rows and reasons. Queryable transcripts
pass a transcription-order GFF3 phase-chain continuity check; phases 0, 1, and
2 are all present and nonempty in every accepted tuple. Genes are queryable
only when a phase-valid coding transcript exists and the complete H1 gene span
is covered by the selected bounded whole-haplotype mapping.

Only overlapping or touching queryable gene spans are merged in
`impg_execution_spans.bed`; `impg_span_feature_map.tsv` maps every merged span
back to its original genes, transcripts, CDS IDs, and deterministic feature
rows. This BED is an annotation target plan, not a SweepGA input partition.
After `impg partition`, the recorded selector intersects these spans with
IMPG's own `partitions.bed` to create the regional focus BED. Current target
inventory is:

| tuple | targeted genes / union bp | queryable genes / union bp | excluded genes / union bp | queryable CDS rows | merged spans |
|---|---:|---:|---:|---:|---:|
| Spinachia | 20,772 / 275,623,286 | 20,491 / 271,066,985 | 281 / 4,603,748 | 707,609 | 18,200 |
| Menidia | 22,548 / 364,299,923 | 20,953 / 330,626,109 | 1,595 / 34,396,197 | 546,831 | 19,459 |
| Tautogolabrus | 21,614 / 441,517,395 | 20,167 / 374,942,909 | 1,447 / 76,876,202 | 380,862 | 17,769 |

These figures are regenerated after the final direct SweepGA PAF and therefore
may change only if a pinned mapping input changes. `query_qc.json` records
targeted, queryable, and excluded genes/bases, phase counts, merge policy,
hashes, and denominator provenance.

## Phase and collapse evidence

Phase identity passes because the two NCBI reports link the exact other
accession with complementary diploid roles and identify the same biological
individual. Labels were not inferred from local filenames. All three provider
assembly methods explicitly include `purge_dups`; the two newer pairs use
Hifiasm and Hi-C scaffolding, while fTauAds1 uses FALCON-Unzip, linked reads,
Bionano, and Hi-C. Provider total spans, contig/scaffold N50s, and H2:H1 span
ratios are preserved in the manifest. No independent deposited k-mer switch
error was found, which is stated rather than replaced with an invented score.

## Bounded SweepGA mapping and callable denominator

The direct command takes H1 first and H2 second. With the pinned WFMASH backend
the emitted PAF has H1 on the query axis and H2 on the target axis; those axes
are detected from the exact FAI dictionaries and recorded rather than assumed.
IMPG indexes the bidirectional mapping and later restores source coordinates.
SweepGA commit
`018e4ce49d2c125820e0ac50dc5feaa02d423683` reads both complete FASTAs. Direct
FastGA was attempted first and reached a reproducible `GIXmake` signal failure
at whole-assembly index scale; `--batch-bytes 32M` does not affect the simple
two-FASTA branch in this revision. The exact failed command and source-based
rationale are in `acquisition_fastga_attempt.md`.

The pinned parser registers `--num-mappings` but no `-n` short option; `-n`
fails as an unknown argument. Commands therefore use the verified long form
even where the task shorthand says “`-n 1:1`”. This distinction is also proven
from help and source in `analysis/sweepga_impg_handoff.md`.

The ranked fallback is SweepGA's documented WFMASH backend. SweepGA still owns
the complete mapping call and the two-axis cap; WFMASH is an internal aligner
companion rather than a separately promoted PAF:

```text
sweepga H1.fna H2.fna --output-file cap1.paf
  --aligner wfmash --map-pct-identity 90 --min-aln-length 25k
  --num-mappings 1:1 --scaffold-jump 0 --overlap 0.95
  --scoring log-length-ani --threads N
```

`N=8` for Spinachia and Tautogolabrus and `N=32` for Menidia; the per-tuple
value and exact expanded command are recorded in the manifest and do not
change mapping semantics.

No annotation interval enters that command. The `M:N` parser and both plane
sweeps were inspected in source: the left value sets query-axis retained-hit
cardinality and the right value sets target-axis cardinality under the stated
overlap policy. At `--overlap 0.95`, raw interval concurrency can exceed one
when records overlap by less than the rejection threshold, so calling raw
concurrency the `1:1` cap would misstate SweepGA's semantics. Validation
instead reapplies the exact pinned 1:1/0.95 policy to each selected PAF and
requires the ordered 12 mandatory PAF fields to be an unchanged fixed point.
Raw query/target interval depths are recorded separately as descriptive QC.
Sequence keys define coordinate domains for the sweep but are not called graph
or annotation partitions.

Direct whole-haplotype results (coverage is H1-query/H2-target for the detected
axes) are:

| tuple | PAF records | H1/H2 union coverage | raw H1/H2 interval depth | exact-policy recheck |
|---|---:|---:|---:|---|
| Spinachia | 753 | 0.982694 / 0.986971 | 12 / 12 | fixed point |
| Menidia | 3,353 | 0.930498 / 0.951193 | 34 / 53 | fixed point |
| Tautogolabrus | 1,972 | 0.949055 / 0.988533 | 12 / 18 | fixed point |

The raw-depth column is intentionally not labeled a SweepGA hit cap: it counts
every geometric half-open interval overlap, including relationships below the
0.95 rejection-overlap threshold. Each exact PAF SHA-256 and recheck hash is in
the manifest.

The pre-existing standalone WFMASH PAF was separately filtered at 1:1, 5:5,
and 10:10 as cap-sensitivity evidence:

| tuple | 1:1 query/target | 5:5 query/target | 10:10 query/target | selected |
|---|---:|---:|---:|---:|
| Spinachia | 0.943811 / 0.939375 | 0.944509 / 0.939375 | 0.944509 / 0.939375 | 1:1 |
| Menidia | 0.859307 / 0.831034 | 0.860614 / 0.831091 | 0.860663 / 0.831091 | 1:1 |
| Tautogolabrus | 0.949826 / 0.910182 | 0.950226 / 0.910186 | 0.950226 / 0.910186 | 1:1 |

All pairs exceeded the predeclared 80% coverage floor at 1:1, and relaxation
added at most 0.14 percentage points. Direct production therefore stays 1:1;
5:5 and 10:10 are recorded sensitivity caps, not IMPG partitions or queries.

`acquisition_build_mask.py` checks the earlier sensitivity PAF and preserves a
preliminary mapping-callability mask. It checks
PAF/CIGAR consumption and FASTA dictionaries, excludes MAPQ 0/255 and secondary
records, and applies exact uniqueness again as a defensive parser gate. The
callable denominator additionally excludes 100 bp at alignment edges, 10 bp
on both sides of insertions/deletions (including the insertion anchor), and
any non-ACGT H1 or H2 base. Reverse-strand PAF records are complemented and
walked in query orientation before allele comparison. BCF is globally
FASTA-normalized, split, sorted, exact-deduplicated, and CSI-indexed. These
objects validate positive mapping denominators but are not production variant
calls. Final callable bases are annotation target spans intersected with the
selected bounded mapping, successfully queried IMPG-native partitions, exact
H1 REF validation, and explicit ambiguity exclusions. Both preliminary and
final-denominator provenance are named separately in the manifest.

The interval parser was also compared exactly with the original per-base
project parser on the 250 kb Spinachia smoke input: both returned 195,113
callable bases and 146 SNVs. All three WFMASH smoke runs used two threads and
produced nonempty PAF, BED, normalized BCF, and CSI objects. fTauAds1's
scaffold-level H2 required a preliminary homologous-scaffold locator; that
finder PAF and exact selected region are staged.

## IMPG native partition and biological VCF proof

IMPG 0.4.1 at superproject commit
`101df81eb28a809c8fac97d297acd9fcfbbfa048` is the production regional
variant extractor, never a concordance statistic. Its embedded
SweepGA dependency is `ddd31d39b6a68fc972025b048076032341b66835`. The binary
was built while the syng checkout was `dd00f52b688c0fb78cb7f25336ef9ac9f6a3e109`;
the superproject gitlink is separately recorded as
`68ac19745201a7d2a17d9bb190671ef7d3ac8c29`, so the difference is explicit.
The immutable executable and sibling gfaffix hashes are in the manifest.

The dependency task's executed *Spinachia* biological gate is incorporated by
hash as `analysis/sweepga_impg_observed.json`. It used exact H1-native coding
genes, ran direct SweepGA `--num-mappings 1:1`, independently observed maximum
query and target depth 1, built an IMPG index, ran `impg partition -w 2000 -d
0`, and selected nine rows from IMPG's own `partitions.bed` by intersection
with the annotation targets. `impg query -b focus.bed -d 0 -o vcf:poa`
emitted regional VCFs; `impg lace` restored source coordinates. bcftools then
normalized, target-trimmed, exact-deduplicated, emitted VCF.gz/BCF, and created
TBI/CSI indexes.

That native biological region had three queryable genes, 50 CDS rows spanning
all phases 0/1/2, 14,507 callable gene-span bases, nine focused native
partitions, and four non-header biological variants (SNP, insertion, and two
deletions). Every REF allele matched H1; applying the ALTs reconstructed the
aligned H2 segment exactly. The separately controlled fixture is supplemental
and explicitly not the completion gate. IMPG 0.4.1 emits VCF, not BCF; claiming
native BCF output would be false.

## Deposited-call search

The accepted NCBI assembly directories contain no exact-H1 deposited VCF/BCF
plus callability mask. NIST GIAB HG002 v5.0q was inspected as the strongest
public dual-modality alternative, but its calls and benchmark BEDs use
CHM13v2.0, GRCh37, or GRCh38 rather than an accessioned HG002 phased H1 with a
native annotation. It fails the exact-reference gate. No sparse callset was
divided by reference length and no cross-reference callset was mixed with an
assembly; a deposited-call modality was therefore not publicly possible for
the accepted exact-reference tuples.

## GNU Guix environment and versions

Acquisition transforms and scientific commands run through
`analysis/slurm/guix_job.sh`, which fail-closes against the recorded pure
profile. Tool build/runtime libraries are also described by the tracked
`acquisition_toolchain_manifest.scm`. Environment provenance is:

- authenticated channel commit `44bbfc24e4bcc48d0e3343cd3d83452721af8c36`;
- runtime manifest `analysis/guix/manifest.scm` and toolchain manifest
  `results/tier3a/acquisition_toolchain_manifest.scm`;
- runtime profile `/gnu/store/z9v2f6faha9cwjz0sm5iphhlzisgi077-profile`;
- SweepGA 0.1.1, binary SHA-256
  `1a5440529f5eff91cb7d82876a83a5282df66fb5e2c4b1a9c6caa0bdb83de7b1`;
- SweepGA WFMASH companion SHA-256
  `0d8a3a72cfda75a30c38e81b90320c5d212d24b8c312ad22fe97d67e553fc0f6`;
- IMPG 0.4.1, binary SHA-256
  `c587dc2326cd24f887b1fcb3938404229ad0f0a27ef0773e90c287b1ade160d4`;
- gfaffix companion SHA-256
  `4bc1c5e236a8fe6aa1dbcff6e6cf515e8a70c808549990a515c3c5212776a627`;
- bcftools 1.14/HTSlib 1.16, samtools 1.14/HTSlib 1.16, Python 3.10.7,
  GNU Guix 1.4.0-7.44bbfc2, and Rust/Cargo front-end 1.94.1 inside the Guix
  build shell.

`acquisition_commands.sh` records immutable URLs, provider checksum checks,
decompression/indexing, contig audits, the Guix build command, direct whole-
FASTA SweepGA mapping, 1:1/5:5/10:10 sensitivity, query-manifest compilation,
IMPG index/native partition/region-query/lacing, bcftools normalization and
indexing, manifest regeneration, and validation.

## Consumer gate

`run-tier3a-biological-recovery` may select any row with
`eligibility_status=eligible_biological`. It should verify declared hashes and
consume the full H1/H2 FASTAs, exact H1-native GFF, annotation provenance,
direct whole-haplotype SweepGA 1:1 PAF, annotation query manifest, execution
span BED, and deterministic span/feature map. It must then run IMPG index,
native partition, selected-partition regional query, and lace before bcftools
normalization/indexing. The preliminary PAF-derived BED/BCF and standalone
WFMASH objects are sensitivity evidence, not the primary numerator. Every row
carries exact releases and labels, contig mapping, phase/collapse evidence,
hit-cap coverage and depth, targeted/queryable/excluded genes and bases,
positive predeclared target denominator, final callable-denominator rule,
tool-stage policies, and Guix provenance; no unstated input is required.
