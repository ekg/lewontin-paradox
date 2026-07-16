# Tier 3B population tuple acquisition QC

Acquisition date: 2026-07-16 UTC

Staging root: `/moosefs/erikg/tier3scratch/tier3b-acquisition`
Status: **PASS — two approved, non-synthetic biological population tuples are staged.** The higher-ranked DGRP candidate was fully staged and smoke-tested, then disqualified on an exact native-annotation/reference dictionary mismatch; it does not count toward this passing floor.

## Approved tuples

| tuple | biological population | n | staged all-sites records | records with a non-reference genotype | callable bases | small-region pi |
|---|---|---:|---:|---:|---:|---:|
| `ag1000g_phase3_ao_coluzzii` | *Anopheles coluzzii*, Luanda, Angola, 2009 | 20 | 994,358 | 84,086 | 808,036 | 0.011535141944890723 |
| `ag1000g_phase3_gm_coluzzii` | *Anopheles coluzzii*, Wali Kunda, The Gambia, 2012 | 20 | 994,358 | 125,073 | 808,665 | 0.015079761117724187 |

Both sample sets are wild, locality/year/species-matched subsets of the public Ag1000G Phase 3 metadata. Selection is deterministic: SHA-256 of `dataset_id + NUL + population_id + NUL + sample_id`, lexicographically ranked, take 20. Selected-list checksums are frozen in the manifest.

## Exact-reference and contig checks

- Exact assembly: AgamP4, INSDC `GCA_000005575.1`; staged uncompressed FASTA SHA-256 `19680ed68a6347f59891ecf0ddc9b54f441bd9b71780db95c02c0dddcd809fe7`.
- Each merged VCF has the same eight-contig dictionary as the FASTA: 273,109,044 total reference bases, including `3R` length 53,200,684.
- `bcftools norm --check-ref exit --fasta-ref ... --regions 3R:10000000-10999999` passed every 994,358 record in both VCFs (`0` REF mismatches, `0` realignments, `0` skipped).
- Both BGZF streams passed `bgzip --test`; both tabix indexes answer `bcftools index --stats` with `3R 53200684 994358`.
- VCF uses 1-based inclusive POS. The selected region is VCF `3R:10000000-10999999`. Callable BED files are 0-based half-open. Native GFF3 is 1-based closed. No liftover or contig rename was performed.

## Native annotation identity

The original byte-for-byte AgamP4.12 GFF3 is retained at the manifest path (SHA-256 `9329b95dc4daff2ee084674bb693bc8e3647dbf84831e758a1c2ce8179632e3b`). Its eight `##sequence-region` entries resolve one-to-one to the exact eight-contig FASTA; 14,979 transcripts were parsed. It is the native VectorBase/MalariaGEN annotation for AgamP4, not a congener, lift, projection, or de novo replacement.

The existing strict all-CDS project auditor was also attempted. It reaches CDS reconstruction after dictionary validation, then rejects native transcript `AGAP000192-RA` as empty or ambiguous. This is a documented biological annotation exclusion for the downstream annotation-stratified calculation, not missing provenance and not grounds to mutate or replace the provider GFF. Downstream must exclude that transcript explicitly if it invokes the strict reconstruction path.

## Callable denominator provenance

The provider mask `3R_sitefilters.gamb_colu.dt_20200416.vcf.gz` has 52,226,568 indexed records on 3R. For each exact 20-sample cohort and the staged interval:

1. retain provider sites with `FILTER="PASS"`;
2. retain cohort sites where at least 18 of 20 samples have non-missing GT, `DP>=5`, and `GQ>=20` (`bcftools` expression `N_PASS(FMT/GT!="mis" && FMT/DP>=5 && FMT/GQ>=20)>=18`);
3. intersect the identical one-base `(CHROM, POS0, END)` records; and
4. merge adjacent bases into a 0-based half-open BED.

The source interval contains 812,045 provider-PASS sites. The AO cohort has 940,367 genotype-complete sites and 808,036 intersected callable bases in 38,658 intervals. The GM cohort has 808,665 intersected callable bases in 38,484 intervals. Thus both denominators are exact-sample, nonzero, coordinate-explicit, and reproducible from retained sources.

## Guix-only environment

Every acquisition, parsing, normalization, indexing, checksum, and smoke calculation ran through `analysis/slurm/guix_job.sh` using `analysis/pilot_results/guix_environment.json`. The record pins:

- Guix channel commit `44bbfc24e4bcc48d0e3343cd3d83452721af8c36`;
- channels `analysis/guix/channels.scm`, SHA-256 `45c055cd1d9010a72eacbb720037a22bccb2d8d6891dbd11b5d66365f29b3a17`;
- manifest `analysis/guix/manifest.scm`, SHA-256 `2fb05e87aa2ac45ce51d4dcf93b232cb98627f525adace98357629ee3f15720a`;
- realized profile `/gnu/store/z9v2f6faha9cwjz0sm5iphhlzisgi077-profile`;
- bcftools 1.14, htslib/bgzip/tabix 1.16, samtools 1.14, bedtools 2.30.0, Python 3.10.7, pandas 1.4.4, pysam 0.20.0, and curl 7.85.0 from its recorded Guix store path.

No host Python, pip/conda environment, synthetic input, or Git-staged large biological object was used.

## Small-region smoke calculations

The first 100 kb (`3R:10000000-10099999`) of each staged VCF and its clipped exact-cohort BED were passed to `analysis/tier3b_popvcf_compute.py` inside the pinned environment:

- AO: 87,138 callable sites, diversity numerator 1,005.1491987938879, finite pi 0.011535141944890723.
- GM: 87,198 callable sites, diversity numerator 1,314.9250099433136, finite pi 0.015079761117724187.

Both assertions require positive callable sites, positive numerator, at least one non-reference genotype record, and a finite estimate. Full 1 Mb analysis is intentionally left to the read-only downstream run task.

## Ranked acquisition ledger

The legacy BCM DGRP Freeze 2 VCF directory returned FTP 550 and was not treated as terminal. Acquisition continued within rank 1 to the authoritative DGN 1.1 DGRP consensus archive, whose provider MD5 was verified. A deterministic 20-line Raleigh cohort was staged from 67 provider-QC-eligible lines. Its one-million-record 2L VCF passes exact r5.57 REF checking and contains 25,652 polymorphic records; its consensus-N plus provider IBD/admixture mask yields 831,671 callable bases; and its 100 kb smoke calculation yields 86,827 callable sites, numerator 426.8169934640482, and finite pi 0.0049157173858828266.

DGRP was nevertheless disqualified rather than approved. The original FlyBase FB2014_03 r5.57 GFF declares every one of its 15 `##sequence-region` ends 1–2 bases beyond the corresponding exact r5.57 FASTA/DGN sequence. Rewriting those native declarations would violate the no-substitution/no-silent-repair requirement and leaving them would fail the downstream exact annotation dictionary gate. The complete staged attempt and checksums remain available for audit, but DGRP is absent from `acquisition_manifest.tsv` and does not satisfy the two-tuple minimum.

Acquisition then continued to rank 2 Ag1000G. Ag1000G's old canonical bucket returned HTTP 403/401; acquisition likewise continued to live MalariaGEN public objects. Exact locations and dispositions are in `acquisition_sources.tsv`. Rank-3 resources were not needed after two approved Ag1000G populations were obtained and are not mislabeled as rejected.
