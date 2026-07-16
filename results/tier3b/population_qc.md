# Tier 3B recovered population analysis QC

Status: **PASS — every approved non-synthetic biological tuple and frozen ratio completed.**

## Biological estimates

| tuple | population | n | polymorphic SNVs | callable sites | population pi | pi_S/pi_W | secondary sample-jackknife SE (pi) |
|---|---|---:|---:|---:|---:|---:|---:|
| `ag1000g_phase3_ao_coluzzii` | `AO_Luanda_2009_coluzzii` | 20 | 906617 | 14972821 | 0.0123500335652 | 0.742656327614 | 6.17517270676e-05 |
| `ag1000g_phase3_gm_coluzzii` | `GM_WaliKunda_2012_coluzzii` | 20 | 1738162 | 14990338 | 0.0158697942084 | 0.723555012744 | 6.47134640517e-05 |

## Coordinate, reference, genotype, and annotation gates

- VCF positions were converted from 1-based POS to zero-based internal coordinates; callable BED remained 0-based half-open; GFF3 was parsed as 1-based closed and converted once.
- The complete eight-contig VCF and AgamP4 FASTA dictionaries matched. Every streamed REF allele was checked against the exact FASTA; duplicates and coordinate disorder were fatal.
- Multiallelic A/C/G/T SNVs were counted once with all allele counts. Indels, symbolic/non-SNV records, filters, ambiguous reference bases, and sites below 36 called chromosomes were removed from numerator and denominator.
- Callable denominators came only from each exact 20-sample cohort BED. No absent sparse-VCF row was treated as callable reference.
- Native VectorBase AgamP4.12 annotation was consumed byte-for-byte. `AGAP000192-RA` was excluded exactly as declared upstream; all additional canonical CDS failures were enumerated with exact reasons in `population_annotation_exclusions.tsv` and the raw audits rather than editing the GFF.
- 4D sites use nuclear code 1 and forward-reference S=G/C, W=A/T. The frozen 10,000-callable-site gate is enforced without relaxation; any underpowered ratio is explicit in the table and failure ledger.

## Uncertainty

The powered 21-Mb tuples meet the frozen minimum of 20 eligible genomic blocks. Population pi and pi_S/pi_W use the deterministic 10,000-replicate, chromosome-stratified 1-Mb block bootstrap with a percentile 95% interval as their reported uncertainty. The pi_S and pi_W component rows use a 20-replicate delete-one-selected-individual jackknife with a normal-approximation 95% interval and SE, conditional on the exact cohort callable mask. Each TSV row labels its method and resampling unit explicitly; raw JSON retains both uncertainty calculations where both apply.

## Reproducible environment and scheduler

- Guix channel commit: `44bbfc24e4bcc48d0e3343cd3d83452721af8c36`
- Guix profile: `/gnu/store/z9v2f6faha9cwjz0sm5iphhlzisgi077-profile`
- Manifest SHA-256: `2fb05e87aa2ac45ce51d4dcf93b232cb98627f525adace98357629ee3f15720a`
- Channels SHA-256: `45c055cd1d9010a72eacbb720037a22bccb2d8d6891dbd11b5d66365f29b3a17`
- Tool versions: `{"bcftools": "bcftools 1.14", "bedtools": "bedtools v2.30.0", "bgzip": "bgzip (htslib) 1.16", "pytest": "pytest 7.1.3", "python3": "Python 3.10.7", "samtools": "samtools 1.14", "tabix": "tabix (htslib) 1.16", "vcftools": "VCFtools (0.1.16)", "wfmash": "wfmash commit e040aa10e87cab44ed5a4db005e784be62b0bd21 (upstream binary emits no version string)"}`
- Heavy tuple analyses and the independent subset reconciliation ran through Slurm. Exact job IDs and `sacct` resource records are in `population_run_manifest.tsv` and `run_logs/`.
- Array concurrency was throttled to `%1` to avoid two simultaneous full scans of shared MooseFS inputs.

## Independent reconciliation

`population_independent_check.tsv` is produced by a separate standard-library VCF/GT parser on the first 100 kb. The production side consumes an exact tabix-indexed `bcftools view --regions` slice; the manual side reads the original BGZF stream and stops at the same coordinate boundary. It independently clips the BED denominator, checks REF, missingness, filters, multiallelic SNVs, and recomputes pairwise diversity, then requires exact denominator/variant-count agreement and numerical agreement within 1e-12 with the production subset calculation.
