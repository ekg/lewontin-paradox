# Exact SweepGA → IMPG handoff

Verification date: 2026-07-16 UTC
Task: `verify-sweepga-impg-handoff`

## Result

The executable handoff is verified, with one necessary correction to the task
wording: SweepGA 0.1.1 at commit
`018e4ce49d2c125820e0ac50dc5feaa02d423683` **does not accept `-n`**.
The literal probe

```text
sweepga -n 1:1 --version
```

exited 2 with `error: unexpected argument '-n' found`.  The supported spelling
is `--num-mappings 1:1`.  The source declares only the long option and describes
it as the query:target mapping-axis limit
([`src/cli.rs:208-212`](https://github.com/pangenome/sweepga/blob/018e4ce49d2c125820e0ac50dc5feaa02d423683/src/cli.rs#L208-L212)).
It would be incorrect and irreproducible to print a successful `-n` command.
The smoke script retains both the failed literal probe and the successful long
form.

With that correction, the complete executed handoff is:

```text
H1/H2 FASTA
  → sweepga --num-mappings 1:1 --scaffold-jump 0 → PAF with cg:Z CIGAR
  → impg index                                           → .impg
  → impg partition (native implicit-graph partitioning)  → partitions.bed
  → impg query -b native-focus.bed -o vcf:poa            → regional VCF
  → impg lace                                            → source-coordinate VCF
  → bcftools norm -Oz/-Ob + bcftools index              → normalized indexed VCF.gz + BCF
```

The staged *Spinachia spinachia* proof is the completion gate.  It has one
SweepGA mapping, query and target overlap
depth 1, nine annotation-intersecting native IMPG partitions, three queryable
H1-native protein-coding genes spanning 14,507 callable bp, and four normalized
biological records (one SNP, one insertion, and two deletions).  Every REF
allele matches H1 and applying all four ALT alleles reconstructs the entire
19,651 bp aligned H2 interval exactly.  The staged whole FASTAs and native GFF
were only read; the proof used newly extracted 20 kb excerpts and never ran
WFMASH.  An additionally executed controlled SNP/indel fixture exactly matches
its three truth records, but is explicitly supplemental and is not part of the
top-level pass predicate.

Machine-readable biological evidence is in
[`sweepga_impg_observed.json`](sweepga_impg_observed.json), expected truth is in
[`sweepga_impg_expected.tsv`](tests/fixtures/sweepga_impg_expected.tsv) for the
optional controlled regression, and the
executable proof is
[`sweepga_impg_smoke.sh`](sweepga_impg_smoke.sh).

## Exact meaning of the SweepGA mapping cap

### Accepted syntax

The parser recognizes `1:1`, unbounded forms, and arbitrary positive `M:N`
pairs.  In particular, `5:5` and `10:10` are accepted; they are not an
assumption.  See
[`src/main.rs:244-277`](https://github.com/pangenome/sweepga/blob/018e4ce49d2c125820e0ac50dc5feaa02d423683/src/main.rs#L244-L277).
The executed dense-overlap fixture retained exactly 1, 5, and 10 records under
`1:1`, `5:5`, and `10:10`, respectively.

The literal successful commands were:

```bash
sweepga raw.paf --output-file cap1.paf  --num-mappings 1:1   --scaffold-jump 0 --overlap 0.95 --scoring log-length-ani --threads 2
sweepga raw.paf --output-file cap5.paf  --num-mappings 5:5   --scaffold-jump 0 --overlap 0.95 --scoring log-length-ani --threads 2
sweepga raw.paf --output-file cap10.paf --num-mappings 10:10 --scaffold-jump 0 --overlap 0.95 --scoring log-length-ani --threads 2
```

### What a query hit and target hit are

A hit is a PAF mapping record whose half-open interval is active at a point in
the relevant coordinate sweep.  It is not a global “at most N records for this
contig” count:

- the query sweep groups records by query sequence inside a genome-pair group;
- the target sweep groups records by target sequence inside that same
  genome-pair group;
- only records passing both sets are retained; and
- a non-overlapping record can survive even when a contig has more than N
  records in total, because the cap constrains concurrent overlap depth.

The grouping and query/target intersection are implemented at
[`src/paf_filter.rs:969-1111`](https://github.com/pangenome/sweepga/blob/018e4ce49d2c125820e0ac50dc5feaa02d423683/src/paf_filter.rs#L969-L1111).
The sweep creates begin/end events and ranks the active set at each coordinate
([query sweep](https://github.com/pangenome/sweepga/blob/018e4ce49d2c125820e0ac50dc5feaa02d423683/src/plane_sweep_exact.rs#L267-L351),
[target sweep](https://github.com/pangenome/sweepga/blob/018e4ce49d2c125820e0ac50dc5feaa02d423683/src/plane_sweep_exact.rs#L354-L432)).
The default score used here is `identity × ln(query length)`
([`src/plane_sweep_exact.rs:24-75`](https://github.com/pangenome/sweepga/blob/018e4ce49d2c125820e0ac50dc5feaa02d423683/src/plane_sweep_exact.rs#L24-L75)).

The proof calculates maximum overlap depth independently from the emitted PAF.
Both controlled and biological `1:1` outputs have query depth 1 and target
depth 1.  The dense cap fixture uses 12 completely overlapping mappings, so its
retained record counts are also its verified query and target depths.

### Ties and repeats

Repetitive hits are ordinary competing mapping records; there is no separate
repeat exemption in this cap.  Among active records, higher score wins.  Equal
scores are ordered by lower start coordinate and then input index, making an
exact score/start tie deterministic
([`src/plane_sweep_exact.rs:161-193`](https://github.com/pangenome/sweepga/blob/018e4ce49d2c125820e0ac50dc5feaa02d423683/src/plane_sweep_exact.rs#L161-L193)).
Records beyond the active cap are tested against the retained records; an
overlap greater than the configured threshold is marked discarded
([`src/plane_sweep_exact.rs:196-258`](https://github.com/pangenome/sweepga/blob/018e4ce49d2c125820e0ac50dc5feaa02d423683/src/plane_sweep_exact.rs#L196-L258)).

The dense fixture deliberately makes its first two records exact ties and tags
them `id:Z:tie_first` and `id:Z:tie_second`.  `1:1` retained only
`tie_first`, directly demonstrating the final input-index tie break.

`--scaffold-jump 0` is essential in this verification.  It ends the SweepGA
pipeline immediately after the mapping-axis sweep
([`src/paf_filter.rs:393-433`](https://github.com/pangenome/sweepga/blob/018e4ce49d2c125820e0ac50dc5feaa02d423683/src/paf_filter.rs#L393-L433)).
Consequently, no SweepGA scaffold construction can be confused with IMPG graph
partitioning.

### Exact handoff format

The handoff is text PAF.  Each line has these 12 required, tab-separated
columns, followed by optional typed tags:

| Column | Field |
|---:|---|
| 1 | query name |
| 2 | query length |
| 3 | query start, 0-based |
| 4 | query end, half-open |
| 5 | strand (`+` or `-`) |
| 6 | target name |
| 7 | target length |
| 8 | target start, 0-based |
| 9 | target end, half-open |
| 10 | matching bases |
| 11 | alignment block length |
| 12 | mapping quality |
| 13+ | optional tags, including `cg:Z:<CIGAR>` |

SweepGA’s parser enforces at least 12 fields and maps them in exactly this order
([`src/paf.rs:135-170`](https://github.com/pangenome/sweepga/blob/018e4ce49d2c125820e0ac50dc5feaa02d423683/src/paf.rs#L135-L170)).
The controlled emitted record is one standard PAF line and contains an
extended `cg:Z:` string with `=`, `X`, `I`, and `D` operations.  IMPG accepts
PAF, 1ALN, or TPA through `-a/--alignment-files`, and `-i` names its reusable
index
([`src/main.rs:4065-4108`](https://github.com/pangenome/impg/blob/101df81eb28a809c8fac97d297acd9fcfbbfa048/src/main.rs#L4065-L4108)).
Its PAF loader consumes the same coordinate fields and records the byte range
of the `cg:Z:` tag
([`src/paf.rs:116-176`](https://github.com/pangenome/impg/blob/101df81eb28a809c8fac97d297acd9fcfbbfa048/src/paf.rs#L116-L176)).

## Exact IMPG interfaces and ownership

| Responsibility | Owner | Verified interface |
|---|---|---|
| H1-versus-H2 alignment and mapping hit limiting | SweepGA | `sweepga H1.fa H2.fa --num-mappings 1:1 --scaffold-jump 0 ...` |
| Mapping ingestion / implicit graph index | IMPG | `impg index -a H1_H2.paf -i H1_H2.impg -t 2` |
| Implicit-graph partitioning | IMPG | `impg partition -a H1_H2.paf -i H1_H2.impg -w ... -d 0 -o bed --output-folder ...` |
| Regional focus | IMPG | `impg query -a H1_H2.paf -i H1_H2.impg -b focus.bed -d 0 ...` where `focus.bed` is selected from IMPG’s own `partitions.bed` |
| Local graph construction and VCF extraction | IMPG | `impg query ... -o vcf:poa --sequence-files H1.fa H2.fa -O calls/` |
| Restore source coordinates / combine regional VCFs | IMPG | `impg lace -l vcf.list --format vcf -o laced.vcf --reference H1.fa --compress none` |
| Annotation intersection and deterministic target trim/deduplication | smoke workflow, after IMPG partition/query | Parse original H1-native GFF; select native partitions intersecting queryable gene spans; normalize, `bcftools view -R target_genes.bed`, then `bcftools norm -d exact` |
| Coding-gene callability | smoke workflow | union length of exact queryable H1-native protein-coding gene spans; explicitly not a genome-wide or CDS-only denominator |
| Normalize and serialize/index VCF and BCF | bcftools in the pinned Guix environment | `bcftools norm -f H1.fa -m -any -Oz -o normalized.vcf.gz laced.vcf && bcftools index -t normalized.vcf.gz`; repeat with `-Ob` and CSI indexing for BCF |

`impg partition` owns the window size and missing-region selection options
([`src/main.rs:4765-4854`](https://github.com/pangenome/impg/blob/101df81eb28a809c8fac97d297acd9fcfbbfa048/src/main.rs#L4765-L4854)).
It writes one BED row per assigned sequence interval with the native partition
number in column 4
([`src/commands/partition.rs:1682-1706`](https://github.com/pangenome/impg/blob/101df81eb28a809c8fac97d297acd9fcfbbfa048/src/commands/partition.rs#L1682-L1706)).
This is why the proof selects focus rows from that output instead of chopping
the input into a preprocessing window scheme.

For a BED-driven graph-like query, `-O` must be a directory and IMPG names one
output per BED column-4 name
([`src/main.rs:10713-10755`](https://github.com/pangenome/impg/blob/101df81eb28a809c8fac97d297acd9fcfbbfa048/src/main.rs#L10713-L10755)).
The `vcf` query branch constructs a local GFA from the query-selected intervals
and converts that graph through POVU
([`src/main.rs:7536-7574`](https://github.com/pangenome/impg/blob/101df81eb28a809c8fac97d297acd9fcfbbfa048/src/main.rs#L7536-L7574),
[`src/main.rs:12148-12189`](https://github.com/pangenome/impg/blob/101df81eb28a809c8fac97d297acd9fcfbbfa048/src/main.rs#L12148-L12189),
[`src/main.rs:10776-10787`](https://github.com/pangenome/impg/blob/101df81eb28a809c8fac97d297acd9fcfbbfa048/src/main.rs#L10776-L10787)).

The `lace` CLI accepts either `-f` inputs or a one-path-per-line `-l` list,
validates `vcf`, accepts a reference and explicit compression mode
([`src/main.rs:4721-4763`](https://github.com/pangenome/impg/blob/101df81eb28a809c8fac97d297acd9fcfbbfa048/src/main.rs#L4721-L4763),
[`src/main.rs:6206-6268`](https://github.com/pangenome/impg/blob/101df81eb28a809c8fac97d297acd9fcfbbfa048/src/main.rs#L6206-L6268)).
For regional VCFs it parses the range suffix on CHROM and adds the range start
to POS while restoring the base contig
([`src/commands/lace.rs:1352-1368`](https://github.com/pangenome/impg/blob/101df81eb28a809c8fac97d297acd9fcfbbfa048/src/commands/lace.rs#L1352-L1368),
[`src/commands/lace.rs:1816-1869`](https://github.com/pangenome/impg/blob/101df81eb28a809c8fac97d297acd9fcfbbfa048/src/commands/lace.rs#L1816-L1869)).

IMPG 0.4.1 emits VCF, not BCF: `query --help` lists `vcf` but no `bcf` output
format.  Therefore the exact non-inferred boundary is “IMPG extracts VCF;
bcftools normalizes, encodes BCF, and creates CSI.”  Describing BCF as a native
IMPG output would be false.  Both normalized encodings are nevertheless
produced and indexed in the handoff proof: BGZF VCF with TBI and BCF with CSI.

## Executed supplemental controlled regression (not a completion gate)

The user correction made native biological evidence the only completion gate.
For extra interface regression coverage, the deterministic H1 is 20,000 bp and
H2 contains:

- SNP `H1#1#chrS:5301 G>T`;
- insertion `H1#1#chrS:5500 G>GTTA`; and
- deletion `H1#1#chrS:5700 GAT>G`.

The literal mapping and IMPG commands executed by the smoke script were:

```bash
sweepga h1.fa h2.fa --output-file h1_h2.1to1.paf \
  --num-mappings 1:1 --scaffold-jump 0 --overlap 0.95 \
  --scoring log-length-ani --threads 2

impg index -a h1_h2.1to1.paf -i h1_h2.impg -t 2

impg partition -a h1_h2.1to1.paf -i h1_h2.impg \
  -w 1000 -d 0 --min-missing-size 1 --min-boundary-distance 0 \
  -o bed --output-folder partitions -t 2

# focus.bed is the H1 row of native partition 5: H1#1#chrS 5000 5999
impg query -a h1_h2.1to1.paf -i h1_h2.impg -b focus.bed \
  -d 0 -o vcf:poa --sequence-files h1.fa h2.fa -O calls -t 2

impg lace -f calls/controlled_truth.vcf --format vcf \
  -o laced.vcf --reference h1.fa --compress none -t 2
bcftools norm -f h1.fa -m -any -Ob -o normalized.bcf laced.vcf
bcftools index -f normalized.bcf
bcftools norm -f h1.fa -m -any -Oz -o normalized.vcf.gz laced.vcf
bcftools index -f -t normalized.vcf.gz
```

SweepGA emitted one PAF record; independently measured maximum query and target
overlap depths were both 1.  IMPG created 20 native partitions and focused a
999 bp partition.  `impg lace` restored the regional VCF to source coordinates;
the normalized, TBI-indexed VCF.gz and CSI-indexed BCF contained exactly:

```text
H1#1#chrS  5301  G    T
H1#1#chrS  5500  G    GTTA
H1#1#chrS  5700  GAT  G
```

Expected and observed `(CHROM, POS, REF, ALT)` arrays are byte-for-byte equal in
the JSON result.

## Executed H1-native annotated biological proof

Input data were read from the staged *Spinachia spinachia* SK-2024b pair:

```text
H1 /moosefs/erikg/tier3data/tier3a-acquisition-20260716/spinachia_spinachia_SK-2024b/h1.fna
   CM106590.1:50001-70000
H2 /moosefs/erikg/tier3data/tier3a-acquisition-20260716/spinachia_spinachia_SK-2024b/h2.fna
   CM106672.1:45001-65000
```

The H1-native annotation is the exact submitted-assembly release
`GCA_048126635.1-GB_2025_08_04`, status `native`, with provenance classification
`native_exact_assembly_submitted_annotation`.  Its audited contig dictionary
maps annotation `CM106590.1` to FASTA `CM106590.1`, length 27,159,969, without a
rename or liftover.  The original GFF SHA-256 recorded by acquisition is
`a984f6e6db964936996d1a0829e3808b4c508c933e21152eea9502dc8ff93c9b`.

The proof parsed that original GFF before deciding which regions were query
targets.  Five protein-coding genes overlap the H1 excerpt.  The deterministic
rule retains genes fully contained in the excerpt and records boundary-truncated
genes as excluded:

| Gene | Original H1 coordinates (1-based, closed) | Strand | Disposition |
|---|---:|:---:|---|
| `gene-AB9W97_009801` | 36,854–52,099 | + | excluded: partial left boundary (2,099 bp in excerpt) |
| `gene-AB9W97_009664` (`tsc22d4`) | 52,959–58,671 | − | targeted and queryable |
| `gene-AB9W97_010084` (`mogat3a`) | 58,914–63,667 | − | targeted and queryable |
| `gene-AB9W97_010351` | 64,450–68,489 | − | targeted and queryable |
| `gene-AB9W97_010352` (`atp1b2a`) | 68,656–76,381 | + | excluded: partial right boundary (1,345 bp in excerpt) |

The three targets total 14,507 non-overlapping gene-span bp.  Their 50 native
CDS rows preserve nine distinct CDS IDs, original transcript parents, strand,
and phase; observed phases are all of `0`, `1`, and `2`.  This denominator is
therefore a **targeted protein-coding gene-span** denominator, not a genome-wide
estimate and not a CDS-only estimate.

The proof extracted the two 20 kb FASTA regions with `samtools faidx`, renamed
only the temporary excerpt headers, and ran the same literal SweepGA/index
sequence.  After SweepGA, all three exact local gene intervals are fully covered
by the retained one-to-one query mapping, so the cap excludes zero target genes
and zero target bp.  The native partition command used `-w 2000`; the workflow
selected the nine rows of IMPG's own `partitions.bed` that intersect at least
one queryable target.  It did not create preprocessing windows or add an
arbitrary padding distance:

```bash
sweepga bio_h1.fa bio_h2.fa --output-file bio.1to1.paf \
  --num-mappings 1:1 --scaffold-jump 0 --overlap 0.95 \
  --scoring log-length-ani --threads 2
impg index -a bio.1to1.paf -i bio.impg -t 2
impg partition -a bio.1to1.paf -i bio.impg \
  -w 2000 -d 0 --min-missing-size 1 --min-boundary-distance 0 \
  -o bed --output-folder partitions -t 2
impg query -a bio.1to1.paf -i bio.impg -b focus.bed \
  -d 0 -o vcf:poa --sequence-files bio_h1.fa bio_h2.fa -O calls -t 2
impg lace -l vcf.list --format vcf -o laced.vcf \
  --reference bio_h1.fa --compress none -t 2
bcftools norm -f bio_h1.fa -m -any -Ob -o normalized.untrimmed.bcf laced.vcf
bcftools index -f normalized.untrimmed.bcf
bcftools view -R target_genes.bed -Ou normalized.untrimmed.bcf \
  | bcftools norm -d exact -Ob -o normalized.bcf
bcftools index -f normalized.bcf
bcftools view -Oz -o normalized.vcf.gz normalized.bcf
bcftools index -f -t normalized.vcf.gz
```

IMPG's selected native partitions cover 18,002 bp.  The 3,495 bp beyond the
exact 14,507 target bp are native-partition query context, not user-specified
padding.  Normalization occurs before target trimming so left-aligned indels
are assigned deterministically.  `view -R` trims to exact gene spans and
`norm -d exact` deduplicates records that could have arisen from adjacent
regional VCFs.  The observed run had four records before and four after this
trim/deduplication, proving that none was spuriously dropped or duplicated.

Observed QC:

| Metric | Value |
|---|---:|
| SweepGA PAF mappings | 1 |
| Maximum query overlap depth | 1 |
| Maximum target overlap depth | 1 |
| H1-native coding genes overlapping excerpt | 5 (17,951 bp within excerpt) |
| Targeted/queryable coding genes | 3 / 3 |
| Boundary-excluded coding genes | 2 (3,444 bp within excerpt) |
| One-to-one-mapping-excluded target genes | 0 (0 bp) |
| Native IMPG partitions focused | 9 |
| Native partition focus context | 18,002 bp |
| Exact targeted/callable denominator | 14,507 gene-span bp |
| Normalized biological records | 4 |
| VCF TBI present | yes |
| BCF CSI present | yes |

The four local-coordinate records (source genomic coordinate is local POS +
50,000) are:

| Local POS | Source H1 POS | REF | ALT | Kind | Exact annotation target |
|---:|---:|---|---|---|---|
| 3,429 | 53,429 | T | TAC | insertion | `gene-AB9W97_009664` |
| 7,967 | 57,967 | TGC | T | deletion | `gene-AB9W97_009664` |
| 14,700 | 64,700 | G | A | SNP | `gene-AB9W97_010351` |
| 15,504 | 65,504 | TTA | T | deletion | `gene-AB9W97_010351` |

Direct sequence validation does not merely trust VCF syntax.  For every record,
the script extracts the indexed H1 bases at POS and asserts equality with REF.
It then applies the four ALT alleles, in coordinate order, to the SweepGA-aligned
H1 interval `[0,19653)` and compares the result base-for-base with the aligned
H2 interval `[349,20000)`.  Both sequences are 19,651 bp after applying the
alleles and are identical.  Thus representative biological SNP and indel calls
are directly confirmed against both haplotype sequences.

The proof does not treat sparse VCF records, raw partition context, or the full
20 kb excerpt as the denominator.  The union of the three exact queryable
H1-native protein-coding gene spans supplies the explicit positive denominator.

## GNU Guix pin and executable builds

The authenticated channel is
[`analysis/guix/channels.scm`](guix/channels.scm), Guix commit
`44bbfc24e4bcc48d0e3343cd3d83452721af8c36`.  The smoke-only manifest is
[`analysis/guix/sweepga_impg_smoke_manifest.scm`](guix/sweepga_impg_smoke_manifest.scm).
The command executed to resolve/build it was:

```bash
guix time-machine -C analysis/guix/channels.scm -- \
  build -m analysis/guix/sweepga_impg_smoke_manifest.scm
```

It resolved Bash 5.1.16, coreutils 9.1, findutils 4.9.0, gawk 5.2.1,
grep 3.8, sed 4.8, Python 3.10.7, samtools 1.14, bcftools 1.14,
git-minimal 2.40.1, GCC 12.3.0, and all native libraries to `/gnu/store`;
`guix shell` materialized profile
`/gnu/store/8x4hx7d9hnv187yprjrzqyg0kxj2z32k-profile`.

The frozen 2023 Guix channel provides Rust 1.67, which cannot compile these
locked 2026 dependency graphs.  The manifest therefore imports
[`analysis/guix/packages/rust-binary.scm`](guix/packages/rust-binary.scm): a
Guix package around the official fixed Rust 1.94.1 release (source SHA-256
`294b3d81fa72e62581276290c60c81eb8b58498d333d422ca1dfc432877d0c40`,
Guix base32 `0h0cgn3k5i6zl4n44g9kim4mi2zbh46cd4324y0jbrkjza0ksjr9`).  That package
rewrites its ELF interpreter and runpaths to the frozen channel's glibc,
GCC-runtime, and zlib store items.  Since the channel's patchelf 0.11 cannot
parse this newer ELF layout, the same Scheme file builds patchelf 0.18.0 from
fixed source (SHA-256
`1952b2a782ba576279c211ee942e341748fdb44997f704dd53def46cd055470b`,
Guix base32 `02s7ap86rx6yagfh9xwp96sgsj0p6hp99vhiq9wn4mxshakv4lhr`).  No rustup,
user-profile compiler, system compiler, or system library participates.
`rustc --version` reported `rustc 1.94.1 (e408947bf 2026-03-25)` and Cargo
reported `cargo 1.94.1 (29ea6fb6a 2026-03-24)`.

The entire proof runs through:

```bash
guix time-machine -C analysis/guix/channels.scm -- \
  shell -m analysis/guix/sweepga_impg_smoke_manifest.scm --pure -- \
  bash analysis/sweepga_impg_smoke.sh --inside-guix [optional-output-directory]
```

The user-facing entry point wraps that command automatically:

```bash
analysis/sweepga_impg_smoke.sh
```

Executable builds actually exercised (not merely version strings) were:

| Program | Source commit | Version | Executable SHA-256 |
|---|---|---|---|
| SweepGA | `018e4ce49d2c125820e0ac50dc5feaa02d423683` | 0.1.1 | `1a5440529f5eff91cb7d82876a83a5282df66fb5e2c4b1a9c6caa0bdb83de7b1` |
| IMPG | `101df81eb28a809c8fac97d297acd9fcfbbfa048` | 0.4.1 | `c587dc2326cd24f887b1fcb3938404229ad0f0a27ef0773e90c287b1ade160d4` |
| gfaffix companion | submodule `460e0dd798a9da7d12aef4f9181419d71489da95` | IMPG companion | `4bc1c5e236a8fe6aa1dbcff6e6cf515e8a70c808549990a515c3c5212776a627` |

These are the literal build settings used for SweepGA (the IMPG invocations
used `/tmp/impg-guix-native-build` and `--bin impg`, then `--bin gfaffix`):

```bash
rm -rf /tmp/sweepga-guix-native-build /tmp/guix-link-shims
mkdir -p /tmp/guix-link-shims
guix time-machine -C analysis/guix/channels.scm -- \
  shell -m analysis/guix/sweepga_impg_smoke_manifest.scm --pure -- bash -c '
    export HOME=/tmp/guix-rust-home
    export CARGO_HOME=/home/erikg/.cargo
    export CARGO_TARGET_DIR=/tmp/sweepga-guix-native-build
    export CARGO_NET_OFFLINE=true CC=gcc CXX=g++
    export CARGO_TARGET_X86_64_UNKNOWN_LINUX_GNU_LINKER=gcc
    export LIBCLANG_PATH="$GUIX_ENVIRONMENT/lib"
    for pair in rt:1 pthread:0 dl:2 util:1; do
      base=${pair%:*}; sonum=${pair#*:}
      ln -sf "$GUIX_ENVIRONMENT/lib/lib${base}.so.${sonum}" \
        "/tmp/guix-link-shims/lib${base}.so"
    done
    export LIBRARY_PATH="/tmp/guix-link-shims${LIBRARY_PATH:+:$LIBRARY_PATH}"
    export RUSTFLAGS="-L native=/tmp/guix-link-shims"
    cargo build --release --locked --offline \
      --manifest-path /moosefs/erikg/sweepga/Cargo.toml --bin sweepga
  '
```

`CARGO_HOME` supplies only Cargo's pre-fetched, content-addressed source cache;
`--locked --offline` prevents dependency resolution or network access.  Every
program executed by the build, including Cargo itself, Rust, GCC, CMake,
pkg-config, patchelf, and link-time libraries, comes from the pinned Guix
profile.  SweepGA's FastGA companion executables and IMPG's gfaffix companion
were staged beside the primary binaries.  The `wfmash` companion is present
only because SweepGA builds it; neither proof invokes it.

The smoke refuses any source/submodule commit or executable/companion hash
drift before running.  It also checks the ELF interpreter and every `ldd`
dependency of SweepGA, IMPG, and gfaffix: all resolve under `/gnu/store` (glibc
2.35 and GCC 12.3 runtime outputs).  This proof therefore tests the newly
recorded Guix builds, not an older development binary.

## Reproduction and interpretation

Run:

```bash
analysis/sweepga_impg_smoke.sh
```

With no argument, the script uses a temporary working directory and publishes
only `analysis/sweepga_impg_observed.json`.  The default executes the biological
native-annotation gate and mapping-cap probes.  To reproduce the separately
reported supplemental controlled regression too, run:

```bash
RUN_SUPPLEMENTAL_CONTROLLED=1 analysis/sweepga_impg_smoke.sh
```

Pass an output directory to retain all PAF, IMPG index, partition BED, VCF,
BCF/CSI, help captures, stderr logs, and multiplicity fixtures.  The published
JSON records `controlled.completion_gate=false`, whether that optional fixture
was executed, and the biological evidence used by the top-level `pass` value.

This verification authorizes only the demonstrated division of labor.  It does
not authorize the former whole-assembly WFMASH-first pipeline, does not move
graph partitioning into SweepGA, and does not replace IMPG native partitioning
with preprocessing windows.
