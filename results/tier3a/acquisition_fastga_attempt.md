# Whole-haplotype SweepGA backend decision

Date: 2026-07-16 UTC. SweepGA commit:
`018e4ce49d2c125820e0ac50dc5feaa02d423683`.

The preferred direct FastGA backend was attempted first on the complete
*Spinachia spinachia* H1/H2 FASTAs, without annotation partitioning:

```bash
sweepga h1.fna h2.fna --output-file cap1.paf \
  --aligner fastga --num-mappings 1:1 --scaffold-jump 0 \
  --overlap 0.95 --scoring log-length-ani --threads 8
```

`FAtoGDB` first exposed a one-argument conversion crash on MooseFS. Prebuilding
the GDB with its explicit output argument allowed SweepGA to continue, but
`GIXmake -T8 ... h1` then terminated by signal (`status.code() = None`) before
any mapping record was produced. Adding `--batch-bytes 32M` did not alter this
two-FASTA CLI branch. SweepGA's own source describes practical FastGA index-size
failures and routes large-batch errors through `IndexCreationError`; the input
H1 is 407,541,755 bp, far beyond the cited small-index range.

The ranked fallback therefore keeps SweepGA as the single whole-FASTA mapping
and cardinality owner, using its documented WFMASH backend:

```bash
sweepga h1.fna h2.fna --output-file cap1.paf \
  --aligner wfmash --map-pct-identity 90 --min-aln-length 25k \
  --num-mappings 1:1 --scaffold-jump 0 --overlap 0.95 \
  --scoring log-length-ani --threads 8
```

The backend resolves to the Guix-built companion whose SHA-256 is
`0d8a3a72cfda75a30c38e81b90320c5d212d24b8c312ad22fe97d67e553fc0f6`.
No annotation region is supplied to SweepGA. Annotation targets are compiled
only after the bounded whole-haplotype PAF exists; IMPG subsequently owns graph
indexing, native partitioning, regional query, and VCF extraction.
