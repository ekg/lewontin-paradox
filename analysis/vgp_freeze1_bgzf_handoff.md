# VGP Freeze 1 BGZF assembly view handoff

## Consumer contract

The immutable upstream mirror remains at
`/moosefs/erikg/vgp/freeze1`.  This task did not modify, rename, or replace
objects below that root.  The derived assembly view is rooted at
`/moosefs/erikg/vgp/derived/freeze1-bgzf`; bulk payloads are deliberately not
stored in Git.

New scale-out work should read `VGP_FREEZE1_ASSEMBLY_BGZF_ROOT` from
`/moosefs/erikg/vgp/derived/freeze1-bgzf/CONFIG.env` and resolve an assembly's
terminal-mirror `source_relative_path` below that value.  Existing pilots keep
their frozen input paths and are not redirected implicitly.

The authoritative machine-readable products are:

- `analysis/vgp_freeze1_bgzf_inventory.tsv`: frozen Slurm input inventory and
  source-to-derived path mapping.
- `analysis/vgp_freeze1_bgzf_manifest.tsv`: one terminal record for every
  mirrored assembly FASTA, including payload/index hashes and provenance.
- `analysis/vgp_freeze1_bgzf_summary.json`: closed-world counts and aggregate
  storage/resource accounting.
- `/moosefs/erikg/vgp/derived/freeze1-bgzf/manifest.tsv` and `summary.json`:
  shared-storage copies for consumers that do not have a repository checkout.

Each accession directory becomes visible through one atomic directory rename
and contains exactly the source-relative `*.fa.gz`, its `.gzi`, its `.fai`,
and `provenance.json`.  The payload and indexes are hardlinks to a
decompressed-SHA-256 content-addressed object only after that object passes
format, digest, and index validation.

## Terminal result

The terminal finalizer accounted for all 581 mirrored assembly FASTAs: 581
were converted, none failed, none were excluded, and none were eligible for
byte-preserving BGZF reuse because every source was ordinary gzip.  Thirty
workers made a reason-coded second attempt after their first allocation found
insufficient space on the nodes' root filesystem; all 30 succeeded on the
explicit node-local `/scratch` mount.  There are 581 complete BGZF/GZI/FAI
triplets and 581 provenance records, with no residual staging object.

Aggregate accounting from `analysis/vgp_freeze1_bgzf_summary.json` is:

- source compressed bytes: 388,785,172,570;
- source decompressed bytes: 1,246,884,198,843;
- derived BGZF bytes: 397,455,489,066;
- GZI plus FAI bytes: 320,131,666;
- BGZF/decompressed compression ratio: 0.31875894283912176;
- worker CPU-hours: 68.16252788583341;
- summed worker elapsed hours: 25.55489021805555;
- maximum observed peak RSS: 46,903,296 bytes;
- maximum observed scratch: 13,154,026,237 bytes;
- terminal failures: 0; and
- retries: 30.

The final shared manifest SHA-256 is
`7aab16fe80156e6f89edc078471ffa8e850a9dc2b160187548da6a78605f0870`.
An independent audit matched the 581 inventory, terminal-mirror, and manifest
source-relative paths; checked every promoted path, size, provenance identity,
CAS hardlink, sequence/name-length record, and probe count; counted exactly
1,743 triplet members; and confirmed that repository and shared-storage copies
of the manifest and summary are identical.

## Validation performed per promoted assembly

Every worker first checked the compressed source SHA-256 against the terminal
mirror manifest.  It streamed ordinary gzip input once through SHA-256, byte
counting, name/length extraction, and pinned htslib `bgzip`, without placing a
full uncompressed assembly on shared storage.  Before promotion it then:

1. ran `bgzip -t` and parsed every BGZF block, including the canonical EOF
   marker;
2. parsed the `.gzi` entry count and required strictly increasing compressed
   and uncompressed offsets within the payload;
3. generated the `.fai` with pinned samtools and compared every contig name and
   length, in order, to the source stream;
4. decompressed the result again and required the exact source-stream SHA-256
   and byte count;
5. ran samtools `faidx` probes at multiple contigs/virtual offsets; and
6. hashed the BGZF, `.gzi`, and `.fai` into accession and CAS provenance.

The triplet was promoted only after all checks passed.  Worker status is
written atomically and records the reason code, attempts, elapsed and CPU time,
peak RSS, peak scratch, sequence count, contig-name/length digest, random probe
count, and BGZF block count.

## Reproduction and operations

The environment is pinned by `analysis/guix/channels.scm` and
`analysis/guix/vgp-freeze1-bgzf-manifest.scm`.  The supported entry point is:

```bash
analysis/run_vgp_freeze1_bgzf_slurm.sh inventory
analysis/run_vgp_freeze1_bgzf_slurm.sh realize
analysis/run_vgp_freeze1_bgzf_slurm.sh submit
analysis/run_vgp_freeze1_bgzf_slurm.sh status
analysis/run_vgp_freeze1_bgzf_slurm.sh retry
analysis/run_vgp_freeze1_bgzf_slurm.sh finalize
```

`submit` and `retry` select only incomplete rows, use bounded arrays by input
size class, and run at low CPU/I/O priority on node-local scratch.  The default
retry ceiling is three attempts.  On this cluster the orchestrator explicitly
sets `VGP_NODE_LOCAL_BASE=/scratch`; sites with a different node-local mount may
override `VGP_BGZF_NODE_LOCAL_BASE` at submission time.  The worker verifies
that the selected filesystem is not a known shared-filesystem type before it
creates scratch.  Invalid or stale destination/CAS objects fail closed instead
of being overwritten.  A failed finalization row remains visible with a reason
code; absence of worker status is reported as `NO_WORKER_STATUS`.

## Repository verification

The task-owned test module passed under the BGZF Guix profile, including a
real-tool streamed conversion/index/promotion integration case.  The complete
repository suite also passed in the pinned comprehensive Guix environment.
Exact commands and terminal counts are recorded in the WG task log.
