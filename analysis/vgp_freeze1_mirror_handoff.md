# VGP Phase 1 Freeze 1 mirror handoff

**Release definition:** VGP/vgp-phase1 commit `dc1b2af5a7741b97d66fb10cb2bce97f41765cdf`, file
`VGPPhase1-freeze-1.0.tsv`, SHA-256 `9c58420484a8b76a2d6175b7c26bf709e68bdc726a67fc7541b8c2b5a2fc13a4`.  The local pinned copy
is exactly 327,466 bytes, 717 physical lines,
and 716 data rows.

## Closed-world reconciliation

All 716 catalog rows have exactly one disposition:
581 rows map one-to-one to unique exact UCSC
browser accession/version roots; 135 catalog
rows have an empty UCSC browser accession and therefore contribute no released
hub object.  The snapshot contains all 581
expected roots, no missing root, no extra root, and no root-level accession or
version drift.  The moving current hub was used only as transport.

## Authoritative metadata-only inventory

The recursive rsync dry-run completed at `2026-07-18T11:19:25Z`
from `rsync://hgdownload.soe.ucsc.edu/hubs/VGP/` (product paths resolve against `rsync://hgdownload.soe.ucsc.edu/hubs/`).
The bound raw listing SHA-256 is `dce0ce8872931c1c8717e2196bc45f9718375413fdf4d26f9187e9639ca6a878`.
It inventoried 43,371 files and
3,916,877,494,936 file bytes.  Exact assembly FASTA:
581 files / 388,785,172,570
bytes.  Exact assembly 2bit: 581 files /
341,066,459,447 bytes.  Their union is
1,162 files /
729,851,632,017 bytes.

The historical approximately 967 GB whole-collection and approximately 520 GB
FASTA-only figures are **unverified historical planning estimates only**.  They
are not checksums, completion criteria, byte ceilings, or evidence.  They have
been replaced by the frozen inventory above.

All 581 official
`md5sum.txt` objects are inventoried.  The worker promotes and validates those
first, binds every checksum that names an exact frozen path, and refuses all
remaining payload if any checksum manifest fails or names an extra/renamed
object.  Objects without a published checksum receive the three-pass local
SHA-256 contract described below.

## Capacity and launch gate

Required durable bytes: 3,916,877,494,936; source-relative staging at
bounded concurrency 2: 119,443,354,696;
checksum/manifests: 67,108,864; one-object quarantine:
59,721,677,348; explicit 20%
operational headroom: 783,375,498,987.  Combined
worst-case requirement: 4,879,485,134,831 bytes and 57,455
inodes.  Filesystem evidence reports 101,082,651,820,032 bytes and
1,004,580,301 inodes available.  No arbitrary global byte or
memory cap is imposed.

Capacity gate: **quota_visibility_unavailable_fail_closed**.  Filesystem capacity passes, but
bulk execution remains fail-closed unless the mounted filesystem's applicable
quota can be queried successfully.  This run did not infer quota absence from
`df`.  The exact failed/unavailable quota probes are preserved in the JSON
summary.  Consequently no bulk payload transfer or Slurm job was launched.

## Harmless transfer fixture

Pinned GNU Guix rsync was deliberately interrupted after
262,144 of
8,388,608 bytes.  Resume added exactly
8,126,464 bytes and produced
SHA-256 `7fabb846ab516f6f68ccd0a750c8cf84649180241a40431fa683562bb2464f54`.  The
fixture also proved checksum-failure quarantine and atomic promotion:
`True`.  Its durable JSON report is
`/moosefs/erikg/lewontin-paradox-data/vgp/freeze1/fixture/fixture-report.json`.

## Durable operation

Run only through the pinned wrapper:

```bash
analysis/run_vgp_freeze1_mirror.sh inventory
analysis/run_vgp_freeze1_mirror.sh build
analysis/run_vgp_freeze1_mirror.sh fixture
analysis/run_vgp_freeze1_mirror.sh worker --concurrency 2
```

`inventory` always creates a new timestamped snapshot and never replaces the
bound snapshot.  To reconcile a deliberate refresh, pass the new raw and
metadata paths explicitly to `build` and review the closed-world result before
allowing the worker to see its capacity evidence.

The worker uses source-relative `.part` staging, `rsync --partial
--append-verify`, bounded concurrency, exponential backoff, size plus published
MD5 verification when available, local SHA-256 before/re-before/after atomic
promotion, and source timestamps as metadata only.  A mismatch is moved to a
source-relative quarantine location and never overwrites or deletes the last
verified destination.  There is no mirror-wide delete operation.  Each object
becomes durable immediately after its own verification.

## Current mutually exclusive accounting

{
  "missing": {
    "bytes": 0,
    "files": 0,
    "objects": 0
  },
  "planned": {
    "bytes": 3916877494936,
    "files": 43371,
    "objects": 47870
  },
  "quarantined": {
    "bytes": 0,
    "files": 0,
    "objects": 0
  },
  "reused": {
    "bytes": 0,
    "files": 0,
    "objects": 0
  },
  "superseded": {
    "bytes": 0,
    "files": 0,
    "objects": 0
  },
  "transferred": {
    "bytes": 0,
    "files": 0,
    "objects": 0
  },
  "verified": {
    "bytes": 0,
    "files": 0,
    "objects": 0
  }
}

No raw sequencing-read archive is in scope.  No bulk payload belongs in git.
