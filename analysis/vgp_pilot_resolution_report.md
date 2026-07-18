# Repaired VGP pilot candidate resolution

Date: 2026-07-18 UTC

Task: `repair-vgp-candidate`

Resolver: `analysis/resolve_vgp_candidates.py` (`vgp-candidate-resolver/2.0`)

## Outcome

Six clade-stratified composition candidates were resolved from immutable,
version-pinned official NCBI responses. All six have an exact versioned GCF H1,
an exact NCBI annotation release whose reference is that same H1, an official
complete sequence-report linkage scoped to that H1, deterministic official
FASTA/GFF locations, and finite exact compressed sizes.

The repaired execution selection is nevertheless **zero**. The explicit
numeric pilot vector passes, but the integrated data-root evidence says that no
user-visible quota interface is available. Free filesystem space is not a
substitute for a quota/allocation limit, and the stricter unknown cap therefore
wins. The irreducible blocker is `QUOTA_UNAVAILABLE`; the resolution decision
is `NO_GO_QUOTA_UNAVAILABLE`.

This is not a backfill. The six repaired rows are composition-only. No GCA H1/H2
pair was transferred to a neighboring GCF accession, and no Tier3A diversity row
was emitted. A later Tier3A row would still require an exact H2 plus affirmative
same-individual and correctly phased H1/H2 evidence.

## Frozen-catalog units

The integrated catalog is
`VGPPhase1-freeze-1.0.commit-dc1b2af5a7741b97d66fb10cb2bce97f41765cdf.tsv`,
SHA-256
`9c58420484a8b76a2d6175b7c26bf709e68bdc726a67fc7541b8c2b5a2fc13a4`.
Its units are deliberately reported with their own denominators:

| Quantity | Count | Meaning |
| --- | ---: | --- |
| Physical catalog lines | 717 | Header plus data rows |
| Header lines | 1 | Schema, not a species record |
| Data rows | 716 | Catalog records after the header |
| Unique species names | 714 | Deduplicated nonempty `Scientific Name` values |
| Data-row excess over unique species | 2 | Duplicate occurrences, not missing content |

The duplicated species are `Lophostoma evotis` (multiplicity 2) and
`Micronycteris microtis` (multiplicity 2). Thus 717, 716, and 714 are all
simultaneously correct and are not compared as if they shared a denominator.

## Prioritized resolution set

| Priority | Clade | Catalog species | Repaired exact H1 | Native annotation release | Disposition |
| ---: | --- | --- | --- | --- | --- |
| 1 | Mammals | *Camelus dromedarius* | `GCF_036321535.1` | `GCF_036321535.1-RS_2024_04` | Metadata eligible; quota blocked |
| 2 | Birds | *Colius striatus* | `GCF_028858725.1` | `GCF_028858725.1-RS_2023_12` | Metadata eligible; quota blocked |
| 3 | Reptiles | *Candoia aspera* | `GCF_035149785.1` | `GCF_035149785.1-RS_2024_02` | Metadata eligible; quota blocked |
| 4 | Amphibians | *Dendropsophus ebraccatus* | `GCF_027789765.1` | `GCF_027789765.1-RS_2024_11` | Metadata eligible; quota blocked |
| 5 | Actinopterygii | *Lepisosteus oculatus* | `GCF_040954835.1` | `GCF_040954835.1-RS_2024_09` | Metadata eligible; quota blocked |
| 6 | Chondrichthyes | *Heterodontus francisci* | `GCF_036365525.1` | `GCF_036365525.1-RS_2024_08` | Metadata eligible; quota blocked |

The repaired manifest contains exactly these six candidates. The repaired
rejection ledger records a precise disposition for all 74 prior candidate
rows: the six above are blocked by the stricter quota cap and the other 68 are
outside the authorized six-species prioritized metadata set.

## Recomputed cap vector

The high resource aggregate for the six metadata-eligible rows is compared to
the authorization vector below. Scratch is the sum of the two largest
per-species high estimates because concurrency is limited to two.

| Dimension | Proposed high value | Maximum | Numeric status |
| --- | ---: | ---: | --- |
| Species | 6 | 6 | Pass |
| Compressed inputs | 3.865305 GiB | 120 GiB | Pass |
| Concurrent scratch | 205.5066 GiB | 750 GiB | Pass |
| Aggregate core-hours | 11.1353 | 1,500 | Pass |
| Concurrent species | 2 | 2 | Pass |
| Peak memory per job | 32 GiB | 256 GiB | Pass |
| Integrated per-user quota | Unknown/unavailable | Must be known and sufficient | **Block** |

Every stricter integrated cap wins, so no row is selected and neither
acquisition nor execution is authorized.

## Cache, retry, and provenance contract

The cache index is `analysis/vgp_resolution_cache/index.json`. It binds 26
official response objects: two batched NCBI Datasets responses, six complete
exact-accession sequence reports (the API exposes no working multi-accession
sequence-report batch), six official checksum catalogs, and twelve HEAD-only
size probes. Each entry contains the normalized request, source/version,
endpoint, retrieval time, response SHA-256, response size, response headers,
and software/environment identity.

Cache keys are SHA-256 over canonical normalized request plus source/version.
Existing objects are immutable. Refreshes honor `Retry-After`, use a respectful
minimum interval, bounded exponential backoff with jitter, at most six
attempts, and the resumable checkpoint
`analysis/vgp_resolution_cache/checkpoint.json`. Offline regeneration reads
only those cached responses.

## Acquisition and run obligations

A remote SHA-256 or MD5 is not a pre-download eligibility condition. For any
future independently gated acquisition, every exact payload must instead be
fully staged, checked against its official expected byte size, locally hashed
with SHA-256, checked against the official MD5 when one exists, have the local
SHA-256 reverified immediately before atomic read-only promotion, and be
quarantined on any mismatch. The per-payload JSON obligations and official MD5
values are recorded in the repaired manifest and size budget.

Callable bases, callable fraction, queryable protein-coding gene count, and
queryable protein-coding CDS bases are not pre-download requirements. They are
post-alignment acceptance measurements under
`vgp_post_alignment_denominators_v1`. The executable contract is:

```sh
python3 analysis/resolve_vgp_candidates.py measure \
  --metrics-json METRICS.json \
  --output-json ACCEPTANCE.json
```

Results are excluded if any measurement is absent or if callable bases are
below 10,000,000, callable fraction is below 0.50, queryable gene count is
below 1,000, or queryable gene bases are below 1,000,000.

## Old versus repaired disposition and authorization boundary

| Artifact state | Manifest candidate rows | Selected | Rejection-ledger rows |
| --- | ---: | ---: | ---: |
| Integrated prior refusal | 74 | 0 | 74 |
| Repaired resolver | 6 | 0 | 74 |

The integrated prior TSVs are preserved byte-for-byte under
`analysis/vgp_resolution_cache/prior_refusal/`; the lineage, review, synthesis,
run, and refusal evidence were not rewritten. This repair downloaded zero
biological payload bytes, submitted zero jobs, and performed no demographic
inference. Full-catalog acquisition, raw population bulk download, and all
demographic inference remain unauthorized.
