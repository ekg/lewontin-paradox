# VGP pilot repair definition quality audit

Date: 2026-07-18 UTC

WG task: `quality-vgp-pilot-2`

## Decision

All seven downstream definitions were reviewed and edited before dispatch.
The resulting graph retains a single fail-closed acquisition/compute path,
preserves the small-pilot authorization, and repairs the observed catalog,
checksum, and denominator semantics without relaxing exact reference,
annotation, haplotype-pair, budget, or storage gates.

This quality audit launched no download, transfer, Slurm submission, cluster
job, or demographic inference. Read-only inspection was limited to the
already frozen catalog, WG state, committed analysis artifacts, and tests.

One material issue was removed from the draft definition. The original
`regate-vgp-pilot` text allowed a GO decision to rely on filesystem free space
when per-user quota was unknown. The integrated audit records
`QUOTA_UNAVAILABLE` as an independent blocker
([integration ledger, lines 292-307](vgp_pilot_lineage_integration.md#known-repair-blockers)).
Free space is not a substitute for a stricter quota/allocation contract, so
the repaired definition now keeps unknown quota unknown and requires every
stricter integrated storage constraint to win.

## Frozen-catalog evidence and corrected units

The exact frozen object was inspected in place at:

`/moosefs/erikg/lewontin-paradox-data/vgp/phase1-freeze-1.0/manifests/VGPPhase1-freeze-1.0.commit-dc1b2af5a7741b97d66fb10cb2bce97f41765cdf.tsv`

Read-only checks reproduced the committed provenance SHA-256
`9c58420484a8b76a2d6175b7c26bf709e68bdc726a67fc7541b8c2b5a2fc13a4`
and the integrated ledger's 327,466-byte/717-line record
([integration ledger, lines 140-146](vgp_pilot_lineage_integration.md)).

The catalog units are:

| Unit | Observed value | Meaning |
| --- | ---: | --- |
| Physical lines | 717 | One header plus all tabular data rows. |
| Header lines | 1 | The schema/header, not a species record. |
| Data rows | 716 | Tabular records after the header. |
| Unique scientific names | 714 | Deduplicated values of `Scientific Name` (column 10). |
| Excess rows over unique names | 2 | Duplicate occurrences, not missing catalog content. |

The two duplicated names each occur twice: `Lophostoma evotis` and
`Micronycteris microtis`. Therefore 717 lines, 716 data rows, and 714 unique
species are simultaneously true. None should be compared directly as if it
were the same denominator. The downstream resolver and gate definitions now
require all units plus duplicate identities and multiplicities.

The prior integrated refusal remains immutable evidence: zero selected rows,
gate `NO_GO`, zero transferred bytes, `NOT_SUBMITTED`, no inferred Ne, review
`FAIL`, and synthesis `stop_repair` are recorded at
[integration ledger lines 270-286](vgp_pilot_lineage_integration.md#refusal-evidence-that-must-not-be-rewritten).
The repair tasks must produce new, provenance-linked artifacts; they must not
rewrite those facts to imply that the earlier biological pilot ran.

## Definition review ledger

### `repair-vgp-candidate`

Edited and logged. The definition now:

- states the 717/1/716/714 units and duplicate multiplicities explicitly;
- requires NCBI/VGP batching, immutable cache keys, rate limits,
  `Retry-After`, bounded retry/backoff, resumable checkpoints, response
  digests, timestamps, endpoints, and environment provenance;
- permits an official exact-version payload without remote SHA-256 or MD5
  only by carrying a staged size-check/local-SHA-256 obligation forward;
- preserves exact accession/version, deterministic official retrieval,
  native H1 annotation/reference linkage, and affirmative same-individual,
  phased H1/H2 evidence for every Tier3A diversity row;
- moves callable/queryable gene/base denominators to a thresholded
  post-alignment measurement contract; and
- states all small-pilot ceilings and forbidden full-catalog, population-bulk,
  and demographic operations.

### `regate-vgp-pilot`

Edited and logged. It independently reproduces catalog units, duplicate
identities, eligibility, exact linkages/pairs, finite sizes, remaining staged
checksum obligations, and the denominator measurement contract. The GO/NO_GO
digest binds the manifest, catalog provenance, data-root/storage contract,
environment, cap vector, retrieval/checksum obligations, pair evidence, and
measurement contract. Unknown quota cannot be promoted to known capacity, and
free-space evidence cannot override a stricter integrated storage gate.

### `acquire-repaired-vgp`

Edited and logged. It refuses `NO_GO`, unknown decisions, altered bound
digests, changed accession/version/URL, or relaxed caps before transferring a
biological byte. Each approved payload must be staged, size-checked, locally
SHA-256-hashed, reverified immediately before atomic read-only promotion, and
quarantined on mismatch. Source checksums, when present, are additional checks
rather than eligibility prerequisites. Exact H1/native-annotation linkage and
Tier3A pair evidence are revalidated from official staged reports.

### `run-repaired-vgp`

Edited and logged. It refuses invalid gates, altered gate/input/local payload
digests, incomplete acquisition, reference/annotation/pair mismatch, or cap
weakening with `NOT_SUBMITTED` and zero use. Callable/queryable denominators
and target totals are measured from post-alignment artifacts before result
acceptance. The definition explicitly prevents treating a VGP H1/H2 assembly
pair as population genotypes or demographic input and repeats every numeric
ceiling, including 256 GiB per job.

### `audit-vgp-demography`

Edited and logged. This fork is metadata/literature only. It distinctly
defines PSMC, MSMC2, SMC++, valid independent Ne, ecological/census evidence,
and circular/excluded estimates. VGP H1/H2 is not assumed to be a diploid
heterozygosity dataset, multiple independent genomes, a population sample, or
otherwise demographic-ready. Repository lookups are batched, cached,
rate-limited, retry-aware, resumable, and provenance-stamped. Bulk genomic
download, compute, and inference are forbidden.

### `review-repaired-vgp`

Edited and logged. The join independently checks both run/refusal evidence and
the metadata-only demography audit. It replays NO_GO/altered-digest refusal,
verifies immutable exact inputs and measured post-alignment denominators,
prevents zero-use refusal from becoming resource calibration, applies all
caps, and keeps PSMC/MSMC2/SMC++ and independent/circular Ne distinct.

### `synthesize-repaired-vgp`

Edited and logged. It preserves refusal outcomes, separates all demographic
method/evidence categories, and prevents an H1/H2 pair from being reported as
demographic genotypes. Any expansion, full-catalog operation, raw population
bulk acquisition, or demographic inference requires a new explicit numeric
authorization. Hypothetical resource estimates do not authorize execution,
and no ready executable tasks may be created.

## Authorization retained

The resolver, independent gate, acquisition, run, review, and synthesis
definitions all preserve these ceilings, with every stricter integrated cap
winning:

| Resource | Maximum without new authorization |
| --- | ---: |
| Species | 6 |
| Compressed inputs | 120 GiB |
| Scratch | 750 GiB |
| Core-hours | 1,500 |
| Concurrent species | 2 |
| Memory per job | 256 GiB |

No definition authorizes the full catalog, raw population bulk acquisition,
or demographic inference. The demography fork may identify and estimate
future data needs only. The synthesis may report non-executable estimates but
must request new numeric authorization before downstream task creation.

## Graph audit

The exact scientific edges observed after editing are:

```text
integrate-vgp-pilot
  -> quality-vgp-pilot-2
       -> repair-vgp-candidate
            |-> regate-vgp-pilot
            |     -> acquire-repaired-vgp
            |          -> run-repaired-vgp --.
            `-> audit-vgp-demography --------+-> review-repaired-vgp
                                                   -> synthesize-repaired-vgp
```

The edge assertions treat WG-generated `.assign-*` predecessors and `.flip-*`
successors as lifecycle edges, not scientific bypasses. Every lifecycle
dependency attached to a controlled task resolves to a real WG task. After
excluding those dot-prefixed lifecycle nodes, each controlled task has exactly
the domain predecessors shown above—no extra edge, missing fork arm, missing
join arm, phantom dependency, or direct path around `regate-vgp-pilot`.
Acquisition depends only on regating; compute depends only on acquisition.

## Automated assertions

The read-only validator
[`assert_vgp_pilot_repair_wg.py`](assert_vgp_pilot_repair_wg.py) checks:

1. presence of all eight controlled tasks (quality plus seven downstream);
2. exact scientific `after` edges and reverse controlled `before` edges;
3. resolution of every dependency, including lifecycle dependencies;
4. the repair fork and review join;
5. exclusive regate-to-acquisition-to-compute ordering with no GO bypass;
6. required catalog/checksum/denominator/reference/pair/refusal language;
7. the complete cap vector and stricter-cap rule; and
8. metadata-only, method-specific, and non-circular demography language.

Validation command:

```sh
python3 analysis/assert_vgp_pilot_repair_wg.py \
  --graph "$WG_PROJECT_ROOT/.wg/graph.jsonl"
```

The expected success sentinel is
`VGP_PILOT_REPAIR_WG_ASSERTIONS_OK`; any mismatch exits nonzero and lists each
failed invariant. The validator performs no external request and mutates no
WG or project state.
