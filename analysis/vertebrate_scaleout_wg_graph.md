# Proposed WG graph: all-vertebrate Tier 3 scale-out

Date: 2026-07-17 UTC  
Status: **proposal only; none of these nodes were created or dispatched by the planning task**

## Authorization rule

Every node that retrieves external content, stages biological assets, submits
Slurm work, raises concurrency, or runs a scientific analysis has a direct
authorization predecessor. An authorization node is a human-decision task,
not an automatic validation task. Its completion record must contain the
individual approver(s), immutable input row IDs, maximum download bytes,
storage/scratch quota, maximum jobs/CPUs/memory, concurrency and bandwidth
caps, expiry, decision timestamp, and the review packet hash. A script existing
in the repository does not satisfy an authorization node.

## Dependency graph

```text
plan-all-vertebrate-tier3 (this task; planning only)
  |
  v
A00 authorize-catalog-metadata-retrieval [HUMAN; <=1 MB catalog object only]
  |
  v
S00 retrieve-and-content-lock-vgp-freeze [network; one small object]
  |
  v
I00 build-row-level-vertebrate-inventory [metadata-only; bounded official API]
  |
  v
R00 review-inventory-and-guidance-reconciliation
  |\
  | +--> E00 capture-combined-pinned-guix-lineage [approved Guix host]
  |        |
  |        v
  |      R01 review-guix-derivations-and-truth-smoke
  |
  v
A10 authorize-pilot-asset-download [HUMAN; exact 8 rows and byte envelope]
  |
  v
S10 stage-read-only-pilot-assets [network; bounded concurrency]
  |
  v
P10 exact-reference-native-annotation-callability-preflight
  |                         
  +------------------------ R01
  |                          |
  v                          v
R10 review-pilot-preflight-and-resource-envelope
  |
  v
A11 authorize-pilot-compute [HUMAN; exact Slurm/I/O envelope]
  |
  v
X10 run-eight-species-stratified-pilot [Slurm; max Tier3A concurrency 2]
  |
  v
C10 collect-pilot-success-failure-and-telemetry [afterany]
  |
  v
Q10 independent-pilot-qc-and-budget-refit
  |
  v
R11 review-pilot-go-no-go [scientific + reproducibility + HPC storage]
  |
  v
A20 authorize-expansion-wave-1 [HUMAN; explicit rows and quota]
  |
  v
S20 stage-wave-1-assets -> P20 preflight-wave-1 -> X20 run-wave-1
  -> C20 collect-wave-1 -> Q20 independent-wave-1-qc -> R20 review-wave-1
  |
  v
A2N authorize-expansion-wave-N [HUMAN; one task per wave, never blanket]
  |
  v
S2N -> P2N -> X2N -> C2N -> Q2N -> R2N
  |
  v
R30 review-full-catalog-readiness
  |
  v
A30 authorize-full-vertebrate-execution [HUMAN JOINT APPROVAL]
  |
  +--> X30 execute-remaining-tier3a-catalog
  |      |
  |      v
  |    C30 collect-tier3a-catalog -> Q30 independent-tier3a-qc
  |
  +--> X31 execute-remaining-tier3c-catalog
         |
         v
       C31 collect-tier3c-catalog -> Q31 independent-tier3c-qc
             \                 /
              v               v
                R40 reconcile-complete-catalog
                          |
                          v
                A40 authorize-final-scientific-analysis [HUMAN]
                          |
                          v
                X40 fit-prespecified-phylogenetic-models
                          |
                          v
                Q40 independent-results-reproduction
                          |
                          v
                A41 authorize-release-and-cleanup [HUMAN]
                          |
                          v
                R41 release-tables-and-checksum-led-cleanup

R40 review-complete-catalog
  |
  v
A50 authorize-demographic-subset-planning [HUMAN; separate program]
  |
  v
D50 method-specific-demographic-inventory-and-plan
  |
  v
A51 authorize-expensive-demographic-execution [HUMAN; if ever justified]
```

`A30` does not authorize `A50` or `A51`. Assembly/composition completion cannot
implicitly start population-genomic or demographic inference.

## Node contracts

| node class | input gate | required output | cannot do |
|---|---|---|---|
| `A*` authorization | signed review packet with quantitative thresholds | immutable decision JSON/message reference and approved envelope | run scripts, infer approval from predecessor success |
| `S00` catalog retrieval | A00 | source URL/revision/hash/UTC record and frozen raw TSV | stage assemblies or follow moving `main` silently |
| `I00` inventory | frozen TSV | complete row ledger, exact reported-count reconciliation, blockers | download full FASTA/GFF catalog |
| `E00/R01` environment | inventory review | channels/manifests/derivations/closure, reproducible binary hashes, truth smoke | introduce Conda/container/module tools |
| `S*` asset staging | corresponding authorization | read-only exact-accession assets, checksums, licenses, failure ledger | replace unavailable version with latest/alternate |
| `P*` preflight | staged manifest | row fingerprints, exact-reference/native-annotation/denominator verdict, resource estimate | waive hard scientific gates |
| `X*` execution | corresponding authorization and passing preflight | attempt-scoped atomic outputs and sentinels | consume legacy/wildcard/partial artifacts |
| `C*` collection | `afterany` on execution | closed-world accepted/rejected/unavailable/failed ledger and retry manifest | drop failed elements |
| `Q*` QC | complete collection | checksum/scientific/telemetry comparison and prediction error | mutate primary artifacts |
| `R*` review/release | passing quantitative packet | go/no-go decision or frozen release | broaden authorized scope |

## Quantitative readiness attached to graph edges

- `R00 -> A10`: every pilot row has a resolved TaxId, exact H1/H2 where
  required, native annotation candidate, evidence URL/date, license review,
  size estimate, and explicit acceptance/rejection; all guidance counts are
  reconciled.
- `R01 -> R10`: exact Guix derivations and closure captured; two-build SweepGA
  identity; IMPG/submodule identity; controlled and H1-native truth smokes
  pass with no ambient executable.
- `R10 -> A11`: exact pilot rows and input hashes frozen; proposed download,
  storage, scratch, memory, concurrency and I/O envelope fits quota.
- `R11 -> A20`: eight-row pilot passes every reference/annotation/denominator
  gate; no more than one first-attempt technical failure; no unexplained
  failure; median prediction error <=35%, p95 <=75%; memory/scratch/I/O below
  80% of allocation/caps; telemetry completeness 100%.
- `R20/R2N -> A2N`: cumulative technical failure <=5%, unexplained <=1%; median
  resource error <=25%, p95 <=50%; lineage/checksum violations zero; closed
  ledger 100%; no storage incident.
- `R30 -> A30`: all blockers resolved; retry backlog zero; 25% quota headroom;
  independent QC passes; full row list and statistical/phylogeny plans frozen.
- `R40 -> A40`: accepted/rejected/unavailable/failed counts exactly reconcile
  to inventory and no partial/legacy/superseded artifact is consumable.

## Instantiation policy

The synthesis task may translate the symbolic graph into `wg add` commands,
but it must preserve every `A*` node as a direct dependency of the external or
expensive node it controls. It must not create an authorization node already
marked done, auto-approve one through a shell command, or place an executable
node on an alternate dependency path that bypasses authorization. Expansion
waves receive distinct IDs; a successful first wave is evidence for asking
permission, not permission for all later waves.
