# Independent review of the repaired VGP three-pair pilot

Review date: 2026-07-29 UTC

Task: `review-vgp-three-pair`

Verdict: **three bounded core assembly results accepted; annotation accepted for P07 only; read sensitivity inconclusive for the repaired results; broad biological scale-out not yet ready**

## Scope and acceptance boundary

The accepted core is the bounded-production lineage for P07, P03, and P02.
The clean P07 whole-genome run is a required reproduction/control, but it is
not an architecture-admissible core result: it enumerated 203,698 regional
queries and hierarchically laced them into a global VCF. Its whole-genome
574,122/267,379,237 result remains historical diagnostic evidence. The
accepted P07 core is the later bounded result, 316,631/270,531,638.

I reran `analysis.audit_vgp_bounded_results` in the pinned Guix environment
against the three durable output trees. The audit rehashed 16,133 promoted
files totaling 12,179,399,820 bytes and reproduced
`analysis/vgp_three_pair_execution_v2.json` byte-for-byte, SHA-256
`f57b8dc2fcfd90eb63a5b1c03fb5e2193f18b87b7b008d4d29749f88ea7b8c15`.
This was not a summary-only check: it reread every range VCF, the callable
BEDs, mask inputs, consensus blocks, primary PSMCFA, 600 bootstrap fits,
mapping PAFs, graph dictionaries, and the exact P07 GFF products.

## Fresh bounded P07 reproduction

I issued a new query from the retained P07 alignment index for the frozen
H1 half-open interval `CM106587.1:0-5000000`. Only the 2,500 native
partitions owned by that range were emitted. The local query created 2,500
regional VCFs, locally laced and H1-normalized them, and produced 3,871
unique `(CHROM,POS,REF,ALT)` keys. Those keys were byte-identical to the
like-for-like clean IMPG subset:

- normalized-key SHA-256:
  `b7f64cd018a693ade59de7daec49b66d9f57e6ea7baa888c6e9e1c7099da1405`;
- callable bp from the independently intersected clean mask: 3,116,900;
- query/lace elapsed: 191/1 seconds;
- peak local calls+temporary state: 4,216,589 bytes;
- calls and temporary directories absent after explicit disposal.

The reproduction used Guix channel
`44bbfc24e4bcc48d0e3343cd3d83452721af8c36`, closure
`8fcdb32021f1cd8eac839509cff47ab6bdd63b656b30e243fdf78d3c4ba24f9d`,
IMPG executable
`c587dc2326cd24f887b1fcb3938404229ad0f0a27ef0773e90c287b1ade160d4`,
and bcftools
`79637872d29b03be83293a56d297519ceb74d861ff61888966d5e17157f57dd4`.
The full record is
`analysis/vgp_three_pair_review_canary_reproduction_v1.json`.

## Recomputed core measurements

| Pair | H1 bp | Ranges (query/nonquery) | Normalized records | Callable bp | Callable SNPs | Recomputed pi |
|---|---:|---:|---:|---:|---:|---:|
| P07 | 407,561,107 | 109 (106/3) | 572,450 | 270,531,638 | 316,631 | 0.0011704028495181033 |
| P03 | 1,212,464,869 | 407 (329/78) | 2,081,450 | 875,683,638 | 1,632,584 | 0.0018643536651303744 |
| P02 | 2,674,608,476 | 930 (925/5) | 5,356,783 | 1,707,746,195 | 4,195,014 | 0.0024564622145154306 |

Each quotient was independently recomputed from callable SNP and denominator
counts. The range plan is disjoint and exhaustive on every H1 contig.
Unaligned sequence is represented by 3, 78, and 5 explicit nonquery ranges,
not silently omitted. Native H1 partition rows have exactly one owner:
203,698 (P07), 603,972 (P03), and 1,375,084 (P02). There are zero
cross-range duplicate normalized keys, zero boundary ownership failures,
zero unowned callable bases, and zero multiply owned callable bases.

The ordered mask reconstruction closes with discrepancy zero for all pairs.
Final callable sequence is exactly 270,531,638, 875,683,638, and
1,707,746,195 bp. Consensus non-`N` sequence equals the callable
intersection in independently selected early, middle, and late ranges.
No noncallable base is encoded as homozygous reference. Heterozygous SNPs use
IUPAC; indel flanks are masked under the frozen 10-bp rule.

Exact CIGAR reconstruction and H1 REF validation pass:

| Pair | Retained PAF rows | Reverse-strand rows checked | Unresolved graph IDs | REF/ALT or coordinate failures |
|---|---:|---:|---:|---:|
| P07 | 4,297 | 967 | 0 | 0 |
| P03 | 394 | 207 | 0 | 0 |
| P02 | 55,907 | 26,584 | 0 | 0 |

The whole-genome VCF/BCF in each accepted tree is a convenience
`bcftools concat` over normalized, nonoverlapping range VCFs. It is not an
IMPG query or lace product.

## IMPG architecture audit

**Pass for all three accepted core results.** The retained whole-assembly
objects are an alignment index and native partition table, which is the
permitted IMPG role. The production runner emits one requested range BED,
invokes `impg query` on that BED, locally laces only that range's VCF list,
normalizes and audits it, and removes the local work directory
(`analysis/slurm/run_vgp_bounded_pair.sh:216`,
`:224`, `:233`, `:308`). Only after every range is closed does it concatenate
the normalized range products (`:398-407`).

The closed ledgers inventory every accepted durable file. They contain the
index, partition table, frozen H1 range plan, range-local normalized
products, consensus blocks, and reductions; they contain no whole-genome
graph, global IMPG-lace output, or global partition-assignment ledger.
Peak local graph state was only 7,570,025, 12,836,941, and 14,389,589 bytes.
The prohibited jobs 1797004-1797008 and 1797029-1797033 were canceled and
are excluded from inference. P02 resume 1804979 reused only the already
bounded 1804024 range products, performed zero IMPG queries, and
strict-H1-revalidated every VCF/BCF.

The clean P07 global-lace products and all earlier exhaustive/global
artifacts are rejected as core inputs. The fresh reproduction used the clean
index and local clean subset only as an equivalence comparator; neither the
clean global laced VCF nor its whole-genome primary callset entered any
accepted result.

## Annotation partitions

**Decision: accept P07; unavailable, not zero, for P02/P03.** The P07 GFF is
the exact H1-native `GCA_048126635.1-GB_2025_08_04` object, SHA-256
`8f640543accd8081d1b7048eda32c9f1eef33b02f321b7b0f8adcf3b01dd6838`.
Its 38-sequence region dictionary equals H1 exactly. The annotation code
parses the GFF features, canonical transcripts, CDS frames, and fourfold
sites and queries them by exact H1 interval against the callable mask and
bounded-range-derived normalized variants
(`analysis/vgp_real_canary_annotation.py:97-148`). It does not query a
global graph or a global-lace output; no additional graph query is needed
because the disjoint bounded tiling already supplies the normalized local
records.

| Partition | Heterozygous records/SNPs | Callable bp | Estimate |
|---|---:|---:|---:|
| CDS | 7,603 | 26,152,171 | 0.0002907215618925098 |
| fourfold | 2,212 | 4,225,276 | 0.0005235160969366262 |
| fourfold W | 987 | 1,452,971 | 0.0006792977974095835 |
| fourfold S | 1,225 | 2,772,305 | 0.0004418705734037200 |
| W-to-S | 868 | 1,452,971 | 0.0005973966445304139 |
| S-to-W | 920 | 2,772,305 | 0.0003318538183929979 |

The 542 frame-discordant overlap positions are explicitly excluded. P02 and
P03 have no cataloged exact-native GFF; no annotation estimate is inferred
or represented as zero.

## PSMC and repaired centering

**Pass as unscaled assembly-derived trajectories.** Primary PSMCFA
populations were independently recounted and contain only `N/K/T`.
Range blocks reconstruct the full contig-aware PSMCFA exactly and never
cross contig boundaries. All 200 fits per pair are nonempty, finite, and
have the expected 64 intervals.

| Pair | Primary theta | Recomputed nearest-index 95% interval | Centered |
|---|---:|---:|---|
| P07 | 0.033275 | [0.030717, 0.038811] | yes |
| P03 | 0.040367 | [0.025469, 0.056461] | yes |
| P02 | 0.061674 | [0.055957, 0.068648] | yes |

The repaired diagnostic uses the final native-iteration theta, requires all
200 finite outputs, and computes bounds as
`sorted[round((n-1)*0.025)]` and `sorted[round((n-1)*0.975)]`. Scaling
scenarios remain separate from the preserved primary trajectory.

## Mapping fallback equivalence

P03's controlled chromosome-1 comparison used byte-identical staged FASTAs.
FastGA and WFMASH both have zero coordinate, strand-transform, or exact
REF/ALT reconstruction failures. Common target coverage is 17,317,077 bp
and target-coverage Jaccard is 0.9706948, but exact raw-variant Jaccard is
only 0.2529223 (3,700 shared variants). WFMASH is therefore accepted as
reproducible fallback infrastructure, not claimed biologically equivalent
to FastGA at the variant level. Backend gap-placement sensitivity remains
an assembly-confidence limitation.

## Read sensitivity

**Decision: inconclusive for the repaired bounded results.** Available P07
Illumina/HiFi work maps reads only to H1 and compares them with the rejected
clean/global-lace P07 callset (574,122 SNPs), not the accepted bounded P07
callset (316,631 SNPs). It has no symmetric H2-reference analysis and no
graph-projected read comparison. Consequently its prior
`concrete_haplotype_reconstruction_failure` label cannot be transferred to
the repaired bounded P07 result.

The old H1-only evidence remains important diagnostic evidence: on the
DP10-80 common mask, read/assembly pi ratio was 0.436332, and
depth-qualified homozygous-reference contradiction fractions were 0.501223
(Illumina) and 0.531662 (HiFi). These numbers may reflect a concrete
assembly or projection error, but the required symmetric representation
test has not localized one. Reads are not treated as a generic truth oracle
and do not override the assembly-primary estimates. P02/P03 have no
accepted-result read sensitivity evidence.

## Provenance and scratch

Input accessions, whole-FASTA digests, sequence dictionaries, mapping PAFs,
GFF identity, graph IDs, normalized outputs, and closed result trees are
checksum-bound. P07 uses H1 `GCA_048126635.1` and H2
`GCA_048127205.1`; P02/P03 whole FASTAs were rehashed against their frozen
input manifests. All graph IDs resolve to the staged H1/H2 dictionaries,
with no alias substitution and no silent omission.

The accepted runners create private `/scratch/vgp-<pair>-bounded-<job>-*`
roots, direct `TMPDIR/TMP/TEMP` and graph temporary directories below them,
copy only closed outputs to durable storage, and prefix-check before
recursive cleanup. The fresh canary independently confirmed local
temporary disposal. However, bounded production did not retain a live
resolved-path/open-file guard or scratch high-water telemetry comparable to
the clean canary. Containment is supported by command paths, cleanup guards,
the common octopus02 node, and absence of durable local graph artifacts, but
broader execution must add per-stage resolved-path and filesystem telemetry.

The scientific environment is pinned to GNU Guix channel
`44bbfc24e4bcc48d0e3343cd3d83452721af8c36`, profile
`/gnu/store/3c2mxm30rbzvnw7qsi235mrkk3m38fym-profile`, and closure
`8fcdb32021f1cd8eac839509cff47ab6bdd63b656b30e243fdf78d3c4ba24f9d`.
bcftools and PSMC are store objects. IMPG, WFMASH, and SweepGA/FastGA are
commit- and digest-pinned executables outside the profile; this is exact
byte provenance, but not a fully Guix-built toolchain.

## Decisions

| Domain | Decision | Meaning |
|---|---|---|
| Core assembly inference | **ACCEPT** | Three bounded H1/H2 assembly-derived pi, consensus, and unscaled primary PSMC results pass all implementation/provenance gates. They are per-individual assembly measurements, not population estimates. |
| Annotation partitions | **ACCEPT P07 ONLY** | Exact native GFF/dictionary and exact H1 interval intersections pass. P02/P03 are unavailable and remain missing. |
| Read sensitivity | **INCONCLUSIVE** | Existing P07 evidence is H1-only and targets the rejected clean callset; no symmetric H1/H2/graph comparison exists for accepted results. |
| Broader scale-out | **NO-GO for biological inference; conditional GO for a bounded technical wave** | The bounded architecture is technically ready for a limited telemetry-gated wave. Broad biological use waits for symmetric representation sensitivity, measured RSS/scratch/CPU telemetry, and explicit assembly-confidence handling. |

## Remaining defect ledger

| ID | Class | Severity | Remaining defect and disposition |
|---|---|---|---|
| I-1 | implementation | blocking broader biological scale-out | No read comparison is symmetric across H1, H2, and graph representations, and none targets the repaired bounded P07 callset. Implement a frozen common-coordinate symmetric analysis; do not transfer the old failure label. |
| I-2 | implementation | blocking resource extrapolation | Slurm `MaxRSS`, disk I/O, and CPU fields are empty for accepted jobs; bounded scratch high-water is absent. Add cgroup/per-stage telemetry before resizing or large fan-out. |
| I-3 | implementation | moderate | Bounded production uses `/scratch` prefix checks but lacks the clean canary's live resolved-path/open-file guard and retained filesystem-type proof. Add `realpath`, filesystem-type, capacity, and managed-open-path telemetry. |
| I-4 | implementation | moderate | Pinned IMPG lace does not progress with one thread; two threads are a required operational minimum. Preserve the guard and report this as implementation behavior. |
| I-5 | implementation | resolved in accepted lineage | P02 exposed graph-allele REF/ALT orientation and a post-reconstruction `ALT=.` record. Accepted resume filtering and strict-H1 validation close the result; failed jobs remain provenance. |
| D-1 | data provenance | moderate | IMPG, WFMASH, and SweepGA/FastGA are digest-pinned external executables, not Guix store builds. Package the accepted binaries/derivations in Guix before claiming a fully Guix-built closure. |
| D-2 | data provenance | low | Some copied P07 graph-ledger metadata records an obsolete `/scratch/...` source path although the durable digest-identical ledger exists. Replace ephemeral source-path fields with durable object identities in future runs. |
| D-3 | data provenance | limiting annotations | P02/P03 lack cataloged exact-native annotations. Preserve missingness; acquire dictionary-identical GFFs before annotation scale-out. |
| A-1 | assembly confidence | high | P07 clean/global-lace pi (0.0021472) differs materially from bounded pi (0.0011704). The bounded architecture is accepted and the clean whole-genome result rejected, but the source of all discarded global differences is not range-by-range localized. |
| A-2 | assembly confidence | high | Old P07 H1-only reads show severe discordance, but do not symmetrically localize an H1, H2, or graph projection error in the bounded result. Retain as diagnostic uncertainty. |
| A-3 | assembly confidence | high | P03 controlled FastGA/WFMASH exact-variant Jaccard is 0.2529 despite high coverage overlap. Treat backend-specific gap placement as method sensitivity. |
| A-4 | assembly confidence | moderate | P02/P03 have no accepted-result raw-read sensitivity evidence. Assembly estimates remain primary but externally uncalibrated. |
| B-1 | biological limitation | fundamental | One phased H1/H2 pair measures within-individual assembly divergence; it is not population heterozygosity or a species mean. |
| B-2 | biological limitation | fundamental | PSMC mutation-rate/generation-time grids are generic sensitivity scenarios, not calibrated absolute demographic histories. |
| B-3 | biological limitation | moderate | Only P07 has exact-native CDS/fourfold partitions, so cross-species annotation contrasts are not supported. |

## Corrected resource handoff

The measured resource record is
`analysis/vgp_three_pair_review_resource_model_v1.json`. Successful
allocations were 16 CPU/96G for P07 (4:31:36), 16 CPU/128G for P03
(14:20:08), and 32 CPU/256G for the P02 resume (17:34:02). P02's accepted
lineage also consumed 5:25:55 in job 1804024; its combined accepted-lineage
cost is 22.999 hours and 735.973 allocated CPU-hours. Requests are not
usage measurements. The only direct production RSS observation is
13,777,560 KiB for serial P02 consensus.

For planning only, successful elapsed time is 11.107 h/H1-Gbp (P07),
11.823 h/H1-Gbp (P03), and 6.568 h/H1-Gbp for the reuse-heavy P02 resume.
Use 11.823 h/H1-Gbp as a conservative observed wall-time coefficient for a
limited wave, but do not extrapolate RAM or scratch from the 7.6-14.4 MB
local graph peaks. Retain successful requests as admission ceilings until
per-stage RSS, scratch, and CPU telemetry are captured. Separate technical
recovery cost from steady-state cost in every budget.
