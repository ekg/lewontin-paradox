# Corrected VGP evidence and claims

Task: `correct-vgp-evidence` · schema `vgp-evidence-correction-v1`

This publication supersedes the prior two-result synthesis (`analysis/vgp_real_synthesis_v1`) with the independently reviewed clean-canary control and bounded three-pair evidence (`analysis/vgp_three_pair_independent_review_v1.md`). Historical artifacts are preserved with explicit supersession links; nothing was silently rewritten.

## Corrected core results

| Pair | Species | Heterozygous SNPs | Callable bp | pi | PSMC theta (95% CI) | Annotation |
|---|---|---:|---:|---:|---|---|
| P07 | Spinachia spinachia | 316,631 | 270,531,638 | `0.0011704028495181033` | 0.033275 [0.030717, 0.038811] | exact native (bounded callset) |
| P03 | Colius striatus | 1,632,584 | 875,683,638 | `0.0018643536651303744` | 0.040367 [0.025469, 0.056461] | missing, not zero |
| P02 | Pseudorca crassidens | 4,195,014 | 1,707,746,195 | `0.0024564622145154306` | 0.061674 [0.055957, 0.068648] | missing, not zero |

The prior P07 whole-genome/global-lace value (574,122/267,379,237 = `0.0021472198306856562`, reproduced byte-identically by the clean canary) is retained as historical diagnostic evidence and reproduction control, not as an admissible core result. The review's fresh bounded reproduction reproduced the P07 range `CM106587.1:0-5000000` key-for-key (3871 keys, SHA-256 `b7f64cd018a693ade59de7daec49b66d9f57e6ea7baa888c6e9e1c7099da1405`, matches clean IMPG subset: true, global graph not materialized) against the like-for-like clean subset.

P04 remains a retained prior-lineage result (pi `0.004604184795871289`) with raw validation still pending; it is not part of the reviewed bounded core.

## P07 reclassification: unresolved read/assembly discordance

P07 is reclassified from *assembly-invalid* (`concrete_haplotype_reconstruction_failure`) to **unresolved read/assembly discordance**. The reclassification is required because:

1. the prior label was derived from H1-only Illumina/HiFi evidence compared against the *rejected* clean/global-lace callset (574,122 SNPs), not the accepted bounded callset (316,631 SNPs);
2. no symmetric H1/H2/graph-representation comparison exists for any accepted result (review decision: read sensitivity **INCONCLUSIVE**; defect I-1);
3. no concrete sequence, provenance, or projection defect has been demonstrated for the bounded result — provenance, graph-ID resolution, REF/ALT reconstruction, and mask accounting all closed with zero failures in the independent audit.

The H1-only diagnostics (read/assembly pi ratio 0.436332 on the DP10-80 common mask; homozygous-reference contradictions 50.122% Illumina / 53.166% HiFi; k-mer heterozygosity 0.0026652671) are retained as diagnostic uncertainty (review A-2). They may reflect a concrete assembly or projection error, but the required symmetric test has not localized one, so they do not exclude the bounded result. The whole-genome/bounded pi discrepancy (review A-1) remains an open assembly-confidence caveat.

## Historical failures are pipeline defects, not biological exclusions

Every historical failure family — FastGA execution, IMPG global architecture and graph/FASTA identifiers, IMPG lace threading, compression/transfer integrity, sequence-dictionary binding, annotation discovery, and scratch containment — is reclassified as a **pipeline defect**. P02 and P03 prove the correction concretely: both were previously `HARD_INVALID_PRIMARY` and both now carry accepted, reviewed bounded results. Species without reviewed reruns (P01, P05, P06, P08, P09, P10) remain eligible; no estimate is imputed for them and none is reported as zero.

| Defect family | Affected | Corrected classification | Resolution |
|---|---|---|---|
| fastga_execution | P03;P08 | PIPELINE_DEFECT_NOT_BIOLOGICAL_EXCLUSION | resolved_in_bounded_lineage |
| impg_global_architecture | P02;P03 | PIPELINE_DEFECT_NOT_BIOLOGICAL_EXCLUSION | resolved_in_bounded_lineage |
| impg_lace_threading | P07;P03;P02 | PIPELINE_DEFECT_NOT_BIOLOGICAL_EXCLUSION | mitigated_two_thread_minimum |
| compression_integrity | P07;annotation_ingest | PIPELINE_DEFECT_NOT_BIOLOGICAL_EXCLUSION | resolved_quarantine_and_decode_fix |
| sequence_dictionary_binding | P02;P03;P04 | PIPELINE_DEFECT_NOT_BIOLOGICAL_EXCLUSION | partially_resolved |
| annotation_discovery | P02;P03 | PIPELINE_DEFECT_NOT_BIOLOGICAL_EXCLUSION | discovery_resolved_partitions_not_computed |
| scratch_containment | P03;P07;P08 | PIPELINE_DEFECT_NOT_BIOLOGICAL_EXCLUSION | resolved_private_scratch_roots_and_guards |

### Pipeline reliability metrics

- Prior scale packet: 658 scheduler allocations — 242 completed, 46 failed, 348 cancelled, remainder nonterminal at freeze (states, not biology).
- Bounded wave: 4 completed allocations yielding 3 accepted results (P02 accepted lineage additionally consumed one failed-but-retained job and one zero-query resume); cancellations of the prohibited global chains are technical provenance.
- Accepted-result quality gates: zero cross-range duplicate keys, zero boundary ownership failures, zero unowned or multiply owned callable bases, zero unresolved graph IDs, zero REF/ALT reconstruction failures, and 600/600 finite centered PSMC bootstraps across the trio.

## Exact-accession annotation catalog

The comprehensive exact-accession annotation catalog is published as the authoritative annotation-path table: 581 Freeze 1 assemblies with preferred exact-dictionary annotations, 1833 accepted parsed annotations (1834 physical objects, 15,879,123,133 bytes), and 10 pilot reference bindings (`analysis/vgp_annotation_catalog.tsv`, `..._assembly_bindings.tsv`, `..._pilot_bindings.tsv`, validation PASS).

**Binding semantics are explicit.** A catalog `EXACT_DICTIONARY` pilot binding binds the annotation to *its own reference assembly* (the annotation's native GCF/GCA target) — it does not by itself certify dictionary equality with a run's H1 assembly. Only P07's annotation is exact-native to an accepted run's H1 with computed biological partitions (P08 also has an exact-native accession binding but no accepted run). P02/P03 bindings are reference-assembly bindings whose sequence-dictionary audit against the run H1 has not been performed; their partitions are *unavailable, not zero* (review D-3). The per-pilot semantics are in `annotation_paths.tsv`.

## What has and has not been scaled

**Scaled (catalog/architecture level):**

- annotation catalog reconciliation: 581 assemblies, 1833 accepted annotations, 10 pilot reference bindings;
- bounded execution architecture across three pairs spanning 0.41-2.67 Gbp H1 assemblies, with range-local IMPG queries, strict H1 validation, and closed ledgers.

**NOT scaled (biological level):**

- biological inference beyond the three accepted bounded pairs (plus the retained prior-lineage P04 result);
- biological annotation partitions beyond P07;
- population sampling of any species;
- any full VGP scale-out. **No scale-out was launched by this correction.** The review's broad biological scale-out NO-GO stands until the symmetric representation test (I-1), resource telemetry (I-2/I-3), and assembly-confidence handling (A-1..A-4) are addressed and an explicit review GO is issued.

**Catalog reconciliation is not a biological scale-out** and must not be described as one.

## Claim classifications

| Claim | Prior | Corrected |
|---|---|---|
| CORE-BOUNDED-TRIO | supported | supported |
| TRIO-OBSERVED-RANGE | supported | bounded |
| P07-DISPOSITION | supported | reclassified_unresolved_discordance |
| P07-DISCORDANCE-EVIDENCE | bounded;supported;suggestive | retained_as_diagnostic_uncertainty |
| P04-DISPOSITION | bounded | unchanged_retained_raw_pending |
| PIPELINE-DEFECTS-NOT-EXCLUSIONS | (failures previously encoded as X hard-invalid primary / X execution error) | supported |
| ANNOT-CATALOG-SCALE | supported | supported_catalog_level_only |
| ANNOT-EXACT-BOUNDED | supported | supported_superseded_values |
| SCALE-STATUS | (new) | supported |
| PSMC-UNSIGNED | supported;bounded | supported_as_unscaled_trajectories |
| POPULATION | unidentifiable | unidentifiable |
| VERTEBRATE-RANGE | unidentifiable | unidentifiable |
| LR-IMPLICATION | unidentifiable | unidentifiable |
| GENE-CONVERSION | unidentifiable | unidentifiable |

Full statements, evidence, and forbidden inferences: `claims_ledger.tsv`.

## Supersession

Historical artifacts are preserved **byte-identically**: the prior synthesis packet and the repair-base status table bind them by SHA-256, and this correction verifies every binding rather than editing any file (adding banners would silently invalidate those frozen evidence digests):

- 33 historical artifacts verified digest-preserved against their frozen bindings.
- `analysis/vgp_real_synthesis_v1/report.md`: digest_preserved_synthesis_manifest
- `analysis/vgp_real_synthesis_v1/paper_pairs.tsv`: digest_preserved_synthesis_manifest
- `analysis/vgp_real_synthesis_v1/claim_ledger.tsv`: digest_preserved_synthesis_manifest
- `analysis/vgp_real_canary_report_v1.md`: preserved_unmodified
- `analysis/vgp_clean_canary_report_v1.md`: preserved_unmodified

The authoritative supersession map is `supersession_ledger.tsv` (superseded artifact -> superseding artifact, with reasons). No historical file content was rewritten, reformatted, or annotated.

## Forbidden inferences (unchanged in kind)

- One phased H1/H2 pair per individual is not population heterozygosity (B-1).
- PSMC scenario grids are not calibrated species demography (B-2).
- No cross-species annotation contrast is supported (B-3).
- No vertebrate diversity range or Lewontin-paradox test is identifiable from this evidence.

## Evidence identities

All input digests are recorded in `manifest.json`.
