# Comprehensive VGP evidence synthesis

**Evidence freeze:** 2026-07-18 UTC

**WG task:** `synthesize-vgp-program`

**Decision:** closed-world reconciliation complete; biological core and all gene-conversion effects remain fail-closed.

## Executive result

The robust result is an audit result, not a vertebrate-diversity effect. The pinned VGP Freeze 1 catalog contains exactly 716 rows (714 unique scientific names; the two duplicate-name occurrences remain separate release rows), of which 581 are released and 135 unreleased in the frozen mirror inventory. It yields 569 catalog-linked entries: 566 distinct non-self candidates and three self-links. None is eligible or complete under the reviewed gates. The mirror has 47,870 planned objects but zero verified or reused objects. Thus there are zero admitted callable-diversity estimates, zero PSMC trajectories, zero exact-annotation subsets, zero biological scale-out jobs, and no estimate of diversity-range compression.

Technical non-execution is not low diversity and is not evidence against heterozygosity. The analysis therefore classifies compression of callable diversity, demographic-shape contrasts, functional partition differences, and phylogenetically adjusted cross-species associations as **not identifiable**.

## Closed-world reconciliation

| Layer | Frozen count/state | Interpretation |
| --- | ---: | --- |
| Freeze 1 catalog | 716 rows; 714 unique names | Release-row membership is the sampling frame. |
| Frozen mirror | 581 released, 135 unreleased; 47,870 planned objects; 0 verified/reused | Inventory and transport evidence only. |
| Linked-pair ledger | 569 entries: 566 non-self, 3 self-links | Every entry has one QC/status row. |
| Pair disposition | 566 failed `UPSTREAM_SCALEOUT_NOT_AUTHORIZED`; 3 excluded `CATALOG_SELF_LINK_NOT_PAIR` | Operational dispositions, not biological classifications. |
| Confidence tier | 566 `UNASSIGNED`; 3 `X` | No Tier A/B/C biological result. |
| Ten-slot review | 0/10 primary passes; 6 alternates retained; 15 PASS, 11 FAIL, 2 FAIL_NOT_ESTIMABLE, 28 NOT_REACHED gates | `CONDITIONAL_GO` means bounded repair/re-pilot only, not scale-out. |
| Core / PSMC / exact annotation | 0 / 0 / 0 | All remain not run/not estimable. |

Figure 1 is `analysis/vgp_comprehensive_figure_closed_world.svg`; the exact paper table is `analysis/vgp_comprehensive_table_core.tsv`.

## Diversity, PSMC, functional partitions, and scaling

The biological sampling unit is one audited diploid individual represented by an exact same-individual H1/H2 pair. Windows and callable blocks quantify within-pair uncertainty; they do not manufacture independent individuals. A future cross-species analysis must use a versioned species tree and a phylogenetic mixed/PGLS or equivalent hierarchical covariance model, with pair/species as the outer unit and technical strata as covariates.

Assembly-derived diversity and PSMC are explicitly non-independent because they reuse the individual, H1/H2 differences, callable mask, and consensus. They must be represented as joint correlated outcomes or descriptive paired summaries with pair-clustered resampling. Same-pair PSMC must never be placed as an independent predictor or validator of same-pair diversity.

Annotation is not an eligibility gate. Exact accession/version and coordinate-dictionary binding may partition an already accepted core result, and partition effects remain nested within the pair. No non-annotated result may be removed retrospectively. Here no core pair passed, so the exact-annotation subset is empty for a pre-result reason.

The scenario ledger preserves both required scales. `UNSCALED_PRIMARY` is absent because core did not run. `SPECIES_SCENARIOS_REQUIRED_AFTER_REVIEW` is also absent: no mutation rate or generation time was selected. Future scenarios must cite species-specific sources, be reported separately with assumption bounds, and leave the unscaled trajectory primary.

## Four non-interchangeable gene-conversion branches

| Evidence branch | Execution state | Sampling unit | Current claim |
| --- | --- | --- | --- |
| Direct pedigree/gamete events | Executed metadata/table preflight; input gate | independent complete meiosis; events clustered within meiosis | D01 reconciles 13 tetrads/52 products and 58 published directional candidates (44 CO-associated, 14 NCO), all excluded. Zero events are admitted, so rates, tracts, CO association, and GC bias are not estimable. |
| Population allele-frequency-spectrum gBGC | **NOT_RUN/DESIGN_ONLY** | unrelated individuals within one population; LD blocks | No conforming execution manifest, frequency spectrum, or B estimate exists. |
| Historical phylogenetic substitution bias | Executed metadata preflight; input gate | branch by callable single-copy partition; chromosome/synteny blocks | H01/H02 freeze ten taxa but have zero verified sequence bytes, callable bases, or substitutions. No historical bias is estimable. |
| Non-allelic conversion among paralogs | **NOT_RUN/DESIGN_ONLY** | resolved copy nested in haplotype/individual and orthogroup | No copy-resolved execution manifest, tract, or homogenization estimate exists. |

The direct pilot's published-marker audit contains 45 S-resolved and 37 W-resolved markers plus one ambiguous/non-SNV marker. This is only a suggestive follow-up signal: linked markers and event candidates are not independent, all candidates fail raw/paralog/error gates, reciprocal detection is uncalibrated, and the registered power threshold is unmet. H1/H2 heterozygous WS/SW states cannot replace parent-of-origin or four-product direction. Historical WS/SW asymmetry, if later estimated, would be long-term substitution evidence only. Because the frozen historical panels may include semi-complete assemblies, any later model must preserve callable-gap accounting and the preregistered fragmentation sensitivity instead of treating missing sequence as absence of substitutions. Population B would remain a model-dependent frequency parameter, and paralog homogenization would remain non-allelic. Figure 2 (`analysis/vgp_comprehensive_figure_evidence.svg`) preserves these boundaries.

## Conclusion classes and claim boundary

The complete claim-by-claim ledger is `analysis/vgp_comprehensive_claim_ledger.tsv`.

- **Supported:** immutable release accounting, pair/QC reconciliation, and explicit non-execution/authorization states.
- **Bounded:** planning/resource and future scenario envelopes, never biological effects.
- **Suggestive:** the unvalidated published direct-candidate marker tally, solely as a follow-up target.
- **Design-only:** population SFS/B and non-allelic copy-aware programs.
- **Not identifiable:** callable-diversity compression, PSMC histories, functional effects, direct rates/bias, historical substitution bias, and cross-species associations.

No row supports extrapolating a direct single-species/plant or human rate across vertebrates. A future sensitivity could introduce an external rate only with provenance, taxonomic relevance, and separately labeled lower/base/upper bounds; it could not become an observed vertebrate result.

## Immutable lineage and reproduction

The final manifest binds 35 reviewed inputs, including core, review, mirror, specialized, design, and environment files, and binds every paper output by SHA-256. The final-analysis environment is GNU Guix channel commit `44bbfc24e4bcc48d0e3343cd3d83452721af8c36`; `analysis/guix/channels.scm` has SHA-256 `45c055cd1d9010a72eacbb720037a22bccb2d8d6891dbd11b5d66365f29b3a17` and `analysis/guix/manifest.scm` has SHA-256 `2fb05e87aa2ac45ce51d4dcf93b232cb98627f525adace98357629ee3f15720a`.

Reproduce without acquisition or biological submission:

```bash
./analysis/run_vgp_comprehensive_synthesis_guix.sh
./analysis/run_vgp_comprehensive_synthesis_guix.sh --tests
./analysis/run_vgp_comprehensive_synthesis_guix.sh --all-tests
```

The generator is intentionally a local, read-only reconciliation over reviewed inputs. It refuses any promoted core result, any same-pair independence claim, any annotation veto, any admitted direct candidate, any phylogenetic estimate without inputs, or any population/non-allelic execution-state promotion.
