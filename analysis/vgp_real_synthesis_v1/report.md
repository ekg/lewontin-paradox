# Real VGP biological evidence synthesis

Generated from the closed-world upstream packets completed 2026-07-21T18:10:15Z
WG task: `synthesize-vgp-real`
Canonical shared VGP root: `/moosefs/erikg/vgp`

## Paper-oriented result

The real VGP execution produced two assembly-derived individual diversity estimates, not zero biological estimates. P07 (*Spinachia spinachia*) has pi = 0.0021472198306856562 from 574,122 differences across 267,379,237 callable bp; P04 (*Falco naumanni*) has pi = 0.004604184795871289 from 3,548,818 differences across 770,780,965 callable bp. Their observed assembly-derived range is therefore 0.002147–0.004604 per bp, or 2.144254-fold.

That computational range is not a validated vertebrate biological range. Exact-individual Illumina and HiFi validation makes P07 a `concrete_haplotype_reconstruction_failure`, so its pi and PSMC are preserved for provenance but excluded from quantitative biological use. P04 remains the sole retained assembly-derived estimate, with its exact CLR validation run pending. Thus there are two completed assembly observations, one retained/raw-pending observation, one read-invalidated observation, and zero fully raw-validated admitted quantitative pairs.

## Execution and confidence accounting

The immutable Freeze 1 accounting covers all 716 catalog rows and all 569 links. Ten links belong to the frozen same-individual, mutually comparable audit roster. Across the pilot and scale packet, all 650 scheduler allocations are retained: 240 completed, 42 failed, 346 cancelled, 21 pending, and one running. The raw-read packet adds eight distinct allocations: two completed, four failed, and two cancelled. Nothing missing or nonterminal is converted to a biological zero.

| Synthesis tier | Pair(s) | Meaning |
|---|---|---|
| T1 retained/raw pending | P04 | Completed assembly pi and PSMC retained; exact raw validation pending |
| T2 invalidated by exact reads | P07 | Artifact preserved, but pi, PSMC, and annotation partitions excluded from quantitative biological synthesis |
| X hard-invalid primary | P01, P02, P03, P05 | Invalid primary execution; no estimate |
| X execution error | P06, P09, P10 | Concrete execution error; no estimate |
| P resumable at freeze | P08 | Running wave at cutoff; no frozen estimate |

The machine-readable pair table is `paper_pairs.tsv`; the complete 658-record scheduler reconciliation is `job_ledger.tsv`.

## Raw-read validation and assembly sensitivity

At the P07 primary DP10-80 common mask, 255,821,332 bp support assembly pi = 0.000872374474229 and mapped-read pi = 0.00038064456642, a read/assembly ratio of 0.436332. Shared calls yield a lower concordant bracket of 0.000353137087098; the union yields an upper bracket of 0.00089988195355. These are paired-caller disagreement bounds, not complete biological confidence intervals.

The primary mask retains 95.677% of inherited callable bases but 38.872% of inherited assembly differences. The excluded assembly-difference density is 34.8067-fold the retained density, a strong mappability/collapse/callability sensitivity signal.

More decisively, depth-qualified homozygous-reference contradictions are the majority in both technologies: Illumina 50.122% (Wilson 95% 49.915%–50.330%) and HiFi 53.166% (Wilson 95% 52.959%–53.373%). These satisfy the predeclared reconstruction-failure rule. The k-mer QV of 35.9217 does not meet the separate severe-sequence-error rule (QV <20), which shows why consensus quality and haplotype reconstruction must remain distinct.

P07's negative-binomial k-mer spectrum gives model-based heterozygosity 0.0026652671, 1.24126 times inherited assembly pi but about 7.002 times stringent common-mask read pi. This disagreement is reported as method sensitivity; estimates are not averaged.

P09's one-cell 0.810x HiFi run is only a mapping compatibility control. P04's 42,344,746,693-byte exact CLR run remains planned. Four canonical raw objects totaling 31,058,137,613 bytes were independently rehashed and reused. One full-size corrupted P07 R2 transfer (4,431,902,981 bytes; SHA-256 `defaed9e929d8acf9d58006a3be51c26dd7c1937079f7d47e9d818420475c965`) remains quarantined; a clean retry produced canonical SHA-256 `c542f6efd9fc1d8f557c89629743fbe4a39584f24002c84185e479057cb443ac` without redownload of verified objects.

## Demographic histories

Each completed pair has a 64-interval assembly-derived PSMC trajectory and 200/200 finite block bootstraps. P04 theta0 is 0.52445 per 100-bp bin (bootstrap 0.435533–0.618109); its unscaled lambda rises to 2.874536 at interval 4 and falls to a broad later minimum of 0.492194 at interval 40. This is a descriptive relative-history shape only.

P07 theta0 is 0.032921 (bootstrap 0.028314–0.038858) and its inherited trajectory reaches lambda 159.317994; exact-read PSMC has a similar theta ratio (0.940281) but strongly discordant lambda shape (Pearson r = -0.144006). The time grids and data sources overlap, so this is sensitivity, not replication, and the inherited P07 history is invalidated for biology.

For both pairs, absolute histories exist only as nine generic mutation-rate × generation-time scenarios. Across those grids, P04 spans 32266.392912500003–753775.20259999996 in scenario Ne and 0–16900317.338 scenario years; P07 spans 2671.7860575–2622453.840237 and 0–3929862.5992360003. No scenario is preferred or species calibrated. Crucially, PSMC and pi use the same H1/H2 pair and are explicitly non-independent.

## Annotation and gene-conversion separation

P07 has an exact-native annotation with equal assembly/annotation sequence dictionaries and six measured partitions: CDS, fourfold, fourfold-W, fourfold-S, WS, and SW. The fourfold estimate is 0.000847594, with W and S estimates 0.001166236 and 0.000680551; WS and SW normalized counts are 0.000925017 and 0.000474837. These are real assembly-derived annotation results, but the parent pair is read-invalidated. Mutation spectrum, polarization, and reconstruction error are unresolved, so the asymmetry is suggestive only.

The direct-pedigree/gamete, population AFS, historical phylogenetic, and non-allelic paralog branches remain separate in `gene_conversion_branches.tsv`. No actual conforming VGP estimate exists in any branch. H1/H2 WS/SW counts are not transmission events, population B, historical substitutions, or copy-resolved paralog tracts.

## Implication for Lewontin's paradox

The measured assembly range demonstrates that the program now contains real nonzero biological estimates. Its stronger lesson is methodological: inferred diversity can change sharply with representation, common-mask definition, and read-backed haplotype validation. Any test of Lewontin's paradox must therefore separate biological between-species variation from assembly/callability variation and must sample populations rather than single assembly individuals.

This pilot cannot itself test the paradox. Two selected assembly observations—one invalidated, one validation-pending—cannot identify a cross-vertebrate diversity distribution, a relationship with census size, or a selection/drift explanation. Reporting 2.144254-fold as the vertebrate diversity range would ignore the exact-read adjudication and nonrandom execution attrition.

## Claim classification and reproducibility

The claim ledger contains 6 supported, 3 bounded, 3 suggestive, and 4 unidentifiable claims. Every claim states its sampling unit, covariance, and forbidden inference. `digest_ledger.tsv` rehashes every repository input and inventories every embedded upstream SHA-256 binding; `manifest.json` binds all emitted outputs.

Both compute environments use GNU Guix channel `44bbfc24e4bcc48d0e3343cd3d83452721af8c36`. The core pi/PSMC workflow profile is `/gnu/store/3c2mxm30rbzvnw7qsi235mrkk3m38fym-profile` (closure SHA-256 `8fcdb32021f1cd8eac839509cff47ab6bdd63b656b30e243fdf78d3c4ba24f9d`); raw-read validation uses `/gnu/store/n3vizxfw5ilggrinaj2mmbmng5ja4d6d-profile` (closure SHA-256 `ac0cb3601e56ef62b9ef99419de3659b2a2ba59b2aead29bc5f1928b50c83da2`). Both completed FastGA mappings passed live node-local `/scratch` contracts and reproduced the frozen PAFs. Canonical data, CAS objects, views, raw reads, run outputs, and promoted products resolve from the single configured VGP root shown above; the legacy Lewontin-paradox path remains migration input only.
