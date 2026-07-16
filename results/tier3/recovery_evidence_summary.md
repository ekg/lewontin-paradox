# Recovered Tier 3 biological evidence

Status: **PASS — corrected Tier 3A, recovered Tier 3B, and Tier 3C composition are synthesized without modality pooling.**

## Dependency and input gate

The WG task has direct dependencies on `rerun-tier3a-origin-sweepga` and `run-tier3b-biological-recovery`; both were confirmed `Done` before this synthesis began. The only numerical Tier 3A input is `results/tier3a/diploid_diversity.tsv`; the old-result and mapping supersession ledgers are negative controls and no listed path or checksum was selected.

Corrected Tier 3A traces all three tuples to fetched SweepGA origin/main commit `018e4ce49d2c125820e0ac50dc5feaa02d423683`, binary SHA-256 `fa7f0edb9b7e275c288db254046020e136d4267dd5ee043379227ef80da0573b`, Guix closure `/gnu/store/z9v2f6faha9cwjz0sm5iphhlzisgi077-profile`, and native `--num-mappings 1:1` commands with observed query/target multiplicity 1/1. Tier 3B uses Guix closure `/gnu/store/z9v2f6faha9cwjz0sm5iphhlzisgi077-profile` at channel commit `44bbfc24e4bcc48d0e3343cd3d83452721af8c36`. Both upstream QC reports pass.

## Assembly modality (Tier 3A)

These are alignment-conditioned H1/H2 estimates over deterministic H1-native coding panels, not population diversity and not genome-wide deposited-individual heterozygosity. Each estimate is variants/callable bases; uncertainty is the 1,000-replicate 50-kb genomic block bootstrap. The ratio intervals below conservatively divide the marginal S and W 95% bounds because paired bootstrap draws were not published; they are intentionally wider than a paired ratio interval.

| assembly pair | coding-gene diversity (95% CI; variants/callable) | CDS diversity (95% CI; variants/callable) | reference-conditioned pi_S/pi_W (conservative interval) |
|---|---|---|---|
| _Menidia menidia_ | 0.0138517 [0.0123545, 0.0153404]; 28233/2038234 | 0.00674737 [0.00589925, 0.00776835]; 2044/302933 | 0.672544 [0.482732, 0.932999] |
| _Spinachia spinachia_ | 0.000397957 [0.000329512, 0.000481235]; 803/2017806 | 0.000232752 [0.0001483, 0.000338608]; 63/270674 | 0.316941 [0.0685349, 1.70202] |
| _Tautogolabrus adspersus_ | 0.00167799 [0.00123835, 0.00212181]; 3551/2116227 | 0.000515957 [0.000344761, 0.000713813]; 115/222887 | 0.760295 [0.27804, 2.12518] |

Coding diversity is heterogeneous across the three pairs (about 0.000398 to 0.01385). All three reference-conditioned ratio point estimates are below one, but the conservative intervals for _Spinachia_ and _Tautogolabrus_ include one; these three selected assembly pairs do not estimate a population-size trend.

## Population modality (Tier 3B)

Both rows use 20 wild diploid biological individuals (40 nominal chromosomes) from the same species and region, with exact population-specific callable masks. Population pi and the ratio use 10,000-replicate chromosome-stratified 1-Mb block bootstrap intervals; component pi_S and pi_W rows use 20-unit delete-one-individual jackknife intervals.

| population | population pi (95% CI; pairwise numerator/callable) | reference-conditioned pi_S/pi_W (95% block-bootstrap CI) |
|---|---|---|
| `AO_Luanda_2009_coluzzii` | 0.01235003 [0.01186339, 0.01281725]; 184914.842/14972821 | 0.7426563 [0.7129171, 0.7708752] |
| `GM_WaliKunda_2012_coluzzii` | 0.01586979 [0.01539494, 0.01633848]; 237893.579/14990338 | 0.723555 [0.6919446, 0.7571679] |

The two populations differ in pi and their bootstrap intervals do not overlap; both reference-conditioned ratio intervals lie below one. This is within-species population heterogeneity, not an across-species effect and not polarized SFS-B.

## Composition modality (Tier 3C)

Composition remains an exact-single-assembly observable: 135 whole-genome GC values and 90 native-annotation GC3 values. The class-fixed GC3 slope remains +0.00892 per Buffalo census-size-proxy unit (10,000-species-bootstrap 95% interval -0.00123 to +0.02149; n=90; BH q=0.217), with positive Aves, uncertain Insecta, negative Mammalia, and null _Drosophila_ estimates. The positive quadratic coefficient remains opposite the predicted concavity. These heterogeneous/null composition findings are not replaced or pooled with diversity.

## Claim boundary and sensitivities

The recovery replaces the old n=0 statements with real biological estimates, but it does not create a cross-species diversity regression: Tier 3A has three selected coding-panel assembly pairs, Tier 3B has two populations from one species, and Tier 3C measures composition. Reference-conditioned pi_S/pi_W is descriptive and is never renamed polarized SFS-B. Callable-mask, sample, panel, native-annotation, and mapping exclusions remain explicit in `recovery_evidence_ledger.tsv`; no causal claim that gBGC resolves Lewontin's paradox follows.

`headline_audit.tsv` independently recomputes all 20 upstream direct estimates plus three composition coefficients; `artifact_index.tsv` records checksums and roles for the committed synthesis. Figures are regenerated from `analysis/tier3_results.tsv` and the manuscript reports the same definitions, denominators, n, and intervals.
