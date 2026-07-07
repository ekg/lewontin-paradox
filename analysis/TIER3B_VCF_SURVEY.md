# Tier 3b VCF availability survey

Status: **survey** (data acquired, not yet downloaded/processed). Written 2026-07-07.

This surveys the availability of population-scale resequencing data
(VCFs or raw WGS reads) for the species where Tier 3b matters: the
**high-Nₑ invertebrates** that define the saturation plateau (the part of
the curve VGP cannot reach), plus the VGP-vertebrate anchors that double
as population-genomics model species.

## Method

- **WGS run/sample counts** from the ENA portal API
  (`https://www.ebi.ac.uk/ena/portal/api/count`, query
  `tax_tree(<taxid>) AND library_strategy=WGS`), which mirrors SRA.
  Counts are *experiments/runs*; for most population projects one run ≈
  one sequenced individual. Run by run, 2026-07-07.
- **Named population-VCF resources** confirmed via web search (resource
  portals + primary papers). These deliver pre-called VCFs — far more
  useful than raw reads — and are the preferred source where they exist.
- NCBI taxids resolved via E-utils esearch on the taxonomy DB.

## Results: WGS run counts (ENA)

| species | taxid | WGS runs | Nₑ (Buffalo log₁₀) | in Buffalo? |
|---|---:|---:|---:|:--|
| **Caenorhabditis elegans** | 6239 | 23,248 | 15.2 | ✓ |
| **Drosophila simulans** | 7240 | 2,766 | 14.7 | ✓ |
| **Drosophila melanogaster** | 7227 | 32,055 | 14.6 | ✓ |
| **Aedes aegypti** | 7159 | 3,127 | 14.4 | ✓ |
| Drosophila ananassae | 7217 | 82 | 14.2 | ✓ |
| Drosophila kikkawai | — | ? | 14.1 | ✓ |
| Nasonia vitripennis | 7425 | 160 | 13.9 | ✓ |
| Ceratitis capitata | 7213 | 200 | 13.9 | ✓ |
| Drosophila subobscura | 7241 | 50 | 13.8 | ✓ |
| **Anopheles gambiae** | 7165 | 71,603 | 13.8 | ✓ |
| Anopheles arabiensis | 7173 | 2,402 | 13.7 | ✓ |
| Drosophila pseudoobscura | 7237 | 1,393 | 13.5 | ✓ |
| **Acyrthosiphon pisum** | 7029 | 209 | 13.4 | ✓ |
| **Drosophila buzzatii** | — | ? | 13.3 | ✓ |
| **Culex pipiens** | 7175 | 378 | 13.0 | ✓ |
| Physa acuta | — | ? | 12.6 | ✓ (mollusk) |
| Anopheles merus | 30066 | 271 | 12.6 | ✓ |
| Drosophila persimilis | 7234 | 60 | 12.6 | ✓ |
| **Daphnia pulex** | 6669 | 4,768 | 12.4 | ✓ |
| Apis mellifera | 7460 | 9,366 | 12.4 | ✓ |
| Daphnia magna | 35525 | 1,387 | 12.3 | ✓ |
| Apis cerana | — | ? | 12.3 | ✓ |
| Pheidole pallidula | — | ? | 12.1 | ✓ |
| Thymelicus lineola | — | ? | 12.1 | ✓ |
| Armadillidium vulgare | — | 32 | 12.1 | ✓ |
| Ixodes ricinus | — | 548 | 12.0 | ✓ |
| Anopheles quadriannulatus | — | ? | 12.0 | ✓ |
| Bombyx mori | 7091 | 6,280 | — | (not in Buffalo core) |
| Tribolium castaneum | 7070 | 128 | — | (not in Buffalo core) |
| **Heliconius melpomene** | 34740 | 809 | — | (not in Buffalo core) |
| Crepidula fornicata | 176853 | 0 | — | ✓ (mollusk, no WGS) |
| Crassostrea gigas | 29159 | 1,768 | — | ✓ (oyster) |
| Ruditapes philippinarum | 129788 | 392 | — | ✓ (clam) |
| Biomphalaria glabrata | 6526 | 725 | — | ✓ (snail) |
| Strongylocentrotus purpuratus | 7668 | 489 | — | (sea urchin) |
| **Vertebrate anchors (VGP + pop data)** | | | | |
| Homo sapiens | 9606 | 1,178,623 | 8.4 | ✓ |
| Mus musculus | 10090 | 76,671 | — | ✓ |
| Bos taurus | 9913 | 25,118 | 6.9 | ✓ |
| Gallus gallus | 9031 | 31,750 | 8.2 | ✓ |
| Danio rerio | 7955 | 17,447 | 10.6 | ✓ |
| Ficedula albicollis | 59894 | 271 | 9.2 | ✓ |
| Gasterosteus aculeatus | 69293 | 10,098 | 10.6 | ✓ |
| Oryzias latipes | 8090 | 4,313 | — | (medaka) |

## Confirmed named VCF resources (pre-called, preferred over raw reads)

| species | resource | n (individuals) | VCF? | source |
|---|---|---:|:--:|---|
| **Anopheles gambiae/arabiensis/merus/coluzzii** | **Ag1000G phase 3 (Ag3.0)** | 2,784 wild + 297 cross-progeny | ✓ SNP VCFs | malariagen.net/data_package/ag1000g-phase3-snp; via `malariagen/data` Python package |
| **Drosophila melanogaster** | **DGRP Freeze 2.0** | 205 inbred lines | ✓ VCF | NCBI bioproject 36679; BCM-HGSC; 4.85M SNPs |
| **Drosophila simulans** | 170-line panel (Rogers et al. 2018) | 170 inbred lines | ✓ VCF (Zenodo) | doi:10.5281/zenodo.154261; PMC5767965 |
| **Drosophila pseudoobscura** | population genomics panel (2026 preprint) | many | ✓ (new, chromosome-scale) | biorxiv 10.64898/2026.02.02.703370 |
| **Caenorhabditis elegans** | **CaeNDR** | many wild isolates, multiple releases | ✓ hard-filtered VCF (latest 20250625) | caendr.org/data/data-release/c-elegans |
| **Daphnia pulex** | Lynch-lab population genomics | 83 → >800 isolates (temporal series) | ✓ (variant calls) | doi:10.1534/genetics.116.190611; PMC5419477, PMC9642108 |
| **Aedes aegypti** | 1206-genome global panel (Science 2025) | 1,206 genomes, 73 locations | ✓ (variant calls) | doi:10.1126/science.ads3732 |
| **Apis mellifera** | population genomics projects | many | (raw reads on ENA: 9,366 runs) | ENA + multiple bioprojects |

## Read of the survey

**The plateau is reachable on data.** Every high-Nₑ species that defines the
saturation zone (Nₑ > 10¹¹) has population-scale resequencing, and the
densest ones have **pre-called VCFs**:

- **Drosophila** — the densest high-Nₑ coverage in Buffalo (~15 species).
  DGRP (melanogaster, 205 lines), the 170-line simulans panel,
  pseudoobscura (2026 preprint), plus raw reads for ananassae, yakuba,
  erecta, mojavensis, persimilis, miranda, subobscura. This is the
  single best-sampled clade for the within-clade B∝Nₑ test.
- **Anopheles** — Ag1000G phase 3 (2,784 mosquitoes, 4 species:
  gambiae, arabiensis, merus, coluzzii) gives a second, independent
  dipteran clade at high Nₑ. Cross-clade Drosophila-vs-Anopheles tests
  the within-clade B∝Nₑ law where the recombination machinery is conserved.
- **C. elegans** — CaeNDR, but **selfing** → π is low and Nₑ estimation
  is subtle; use for **composition (GC3)**, not the π fit. Flag separately.
- **Aedes aegypti** — 1,206-genome global panel (2025) is large, but
  Aedes has huge structural variation and transposable-element activity;
  π handling needs care.
- **Daphnia pulex** — cyclical parthenogen, >800 isolates; π and Nₑ
  estimation are non-standard but the data is deep.
- **Nasonia vitripennis** — haplodiploid; π handling differs (males
  haploid). Flag.
- **Acyrthosiphon pisum** (pea aphid) — cyclical parthenogen;
  low-coverage WGS exists (209 runs).
- **Mollusks** (Crepidula, Physa, Bostrycapulus) — the high-Nₑ mollusks
  in Buffalo have **little to no WGS** (Crepidula fornicata: 0 runs).
  These drop to Tier 3c (assembly GC3 only).

**Vertebrate anchors** (human, mouse, cow, chicken, zebrafish, flycatcher,
stickleback) all have deep population data — but they sit on the **rising
limb / low plateau edge** (Nₑ 6–11), where VGP phased assemblies already
serve Tier 3a. Their population VCFs are a cross-check on the
individual-heterozygosity estimate, not the plateau test.

## Practical consequence for the plan

- **Tier 3b is feasible at meaningful n for the plateau**: Drosophila
  (~6–10 species with VCFs or enough raw reads to call) + Anopheles (4
  species, one VCF resource) + Aedes (1) + Daphnia (1) + a few others =
  **~12–18 high-Nₑ species with true population π**. This is the
  discriminating subset, and it is exactly where the model predicts
  saturation.
- **Within-clade Drosophila** is the strongest single test: ~6 species
  spanning Nₑ 10¹²·⁵–10¹⁴·⁷, all with population data, same recombination
  machinery → B∝Nₑ should hold → GC3 and π_S/π_W should saturate together.
- **Cross-clade Drosophila vs Anopheles** (both Diptera, conserved
  machinery, independent population resources) tests the
  within-recombination-class prediction the §4 mechanism makes.
- **Selfing/haplo-diploid/parthenogenetic species** (C. elegans, Nasonia,
  Daphnia, Acyrthosiphon, some aphids) need separate handling: keep for
  **composition (GC3)**, exclude from or flag in the **π fit**. This is a
  known, manageable confound, not a blocker.

## What is NOT available

- **High-Nₑ mollusks** (Crepidula plana, Physa acuta, Bostrycapulus
  aculeatus — all Nₑ > 10¹² in Buffalo): essentially no population WGS.
  These can only contribute GC3 (Tier 3c, from assemblies), not π. They
  widen the composition picture but cannot test the plateau in π.
- **Most non-model high-Nₑ invertebrates** (termites, isopods, some
  ants, butterflies other than Heliconius): sparse WGS, no pre-called VCFs.
  Tier 3c (assembly GC3) only.
- **Pre-called VCFs with 4D-site annotations** ready to use: none ship
  ready-made. We must compute π and the W/S stratification ourselves from
  the VCF + the species' reference GFF. Tooling gap: `pysam`, `pyfaidx`,
  `biopython`, `cyvcf2` (none installed; all pip-installable).

## Next concrete step

Install the tooling (`pip install pysam pyfaidx biopython cyvcf2`),
then fetch the two flagship VCFs that are easiest and span the plateau
extremes: **DGRP** (D. melanogaster, Nₑ 10¹⁴·⁶, the top) and **Ag1000G
phase 3** (Anopheles, Nₑ 10¹³·⁸, a second clade). Compute population π
+ W/S-stratified π at 4D for both, as the first real Tier-3b data points
and a pipeline validation. Then extend across the Drosophila clade.
