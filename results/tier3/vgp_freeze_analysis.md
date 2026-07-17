# VGP Phase 1 Freeze Analysis

## Data

The raw manifest is committed alongside this document:

- **TSV path:** `results/tier3/VGPPhase1-freeze-1.0.tsv`
- **Source URL:** <https://github.com/VGP/vgp-phase1/blob/main/VGPPhase1-freeze-1.0.tsv>
- **Format:** 717 lines (1 header + 716 data rows), tab-separated
- **Key columns for Tier 3A work:**
  - Col 10: Scientific Name (species)
  - Col 13: Assembly status ("4" = completed)
  - Col 16: Main haplotype GCA accession
  - Col 17: RefSeq annotation GCF accession (GFF available when present)
  - Col 21: Paired haplotype assembly ID (non-empty = reciprocal pair exists)
  - Col 26: Annotation status ("Completed NCBI", "Ready", "Wait", empty)

To fetch the latest version:
```bash
wget -q https://raw.githubusercontent.com/VGP/vgp-phase1/main/VGPPhase1-freeze-1.0.tsv \
     -O results/tier3/VGPPhase1-freeze-1.0.tsv
```

## Source
- File: `VGPPhase1-freeze-1.0.tsv`
- URL: https://github.com/VGP/vgp-phase1/blob/main/VGPPhase1-freeze-1.0.tsv
- Date downloaded: 2025-01-XX

## Species-Level Summary

| Metric | Count |
|--------|-------|
| Total unique species in manifest | 714 |
| Completed assemblies (status=4) | 223 |
| With annotation (Completed NCBI/Ready) | 248 |
| With paired haplotype | 271 |
| **Assembly + Annotation** | **120** |
| **Assembly + Annotation + Paired (Tier 3A-ready)** | **40** |

## Taxonomic Breakdown of Tier 3A Candidates (Assembly + Annotation + Paired)

| Group | Count |
|-------|-------|
| Fish (Actinopterygii/Sarcopterygii) | 13 |
| Amphibia | 4 |
| Reptilia | 3 |
| Mammalia | 9 |
| Aves | 6 |
| Other (singletons) | 5 |

## Fish Species with Full Pipeline Requirements (Assembly + RefSeq Annotation + Paired Haplotype)

| Species | Order | H1 GCA | H2 Haplotype |
|---------|-------|--------|-------------|
| Acipenser ruthenus | Acipenseriformes | GCA_902713425.2 | fAciRut3.pat |
| Amia calva | Amiiformes | GCA_036373705.1 | fAmiCal2.hap2 |
| Lepisosteus oculatus | Lepisosteiformes | GCA_040954835.1 | fLepOcu1.hap1 |
| Lampris incognitus | Lampridiformes | GCA_029633865.1 | fLamInc1.hap1 |
| Syngnathus acus | Syngnathiformes | GCA_901709675.2 | fSynAcu2.pri |
| Enoplosus armatus | Centrarchiformes | GCA_043641665.1 | fEnoArm2.hap2 |
| Pempheris klunzingeri | Pempheriformes | GCA_042242105.1 | fPemKlu1.hap2 |
| Cyclopterus lumpus | Scorpaeniformes | GCA_009769545.1 | fCycLum2.pri |
| Heterodontus francisci | Heterodontiformes | GCA_036365525.1 | sHetFra1.hap1 |
| Hemiscyllium ocellatum | Orectolobiformes | GCA_020745735.1 | sHemOce1.mat |
| Heptranchias perlo | Hexanchiformes | GCA_035084215.1 | sHepPer1.hap2 |
| Pristiophorus japonicus | Pristiophoriformes | GCA_044704955.1 | sPriJap1.hap2 |
| Narcine bancroftii | Torpediniformes | GCA_036971445.1 | sNarBan1.hap1 |

## Current Tier 3A Species vs VGP Freeze

| Species | In VGP? | Assembly? | Annotation? | Paired? |
|---------|---------|-----------|-------------|---------|
| Spinachia spinachia | Yes (fSpiSpi1) | Yes (GCA_048126635.1) | **No (VGP-deposited only)** | Yes (fSpiSpi1.hap2) |
| Menidia menidia | Yes (fMenAtl1) | Yes (GCA_048628825.1) | **No (VGP-deposited only)** | Yes (fMenAtl1.hap2) |
| Tautogolabrus adspersus | Yes (fTauAds1) | Yes (GCA_020745685.1) | **No** | No (only alt haplotype) |

## Key Finding: Current Tier 3A vs Annotated Tier 3A Candidates

The **current 3 Tier 3A species** were acquired from VGP assemblies with **VGP-deposited** GFFs (not RefSeq-annotated). The **13 fish candidates** listed above have **RefSeq (GCF) annotations** which represent the highest quality, version-coupled annotation standard.

For Tier 3A (diploid diversity), the pipeline needs:
1. Paired H1/H2 assemblies ✓ (all candidates have this)
2. Native annotation on H1 ✓ (all candidates have this via VGP/EGAPx)

The **RefSeq annotation** is more critical for Tier 3C (GC3 composition).

## Tier 3C Candidates (Assembly + RefSeq Annotation, any haplotype status)

46 fish species have completed assemblies with RefSeq annotations. Of these:
- 13 have paired haplotypes (strongest candidates)
- 46 total (can still run Tier 3C composition analysis)

## Data Source Clarification

- **Tier 3A** uses GCA assemblies downloaded from NCBI Datasets API for PRJNA489243
- Annotations come from VGP-deposited GFFs (NCBI EGAPx 0.4.1-alpha; Gnomon)
- RefSeq annotation status is tracked separately in the manifest column "Annotation status"
- VGP assemblies can have annotations before they enter the RefSeq pipeline
