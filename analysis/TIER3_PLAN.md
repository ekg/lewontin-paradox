# Tier 3 data-driven plan: collecting GC, π, and W/S-stratified diversity to test gBGC saturation

Status: **plan** (not yet executed). Written 2026-07-07.

## 0. The question and the discriminating predictions

Tier 1 (Buffalo's total π, n=172) showed a rising-then-plateau but is
**under-identified**: the scatter is too large for curvature to beat a
straight power law (ΔAIC≈0), and pure saturation is rejected (ΔAIC=+179).
Total π *cannot* statistically establish gBGC saturation. Tier 3 is the
discriminating test — gBGC's two cross-species signatures that total-π
scatter cannot obscure, because they are different observables:

1. **Composition.** The equilibrium GC content at W/S-accessible sites,
   GC*, rises with Nₑ and **saturates** at high Nₑ. Background selection
   predicts no composition change. (Lean: `Composition.lean`.)
2. **W/S-stratified diversity.** At fourfold-degenerate (4D) codon sites,
   gBGC drives **S→ (G/C) alleles toward fixation**, so heterozygosity at
   S-ending sites should be **lower** than at W-ending (A/T) sites in
   high-Nₑ species, and **equal** in low-Nₑ species. The gap grows with
   Nₑ then **saturates**. BGS predicts no W/S asymmetry.
   (Lean: `MutationSelectionDrift.lean`, `GenomeMix.lean`.)

Both are cross-species (across Nₑ), both qualitative (shape, not just
slope), both impossible to explain with BGS alone. **That** is Tier 3.

## 1. The central design constraint: no single data source spans the Nₑ range

Buffalo's 173 core species span log₁₀ Nₑ ≈ 2.8 to 15.2. The **plateau** —
where saturation should bite — is at the top of that range (Nₑ > 10¹¹,
log₁₀ > 11). Those species are almost all **invertebrates**:

| clade (Buffalo core) | n | log₁₀ Nₑ range | where in the curve |
|---|---|---|---|
| Mammalia | 40 | 2.8 – 10.7 | rising limb → low plateau edge |
| Aves | 13 | 5.0 – 9.6 | rising limb |
| Teleostei | 10 | 7.0 – 10.6 | rising limb |
| Reptilia | 4 | 4.5 – 9.6 | rising limb |
| **Insecta (Drosophila, Anopheles, Aedes…)** | **58** | **~10 – 14.7** | **PLATEAU (saturation zone)** |
| Nematoda (C. elegans) | 1 | 15.2 | PLATEAU |
| Mollusca, Echinodermata, Cnidaria, Porifera, Annelida, Nemertea | ~30 | 7 – 13 | mixed; sparse assemblies |

**VGP is vertebrate-only.** The Vertebrate Genomes Project has produced
~hundreds of reference-grade phased diploid assemblies (NCBI bioproject
PRJNA538202 and sub-projects), covering mammals, birds, reptiles, amphibians,
fish — i.e. log₁₀ Nₑ ≈ 2.8–13.4, the **rising limb and the low edge of the
plateau**, but **not** the high-Nₑ invertebrates (Drosophila at 10¹⁴·⁵,
C. elegans at 10¹⁵·²) where saturation is sharpest.

**Consequence:** VGP alone cannot test the plateau. It tests the
**within-clade rising limb** densely (where the literature already has
support — Romiguier 2010 mammals, Weber 2014 birds, Barton & Zeng 2021
passerines all live here). To see **saturation**, we must add the
high-Nₑ invertebrates from **population resequencing** (DGRP, C. elegans
wild isolates, Anopheles population projects) — the one source with true
π *and* the Nₑ range to reach the plateau.

This is the honest core of the plan: **a hybrid of three sources**, each
chosen for the part of the curve and the observable it can deliver.

## 2. wdyt: VGP diploid → π? Yes, with two honest caveats

The instinct is right and elegant. A VGP assembly is a **phased diploid**:
two haplotypes (H1, H2) per chromosome, assembled from a trio (offspring +
two unrelated parents) or from long-read + HiC phasing of a single
individual. Comparing H1 vs H2 gives the **individual heterozygosity**,
which for unrelated-parent trios ≈ population heterozygosity ≈ π. So a
single VGP assembly delivers **all three Tier-3 observables from the same
genome**:

- **GC3** from the CDS (composition)
- **individual π** from H1-vs-H2 variants (the diversity point)
- **W/S-stratified π at 4D sites** from H1-vs-H2 variants partitioned by
  the current codon's W/S class (the clean within-species signature)

This is the strongest single-assembly design available. Two caveats:

1. **VGP doesn't reach the plateau** (above). It tests the rising limb +
   low-plateau edge within vertebrate clades. The saturation test needs
   the invertebrate population data (Tier 3b). Don't claim VGP alone
   settles saturation.
2. **Individual heterozygosity ≠ population π.** For unrelated-parent
   trios it's a good estimate (≈ population het), but trio assemblies from
   related parents are biased low (IBD), and a single individual gives no
   population variance. For the **shape** of the Nₑ–π curve (does it
   plateau?) this is fine — individual het tracks population π. For the
   **magnitude** of the plateau and for SFS-based gBGC estimates (B = 4Nₑb),
   a single individual is under-powered; those need population samples
   (Tier 3b). We'll report individual het as individual het and not
   conflate it with population π.

Net: **use VGP for the vertebrate within-clade composition + individual π +
within-individual W/S stratification (the rising limb, densely), and add
population VCFs for the high-Nₑ invertebrates (the plateau).** This is the
design.

## 3. The three collection tiers

### Tier 3a — VGP phased vertebrate assemblies (composition + individual π + W/S-strat)

- **Source.** NCBI assembly DB, bioproject PRJNA538202 (+ VGP sub-projects;
  enumerate via E-utils `esearch` on the assembly DB filtered by the VGP
  bioproject link, then `esummary` for accession + FTP). Cross-match to
  Buffalo's 71 vertebrate core species by scientific name (congener
  fallback where the exact species isn't in VGP — GC3 is conserved within
  genera well enough for a within-clade fit).
- **What we download.** Per species: the **primary haplotype FASTA** + the
  **alternate haplotype FASTA** (or the assembly variants file, `.vcf`/`.bed`
  of H1-vs-H2 variants, which VGP deposits) + the **GFF/GTF annotation**
  (needed to locate CDS / 4D sites). VGP deposits these on the NCBI FTP.
- **What we compute.**
  1. **GC3**: extract CDS, take the third codon position, GC fraction. (Whole-genome GC% as a cross-check.)
  2. **individual π**: density of H1-vs-H2 variants genome-wide (or from
     the deposited variant file). Convert to per-site heterozygosity.
  3. **W/S-stratified π at 4D sites**: restrict H1-vs-H2 variants to
     fourfold-degenerate codon sites (from the GFF + a genetic code table),
     partition by the site's W/S class (S = G|C, W = A|T), compute
     het-density per class. → the within-species signature.
- **The test.** Within each vertebrate clade (mammals, birds, teleosts),
  regress GC3 and the S/W het-gap on Nₑ (Buffalo's `pred_log10_N`).
  Prediction: GC3 rises then begins to plateau; the S/W het-gap grows
  then saturates. **Across** vertebrate clades the signal should attenuate
  (Galtier 2018) — the within-vs-across distinction the §4 mechanism predicts.
- **Expected n.** ~40–60 species with full VGP + annotation (after congener
  fallback and annotation-availability filtering). Enough for within-clade
  fits; not the plateau.
- **Confounds.** Annotation quality (4D-site calls need a clean GFF —
  VGP annotations lag assemblies, so some species will be GC3-only).
  Individual-vs-population π (caveat 2). 4D-site saturation at long
  within-clade divergence (use close outgroups / restrict to the rising
  limb where VGP sits). Controlled by clade and by annotation status.

### Tier 3b — Population resequencing for the high-Nₑ invertebrates (the plateau)

- **Source.** Public population VCFs + reference assemblies for the
  plateau species. The ones that matter (Buffalo's high-Nₑ core):
  - **Drosophila** (Buffalo has ~15 Drosophila species, the densest
    high-Nₑ coverage): DGRP (D. melanogaster, ~205 inbred lines), Drosophila
    Population Genomics Project (D. pseudoobscura, D. simulans, D.
    miranda, …). Reference assemblies on NCBI/FlyBase.
  - **Anopheles gambiae / arabiensis / merus** (Ag1000G, ~1000+ genomes).
  - **Aedes aegypti** (AaegL5 + population VCFs).
  - **Culex pipiens** (population data exists).
  - **Nasonia vitripennis** (haplodiploid — π handling differs; flag it).
  - **Caenorhabditis elegans** (wild-isolate VCFs, ~ hundreds of isolates;
    selfing — π is low and Nₑ estimation is subtle; the C. elegans point is
    composition, not π).
  - **Apis mellifera / Acyrthosiphon pisum / Ceratitis capitata** (population
    data exists for some).
- **What we download.** Per species: the reference assembly FASTA + GFF, and
  the population VCF (or a downsampled subset — n≈20 lines is plenty for a
  π estimate and a W/S-stratified SFS).
- **What we compute.**
  1. **GC3** (composition) from the reference CDS.
  2. **population π** (true π, not individual het) from the VCF.
  3. **W/S-stratified π at 4D sites** from the VCF, partitioned by W/S class.
  4. (If the VCF supports it) **B = 4Nₑb** via the SFS skew at W/S sites —
     the strength signature (Galtier 2018's method). This is the gold
     standard; only feasible here, not in VGP.
- **The test.** This is where the **plateau** lives. Regress GC3, π_S/π_W,
  and B on Nₑ across these high-Nₑ species. Prediction: GC3 plateaued
  (flat above Nₑ~10¹²), π_S ≪ π_W at 4D, B saturating. **Cross-clade**
  among insects (Drosophila vs Anopheles vs Aedes) tests the within-clade
  B∝Nₑ law where the machinery is conserved.
- **Expected n.** ~10–20 species with population VCFs after availability
  filtering. Small but it spans Nₑ ~10¹²–10¹⁵ — the saturation zone — which
  VGP cannot reach. This is the discriminating subset.
- **Confounds.** Sampling depth differences across VCFs (standardize π to
  a fixed sample size, n=20, by downsampling). Reference bias (use the
  focal-species reference, not a congener, for π). Selfing/haplo-diploidy
  (C. elegans, Nasonia — handle separately or exclude from the π fit, keep
  for composition). 4D-site polarization needs an outgroup for the SFS-B
  estimate (use a sister species); for the simple π_S vs π_W within-species
  comparison, no outgroup is needed (classify by current codon).

### Tier 3c — NCBI assemblies, GC3-only (the cross-species composition fill)

- **Source.** NCBI assembly DB via E-utils (already verified working: `esearch`
  + `esummary` returns accession + FTP path; tested on D. melanogaster).
  For every Buffalo core species (or congener) without a VGP assembly and
  without a population VCF — the non-model invertebrates (mollusks,
  echinoderms, cnidarians, sponges, annelids, nemerteans).
- **What we download.** The best-available RefSeq/GenBank assembly FASTA +
  GFF (GC3 needs annotation; whole-genome GC% needs only the FASTA). Single
  haplotype is fine — **composition doesn't need phasing**.
- **What we compute.** GC3 (where annotated) + whole-genome GC%. **No π**
  (single assembly) and **no W/S stratification** (no individual variation).
  This tier delivers composition only.
- **The test.** Adds the non-model invertebrates to the cross-species GC3-vs-Nₑ
  fit, widening the Nₑ range and the phylogenetic spread. Tests the
  **across-clade attenuation** prediction (Galtier 2018): GC3 should track
  Nₑ within clades but the relationship should be noisier/attenuated across
  these deeply divergent invertebrates — the §4 mechanism's signature.
- **Expected n.** ~40–80 additional species (after assembly + annotation
  availability; many non-model invertebrates have assemblies but poor CDS
  annotation → GC%-only). Whole-genome GC% covers more.
- **Confounds.** Annotation absence (whole-genome GC% is a weaker proxy for
  gBGC than GC3 because noncoding GC is dominated by other forces — report
  both, weight GC3). Assembly quality (contigs, contamination). Congener
  substitution (GC3 conserved enough within genera; flag the substitution).

## 4. The W/S-stratification method (the clever part)

For Tiers 3a and 3b, the within-species signature is computed as:

1. From the annotation (GFF/GTF), locate fourfold-degenerate codon sites in
   CDS (4D = codons where all changes at position 3 are synonymous — the
   standard vertebrate/insect genetic code table).
2. For each 4D site, classify the **current (reference) third base** as
   **S** (G or C) or **W** (A or T).
3. Compute heterozygosity (Tier 3a: H1-vs-H2 variant density; Tier 3b:
   population π from the VCF) **separately** for S-class and W-class 4D sites.
4. The signature is **π_W − π_S** (or the ratio π_S/π_W), as a function of
   Nₑ. Prediction: ≈0 at low Nₑ; positive and growing at mid Nₑ; **saturating**
   (not growing further) at high Nₑ. BGS predicts ≈0 everywhere.

This needs **no outgroup** (classification is by the current codon, not the
ancestral state). An outgroup is only needed for the SFS-based B estimate
(Tier 3b, optional gold standard). This makes Tier 3a (VGP, no outgroups)
fully capable of the within-species W/S test — the user's instinct that a
single diploid assembly suffices is correct.

## 5. Tooling status and gaps

Verified present: `curl` + NCBI E-utils (assembly discovery + FTP paths
work), `samtools`, `bedtools`, python3 (numpy, scipy, pandas, matplotlib).

**Missing (need to `pip install` before execution):**
- `pysam` (VCF/FASTA parsing) — Tier 3b VCF handling
- `pyfaidx` or `biopython` (FASTA indexing, CDS extraction) — all tiers
- `cyvcf2` (faster VCF; optional) — Tier 3b
- `bcftools` / `vcftools` (CLI; optional, pysam covers it) — Tier 3b
- `samtools` is present but `bcftools` is not (variant calling from
  H1-vs-H2 alignments if the VGP deposited variant file is unavailable) —
  Tier 3a fallback; install `bcftools` or use pysam to call H1-vs-H2 diffs
  directly from the two haplotype FASTAs.

All pip-installable; no admin needed. First execution step is to install
these and re-verify.

## 6. Deliverables and order of operations

1. **Tooling.** `pip install pysam pyfaidx biopython cyvcf2`; install or
   alias `bcftools`. Verify on a test VCF/FASTA.
2. **VGP enumeration (Tier 3a).** Script `analysis/tier3a_vgp_collect.py`:
   E-utils query of VGP bioproject → assembly accessions + FTP →
   cross-match to Buffalo's 71 vertebrate core (name + congener) →
   download manifest (species, accession, hap1 FASTA, hap2 FASTA/variants,
   GFF). Log which species have annotations (gate the 4D analysis).
3. **Tier 3a compute.** Script `analysis/tier3a_vgp_compute.py`:
   per species → GC3, whole-genome GC%, individual π (H1-vs-H2),
   W/S-stratified π at 4D. Output `analysis/tier3a_data.tsv`.
4. **Tier 3b enumeration + compute.** Script `analysis/tier3b_popvcf_collect.py`:
   hard-code the known population-VCF projects (DGRP, Ag1000G, Aedes, C. elegans
   wild isolates, …), fetch a downsampled VCF (n=20) + reference + GFF per
   species. Script `analysis/tier3b_popvcf_compute.py`: population π,
   W/S-stratified π at 4D, (SFS-B where outgroup available). Output
   `analysis/tier3b_data.tsv`.
5. **Tier 3c fill.** Script `analysis/tier3c_ncbi_gc.py`: E-utils query for
   the remaining core species (or congeners) → best assembly → GC3 + GC%.
   Output `analysis/tier3c_data.tsv`.
6. **Merge + fit.** Script `analysis/tier3_fit.py`: join all three tiers to
   Buffalo's `pred_log10_N` (+ clade, + class). Fits: GC3-vs-Nₑ
   (within-clade and across-clade), π_S/π_W-vs-Nₑ (the saturation shape),
   B-vs-Nₑ (Tier 3b). Figure `analysis/fig_tier3.{pdf,png}`: panels for
   composition, W/S-stratified diversity, within-vs-across-clade.
7. **Manuscript.** Replace the literature-citation Tier-3 section with the
   data-driven results, keeping the within-vs-across-clade framing and the
   §4-mechanism reconciliation. Honest reporting of n, confounds, and what
   each tier can/cannot establish.

## 7. What we will be able to claim (and not)

- **Tier 3a (VGP) can claim:** within vertebrate clades, GC3 rises with Nₑ
  and the within-individual W/S het-gap grows — the composition + rising-limb
  signature, on our own data, for the clades the literature already studied.
  It **cannot** claim the plateau (doesn't reach the Nₑ range) and cannot
  estimate B (single individual).
- **Tier 3b (pop VCFs) can claim:** at high Nₑ, π_S/π_W saturates and (if
  SFS-B computable) B saturates — the **plateau**, on the species where it
  should bite. Small n (~10–20) but it is the discriminating subset. It
  **cannot** claim dense cross-clade coverage.
- **Tier 3c (NCBI GC3) can claim:** the cross-species composition
  relationship across deeply divergent invertebrates, testing the
  across-clade attenuation. It **cannot** say anything about π or W/S.
- **Together** the three tiers test the full prediction: gBGC scales with
  Nₑ **within** clades (3a + 3b), saturates at high Nₑ **within** the
  saturation zone (3b), and attenuates **across** deep clades (3c) — exactly
  the within-vs-across-clade pattern the §4 homology-maintenance mechanism
  predicts, now on data we collect rather than cite.

**Not the whole story (held):** the low-Nₑ rising limb is Nₑ/Nₑ reduction
(Buffalo's territory); gBGC saturation is the high-Nₑ plateau and the
within-clade composition/strength signatures. The plan is designed to show
that distinction empirically, not to claim gBGC explains all of Lewontin's
Paradox.

## 8. Honest risks

- **Annotation lag.** VGP and many NCBI assemblies ship FASTAs before clean
  GFFs. The 4D analysis (the W/S test) is gated on annotation; GC3 is gated
  on CDS annotation; whole-genome GC% is not. Expect the W/S test on fewer
  species than the GC test. Report n per observable.
- **Population-VCF sampling heterogeneity.** π is sensitive to sample size
  and structure; standardize by downsampling to a fixed n and report
  per-species n. Selfing/haplo-diploid species (C. elegans, Nasonia, some
  aphids) need separate handling or exclusion from the π fit.
- **VGP phasing completeness.** Some VGP assemblies are primary-only or
  have low alternate-haplotype coverage. The individual-π and W/S tests
   need a real H2 (or a deposited variant file). Filter on phased-ness.
- **Congener substitution.** GC3 is conserved within genera, but π is not —
  congener fallback is valid for composition (Tier 3c) and questionable for
  π (Tier 3a/b use the focal species or skip π). Log every substitution.
- **4D-site definition.** The standard code's 4D sites are well-defined for
  the species here (vertebrates, insects, nematodes); mitochondrial and
  non-standard codes are excluded. C. elegans uses the standard nuclear
  code; fine.
