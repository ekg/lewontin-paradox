# Tier 3 execution plan: moving the empirical work to a cluster

Status: **plan** (not yet executing). Written 2026-07-07.

This document *operationalizes* `analysis/TIER3_PLAN.md` (the science) and
`analysis/TIER3B_VCF_SURVEY.md` (the data availability) onto a cluster.
It does **not** pick a specific cluster — see §2 (cluster selection).

## 1. Why move off the laptop

Concrete limits we are about to hit:

- **Uptime.** VGP assembly downloads and Ag1000G Zarr pulls take hours;
  the laptop sleeps, networks drop, jobs don't resume. A cluster job
  queues and runs unattended.
- **Disk.** Phased vertebrate assemblies (Tier 3a) are ~1–3 GB/haplotype
  ×2 ×~60 species ≈ a few hundred GB; population VCFs add tens of GB;
  NCBI assemblies (Tier 3c) tens more. Laptop scratch is fragile for this.
- **Parallelism.** Tier 3c is ~40–80 independent per-species GC3 jobs; Tier
  3a is ~60 independent per-species compute jobs. These are embarrassingly
  parallel — perfect for an array job, wasteful to serialize on a laptop.
- **Memory.** Whole-genome VCF iteration (Ag1000G, DGRP) and 4D-site
  stratification are memory-happier with ≥16–32 GB; some laptops can't.

The analysis *outputs* (TSVs, figures) are small — KB to MB — so the
laptop remains the place to edit the manuscript and read results. The
cluster is where the data lives and the compute runs.

## 2. Cluster selection (open decision — needs your input)

I don't know your compute environment, so this is parameterized. Pick one:

| option | scheduler | notes | fits? |
|---|---|---|---|
| **Institutional HPC** | SLURM (most common) / PBS-Torque / LSF | Usually free for affiliates; scratch + module system + conda/Singularity. **Best default if you have access.** | ✓ if you have an account |
| **ACCESS (US national)** | SLURM | Bridges-2 (PSC), Expanse (SDSC), Delta (NCSA) — big-mem nodes, good for VCF iteration. Needs an allocation (startup/beginner allocation is easy to get). | ✓ if US-affiliated |
| **Commercial cloud (GCP/AWS spot)** | none (or Slurm-on-cloud) | Strongest for Ag1000G: the data is already on GCS (same region = free egress). Spot VMs are cheap for batch. More setup. | ✓ if you want GCS-colocated |
| **A lab server / collaborator's machine** | varies | Often the path of least resistance if one exists. | ✓ if available |

**What the cluster must provide (real requirements):**

- A scheduler (SLURM assumed in the scripts below; trivial to adapt).
- **Internet egress on compute nodes** — non-negotiable: we stream from
  NCBI, ENA, VGP FTPs, and Google Cloud Storage (Ag1000G Zarr). Some HPCs
  only allow internet from login nodes; that's fine for *staging*, but the
  Zarr streaming path needs it on compute (or pre-stage to scratch).
- **Scratch storage**, ≥200 GB transient, ideally ≥500 GB for headroom.
- **≥16 GB RAM/node** for VCF iteration; 32–64 GB nice for the big VCFs.
- A way to install Python packages (`conda`/`mamba`, or `pip --user`, or
  Singularity/Apptainer if conda is blocked).
- `git` (to keep the repo on the cluster and push results back).

**Open question for you:** which cluster, and is it SLURM? I'll finalize
the job scripts once you confirm. Everything below is SLURM-shaped but
the per-species logic is scheduler-agnostic Python.

## 3. Software environment (one-time setup)

The survey identified the gaps (laptop has `samtools`, `bedtools`, numpy/
scipy/pandas/matplotlib; missing `pysam`, `pyfaidx`, `biopython`,
`cyvcf2`, `bcftools`). Fix once on the cluster:

```yaml
# analysis/envs/tier3.yaml
name: tier3
channels: [conda-forge, bioconda]
dependencies:
  - python=3.12
  - numpy, scipy, pandas, matplotlib
  - pysam          # VCF/FASTA parsing (Tier 3b)
  - pyfaidx        # FASTA indexing, CDS extraction (all tiers)
  - biopython      # GenBank/EML parsing for annotation
  - cyvcf2         # fast VCF iteration (Tier 3b)
  - htslib         # bgzip/tabix, bgzipped VCF
  - bcftools       # VCF subsetting/downsampling (Tier 3b)
  - samtools       # already present, but pin in the env
  - bedtools       # 4D-site BED intersection
  - wget, curl     # NCBI/ENA/VGP downloads
  - malariagen-data # Ag1000G Zarr-from-GCS (Tier 3b)
```

Create with `mamba env create -f analysis/envs/tier3.yaml` (or `conda`).
If conda is blocked on the cluster, build the same env as a
**Singularity/Apptainer** image from a conda-pack or a Dockerfile — the
job scripts then `apptainer run ...`.

## 4. Repo layout on the cluster (and what stays off the laptop)

Clone the repo **on the cluster** (login node or scratch checkout):

```
<scratch>/lewontin-paradox/        # git repo (scripts, manuscript, lean)
<scratch>/tier3data/              # staged assemblies + VCFs (gitignored)
  ├── assemblies/<species>/       # Tier 3a/3c FASTAs + GFF
  ├── vcfs/<species>/             # Tier 3b VCFs (or symlinks to Zarr cache)
  └── manifest.tsv                # accession + path + checksum per item
<scratch>/tier3out/               # per-species results (gitignored)
  └── <tier>/<species>.tsv
```

Principles:

- **Raw data is never committed.** It lives on scratch with a `manifest.tsv`
  recording accession + source URL + MD5 (reproducible; re-fetchable).
  Add `tier3data/` and `tier3out/` to `.gitignore`.
- **Small results are committed.** Each per-species job writes a tiny TSV
  (species, Nc, GC3, GC%, pi, pi_S, pi_W, n_samples, notes). The merge
  step concatenates these into `analysis/tier3_results.tsv` → committed.
  Figures (`analysis/fig_tier3.{pdf,png}`) → committed.
- The laptop keeps editing `manuscript.typ`, `lean/`, this plan. The
  cluster pushes results; you pull.

## 5. Storage & compute budget (estimates)

| item | size | notes |
|---|---|---|
| Tier 3a VGP assemblies | ~200–350 GB | ~60 species × 2 haplotypes × ~1–3 GB + GFFs. Stream-or-store. |
| Tier 3c NCBI assemblies | ~30–60 GB | ~60 species, single haplotype, mostly small. |
| Tier 3b DGRP VCF | ~5–10 GB | 4.85M SNPs × 205 lines. Small. |
| Tier 3b Ag1000G | ~3 GB/arm if staged; **~0 if streamed** | Zarr on GCS — `malariagen-data` reads regions/samples on demand. Downsample n=20. **Prefer streaming.** |
| Tier 3b other VCFs (simulans, Daphnia, Aedes) | ~5–40 GB each | Aedes 1206-genome is the big one; downsample. |
| Per-species result TSVs | KB | negligible |
| Total scratch | **~300–500 GB** headroom | one moderate scratch allocation |

Compute: each per-species GC3/π/4D job is minutes (3c) to ~1 h (3b big VCF).
~150 species-jobs total → a few hours of wall-time on an array, dominated by
download bandwidth, not CPU.

## 6. Execution sequence (phases)

The order is chosen to **validate the pipeline cheaply before the big
downloads**, and to hit the plateau first (where the signal is).

### Phase 0 — environment + manifest (1 job, login node)
- `mamba env create -f analysis/envs/tier3.yaml`.
- Generate `tier3data/manifest.tsv` by joining Buffalo's 173 core species
  to: VGP bioproject accessions (Tier 3a), population-VCF accessions
  (Tier 3b, from the survey), and NCBI assembly accessions (Tier 3c). This
  script (`analysis/tier3_manifest.py`) is the single source of truth for
  what gets fetched; it already exists conceptually as the VCF survey +
  the NCBI E-utils calls we verified.

### Phase 1 — Tier 3c NCBI GC3 (array job, warm-up + validation)
- *Why first:* highest n (~60–80 species), cheapest (1 assembly each, no
  phasing, no VCF), pure composition. Validates the GC3 pipeline on the
  easy case before the phased/VCF variants.
- *Job:* SLURM array over species; each task: `esearch`/`esummary` → best
  assembly → `wget` FASTA + GFF → compute GC3 + whole-genome GC% → write
  `tier3out/3c/<species>.tsv`. Idempotent (skip if output exists).
- *Test:* spot-check GC3 on *D. melanogaster* and *H. sapiens* against
  published values (dm6 GC3 ≈ 0.55; hg38 GC3 ≈ 0.52). Confirms the CDS
  extraction + third-position code before trusting it across 60 species.

  **Post-run audit (2026-07-14):** the two approximate values above had no
  recoverable citation or pooled-versus-gene-weighted definition and generated
  invalid upper bounds.  They remain historical failed anchors and MUST NOT be
  silently widened.  `TIER3C_CONTROL_AUDIT.md` replaces them for promotion
  with exact reproduction by an independent implementation, unchanged native
  provenance/CDS hard gates, and definition-aware published context.

### Phase 2 — Tier 3b population VCFs (the plateau)
- **2a — DGRP + Ag1000G first** (pipeline validation on the plateau
  extremes): *D. melanogaster* (Nₑ 10¹⁴·⁶, DGRP, the top) and *Anopheles
  gambiae* (Nₑ 10¹³·⁸, Ag1000G, second clade). Compute population π +
  W/S-stratified π at 4D for both. This is the **first real Tier-3b data
  point** and a pipeline check on the two biggest resources.
- **2b — the Drosophila clade** (the strongest single test, per the
  survey): extend to *D. simulans* (170-line VCF), *D. pseudoobscura*,
  and the raw-read species (*ananassae, yakuba, erecta, mojavensis,
  persimilis, miranda, subobscura* — call variants with `bcftools` from
  ENA WGS reads against the species reference). ~6–10 species spanning
  Nₑ 10¹²·⁵–10¹⁴·⁷, all Diptera → within-clade B∝Nₑ test.
- **2c — cross-clade check**: *Anopheles* (Ag1000G, 4 species) as the
  second dipteran clade. Drosophila-vs-Anopheles tests the
  within-recombination-class prediction.
- **2d — other high-Nₑ**: Aedes (1206-genome, downsample n=20), Daphnia
  pulex (cyclical parthenogen — flag), C. elegans (selfing → **composition
  only**, no π fit), Nasonia (haplodiploid → flag).
- *Selfing/haplo-diploid rule* (from the survey, carried through): these
  species contribute **GC3 to the composition fit**, are **excluded from /
  flagged in the π fit**. Documented in each result TSV's `notes` column.

### Phase 3 — Tier 3a VGP phased vertebrate assemblies (the rising limb)
- Array over the ~71 vertebrate core species (congener fallback where the
  exact species isn't in VGP — log every substitution). Each task: fetch H1
  + H2 (or deposited variant file) + GFF → GC3 + individual π (H1-vs-H2
  het) + W/S-stratified π at 4D (classify by current codon — **no outgroup
  needed**, per the plan). Gated on annotation availability: 4D analysis
  only where the GFF is clean; GC3 where CDS annotated; whole-genome GC%
  always.
- *Note:* VGP annotations lag assemblies, so expect 4D on fewer species
  than GC3. Report n per observable, honestly.

### Phase 4 — merge + fit + figure + manuscript (1 job, login or interactive)
- `analysis/tier3_fit.py`: join all three tiers' TSVs to Buffalo's
  `pred_log10_N` + clade/class. Fits: GC3-vs-Nₑ (within-clade and
  across-clade), π_S/π_W-vs-Nₑ (saturation shape), B-vs-Nₑ (Tier 3b where
  SFS computable). `analysis/fig_tier3.{pdf,png}`: composition / W/S
  stratified / within-vs-across-clade panels.
- Commit results + figure. Pull on laptop. Rewrite the manuscript Tier-3
  section from literature-citation to data-driven, keeping the
  within-vs-across-clade framing and the §4-mechanism reconciliation.

## 7. SLURM job template (concrete, short)

A real array-job skeleton for Phase 1 (Tier 3c). The other phases share
this shape with a different per-task script.

```bash
#!/bin/bash
#SBATCH --job-name=tier3c_gc3
#SBATCH --array=1-80          # one per species (manifest rows)
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --partition=YOUR_PARTITION   # fill in
#SBATCH --output=slurm/tier3c_%A_%a.out

set -euo pipefail
source activate tier3          # or: apptainer run tier3.sif
cd $SLURM_SUBMIT_DIR

SPECIES=$(sed -n "${SLURM_ARRAY_TASK_ID}p" tier3data/manifest_3c.txt)
OUT="tier3out/3c/${SPECIES// /_}.tsv"

# idempotent: skip if done
[ -s "$OUT" ] && { echo "skip $SPECIES (exists)"; exit 0; }

python analysis/tier3c_ncbi_gc.py "$SPECIES" "$OUT"
```

Run: `sbatch analysis/slurm/tier3c.sh`. Same pattern for `tier3a.sh`,
`tier3b.sh` with different `--mem`/`--time` (3b big VCFs: `--mem=32G
--time=04:00:00`). Submit dependent merge with `--dependency=afterok`.

## 8. Workflow manager: start with arrays, graduate to Snakemake

SLURM array scripts (above) are enough for ~150 embarrassingly-parallel
species-jobs and are easy to read. If the per-species logic fans into
multi-step dependencies (download → index → call → stratify → summarize)
with resumability needs, switch to **Snakemake** (Python-based, fits the
repo's style) with a SLURM profile — it gives automatic resumption on the
partial outputs in `tier3out/`. Recommendation: **start with arrays for
Phases 1–2**; adopt Snakemake if Phase 3's VGP chain gets fiddly. Don't
build the workflow manager before you need it.

## 9. Reproducibility & bring-it-back

- Every per-species script is deterministic given the manifest; the
  manifest records accessions + URLs; outputs are keyed by species. Anyone
  can re-run on a cluster from the repo + the manifest.
- `manifest.tsv` is committed (it's small text). Raw data is not.
- Result TSVs and figures are committed → they appear in the repo, the
  laptop pulls them, the manuscript cites them.
- The compute environment is pinned by `analysis/envs/tier3.yaml`; if
  containerizing, the Apptainer recipe is committed too.

## 10. Failure handling & checkpoints

- **Idempotent per-species outputs**: each job writes one TSV and skips if
  it exists (`[ -s "$OUT" ] && exit 0`). Re-running the array after a
  partial failure only redoes the missing species. No wasted compute.
- **Download robustness**: `wget -c` (resume); retry wrapper for NCBI/ENA
  (they rate-limit). Stage to scratch, verify MD5 against manifest.
- **Annotation gaps** (the main real failure mode for 4D): the script
  degrades gracefully — emits GC3-only rows with `notes="no_CDS"` so the
  merge still includes the species for composition. Nothing crashes the
  array; one species's bad annotation never blocks another.
- **VCF-too-big**: downsample to fixed n=20 (`bcftools view -S
  .,<comma-separated-20-samples>` or random-sample) and record n per
  species. π is stable at n=20; SFS-B needs more but is optional.

## 11. Open decisions (need your call before Phase 0)

1. **Which cluster, and is it SLURM?** (§2) — determines the job-script
   flavor and whether to use conda or Apptainer.
2. **Internet on compute nodes?** If only login nodes have egress, we
   *stage* everything to scratch first (fine for 3a/3c; for Ag1000G we
   then pre-download the Zarr subset rather than stream).
3. **Downsample n for population VCFs** — propose 20 (enough for π; SFS-B
   optional and needs more — flag if you want it).
4. **Do we build Typst + Lean on the cluster too?** Not needed for Tier 3
   (the manuscript/Lean edits stay laptop-side), but cheap if you want a
   fully-on-cluster workflow. Default: no.
5. **Congener substitution policy** (Tier 3a/3c when the exact species
   isn't assembled): allow, log every substitution; valid for GC3, not for
  π. Confirm this is the policy you want.

## 12. What this gets us, honestly

- Phase 1 alone delivers the **cross-species composition** picture (GC3
  vs Nₑ, ~60–80 species) on our own data — replacing the literature
  citation in the manuscript with a fit, and testing the across-clade
  attenuation.
- Phases 2a–c deliver the **plateau** on data: within-clade Drosophila
  (B∝Nₑ should hold → GC3 and π_S/π_W saturate together) and cross-clade
  Drosophila-vs-Anopheles (the within-recombination-class test). This is
  the discriminating Tier-3 signature.
- Phase 3 delivers the **rising limb** (VGP vertebrates) with individual
  π + within-individual W/S stratification, the part VGP can actually do.
- Together: the full within-vs-across-clade prediction the §4
  homology-maintenance mechanism makes, on data we collect rather than
  cite. Not the whole story (low-Nₑ rising limb stays Buffalo's
  Nₑ/Nₑ-reduction territory) — held.
