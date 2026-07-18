# Independent review of the repaired bounded VGP pilot

Date: 2026-07-18 UTC

## Audited outcome

**Review decision: PASS (correctly refused; `NOT_SUBMITTED`).** The exact repaired gate is `NO_GO`, with decision SHA-256 `6ba08bae6be0903b04cbf2b7a73c5201554ade82d2ea10a22280eff63289eb16`. Acquisition stopped before a provider request or biological byte, and compute stopped before `sbatch`. This is a valid audited control outcome, not an executed pilot row and not a biological result. QC totals are PASS=57, FAIL=0.

The three exact blockers are `CAP_MOOSEFS_READ_GB_EXCEEDED; CAP_SCRATCH_GIB_EXCEEDED; QUOTA_UNAVAILABLE` (`analysis/vgp_pilot_gate.json`; promoted run summary at `analysis/vgp_pilot_run_manifest.tsv:2-5`). No download, expansion, job, SweepGA/IMPG analysis, VCF/BCF generation, or demographic inference was launched by this review.

## Exact authorization boundary and refusal reproduction

I rebuilt the gate locally and matched its stable decision, authorization tuple, inputs, row audit, retrieval obligations, pair evidence, measurement contract, storage audit, environment, and strict cap vector. I then called both refusal entrypoints with downloader and submitter spies for current `NO_GO`, an unknown decision, an un-rehashed gate alteration, every bound manifest/root/environment/cap/retrieval-obligation/input/pair/measurement digest alteration, a relaxed species cap, and an altered approved retrieval URL. Every case returned refusal evidence with zero spy calls. The current refusal also reproduced all promoted acquisition, run, telemetry, result, exclusion, and refusal rows after removing timestamps, run IDs, and worktree prefixes.

The refusal matrix is recorded row-by-row in `analysis/repaired_vgp_pilot_qc.tsv`. The executable branch remains restricted to literal `GO` plus exact recomputation of all bound digests. A `NO_GO`, unknown token, altered gate, or changed bound contract cannot be reinterpreted as execution.

## Acquisition, identity, and immutable objects

The six metadata candidates in `analysis/vgp_pilot_manifest.tsv:2-7` have exact current taxon names/TaxIds, exact-version H1 RefSeq accessions, and official native annotation references equal to those H1 versions. The gate has two finite obligations per candidate (`h1_fasta` and `native_h1_annotation`) with official MD5 values, expected compressed sizes, mandatory staged local SHA-256, immediate local re-verification, and atomic read-only promotion.

However, `analysis/vgp_pilot_acquisition_manifest.tsv:2-5` contains only refusal/blocker rows and `analysis/vgp_pilot_immutable_object_inventory.tsv` contains zero objects. Therefore there are **no staged or promoted biological local SHA-256 values to claim or verify**. The review verifies the empty inventory and refusal SHA bindings; it does not substitute metadata checksums for nonexistent promoted-object hashes.

All six repaired rows resolve to `tier3c_composition`. Their linked H2 accessions share assembly BioSample/isolate metadata, but `h2_accession_version` is intentionally blank, same-individual and phase statuses are `not_applicable_composition_only`, and no Tier3A row is authorized. Thus the linked H2 records are discovery linkage, not validated phased H1/H2 evidence and not a diploid or population genotype dataset.

## Scientific outputs, denominators, and job evidence

There are no executed candidate rows. `analysis/vgp_pilot_results.tsv:2` is an excluded run summary, not a diversity/composition measurement; its run-level validated-species numerator, denominator, target, and value are all zero, while its measurement method and artifact SHA-256 are blank. `analysis/vgp_pilot_exclusions.tsv:2-5` reproduces the gate and three blocker exclusions with `imputed=false` and `demographic_input_used=false`.

Consequently SweepGA mappings, IMPG partitions/queries, callable/queryable gene/base denominators, target totals, VCF/BCF validity, and primary calculations are **not applicable for this refused run**. The dormant toolchain contracts remain locally pinned and digest-bound, but the pre-existing SweepGA/IMPG smoke artifacts are not misreported as outputs of this pilot. There are no promoted immutable biological artifacts from which a primary output could be recomputed.

The sole telemetry row (`analysis/vgp_pilot_slurm_telemetry.tsv:2`) is `NOT_SUBMITTED`, has no job/array/dependency/command, and records zero elapsed/CPU/scratch/I/O/metadata/network use. Thus every attributable job is terminal vacuously—there are zero job IDs—and compute-node retrieval was impossible on this path. This is evidence of refusal, not performance calibration.

## Strict caps, allocation, and headroom

| dimension | strict gate limit | six-row proposal | disposition |
| --- | ---: | ---: | --- |
| aggregate_wall_hours | 17.5 hours | 8.7311 | within proposal |
| compressed_inputs_gib | 120 GiB | 3.865305 | within proposal |
| concurrent_species | 2 count | 2 | within proposal |
| core_hours | 280 core-hours | 11.1353 | within proposal |
| cpus_per_element | 8 count | 8 | within proposal |
| file_inodes | 60000 count | 22596 | within proposal |
| memory_per_job_gib | 96 GiB | 32 | within proposal |
| metadata_operations | 500000 count | 271151 | within proposal |
| moosefs_read_gb | 180 GB | 406.7274 | proposal exceeded; refusal required |
| moosefs_write_gb | 200 GB | 3.1634 | within proposal |
| peak_bandwidth_mib_s | 120 MiB/s | 120 | within proposal |
| persistent_input_gb | 160 GB | 4.150340415 | within proposal |
| persistent_output_gb | 32 GB | 0 | within proposal |
| scratch_gib | 139.698386192 GiB | 205.5066 | proposal exceeded; refusal required |
| species | 6 count | 6 | within proposal |

The strict integrated limits are stronger than the task-wide ceilings where applicable: 6 species, 120 GiB compressed input, 139.698386192322 GiB scratch (stronger than 750 GiB), 280 core-hours (stronger than 1,500), 2 concurrent species, and 96 GiB per job (stronger than 256). The six-row proposal exceeded strict MooseFS read and scratch limits and lacked an enforceable quota/allocation report with required 25% headroom, so no executable selection was permitted. Filesystem free space is not treated as a user quota.

`analysis/repaired_vgp_resource_calibration.tsv` leaves low/base/high empty for every metric. Zero-use refusal is excluded from calibration, and an unsupported full-eligible-catalog projection is labeled `REQUIRES_NEW_AUTHORIZATION`.

## Exact demography join

The join uses candidate ID plus exact scientific name, TaxId, H1 reference version, BioSample, and individual/isolate. Literature populations remain separate from the VGP individual and from one another.

| taxon | exact VGP reference | VGP BioSample / individual | PSMC | MSMC2 | SMC++ | independent-Ne disposition |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Camelus dromedarius | `GCF_036321535.1` | `SAMN39296380` / `mCamDro1` | no | no | no | 6 independent LD-Ne populations |
| Colius striatus | `GCF_028858725.1` | `SAMN33339572` / `bColStr4` | no | no | no | none valid as absolute independent Ne |
| Candoia aspera | `GCF_035149785.1` | `SAMN37159891` / `rCanAsp1` | no | no | no | none valid as absolute independent Ne |
| Dendropsophus ebraccatus | `GCF_027789765.1` | `SAMN32145295` / `aDenEbr1` | no | no | no | none valid as absolute independent Ne |
| Lepisosteus oculatus | `GCF_040954835.1` | `SAMN41155427` / `fLepOcu1` | no | no | no | none valid as absolute independent Ne |
| Heterodontus francisci | `GCF_036365525.1` | `SAMN39432692` / `sHetFra1` | no | no | no | none valid as absolute independent Ne |

All three method decisions are independently `no` for all six candidates (`analysis/vgp_demography_input_audit.tsv:2-7`). PSMC lacks a heterozygosity-retaining callable diploid consensus/mask; MSMC2 lacks validated phased comparable haplotypes, masks, and relationships; SMC++ lacks an exact-reference population genotype set, population definition, masks, and QC. VGP H1/H2 assembly linkage satisfies none of those contracts.

Only *Camelus dromedarius* has valid independent numeric Ne observations in this bounded audit: six LD-Ne values for six named Saudi breed populations, from different animals/project, with an explicit missing-time/interval caveat. Historical camel PSMC is separate; spotted-gar secondary Nb lacks value-to-population/time/uncertainty mapping; horn-shark mitochondrial theta is coalescent-scaled, not absolute nuclear Ne; census and frog population-structure records are different estimands. One `pi/(4mu)` policy row per candidate is `circular_excluded`; no same-response pi-derived Ne is admitted as an independent predictor (`analysis/vgp_independent_ne_sources.tsv:2-22`).

## Reproducibility and authorization boundary

The pinned environment is `/gnu/store/x7vw3qvf5yxsff7x5cjhxs713m90ni6n-profile.drv` -> `/gnu/store/z9v2f6faha9cwjz0sm5iphhlzisgi077-profile`, channel commit `44bbfc24e4bcc48d0e3343cd3d83452721af8c36`, channels SHA-256 `45c055cd1d9010a72eacbb720037a22bccb2d8d6891dbd11b5d66365f29b3a17`, and manifest SHA-256 `2fb05e87aa2ac45ce51d4dcf93b232cb98627f525adace98357629ee3f15720a`. The full analysis test suite was run in the recorded profile realized by the pinned GNU Guix time-machine, and its result is recorded in QC.

Full catalog acquisition, raw population bulk download, population genotype construction, expansion, and all PSMC/MSMC2/SMC++ inference remain unauthorized. Any future resource projection or execution requires a new exact GO after the strict cap, allocation/headroom, immutable acquisition, and method-specific input contracts all pass.
