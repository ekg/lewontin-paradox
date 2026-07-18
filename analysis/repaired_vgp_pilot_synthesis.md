# Repaired VGP pilot: paper-oriented synthesis

Date: 2026-07-18 UTC

## Bottom line

**Audited outcome: `NO_GO`, correctly refused, `NOT_SUBMITTED`. No biological pilot ran.** The immutable decision SHA-256 is `9f39b13be5e0b1999c4cd98498399aee8700a7487d57fa81dfb7c59c29ff867d`. Acquisition stopped before a provider request and compute stopped before `sbatch`. The zero-byte, zero-job, empty-result outcome is evidence that the authorization boundary worked; it is not diversity, composition, or performance evidence.

The reviewed scope comprised **6 metadata candidates / 6 species / 6 classes**, one each from Mammalia, Aves, Lepidosauria, Amphibia, Actinopteri, Chondrichthyes. Metadata repair made all 6/6 composition candidates pre-download eligible, but the stricter integrated gate selected and executed 0/6. Diversity eligibility was 0/6, population-genomic eligibility 0/6, and demographic eligibility 0/6. Review QC was 57 PASS, 0 FAIL.

This synthesis uses only the independently reviewed evidence packet: exact candidate metadata in `analysis/vgp_pilot_manifest.tsv:2-7`; refusal/run/acquisition evidence in `analysis/vgp_pilot_run_manifest.tsv:2-5`, `analysis/vgp_pilot_acquisition_manifest.tsv:2-5`, `analysis/vgp_pilot_refusals.tsv:2`, and `analysis/vgp_pilot_slurm_telemetry.tsv:2`; the header-only immutable inventory; the excluded result and exclusions in `analysis/vgp_pilot_results.tsv:2` and `analysis/vgp_pilot_exclusions.tsv:2-5`; review decisions in `analysis/repaired_vgp_pilot_qc.tsv:2-58`; and the metadata-only demography audit in `analysis/vgp_demography_input_audit.tsv:2-7` plus its classified source ledger `analysis/vgp_independent_ne_sources.tsv:2-22`. It does not join older unreviewed pilot tables or treat metadata search results as executed outputs.

## Exact evidence obtained and measured denominators

| evidence layer | reviewed denominator | observed outcome | paper disposition |
| --- | ---: | --- | --- |
| repaired candidate manifest | 6 candidate rows | 6 exact current TaxId/name + exact-version H1 RefSeq + native exact-H1 annotation locations | metadata/provenance evidence only |
| taxonomic strata | 6 classes | one candidate in each listed class | no class was biologically analyzed |
| acquisition obligations | 12 finite obligations (H1 FASTA + native H1 annotation for each candidate) | 0 provider requests; 0 transferred biological bytes; 0 staged, quarantined, or promoted paths | zero-byte refusal evidence |
| immutable object inventory | 0 object rows | header-only inventory | empty inventory is retained as evidence; no local biological SHA-256 exists |
| executed candidate rows | 6 proposed; 0 selected/executed | no candidate result rows | do not report a biological sample size of six |
| run result summary | numerator 0; denominator 0; target 0 | one excluded `validated_species_count=0` summary with blank method and artifact hash | empty-result control row, not a biological estimate |
| exclusions | 4 rows | gate exclusion plus `CAP_MOOSEFS_READ_GB_EXCEEDED; CAP_SCRATCH_GIB_EXCEEDED; QUOTA_UNAVAILABLE` | all `imputed=false`, `demographic_input_used=false` |
| scheduler | 0 job IDs | `NOT_SUBMITTED`; no command, array, dependency, elapsed time, CPU, scratch, I/O, metadata operations, or network use | terminal refusal, not performance telemetry |

All candidate post-alignment fields—callable bases/fraction and queryable gene/base denominators—remain `POST_ALIGNMENT_REQUIRED`; they are **not zeros** and were never measured. No SweepGA mapping, IMPG partition/query, VCF/BCF, coding-diversity numerator, CDS/fourfold composition target, uncertainty interval, or candidate artifact SHA-256 was produced. The result uncertainty is therefore “not estimable because no biological measurement exists,” not a zero-width interval.

## Diversity and composition claims

Supported paper claims:

- The bounded repaired proposal produced no promoted diversity or composition estimate and no biological candidate result.
- The six exact metadata rows were composition-only candidates whose required denominators could only be measured after an authorized acquisition and alignment; that boundary was never crossed.
- The gate failed for exactly three reasons: `CAP_MOOSEFS_READ_GB_EXCEEDED; CAP_SCRATCH_GIB_EXCEEDED; QUOTA_UNAVAILABLE`. The proposed worst-case MooseFS read and scratch loads exceeded their strict limits, and enforceable quota/headroom evidence was absent.
- The refusal ledgers, zero-byte acquisition summary, header-only immutable inventory, excluded result summary, exclusion ledger, and `NOT_SUBMITTED` telemetry are valid control/provenance outcomes.

Unsupported paper claims:

- Any estimate, range, rank, class summary, or uncertainty interval for nucleotide diversity, coding diversity, CDS/fourfold composition, callability, queryability, or mapping multiplicity from this repaired pilot.
- Any statement that six species were biologically sampled, processed, or validated. Six is the metadata-candidate denominator; the executed denominator is zero.
- Any reuse of dormant toolchain smoke fixtures as a repaired-pilot biological result.
- Any inference that the zero resource row calibrates runtime, memory, scratch, I/O, storage, bandwidth, or throughput.

## Exact demography and independent-Ne audit

The join in `analysis/repaired_vgp_paper_table.tsv` is exact on candidate ID, scientific name, TaxId, H1 reference accession.version, BioSample, and individual/isolate. Literature populations remain separately identified and are never collapsed onto the VGP individual.

| species | exact VGP reference | VGP BioSample / individual | future PSMC eligibility | future MSMC2 eligibility | future SMC++ eligibility | valid independent absolute-Ne evidence |
| --- | --- | --- | --- | --- | --- | --- |
| *Camelus dromedarius* | `GCF_036321535.1` | `SAMN39296380` / `mCamDro1` | no | no | no | 6 LD-Ne population records |
| *Colius striatus* | `GCF_028858725.1` | `SAMN33339572` / `bColStr4` | no | no | no | 0 valid independent absolute-Ne records |
| *Candoia aspera* | `GCF_035149785.1` | `SAMN37159891` / `rCanAsp1` | no | no | no | 0 valid independent absolute-Ne records |
| *Dendropsophus ebraccatus* | `GCF_027789765.1` | `SAMN32145295` / `aDenEbr1` | no | no | no | 0 valid independent absolute-Ne records |
| *Lepisosteus oculatus* | `GCF_040954835.1` | `SAMN41155427` / `fLepOcu1` | no | no | no | 0 valid independent absolute-Ne records |
| *Heterodontus francisci* | `GCF_036365525.1` | `SAMN39432692` / `sHetFra1` | no | no | no | 0 valid independent absolute-Ne records |

All method decisions are currently `no`, but their blockers remain method-specific:

- **PSMC:** all 6/6 lack a heterozygosity-retaining callable diploid consensus, a compatible callable mask, and established coverage. H1 is explicitly haploid.
- **MSMC2:** all 6/6 lack validated accurate mutually comparable phased genomes/haplotypes, compatible masks, and audited population/individual relationships.
- **SMC++:** all 6/6 lack an exact-reference population genotype set plus compatible masks, population definitions, and method-specific QC. The camel candidate has a population VCF only on the incompatible `GCA/GCF_000803125` lineage; it is not ready for `GCF_036321535.1`.

VGP H1/H2 linkage is assembly discovery metadata only. It is never assumed to be a heterozygosity-retaining demographic genotype, two independent genomes, accurate comparable phasing, or a population dataset.

Only *Camelus dromedarius* has valid independent numeric Ne evidence in this bounded audit: 6 non-circular LD-Ne records across Awarik, Haddana, Majaheem, Sahliah, Shul, and Sofor (values 15, 11, 37, 24, 17, and 23; final per-breed sample sizes 5, 4, 9, 7, 4, and 5 diploid individuals). These are different animals/project from VGP `SAMN39296380/mCamDro1`; the exact LD time slice and intervals were not reported, so they are 6 population observations for 1 species, not six species-level replicates and not VGP-sample Ne.

Distinct non-promoted evidence classes are retained:

- `camel_psmc_fitak2020` is a published historical PSMC scenario on other animals/reference. Its approximate absolute Ne/time trajectory depends on the stated mutation-rate and generation-time scenario; it is not a new repaired-pilot inference and is not the contemporary independent LD-Ne field.
- `horn_theta_ima3_2022` is independent mitochondrial `theta=4Ne-mu`, retained only in a coalescent-scaled field. Without the appropriate locus mutation scenario—and because the estimand is mitochondrial rather than diploid nuclear—it is not absolute/contemporary Ne.
- `gar_nb_cosewic2015` is Nb, not Ne, and lacks value-to-population/time/uncertainty mapping. Camel and gar census measures and frog population structure are separate ecological/population fields, not Ne.
- Exactly one same-response `pi/(4mu)` policy row per candidate (6/6) is `circular_excluded`; no value was calculated. A predictor algebraically derived from the response pi cannot be used to explain that response.

Coalescent-scaled quantities, absolute scenario-scaled Ne/time, valid independent LD-Ne, census/Nb/structure evidence, and circular exclusions therefore remain separate in both the table and conclusions.

## Lewontin-paradox implications

The repaired refusal supports an infrastructure conclusion only: exact-reference metadata, measured denominators, immutable acquisition, storage headroom, and method-specific population inputs are material constraints on a comparative analysis. It supplies **no new cross-species pi or composition values**, so it cannot confirm, refute, narrow, or quantify Lewontin's paradox. The single-species independent LD-Ne inventory is insufficient for a cross-species Ne–pi relationship, and the circular `pi/(4mu)` rows are expressly unusable. Any biological statement about compressed diversity range, census size, linked selection, mutation, life history, or demographic history would exceed this evidence.

## Resource use, proposal, and uncertainty

Actual attributable use was zero: 0 provider requests, 0 biological bytes, 0 promoted objects, 0 Slurm submissions, 0 compute jobs, 0 core-seconds, 0 scratch bytes, 0 read/write/network bytes, and 0 demographic inferences. These zeros prove refusal/cap compliance only; they are **not performance telemetry**.

The following values are the gate's finite six-row **pre-run proposal**, not observed use:

| dimension | strict current ceiling | proposed finite value | gate disposition |
| --- | ---: | ---: | --- |
| `species` | 6 count | 6 count | within |
| `compressed_inputs_gib` | 120 GiB | 3.865305 GiB | within |
| `scratch_gib` | 139.698386192322 GiB | 205.5066 GiB | exceeded: refusal required |
| `core_hours` | 280 core-hours | 11.1353 core-hours | within |
| `concurrent_species` | 2 count | 2 count | within |
| `memory_per_job_gib` | 96 GiB | 32 GiB | within |
| `aggregate_wall_hours` | 17.5 hours | 8.7311 hours | within |
| `cpus_per_element` | 8 count | 8 count | within |
| `file_inodes` | 60000 count | 22596 count | within |
| `metadata_operations` | 500000 count | 271151 count | within |
| `moosefs_read_gb` | 180 GB | 406.7274 GB | exceeded: refusal required |
| `moosefs_write_gb` | 200 GB | 3.1634 GB | within |
| `persistent_input_gb` | 160 GB | 4.150340415 GB | within |
| `persistent_output_gb` | 32 GB | 0 GB | within |
| `peak_bandwidth_mib_s` | 120 MiB/s | 120 MiB/s | within |

Every stricter integrated cap wins. Thus the current ceiling is no more than 6 species, 120 GiB compressed inputs, **139.698386192 GiB scratch** (stricter than the outer 750 GiB bound), **280 core-hours** (stricter than 1,500), 2 concurrent species, and **96 GiB/job** (stricter than 256), with the additional recorded ceilings shown above. A ceiling is not a `GO`.

Observed-calibrated low/base/high projections for both a next wave and the full eligible catalog are deliberately blank in `analysis/repaired_vgp_next_decision.tsv`: there was no successful observation from which to estimate them. The zero-use refusal is excluded from calibration. “Not estimable” is the reviewed answer; no historical or planning proxy is silently relabeled as repaired-pilot telemetry. Any future numeric estimate is a planning artifact, not authorization.

## Active-state, retention, and quarantine verification

At synthesis time (2026-07-18 UTC), a read-only `squeue -h -u $USER` snapshot returned zero jobs, and a process-table search returned zero matching repaired VGP acquisition/run processes. This agrees with the immutable ledgers: zero job IDs, zero provider requests, zero network bytes, `NOT_SUBMITTED`, and no acquisition paths. **No attributable active jobs or download processes were found.**

The retained evidence is metadata/refusal material: gate and decision hashes, manifests, refusal/acquisition/run summaries, header-only immutable inventory, empty-result summary, exclusions, telemetry, QC, and the metadata-only demography audit. There are no attributable biological objects or failed partial paths to retain, quarantine, resume, or delete; acquisition rows have blank staging, quarantine, and promoted paths. No cleanup or deletion was performed by this synthesis.

For any future authorized work, the reviewed policy keeps failed partials for at most 14 days, reproducible intermediates through independent QC plus 90 days, and provenance/final results/failure ledgers/sentinels/telemetry through publication plus seven years, subject to institutional policy. Compressed sources persist through final reproducibility review and can be evicted only by checksum-listed human approval with tested rehydration. A new authorization must replace those general defaults with an explicit numeric storage and retention scope.

## Decision boundary: options only

This handoff defines exactly four options and selects or authorizes none:

1. **Repair remaining candidate metadata**—metadata-only identity, quota, checksum, annotation, and finite resource-envelope work; no biological acquisition, job, or inference.
2. **Request authorization for a numerically bounded expansion wave**—a new request must name exact species/accessions and numeric compressed-input, scratch, core-hour, concurrency, per-job-memory, storage, I/O, retention/expiry, and rollback limits. Full-catalog acquisition/execution is not covered.
3. **Request authorization for a population-data subset**—a new request must name exact taxon/reference/population/sample/assets and the same numeric resource/retention limits. Raw population bulk download, genotype construction, and each PSMC/MSMC2/SMC++ inference require separately explicit authorization.
4. **Stop**—retain this refusal and metadata audit as terminal evidence.

No expansion, full-catalog acquisition/execution, population bulk download, demographic inference, or ready executable task was created, authorized, or launched. Reporting a blank or hypothetical estimate grants no authority.
