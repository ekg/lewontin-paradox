# VGP demography input audit

Date: 2026-07-18 UTC

## Outcome

None of the six repaired, prioritized metadata-eligible pilot candidates is currently eligible for PSMC, MSMC2, or SMC++. Eligibility is deliberately method-specific. A VGP H1 (and a linked H2 accession) is an assembly representation, not evidence for a heterozygosity-retaining diploid consensus, multiple independent genomes, a population sample, or a callable mask.

The audit denominator is the six rows in `analysis/vgp_pilot_manifest.tsv` produced by the repaired resolver. They remain `pilot_selected=no` only because `QUOTA_UNAVAILABLE` is stricter than their metadata eligibility; treating the literal execution flag as an empty denominator would repeat the obsolete fail-closed inventory rather than audit the repaired candidates.

| species | PSMC | MSMC2 | SMC++ | independent evidence |
| --- | --- | --- | --- | --- |
| *Camelus dromedarius* | no | no | no | six independent LD Ne points usable with missing-time/interval caveat; historical PSMC scenario separate |
| *Colius striatus* | no | no | no | missing after bounded exact-species search; absence not proof of nonexistence |
| *Candoia aspera* | no | no | no | missing after bounded exact-species search; absence not proof of nonexistence |
| *Dendropsophus ebraccatus* | no | no | no | no Ne located; population-structure evidence is not Ne |
| *Lepisosteus oculatus* | no | no | no | secondary Nb series located; value-to-population/time/uncertainty mapping incomplete |
| *Heterodontus francisci* | no | no | no | independent mitochondrial theta located; coalescent-scaled only, not absolute/contemporary Ne |

## Method contracts and blockers

PSMC requires one demonstrably diploid callable genome/consensus retaining heterozygosity, a compatible mask, and individual/read/reference/coverage provenance. Every H1 here is explicitly `haploid`; no compatible callable consensus/mask was located. Five exact BioSamples have public WGS metadata, but SRA experiment totals include alternate or duplicate representations and therefore do not establish deduplicated depth. *Colius striatus* has only RNA-seq linked to its exact BioSample.

MSMC2 requires validated accurate and mutually comparable phased haplotypes/genomes, masks, and individual/population relationships. Linked H2 accessions are recorded as discovery leads only. The composition-only resolver did not validate demographic phasing, masks, or the relationship set, and two assembly haplotypes are not multiple independent individuals.

SMC++ requires multi-sample population genotypes with sample identities, a population definition, exact reference compatibility, masks, and QC. No candidate meets all fields. For *Camelus dromedarius*, a 159.95 MB Dryad VCF and a 63-genome raw-read BioProject are real discovery leads, but the VCF uses the older GCA/GCF_000803125 reference lineage, not GCF_036321535.1; a compatible mask and SMC++-specific QC remain missing.

Every required status and blocker is explicit in `analysis/vgp_demography_input_audit.tsv`; blanks are not used to imply availability.

## Coalescent scale versus absolute scale

No inference was run. If eligible data are authorized later, coalescent-scaled output must be retained before any conversion. Absolute Ne and calendar time require an explicit mutation-rate and generation-time scenario. Only dromedary has a curated literature scenario here (1.1e-8 substitutions/site/generation and 5 years); it is a published assumption, not a direct estimate for mCamDro1. The horn-shark literature supplies mitochondrial theta=4Ne-mu and some HPDs, but no accepted absolute conversion; the mitochondrial estimand must not be silently treated as diploid nuclear Ne.

## Independent Ne, census, and circularity

The source table contains 6 independent dromedary LD-Ne point records (six named Saudi breeds). Their animals/project do not overlap the repaired VGP sample, but the exact recent generation bin and Ne intervals are not reported in the paper's Table 2, so those limitations travel with the values. A historical dromedary PSMC scenario is separate from contemporary Ne.

For spotted gar, the authoritative COSEWIC report gives six Nb values (26.9, 37.8, 50.1, 58.8, 61.4, 567.5), but it is a secondary citation and the passage does not map each value to a population, sampling year, sample size, or interval. Those values remain candidates, not promoted primary covariates. The same report supplies independent census estimates with 95% CIs for Point Pelee and Rondeau Bay. Horn-shark theta is coalescent-scaled; the frog source establishes population structure but no Ne. No exact-species Ne was located for *Colius striatus* or *Candoia aspera* in this bounded search, which is explicit missingness rather than proof of absence.

For every species, an audit-policy row excludes `pi/(4mu)` derived from the same response data. No such value was calculated. Estimates sharing response-generating sites or samples must also remain excluded unless an independent estimand is justified and provenance is recorded.

## Repository evidence and immutable cache

`analysis/vgp_demography_cache/index.json` records four metadata-only NCBI calls totaling 144272 response bytes: one batched OR search plus one batched summary for each of Taxonomy and SRA. Requests are keyed by SHA-256 of canonical normalized requests; every response has a SHA-256 and retrieval timestamp. The client enforces a 0.34-second minimum interval, bounded exponential backoff, Retry-After parsing, atomic writes, immutable conflicts, and a complete resume checkpoint with no pending keys. An offline rerun re-verifies every digest. The earlier exact-assembly cache remains separately cited for exact reference identity.

The cache authorization record states zero biological payloads, zero jobs, and no inference. Search-result HTML and literature pages were inspected, but no bulk reads, VCFs, assemblies, or masks were downloaded.

## Additional acquisition requiring new authorization

Five exact VGP BioSamples have all-WGS-run SRA `total_size` aggregates ranging from 50,238,621,306 to 478,930,410,963 bytes; per-species values are in the TSV. These are sizing disclosures, not approved downloads or reliable unique-read depth, because alternate representations can be present. *Colius striatus* WGS size remains unknown. Dromedary's candidate Dryad VCF is reported as 159.95 MB (exact bytes were not fetched). Any of these payloads, any population raw-read acquisition, reference remapping, genotype calling, mask construction, phasing validation, or demographic inference needs a separate authorization and resource plan.

## Bounded search ledger

- NCBI exact BioSample/TaxId batch: all six species; public SRA experiment metadata only.
- Exact-name repository/literature queries combined each species name with effective population size, population genomics, LD/temporal/pedigree/inbreeding, PSMC/MSMC/SMC++, VCF, and raw reads.
- Primary/authoritative sources retained: Ibrahim et al. 2025 (dromedary LD Ne and PRJNA1219399); Fitak et al. 2020 (independent historical PSMC scenario and Dryad VCF); Robertson et al. 2009 (frog population structure); Canfield et al. 2022 (horn-shark theta); COSEWIC 2015 (spotted-gar Nb candidates and census).
- No qualifying exact-species estimate was located for mousebird or viper boa. A future systematic review could expand databases and synonym searches but is outside this bounded pilot authorization.

## Reproduction and validation

```sh
analysis/slurm/guix_job.sh "$PWD/analysis/pilot_results/guix_environment.json" \
  python3 analysis/refresh_vgp_demography_metadata.py
analysis/slurm/guix_job.sh "$PWD/analysis/pilot_results/guix_environment.json" \
  python3 analysis/build_vgp_demography_audit.py
analysis/slurm/guix_job.sh "$PWD/analysis/pilot_results/guix_environment.json" \
  python3 analysis/validate_vgp_demography_audit.py
```

The first command is offline by default. `--refresh-cache` is only for a separately deliberate metadata refresh and never requests biological payloads.
