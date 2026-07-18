#!/usr/bin/env python3
"""Build the method-specific VGP demography metadata/literature audit."""

from __future__ import annotations

import csv
import json
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any, Iterable

from analysis.refresh_vgp_demography_metadata import PRIORITIZED_IDS, file_digest


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "analysis/vgp_pilot_manifest.tsv"
CACHE_INDEX = ROOT / "analysis/vgp_demography_cache/index.json"
RESOLUTION_INDEX = ROOT / "analysis/vgp_resolution_cache/index.json"
AUDIT = ROOT / "analysis/vgp_demography_input_audit.tsv"
SOURCES = ROOT / "analysis/vgp_independent_ne_sources.tsv"
REPORT = ROOT / "analysis/vgp_demography_audit.md"

AUDIT_HEADER = [
    "candidate_id", "scientific_name", "ncbi_taxid", "taxonomy_status", "taxonomy_source",
    "audit_denominator_basis", "exact_reference_accession_version", "exact_reference_status",
    "exact_reference_evidence", "biosample_accession", "individual_or_isolate_id",
    "assembly_haplotype_role", "linked_h2_accession_version", "h1_h2_demography_disposition",
    "raw_read_status", "raw_read_bioprojects", "raw_wgs_run_count",
    "raw_wgs_total_bases_metadata", "raw_wgs_total_size_bytes_metadata",
    "coverage_provenance_status", "callable_diploid_genome_status", "compatible_mask_status",
    "phasing_accuracy_status", "individual_population_relationship_status",
    "population_genotype_status", "population_dataset_id", "population_sample_size",
    "population_definition_status", "population_reference_status", "population_mask_status",
    "population_qc_status", "mutation_rate_status", "mutation_rate_scenario",
    "generation_time_status", "generation_time_scenario", "recombination_status",
    "psmc_eligible", "psmc_blockers", "msmc2_eligible", "msmc2_blockers",
    "smcpp_eligible", "smcpp_blockers", "coalescent_scaled_output_status",
    "absolute_ne_time_status", "independent_ne_status", "independent_ne_source_ids",
    "ecological_covariate_status", "census_measure_status", "circularity_disposition",
    "additional_acquisition_requires_authorization", "additional_acquisition_scope",
    "additional_acquisition_size", "metadata_cache_request_keys", "audit_provenance",
]

SOURCE_HEADER = [
    "record_id", "candidate_id", "scientific_name", "ncbi_taxid", "record_status",
    "classification", "independence_status", "circularity_status", "exclusion_reason",
    "estimand", "estimand_definition", "method", "population", "geography",
    "measurement_time", "sample_size", "value", "unit", "interval_lower", "interval_upper",
    "interval_type", "uncertainty_status", "mutation_rate", "generation_time",
    "response_dataset_overlap", "source_kind", "source_title", "source_authors", "source_year",
    "doi", "source_url", "source_locator", "source_retrieved_utc", "provenance_notes",
]

LITERATURE = {
    "camel_ld": {
        "source_kind": "primary_literature", "source_title": "Geographical distribution, genetic diversity, and environmental adaptations of dromedary camel breeds in Saudi Arabia",
        "source_authors": "Ibrahim et al.", "source_year": "2025", "doi": "10.3389/fvets.2024.1490186",
        "source_url": "https://www.frontiersin.org/journals/veterinary-science/articles/10.3389/fvets.2024.1490186/full",
        "source_locator": "Methods lines 398-400; Results Table 2 lines 428-437; Data availability lines 524-526",
    },
    "camel_psmc": {
        "source_kind": "primary_literature", "source_title": "Genomic signatures of domestication in Old World camels",
        "source_authors": "Fitak et al.", "source_year": "2020", "doi": "10.1038/s42003-020-1039-5",
        "source_url": "https://doi.org/10.1038/s42003-020-1039-5",
        "source_locator": "Figure 2 caption and demographic methods: generation time 5 y, mutation rate 1.1e-8; PSMC consensus/mask description and 100 bootstraps",
    },
    "gar_report": {
        "source_kind": "authoritative_government_status_report", "source_title": "COSEWIC assessment and status report on the Spotted Gar Lepisosteus oculatus in Canada",
        "source_authors": "COSEWIC", "source_year": "2015", "doi": "not_assigned", "source_url": "https://www.canada.ca/en/environment-climate-change/services/species-risk-public-registry/cosewic-assessments-status-reports/spotted-gar-2015.html",
        "source_locator": "Abundance lines 560-568 (census CI and secondary Nb values)",
    },
    "horn": {
        "source_kind": "primary_literature", "source_title": "Little Sharks in a Big World: Mitochondrial DNA Reveals Small-scale Population Structure in the California Horn Shark (Heterodontus francisci)",
        "source_authors": "Canfield, Galvan-Magana, and Bowen", "source_year": "2022", "doi": "10.1093/jhered/esac008",
        "source_url": "https://academic.oup.com/jhered/article/113/3/298/6570621",
        "source_locator": "Materials: N=318, 14 sites, 2004-2019; Table 4 IMa3 theta=4Ne-mu and HPD notes",
    },
    "frog": {
        "source_kind": "primary_literature", "source_title": "Discordant patterns of evolutionary differentiation in two Neotropical treefrogs",
        "source_authors": "Robertson, Duryea, and Zamudio", "source_year": "2009", "doi": "10.1111/j.1365-294X.2009.04126.x",
        "source_url": "https://www.csun.edu/~jrobertso/Publications_files/Robertson_Duryea_Zamudio_2009.pdf",
        "source_locator": "Figure 1 and Tables 2-3: 15 D. ebraccatus populations; mtDNA and microsatellite population structure; no Ne estimate",
    },
}


def tsv(path: Path, header: list[str], rows: Iterable[dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_manifest() -> dict[str, dict[str, str]]:
    with MANIFEST.open(newline="", encoding="utf-8") as handle:
        rows = {row["candidate_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    if set(rows) != set(PRIORITIZED_IDS):
        raise RuntimeError("audit denominator drifted from repaired bounded six")
    return rows


def sra_summary(index: dict[str, Any]) -> tuple[dict[str, dict[str, Any]], str]:
    entry = next(
        item for item in index["responses"]
        if item["request"]["parameters"].get("db") == "sra" and "id" in item["request"]["parameters"]
    )
    payload = json.loads((ROOT / entry["response_path"]).read_text(encoding="utf-8"))["result"]
    result: dict[str, dict[str, Any]] = {}
    for uid in payload["uids"]:
        item = payload[uid]
        root = ET.fromstring("<ROOT>" + item["expxml"] + item["runs"] + "</ROOT>")
        biosample = root.findtext("Biosample") or ""
        strategy = root.findtext(".//LIBRARY_STRATEGY") or "missing"
        stat = root.find(".//Statistics")
        assert stat is not None
        aggregate = result.setdefault(biosample, {"projects": set(), "strategies": set(), "wgs_runs": [], "wgs_bases": 0, "wgs_bytes": 0})
        aggregate["projects"].add(root.findtext("Bioproject") or "")
        aggregate["strategies"].add(strategy)
        if strategy == "WGS":
            aggregate["wgs_runs"].extend(run.get("acc") or "" for run in root.findall("Run"))
            aggregate["wgs_bases"] += int(stat.get("total_bases") or 0)
            aggregate["wgs_bytes"] += int(stat.get("total_size") or 0)
    return result, entry["request_key"]


def source(record_id: str, row: dict[str, str], **values: str) -> dict[str, str]:
    # TSV missingness is lexical and explicit; empty cells are never overloaded
    # to mean unavailable, not reported, or not applicable.
    result = {key: "not_applicable" for key in SOURCE_HEADER}
    result.update({"record_id": record_id, "candidate_id": row["candidate_id"], "scientific_name": row["scientific_name_source"], "ncbi_taxid": row["ncbi_taxid"], "source_retrieved_utc": "2026-07-18"})
    result.update(values)
    return result


def literature_sources(manifest: dict[str, dict[str, str]]) -> list[dict[str, str]]:
    by_name = {row["scientific_name_source"]: row for row in manifest.values()}
    camel = by_name["Camelus dromedarius"]
    rows: list[dict[str, str]] = []
    camel_populations = [
        ("awa", "Awarik", "coastal regions near Jeddah and South Jazan, Saudi Arabia", "5", "15"),
        ("had", "Haddana", "Hijaz Mountains and southwestern Saudi Arabia", "4", "11"),
        ("maj", "Majaheem", "Najd and Riyadh, Saudi Arabia", "9", "37"),
        ("sah", "Sahliah", "coastal regions near Jeddah and South Jazan, Saudi Arabia", "7", "24"),
        ("shu", "Shul", "Najd and Riyadh, Saudi Arabia", "4", "17"),
        ("sof", "Sofor", "Najd and Riyadh, Saudi Arabia", "5", "23"),
    ]
    for code, population, geography, n, value in camel_populations:
        rows.append(source(
            "camel_ld_" + code, camel, record_status="accepted_independent_with_time_caveat",
            classification="independent_literature_ne", independence_status="independent_different_animals_and_project",
            circularity_status="non_circular", estimand="Ne_LD", estimand_definition="LD-based effective population size reported by SNeP v1.11",
            method="linkage_disequilibrium_SNeP_1.11", population=population, geography=geography,
            measurement_time="recent LD time slice; exact generation bin not reported for Table 2 point",
            sample_size=n, value=value, unit="diploid individuals", uncertainty_status="interval not reported",
            mutation_rate="SNeP mutation correction enabled; numeric value not reported", generation_time="not required for reported generation-scaled point",
            response_dataset_overlap="none; PRJNA1219399 versus repaired VGP BioSample SAMN39296380",
            provenance_notes="Final per-population N after removal of related individuals; article reports 34 total final samples.", **LITERATURE["camel_ld"],
        ))
    rows.append(source(
        "camel_psmc_fitak2020", camel, record_status="accepted_historical_scenario_not_contemporary",
        classification="historical_coalescent_absolute_scenario", independence_status="independent_different_animals_and_project",
        circularity_status="non_circular", estimand="historical_Ne_PSMC", estimand_definition="PSMC historical diploid effective population size trajectory",
        method="PSMC_v0.6.4", population="sampled domestic dromedaries in PRJNA276064", geography="multiple origins; see source supplements",
        measurement_time="~700,000 to 16,000 years before present (reported narrative)", sample_size="study has 25 camel genomes across three species; exact dromedary n in supplement",
        value="~40000 declining to ~15000", unit="diploid individuals", interval_type="100 bootstrap trajectories; numeric CI not transcribed",
        uncertainty_status="bootstrap performed; pointwise bounds not extracted", mutation_rate="1.1e-8 changes/site/generation",
        generation_time="5 years", response_dataset_overlap="none known; older PRJNA276064 and GCA_000803125.2, not repaired sample/reference",
        provenance_notes="Absolute Ne/time are scenario-scaled. Coalescent output must remain separate if assumptions change.", **LITERATURE["camel_psmc"],
    ))
    rows.append(source(
        "camel_global_census_2025", camel, record_status="accepted_ecological_covariate", classification="census_measure",
        independence_status="independent_non_genomic", circularity_status="non_circular", estimand="global_census_size",
        estimand_definition="Approximate worldwide number of dromedary camels", method="compiled livestock census cited by primary article",
        population="global domestic dromedary", geography="global", measurement_time="source article current context (2025); underlying census year not stated",
        value="35000000", unit="animals", uncertainty_status="approximate; no interval", response_dataset_overlap="none",
        provenance_notes="Ecological/census covariate, not Ne.", **LITERATURE["camel_ld"],
    ))

    gar = by_name["Lepisosteus oculatus"]
    rows.append(source(
        "gar_nb_cosewic2015", gar, record_status="candidate_secondary_mapping_incomplete", classification="independent_literature_nb",
        independence_status="independent_different_animals", circularity_status="non_circular", estimand="Nb",
        estimand_definition="effective number of breeders; not interchangeable with generational Ne",
        method="microsatellite population-genetic estimate reported secondarily as Glass et al. 2015",
        population="six genetically distinct Canadian populations; value-to-population mapping unavailable in cited status-report text",
        geography="Point Pelee, Rondeau Bay, Long Point Bay / Lake Erie, Ontario, Canada",
        measurement_time="samples collected before 2015; exact cohort/year per estimate unavailable",
        sample_size="not reported in status-report passage", value="26.9;37.8;50.1;58.8;61.4;567.5", unit="breeders",
        uncertainty_status="intervals not reported in status-report passage", response_dataset_overlap="none known; Canadian microsatellite samples versus VGP fLepOcu1",
        provenance_notes="Authoritative secondary source. Do not promote to primary covariate until Glass et al. mapping, sample years, n, and uncertainty are recovered.", **LITERATURE["gar_report"],
    ))
    for record_id, population, value, low, high in [
        ("gar_census_point_pelee", "Point Pelee", "483", "433", "519"),
        ("gar_census_rondeau", "Rondeau Bay", "8121", "7281", "8278"),
    ]:
        rows.append(source(
            record_id, gar, record_status="accepted_ecological_covariate", classification="census_measure",
            independence_status="independent_non_genomic", circularity_status="non_circular", estimand="mature_adult_census_size",
            estimand_definition="mark-recapture abundance (Point Pelee) or habitat-density extrapolation (Rondeau Bay)",
            method="mark-recapture" if population == "Point Pelee" else "extrapolation from Point Pelee density",
            population=population, geography="Ontario, Canada", measurement_time="2007-2009 sampling context; report published 2015",
            value=value, unit="mature adults", interval_lower=low, interval_upper=high, interval_type="95% CI",
            uncertainty_status="reported 95% CI", response_dataset_overlap="none", provenance_notes="Census covariate, not Ne or Nb.", **LITERATURE["gar_report"],
        ))

    horn = by_name["Heterodontus francisci"]
    rows.append(source(
        "horn_theta_ima3_2022", horn, record_status="accepted_coalescent_scaled_only", classification="coalescent_scaled_not_absolute_ne",
        independence_status="independent_different_animals", circularity_status="non_circular", estimand="theta=4Ne_mu",
        estimand_definition="IMa3 mitochondrial-locus population mutation parameter; not Ne without a justified mutation scenario",
        method="IMa3 on 724-bp mitochondrial control region", population="NCI, MLCA, CAT, BT, and LSIBM pairwise models",
        geography="California Channel Islands/mainland and Baja California", measurement_time="318 tissues collected 2004-2019; coalescent time depth",
        sample_size="318 total across 14 localities; model uses population subsets", value="NCI 1.18; MLCA 3.74/4.26; CAT 1.64/2.67/2.41; BT 3.16; LSIBM 5.12",
        unit="theta (4Ne-mu)", interval_type="95% HPD where bounded", uncertainty_status="some HPDs reported; several right tails unbounded",
        mutation_rate="not supplied as an absolute scaling scenario in Table 4", generation_time="not applicable to unscaled theta",
        response_dataset_overlap="none; mtDNA specimens differ from VGP sHetFra1", provenance_notes="Do not label these values absolute or contemporary Ne.", **LITERATURE["horn"],
    ))

    frog = by_name["Dendropsophus ebraccatus"]
    rows.append(source(
        "frog_population_structure_2009", frog, record_status="accepted_ecological_covariate_no_ne", classification="population_structure_covariate",
        independence_status="independent_different_animals", circularity_status="non_circular", estimand="population_structure",
        estimand_definition="mtDNA and microsatellite differentiation; no Ne was estimated", method="mtDNA plus 12 microsatellite loci",
        population="15 populations", geography="Costa Rica and Panama", measurement_time="collection dates not recovered in bounded audit",
        sample_size="mtDNA per-site n=4-40; microsatellite multilocus sampling described in source", value="six inferred demes", unit="demes",
        uncertainty_status="not an Ne estimate", response_dataset_overlap="none known", provenance_notes="Useful population-definition evidence only; not a genome-wide genotype set for SMC++.", **LITERATURE["frog"],
    ))

    for name in ["Colius striatus", "Candoia aspera"]:
        row = by_name[name]
        rows.append(source(
            "no_independent_ne_" + row["candidate_id"].split("_gca_")[0], row,
            record_status="missing_after_bounded_search", classification="missing_source_sentinel",
            independence_status="not_applicable", circularity_status="not_applicable", estimand="none located",
            estimand_definition="No exact-species LD, temporal, pedigree, variance, inbreeding, or well-provenanced literature Ne located in bounded audit",
            method="repository and exact-name literature search", population="missing", geography="missing", measurement_time="missing",
            sample_size="missing", value="missing", unit="missing", uncertainty_status="missing", response_dataset_overlap="not_applicable",
            source_kind="audit_search_record", source_title="Bounded exact-species source audit", source_year="2026",
            source_url="not_located", source_locator="search terms and disposition documented in analysis/vgp_demography_audit.md",
            provenance_notes="Explicit missingness; absence is not evidence no estimate exists.",
        ))

    for row in manifest.values():
        rows.append(source(
            "excluded_pi_div_4mu_" + row["candidate_id"], row, record_status="excluded_if_proposed",
            classification="excluded_circular", independence_status="not_independent", circularity_status="circular_excluded",
            exclusion_reason="algebraically derived from the same pi response and assumed mutation rate",
            estimand="pi/(4mu)", estimand_definition="hypothetical algebraic transformation of response-generating diversity; not an observed independent Ne",
            method="prohibited algebraic derivation", population="same as response", geography="same as response", measurement_time="undefined",
            sample_size="same response data", value="not calculated", unit="not applicable", uncertainty_status="not calculated",
            response_dataset_overlap="derived_from_response", source_kind="audit_exclusion_rule", source_title="Pre-specified circularity exclusion",
            source_authors="project audit", source_year="2026", source_url="not_applicable", source_locator="task circularity contract",
            provenance_notes="Policy row, not a literature estimate; no calculation was performed.",
        ))
    return rows


def method_rows(manifest: dict[str, dict[str, str]], cache: dict[str, Any], sources: list[dict[str, str]]) -> list[dict[str, str]]:
    sra, sra_key = sra_summary(cache)
    resolution = json.loads(RESOLUTION_INDEX.read_text(encoding="utf-8"))
    common_keys = [item["request_key"] for item in cache["responses"]]
    accepted_by_candidate: dict[str, list[str]] = {}
    independent_classes = {
        "independent_literature_ne", "independent_literature_nb",
        "historical_coalescent_absolute_scenario", "coalescent_scaled_not_absolute_ne",
    }
    for item in sources:
        if item["classification"] in independent_classes:
            accepted_by_candidate.setdefault(item["candidate_id"], []).append(item["record_id"])
    rows: list[dict[str, str]] = []
    for candidate_id in PRIORITIZED_IDS:
        source_row = manifest[candidate_id]
        sample = sra.get(source_row["biosample_accession"], {"projects": set(), "strategies": set(), "wgs_runs": [], "wgs_bases": 0, "wgs_bytes": 0})
        wgs = bool(sample["wgs_runs"])
        name = source_row["scientific_name_source"]
        row = {key: "" for key in AUDIT_HEADER}
        row.update({
            "candidate_id": candidate_id, "scientific_name": name, "ncbi_taxid": source_row["ncbi_taxid"],
            "taxonomy_status": "exact current NCBI scientific name and TaxId match",
            "taxonomy_source": "batched NCBI Taxonomy ESummary plus repaired manifest",
            "audit_denominator_basis": "repaired prioritized metadata-eligible six; pilot_selected remains no only because QUOTA_UNAVAILABLE",
            "exact_reference_accession_version": source_row["h1_accession_version"], "exact_reference_status": "exact version verified",
            "exact_reference_evidence": "analysis/vgp_resolution_cache/index.json and manifest h1_exact_version_status",
            "biosample_accession": source_row["biosample_accession"], "individual_or_isolate_id": source_row["individual_or_isolate_id"],
            "assembly_haplotype_role": source_row["h1_haplotype_role"], "linked_h2_accession_version": source_row["linked_h2_accessions_ncbi"],
            "h1_h2_demography_disposition": "assembly haplotypes only; not assumed diploid consensus, independent genomes, or population sample",
            "raw_read_status": ("public genomic WGS metadata tied to exact BioSample" if wgs else "no genomic WGS located for exact BioSample; RNA-Seq only"),
            "raw_read_bioprojects": ";".join(sorted(value for value in sample["projects"] if value)),
            "raw_wgs_run_count": str(len(sample["wgs_runs"])), "raw_wgs_total_bases_metadata": str(sample["wgs_bases"]) if wgs else "0",
            "raw_wgs_total_size_bytes_metadata": str(sample["wgs_bytes"]) if wgs else "0",
            "coverage_provenance_status": ("missing usable coverage: SRA totals include alternate/duplicate read representations and were not deduplicated" if wgs else "missing genomic reads and coverage"),
            "callable_diploid_genome_status": "missing; H1 is explicitly haploid and no heterozygosity-retaining consensus was located",
            "compatible_mask_status": "missing; no callable mask tied to an audited diploid consensus/reference",
            "phasing_accuracy_status": "missing for demography; linked H2 accession is not validated here as an accurate comparable phased haplotype",
            "individual_population_relationship_status": "single assembly BioSample/individual known; population membership and relatedness missing",
            "population_genotype_status": "missing multi-sample exact-reference VCF/genotype dataset",
            "population_dataset_id": "missing", "population_sample_size": "missing", "population_definition_status": "missing",
            "population_reference_status": "missing", "population_mask_status": "missing", "population_qc_status": "missing",
            "mutation_rate_status": "missing exact-species direct estimate; no absolute scaling authorized",
            "mutation_rate_scenario": "missing", "generation_time_status": "missing exact audit scenario",
            "generation_time_scenario": "missing", "recombination_status": "missing map/rate compatible with exact reference",
            "psmc_eligible": "no", "psmc_blockers": "no diploid heterozygosity-retaining callable consensus; no compatible mask; coverage not established",
            "msmc2_eligible": "no", "msmc2_blockers": "no validated accurate mutually comparable phased genomes/haplotypes; no masks; population/individual relationships missing",
            "smcpp_eligible": "no", "smcpp_blockers": "no exact-reference multi-sample population VCF; sample/population definitions, mask, and required QC missing",
            "coalescent_scaled_output_status": "not generated; inputs ineligible and inference unauthorized",
            "absolute_ne_time_status": "not generated; requires explicit mutation-rate and generation-time scenario after input eligibility and authorization",
            "independent_ne_status": "no eligible exact estimate located", "independent_ne_source_ids": ";".join(sorted(accepted_by_candidate.get(candidate_id, []))) or "none",
            "ecological_covariate_status": "no curated numeric covariate located in bounded audit", "census_measure_status": "no curated census measure located in bounded audit",
            "circularity_disposition": "same-response pi/(4mu) and shared response-generating data excluded; see excluded policy row",
            "additional_acquisition_requires_authorization": "yes",
            "additional_acquisition_scope": ("all public WGS runs tied to exact BioSample, followed by separately authorized diploid calling/masking validation" if wgs else "locate genomic raw reads or generate a diploid dataset, then size before authorization"),
            "additional_acquisition_size": (str(sample["wgs_bytes"]) + " bytes total_size reported by SRA metadata; upper-bound-like aggregate may include alternate/duplicate representations" if wgs else "unknown; no WGS payload located"),
            "metadata_cache_request_keys": ";".join(common_keys),
            "audit_provenance": "analysis/vgp_demography_cache/index.json; analysis/vgp_resolution_cache/index.json",
        })
        if name == "Camelus dromedarius":
            row.update({
                "population_genotype_status": "candidate independent population VCF exists but is not exact-reference ready",
                "population_dataset_id": "Dryad 10.5061/dryad.prr4xgxj2 Drom.SNPs.filtered.vcf.gz; PRJNA1219399 raw 63-genome study",
                "population_sample_size": "25-genome Dryad study across camel species; 63 initial/34 post-QC in PRJNA1219399 study",
                "population_definition_status": "available in study supplements", "population_reference_status": "mismatch: GCA/GCF_000803125 lineage, not GCF_036321535.1",
                "population_mask_status": "missing compatible callable mask", "population_qc_status": "study QC reported; SMC++-specific missingness/callability validation not audited",
                "mutation_rate_status": "literature scenario available, not a direct estimate for repaired sample", "mutation_rate_scenario": "1.1e-8/site/generation (Fitak et al. 2020 PSMC scenario)",
                "generation_time_status": "literature scenario available", "generation_time_scenario": "5 years (Fitak et al. 2020 PSMC scenario)",
                "independent_ne_status": "six independent LD Ne points usable with missing-time/interval caveat; historical PSMC scenario separate",
                "ecological_covariate_status": "global livestock census available", "census_measure_status": "~35,000,000 global, census year/interval missing",
                "smcpp_blockers": "available VCF uses incompatible older reference; compatible mask and SMC++ QC missing; exact sample/population join requires verification",
                "additional_acquisition_scope": "option A: Dryad Drom.SNPs.filtered.vcf.gz only (metadata review/remap plan); option B: exact BioSample WGS; neither is authorized",
                "additional_acquisition_size": "Dryad reports 159.95 MB (exact bytes not retrieved); exact VGP BioSample WGS SRA total_size=%s bytes" % sample["wgs_bytes"],
            })
        elif name == "Lepisosteus oculatus":
            row.update({
                "independent_ne_status": "secondary Nb series located; value-to-population/time/uncertainty mapping incomplete",
                "ecological_covariate_status": "mark-recapture abundance available", "census_measure_status": "Point Pelee 483 (95% CI 433-519); Rondeau Bay extrapolated 8121 (7281-8278)",
            })
        elif name == "Heterodontus francisci":
            row.update({
                "independent_ne_status": "independent mitochondrial theta located; coalescent-scaled only, not absolute/contemporary Ne",
                "coalescent_scaled_output_status": "published independent mtDNA theta=4Ne-mu available; no new output generated",
                "absolute_ne_time_status": "published theta not converted: locus mutation scenario absent and mtDNA estimand differs from diploid nuclear Ne",
            })
        elif name == "Dendropsophus ebraccatus":
            row.update({"ecological_covariate_status": "15-population mtDNA/microsatellite structure source available; no Ne", "independent_ne_status": "no Ne located; population-structure evidence is not Ne"})
        elif name in {"Colius striatus", "Candoia aspera"}:
            row["independent_ne_status"] = "missing after bounded exact-species search; absence not proof of nonexistence"
        # Preserve the inherited resolution request keys in addition to this audit cache.
        row["metadata_cache_request_keys"] += ";" + source_row["metadata_cache_request_keys"]
        rows.append(row)
    return rows


def render_report(rows: list[dict[str, str]], source_rows: list[dict[str, str]], cache: dict[str, Any]) -> str:
    cache_bytes = sum(item["response_size_bytes"] for item in cache["responses"])
    independent_count = sum(item["classification"] == "independent_literature_ne" for item in source_rows)
    return f"""# VGP demography input audit

Date: 2026-07-18 UTC

## Outcome

None of the six repaired, prioritized metadata-eligible pilot candidates is currently eligible for PSMC, MSMC2, or SMC++. Eligibility is deliberately method-specific. A VGP H1 (and a linked H2 accession) is an assembly representation, not evidence for a heterozygosity-retaining diploid consensus, multiple independent genomes, a population sample, or a callable mask.

The audit denominator is the six rows in `analysis/vgp_pilot_manifest.tsv` produced by the repaired resolver. They remain `pilot_selected=no` only because `QUOTA_UNAVAILABLE` is stricter than their metadata eligibility; treating the literal execution flag as an empty denominator would repeat the obsolete fail-closed inventory rather than audit the repaired candidates.

| species | PSMC | MSMC2 | SMC++ | independent evidence |
| --- | --- | --- | --- | --- |
""" + "\n".join(
        f"| *{r['scientific_name']}* | {r['psmc_eligible']} | {r['msmc2_eligible']} | {r['smcpp_eligible']} | {r['independent_ne_status']} |" for r in rows
    ) + f"""

## Method contracts and blockers

PSMC requires one demonstrably diploid callable genome/consensus retaining heterozygosity, a compatible mask, and individual/read/reference/coverage provenance. Every H1 here is explicitly `haploid`; no compatible callable consensus/mask was located. Five exact BioSamples have public WGS metadata, but SRA experiment totals include alternate or duplicate representations and therefore do not establish deduplicated depth. *Colius striatus* has only RNA-seq linked to its exact BioSample.

MSMC2 requires validated accurate and mutually comparable phased haplotypes/genomes, masks, and individual/population relationships. Linked H2 accessions are recorded as discovery leads only. The composition-only resolver did not validate demographic phasing, masks, or the relationship set, and two assembly haplotypes are not multiple independent individuals.

SMC++ requires multi-sample population genotypes with sample identities, a population definition, exact reference compatibility, masks, and QC. No candidate meets all fields. For *Camelus dromedarius*, a 159.95 MB Dryad VCF and a 63-genome raw-read BioProject are real discovery leads, but the VCF uses the older GCA/GCF_000803125 reference lineage, not GCF_036321535.1; a compatible mask and SMC++-specific QC remain missing.

Every required status and blocker is explicit in `analysis/vgp_demography_input_audit.tsv`; blanks are not used to imply availability.

## Coalescent scale versus absolute scale

No inference was run. If eligible data are authorized later, coalescent-scaled output must be retained before any conversion. Absolute Ne and calendar time require an explicit mutation-rate and generation-time scenario. Only dromedary has a curated literature scenario here (1.1e-8 substitutions/site/generation and 5 years); it is a published assumption, not a direct estimate for mCamDro1. The horn-shark literature supplies mitochondrial theta=4Ne-mu and some HPDs, but no accepted absolute conversion; the mitochondrial estimand must not be silently treated as diploid nuclear Ne.

## Independent Ne, census, and circularity

The source table contains {independent_count} independent dromedary LD-Ne point records (six named Saudi breeds). Their animals/project do not overlap the repaired VGP sample, but the exact recent generation bin and Ne intervals are not reported in the paper's Table 2, so those limitations travel with the values. A historical dromedary PSMC scenario is separate from contemporary Ne.

For spotted gar, the authoritative COSEWIC report gives six Nb values (26.9, 37.8, 50.1, 58.8, 61.4, 567.5), but it is a secondary citation and the passage does not map each value to a population, sampling year, sample size, or interval. Those values remain candidates, not promoted primary covariates. The same report supplies independent census estimates with 95% CIs for Point Pelee and Rondeau Bay. Horn-shark theta is coalescent-scaled; the frog source establishes population structure but no Ne. No exact-species Ne was located for *Colius striatus* or *Candoia aspera* in this bounded search, which is explicit missingness rather than proof of absence.

For every species, an audit-policy row excludes `pi/(4mu)` derived from the same response data. No such value was calculated. Estimates sharing response-generating sites or samples must also remain excluded unless an independent estimand is justified and provenance is recorded.

## Repository evidence and immutable cache

`analysis/vgp_demography_cache/index.json` records four metadata-only NCBI calls totaling {cache_bytes} response bytes: one batched OR search plus one batched summary for each of Taxonomy and SRA. Requests are keyed by SHA-256 of canonical normalized requests; every response has a SHA-256 and retrieval timestamp. The client enforces a 0.34-second minimum interval, bounded exponential backoff, Retry-After parsing, atomic writes, immutable conflicts, and a complete resume checkpoint with no pending keys. An offline rerun re-verifies every digest. The earlier exact-assembly cache remains separately cited for exact reference identity.

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
analysis/slurm/guix_job.sh "$PWD/analysis/pilot_results/guix_environment.json" \\
  python3 analysis/refresh_vgp_demography_metadata.py
analysis/slurm/guix_job.sh "$PWD/analysis/pilot_results/guix_environment.json" \\
  python3 analysis/build_vgp_demography_audit.py
analysis/slurm/guix_job.sh "$PWD/analysis/pilot_results/guix_environment.json" \\
  python3 analysis/validate_vgp_demography_audit.py
```

The first command is offline by default. `--refresh-cache` is only for a separately deliberate metadata refresh and never requests biological payloads.
"""


def build() -> None:
    manifest = read_manifest()
    cache = json.loads(CACHE_INDEX.read_text(encoding="utf-8"))
    sources = literature_sources(manifest)
    rows = method_rows(manifest, cache, sources)
    tsv(AUDIT, AUDIT_HEADER, rows)
    tsv(SOURCES, SOURCE_HEADER, sources)
    REPORT.write_text(render_report(rows, sources, cache), encoding="utf-8")
    print(f"VGP_DEMOGRAPHY_AUDIT_BUILT species={len(rows)} sources={len(sources)}")


if __name__ == "__main__":
    build()
