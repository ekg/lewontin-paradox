import csv
import math
from pathlib import Path

import numpy as np
import pytest

from analysis import tier3_fit as fit


def write_tsv(path: Path, rows):
    fields = list(rows[0])
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def composition(name, x, gc3="0.5", status="native"):
    return {
        "dataset_id": name.lower().replace(" ", ".") + ".tier3c",
        "scientific_name": name,
        "taxon_id": "1",
        "buffalo_diversity": "0.01",
        "buffalo_pred_log10_N": str(x),
        "assembly_accession_version": "GCF_000000001.1",
        "whole_genome_gc": "0.4",
        "whole_genome_callable_bases": "1000000",
        "gc3": gc3,
        "gc3_callable_third_positions": "20000" if gc3 else "",
        "callable_genes": "100" if gc3 else "",
        "annotation_status": status,
    }


def buffalo(name, x, clade="Insecta"):
    return {
        "species": name,
        "diversity": "0.01",
        "pred_log10_N": str(x),
        "kingdom": "Animalia",
        "phylum": "Arthropoda",
        "class": clade,
        "order": "Diptera",
        "family": "Exampleidae",
        "genus": name.split()[0],
        "ave_rec": "NA",
        "map_length": "NA",
    }


def test_exact_join_rejects_duplicate_and_congener(tmp_path):
    c_path, b_path = tmp_path / "c.tsv", tmp_path / "b.tsv"
    write_tsv(c_path, [composition("Drosophila melanogaster", 12)])
    write_tsv(b_path, [buffalo("Drosophila melanogaster", 12)])
    joined = fit.join_composition(fit.read_tsv(c_path), fit.read_buffalo(b_path))
    assert joined[0]["scientific_name"] == "Drosophila melanogaster"

    with pytest.raises(fit.SynthesisError, match="duplicate"):
        fit.join_composition(
            fit.read_tsv(c_path),
            fit.read_buffalo_rows(
                [buffalo("Drosophila melanogaster", 12)] * 2
            ),
        )

    write_tsv(b_path, [buffalo("Drosophila simulans", 12)])
    with pytest.raises(fit.SynthesisError, match="exact-species"):
        fit.join_composition(fit.read_tsv(c_path), fit.read_buffalo(b_path))


def test_join_checks_embedded_buffalo_covariate(tmp_path):
    c_path, b_path = tmp_path / "c.tsv", tmp_path / "b.tsv"
    write_tsv(c_path, [composition("Drosophila melanogaster", 11.0)])
    write_tsv(b_path, [buffalo("Drosophila melanogaster", 12.0)])
    with pytest.raises(fit.SynthesisError, match="predictor mismatch"):
        fit.join_composition(fit.read_tsv(c_path), fit.read_buffalo(b_path))


def test_native_annotation_provenance_is_carried_and_mismatch_rejected():
    source = composition("Drosophila melanogaster", 12)
    dataset_id = source["dataset_id"]
    annotation = {
        "provider": "NCBI RefSeq",
        "release": "2024-01-01",
        "assembly_accession": source["assembly_accession_version"],
        "fasta_sha256": "a" * 64,
        "gff_sha256": "b" * 64,
        "contig_mapping_sha256": "c" * 64,
        "sequence_regions_sha256": "d" * 64,
        "genetic_code": 1,
        "native_vs_projected": "native",
        "contig_dictionary_validated": True,
        "cds_audit": {
            "all_retained_cds_phase_translation_passed": True,
            "sampled_cds_mismatches": 0,
        },
    }
    provenance = {
        dataset_id: {
            "dataset_id": dataset_id,
            "reference": {
                "accession": source["assembly_accession_version"],
                "fasta_sha256": "a" * 64,
                "contig_dictionary_sha256": "d" * 64,
            },
            "annotation_provenance": annotation,
        }
    }
    joined = fit.join_composition([source], fit.read_buffalo_rows([buffalo("Drosophila melanogaster", 12)]), provenance)
    assert joined[0]["annotation_provider"] == "NCBI RefSeq"
    assert joined[0]["annotation_contig_mapping_sha256"] == "c" * 64
    assert joined[0]["annotation_genetic_code"] == 1

    provenance[dataset_id]["annotation_provenance"]["assembly_accession"] = "GCF_999999999.1"
    with pytest.raises(fit.SynthesisError, match="FASTA/GFF accession mismatch"):
        fit.join_composition([source], fit.read_buffalo_rows([buffalo("Drosophila melanogaster", 12)]), provenance)


def test_modality_separation_never_promotes_individual_het_to_population_pi(tmp_path):
    a = tmp_path / "a.tsv"
    row = {
        "dataset_id": "x",
        "scientific_name": "Species one",
        "modality": "direct_wfmash_h2_to_h1",
        "statistic_label": "individual_snv_heterozygosity",
        "individual_snv_heterozygosity": "0.01",
        "total_denominator": "10000000",
        "pi_S_over_pi_W": "0.8",
    }
    write_tsv(a, [row])
    observables = fit.load_diversity_observations(a, None)
    assert observables[0]["observable"] == "individual_snv_heterozygosity"
    assert observables[0]["observable_tier"] == "alignment_conditioned_individual"
    assert all(item["observable"] != "population_pi" for item in observables)


def test_empty_diversity_inputs_become_structured_missingness(tmp_path):
    a, b = tmp_path / "a.tsv", tmp_path / "b.tsv"
    a.write_text("dataset_id\tindividual_snv_heterozygosity\n", encoding="utf-8")
    b.write_text("dataset_id\tpopulation_pi\n", encoding="utf-8")
    claims = fit.diversity_claim_rows(fit.load_diversity_observations(a, b))
    by_observable = {row["observable"]: row for row in claims}
    assert by_observable["population_pi"]["status"] == "unavailable"
    assert by_observable["individual_snv_heterozygosity"]["status"] == "unavailable"
    assert by_observable["pi_S_over_pi_W"]["status"] == "unavailable"
    assert by_observable["polarized_sfs_B"]["status"] == "deferred"
    assert all(row["n"] == "0" for row in claims)


def test_non_native_gc3_is_missing_but_whole_genome_gc_remains(tmp_path):
    c_path, b_path = tmp_path / "c.tsv", tmp_path / "b.tsv"
    write_tsv(c_path, [composition("Species one", 4, "0.7", "projected")])
    write_tsv(b_path, [buffalo("Species one", 4)])
    joined = fit.join_composition(fit.read_tsv(c_path), fit.read_buffalo(b_path))
    assert joined[0]["gc3"] is None
    assert joined[0]["whole_genome_gc"] == pytest.approx(0.4)
    assert joined[0]["gc3_missing_reason"] == "annotation_not_native"


def test_clade_split_and_synthetic_slope_recovery():
    rows = []
    for clade, offset in [("Insecta", 0.1), ("Mammalia", 0.3)]:
        for i in range(8):
            x = float(i)
            rows.append(
                {
                    "scientific_name": f"{clade} species{i}",
                    "class": clade,
                    "genus": f"Genus{i // 2}",
                    "pred_log10_N": x,
                    "gc3": offset + 0.025 * x,
                    "whole_genome_gc": 0.4,
                    "gc3_denominator": 20000,
                }
            )
    splits = fit.split_clades(rows, min_n=8)
    assert set(splits) == {"Insecta", "Mammalia"}
    model = fit.fit_ols(rows, "gc3", fixed_effect="class", bootstrap_replicates=200)
    assert model.n == 16
    assert model.effect == pytest.approx(0.025, abs=1e-10)
    assert model.ci_low <= model.effect <= model.ci_high


def test_species_and_clade_leave_one_out_retain_negative_results():
    rows = [
        {
            "scientific_name": f"Species {i}",
            "class": "A" if i < 5 else "B",
            "genus": "G",
            "pred_log10_N": float(i),
            "gc3": 0.8 - 0.02 * i,
        }
        for i in range(10)
    ]
    species = fit.leave_one_out(rows, "gc3", unit="species")
    clades = fit.leave_one_out(rows, "gc3", unit="class")
    assert species["max"] < 0
    assert clades["max"] < 0


def test_synthetic_quadratic_concavity_recovery():
    rows = []
    for clade, offset in [("A", 0.0), ("B", 0.2)]:
        for i in range(10):
            x = i - 4.5
            rows.append(
                {
                    "scientific_name": f"{clade} species{i}",
                    "class": clade,
                    "pred_log10_N": x,
                    "gc3": 0.5 + offset + 0.02 * x - 0.004 * x * x,
                }
            )
    model = fit.fit_quadratic(rows, "gc3", fixed_effect="class", bootstrap_replicates=200)
    assert model.effect == pytest.approx(-0.004, abs=1e-10)
    assert model.ci_high < 0


def test_bh_policy_preserves_null_and_negative_claims():
    rows = [
        {"status": "estimated", "p_value": "0.01", "analysis_family": "primary_composition"},
        {"status": "estimated", "p_value": "0.9", "analysis_family": "primary_composition"},
        {"status": "estimated", "p_value": "0.04", "analysis_family": "primary_composition"},
    ]
    fit.apply_bh(rows)
    assert [float(row["q_value"]) for row in rows] == pytest.approx([0.03, 0.9, 0.06])
    assert len(rows) == 3


def test_results_round_trip_and_png_pdf_generation(tmp_path):
    rows = [fit.empty_result_row()]
    rows[0].update(
        row_kind="point",
        status="estimated",
        analysis_id="point.gc3.one",
        observable="gc3",
        observable_tier="exact_assembly_composition",
        scientific_name="Species one",
        clade="Insecta",
        predictor="pred_log10_N",
        predictor_value="12",
        estimate="0.55",
        n="1",
        eligible_n="1",
        missing_n="0",
        denominator="20000 callable CDS third positions",
    )
    result = tmp_path / "result.tsv"
    fit.write_results(rows, result)
    png, pdf = tmp_path / "figure.png", tmp_path / "figure.pdf"
    fit.render_figure(fit.read_tsv(result), png, pdf)
    assert png.read_bytes().startswith(b"\x89PNG\r\n\x1a\n")
    assert pdf.read_bytes().startswith(b"%PDF-")
    assert len(png.read_bytes()) > 1000
    assert len(pdf.read_bytes()) > 1000
