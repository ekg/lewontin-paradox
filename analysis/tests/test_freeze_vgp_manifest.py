from analysis import freeze_vgp_manifest as freeze


def test_compute_counts_uses_frozen_column_rules():
    rows = [
        {
            "Scientific Name": "Species one",
            "Status": "4",
            "Annotation status": "Completed NCBI",
            "Assembly IDs other high-quality haplotypes": "x1.hap2",
            "Lineage": "Fishes",
            "RefSeq annotation main haplotype": "GCF_000001.1",
        },
        {
            "Scientific Name": "Species two",
            "Status": "4",
            "Annotation status": "Ready",
            "Assembly IDs other high-quality haplotypes": "",
            "Lineage": "Fishes",
            "RefSeq annotation main haplotype": "GCF_000002.1",
        },
        {
            "Scientific Name": "Species three",
            "Status": "4",
            "Annotation status": "Wait",
            "Assembly IDs other high-quality haplotypes": "x3.hap2",
            "Lineage": "Mammals",
            "RefSeq annotation main haplotype": "",
        },
    ]
    counts = freeze.compute_counts(rows)
    assert counts["observed"] == {
        "unique_species": 3,
        "completed": 3,
        "completed_annotated": 2,
        "triple_eligible": 1,
        "triple_eligible_fish": 1,
        "completed_refseq_fish": 2,
    }


def test_build_seed_rows_merges_modalities_on_h1_accession():
    rows = [
        {
            "Status": "4",
            "Scientific Name": "Acipenser ruthenus",
            "Lineage": "Fishes",
            "Order": "f2",
            "Family Scientific Name": "Acipenseridae",
            "Annotation status": "Completed NCBI",
            "NCBI taxon ID": "7906",
            "Accession # for main haplotype": "GCA_902713425.2",
            "Accession #s other high-quality haplotypes": "GCA_902713435.2",
            "Assembly IDs other high-quality haplotypes": "fAciRut3.pat",
            "RefSeq annotation main haplotype": "GCF_902713425.2",
        },
        {
            "Status": "4",
            "Scientific Name": "Acipenser ruthenus",
            "Lineage": "Fishes",
            "Order": "f2",
            "Family Scientific Name": "Acipenseridae",
            "Annotation status": "Completed NCBI",
            "NCBI taxon ID": "7906",
            "Accession # for main haplotype": "GCA_902713425.2",
            "Accession #s other high-quality haplotypes": "",
            "Assembly IDs other high-quality haplotypes": "",
            "RefSeq annotation main haplotype": "GCF_902713425.2",
        },
    ]
    seeds = freeze.build_seed_rows(rows)
    assert len(seeds) == 1
    assert seeds[0].seed_modalities == (freeze.SEED_MODALITY_TIER3A, freeze.SEED_MODALITY_TIER3C)


def test_normalize_group_prefers_ncbi_class_then_manifest_lineage():
    assert freeze.normalize_group("Actinopterygii", "Other") == "Fishes"
    assert freeze.normalize_group("Mammalia", "Fishes") == "Mammals"
    assert freeze.normalize_group(None, "Birds") == "Birds"
    assert freeze.normalize_group(None, "Unknown") == "Other"


def test_select_pilot_spreads_categories_and_marks_selected():
    rows = [
        {
            "candidate_id": "fish_a",
            "catalog_row_number": 10,
            "assembly_composition_eligible": "yes",
            "assembly_diversity_eligible": "yes",
            "lineage_group": "Fishes",
            "class": "Actinopterygii",
            "combined_genome_bp": 1,
            "combined_contig_count": 10,
            "annotation_gff_uncompressed_bytes": 10,
            "predicted_core_hours_high": 1.0,
            "predicted_download_bytes_exact": 10,
        },
        {
            "candidate_id": "bird_b",
            "catalog_row_number": 20,
            "assembly_composition_eligible": "yes",
            "assembly_diversity_eligible": "no",
            "lineage_group": "Birds",
            "class": "Aves",
            "combined_genome_bp": 2,
            "combined_contig_count": 20,
            "annotation_gff_uncompressed_bytes": 20,
            "predicted_core_hours_high": 2.0,
            "predicted_download_bytes_exact": 20,
        },
        {
            "candidate_id": "mammal_c",
            "catalog_row_number": 30,
            "assembly_composition_eligible": "yes",
            "assembly_diversity_eligible": "yes",
            "lineage_group": "Mammals",
            "class": "Mammalia",
            "combined_genome_bp": 3,
            "combined_contig_count": 30,
            "annotation_gff_uncompressed_bytes": 30,
            "predicted_core_hours_high": 3.0,
            "predicted_download_bytes_exact": 30,
        },
    ]
    freeze.select_pilot(rows, limit=2)
    selected = [row["candidate_id"] for row in rows if row["pilot_selected"] == "yes"]
    assert len(selected) == 2
    assert "fish_a" in selected
    assert "mammal_c" in selected or "bird_b" in selected


def test_predicted_resources_uses_diversity_vs_composition_paths():
    diversity = {
        "assembly_diversity_eligible": "yes",
        "h1_length_bp": 1_200_000_000,
        "h2_length_bp": 1_200_000_000,
        "h1_contig_count": 500,
        "h2_contig_count": 500,
        "annotation_gff_uncompressed_bytes": 600_000_000,
        "h1_fasta_compressed_bytes": 100,
        "h2_fasta_compressed_bytes": 100,
        "annotation_gff_compressed_bytes": 50,
    }
    freeze.predicted_resources(diversity)
    assert diversity["predicted_core_hours_base"] > 0.0
    assert diversity["predicted_scratch_gb_high"] > diversity["predicted_scratch_gb_base"]

    composition = {
        "assembly_diversity_eligible": "no",
        "h1_length_bp": 1_200_000_000,
        "h2_length_bp": 0,
        "h1_contig_count": 500,
        "h2_contig_count": 0,
        "annotation_gff_uncompressed_bytes": 600_000_000,
        "h1_fasta_compressed_bytes": 100,
        "h2_fasta_compressed_bytes": 0,
        "annotation_gff_compressed_bytes": 50,
    }
    freeze.predicted_resources(composition)
    assert composition["predicted_core_hours_base"] < diversity["predicted_core_hours_base"]
    assert composition["predicted_download_bytes_exact"] == 150

