#!/usr/bin/env python3
"""Focused regression tests for the VGP annotation cataloger."""

from __future__ import annotations

import gzip
import tempfile
import unittest
from pathlib import Path

import catalog_vgp_annotations as catalog


class AnnotationParserTests(unittest.TestCase):
    def test_gff3_parser_validates_dictionary_phase_and_parents(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "fixture.gff3.gz"
            payload = (
                "##gff-version 3\n"
                "#!genome-build-accession NCBI_Assembly:GCA_000000001.1\n"
                "##sequence-region chr1 1 100\n"
                "chr1\tFixture\tgene\t1\t90\t.\t+\t.\tID=g1\n"
                "chr1\tFixture\tmRNA\t1\t90\t.\t+\t.\tID=t1;Parent=g1\n"
                "chr1\tFixture\tCDS\t1\t30\t.\t+\t0\tID=c1;Parent=t1\n"
            ).encode()
            with gzip.open(path, "wb") as handle:
                handle.write(payload)
            digest = catalog.sha256_file(path)
            result = catalog.parse_annotation_task(
                {
                    "physical_path": str(path),
                    "format": "GFF3",
                    "compression": "GZIP",
                    "expected_sha256": digest,
                }
            )
            self.assertEqual(result["feature_rows"], 3)
            self.assertEqual(result["feature_counts"], {"CDS": 1, "gene": 1, "mRNA": 1})
            self.assertEqual(result["declared_sequence_dictionary"], {"chr1": 100})
            self.assertEqual(result["missing_parent_references"], 0)
            self.assertEqual(result["cds_phase_errors"], 0)
            self.assertEqual(result["actual_sha256"], digest)

    def test_parser_reports_invalid_cds_phase_and_missing_parent(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "bad.gff"
            path.write_text("chr1\tFixture\tCDS\t1\t3\t.\t+\t.\tID=c1;Parent=absent\n", encoding="utf-8")
            result = catalog.parse_annotation_task(
                {
                    "physical_path": str(path),
                    "format": "GFF3",
                    "compression": "NONE",
                    "expected_sha256": catalog.sha256_file(path),
                }
            )
            self.assertEqual(result["cds_phase_errors"], 1)
            self.assertEqual(result["missing_parent_references"], 1)


class BindingTests(unittest.TestCase):
    def setUp(self) -> None:
        self.assembly = {
            "accession_version": "GCA_000000001.1",
            "derived_path": "/assembly.fa.gz",
            "fai_path": "/assembly.fa.gz.fai",
            "derived_bgzf_sha256": "a" * 64,
            "fai_sha256": "b" * 64,
        }

    def test_exact_dictionary(self) -> None:
        annotation = {
            "embedded_assembly_accession_version": "GCA_000000001.1",
            "observed_sequence_max_end": {"chr1": 90},
            "declared_sequence_dictionary": {"chr1": 100},
        }
        result = catalog.bind_annotation(annotation, self.assembly, {"chr1": 100}, {})
        self.assertEqual(result["binding_status"], "EXACT_DICTIONARY")

    def test_alias_requires_one_to_one_equal_length_mapping(self) -> None:
        annotation = {
            "embedded_assembly_accession_version": "GCA_000000001.1",
            "observed_sequence_max_end": {"GB1.1": 90},
            "declared_sequence_dictionary": {"GB1.1": 100},
        }
        result = catalog.bind_annotation(
            annotation, self.assembly, {"chr1": 100}, {"GB1.1": ("chr1", 100)}
        )
        self.assertEqual(result["binding_status"], "VALIDATED_ALIAS")
        self.assertEqual(result["alias_map"], {"GB1.1": "chr1"})

    def test_dictionary_length_mismatch_is_rejected(self) -> None:
        annotation = {
            "embedded_assembly_accession_version": "GCA_000000001.1",
            "observed_sequence_max_end": {"chr1": 90},
            "declared_sequence_dictionary": {"chr1": 101},
        }
        result = catalog.bind_annotation(annotation, self.assembly, {"chr1": 100}, {})
        self.assertEqual(result["binding_status"], "DICTIONARY_MISMATCH")
        self.assertEqual(result["length_mismatch_count"], 1)

    def test_explicit_circular_sequence_allows_origin_spanning_coordinate(self) -> None:
        annotation = {
            "embedded_assembly_accession_version": "GCA_000000001.1",
            "observed_sequence_max_end": {"chr1": 105},
            "declared_sequence_dictionary": {"chr1": 100},
            "circular_sequences": ["chr1"],
        }
        result = catalog.bind_annotation(annotation, self.assembly, {"chr1": 100}, {})
        self.assertEqual(result["binding_status"], "EXACT_DICTIONARY")
        self.assertEqual(result["coordinate_overrun_count"], 0)
        self.assertEqual(result["circular_coordinate_wrap_count"], 1)

    def test_embedded_accession_mismatch_precedes_dictionary(self) -> None:
        annotation = {
            "embedded_assembly_accession_version": "GCF_000000001.1",
            "observed_sequence_max_end": {"chr1": 90},
            "declared_sequence_dictionary": {"chr1": 100},
        }
        result = catalog.bind_annotation(annotation, self.assembly, {"chr1": 100}, {})
        self.assertEqual(result["binding_status"], "ACCESSION_MISMATCH")


if __name__ == "__main__":
    unittest.main()
