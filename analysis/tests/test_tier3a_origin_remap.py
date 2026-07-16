import argparse
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "tier3a_origin_remap", ROOT / "analysis/tier3a_origin_remap.py"
)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class Tier3AOriginRemapTests(unittest.TestCase):
    def test_native_mapping_audit_observes_one_to_one_without_rewriting(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            h1_fai, h2_fai = root / "h1.fai", root / "h2.fai"
            paf, output, contig_map = root / "map.paf", root / "audit.json", root / "contigs.tsv"
            h1_fai.write_text("h1\t1000\t0\t80\t81\n", encoding="utf-8")
            h2_fai.write_text("h2\t1000\t0\t80\t81\n", encoding="utf-8")
            paf.write_text(
                "h1\t1000\t0\t400\t+\th2\t1000\t0\t400\t390\t400\t60\n"
                "h1\t1000\t500\t900\t+\th2\t1000\t500\t900\t390\t400\t60\n",
                encoding="utf-8",
            )
            original = paf.read_bytes()
            recheck = root / "recheck.paf"
            recheck.write_bytes(original)
            MODULE.audit_mapping(argparse.Namespace(
                paf=paf, h1_fai=h1_fai, h2_fai=h2_fai, output=output,
                contig_map=contig_map, native_recheck_paf=recheck, overlap_threshold=0.95,
            ))
            audit = json.loads(output.read_text(encoding="utf-8"))
            self.assertEqual(audit["observed_native_query_multiplicity_cap"], 1)
            self.assertEqual(audit["observed_native_target_multiplicity_cap"], 1)
            self.assertEqual(audit["audit_mode"], "pinned_native_1to1_recheck_read_only_no_production_replacement")
            self.assertEqual(paf.read_bytes(), original)
            self.assertIn("origin_main_native_1to1_mapping", contig_map.read_text(encoding="utf-8"))

    def test_native_mapping_audit_rejects_a_changed_native_recheck(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "h1.fai").write_text("h1\t1000\t0\t80\t81\n", encoding="utf-8")
            (root / "h2.fai").write_text("h2\t1000\t0\t80\t81\n", encoding="utf-8")
            (root / "map.paf").write_text(
                "h1\t1000\t0\t400\t+\th2\t1000\t0\t400\t390\t400\t60\n"
                "h1\t1000\t5\t405\t+\th2\t1000\t5\t405\t390\t400\t60\n",
                encoding="utf-8",
            )
            (root / "recheck.paf").write_text(
                "h1\t1000\t0\t400\t+\th2\t1000\t0\t400\t390\t400\t60\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "recheck changed"):
                MODULE.audit_mapping(argparse.Namespace(
                    paf=root / "map.paf", h1_fai=root / "h1.fai", h2_fai=root / "h2.fai",
                    output=root / "audit.json", contig_map=root / "contigs.tsv",
                    native_recheck_paf=root / "recheck.paf", overlap_threshold=0.95,
                ))

    def test_corrected_selector_rejects_legacy_binary(self):
        rows = []
        for dataset in MODULE.EXPECTED_DATASETS:
            rows.append({
                "dataset_id": dataset,
                "tier3a_correction_status": "current_production",
                "sweepga_origin_main_sha256": MODULE.LEGACY_SWEEPGA_SHA256,
                "sweepga_binary_sha256": MODULE.LEGACY_SWEEPGA_SHA256,
                "sweepga_origin_main_commit": MODULE.SWEEPGA_COMMIT,
                "sweepga_direct_command": f"sweepga H1 H2 {MODULE.REQUIRED_COMMAND_TOKEN}",
            })
        with self.assertRaisesRegex(ValueError, "unapproved SweepGA"):
            MODULE.validate_current_rows(rows)


if __name__ == "__main__":
    unittest.main()
