import copy
import importlib.util
import os
from pathlib import Path
import unittest


MODULE_PATH = Path(__file__).parents[1] / "assert_vgp_pilot_repair_wg.py"
SPEC = importlib.util.spec_from_file_location("assert_vgp_pilot_repair_wg", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
assertions = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(assertions)


def live_graph_path() -> Path:
    project_root = os.environ.get("WG_PROJECT_ROOT")
    if project_root:
        return Path(project_root) / ".wg" / "graph.jsonl"
    return Path(__file__).parents[2] / ".wg" / "graph.jsonl"


class VgpPilotRepairWgAssertionsTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = live_graph_path()
        if not path.is_file():
            raise unittest.SkipTest(f"live WG graph is unavailable: {path}")
        cls.live_tasks = assertions.load_tasks(path)

    def test_live_reviewed_graph_and_definitions_pass(self):
        self.assertEqual(assertions.validate_tasks(self.live_tasks), [])

    def test_phantom_dependency_is_rejected(self):
        tasks = copy.deepcopy(self.live_tasks)
        tasks["repair-vgp-candidate"]["after"].append("phantom-vgp-gate")

        errors = assertions.validate_tasks(tasks)

        self.assertTrue(any("phantom dependencies" in error for error in errors), errors)

    def test_direct_acquisition_bypass_is_rejected(self):
        tasks = copy.deepcopy(self.live_tasks)
        assignment_edges = [
            dependency
            for dependency in tasks["acquire-repaired-vgp"]["after"]
            if dependency.startswith(".")
        ]
        tasks["acquire-repaired-vgp"]["after"] = assignment_edges + [
            "repair-vgp-candidate"
        ]

        errors = assertions.validate_tasks(tasks)

        self.assertTrue(any("acquisition bypasses" in error for error in errors), errors)

    def test_compute_altered_digest_refusal_text_is_required(self):
        tasks = copy.deepcopy(self.live_tasks)
        description = tasks["run-repaired-vgp"]["description"]
        tasks["run-repaired-vgp"]["description"] = description.replace(
            "altered local payload SHA-256", "changed payload"
        )

        errors = assertions.validate_tasks(tasks)

        self.assertTrue(
            any("altered local payload SHA-256" in error for error in errors), errors
        )

    def test_method_specific_demography_language_is_required(self):
        tasks = copy.deepcopy(self.live_tasks)
        description = tasks["audit-vgp-demography"]["description"]
        tasks["audit-vgp-demography"]["description"] = description.replace(
            "MSMC2:", "generic multiple-sequence method:"
        )

        errors = assertions.validate_tasks(tasks)

        self.assertTrue(any("MSMC2:" in error for error in errors), errors)


if __name__ == "__main__":
    unittest.main()
