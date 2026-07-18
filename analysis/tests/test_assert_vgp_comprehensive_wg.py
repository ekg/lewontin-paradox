import copy
import importlib.util
import os
from pathlib import Path
import unittest


MODULE_PATH = Path(__file__).parents[1] / "assert_vgp_comprehensive_wg.py"
SPEC = importlib.util.spec_from_file_location("assert_vgp_comprehensive_wg", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
assertions = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(assertions)


def live_graph_path() -> Path:
    project_root = os.environ.get("WG_PROJECT_ROOT")
    if project_root:
        return Path(project_root) / ".wg" / "graph.jsonl"
    return Path(__file__).parents[2] / ".wg" / "graph.jsonl"


class VgpComprehensiveWgAssertionsTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = live_graph_path()
        if not path.is_file():
            raise unittest.SkipTest(f"live WG graph is unavailable: {path}")
        cls.live_tasks = assertions.load_tasks(path)

    def test_live_quality_reviewed_graph_passes(self):
        self.assertEqual(assertions.validate_tasks(self.live_tasks), [])

    def test_pilot_run_cannot_bypass_acquisition(self):
        tasks = copy.deepcopy(self.live_tasks)
        tasks["run-vgp-10-pilot"]["after"].remove("acquire-vgp-10-pilot")

        errors = assertions.validate_tasks(tasks)

        self.assertTrue(any("implementation/acquisition join" in e for e in errors), errors)

    def test_mirror_timeout_regression_is_rejected(self):
        tasks = copy.deepcopy(self.live_tasks)
        tasks["mirror-vgp-freeze1"]["timeout"] = "4d"

        errors = assertions.validate_tasks(tasks)

        self.assertTrue(any("mirror-vgp-freeze1: timeout" in e for e in errors), errors)

    def test_annotation_cannot_become_a_core_gate(self):
        tasks = copy.deepcopy(self.live_tasks)
        description = tasks["scale-vgp-core"]["description"]
        tasks["scale-vgp-core"]["description"] = description.replace(
            "Annotation absence is never a core veto", "Annotation is required for core"
        )

        errors = assertions.validate_tasks(tasks)

        self.assertTrue(any("Annotation absence" in e for e in errors), errors)

    def test_direct_transmission_boundary_is_required(self):
        tasks = copy.deepcopy(self.live_tasks)
        description = tasks["pilot-pedigree-gbgc"]["description"]
        tasks["pilot-pedigree-gbgc"]["description"] = description.replace(
            "establish parent-of-origin, direction, and transmission",
            "use unpolarized H1/H2 states",
        )

        errors = assertions.validate_tasks(tasks)

        self.assertTrue(any("parent-of-origin" in e for e in errors), errors)

    def test_legacy_global_cap_is_rejected(self):
        tasks = copy.deepcopy(self.live_tasks)
        tasks["scale-vgp-core"]["description"] += " Hard ceiling: 6 species."

        errors = assertions.validate_tasks(tasks)

        self.assertTrue(any("forbidden legacy global cap" in e for e in errors), errors)

    def test_missing_machine_readable_deliverable_is_rejected(self):
        tasks = copy.deepcopy(self.live_tasks)
        description = tasks["review-vgp-10-pilot"]["description"]
        tasks["review-vgp-10-pilot"]["description"] = description.replace(
            "analysis/vgp_10_pilot_review_decision.json", "a review decision"
        )

        errors = assertions.validate_tasks(tasks)

        self.assertTrue(any("missing deliverable" in e for e in errors), errors)


if __name__ == "__main__":
    unittest.main()
