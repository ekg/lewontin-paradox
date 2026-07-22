#!/usr/bin/env python3
"""Build and validate the corrective VGP implementation-base manifest.

This module is deliberately repository-only: it never submits work, reads the
shared VGP payload tree, or mutates a prior result packet.  ``build`` freezes
the audit of the source branches while those refs are available.  ``check``
validates the frozen audit, the sidecar status ledger, and the implementation
contracts needed by the next clean canary.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import subprocess
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = ROOT / "analysis"
MANIFEST = ANALYSIS / "vgp_repair_base_manifest.json"
STATUS_LEDGER = ANALYSIS / "vgp_repair_base_artifact_status.tsv"

AUDITED_MAIN_HEAD = "0c20b121b25a8ff139e8c9704b6fa9f4f0de743f"
CHANNEL_COMMIT = "44bbfc24e4bcc48d0e3343cd3d83452721af8c36"

SELECTED = "SELECTED_IMPLEMENTATION"
SELECTED_HISTORICAL = "SELECTED_IMPLEMENTATION_WITH_HISTORICAL_OUTPUT"
SELECTED_SUPERSEDED = "SELECTED_IMPLEMENTATION_WITH_SUPERSEDED_OUTPUT"
HISTORICAL_ONLY = "HISTORICAL_OUTPUT_ONLY"
SKIPPED = "SKIPPED_CONCLUSION_ONLY"
DEPENDENCY = "DEPENDENCY_ALREADY_INTEGRATED"
MERGE = "MERGE_METADATA_NOT_REPLAYED"
DISPOSITIONS = {
    SELECTED, SELECTED_HISTORICAL, SELECTED_SUPERSEDED,
    HISTORICAL_ONLY, SKIPPED, DEPENDENCY, MERGE,
}


def _commits(text: str) -> list[str]:
    return text.split()


# The order in each list is the source-branch first-parent order.  These lists
# are intentionally explicit rather than discovered during validation: branch
# cleanup must not erase the integration decision.
LINEAGES = [
    {
        "task": "mirror-vgp-freeze1",
        "source_branch": "wg/agent-258/mirror-vgp-freeze1",
        "main_integration_commit": "7a727349b657b61cbf9a2fb29ef3cb76316c8c5a",
        "commits": _commits("66a3aaeee3f99fd9ef9cb20f31152fea535e9058"),
    },
    {
        "task": "resume-vgp-freeze1-mirror-real",
        "source_branch": "wg/agent-306/resume-vgp-freeze1-mirror-real",
        "main_integration_commit": "8e62848e11d7974f674a2bb24043a4fcb116ba2e",
        "commits": _commits("""
            7f0740b4fee5a09c37969b02831d503d42488039
            561adfec6dd57d7f2a89423e31763796aa4e06eb
            f865c692e88f33ba2f2574acaa5c1772b806edcb
            1465befaa79e47b4c79b6cf1f774610538a54d7d
            260a7001eb47df735bd472f93a9da97a8f388854
            603b825d7ff2aca0a65a9ee38afc52fb683e273e
            271da58532ff361c2f3039725cadcec357188ae8
        """),
    },
    {
        "task": "bgzip-vgp-freeze1-assemblies",
        "source_branch": "wg/agent-323/bgzip-vgp-freeze1-assemblies",
        "main_integration_commit": "47c9d91953dd6fcbd491c72973166251ee2dea83",
        "commits": _commits("""
            5c10bd0204a9a739ef9a8fefcaf4aff72e099da6
            3a99c9fb2cf7be2b3909d45a55089cca673a497a
            c86758731b0931a870b633ef0b50fb7a4328cf5d
        """),
    },
    {
        "task": "run-vgp-real-pilot",
        "source_branch": "wg/agent-318/run-vgp-real-pilot",
        "main_integration_commit": "a017809f44814f660093bec55b986b91daf9bcce",
        "commits": _commits("""
            60bc724bc86ccac2faf5dcaa6ad2f3798d4e9087
            23764c40fcccf5c99fabca958ba426c67b1f70e1
            9ca6df91abb9833c7466e491a0e878feaf9d5da8
            0c28617932aec60e4241cf63110d3c819173b628
            099a26c3e35c631a5af3c794354651e4b0c9b3d8
            7d8c9261d726f0c0ea5abcf9cce5ec57dbc1cdba
            d84506ea7f5715e6f3b58897ba82a3e765c34652
            eeb076241f3d7f86ba4ae4e10df3ef885ab4bb61
            17b7f5941de649d8b81df352dd226abcdd59c979
            d668a2af7dac75cfde4c7459409dbde1df2e129c
            010619e1b0123f3bdb3a8097d3c3d8dbed6ec560
            d1cb15c93a3b231e989e27bedf995df05385affd
            069cae8c87f210b1a28484bc05c471b69ebd44eb
            866e36c89f6c5a99c92e09ab11ed76e4dcae1190
            c6026eb33f5bcf4b383c00f4394ef53f433fd41d
            29f690b2e69a43cb7e7e5e655e050c89747ae30c
            6205dac2040fcd9386a11c204d137f1c99546d4f
            79a0297b9cdb231a5c9ce7467755338cef5e4f6f
            79e958e0ef2550e14120e478981f721e09b1097e
            d67e6cc31336c328a693d8a84e799a87f6ed67ce
            7eb0d31f78dc933b77384c6a035816cca88a8a91
            6684774317d645380cf1be39b8cbad888601c419
            4dff358fc306bd97da5c67a73850266740cb707e
            52e18707b3a0308de8c4b4e828039e484ca8686b
            9fee444da3ffac48ac4ee3d7ef46479fee925780
            357e816bb4085564cb3321a5d5b974936f6bf773
            6e6a0727aa331469ceba826e35c8c28cf5f8f401
            0c622ce32c5850ff9266eaf4cb9faee57f85a3fb
            642ae186936a6708f4069345d3ce9f48756b64ed
            306b5b0acab6db474d8bb9a1a6b1a1ebc8384e64
            e0f20710f2cf8124fa00fab5ce963f5fe54b1abf
            48aafe61c8525bf9d9115907b10ffe8dce81f50d
            910b402765bdab29d749c0742408c3f11d3ba2bb
            9f0cb7b5741a623118fa208a850a6117b926e2a3
            5c0a6c870c9a62e20c031f2eb5f38d205891724a
            0b3272e612ea20fc45b1d7ae1be04e24a849ccc8
            8d29d53482a0ebab73ea56e8fe6123f96cbb3f7e
            12337d73b1c674f8f3bd859bbdd5c0ec2dd6fd5f
            636575319d6d9ff2db6352e5772b5059d186411d
            f9666613d7d75e2b4ecd8bf9b2e78e2e39290e54
            629b6e863e611fd1c84e97a12444b9ea34a6d023
            67c135348ae59d0fdabb4f9029088a944c573774
            593d8722ced588b4a533f4e79bb14d89e321c364
            460f1ba76482957b10737310634c9df117b70204
            f19132baf637d9bb5251cfaed6485879855be0f1
            9f26e0c4e3354f424858c8ca4c8e7a05544930f0
            870715e509fe9ea01afaed6fba18ed103c025629
            82ce5838382ce88c35f7856c4f366bc2f3503bd3
            2acce079e78d114943b1735df61099deaaf53158
            ab1fa8cf7d82aeeee91ceedb5afc40fe8c6e1c8c
            fdd2e8a4eff270dc6dff697e7214a2e6af9a8833
        """),
    },
    {
        "task": "repair-vgp-psmc",
        "source_branch": "wg/agent-368/repair-vgp-psmc",
        "main_integration_commit": "ee451c91bdcabebdfcfd17d77bcce91f428f6361",
        "commits": _commits("""
            7ed1967de41ab02723fc2cdaf88dfad9c730619b
            aeadf1442b27f29b298873b19c5f5bc2df0e5015
            1e386d658f951d220584325dec9eba0802542b3a
        """),
    },
    {
        "task": "validate-vgp-pilot-reads",
        "source_branch": "wg/agent-365/validate-vgp-pilot-reads",
        "main_integration_commit": "d3ca516541f37698bd8afbb8bad64d6de54ba160",
        "commits": _commits("""
            811d8d01faa4e0df4753ae2fa62966f531d5e400
            02fd12121cab661b66f7396206a64f7960aa105a
            35f0fc178ac860df77fe2c3450a1f4b2625f042b
            273904bffeee69f8e213b5630d5f9472dab8703f
            cf4197f88157b8640f6d940d12f61f13181ce63c
            01a5a372b7f3c470b63c7da70bed79be2eacc62c
            92eeddbac565391f943fe9bc8531027a1dc5faec
            920f028efaa0f9ffb880da0761576021b3dac18a
            69f39d34e900bb2939fa325f4b9b097efaa679a4
            445efc179f8220fbdebcd41c41af1e8c4dbd9c53
            ba1ba9a999e0abd0c4b41744585350c6cab34910
            095c2fdd58a9e6758c472f6773eec9c379c3d3e3
            74cf49ccf4c88bc016a35d7c821d8b353efb0682
            b86bcce00e5709b04da61973a1f85da3f94584d8
        """),
    },
    {
        "task": "scale-vgp-real",
        "source_branch": "wg/agent-371/scale-vgp-real",
        "main_integration_commit": "2044765ca50fc6526c7656da161bcc3f113d1276",
        "commits": _commits("""
            769aa5b5243ef8d32c1b8a5ffe79a954caa30575
            7b1e6c54ae7d590615a5ef105ea2c90eab50df49
            4e4ddae1ba6f4b0b615265bd83f32d8994e12b49
            277e87d27447f02fd9364ff0727aa9004708e140
            e384dc8ee34a95a665e7202d8a5c133a396a8973
            ae818339de304d4edd8853d8bd7e1715ce3092b8
            bfbe0d555db8b2d51bca9f37ba5d915466b9832f
            787d402f9d9c073087dcbba600e3253cff9d78cd
            73c9a955df4b3625d551835ee47e4a682141b0c7
        """),
    },
    {
        "task": "synthesize-vgp-real",
        "source_branch": "wg/agent-376/synthesize-vgp-real",
        "main_integration_commit": "0c20b121b25a8ff139e8c9704b6fa9f4f0de743f",
        "commits": _commits("""
            296ac54d96d62d2bd31e33e011024ec163c027eb
            36c0212c9e3efe90ed800eb1e6e92c387fb36b54
        """),
    },
]


DISPOSITION_OVERRIDES = {
    "7f0740b4fee5a09c37969b02831d503d42488039": DEPENDENCY,
    "561adfec6dd57d7f2a89423e31763796aa4e06eb": MERGE,
    "271da58532ff361c2f3039725cadcec357188ae8": MERGE,
    "66a3aaeee3f99fd9ef9cb20f31152fea535e9058": SELECTED_HISTORICAL,
    "f865c692e88f33ba2f2574acaa5c1772b806edcb": SELECTED_HISTORICAL,
    "1465befaa79e47b4c79b6cf1f774610538a54d7d": SELECTED_HISTORICAL,
    "260a7001eb47df735bd472f93a9da97a8f388854": SELECTED_HISTORICAL,
    "603b825d7ff2aca0a65a9ee38afc52fb683e273e": SELECTED_HISTORICAL,
    "c86758731b0931a870b633ef0b50fb7a4328cf5d": SELECTED_HISTORICAL,
    "fdd2e8a4eff270dc6dff697e7214a2e6af9a8833": HISTORICAL_ONLY,
    "aeadf1442b27f29b298873b19c5f5bc2df0e5015": HISTORICAL_ONLY,
    "445efc179f8220fbdebcd41c41af1e8c4dbd9c53": SKIPPED,
    "74cf49ccf4c88bc016a35d7c821d8b353efb0682": SKIPPED,
    "b86bcce00e5709b04da61973a1f85da3f94584d8": HISTORICAL_ONLY,
    "73c9a955df4b3625d551835ee47e4a682141b0c7": SELECTED_SUPERSEDED,
    "296ac54d96d62d2bd31e33e011024ec163c027eb": SKIPPED,
    "36c0212c9e3efe90ed800eb1e6e92c387fb36b54": SKIPPED,
}


ARTIFACTS = [
    ("analysis/vgp_freeze1_source_inventory.tsv", "HISTORICAL", "frozen mirror inventory"),
    ("analysis/vgp_freeze1_mirror_manifest.tsv", "HISTORICAL", "completed mirror output manifest"),
    ("analysis/vgp_freeze1_mirror_summary.json", "HISTORICAL", "completed mirror output summary"),
    ("analysis/vgp_freeze1_exception_ledger.json", "HISTORICAL", "mirror exception evidence"),
    ("analysis/vgp_freeze1_bgzf_inventory.tsv", "HISTORICAL", "BGZF conversion input snapshot"),
    ("analysis/vgp_freeze1_bgzf_manifest.tsv", "HISTORICAL", "completed BGZF output manifest"),
    ("analysis/vgp_freeze1_bgzf_summary.json", "HISTORICAL", "completed BGZF output summary"),
    ("analysis/vgp_real_canary_execution_v1.json", "HISTORICAL", "prior canary execution packet"),
    ("analysis/vgp_real_canary_promotion_v1.json", "HISTORICAL", "prior canary promotion packet"),
    ("analysis/vgp_real_pilot_P04_execution_v1.json", "HISTORICAL", "prior pilot execution packet"),
    ("analysis/vgp_real_pilot_closed_world_v1.json", "HISTORICAL", "prior ten-pair accounting snapshot"),
    ("analysis/vgp_real_pilot_fastga_scratch_v1.json", "HISTORICAL", "prior live scratch evidence"),
    ("analysis/vgp_real_pilot_sacct_v1.tsv", "HISTORICAL", "prior pilot scheduler telemetry"),
    ("analysis/vgp_psmc_bootstrap_repair_v1.json", "HISTORICAL", "repair diagnostic evidence, not new run authority"),
    ("analysis/vgp_read_validation_environment_v1.json", "HISTORICAL", "realized raw-validation environment"),
    ("analysis/vgp_read_validation_evidence_manifest_v1.json", "HISTORICAL", "raw measurement evidence index"),
    ("analysis/vgp_read_validation_results_v1.json", "HISTORICAL", "prior paired sensitivity measurements"),
    ("analysis/vgp_read_validation_per_pair_v1.tsv", "SUPERSEDED", "prior pair disposition is not a biological exclusion"),
    ("analysis/vgp_read_validation_report_v1.md", "SUPERSEDED", "prior interpretive decision is non-authoritative"),
    ("analysis/vgp_real_scaleout_v1/summary.json", "SUPERSEDED", "two-result accounting is not VGP scale-out"),
    ("analysis/vgp_real_scaleout_v1/pair_accounting.tsv", "SUPERSEDED", "technical outcomes are not biological exclusions"),
    ("analysis/vgp_real_scaleout_v1/results.md", "SUPERSEDED", "two-result scale-out conclusion is non-authoritative"),
    ("analysis/vgp_real_synthesis_v1/manifest.json", "SUPERSEDED", "conclusion packet excluded from corrective base"),
    ("analysis/vgp_real_synthesis_v1/paper_pairs.tsv", "SUPERSEDED", "P07 biological exclusion is not adopted"),
    ("analysis/vgp_real_synthesis_v1/claim_ledger.tsv", "SUPERSEDED", "conclusion claims are not adopted"),
    ("analysis/vgp_real_synthesis_v1/report.md", "SUPERSEDED", "conclusion report is not adopted"),
]


REQUIRED_CONTRACTS = {
    "analysis/fastga_scratch_guard.py": (
        "/proc/{pid}/cwd", "temp_environment_resolved", "managed_open_paths_resolved",
    ),
    "analysis/slurm/vgp_10_pilot/mapping_stage.sh": (
        "VGP_NODE_LOCAL_BASE:=/scratch", "mktemp -d -- \"$VGP_NODE_LOCAL_BASE/", "export TMPDIR=", "export TMP=", "export TEMP=", "fastga_scratch_guard", "sweepga",
    ),
    "analysis/slurm/scale_vgp_real/revalidate_p07_fastga.sh": (
        "node_local_base_resolved", "scratch_resolved", "refusing cleanup outside validated requested/resolved scratch roots",
    ),
    "analysis/vgp_10_pilot.py": (
        "freeze_psmcfa_bootstrap_units", "bootstrap_psmcfa", "primary_psmcfa_NKT_bins",
    ),
    "analysis/repair_vgp_psmc_bootstrap.py": (
        "bootstrap_manifest.tsv", "sampled_unit_indices", "primary_psmcfa_NKT_bins",
    ),
    "analysis/build_vgp_freeze1_bgzf.py": (
        "VGP_FREEZE1_ASSEMBLY_BGZF_ROOT", "RESOURCE_CLASSES", "scratch_multiplier",
    ),
    "analysis/guix/vgp-freeze1-bgzf-manifest.scm": ("samtools", "htslib", "python-pytest"),
    "analysis/slurm/vgp_10_pilot/pair_stage.sh": ("IMPG", "impg_hierarchical_lace.sh"),
    "analysis/slurm/vgp_10_pilot/impg_hierarchical_lace.sh": ("VGP_IMPG_LACE_THREADS", "IMPG"),
    "analysis/vgp_read_validation.py": ("paired sensitivity", "stream_depth_masks", "assembly_evidence"),
    "analysis/compile_vgp_read_validation_results.py": ("candidate_false_positive", "mask_sensitivity"),
    "analysis/guix/vgp_read_validation/manifest.scm": ("psmc-vgp-pinned", "minimap2", "jellyfish"),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git(*args: str) -> str:
    return subprocess.check_output(["git", *args], cwd=ROOT, text=True).strip()


def disposition(commit: str) -> str:
    return DISPOSITION_OVERRIDES.get(commit, SELECTED)


def reason_for(value: str) -> str:
    return {
        SELECTED: "implementation, environment, worker, test, or provenance patch retained",
        SELECTED_HISTORICAL: "implementation retained; emitted run packet is historical evidence only",
        SELECTED_SUPERSEDED: "implementation fix retained; emitted accounting/conclusion is superseded",
        HISTORICAL_ONLY: "output packet preserved without granting execution or conclusion authority",
        SKIPPED: "presentation or conclusion patch excluded from the corrective implementation base",
        DEPENDENCY: "prerequisite already represented by its own main integration",
        MERGE: "source-branch merge topology recorded but not replayed as a patch",
    }[value]


def build() -> None:
    lineages = []
    all_commits: set[str] = set()
    for spec in LINEAGES:
        patches = []
        for sequence, commit in enumerate(spec["commits"], 1):
            if commit in all_commits:
                raise ValueError(f"duplicate source commit: {commit}")
            all_commits.add(commit)
            value = disposition(commit)
            subject = git("show", "-s", "--format=%s", commit)
            patches.append({
                "sequence": sequence,
                "commit": commit,
                "subject": subject,
                "disposition": value,
                "reason": reason_for(value),
            })
        lineage = {key: value for key, value in spec.items() if key != "commits"}
        lineage["patches"] = patches
        lineages.append(lineage)

    payload = {
        "schema_version": "vgp-corrective-repair-base-v1",
        "task_id": "integrate-vgp-repair-base",
        "status": "CURRENT_IMPLEMENTATION_BASE",
        "audited_main_head": AUDITED_MAIN_HEAD,
        "pinned_guix_channel_commit": CHANNEL_COMMIT,
        "job_submission_count": 0,
        "selection_policy": {
            "implementation": "retain dependency-correct code, pinned environments, tests, workers, manifests, and provenance tooling",
            "historical_outputs": "preserve bytes and bind their status in vgp_repair_base_artifact_status.tsv",
            "forbidden_inferences": [
                "P07 is assembly-invalid",
                "technical pipeline failure is a biological exclusion",
                "two completed result packets constitute full VGP scale-out",
            ],
        },
        "logical_integration_order": [
            {"order": 1, "commit": "f83edbdc8c84c0993ac7b6a3e0239d734ea5e09f", "role": "canonical root and authorization contracts"},
            {"order": 2, "commit": "8e62848e11d7974f674a2bb24043a4fcb116ba2e", "role": "completed resumable mirror tooling and provenance"},
            {"order": 3, "commit": "47c9d91953dd6fcbd491c72973166251ee2dea83", "role": "BGZF resource classes, worker, indexes, and tests"},
            {"order": 4, "commit": "3aea32edb52a2bdf9b2eef7f1329dba743d5ff26", "role": "SweepGA/IMPG canary implementation prerequisite"},
            {"order": 5, "commit": "a017809f44814f660093bec55b986b91daf9bcce", "role": "general pilot, annotation, provenance, and scratch workflow"},
            {"order": 6, "commit": "ee451c91bdcabebdfcfd17d77bcce91f428f6361", "role": "PSMCFA bootstrap constructor repair"},
            {"order": 7, "commit": "2044765ca50fc6526c7656da161bcc3f113d1276", "role": "final /scratch alias, live guard, cleanup, and resource fixes; results superseded"},
            {"order": 8, "commit": "8a605e310d3ff55fc4960a0a667335d558b6b7a7", "role": "raw-read acquisition and immutable evidence inputs"},
            {"order": 9, "commit": "d3ca516541f37698bd8afbb8bad64d6de54ba160", "role": "raw-read paired sensitivity tooling and pinned environment; decisions superseded"},
        ],
        "observed_main_integration_order": [
            "f83edbdc8c84c0993ac7b6a3e0239d734ea5e09f",
            "8a605e310d3ff55fc4960a0a667335d558b6b7a7",
            "3aea32edb52a2bdf9b2eef7f1329dba743d5ff26",
            "8e62848e11d7974f674a2bb24043a4fcb116ba2e",
            "47c9d91953dd6fcbd491c72973166251ee2dea83",
            "a017809f44814f660093bec55b986b91daf9bcce",
            "f7e48cba9294ebf7a233dd4798e443f347bdf240",
            "ee451c91bdcabebdfcfd17d77bcce91f428f6361",
            "2044765ca50fc6526c7656da161bcc3f113d1276",
            "d3ca516541f37698bd8afbb8bad64d6de54ba160",
            "0c20b121b25a8ff139e8c9704b6fa9f4f0de743f",
        ],
        "integration_delta": {
            "analysis/build_vgp_freeze1_bgzf.py": "remove undeclared mkfifo/coreutils runtime while retaining one-pass streaming conversion",
            "analysis/build_vgp_repair_base.py": "freeze lineage selection, status binding, and fail-closed contract validation",
            "analysis/tests/test_vgp_repair_base.py": "regress claim boundaries, statuses, source enumeration, and required paths",
            "analysis/vgp_repair_base_artifact_status.tsv": "bind unchanged historical packet bytes to explicit status",
            "analysis/vgp_repair_base_handoff.md": "human-readable integration decision and order",
        },
        "historical_review_dependency": {
            "commit": "f7e48cba9294ebf7a233dd4798e443f347bdf240",
            "role": "identified bootstrap defect; review conclusions remain historical",
        },
        "skipped_main_conclusion": {
            "commit": "0c20b121b25a8ff139e8c9704b6fa9f4f0de743f",
            "role": "files preserved and status-bound, but synthesis claims are excluded from this base",
        },
        "lineages": lineages,
        "required_contracts": {key: list(value) for key, value in REQUIRED_CONTRACTS.items()},
        "artifact_status_ledger": str(STATUS_LEDGER.relative_to(ROOT)),
    }
    MANIFEST.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

    with STATUS_LEDGER.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("schema_version", "path", "status", "sha256", "reason"))
        for path_text, status, reason in ARTIFACTS:
            path = ROOT / path_text
            if not path.is_file():
                raise FileNotFoundError(path)
            writer.writerow(("vgp-repair-base-artifact-status-v1", path_text, status, sha256(path), reason))


def _require_tokens(path_text: str, tokens: Iterable[str]) -> None:
    path = ROOT / path_text
    if not path.is_file():
        raise AssertionError(f"missing required implementation path: {path_text}")
    content = path.read_text(errors="replace")
    missing = [token for token in tokens if token not in content]
    if missing:
        raise AssertionError(f"{path_text} is missing contract tokens: {missing}")


def check() -> dict[str, int]:
    data = json.loads(MANIFEST.read_text())
    assert data["schema_version"] == "vgp-corrective-repair-base-v1"
    assert data["status"] == "CURRENT_IMPLEMENTATION_BASE"
    assert data["job_submission_count"] == 0
    assert data["pinned_guix_channel_commit"] == CHANNEL_COMMIT
    assert len(data["selection_policy"]["forbidden_inferences"]) == 3

    seen: set[str] = set()
    disposition_counts: dict[str, int] = {value: 0 for value in DISPOSITIONS}
    for lineage in data["lineages"]:
        assert len(lineage["main_integration_commit"]) == 40
        for expected_sequence, patch in enumerate(lineage["patches"], 1):
            assert patch["sequence"] == expected_sequence
            assert len(patch["commit"]) == 40 and patch["commit"] not in seen
            assert patch["disposition"] in DISPOSITIONS
            seen.add(patch["commit"])
            disposition_counts[patch["disposition"]] += 1
    assert seen == {commit for lineage in LINEAGES for commit in lineage["commits"]}
    assert disposition_counts[SELECTED] > 0
    assert disposition_counts[SKIPPED] > 0
    assert disposition_counts[HISTORICAL_ONLY] > 0

    order = data["logical_integration_order"]
    assert [item["order"] for item in order] == list(range(1, len(order) + 1))
    observed_order = data["observed_main_integration_order"]
    assert len(observed_order) == 11 and observed_order[-1] == AUDITED_MAIN_HEAD
    git_checks = 0
    if shutil.which("git"):
        for item in order:
            if subprocess.run(
                ["git", "merge-base", "--is-ancestor", item["commit"], "HEAD"],
                cwd=ROOT, check=False,
            ).returncode:
                raise AssertionError(f"selected main integration is not on HEAD: {item['commit']}")
            git_checks += 1

    for path_text, tokens in REQUIRED_CONTRACTS.items():
        _require_tokens(path_text, tokens)

    rows = list(csv.DictReader(STATUS_LEDGER.open(newline=""), delimiter="\t"))
    assert len(rows) == len(ARTIFACTS)
    assert {row["status"] for row in rows} == {"HISTORICAL", "SUPERSEDED"}
    assert len({row["path"] for row in rows}) == len(rows)
    for row in rows:
        assert row["schema_version"] == "vgp-repair-base-artifact-status-v1"
        path = ROOT / row["path"]
        assert path.is_file(), row["path"]
        assert sha256(path) == row["sha256"], row["path"]

    repair = json.loads((ANALYSIS / "vgp_psmc_bootstrap_repair_v1.json").read_text())
    assert repair["sampling_population"] == "primary_psmcfa_NKT_bins"
    assert repair["passed"] is True
    bgzf = json.loads((ANALYSIS / "vgp_freeze1_bgzf_summary.json").read_text())
    assert bgzf["config_variable"] == "VGP_FREEZE1_ASSEMBLY_BGZF_ROOT"
    assert bgzf["closed_world"] is True

    return {
        "source_commits": len(seen),
        "selected_patches": sum(disposition_counts[value] for value in (SELECTED, SELECTED_HISTORICAL, SELECTED_SUPERSEDED)),
        "skipped_conclusion_patches": disposition_counts[SKIPPED],
        "historical_only_patches": disposition_counts[HISTORICAL_ONLY],
        "status_bound_artifacts": len(rows),
        "main_ancestry_checks": git_checks,
        "submitted_biological_jobs": 0,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("command", choices=("build", "check"))
    args = parser.parse_args()
    if args.command == "build":
        build()
    result = check()
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
