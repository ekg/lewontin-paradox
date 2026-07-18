#!/usr/bin/env python3
"""Fail-closed validation for the gene-conversion evidence design.

The validator reads repository metadata only.  It performs no network access,
biological download, environment realization, or scheduler submission.
"""

from __future__ import annotations

import csv
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
PLAN = ROOT / "analysis/gene_conversion_evidence_plan.md"
DATASETS = ROOT / "analysis/gene_conversion_dataset_manifest.tsv"
ESTIMANDS = ROOT / "analysis/gene_conversion_estimand_manifest.tsv"
CLAIMS = ROOT / "analysis/gene_conversion_claim_matrix.tsv"

BRANCHES = {"direct", "population", "historical_phylogenetic", "non_allelic"}
EXPECTED_DATASETS = {
    "D01", "D02", "P01", "P02", "H01", "H02", "H03", "H04", "N01", "N02"
}
EXPECTED_ESTIMANDS = {
    "D_EVT", "D_TRACT", "D_CO", "D_GCBIAS", "P_SFS", "P_B",
    "H_SUB", "H_GBGC", "H_CLUSTER", "N_TRACT", "N_HOM",
}
EXPECTED_CLAIMS = {f"C{i:02d}" for i in range(1, 11)}

REQUIRED_ACCESSIONS = {
    "D01": ("PRJEB4500", "ERP003793", "GCF_000001735.4"),
    "D02": ("phs001224", "PRJEB3381", "PRJEB3246"),
    "P01": ("PRJEB31736", "PRJEB36890"),
    "P02": ("PRJNA273563",),
    "H01": (
        "GCA_048126635.1", "GCA_949316345.1", "GCA_964276395.1",
        "GCA_948146105.1", "GCA_048301445.1",
    ),
    "H02": (
        "GCA_017639655.1", "GCA_976974455.1", "GCA_965282525.1",
        "GCA_023634085.1", "GCA_963210335.1",
    ),
    "H03": (
        "GCA_039906515.1", "GCA_963455315.2", "GCA_937001465.1",
        "GCA_011762595.2", "GCA_949987515.2",
    ),
    "H04": (
        "GCA_048301445.1", "GCA_965643225.1", "GCA_901709675.2",
        "GCA_965263715.1", "GCA_949316345.1",
    ),
    "N01": ("PRJNA730822", "GCA_009914755.4"),
    "N02": ("PRJEB36100",),
}

# These identifiers must occur in the collection/assembly accession field, not
# merely in a citation URL or annotation description.  D01's GCF identifier is
# a reference rather than a study accession and is intentionally checked only
# by REQUIRED_ACCESSIONS above.
REQUIRED_EXACT_ACCESSION_FIELD = {
    dataset_id: tuple(
        accession for accession in accessions
        if not (dataset_id == "D01" and accession == "GCF_000001735.4")
    )
    for dataset_id, accessions in REQUIRED_ACCESSIONS.items()
}

REQUIRED_PLAN_PHRASES = (
    "four separate evidence branches",
    "Direct allelic evidence",
    "Population evidence",
    "Historical phylogenetic evidence",
    "Non-allelic evidence",
    "H1/H2 assemblies are state observations, not conversion observations",
    "Historical substitution bias is not a count of direct events",
    "Population B is model-dependent and is not a direct tract rate",
    "Paralog homogenization is not allelic gBGC",
    "bounded sensitivity prior",
    "gBGC can raise S alleles in frequency and fixation in a way that resembles positive selection",
    "population gBGC not measured",
    "non-allelic conversion not measured",
    "NOT_RUN_DESIGN_ONLY",
    "at least three ingroup taxa and two taxa",
    "at least three distinguishable copies",
    "downloaded no biological data",
    "submitted no local or scheduler analysis",
)


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{path}: missing header")
        return reader.fieldnames, list(reader)


def _check_rectangular_nonempty(
    path: Path, header: list[str], rows: list[dict[str, str]], errors: list[str]
) -> None:
    if len(header) != len(set(header)):
        errors.append(f"{path}: duplicate column name")
    for line, row in enumerate(rows, start=2):
        if None in row:
            errors.append(f"{path}:{line}: too many fields")
        missing = [name for name in header if not row.get(name, "").strip()]
        if missing:
            errors.append(f"{path}:{line}: empty fields {missing}")
        whitespace = [
            name for name in header
            if row.get(name, "") != row.get(name, "").strip()
        ]
        if whitespace:
            errors.append(f"{path}:{line}: outer whitespace in {whitespace}")


def validate_design() -> list[str]:
    errors: list[str] = []
    for path in (PLAN, DATASETS, ESTIMANDS, CLAIMS):
        if not path.is_file():
            errors.append(f"missing artifact: {path}")
    if errors:
        return errors

    dataset_header, datasets = read_tsv(DATASETS)
    estimand_header, estimands = read_tsv(ESTIMANDS)
    claim_header, claims = read_tsv(CLAIMS)
    for path, header, rows in (
        (DATASETS, dataset_header, datasets),
        (ESTIMANDS, estimand_header, estimands),
        (CLAIMS, claim_header, claims),
    ):
        _check_rectangular_nonempty(path, header, rows, errors)

    dataset_by_id = {row["dataset_id"]: row for row in datasets}
    if set(dataset_by_id) != EXPECTED_DATASETS:
        errors.append(f"dataset IDs differ: {sorted(dataset_by_id)}")
    if {row["branch"] for row in datasets} != BRANCHES:
        errors.append("dataset manifest does not contain exactly four evidence branches")
    expected_counts = {"direct": 2, "population": 2, "historical_phylogenetic": 4, "non_allelic": 2}
    for branch, expected in expected_counts.items():
        observed = sum(row["branch"] == branch for row in datasets)
        if observed != expected:
            errors.append(f"dataset branch {branch}: observed {observed}, expected {expected}")

    if dataset_by_id.get("D01", {}).get("execution_state") != "AUTHORIZED_DOWNSTREAM_PILOT":
        errors.append("D01 is not the sole selected direct executable pilot")
    if [row["dataset_id"] for row in datasets if row["branch"] == "direct" and row["execution_state"] == "AUTHORIZED_DOWNSTREAM_PILOT"] != ["D01"]:
        errors.append("direct executable selection must be exactly D01")
    phylo_executable = {
        row["dataset_id"] for row in datasets
        if row["branch"] == "historical_phylogenetic"
        and row["execution_state"] == "AUTHORIZED_DOWNSTREAM_PILOT"
    }
    if phylo_executable != {"H01", "H02"}:
        errors.append(f"phylogenetic executables differ: {sorted(phylo_executable)}")
    for row in datasets:
        if row["branch"] in {"population", "non_allelic"} and row["execution_state"] != "NOT_RUN_DESIGN_ONLY":
            errors.append(f"{row['dataset_id']}: unauthorized design-only branch state")

    for dataset_id, needles in REQUIRED_ACCESSIONS.items():
        text = "\t".join(dataset_by_id.get(dataset_id, {}).values())
        for needle in needles:
            if needle not in text:
                errors.append(f"{dataset_id}: missing accession {needle}")
    for dataset_id, needles in REQUIRED_EXACT_ACCESSION_FIELD.items():
        accession_field = dataset_by_id.get(dataset_id, {}).get("exact_accessions", "")
        for needle in needles:
            if needle not in accession_field:
                errors.append(f"{dataset_id}: missing accession {needle} from exact_accessions")
    for row in datasets:
        urls = row["source_urls"].split(";")
        if not urls or any(not url.startswith("https://") for url in urls):
            errors.append(f"{row['dataset_id']}: source URLs must be explicit HTTPS URLs")

    estimand_by_id = {row["estimand_id"]: row for row in estimands}
    if set(estimand_by_id) != EXPECTED_ESTIMANDS:
        errors.append(f"estimand IDs differ: {sorted(estimand_by_id)}")
    if {row["branch"] for row in estimands} != BRANCHES:
        errors.append("estimand manifest does not contain exactly four evidence branches")
    for estimand_id in ("P_SFS", "P_B", "N_TRACT", "N_HOM"):
        if estimand_by_id.get(estimand_id, {}).get("execution_state") != "NOT_RUN_DESIGN_ONLY":
            errors.append(f"{estimand_id}: design-only state missing")
    required_control_columns = {
        "callable_opportunity", "primary_model_and_null_simulation",
        "polarization_error_control", "mutation_bias_and_demography_control",
        "linkage_recombination_or_phylogeny_control", "multiple_testing",
    }
    if not required_control_columns.issubset(estimand_header):
        errors.append("estimand manifest lacks branch-specific control columns")

    claim_by_id = {row["claim_id"]: row for row in claims}
    if set(claim_by_id) != EXPECTED_CLAIMS:
        errors.append(f"claim IDs differ: {sorted(claim_by_id)}")
    claim_text = "\n".join("\t".join(row.values()) for row in claims)
    for phrase in (
        "without parents/gametes, population frequencies, or outgroups",
        "Population-frequency gBGC was not measured",
        "Non-allelic conversion was not measured",
        "Historical substitution bias is a long-term model-dependent branch signal",
        "population B to event rate",
    ):
        if phrase not in claim_text:
            errors.append(f"claim matrix missing phrase: {phrase}")

    plan = PLAN.read_text(encoding="utf-8")
    for phrase in REQUIRED_PLAN_PHRASES:
        if phrase not in plan:
            errors.append(f"plan missing required phrase: {phrase}")
    if len(plan.splitlines()) < 200:
        errors.append("plan is unexpectedly short")

    return errors


def main() -> int:
    errors = validate_design()
    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        return 1
    print(
        "gene-conversion evidence design valid: 1 direct + 2 phylogenetic "
        "pilots; population/non-allelic design-only"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
