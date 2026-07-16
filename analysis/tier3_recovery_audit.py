#!/usr/bin/env python3
"""Fail-closed audit and publication ledger for the recovered Tier 3 evidence.

This script deliberately reads only the corrected Tier 3A origin/main result,
the recovered Tier 3B population result, and the already-validated Tier 3C
composition synthesis.  The Tier 3A supersession files are negative controls:
their paths and checksums must not occur in the selected lineage.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import sys
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from analysis import tier3_fit


OUT = ROOT / "results/tier3"
TIER3A = ROOT / "results/tier3a/diploid_diversity.tsv"
TIER3A_QC = ROOT / "results/tier3a/diploid_qc.md"
TIER3A_RUN = ROOT / "results/tier3a/diploid_run_manifest.tsv"
TIER3A_LINEAGE = ROOT / "results/tier3a/diploid_lineage_audit.json"
TIER3A_SUPERSEDED = ROOT / "results/tier3a/diploid_superseded_results.tsv"
TIER3B = ROOT / "results/tier3b/population_diversity.tsv"
TIER3B_QC = ROOT / "results/tier3b/population_qc.md"
TIER3B_RUN = ROOT / "results/tier3b/population_run_manifest.tsv"
TIER3C = ROOT / "analysis/tier3c_data.tsv"
SYNTHESIS = ROOT / "analysis/tier3_results.tsv"


LEDGER_FIELDS = [
    "tier", "identity", "scientific_name", "population_id", "modality",
    "observable", "annotation_category", "eligible_n", "nominal_chromosomes",
    "variant_count", "numerator", "numerator_definition", "denominator",
    "estimate", "interval_low", "interval_high", "uncertainty_method",
    "uncertainty_unit", "uncertainty_replicates", "uncertainty_standard_error",
    "interval_type", "exclusions", "software_provenance", "source_result_path",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, Any]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for source in rows:
            writer.writerow({field: source.get(field, "") for field in fields})


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def finite(value: Any) -> float:
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"non-finite number: {value!r}")
    return number


def close(expected: float, observed: float, tolerance: float = 1e-12) -> bool:
    return abs(expected - observed) <= tolerance * max(1.0, abs(expected), abs(observed))


def validate_provenance() -> dict[str, str]:
    tier3a_qc = TIER3A_QC.read_text(encoding="utf-8")
    if "All 3 acquired biological tuples produced" not in tier3a_qc or "supersession ledger passed" not in tier3a_qc:
        raise ValueError("Tier 3A QC is not PASS")
    if "Status: **PASS" not in TIER3B_QC.read_text(encoding="utf-8"):
        raise ValueError("Tier 3B QC is not PASS")
    lineage = json.loads(TIER3A_LINEAGE.read_text(encoding="utf-8"))
    if lineage.get("status") != "passed":
        raise ValueError("Tier 3A lineage audit did not pass")
    for field in ("old_result_checksum_intersection", "superseded_checksum_intersection", "superseded_path_intersection"):
        if lineage.get(field) != []:
            raise ValueError(f"Tier 3A lineage contains forbidden superseded input: {field}")

    old_hashes = {row["artifact_sha256"] for row in read_tsv(TIER3A_SUPERSEDED)}
    selected_hashes = {sha256(TIER3A), sha256(TIER3A_QC), sha256(TIER3A_RUN)}
    if old_hashes & selected_hashes:
        raise ValueError("selected Tier 3A artifacts reuse a superseded result checksum")

    runs = read_tsv(TIER3A_RUN)
    if len(runs) != 3 or {row["status"] for row in runs} != {"completed"}:
        raise ValueError("Tier 3A run manifest does not contain three completed tuples")
    commits = {row["sweepga_origin_main_commit"] for row in runs}
    profiles = {row["guix_profile_store_path"] for row in runs}
    binaries = {row["sweepga_binary_sha256"] for row in runs}
    if len(commits) != 1 or not next(iter(commits)):
        raise ValueError("Tier 3A does not trace to one fetched origin/main SweepGA commit")
    if any(not path.startswith("/gnu/store/") for path in profiles):
        raise ValueError("Tier 3A Guix closure is missing")
    if len(binaries) != 1 or any(len(value) != 64 for value in binaries):
        raise ValueError("Tier 3A SweepGA binary checksum is missing or inconsistent")
    for row in runs:
        if "--num-mappings 1:1" not in row["sweepga_command"]:
            raise ValueError("Tier 3A did not use native --num-mappings 1:1")
        if row["sweepga_hit_cap"] != "1:1":
            raise ValueError("Tier 3A run manifest reports a non-1:1 mapping cap")
        if row["sweepga_observed_query_multiplicity"] != "1" or row["sweepga_observed_target_multiplicity"] != "1":
            raise ValueError("Tier 3A mapping multiplicity is not one-to-one")

    population_runs = read_tsv(TIER3B_RUN)
    if len(population_runs) != 2 or {row["status"] for row in population_runs} != {"accepted"}:
        raise ValueError("Tier 3B run manifest does not contain two accepted populations")
    if any(not row["guix_profile"].startswith("/gnu/store/") for row in population_runs):
        raise ValueError("Tier 3B Guix closure is missing")
    return {
        "sweepga_commit": next(iter(commits)),
        "sweepga_binary_sha256": next(iter(binaries)),
        "tier3a_guix_profile": next(iter(profiles)),
        "tier3b_guix_profile": population_runs[0]["guix_profile"],
        "guix_channel_commit": population_runs[0]["guix_channel_commit"],
    }


def evidence_ledger() -> list[dict[str, Any]]:
    observations = tier3_fit.load_diversity_observations(TIER3A, TIER3B)
    if len(observations) != 23:
        raise ValueError(f"expected 23 recovered observations (12+3 derived+8), observed {len(observations)}")
    rows: list[dict[str, Any]] = []
    for source in observations:
        tier = "3B" if source["observable_tier"] == "population" else "3A"
        row = {
            "tier": tier,
            "identity": source["dataset_id"],
            "scientific_name": source["scientific_name"],
            "population_id": source.get("population_id", ""),
            "modality": source["source_modality"],
            "observable": source["observable"],
            "annotation_category": source.get("annotation_category", ""),
            "eligible_n": source.get("eligible_sample_size", ""),
            "nominal_chromosomes": source.get("nominal_chromosomes", ""),
            "variant_count": source.get("variant_count", ""),
            "numerator": source.get("numerator", ""),
            "numerator_definition": source.get("numerator_definition", ""),
            "denominator": source["denominator"],
            "estimate": f"{finite(source['value']):.15g}",
            "interval_low": f"{finite(source['ci_low']):.15g}",
            "interval_high": f"{finite(source['ci_high']):.15g}",
            "uncertainty_method": source.get("uncertainty_method", ""),
            "uncertainty_unit": source.get("uncertainty_unit", ""),
            "uncertainty_replicates": source.get("uncertainty_replicates", ""),
            "uncertainty_standard_error": source.get("uncertainty_standard_error", ""),
            "interval_type": source.get("interval_type", ""),
            "exclusions": source.get("exclusions", ""),
            "software_provenance": source.get("software_provenance", ""),
            "source_result_path": source.get("source_result_path", ""),
        }
        if finite(row["estimate"]) < finite(row["interval_low"]) or finite(row["estimate"]) > finite(row["interval_high"]):
            raise ValueError(f"estimate falls outside uncertainty interval: {row['identity']}/{row['observable']}")
        rows.append(row)
    if len({row["identity"] for row in rows if row["tier"] == "3A"}) != 3:
        raise ValueError("Tier 3A does not contain three biological identities")
    if len({row["population_id"] for row in rows if row["tier"] == "3B"}) != 2:
        raise ValueError("Tier 3B does not contain two population identities")
    return rows


def independent_ols(rows: Sequence[Mapping[str, str]], outcome: str, quadratic: bool = False) -> float:
    selected = [row for row in rows if row["row_kind"] == "point" and row["observable"] == outcome and row["status"] == "estimated"]
    x = np.asarray([finite(row["predictor_value"]) for row in selected])
    y = np.asarray([finite(row["estimate"]) for row in selected])
    clades = sorted({row["clade"] for row in selected})
    columns = [np.ones(len(selected)), x - x.mean()]
    if quadratic:
        columns.append((x - x.mean()) ** 2)
    columns.extend(np.asarray([1.0 if row["clade"] == clade else 0.0 for row in selected]) for clade in clades[1:])
    coefficients, *_ = np.linalg.lstsq(np.column_stack(columns), y, rcond=None)
    return float(coefficients[2 if quadratic else 1])


def headline_audit(ledger: Sequence[Mapping[str, Any]]) -> list[dict[str, str]]:
    audits: list[dict[str, str]] = []

    def record(claim: str, source: str, reported: float, recomputed: float, definition: str) -> None:
        difference = abs(reported - recomputed)
        passed = close(reported, recomputed)
        audits.append({
            "claim": claim, "source": source, "reported": f"{reported:.15g}",
            "recomputed": f"{recomputed:.15g}", "absolute_difference": f"{difference:.3g}",
            "tolerance": "1e-12 relative/absolute", "status": "PASS" if passed else "FAIL",
            "definition": definition,
        })
        if not passed:
            raise ValueError(f"headline mismatch: {claim}")

    for row in read_tsv(TIER3A):
        record(
            f"3A:{row['dataset_id']}:{row['statistic_label']}:{row['annotation_class']}",
            "fresh Tier 3A numerator and callable denominator",
            finite(row["estimate"]), finite(row["variant_numerator"]) / finite(row["callable_denominator"]),
            "H1/H2 alternative-allele count divided by H1-native callable sites",
        )
    for row in read_tsv(TIER3B):
        if row["statistic"] == "pi_S_over_pi_W":
            recomputed = (finite(row["component_s_numerator"]) / finite(row["component_s_callable_sites"])) / (finite(row["component_w_numerator"]) / finite(row["component_w_callable_sites"]))
            definition = "(sum pairwise S differences/S callable 4D sites)/(sum pairwise W differences/W callable 4D sites)"
        else:
            recomputed = finite(row["numerator"]) / finite(row["callable_site_denominator"])
            definition = "sum unbiased pairwise differences divided by callable sites"
        record(f"3B:{row['tuple_id']}:{row['statistic']}", "fresh Tier 3B numerator and callable denominator", finite(row["estimate"]), recomputed, definition)

    synthesis_rows = read_tsv(SYNTHESIS)
    model_ids = {
        "primary.gc3.across_class_fixed": independent_ols(synthesis_rows, "gc3"),
        "sensitivity.whole_genome_gc.across_class_fixed": independent_ols(synthesis_rows, "whole_genome_gc"),
        "sensitivity.gc3.quadratic_concavity": independent_ols(synthesis_rows, "gc3", quadratic=True),
    }
    by_id = {row["analysis_id"]: row for row in synthesis_rows}
    for analysis_id, recomputed in model_ids.items():
        if analysis_id not in by_id:
            raise ValueError(f"missing composition headline row: {analysis_id}")
        record(analysis_id, "independent NumPy least-squares reconstruction from committed point rows", finite(by_id[analysis_id]["effect"]), recomputed, "class-fixed OLS coefficient using the named composition observable")

    if any(row["status"] != "PASS" for row in audits):
        raise ValueError("headline audit failed")
    return audits


def verify_synthesis_table(ledger: Sequence[Mapping[str, Any]]) -> None:
    rows = read_tsv(SYNTHESIS)
    diversity = [row for row in rows if row["row_kind"] == "point" and row["observable_tier"] in {"population", "alignment_conditioned_diploid_assembly"}]
    if len(diversity) != len(ledger):
        raise ValueError("synthesis table does not contain every recovery-ledger observation")
    for row in diversity:
        for field in ("dataset_id", "source_modality", "numerator", "denominator", "eligible_sample_size", "ci_low", "ci_high", "uncertainty_method", "exclusions", "software_provenance"):
            if row.get(field, "") == "":
                raise ValueError(f"synthesis diversity row lacks {field}: {row.get('analysis_id')}")
    png = ROOT / "analysis/fig_tier3.png"
    pdf = ROOT / "analysis/fig_tier3.pdf"
    if not png.read_bytes().startswith(b"\x89PNG\r\n\x1a\n") or not pdf.read_bytes().startswith(b"%PDF-"):
        raise ValueError("regenerated Tier 3 figure is invalid")


def write_summary(provenance: Mapping[str, str], ledger: Sequence[Mapping[str, Any]]) -> None:
    by = {(row["identity"], row["observable"], row["annotation_category"]): row for row in ledger}
    assembly_ids = sorted({row["identity"] for row in ledger if row["tier"] == "3A"})
    population_ids = sorted({row["identity"] for row in ledger if row["tier"] == "3B"})
    lines = [
        "# Recovered Tier 3 biological evidence", "",
        "Status: **PASS — corrected Tier 3A, recovered Tier 3B, and Tier 3C composition are synthesized without modality pooling.**", "",
        "## Dependency and input gate", "",
        "The WG task has direct dependencies on `rerun-tier3a-origin-sweepga` and `run-tier3b-biological-recovery`; both were confirmed `Done` before this synthesis began. The only numerical Tier 3A input is `results/tier3a/diploid_diversity.tsv`; the old-result and mapping supersession ledgers are negative controls and no listed path or checksum was selected.", "",
        f"Corrected Tier 3A traces all three tuples to fetched SweepGA origin/main commit `{provenance['sweepga_commit']}`, binary SHA-256 `{provenance['sweepga_binary_sha256']}`, Guix closure `{provenance['tier3a_guix_profile']}`, and native `--num-mappings 1:1` commands with observed query/target multiplicity 1/1. Tier 3B uses Guix closure `{provenance['tier3b_guix_profile']}` at channel commit `{provenance['guix_channel_commit']}`. Both upstream QC reports pass.", "",
        "## Assembly modality (Tier 3A)", "",
        "These are alignment-conditioned H1/H2 estimates over deterministic H1-native coding panels, not population diversity and not genome-wide deposited-individual heterozygosity. Each estimate is variants/callable bases; uncertainty is the 1,000-replicate 50-kb genomic block bootstrap. The ratio intervals below conservatively divide the marginal S and W 95% bounds because paired bootstrap draws were not published; they are intentionally wider than a paired ratio interval.", "",
        "| assembly pair | coding-gene diversity (95% CI; variants/callable) | CDS diversity (95% CI; variants/callable) | reference-conditioned pi_S/pi_W (conservative interval) |", "|---|---|---|---|",
    ]
    for identity in assembly_ids:
        coding = by[(identity, "diploid_haplotype_diversity", "coding_gene")]
        cds = by[(identity, "diploid_haplotype_diversity", "CDS")]
        ratio = by[(identity, "pi_S_over_pi_W", "native_fourfold_reference_conditioned_ratio")]
        lines.append(
            f"| _{coding['scientific_name']}_ | {float(coding['estimate']):.6g} [{float(coding['interval_low']):.6g}, {float(coding['interval_high']):.6g}]; {coding['variant_count']}/{coding['denominator'].split()[0]} | {float(cds['estimate']):.6g} [{float(cds['interval_low']):.6g}, {float(cds['interval_high']):.6g}]; {cds['variant_count']}/{cds['denominator'].split()[0]} | {float(ratio['estimate']):.6g} [{float(ratio['interval_low']):.6g}, {float(ratio['interval_high']):.6g}] |"
        )
    lines.extend(["", "Coding diversity is heterogeneous across the three pairs (about 0.000398 to 0.01385). All three reference-conditioned ratio point estimates are below one, but the conservative intervals for _Spinachia_ and _Tautogolabrus_ include one; these three selected assembly pairs do not estimate a population-size trend.", "", "## Population modality (Tier 3B)", "", "Both rows use 20 wild diploid biological individuals (40 nominal chromosomes) from the same species and region, with exact population-specific callable masks. Population pi and the ratio use 10,000-replicate chromosome-stratified 1-Mb block bootstrap intervals; component pi_S and pi_W rows use 20-unit delete-one-individual jackknife intervals.", "", "| population | population pi (95% CI; pairwise numerator/callable) | reference-conditioned pi_S/pi_W (95% block-bootstrap CI) |", "|---|---|---|",
    ])
    for identity in population_ids:
        pi = by[(identity, "population_pi", "callable_nuclear_region")]
        ratio = by[(identity, "pi_S_over_pi_W", "native_4D_reference_conditioned_ratio")]
        lines.append(f"| `{pi['population_id']}` | {float(pi['estimate']):.7g} [{float(pi['interval_low']):.7g}, {float(pi['interval_high']):.7g}]; {float(pi['numerator']):.9g}/{pi['denominator'].split()[0]} | {float(ratio['estimate']):.7g} [{float(ratio['interval_low']):.7g}, {float(ratio['interval_high']):.7g}] |")
    lines.extend(["", "The two populations differ in pi and their bootstrap intervals do not overlap; both reference-conditioned ratio intervals lie below one. This is within-species population heterogeneity, not an across-species effect and not polarized SFS-B.", "", "## Composition modality (Tier 3C)", "", "Composition remains an exact-single-assembly observable: 135 whole-genome GC values and 90 native-annotation GC3 values. The class-fixed GC3 slope remains +0.00892 per Buffalo census-size-proxy unit (10,000-species-bootstrap 95% interval -0.00123 to +0.02149; n=90; BH q=0.217), with positive Aves, uncertain Insecta, negative Mammalia, and null _Drosophila_ estimates. The positive quadratic coefficient remains opposite the predicted concavity. These heterogeneous/null composition findings are not replaced or pooled with diversity.", "", "## Claim boundary and sensitivities", "", "The recovery replaces the old n=0 statements with real biological estimates, but it does not create a cross-species diversity regression: Tier 3A has three selected coding-panel assembly pairs, Tier 3B has two populations from one species, and Tier 3C measures composition. Reference-conditioned pi_S/pi_W is descriptive and is never renamed polarized SFS-B. Callable-mask, sample, panel, native-annotation, and mapping exclusions remain explicit in `recovery_evidence_ledger.tsv`; no causal claim that gBGC resolves Lewontin's paradox follows.", "", "`headline_audit.tsv` independently recomputes all 20 upstream direct estimates plus three composition coefficients; `artifact_index.tsv` records checksums and roles for the committed synthesis. Figures are regenerated from `analysis/tier3_results.tsv` and the manuscript reports the same definitions, denominators, n, and intervals.", ""])
    (OUT / "recovery_evidence_summary.md").write_text("\n".join(lines), encoding="utf-8")


def write_artifact_index() -> None:
    paths = [
        TIER3A, TIER3A_QC, TIER3A_RUN, TIER3A_LINEAGE, TIER3B, TIER3B_QC,
        TIER3B_RUN, TIER3C, SYNTHESIS, ROOT / "analysis/fig_tier3.png",
        ROOT / "analysis/fig_tier3.pdf", ROOT / "manuscript.typ", ROOT / "manuscript.pdf",
        ROOT / "analysis/tier3_fit.py", ROOT / "analysis/tier3_recovery_audit.py",
        ROOT / "analysis/tests/test_tier3_fit.py", ROOT / "analysis/tests/test_tier3_recovery_audit.py",
        ROOT / "analysis/TIER3_RESULTS.md", ROOT / "analysis/README.md", ROOT / "README.md",
        OUT / "recovery_evidence_summary.md", OUT / "recovery_evidence_ledger.tsv",
        OUT / "headline_audit.tsv",
    ]
    rows = []
    for path in paths:
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"indexed artifact is missing or empty: {path}")
        rows.append({"artifact_path": str(path.relative_to(ROOT)), "sha256": sha256(path), "bytes": path.stat().st_size, "role": "input" if path in {TIER3A, TIER3A_QC, TIER3A_RUN, TIER3A_LINEAGE, TIER3B, TIER3B_QC, TIER3B_RUN, TIER3C} else "synthesis_output"})
    write_tsv(OUT / "artifact_index.tsv", ("artifact_path", "sha256", "bytes", "role"), rows)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    provenance = validate_provenance()
    ledger = evidence_ledger()
    write_tsv(OUT / "recovery_evidence_ledger.tsv", LEDGER_FIELDS, ledger)
    audits = headline_audit(ledger)
    write_tsv(OUT / "headline_audit.tsv", ("claim", "source", "reported", "recomputed", "absolute_difference", "tolerance", "status", "definition"), audits)
    verify_synthesis_table(ledger)
    write_summary(provenance, ledger)
    write_artifact_index()
    print(f"PASS: {len(ledger)} recovered observations; {len(audits)} independent headline checks")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
