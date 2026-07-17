# VGP pilot independent Ne inventory

Date: 2026-07-17 UTC

## Scope

This inventory is keyed exactly to `analysis/vgp_pilot_manifest.tsv` and the
upstream freeze provenance in `analysis/vgp_phase1_freeze_provenance.json`.
On 2026-07-17 UTC the frozen manifest contained `pilot_selected = no` for every
row, and the provenance summary reported `selected_count = 0`.

Because no pilot species were selected, this task did not retrieve
species-specific literature or repository metadata. Instead it emits a
fail-closed empty inventory whose schema is ready for downstream joins and
whose validator will stop if a later manifest selects one or more species
without corresponding manual curation.

## Manifest status

| field | value |
| --- | --- |
| selected rows in manifest | 0 |
| selected rows in provenance summary | 0 |
| freeze provenance generated at | 2026-07-17T18:24:11Z |
| species-specific Ne/life-history records retrieved | 0 |
| population-data availability rows emitted | 0 |

## Missing data disposition

- No selected pilot species were available for inventory.
- `analysis/vgp_pilot_ne_sources.tsv` is intentionally header-only because the
  denominator of selected taxa is empty.
- `analysis/vgp_pilot_population_data_availability.tsv` is intentionally
  header-only for the same reason.
- If a future manifest changes `pilot_selected` away from all-`no`, rerun this
  task with manual source curation; the builder and validator fail closed in
  that state.

## Circularity and authorization notes

- No Ne estimate was algebraically derived from `pi/(4mu)` or reused as an
  independent predictor.
- No raw reads, assemblies, VCFs, or other population payloads were downloaded.
- No demographic inference, PSMC, MSMC2, or SMC++ execution was attempted.
- No claim is made that VGP H1/H2 assembly pairs are callable diploid genotypes.

## Outputs

- `analysis/vgp_pilot_ne_sources.tsv`
- `analysis/vgp_pilot_population_data_availability.tsv`
- `analysis/vgp_pilot_ne_inventory.md`
