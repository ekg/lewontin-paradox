# Ten-pair VGP pilot run result

Run ID: `vgp10-20260718-preflight-v1`
Terminal state: **10/10 primary slots failed mandatory preflight; 0 biological jobs submitted**

## Outcome

The closed-world acquisition and pinned GNU Guix capture were revalidated successfully. Exact pair, BioSample/individual, accession.version, accepted core-object size, and SHA-256 identity are resolved for all ten immutable primaries. Those successes do not authorize biological compute: every primary lacks exact-final-sequence QV, H2 BUSCO, a manifest-bound k-mer/copy-number collapse audit, and a repeat/low-complexity mask. P07 and P08 also lack H1 BUSCO; P09 retains unresolved exact read chemistry.

The approved preflight validator requires QV >= 40, BUSCO completeness/missing/duplication measurements for both haplotypes, and a passing copy-number/k-mer audit. Missing values were not imputed from technology labels, H1-only annotation BUSCO, assembly length, or published plots. Consequently SweepGA, multiplicity auditing, IMPG, VCF/BCF, masks, consensus, diversity, PSMC, 200 bootstraps, scaling scenarios, and validation comparisons are explicitly `not_run_preflight_failed` for every pair. This is a failed pilot execution, not a zero-diversity result.

## Primary-slot accounting

| state | count |
|---|---:|
| completed core analysis | 0 |
| explicitly failed at preflight | 10 |
| silently dropped or replaced | 0 |
| alternate activations | 0 |
| Slurm jobs / dependency edges | 0 / 0 |

No alternate was activated: there is no versioned amendment authorizing one, and the six declared alternates remain `standby_not_triggered`. Primary failure rows are retained in the result manifest.

## Machine-readable failure reasons

- `MISSING_EXACT_FINAL_SEQUENCE_QV`: 10/10
- `MISSING_H1_BUSCO`: 2/10
- `MISSING_H2_BUSCO`: 10/10
- `MISSING_KMER_COPY_NUMBER_AUDIT`: 10/10
- `UNRESOLVED_EXACT_READ_CHEMISTRY`: 1/10
- `MISSING_REPEAT_OR_LOW_COMPLEXITY_MASK`: 10/10

Annotation status remains branch-local. P03 and P04 retain their acquisition-time annotation mismatch; all other annotation branches were not reached because core preflight failed. No annotation disposition was used to create a core failure.

## Validation interpretation

There are zero unresolved pair/accession identities and zero checksum failures. There are also zero retained mappings with multiplicity greater than one because no mapping was retained; overlap depth, H1 REF, direct H2 reconstruction, callable-mask totals, VCF/consensus/PSMC consistency, bootstrap success fraction, and scaled trajectories are **not measured**, not passing zeros. D01 and raw-read/k-mer/published comparisons were not run because no primary produced a core result to compare. The per-pair artifacts preserve these distinctions explicitly.

## Reproducibility and telemetry

Execution used captured Guix channel `44bbfc24e4bcc48d0e3343cd3d83452721af8c36`, manifest `afd14da446d29d4e66b64837194dfae524342dc5b8de3c2778ec85c1cca09a3d`, profile `/gnu/store/3c2mxm30rbzvnw7qsi235mrkk3m38fym-profile`, derivation `/gnu/store/g05w4h08h891dcaf7v4gy3irf2zwbwa7-profile.drv`, and closure `8fcdb32021f1cd8eac839509cff47ab6bdd63b656b30e243fdf78d3c4ba24f9d`. Live CAS revalidation read 12567760437 logical bytes across exactly 100 primary core objects. No node-local scratch was allocated because the submission boundary was not crossed. Per-pair and aggregate elapsed time, CPU, peak RSS, filesystem/logical reads and writes, metadata operations, retries, scheduler IDs, and dispositions are recorded in `analysis/vgp_10_pilot_resource_telemetry.tsv`.

## Delivered artifacts

- `analysis/vgp_10_pilot_result_manifest.tsv`: exactly ten retained primary terminal rows and per-pair artifact paths.
- `analysis/vgp_10_pilot_qc.tsv`: gate-by-gate measured/missing/not-reached distinctions.
- `analysis/vgp_10_pilot_resource_telemetry.tsv`: ten preflight revalidation rows plus aggregate telemetry.
- `analysis/vgp_10_pilot_pair_artifacts/P01..P10/`: QC, diversity, PSMC trajectory, bootstrap, scenario, validation, annotation, telemetry, and failure artifacts.
- `analysis/vgp_10_pilot_run_summary.json` and `analysis/vgp_10_pilot_command_log.json`: independently recomputable counts, identities, command dispositions, and zero-job ledger.
