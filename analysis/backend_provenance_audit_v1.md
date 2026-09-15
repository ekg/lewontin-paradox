# Backend provenance audit — three-pair bounded-production PAFs

Schema: `vgp-backend-provenance-audit-v1` · date: 2026-09-07 (UTC)

Question: which aligner backend produced each bounded-production mapping PAF
(P02/P03/P07), from recorded evidence only.

## Verdicts

| Pair | production PAF sha256 (bound in `analysis/vgp_three_pair_execution_v2.json` → `pairs[].coordinate_and_strand_audit.paf_sha256`) | backend | verdict |
|---|---|---|---|
| P02 | `17c6cefbf39f35fac330f5917e9f077469f0f927fd715335d4ab31c0785ab61c` | **fastga** | determined |
| P03 | `1c062bc4fa566c18db904984920c14efc8169649aa88754ab2c4b273df18860a` | **wfmash** (documented fallback) | determined |
| P07 | `0bcf45366bca9217cd012242cedf79921a6cafd6b1fb1f5247db7e54512d00fa` | **fastga** | determined |

**The accepted three-pair core is backend-mixed: P02 fastga, P03 wfmash, P07 fastga.**

## Evidence chains

### P02 → fastga

1. `analysis/vgp_three_pair_execution_v2.json` → `pairs[P02].coordinate_and_strand_audit.paf_sha256`
   = `17c6cefb…ab61c` (55,907 rows, H2_query_to_H1_reference) — the digest the
   audited production consumed.
2. That digest equals `h2_to_h1.1to1.paf` in
   `/moosefs/erikg/vgp/pilot/runs/vgp10-auth-20260718-v2-pilot-v1/P02/mapping/.complete.json`
   (= `17c6cefbf39f35fac330f5917e9f077469f0f927fd715335d4ab31c0785ab61c`).
3. Same directory, `sweepga.stderr` (job 1782140):
   `[INFO ] [sweepga::align 5.6s] Running fastga alignment: …/h2.fa -> …/h1.fa`
   followed by `[FastGA] prepare_gdb: Converting …` and
   `[FastGA] Calling FAtoGDB: /moosefs/erikg/tier3scratch/sweepga-origin-main-018e4ce/bin-1/FAtoGDB …`.
4. `mapping_execution.json` records `fastga_amendment` (authorization
   `vgp10-auth-20260718-v2`, `required_option = --num-mappings 1:1`).

### P03 → wfmash

1. `analysis/vgp_three_pair_execution_v2.json` → `pairs[P03].coordinate_and_strand_audit.paf_sha256`
   = `1c062bc4…18860a` (394 rows).
2. That digest equals `h2_to_h1.1to1.paf` in
   `/moosefs/erikg/vgp/pilot/runs/vgp-three-pair-20260722-v1/P03/mapping/.complete.json`
   and `multiplicity.json` → `paf_sha256`.
3. Same directory, `wfmash_fallback_contract.json`
   (schema `vgp-three-pair-wfmash-fallback-v1`):
   `"fallback_trigger": "reproducible corrected FastGA failure in frozen selection evidence"`,
   `"same_staged_fasta_bytes_as_corrected_retry": true`, with staged H1/H2 sha256s
   `92697b85…` / `8270fea2…`.
4. Consistent with `analysis/vgp_three_pair_execution_v2.json`
   `pairs[P03].failure_class = "prior_fastga_execution_failure"`.

### P07 → fastga

1. `analysis/vgp_three_pair_execution_v2.json` → `pairs[P07].coordinate_and_strand_audit.paf_sha256`
   = `0bcf4536…2d00fa` (4,297 rows).
2. That digest equals `mapping.paf_sha256` in
   `/moosefs/erikg/vgp/pilot/clean-canary/vgp-clean-canary-20260722-v1/P07/execution.json`,
   whose `mapping.fastga_scratch_contract` records:
   `observed_cwds: ["/mnt/sdb1/scratch/vgp-map-P07-1791510-sXWE6K/inputs"]`,
   `fastga_pids: [3857]`, `fastga_snapshot_count: 889`, and 40+ observed managed
   open paths `…/_algn.3857.N.las` (FastGA native LAS intermediates), with
   `contract_valid: true`.
3. The P07 bounded production reused the clean-canary PAF/IMPG index (recorded
   in the three-pair evidence; the bounded `execution.json` itself does not
   restate the input digest — a provenance gap noted below).

## Other lineages (context)

- Tier 3A corrected remaps (all three tuples): `--aligner wfmash` explicitly —
  `/moosefs/erikg/tier3data/tier3a-origin-remap-20260716/spinachia_spinachia_SK-2024b/mapping/command.txt`
  (`sweepga … --aligner wfmash --num-mappings 1:1 --scaffold-jump 0 …`).
- Wave P09 mapping (job 2833615): same remap template, wfmash backend.
- The P07 paralog-triage input PAF (`bfd2b8c3…`, tier3a remap) is therefore
  wfmash-backend while the bounded P07 core PAF is fastga-backend.

## Consequences and recommendations

1. Cross-leg π comparability (bounded core vs tier3a coding panels vs wave)
   currently mixes backends; review defect A-3 measured 0.25 variant Jaccard
   between backends on P03 chr1 despite 0.97 coverage Jaccard.
2. Standardizing on **wfmash** (majority lineage: tier3a remaps, P03, wave;
   byte-identical reproducible Guix builds; proven at scale) would require
   re-mapping + re-running bounded production for P02 and P07 (small genomes).
   Standardizing on fastga would require remapping the tier3a panels, P03, and
   the wave, with FastGA robustness unproven at 6 Gbp scale.
3. The P09 both-backend calibration (job 2833615 wfmash + fastga companion)
   quantifies A-3 on the largest genome before the choice is finalized.
4. Provenance gap to close regardless of choice: bounded `execution.json`
   should record the input PAF path+digest (it is currently only bound via the
   three-pair execution evidence and `.complete.json` of the producing run).
