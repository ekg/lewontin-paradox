# Wave mapping reuse audit — P01/P05 pilot mappings (v1)

Schema: `wave-mapping-reuse-audit-v1` · 2026-09-18 · actor: wave-launch agent

## Question

Can the complete pilot mappings under
`/moosefs/erikg/vgp/pilot/runs/vgp10-auth-20260718-v2-pilot-v1/{P01,P05}/mapping/`
be reused as production wave inputs, or must they be remapped?

## (a) Binary identity — MATCH

The first `sweepga.stderr` line of each run records the exact binary invoked:

- Path: `/moosefs/erikg/tier3scratch/sweepga-origin-main-018e4ce/bin-1/sweepga`
  (identical for P01 and P05)
- sha256 of that path, computed 2026-09-18:
  `fa7f0edb9b7e275c288db254046020e136d4267dd5ee043379227ef80da0573b`
- Equals the pinned origin/main build digest recorded in
  `analysis/sweepga_origin_main_build.json` (byte-identical double-build) — MATCH.

## (b) Flag set — MATCH (accepted core, default backend)

Both stderr first lines carry exactly:

`--num-mappings 1:1 --scaffold-jump 0 --overlap 0 --scoring log-length-ani --threads 32`

No `--aligner` flag → sweepga default backend (fastga). Input order `h2.fa h1.fa`
(query = h2, target = h1), consistent with the accepted production orientation.
No deviation from the accepted-core parameterization.

## (c) PAF digests and row counts

| Pair | file | sha256 | rows |
|---|---|---|---:|
| P01 | h2_to_h1.1to1.paf | `eff2b4dddcb8d404256b6fdea566c99dcf9727eff319690828dc412d943341b6` | 12,602 |
| P01 | h2_to_h1.native.1to1.paf | `38c64cb4d13b8af7608d35fda566ae1cdea5d6cc7f3f60086b5629006673ebec` | — |
| P05 | h2_to_h1.1to1.paf | `c91646896257a22f7f54e39e303dd4fc5afcae73958df08c6851f98e12587118` | 38,555 |
| P05 | h2_to_h1.native.1to1.paf | `dbbded38b745eb4c8799c37b1652d2e4f30c32279d9e8f9f6ed3c5ced7241aa6` | — |

## (d) Input binding — MATCH at four levels

1. `pilot/inputs/{P01,P05}/input-manifest.json` sha256 fields:

| Pair | file | manifest sha256 |
|---|---|---|
| P01 | h1.fa | `4a11263679303b8bc398c2906d192bfdde5efb18e10d38b9d142d6d38c6d57f8` |
| P01 | h2.fa | `67c266f5baf39accde67156da0876de47ad586dc999e9437558e015435e9edd0` |
| P05 | h1.fa | `e613177be408823e216a61b568dc4e66ac0e54e50ac3c64a4cf0967800405f68` |
| P05 | h2.fa | `ffe0e8582d655e10882e1d204bfdc374b79249f2e0daf8d9220b281e89547975` |

2. Recomputed sha256 of each file on 2026-09-18 — all four MATCH the manifest.
3. Sequence-dictionary binding: every PAF query contig (P01: 170 unique;
   P05: 160 unique) is in the manifest h2 dictionary, and every target contig
   (P01: 153; P05: 134 unique) is in the manifest h1 dictionary — 0 unbound
   names in all four checks.
4. Transient-evidence note (same standard as accepted P02): the staged scratch
   copies (`/scratch/vgp-map-P01-1782133-AIO9Gc/inputs/…`) were deleted at
   completion, so byte-level identity of the consumed copies rests on the
   staging provenance system shared with the verified P02 accepted run plus
   checks 1–3. `mapping_execution.json` does not restate input digests — the
   gap already documented in `analysis/backend_provenance_audit_v1.md`.

## REUSE_VERDICT

| Pair | verdict |
|---|---|
| P01 | **REUSE** — production mapping source: `pilot/runs/.../P01/mapping/h2_to_h1.1to1.paf` (`eff2b4dd…`) |
| P05 | **REUSE** — production mapping source: `pilot/runs/.../P05/mapping/h2_to_h1.1to1.paf` (`c9164689…`) |

Downstream consumers must reference the `h2_to_h1.1to1.paf` digests above;
the `.native.1to1.paf` files are the pre-filter raws (retained, larger).
