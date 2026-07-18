# VGP ten-pair pilot acquisition handoff

**Manifest version:** `vgp-10-pilot-acquisition-v1.0.0`

**Acquisition UTC:** `2026-07-18T12:18:01Z`

**Execution:** pinned GNU Guix commit `44bbfc24e4bcc48d0e3343cd3d83452721af8c36`; local process only; zero Slurm jobs.

## Outcome

The closed world contains exactly 10 approved primaries and 6 pre-ranked alternates. Exact species, TaxId, catalog row, BioSample, isolate, accession.version, reciprocal linked-assembly roles, VGP BioProject lineage, technologies, and annotation disposition are frozen in `vgp_10_pilot_acquisition_manifest.tsv`. 10/10 primary core acquisition rows have every required core object verified or reused. An alternate was not activated: all alternate biological payload rows remain `superseded` by the retained primary, and no manifest amendment exists.

The direct-control disposition is `selected` / `reused` for `D01`. D01 is bound to 13 complete four-product tetrads and exact Col/Ler/Cvi parents. No trio label or alternate dataset is treated as a complete pedigree, and the 357-GB ENA FASTQ archive is explicitly superseded by the selected immutable event/filter supplements and run metadata.

## Exact accounting

| category | objects | bytes |
|---|---:|---:|
| planned | 182 | 375689532447 |
| transferred (I/O flow; orthogonal) | 0 | 0 |
| verified by local SHA-256 (validation state; orthogonal) | 149 | 12778465073 |
| newly promoted | 0 | 0 |
| reused | 149 | 12778465073 |
| missing | 20 | 0 |
| superseded/not activated | 13 | 362911067374 |
| quarantined | 0 | 0 |

Newly promoted/reused/missing/superseded/quarantined are mutually exclusive and reconcile to planned objects and bytes: `True`. Transferred bytes are an orthogonal physical-I/O measure, and verified bytes are an orthogonal validation-state measure; neither is added to terminal-disposition bytes.

## Fail-closed branch status

Same-individual or H1/H2 ambiguity, accession drift, size mismatch, provider-MD5 mismatch, and local-SHA mismatch are hard core refusals. Annotation absence, a paired-RefSeq dictionary difference, or annotation download failure is branch-local and cannot veto a valid core. Published H1 annotation BUSCO values are retained where exact; missing H2 BUSCO and missing exact-final-sequence Merqury QV are explicitly missing, never imputed. Hi-C absence alone is not used as a refusal.

No raw-read or k-mer payload is in the object inventory. The resolver never lists SRA runs or GenomeArk raw-data prefixes. This proves that zero unmanifested bulk raw-read objects/bytes were acquired; the choice deliberately leaves selective validation/QV evidence incomplete rather than expanding scope silently.

## Storage protocol demonstrated

`summary.json` records a same-code-path demonstration of partial resume, size/provider digest checks, pre/post-promotion SHA-256, atomic `os.replace`, read-only promotion, CAS reuse with local SHA-256 revalidation, and mismatch quarantine. The demonstration accounting is one verified/reused 16-byte object plus one quarantined 16-byte candidate and is deliberately separate from biological-object accounting. Real objects use the same `objects/sha256/<2>/<2>/<sha256>` contract and are directly reusable by the official Freeze 1 mirror.
