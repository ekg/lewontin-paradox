# VGP pilot raw-read validation handoff v1

Date: 2026-07-18  
Task: `acquire-vgp-validation-reads`  
Canonical shared root: `/moosefs/erikg/vgp`  
Machine-readable plan: `analysis/vgp_validation_read_plan_v1.json`  
Machine-readable outcome: `analysis/vgp_validation_reads_manifest_v1.json`

This acquisition is independent of and does not gate the ten-pair
assembly-derived run. All biological payload, staging, CAS, accession-view,
quarantine, and downstream-manifest paths resolve from
`analysis/vgp_data_root_config.json`. The former
`/moosefs/erikg/lewontin-paradox-data/vgp/phase1-freeze-1.0` root was inventoried
as migration input only; it contained zero matching validation-read objects or
bytes, so nothing was redownloaded in place or newly written there.

## Stratified selection

| pilot | stratum | exact individual | exact diploid assembly pair | raw-read scope |
|---|---|---|---|---|
| P07, *Spinachia spinachia* | small (813,562,223 bp across both haplotypes), later HiFi/Hi-C, low expected diversity | `SAMN36735485`, `fSpiSpi1` | `GCA_048126635.1` / `GCA_048127205.1` | complete PacBio run plus both useful Illumina WGS mates planned |
| P09, *Heterodontus francisci* | very large/repeat-sensitive (11,209,934,419 bp across both haplotypes), later HiFi/Hi-C | `SAMN39432692`, `sHetFra1` | `GCA_036365525.1` / `GCA_036365495.1` | one complete HiFi validation run from the seven-run library acquired; not represented as the full library |
| P04, *Falco naumanni* | early CLR/TrioCanu parental generation | `SAMN16870685`, `bFalNau1` | `GCA_017639655.1` / `GCA_017639645.1` | complete public primary CLR run planned |

Every object row repeats its pilot ID, species/TaxId, BioSample,
individual/isolate, exact H1/H2 accession.version, run, experiment, BioProject,
platform, instrument, library fields, read/base counts, scope, authoritative ENA
URL, expected bytes, and ENA MD5. NCBI Datasets exact assembly reports and the
frozen authorization roster bind both assemblies to the BioSample/individual;
the ENA and NCBI SRA exact-BioSample records bind each library/run to the same
individual.

## Verified real payload

`SRR29944135` / `SRX25437651` is PacBio Sequel II genomic WGS from library
`sHetFra1_PacBio_HiFi_fastq_4`. The submitter's exact library name resolves the
generic “PacBio Sequel” label in the P09 assembly reports to HiFi for this run.
It is a complete archived run/SMRT-cell validation subset, not all seven HiFi
runs used by the assembly.

| check | observed |
|---|---:|
| ENA object | `SRR29944135_subreads.fastq.gz` |
| expected and observed bytes | 8,082,600,971 |
| ENA MD5 | `835004b084811fa992fa79bfaa0de5f0` (match) |
| local SHA-256 | `0f1e4794dce9552852dc28042c2c31e58f4a11f07cc2043ed3825cf8facada3f` |
| gzip magic | pass |
| post-promotion full rehash | pass |
| CAS mode | `0440` |
| CAS path | `/moosefs/erikg/vgp/objects/sha256/0f/1e/0f1e4794dce9552852dc28042c2c31e58f4a11f07cc2043ed3825cf8facada3f` |
| accession view | `/moosefs/erikg/vgp/views/accession/P09/SRR29944135/SRR29944135_subreads.fastq.gz` |

The interrupted single-stream transfer was resumed by four bounded HTTPS
ranges using the cluster CA bundle, without disabling certificate checks. A
host-Python compatibility error happened only after atomic CAS promotion and
before view publication. The next invocation discovered the exact-size
unmanifested CAS object, revalidated ENA MD5 and local SHA-256, created the
view, and published the ledger without redownload. This recovery is covered by
`test_recover_unmanifested_promotion_without_redownload`.

## Complete object/byte ledger

| disposition | objects | bytes |
|---|---:|---:|
| planned | 5 | 73,402,884,306 |
| transferred | 1 | 8,082,600,971 |
| verified | 1 | 8,082,600,971 |
| newly promoted | 1 | 8,082,600,971 |
| reused | 0 | 0 |
| missing | 0 | 0 |
| quarantined | 0 | 0 |
| pending planned transfer | 4 | 65,320,283,335 |

The status partition reconciles exactly to all 5 planned objects and all
73,402,884,306 planned bytes. Live preflight measured 101,120,211,419,136
available bytes on the canonical MooseFS before the selected transfer, versus
8,082,600,971 selected bytes plus a 5-GiB free-space reserve. No former
quota-unavailable or global-byte policy was used.

The repository and canonical outcome manifests have identical SHA-256
`2c6b9ff843a98fc309aa124ce0337519f95cbfeb760ed497041006a5966d4a3e`.
Downstream work should consume the canonical copy at
`/moosefs/erikg/vgp/pilot/manifests/vgp_validation_reads_manifest_v1.json` and
the verified accession view, not a staging path.

## Revalidation and resume

Offline full-payload revalidation (no network):

```sh
python3 analysis/acquire_vgp_validation_reads.py --verify-only
```

Resume every pending planned P04/P07 object sequentially:

```sh
python3 analysis/acquire_vgp_validation_reads.py
```

The pure Guix Python trust store does not include the cluster's self-signed TLS
chain and correctly refused ENA at zero bytes. Acquisition on this cluster must
therefore use host Python/site trust as shown; TLS verification must not be
disabled. Interrupted transfers remain under the configured canonical staging
tree, and checksum/size/gzip failures are atomically quarantined rather than
promoted.
