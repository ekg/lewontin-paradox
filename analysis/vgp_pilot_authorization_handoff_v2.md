# VGP ten-pair execution authorization handoff v2

Authorization ID: `vgp10-auth-20260718-v2`

Schema: `vgp-pilot-execution-authorization-v2.0.0`

Current decision: **GO for core execution on all 10 frozen primaries**

Canary: **P07, *Spinachia spinachia***

Biological jobs submitted by this authorization task: **0**

The machine-readable authority is `analysis/vgp_pilot_authorization_v2.json`.
The read-only live preflight is
`analysis/vgp_pilot_authorization_preflight_v2.json`. Historical NO_GO,
CONDITIONAL_GO, zero-job, and review artifacts remain unchanged evidence of
earlier runs. They do not supply or override this authorization.

The active shared root is `/moosefs/erikg/vgp`. The former
`lewontin-paradox-data/vgp/phase1-freeze-1.0` tree is migration input only.
`vgp_data_root_migration_v1.json` records 140 unique verified objects
(12,692,829,704 logical bytes) hard-linked into the canonical digest paths,
with zero downloads and zero source removals. No active authorization,
compute, mirror, raw/validation, or derived-output path resolves to the legacy
tree.

## Exact authorized roster

Every row below passed exact same-individual provenance, accession-version,
accepted CAS size/SHA-256, complete gzip/FASTA readability, and mutual assembly
span comparability. Confidence is currently C because optional validation is
incomplete; C is not a refusal state.

| pair | species | exact H1 / H2 | measured compressed CAS bytes | frozen assembly bp | confidence | core |
|---|---|---|---:|---:|---|---|
| P01 | Camelus dromedarius | `GCA_036321535.1` / `GCA_036321565.1` | 1,305,801,004 | 4,393,098,473 | C | **AUTHORIZED** |
| P02 | Pseudorca crassidens | `GCA_039906515.1` / `GCA_039906525.1` | 1,508,064,765 | 5,353,123,058 | C | **AUTHORIZED** |
| P03 | Colius striatus | `GCA_028858725.2` / `GCA_028858625.2` | 713,759,958 | 2,365,968,587 | C | **AUTHORIZED** |
| P04 | Falco naumanni | `GCA_017639655.1` / `GCA_017639645.1` | 745,921,893 | 2,372,467,332 | C | **AUTHORIZED** |
| P05 | Candoia aspera | `GCA_035149785.1` / `GCA_035125265.1` | 966,840,519 | 3,060,141,947 | C | **AUTHORIZED** |
| P06 | Dendropsophus ebraccatus | `GCA_027789765.1` / `GCA_027789725.1` | 1,317,901,491 | 4,567,738,964 | C | **AUTHORIZED** |
| P07 | Spinachia spinachia | `GCA_048126635.1` / `GCA_048127205.1` | 245,511,792 | 813,562,223 | C | **AUTHORIZED** |
| P08 | Menidia menidia | `GCA_048628825.1` / `GCA_048544195.1` | 338,655,435 | 1,124,806,929 | C | **AUTHORIZED** |
| P09 | Heterodontus francisci | `GCA_036365525.1` / `GCA_036365495.1` | 2,948,092,316 | 11,209,934,419 | C | **AUTHORIZED** |
| P10 | Hemiscyllium ocellatum | `GCA_020745735.1` / `GCA_020745765.1` | 2,474,337,026 | 8,132,928,280 | C | **AUTHORIZED** |

P03 and P04 lack an accepted exact annotation binding. That disables only their
annotation-derived outputs. It does not affect their core diversity or
unscaled PSMC authorization. Missing exact-final QV, complete two-haplotype
BUSCO, k-mer/copy-number audit, standalone repeat report, raw-read validation,
exact P09 chemistry, and independent Ne evidence are explicit confidence
covariates. Even when all optional fields are set missing, the deterministic
fixture retains ten authorized jobs (`missing_optional_qc_global_job_count_effect=0`).

## Authorization boundary

The only hard pre-execution gates are:

1. exact same-individual H1/H2 identity and accession provenance;
2. accepted CAS byte sizes and SHA-256 digests;
3. readable, mutually comparable assemblies;
4. the captured accepted SweepGA and IMPG executable bytes;
5. writable storage with headroom; and
6. a live Slurm controller, suitable partition, and scheduler-accepted packet.

The hard per-pair result gates remain alignment success, retained query and
target multiplicity no greater than one, exact H1 REF checks, manifest-bound
H2 reconstruction, exact ordered mask accounting, at least 100,000,000
callable bases, and at least 60% callable H1 sequence. A pair that fails one of
those result gates stops and cannot be reported as a diversity/PSMC success;
it does not retroactively deauthorize the other pairs.

When no standalone repeat report exists, mapping generates a deterministic
assembly-derived simple-repeat/low-complexity mask (long homopolymers and exact
di-/trinucleotide tandems) and combines it with sequence gap, mapping, and any
manifested masks under the fixed first-reason-wins accounting order. Absolute
PSMC scaling and annotation are optional branches. Their absence preserves the
unscaled PSMC output and never creates a zero estimate.

## Canary selection and exact packet

P07 is the canary because its immutable summed compressed CAS input is
245,511,792 bytes, the minimum across the ten primaries; selection ID is the
predeclared tie-break. This uses no mutable benchmark or outcome.

The initial packet is:

```sh
sbatch analysis/slurm/vgp_10_pilot/authorized/v2.0.0/P07.canary.sbatch
```

It requests the live `highmem` partition, one node, 32 CPUs, 128 GiB RAM, and
72 hours. The worker authenticates and decompresses both CAS objects into
`$SLURM_TMPDIR`, performs high-I/O mapping/IMPG/normalization/consensus/PSMC
there, and checkpoints every completed stage plus every ten PSMC replicates to:

```text
/moosefs/erikg/vgp/pilot/authorized-checkpoints/vgp10-auth-20260718-v2/P07
```

On an observed OOM, submit the next immutable packet; it resumes from that
checkpoint without relaxing a scientific gate:

```sh
sbatch analysis/slurm/vgp_10_pilot/authorized/v2.0.0/P07.canary.mem256.sbatch
sbatch analysis/slurm/vgp_10_pilot/authorized/v2.0.0/P07.canary.mem512.sbatch
```

Do not submit both retries speculatively. Use 256 GiB after a 128-GiB OOM and
512 GiB only after a 256-GiB OOM. The worker records retry index, Slurm job ID,
scratch staging, and checkpoint time.

## Completed non-biological preflight

The live preflight rehashed every accepted canonical CAS object, performed a
complete gzip-stream integrity read of all twenty FASTAs, and verified a
nonempty FASTA header and sequence row for each. The exact authorization JSON
is bound into that report by payload SHA-256
`f8145b384b5abe4282fe9d9732a1858b6e668bdf8e0863005f5ae512154c2b15`.
It authenticated SweepGA SHA-256
`fa7f0edb9b7e275c288db254046020e136d4267dd5ee043379227ef80da0573b`
and IMPG SHA-256
`c587dc2326cd24f887b1fcb3938404229ad0f0a27ef0773e90c287b1ade160d4`
through the pinned Guix capture. The Slurm controller reported UP; `highmem`
reported more than 128 GiB/node; storage was writable with
101,119,647,088,640 free bytes versus a conservative 687,194,767,360-byte
ten-pair scratch envelope.

All 128-, 256-, and 512-GiB packets passed `sbatch --test-only`. Slurm printed
hypothetical test-only job numbers 1781493–1781495, but none appeared in
`squeue`; `slurm_jobs_submitted=0` and `biological_execution=false` are recorded
in the preflight JSON.

## Reproduction

Regenerate and validate the immutable packet:

```sh
python3 -m analysis.vgp_pilot_authorization generate
python3 -m analysis.vgp_pilot_authorization validate
```

Repeat the read-only CAS, environment, storage, and scheduler preflight before
submission:

```sh
python3 -m analysis.vgp_pilot_authorization preflight
```

The deterministic authorization fixtures and full analysis suite run through
the repository's pinned GNU Guix time machine; no host Python environment is
accepted as substitute evidence.
