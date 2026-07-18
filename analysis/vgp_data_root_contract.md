# VGP data root contract

Date: 2026-07-18T17:38:47Z

Versioned config artifact: `analysis/vgp_data_root_config.json`

## Root status

- Authorized root: `/moosefs/erikg/vgp`
- Owner/group: `erikg:erikg`
- Root mode: `0o2770`
- World-writable: `false`
- Downstream acquisition ready: `true`

## Durable layout

| key | relative path | mode | purpose |
|---|---|---|---|
| `manifests` | `manifests` | `0o2770` | metadata-only release and pilot manifests |
| `immutable_objects` | `objects/sha256` | `0o2770` | content-addressed verified immutable objects |
| `accession_views` | `views/accession` | `0o2770` | accession.version views that resolve only to verified immutable objects |
| `version_views` | `views/version` | `0o2770` | inventory or release version views that resolve only to verified immutable objects |
| `staging` | `staging` | `0o2770` | bounded mutable staging root |
| `staging_acquisition` | `staging/acquisition` | `0o2770` | temporary acquisition workspace before verification |
| `staging_partials` | `staging/partials` | `0o2770` | the only location where partial file bytes may exist |
| `staging_outputs` | `staging/outputs` | `0o2770` | temporary validated-output promotion workspace |
| `quarantine` | `quarantine` | `0o2770` | failed checksum, size, provenance, or format candidates |
| `logs` | `logs` | `0o2770` | transfer, validation, and pilot logs |
| `locks` | `locks` | `0o2770` | advisory lock files for acquisition and promotion coordination |
| `pilot` | `pilot` | `0o2770` | bounded pilot-only outputs and run packets |
| `pilot_manifests` | `pilot/manifests` | `0o2770` | pilot manifest snapshots used by later tasks |
| `pilot_runs` | `pilot/runs` | `0o2770` | per-run telemetry packets and run-local metadata |
| `pilot_outputs` | `pilot/outputs` | `0o2770` | validated pilot outputs promoted after checks pass |

## Transfer contract

- Acquisition writes mutable bytes only under `staging/acquisition` and `staging/partials`; partial files may exist only under `staging/partials`.
- Promotion target for immutable verified content is `objects/sha256`.
- Accession views live under `views/accession` and version views live under `views/version`; both must resolve only to verified immutable objects.
- Before promotion, acquisition must verify source locator, accession.version, byte size, and checksum.
- Mount or read only verified immutable inputs or accession views that resolve to verified immutable inputs.
- Copy mutable working sets to SLURM_TMPDIR before compute stages write temporary data.
- Compute arrays never download biological payloads or reach out to providers.
- Validate outputs in scratch or staging, then atomically promote only validated final outputs into pilot/outputs.
- Directory scans and cleanup must remain inside the configured VGP_DATA_ROOT and may not traverse unrelated project data.

## Validation evidence

### Small-file, fsync, checksum, atomic rename, lock, cleanup

- File fsync: `pass`
- Staging dir fsync: `pass`
- Atomic promotion: `pass`
- Same filesystem: `true`
- Checksum verification: `pass`
- Lock behavior: `pass`
- Cleanup: `pass`

### Filesystem and MooseFS evidence

- `2026-07-18T17:38:47Z` `/bin/ls -ld /moosefs/erikg/vgp` exit=`0`

```text
drwxrws--- 10 erikg erikg 3001182 Jul 18 17:30 /moosefs/erikg/vgp
```
- `2026-07-18T17:38:47Z` `/bin/findmnt -T /moosefs/erikg/vgp` exit=`0`

```text
TARGET   SOURCE             FSTYPE OPTIONS
/moosefs mfs#octopus04:9521 fuse   rw,nosuid,nodev,noatime,user_id=0,group_id=0,allow_other
```
- `2026-07-18T17:38:47Z` `/usr/bin/stat -f -c 'fstype=%T block_size=%S blocks=%b blocks_available=%a files=%c files_free=%d' /moosefs/erikg/vgp` exit=`0`

```text
fstype=fuseblk block_size=65536 blocks=3163819648 blocks_available=1542289743 files=1014444374 files_free=1004580301
```
- `2026-07-18T17:38:47Z` `/bin/df -hP /moosefs/erikg/vgp` exit=`0`

```text
Filesystem          Size  Used Avail Use% Mounted on
mfs#octopus04:9521  189T   97T   92T  52% /moosefs
```
- `2026-07-18T17:38:47Z` `/bin/df -iP /moosefs/erikg/vgp` exit=`0`

```text
Filesystem             Inodes   IUsed      IFree IUse% Mounted on
mfs#octopus04:9521 1014444374 9864073 1004580301    1% /moosefs
```
- `2026-07-18T17:38:47Z` `/usr/bin/mfsmount -V` exit=`0`

```text
fusermount version: 2.9.9
```
```text
LizardFS version 3.12.0-devel
FUSE library version: 2.9.9
```

## Blockers

- None.

## Confidence warnings

- `QUOTA_INTERFACE_UNAVAILABLE`: No user-visible quota helper was found; direct filesystem headroom remains authoritative.

## Git hygiene note

- No biological bulk data, symlinked assembly tree, or generated biological asset was added to Git by this task.
