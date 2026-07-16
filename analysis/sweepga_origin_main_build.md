# SweepGA `origin/main` reproducible-build audit

## Outcome: blocked by the native CLI

The fetched, unmodified `origin/main` source builds reproducibly, but its native
binary rejects the mandatory literal `-n 1:1` option. The task therefore cannot
be completed and must not unblock Tier 3A remapping.

The decisive invocation used the absolute path of the newly built binary on
real Spinachia H1/H2 20-kb excerpts:

```text
/moosefs/erikg/tier3scratch/sweepga-origin-main-018e4ce/bin-1/sweepga \
  h1.fa h2.fa --output-file literal-n-1to1.paf \
  -n 1:1 --scaffold-jump 0 --threads 2
```

It exited 2 before mapping and printed:

```text
error: unexpected argument '-n' found
```

No SweepGA source patch, wrapper, alias, argument rewrite, or compatibility
shim was used. The task's build shims supply only Guix library search paths.

## Repository and source isolation

- Parent repository: `/moosefs/erikg/sweepga`
- `origin` fetch/push URL: `git@github.com:pangenome/sweepga.git`
- Command: `git -C /moosefs/erikg/sweepga fetch --no-tags origin main`
- Fetched `FETCH_HEAD` and `refs/remotes/origin/main`:
  `018e4ce49d2c125820e0ac50dc5feaa02d423683`
- Existing checkout: branch `main`, HEAD at the same commit, divergence `0 0`
  relative to `origin/main`.
- Existing uncommitted paths, preserved without reset/clean/checkout:
  `README.md`, `docs/BUILD-NOTES.md`, `install.sh`, and
  `scripts/build-clean.sh` (6 insertions and 6 deletions in total).

Two independent clean archives of that commit were extracted at:

- `/moosefs/erikg/tier3scratch/sweepga-origin-main-018e4ce/checkout-1`
- `/moosefs/erikg/tier3scratch/sweepga-origin-main-018e4ce/checkout-2`

Both `git archive` streams have SHA-256
`a87575b6c8fb6a07a5ceac8c8b17c1cea2dfcc3318b4a10ee86878c4129347ed`.
Both copies have the locked Cargo graph SHA-256
`b1aee267fbe9db8c1ed72bcf2be8b9f7b3ba2a9a71778ce35ee1fb8ac067636a`.

## Guix build and reproducibility

The tracked build definition is
`analysis/guix/sweepga_origin_main_manifest.scm`; the authenticated channel is
`analysis/guix/channels.scm` at Guix commit
`44bbfc24e4bcc48d0e3343cd3d83452721af8c36`. The resolved profile is
`/gnu/store/yfffyhdm3a9bsah4gzw9dzri623af3f6-profile`.

Important pinned tools and libraries were Rust/Cargo 1.94.1, GCC 12.3.0,
Clang 15.0.7, CMake 3.25.1, pkg-config 0.29.2, Git 2.40.1, binutils 2.38,
glibc 2.35, GSL 2.7.1, jemalloc 5.3.0, and HTSlib 1.16. The complete package
list and `guix describe` output are retained in the scratch evidence.

Run the complete two-build recipe with:

```bash
analysis/sweepga_origin_main_rebuild.sh
```

Each pass uses a separate clean archive, Cargo home, target directory, HOME,
and cache. It executes `cargo fetch --locked`, then
`cargo build --release --locked --bin sweepga`, entirely inside
`guix time-machine ... shell --pure`. Source/build prefixes are remapped and
the deployment executable is deterministically stripped. `SOURCE_DATE_EPOCH`
is the commit time, `1776293668`.

The primary executable is byte-identical across both builds:

| Build | Absolute realpath | Size | SHA-256 |
|---|---|---:|---|
| 1 | `/moosefs/erikg/tier3scratch/sweepga-origin-main-018e4ce/bin-1/sweepga` | 5,193,792 | `fa7f0edb9b7e275c288db254046020e136d4267dd5ee043379227ef80da0573b` |
| 2 | `/moosefs/erikg/tier3scratch/sweepga-origin-main-018e4ce/bin-2/sweepga` | 5,193,792 | `fa7f0edb9b7e275c288db254046020e136d4267dd5ee043379227ef80da0573b` |

The Guix ELF interpreter is
`/gnu/store/gsjczqir1wbz8p770zndrpw4rnppmxi3-glibc-2.35/lib/ld-linux-x86-64.so.2`,
and every resolved `ldd` library is in `/gnu/store`. Seven FastGA/ONE helper
executables were also byte-identical. The independently compiled `wfmash`
helper retained native C++ build metadata and differed; it is not the SweepGA
binary being compared and was not invoked by the mandatory failing probe.

## Exact-binary and CLI evidence

The executable gate is:

```bash
analysis/sweepga_origin_main_smoke.sh
```

It enters a pinned pure Guix shell, verifies the expected SHA-256, prepends only
the rebuilt binary directory to `PATH`, records `type -a`, `command -v`, and
`realpath`, extracts the real biological regions, and invokes the binary by
absolute path. The recorded values all resolve to build 1 above. `--version`
prints `sweepga 0.1.1`.

The generated help exposes only:

```text
--num-mappings <NUM_MAPPINGS>
```

It contains no `-n` alias. This is confirmed by the unmodified source in
`src/cli.rs`: its module-level design note says that no short flags live on
`AlnArgs`, and the field is declared only with
`#[clap(long = "num-mappings", ...)]`.

Before the hard-gate clarification arrived, the exact binary was also tested
with the actually supported long syntax on a dense overlapping PAF fixture.
`--num-mappings 1:1`, `5:5`, and `10:10` retained exactly 1, 5, and 10 records,
respectively, with the deterministic first input record winning the exact tie.
This validates the query:target multiplicity implementation but does not
satisfy or replace the mandatory literal short-option gate.

## IMPG disposition

The IMPG native partition, regional query, and normalized/indexed VCF/BCF proof
was not allowed to proceed after the native `-n` failure was confirmed. No
previous SweepGA binary and no `--num-mappings` substitution was used to claim
completion. Consequently all IMPG completion fields are false in
`analysis/sweepga_origin_main_build.json`.

The scratch evidence is under
`/moosefs/erikg/tier3scratch/sweepga-origin-main-018e4ce`, including build logs,
tool versions, profile package inventory, source hashes, ELF/`ldd` evidence,
help/version output, type/path/realpath evidence, and the literal-option stderr.

Resolution requires an upstream `origin/main` commit whose unmodified native
binary actually declares and accepts `-n` for `--num-mappings`, followed by a
fresh run of both tracked scripts and the IMPG end-to-end handoff. Changing the
task to accept only `--num-mappings` would be a specification decision, not a
build workaround.
