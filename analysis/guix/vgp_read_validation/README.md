# Pinned VGP raw-read validation environment

This manifest is evaluated only through
`analysis/guix/vgp_10_pilot/channels.scm`, whose Guix channel is frozen at
commit `44bbfc24e4bcc48d0e3343cd3d83452721af8c36`. It adds Minimap2 2.24 and
Jellyfish 2.3.0 to the already pinned PSMC 0.6.5, HTSlib/SAMtools/BCFtools,
BEDTools, Python, SciPy, and test stack.

The realized profile, derivation, closure digest, package versions, and
executable digests are captured beneath the canonical shared VGP root at
`/moosefs/erikg/vgp/derived/read-validation/environment/`. Biological inputs
and derived products are never placed beneath a Lewontin-paradox-named data
directory.
