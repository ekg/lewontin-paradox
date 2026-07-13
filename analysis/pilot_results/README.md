# Tier 3 pilot records

`pilot_registry.json` predeclares the five scheduler pilots without promoting
survey-level availability into analysis eligibility.  In particular, DGRP,
Ag1000G, and the VGP individual remain fail-closed until their exact reference,
sample, denominator, native-annotation, checksum, and access tuples are locked.

Runtime records written here are small provenance/QC products only:

- `compute_smoke.json` is the real compute-node Guix/store and synthetic truth
  smoke record;
- `<tier>/<dataset>.run.json` is the atomic stage/checksum run record when the
  workflow `output_root` is this directory;
- `logs/` receives Slurm stdout/stderr.

Raw FASTA, GFF, VCF/BCF, PAF, masks, and scratch-stage directories must remain
under the declared Tier 3 data/scratch roots and must not be committed.
