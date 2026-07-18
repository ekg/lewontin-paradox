;; Pinned metadata-only preflight environment for pilot-vgp-phylo-gbgc.
;; This is deliberately not a biological analysis environment: the sequence,
;; checksum, coordinate, and orthology gates failed before alignment/modeling.
(use-modules (guix profiles)
             (gnu packages base)
             (gnu packages bash)
             (gnu packages check)
             (gnu packages python)
             (gnu packages python-xyz))

(packages->manifest
 (list bash-minimal
       coreutils
       python
       python-pytest))
