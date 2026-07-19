;; Pinned execution environment for the VGP Freeze 1 BGZF derived view.
;; Evaluate only through analysis/guix/channels.scm at commit
;; 44bbfc24e4bcc48d0e3343cd3d83452721af8c36.
(use-modules (guix profiles)
             (gnu packages base)
             (gnu packages bash)
             (gnu packages bioinformatics)
             (gnu packages check)
             (gnu packages compression)
             (gnu packages gawk)
             (gnu packages linux)
             (gnu packages python)
             (gnu packages python-xyz))

(packages->manifest
 (list bash-minimal
       coreutils
       gzip
       gawk
       htslib
       util-linux
       samtools
       python
       python-pytest))
