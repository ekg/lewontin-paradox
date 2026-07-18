;; Minimal, pinned acquisition environment for the VGP Freeze 1 mirror.
;; Evaluate only through analysis/guix/channels.scm.
(use-modules (guix profiles)
             (gnu packages base)
             (gnu packages bash)
             (gnu packages python)
             (gnu packages rsync))

(packages->manifest
 (list bash-minimal
       coreutils
       python
       rsync))
