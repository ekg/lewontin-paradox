;; Raw-read validation environment for the VGP pilot.  Resolve only through
;; ../vgp_10_pilot/channels.scm (frozen Guix commit 44bbfc24e4bcc48...).
(set! %load-path (cons "analysis/guix" %load-path))
(set! %load-path (cons "analysis/guix/vgp_10_pilot" %load-path))

(use-modules (gnu packages)
             (gnu packages python)
             (gnu packages python-build)
             (gnu packages python-science)
             (gnu packages python-xyz)
             (guix packages)
             (guix profiles)
             (vgp_10_pilot packages psmc))

;; The frozen SciPy propagates Matplotlib and an obsolete GUI/WebKit stack even
;; though validation imports only SciPy/NumPy.  Retain the exact frozen recipe
;; and its real NumPy runtime, matching the established Tier 3 environment.
(define python-scipy-vgp-read-validation
  (package/inherit python-scipy
    ;; Keep the established package identity so this exact slim frozen SciPy
    ;; realization is reused across Tier 3 and VGP validation profiles.
    (name "python-scipy-tier3")
    (propagated-inputs (list python-numpy))
    (native-inputs
     (modify-inputs (package-native-inputs python-scipy)
       (append python-pythran python-pyparsing)))))

(packages->manifest
 (cons* psmc-vgp-pinned
        python-scipy-vgp-read-validation
        (map specification->package+output
             '("bash-minimal"
               "coreutils"
               "findutils"
               "gawk"
               "grep"
               "sed"
               "gzip"
               "time"
               "python"
               "python-pytest"
               "python-jsonschema"
               "python-pysam"
               "minimap2"
               "jellyfish"
               "samtools"
               "bcftools"
               "htslib"
               "bedtools"))))
