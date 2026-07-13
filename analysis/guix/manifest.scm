;; Reproducible Tier 3 environment.  Evaluate only through channels.scm.
(set! %load-path
      (cons "analysis/guix" %load-path))

(use-modules (guix packages)
             (guix profiles)
             (gnu packages bash)
             (gnu packages bioinformatics)
             (gnu packages certs)
             (gnu packages check)
             (gnu packages maths)
             (gnu packages python)
             (gnu packages python-build)
             (gnu packages python-science)
             (gnu packages python-xyz)
             (gnu packages statistics)
             (packages impg)
             (packages wfmash))

;; SciPy 1.10.1 in the frozen channel propagates Matplotlib even though SciPy
;; itself does not import it.  That pulls an obsolete GUI/WebKit stack into a
;; headless analysis profile.  Keep the exact frozen SciPy source/build recipe,
;; retain Pythran and pyparsing while building, and expose only its true NumPy
;; runtime dependency.
(define python-scipy-tier3
  (package/inherit python-scipy
    (name "python-scipy-tier3")
    (propagated-inputs (list python-numpy))
    (native-inputs
     (modify-inputs (package-native-inputs python-scipy)
       (append python-pythran python-pyparsing)))))

(packages->manifest
 (list python
       bash-minimal
       nss-certs
       python-pytest
       python-pysam
       python-pyfaidx
       python-biopython
       python-numpy
       python-pandas
       python-scipy-tier3
       python-jsonschema
       bcftools
       samtools
       htslib
       bedtools
       vcftools
       wfmash-tier3
       ;; Source-only and deliberately has no `impg` executable.  See package.
       impg-source-pinned))
