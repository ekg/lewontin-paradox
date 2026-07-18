(define-module (vgp_10_pilot packages psmc)
  #:use-module (guix build-system gnu)
  #:use-module (guix download)
  #:use-module (guix gexp)
  #:use-module ((guix licenses) #:prefix license:)
  #:use-module (guix packages)
  #:use-module (gnu packages compression)
  #:use-module (gnu packages perl))

(define-public psmc-vgp-pinned
  (package
    (name "psmc-vgp-pinned")
    (version "0.6.5")
    (source
     (origin
       (method url-fetch)
       ;; The immutable tag resolves to this exact commit.  Commit and archive
       ;; digests are duplicated in environment-lock.json for audit tooling.
       (uri "https://github.com/lh3/psmc/archive/b37b1cfa05b89c67c2ad1b63c699a27600d5516e.tar.gz")
       (sha256
        (base32 "0h40qiq6xx4mczrs0v8qv9bmdv1ggm48f30h524zjk2s9n5gyvag"))))
    (build-system gnu-build-system)
    (arguments
     (list
      #:tests? #f                    ; upstream has no test target
      #:make-flags #~(list "CC=gcc")
      #:phases
      #~(modify-phases %standard-phases
          (delete 'configure)
          (add-after 'build 'build-utils
            (lambda _
              (invoke "make" "-C" "utils" "fq2psmcfa" "splitfa")))
          (replace 'install
            (lambda _
              (let ((bin (string-append #$output "/bin"))
                    (share (string-append #$output "/share/psmc")))
                (mkdir-p bin)
                (mkdir-p share)
                (install-file "psmc" bin)
                (for-each (lambda (file) (install-file file bin))
                          (find-files "utils" "^(fq2psmcfa|splitfa)$"))
                (for-each (lambda (file) (install-file file share))
                          (find-files "utils" "\\.(pl|py|R)$"))))))))
    (inputs (list zlib))
    (native-inputs (list perl))
    (home-page "https://github.com/lh3/psmc")
    (synopsis "Pairwise sequentially Markovian coalescent model")
    (description
     "This source-pinned PSMC build is used only for unscaled inference in the
VGP ten-pair pilot.  Scenario scaling remains a separate workflow stage.")
    (license license:expat)))
