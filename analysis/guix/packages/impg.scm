(define-module (packages impg)
  #:use-module ((guix licenses) #:prefix license:)
  #:use-module (guix build-system copy)
  #:use-module (guix gexp)
  #:use-module (guix git-download)
  #:use-module (guix packages)
  #:export (impg-source-pinned))

(define %impg-commit "101df81eb28a809c8fac97d297acd9fcfbbfa048")

(define-public impg-source-pinned
  (package
    (name "impg-source-pinned")
    (version (string-append "0.4.1-" (string-take %impg-commit 7)))
    (source
     (origin
       (method git-fetch)
       (uri (git-reference
             (url "https://github.com/pangenome/impg.git")
             (commit %impg-commit)
             (recursive? #t)))
       (file-name (git-file-name name version))
       ;; Recursive submodules: vendor/gfaffix=460e0dd... and
       ;; vendor/syng=68ac197....
       (sha256
        (base32 "1g7z0xwaryz6cbkc2mhsrysi9cvg8599k50n3qh6ny4i0rx1z4di"))))
    (build-system copy-build-system)
    (arguments
     (list
      #:install-plan
      #~'(("." "share/impg-source"))
      #:phases
      #~(modify-phases %standard-phases
          (add-before 'install 'verify-immutable-lock
            (lambda _
              (unless (and (file-exists? "Cargo.lock")
                           (file-exists? "vendor/gfaffix/Cargo.toml")
                           (file-exists? "vendor/syng/Makefile"))
                (error "recursive IMPG source is incomplete")))))))
    (home-page "https://github.com/pangenome/impg")
    (synopsis "Immutable, execution-disabled Tier 3 IMPG source snapshot")
    (description
     "This source-only package preserves the recursively pinned IMPG candidate
for optional concordance work.  It intentionally installs no executable:
Cargo.lock still contains branch/unqualified Git sources, so
tier3-decisions-v1 forbids execution until every dependency is Guix-vendored
and hash-locked and the 10 kb truth test passes.  This fail-closed package
prevents an unlocked Cargo fetch or user-profile binary from enabling IMPG.")
    (license license:expat)))

impg-source-pinned
