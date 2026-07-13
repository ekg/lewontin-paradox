(define-module (packages wfmash)
  #:use-module ((guix licenses) #:prefix license:)
  #:use-module (guix build-system cmake)
  #:use-module (guix gexp)
  #:use-module (guix git-download)
  #:use-module (guix packages)
  #:use-module (gnu packages algebra)
  #:use-module (gnu packages bioinformatics)
  #:use-module (gnu packages compression)
  #:use-module (gnu packages maths)
  #:use-module (gnu packages pkg-config)
  #:export (wfmash-tier3))

(define %wfmash-commit "e040aa10e87cab44ed5a4db005e784be62b0bd21")

(define-public wfmash-tier3
  (package
    (name "wfmash-tier3")
    (version (string-append "0.24.2-12." (string-take %wfmash-commit 7)))
    (source
     (origin
       (method git-fetch)
       (uri (git-reference
             (url "https://github.com/waveygang/wfmash.git")
             (commit %wfmash-commit)
             (recursive? #t)))
       (file-name (git-file-name name version))
       ;; Recursive source includes deps/WFA2-lib at
       ;; 49c255df126ee536fe92caff7a9f7c183ec3ff29.
       (sha256
        (base32 "0xvrgb7aimsdyq7kn3lzf7rrxxljs54dv81kv3apsfpcyp6nxnwj"))))
    (build-system cmake-build-system)
    (arguments
     (list
      #:configure-flags
      #~(list "-DCMAKE_BUILD_TYPE=Generic"
              "-DDISABLE_LTO=ON"
              "-DBUILD_DEPS=OFF"
              "-DVENDOR_HTSLIB=OFF"
              "-DVENDOR_EVERYTHING=OFF")
      #:phases
      #~(modify-phases %standard-phases
          (add-after 'unpack 'remove-host-specific-flags
            (lambda _
              (substitute* (find-files "." "CMakeLists\\.txt$")
                (((string-append "-march" "=native")) "")
                (((string-append "-march" "=x86-64-v3")) ""))
              ;; This revision always emits extended CIGAR, but upstream
              ;; removed the historical -4 spelling.  Preserve the frozen
              ;; command line as an explicitly documented compatibility flag.
              (substitute* "src/interface/parse_args.hpp"
                (("    // SYSTEM")
                 (string-append
                  "    args::Flag extended_cigar(alignment_opts, \"\", "
                  "\"emit extended =/X/I/D CIGAR (always enabled)\", "
                  "{'4', \"extended-cigar\"});\n\n    // SYSTEM")))))
          (replace 'check
            (lambda* (#:key tests? #:allow-other-keys)
              (when tests?
                ;; Exercise an upstream CTest with its committed data.
                (invoke "ctest" "--output-on-failure" "-R" "^wfmash-time-LPA$")
                ;; A small base-alignment smoke.  It must emit an extended
                ;; CIGAR and must never collapse it to M.
                (invoke "bash" "-c"
                        (string-append
                         "set -euo pipefail; "
                         "bin/wfmash ../source/data/LPA.subset.fa.gz "
                         "-p 80 -n 5 -t 4 -4 > tiny-extended.paf; "
                         "test -s tiny-extended.paf; "
                         "grep -Eq 'cg:Z:.*=' tiny-extended.paf; "
                         "grep -Eq 'cg:Z:.*X' tiny-extended.paf; "
                         "grep -Eq 'cg:Z:.*I' tiny-extended.paf; "
                         "grep -Eq 'cg:Z:.*D' tiny-extended.paf; "
                         "! grep -Eq 'cg:Z:[^[:space:]]*M' tiny-extended.paf"))))))))
    (inputs
     (list bzip2 gsl htslib libdeflate xz zlib))
    (native-inputs
     (list pkg-config samtools))
    (home-page "https://github.com/waveygang/wfmash")
    (synopsis "Immutable Tier 3 WFMASH base-accurate aligner")
    (description
     "This package pins the recursively fetched WFMASH source used by
tier3-decisions-v1, removes host-specific compiler tuning, runs an upstream
CTest, restores the frozen @code{-4} compatibility spelling, and requires an
extended-CIGAR truth smoke during every build.")
    (license license:expat)))

wfmash-tier3
