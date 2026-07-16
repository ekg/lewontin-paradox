(define-module (packages rust-binary)
  #:use-module (guix packages)
  #:use-module (guix download)
  #:use-module (guix gexp)
  #:use-module (guix build-system gnu)
  #:use-module ((guix licenses) #:prefix license:)
  #:use-module (gnu packages base)
  #:use-module (gnu packages gcc)
  #:use-module (gnu packages elf)
  #:use-module (gnu packages compression)
  #:export (rust-toolchain-1.94.1-binary))

;; Guix commit 44bbfc2 carries patchelf 0.11, which predates the ELF layout in
;; the official Rust 1.94 release.  Build the needed newer release from its
;; authenticated source rather than reaching outside the pinned environment.
(define patchelf-for-rust
  (package
    (name "patchelf-for-rust")
    (version "0.18.0")
    (source
     (origin
       (method url-fetch)
       (uri (string-append
             "https://github.com/NixOS/patchelf/releases/download/" version
             "/patchelf-" version ".tar.bz2"))
       (sha256
        (base32 "02s7ap86rx6yagfh9xwp96sgsj0p6hp99vhiq9wn4mxshakv4lhr"))))
    (build-system gnu-build-system)
    (arguments (list #:tests? #f))
    (home-page "https://github.com/NixOS/patchelf")
    (synopsis "Modify ELF executables and libraries")
    (description "Patchelf 0.18 is used while packaging the fixed Rust binary release.")
    (license license:gpl3+)))

;; The frozen 2023 channel provides Rust 1.67, while the locked SweepGA and
;; IMPG dependency graphs require a current compiler.  Import the official
;; release as a fixed-output Guix package and rewrite every relevant ELF to the
;; frozen channel's glibc/GCC/zlib closure.  No user-profile or rustup program
;; participates in the build.
(define-public rust-toolchain-1.94.1-binary
  (package
    (name "rust-toolchain-binary")
    (version "1.94.1")
    (source
     (origin
       (method url-fetch)
       (uri (string-append
             "https://static.rust-lang.org/dist/rust-" version
             "-x86_64-unknown-linux-gnu.tar.xz"))
       (sha256
        (base32 "0h0cgn3k5i6zl4n44g9kim4mi2zbh46cd4324y0jbrkjza0ksjr9"))))
    (build-system gnu-build-system)
    (arguments
     (list
      #:tests? #f
      #:strip-binaries? #f
      ;; The channel's validator is older than this upstream ELF set.  Runtime
      ;; linkage is instead asserted by the smoke script with ldd plus actual
      ;; locked Cargo builds.
      #:validate-runpath? #f
      #:phases
      #~(modify-phases %standard-phases
          (delete 'configure)
          (delete 'build)
          (delete 'check)
          (replace 'install
            (lambda* (#:key inputs outputs #:allow-other-keys)
              (let* ((out (assoc-ref outputs "out"))
                     (glibc (assoc-ref inputs "glibc"))
                     (gcc-lib (assoc-ref inputs "gcc-lib"))
                     (zlib (assoc-ref inputs "zlib"))
                     (interpreter
                      (string-append glibc "/lib/ld-linux-x86-64.so.2"))
                     (rpath
                      (string-append out "/lib:" gcc-lib "/lib:"
                                     glibc "/lib:" zlib "/lib")))
                (invoke "bash" "install.sh"
                        (string-append "--prefix=" out)
                        "--disable-ldconfig"
                        "--components=rustc,cargo,rust-std-x86_64-unknown-linux-gnu")
                ;; The installer uses both regular files and symlinks, so ask
                ;; patchelf itself whether each regular file is an ELF.  Its
                ;; --print-interpreter exit status reliably distinguishes
                ;; executables from shared objects.
                (invoke
                 "bash" "-c"
                 (string-append
                  "while IFS= read -r -d '' f; do "
                  "  if patchelf --print-rpath \"$f\" >/dev/null 2>&1; then "
                  "    patchelf --force-rpath --set-rpath \"$1\" \"$f\"; "
                  "    if patchelf --print-interpreter \"$f\" >/dev/null 2>&1; then "
                  "      patchelf --set-interpreter \"$2\" \"$f\"; "
                  "    fi; "
                  "  fi; "
                  "done < <(find \"$3\" -type f -print0)")
                 "patch-rust-elfs" rpath interpreter out)))))))
    (native-inputs
     (list patchelf-for-rust))
    (inputs
     `(("glibc" ,glibc)
       ("gcc-lib" ,gcc-12 "lib")
       ("zlib" ,zlib)))
    (home-page "https://www.rust-lang.org")
    (synopsis "Fixed official Rust compiler and Cargo toolchain")
    (description
     "This package installs the fixed official Rust 1.94.1 binary release and
rewrites its ELF interpreter and runpaths to the authenticated Guix channel's
runtime closure.  It exists solely to build the pinned SweepGA/IMPG proof.")
    (license (list license:asl2.0 license:expat))))
