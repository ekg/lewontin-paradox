;; Hermetic build and smoke environment for SweepGA origin/main.
;; Resolve only through analysis/guix/channels.scm.
(set! %load-path
      (cons "analysis/guix" %load-path))

(use-modules (guix profiles)
             (packages rust-binary))

(packages->manifest
 (cons rust-toolchain-1.94.1-binary
       (map specification->package+output
            '("bash-minimal"
              "coreutils"
              "diffutils"
              "findutils"
              "gawk"
              "grep"
              "sed"
              "tar"
              "git-minimal"
              "make"
              "gcc-toolchain"
              "pkg-config"
              "cmake"
              "clang"
              "htslib"
              "gsl"
              "jemalloc"
              "glibc:static"
              "libdeflate"
              "zlib"
              "zstd"
              "bzip2"
              "xz"
              "openssl"
              "nss-certs"))))
