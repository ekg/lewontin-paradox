;; Build and runtime tools for the SweepGA -> IMPG executable handoff proof.
;; Resolve this manifest only through ../channels.scm.
(set! %load-path
      (cons "analysis/guix" %load-path))

(use-modules (guix profiles)
             (packages rust-binary))

(packages->manifest
 (cons rust-toolchain-1.94.1-binary
       (map specification->package+output
            '("bash-minimal"
              "coreutils"
              "findutils"
              "gawk"
              "grep"
              "sed"
              "python"
              "samtools"
              "bcftools"
              "git-minimal"
              "make"
              "gcc-toolchain"
              "gsl"
              "jemalloc"
              "pkg-config"
              "cmake"
              "clang"
              "htslib"
              "glibc:static"
              "libdeflate"
              "zlib"
              "zstd"
              "bzip2"
              "xz"
              "openssl"
              "nss-certs"))))
