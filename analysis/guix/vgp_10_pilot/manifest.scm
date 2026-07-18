;; Production VGP ten-pair environment.  Resolve only through channels.scm.
(set! %load-path (cons "analysis/guix" %load-path))
(set! %load-path (cons "analysis/guix/vgp_10_pilot" %load-path))

(use-modules (gnu packages)
             (guix profiles)
             (vgp_10_pilot packages psmc)
             (packages rust-binary))

(packages->manifest
 (cons* psmc-vgp-pinned
        rust-toolchain-1.94.1-binary
        (map specification->package+output
             '("bash-minimal"
               "coreutils"
               "findutils"
               "gawk"
               "grep"
               "sed"
               "python"
               "python-pytest"
               "python-jsonschema"
               "samtools"
               "bcftools"
               "htslib"
               "bedtools"
               "git-minimal"
               "make"
               "gcc-toolchain"
               "gsl"
               "jemalloc"
               "pkg-config"
               "cmake"
               "clang"
               "glibc:static"
               "libdeflate"
               "zlib"
               "zstd"
               "bzip2"
               "xz"
               "openssl"
               "nss-certs"))))
