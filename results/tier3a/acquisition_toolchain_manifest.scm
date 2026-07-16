;;; Guix build/runtime libraries for the staged SweepGA and IMPG executables.
;;; Evaluate only with analysis/guix/channels.scm (Guix 44bbfc24e4bcc48d).
(use-modules (guix profiles))

(specifications->manifest
 '("bash-minimal"
   "coreutils"
   "findutils"
   "git"
   "make"
   "nss-certs"
   "sed"
   "gcc-toolchain"
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
   "gsl"
   "jemalloc"))
