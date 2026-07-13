;; Sole authenticated channel for tier3-decisions-v1.  Do not append %default-channels.
(list
 (channel
  (name 'guix)
  (url "https://git.savannah.gnu.org/git/guix.git")
  (commit "44bbfc24e4bcc48d0e3343cd3d83452721af8c36")
  (introduction
   (make-channel-introduction
    "9edb3f66fd807b096b48283debdcddccfea34bad"
    (openpgp-fingerprint
     "BBB0 2DDF 2E5F 5EAA C74C 3F7F 27D5 86A4 552F 384F")))))
