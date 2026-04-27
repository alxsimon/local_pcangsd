(use-modules (guix git)
             (guix git-download)
             (guix packages)
             ((guix licenses) #:prefix license:)
             (guix build-system pyproject)
             (gnu packages)
             (gnu packages python-build)
             (gnu packages python-xyz)
             (gnu packages python-science)
             (guix-arg packages python-extra))
             ; (guix-science packages bioinformatics))

(define-public python-pcangsd
  (package
    (name "python-pcangsd")
    (version "1.36.4")
    (source
     (origin
       (method git-fetch)
       (uri (git-reference
             (url "https://github.com/Rosemeis/pcangsd")
             (commit (string-append "v" version))))
       (file-name (git-file-name name version))
       (sha256
        (base32 "04j5545z8lz9fjac9smx095m7zdjb6gli1nmsm28r4ha4pp17g9g"))))
    (build-system pyproject-build-system)
    (arguments
      (list
        ;; No tests in package.
        #:tests? #f))
    (propagated-inputs (list python-numpy
                             python-scipy))
    (native-inputs (list python-cython-next
                         python-setuptools))
    (home-page "https://github.com/Rosemeis/pcangsd")
    (synopsis
     "Framework for analyzing low depth @acronym{NGS, Next-Generation Sequencing} data using @acronym{PCA, Principal Component Analysis}")
    (description
     "Framework for analyzing low-depth @acronym{NGS, Next-Generation Sequencing} data
in heterogeneous/structured populations using @acronym{PCA, Principal Component Analysis}.")
    (license license:gpl3)))

(package
  (name "python-local-pcangsd")
  (version "0.0.3")
  (source (git-checkout (url (dirname (current-filename)))))
  (build-system pyproject-build-system)
  (arguments '(#:tests? #f))
  (propagated-inputs (list python-numpy
                           python-dask
                           python-pandas
                           python-xarray
                           python-sgkit
                           python-pcangsd
                           python-scipy))
  (native-inputs (list python-cython-next python-setuptools))
  (home-page "https://github.com/alxsimon/local_pcangsd")
  (synopsis
   "")
  (description
   "")
  (license license:gpl3+))
