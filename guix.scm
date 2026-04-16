(use-modules (guix git)
             (guix packages)
             (guix licenses)
             (guix build-system pyproject)
             (gnu packages)
             (gnu packages python-build)
             (gnu packages python-xyz)
             (gnu packages python-science)
             (guix-arg packages python-extra))

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
                           ; python-pcangsd
                           python-scipy))
  (native-inputs (list python-setuptools))
  (home-page "https://github.com/alxsimon/local_pcangsd")
  (synopsis
   "")
  (description
   "")
  (license gpl3+))
