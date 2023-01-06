## Resubmission
This is a resubmission. In this version we have:

* Updated the source URL of funGp's logo in README.md

* Edited the line type in one plot to improve its visualization

* Added trace, pbars and control.optim arguments to some functions to improve
  the user control over displays

* Added two slots to the fgpm class to make it more informative

* Enabled forwarding of the control argument to optim through fgpm calls

* Protected pre-existent future backend registers as suggested by the
  author of future




## Test environments
* local Windows install, R 4.2.0
* rhub, Fedora Linux, R-devel, clang, gfortran
* rhub, Ubuntu Linux 20.04.1 LTS, R-release, GCC
* rhub, Windows Server 2022, R-devel, 64 bit

Also checked with --run-donttest

## R CMD check results
There were no ERRORs or WARNINGs.
There were some NOTEs indicating:
  - Possibly invalid URLs:
    URL: https://www.sciencedirect.com/science/article/abs/pii/S0951832019301693
    From: man/fgpm.Rd
          man/fgpm_factory.Rd
          man/funGp-package.Rd
    Status: 403
    Message: Forbidden
    - Confirmed directly from the NOTEs that the URL is working fine
  - checking for detritus in the temp directory ... NOTE
    - Only on Windows Server 2022, R-devel, 64 bit
    - As noted in R-hub issue #503, this could be due to a bug/crash in MiKTeX
      and can likely be ignored

## Downstream dependencies
There are currently no downstream dependencies for this package.
