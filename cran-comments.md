## Resubmission
This is a resubmission. In this version we have:

* Implemented plot and summary methods for our most important
  objects to adjust to the standard use in R packages
  
* Removed our dependency on the qdapRegex package which is
  already scheduled for archival by CRAN
  
* Improved the documentation to simplify help pages indexing

## Test environments
* local Windows install, R 4.1.0
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
    - might be related to the first NOTE

## Downstream dependencies
There are currently no downstream dependencies for this package.
