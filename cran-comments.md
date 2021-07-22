## Resubmission
This is a resubmission. Previous version (0.2.1) was archived on
CRAN due to a check ERROR not corrected in time. In this version
we have:

* Resolved the line that was causing the check ERROR displayed for
  version 0.2.1.
  
* Replaced the \donttest command around long examples by \dontrun
  to prevent issues during CRAN revision.
  
* Improved the documentation to simplify help pages indexing

## Test environments
* local Windows install, R 4.1.0
* rhub, Fedora Linux, R-devel, clang, gfortran
* rhub, Ubuntu Linux 20.04.1 LTS, R-release, GCC
* rhub, Windows Server 2008 R2 SP1, R-devel, 32/64 bit

Also checked with --run-donttest

## R CMD check results
There were no ERRORs or WARNINGs.
There were some NOTEs indicating:
  - Package was archived on CRAN on 2021-04-29 as check problems
  were not corrected in time
  - Possibly mis-spelled words in DESCRIPTION:
    - Betancourt (my last name)
    - Metamodeling (valid technical word)
    - al (standard part of a reference)
    - et (standard part of a reference)

## Downstream dependencies
There are currently no downstream dependencies for this package.
