## Resubmission
This is a resubmission. In this version we have:

* Updated the source URL of funGp's logo in README.md

* Added a BugReports section in DESCRIPTION

* Edited the line type in one plot to improve its visualization

* Added trace, pbars and control.optim arguments to some functions to improve
  the user control over displays

* Added two slots to the fgpm class to make it more informative

* Enabled forwarding of the control argument to optim through fgpm calls

* Protected pre-existent future backend registers as suggested by the
  author of future

* Also protected pre-existent foreach adapters as suggested by the
  author of future

* Corrected an issue with an URL from the previous submission attempt of this
  version

## Test environments
* local Windows install, R 4.2.0
* rhub, Fedora Linux, R-devel, clang, gfortran
* rhub, Ubuntu Linux 20.04.1 LTS, R-release, GCC
* rhub, Windows Server 2022, R-devel, 64 bit

Also checked with --run-donttest

## R CMD check results
There were no ERRORs, NOTEs or WARNINGs.

## Downstream dependencies
There are currently no downstream dependencies for this package.
