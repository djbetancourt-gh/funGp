## Resubmission
This is a resubmission. In this version we have:

* Improved the display of progress bars for the fgpm() and
  fgpm_factory() functions.

* Corrected two functionality bugs concerning particular uses of
  fgpm() and fgpm_factory().
  
* Corrected a mistake in the ordering of a data structure delivered
  by fgpm_factory().
  
* Improved the display of the legend in plotX().

## Test environments
* local Linux install, R 3.6.3
* rhub, Fedora Linux, R-devel, clang, gfortran
* rhub, Ubuntu Linux 16.04 LTS, R-release, GCC
* rhub, Windows Server 2008 R2 SP1, R-devel, 32/64 bit

Also checked with --run-donttest

## R CMD check results
There were no ERRORs, WARNINGs or NOTEs.

## Downstream dependencies
There are currently no downstream dependencies for this package.
