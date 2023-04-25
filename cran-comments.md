## Resubmission
This is a resubmission. In this version we have:

* Added extend argument to update method enabling a new way of model updating
  where current hyperparameters are used as starting points for the new
  optimization.

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
