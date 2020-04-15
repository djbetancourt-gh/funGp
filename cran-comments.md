## Resubmission
This is a resubmission. In this version we have:

* Added the Value field previously missing in the documentation
  of some exported methods, e.g., in decay.Rd.
  
* Replaced cat(..) with message(..) and also added the boolean
  param 'trace' in some functions in order to allow the user to
  easily suppress messages printed to console,
  e.g., in 3_training_S.R.
  
* Added the boolean pbars param in relevant functions to allow
  easy suppression of progress bars.

## Test environments
* local Linux install, R 3.6.3
* rhub, Fedora Linux, R-devel, clang, gfortran
* rhub, Ubuntu Linux 16.04 LTS, R-release, GCC
* rhub, Windows Server 2008 R2 SP1, R-devel, 32/64 bit

Also checked with --run-donttest

## R CMD check results
There were no ERRORs or WARNINGs. 

There were 2 NOTEs:

* New submission

This is indeed the first submission of this package.

* Possibly mis-spelled words in DESCRIPTION:
  Betancourt (20:2)
  Metamodeling (19:2)
  al (20:16)
  et (20:13)

"Betancourt"" is a name, "Metamodeling"" is a valid word used
in computer experiments, "et" and "al" make part of the
reference.

## Downstream dependencies
There are currently no downstream dependencies for this package.
