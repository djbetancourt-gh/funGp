## Resubmission
This is a resubmission. In this version I have:

* Extended the description of the package at the Description field
  of the DESCRIPTION file. Added relevant references describing
  the techniques behind the methods in the package.

* Double checked all the files for incomplete documentation and
  unnecessary developer notes. In particular, "Fill!!!!!!!'" and
  "fill this!!!" were removed everywhere.

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
