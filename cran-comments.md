## Resubmission
This is a resubmission. In this version we have:

* Replaced the `djbetancourt@uninorte.edu.co` belonging to the
  package maintainer with the new `fungp.rpack@gmail.com` one.
  The former is about to turn inactive. The later, being a personal
  `@gmail` account has granted long-term functionality.

* Replaced the depreacted `@docType "package"` in `0_funGp_Doc.R`
  with the new `"_PACKAGE"` sentinel.

* Updated the documentation of `fgpm-class`, `kernel-class` and
  `proj-Class` to escape some textual braces.

* Resolved some NOTEs listed in CRAN related to RMarkdown and
  Rozygen syntax.

* Updated a few links to reference files which now point to the
  general HAL site instead of the French ramification.


## Test environments
* local Windows install, R 4.3.1
  - 0 errors ✔ | 0 warnings ✔ | 0 notes ✔
* rhub, Fedora Linux, R-devel, clang, gfortran
  - 0 errors ✔ | 0 warnings ✔
  - Gave the following note, which is unrelated to funGp:
  ```
  * checking HTML version of manual ... NOTE
  Skipping checking HTML validation: no command 'tidy' found
  Skipping checking math rendering: package 'V8' unavailable
  ```
* rhub, Ubuntu Linux 20.04.1 LTS, R-release, GCC
  - 0 errors ✔ | 0 warnings ✔
  - Gave the following note, which is unrelated to funGp:
  ```
  * checking HTML version of manual ... NOTE
  Skipping checking HTML validation: no command 'tidy' found
  Skipping checking math rendering: package 'V8' unavailable
  ```
* rhub, Windows Server 2022, R-devel, 64 bit
  - 0 errors ✔ | 0 warnings ✔
  - Gave following notes, which are unrelated to funGp:
    ```
    * checking HTML version of manual ... NOTE
    Skipping checking math rendering: package 'V8' unavailable
    
    * checking for non-standard things in the check directory ... NOTE
    Found the following files/directories:
    ''NULL''
    
    * checking for detritus in the temp directory ... NOTE
    Found the following files/directories:
    'lastMiKTeXException'
    ```
    
  - All environments gave the following note related to the
    maintainer's email update pointed out above in the
    resubmission notes:
    ```
    * checking CRAN incoming feasibility ... [17s] NOTE
    New maintainer:
      Jose Betancourt <fungp.rpack@gmail.com>
    Old maintainer(s):
      Jose Betancourt <djbetancourt@uninorte.edu.co>
    ```

Also checked with --run-donttest


## R CMD check results
There were no ERRORs or WARNINGs. Only the following NOTE
related to the maintainer's email update pointed out above
in the resubmission notes:
  ```
  New maintainer:
    Jose Betancourt <fungp.rpack@gmail.com>
  Old maintainer(s):
    Jose Betancourt <djbetancourt@uninorte.edu.co>
  ```


## Downstream dependencies
There are currently no downstream dependencies for this package.
