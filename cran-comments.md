## Resubmission
This is a resubmission. In this version we have:

* Fixed an invalid `\doi` specification that was present in multiple
  md files. We corrected it as requested by the CRAN team.

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

* Replaced the `\url` and `\href` references of the funGp manual
  with the `\doi` one pointing to the upcoming JSS publication.
  
* The DOI in the CITATION is for a new JSS publication that will
  be registered after publication on CRAN.


## Test environments
* local Windows install, R 4.3.1
  - 0 errors ✔ | 0 warnings ✔ | 0 notes ✔
  - Status: OK
* rhub v2, linux, R-* (any version), ubuntu-latest on GitHub
  - 0 errors ✔ | 0 warnings ✔ | 0 notes ✔
  - Status: OK
* rhub v2, macos, R-* (any version), macos-13 on GitHub
  - 0 errors ✔ | 0 warnings ✔ | 0 notes ✔
  - Status: OK
* rhub v2, windows, R-* (any version), windows-latest on GitHub
  - 0 errors ✔ | 0 warnings ✔ | 0 notes ✔
  - Status: OK
* rhub v2, atlas, R-devel (2024-05-07 r86527), Fedora Linux 38 (Container Image)
  - 0 errors ✔ | 0 warnings ✔ | 0 notes ✔
  - Status: OK
* rhub v2, ubuntu-release, R-4.4.0 (2024-04-24), Ubuntu 22.04.4 LTS
  - 0 errors ✔ | 0 warnings ✔ | 0 notes ✔
  - Status: OK
    
  - Previous checks with Rhub v1 gave the following note related
    to the maintainer's email update pointed out above in the
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
  - 0 errors ✔ | 0 warnings ✔ | 0 notes ✔
  - Status: OK


## Downstream dependencies
There are currently no downstream dependencies for this package.
