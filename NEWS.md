# funGp 0.3.1 (2023-01-02)

## Minor display changes

* Line type of the mean function in `plotSims.fgpm` changed to dashed
  (`lty` = 2) for easier inspection.

* A `trace` argument was added to `update.fgpm`.

* Convergence code and negated log-likelihood value added as slots of
  fgpm objects to make them more informative.

* Forwarding of `optim` `control` argument enabled in fgpm through
  `control.optim`.

* Following the addition of control.optim, we improved the functionality
  of the `fgpm` and `fgpm_factory` arguments to control the display of
  `funGp` native progress messages

* Pre-existent future backend registers protected as suggested by the
  author of future in a previous issue



# funGp 0.3.0 (2022-05-30)

## Main functionality corrections

* The `summary` method has been added for `Xfgpm` and `fgpm` objects.

* New classes `"predict.Xfgpm"` and `"simulate.Xfgpm"`.

* All plotters are now `plot` methods for the classes `"Xfgpm"`, `"fgpm"`,
  `"predict.Xfgpm"` or `"simulate.Xfgpm"`.

* New functions `modelDef` and `[[` to refit a `fgpm` object from a `Xfgpm` one.

## Main flow changes

* `format4pred()` removed and replaced by `get_active_in()` as anticipated in 0.2.0 release.

## Minor flow changes

* Own `rm_between()` function implemented to remove dependency on the `qdapRegex` package.
  
* Minor change to resolve a `%dopar%` issue during parallelized processes based on `%dorng%`.



# funGp 0.2.2 (2021-07-21)

## Minor flow changes

* Minor change to resolve a `check` issue evoked by CRAN.



# funGp 0.2.1 (2020-11-24)

## Minor flow changes

* Minor change to resolve a `check` `NOTE` related to the declaration in `Imports` of the `plyr`
  package.



# funGp 0.2.0 (2020-11-17)

## Main functionality corrections

* Fixed incoherence between the sorting order of the list of explored models and their performance
  statistic, stored at the `xm@log.success@sols` and `xm@log.success@fitness` slots of the `Xfgpm`
  object delivered by `fgpm_factory()`. Now the models in `xm@log.success@sols` are correctly
  sorted in decreasing order by performance.
  
* Resolved bug causing crash in `fgpm_factory()` during structural optimization with only functional
  inputs, e.g., for the call `fgpm_factory(fIn = fIn, sOut = sOut)`. Thanks to Jérémy Rohmer for
  pointing out the bug.

* Resolved bug causing crash in `fgpm_factory()` during structural optimization with only scalar
  inputs, e.g., for the call `fgpm_factory(sIn = sIn, sOut = sOut)`.

## Minor display changes

* Improved display of progress bars for `fgpm()` and `fgpm_factory()`.

* Improved placing of the legend in `plotX()` to prevent hiding important information.

## Minor flow changes

* `format4pred()` deprecated for `get_active_in()`.



# funGp 0.1.0 (2020-04-22)

* This is the first release of funGp.
