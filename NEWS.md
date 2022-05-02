# funGp 0.3.0 (2022-05-02)

## Main functionality corrections

* `Summary` method added for `Xfgpm`

* New classes `predict.Xfgpm` and `simulate.Xfgpm`

* All plotters are `plot` methods

* New functions `modelDef` and `[[` to refit a `fgpm` object from a `Xfgpm` one

## Main flow changes

* `format4pred()` removed and replaced by `get_active_in()` as  anticipated in 0.2.0 release

## Minor flow changes

* Own `rm_between()` and function implemented to remove dependency on the `qdapRegex` package
  already scheduled by CRAN for archival
  
* Minor change to resolve a `%dopar%` issue during parallelized processes based on `%dorng%`



# funGp 0.2.2 (2021-07-21)

## Minor flow changes

* Minor change to resolve a `check` issue evoked by CRAN



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

* Improved display of progress bars for `fgpm()` and `fgpm_factory()`

* Improved placing of the legend in `plotX()` to prevent hiding important information

## Minor flow changes

* `format4pred()` deprecated for `get_active_in()`



# funGp 0.1.0 (2020-04-22)

* This is the first release of funGp.
