# funGp 0.1.0 (2020-04-22)

* This is the first release of funGp.

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

* Iproved display of progress bars for `fgpm()` and `fgpm_factory()`

* Improved placing of the legend in `plotX()` to prevent hiding important information

## Minor flow changes

* `format4pred()` deprecated for `get_active_in()`
