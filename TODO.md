# TODO or not TODO

-   **(x)** is used for tasks that are completed,
-   **(WIP)** is used for tasks that are in progress,
-   **(?)** Idea to be discussed.

## Checks

-   **(x)** Check that the dimensions/the content of the data given as
    inputs for the `predict` and the `simulate` methods are consistent
    with the content of the object used to predict or simulate from.

## `summary` methods

-   **(x)** Add a `summary` method for the class `"fgpm"`.

-   **(x)** Add a `summary` method for the class `"Xfgpm"`.

## `plot` methods

-   **(x)** Add a `plot` method for the class `"fgpm"`.

-   **(x)** Add a `plot` method for the class `"Xfgpm"`.

-   **(WIP)** For the plot methods that produce two plots on a same
    graphic, add an extra argument `horiz` to allow the user to chose
    the the way in which the plots are set. This will be useful for
    reports or articles such as the JSS.

-   Add an `horiz` argument in the `plot`method for the S3 class
    `"predict.fgpm"`and in the `plotPred` method/function.

-   **(?)** Add `autoplot` methods allowing to get \*\*ggplot2\*
    graphics.

## Class `"Xfgpm"`

-   **(x)** Store a copy of `sIn`, `fIn` and `sOut` as three new slots
    of the class. This will help at recreating any `fgpm` model visited
    from the data used in the fit. Maybe write a new function for that
    aim? How to specify the model to be re-created: by its number?

-   **(x)** Improve the `summary` method. Should work correctly even
    when there are many inputs. Maybe add options to provide details on
    the scalar and the functional inputs separately in order to avoid
    having a table which is too wide.

-   **(x)** When an object is created by using `fgpm_factory`, the
    progression bars whould be displayed only in an interactive session
    i.e., where `interactive()` is `TRUE`.

## S3 classes `"predict.fgpm"` and `"simulate.fgpm"`

-   **(x)** Set the S3 class ot the output of the `predict` and
    `simulate` methods to `"predict.fgpm"` and `"simulate.fgpm"` both
    inheriting from `"list"`. This allows to write (S3) plot methods
    without breaking the inherited standard behaviours (in `print`, …).

## Simplifications

-   **(?)** Should we still export the methods `plotSims`, `plotPred`, …
    ? As a compromise we could continue to export the functions, but no
    longer turn them into S4 generics and methods because e.g.,
    `plotEvol` is unlikely to be mnemonic for the user or to be used in
    some other package.

-   **(?)** Rename the formal arguments of the `predict` and `simulate`
    methods for the class `"fgpm"` to use the same formal names? Could
    be for instance `sIn.new` and `fIn.new` or even simply `sIn`and
    `fIn`, since it is clear that this must be “new” values.

-   **(?)** Rename `fgpm_factory` into `Xfgpm` to conform to the R rule:
    the creator of a class has the same name as the class itself.

## Tests

-   **(?)** Add tests to make sure that the checks on the size of
    objects actually do what we expect.

## New functions

-   Add a function `model_def` (or other name) that takes `Xfgpm` object
    and returns a parsed code for one or several model. Bys using `eval`
    on this It could also return a list of model definitions. Then
    `eval(model_def[[4]])` would rebuild the 4-th model. Using `fIn` and
    `sIn` in the current env?

## Package documentation

-   Read again carefully the doc for the new functions and methods. Find
    the `"XXXY"` code that shows something that still requires work.

-   In the examples the output provided as comments should be removed,
    because this is an unusual practice in R. This behaviour is for
    vignettes. Some examples can be made simpler.

-   In the examples the `plot` methods should be used in place of
    `plotLOO`, `plotPreds`, … Similarly `summary` can be useful.

## modelDef

-   Find out why the class `"modelCall"` (with no method other than
    `show` exists. Maybe then rename `modelDef` into `modelCall`?

-   Make the **Examples** section faster to run (using `dontrun`?) and
    improve the example which leads to a much smaller Q2 when evaluated
    on the new data.

## Misc

-   No longer export `formatShort`.

## NEWS

Try

    library(funGp)
    example(plot.predict.fgpm)

    ## 
    ## plt.p.> # plotting predictions without the true output values_______________________
    ## plt.p.> # building the model
    ## plt.p.> set.seed(100)
    ## 
    ## plt.p.> n.tr <- 25
    ## 
    ## plt.p.> sIn <- expand.grid(x1 = seq(0, 1, length = sqrt(n.tr)),
    ## plt.p.+                    x2 = seq(0, 1, length = sqrt(n.tr)))
    ## 
    ## plt.p.> fIn <- list(f1 = matrix(runif(n.tr * 10), ncol = 10),
    ## plt.p.+             f2 = matrix(runif(n.tr * 22), ncol = 22))
    ## 
    ## plt.p.> sOut <- fgp_BB3(sIn, fIn, n.tr)
    ## 
    ## plt.p.> m1 <- fgpm(sIn = sIn, fIn = fIn, sOut = sOut)

    ## ** Presampling...

    ## ** Optimising hyperparameters...

    ## final  value 2.841058 
    ## converged

    ## ** Hyperparameters done!

    ## 
    ## plt.p.> # making predictions
    ## plt.p.> n.pr <- 100
    ## 
    ## plt.p.> sIn.pr <- as.matrix(expand.grid(x1 = seq(0,1,length = sqrt(n.pr)),
    ## plt.p.+                                 x2 = seq(0,1,length = sqrt(n.pr))))
    ## 
    ## plt.p.> fIn.pr <- list(f1 = matrix(runif(n.pr * 10), ncol = 10),
    ## plt.p.+                f2 = matrix(runif(n.pr * 22), ncol = 22))
    ## 
    ## plt.p.> m1.preds <- predict(m1, sIn.pr = sIn.pr, fIn.pr = fIn.pr)
    ## 
    ## plt.p.> # plotting predictions
    ## plt.p.> plotPreds(m1, preds = m1.preds)

    ## 
    ## plt.p.> # plotting predictions and true output values_______________________________
    ## plt.p.> # building the model
    ## plt.p.> set.seed(100)
    ## 
    ## plt.p.> n.tr <- 25
    ## 
    ## plt.p.> sIn <- expand.grid(x1 = seq(0, 1, length = sqrt(n.tr)),
    ## plt.p.+                    x2 = seq(0, 1, length = sqrt(n.tr)))
    ## 
    ## plt.p.> fIn <- list(f1 = matrix(runif(n.tr * 10), ncol = 10),
    ## plt.p.+             f2 = matrix(runif(n.tr * 22), ncol = 22))
    ## 
    ## plt.p.> sOut <- fgp_BB3(sIn, fIn, n.tr)
    ## 
    ## plt.p.> m1 <- fgpm(sIn = sIn, fIn = fIn, sOut = sOut)

    ## ** Presampling...
    ## ** Optimising hyperparameters...

    ## final  value 2.841058 
    ## converged

    ## ** Hyperparameters done!

![](TODO_files/figure-markdown_strict/unnamed-chunk-1-1.png)

    ## 
    ## plt.p.> # making predictions
    ## plt.p.> n.pr <- 100
    ## 
    ## plt.p.> sIn.pr <- as.matrix(expand.grid(x1 = seq(0,1,length = sqrt(n.pr)),
    ## plt.p.+                                 x2 = seq(0,1,length = sqrt(n.pr))))
    ## 
    ## plt.p.> fIn.pr <- list(f1 = matrix(runif(n.pr*10), ncol = 10),
    ## plt.p.+                f2 = matrix(runif(n.pr*22), ncol = 22))
    ## 
    ## plt.p.> m1.preds <- predict(m1, sIn.pr = sIn.pr, fIn.pr = fIn.pr)
    ## 
    ## plt.p.> # generating output data for validation
    ## plt.p.> sOut.pr <- fgp_BB3(sIn.pr, fIn.pr, n.pr)
    ## 
    ## plt.p.> # plotting predictions. Note that the 2-nd argument is 'y'
    ## plt.p.> plot(m1.preds, sOut.pr)

![](TODO_files/figure-markdown_strict/unnamed-chunk-1-2.png)

    ## 
    ## plt.p.> # only calibration plot
    ## plt.p.> plot(m1.preds, sOut.pr = sOut.pr, sortp = FALSE)

![](TODO_files/figure-markdown_strict/unnamed-chunk-1-3.png)

    ## 
    ## plt.p.> # only sorted output plot
    ## plt.p.> plot(m1.preds, sOut.pr = sOut.pr, calib = FALSE)

![](TODO_files/figure-markdown_strict/unnamed-chunk-1-4.png)

    example(plt.simulate.fgpm)

    ## Warning in example(plt.simulate.fgpm): no help found for 'plt.simulate.fgpm'
    ##FB test git 24/02/2022