## ******************************************************************************
## 'New' plot methods. The are simply standard interfaces to the
## existing methods 'plotLOO', 'plotEvol", 'plotX', 'plotPreds' and
## 'plotSim', but the user can access to these by following in a
## standard R route.
##
## ******************************************************************************

## ==============================================================================
## plot generic. Left undocumented by omitting the title
## ==============================================================================
##' @name plot
##' @rdname plot-methods
##' @exportMethod plot
setGeneric(name = "plot", def = function(x, y, ...) standardGeneric("plot"))

## ==============================================================================
## plot method, class "fgpm"
## ==============================================================================
##' Plot the Leave-One-Out (LOO) calibration.
##'
##' @description This method provides a diagnostic plot for the
##'     validation of regression models. It displays a calibration plot
##'     based on the leave-one-out predictions of the output at the
##'     points used to train the model.
##'
##' @title Plot method for the class \code{"fgpm"}
##'
##' @param x A \code{fgpm} object.
##' @param y Not used.
##' @param ... Graphical parameters. These currently include
##' \itemize{
##'    \item{\code{xlim}, \code{ylim}}{ to set the limits of the axes.}
##'    \item{\code{pch}, \code{pt.col}, \code{pt.bg}, \code{pt.cex}}{ to set
##'        the symbol used for the points and the related properties.}
##'    \item{\code{line}}{ to set the color used for the line.}
##'    \item{\code{xlab}, \code{ylab}, \code{main}}{ to set
##'        the labels of the axes and the main title.} See
##'        \bold{Examples}.
##' }
##' @return Nothing.
##' @export
##' @method plot fgpm
##' @examples
##' # generating input and output data for training
##' set.seed(100)
##' n.tr <- 25
##' sIn <- expand.grid(x1 = seq(0,1,length = sqrt(n.tr)),
##'                    x2 = seq(0, 1, length = sqrt(n.tr)))
##' fIn <- list(f1 = matrix(runif(n.tr*10), ncol = 10),
##'             f2 = matrix(runif(n.tr*22), ncol = 22))
##' sOut <- fgp_BB3(sIn, fIn, n.tr)
##'
##' # building the model
##' m1 <- fgpm(sIn = sIn, fIn = fIn, sOut = sOut)
##'
##' # plotting the model
##' plot(m1)
##' # change some graphical parameters if wanted
##' plot(m1, line = "SpringGreen3" ,
##'      pch = 21, pt.col = "orangered", pt.bg = "gold",
##'      main = "LOO cross-validation")
##'
setMethod("plot", "fgpm",
          function(x, y = NULL, ...) {
               plotLOO(model = x, ...)
          })

## ==============================================================================
## plot method, class "Xfgpm"
## ==============================================================================
##' @description Plot an object with class \code{"Xfgpm"} representing
##'     a collection of functional GP models corresponding to
##'     different structural parameters.
##'
##' Two types of graphics can be shown depending on the choice of
##' \code{which}. The choice \code{which = "evolution"} is used to
##' assess the quality of the fitted \code{fgpm} models on the basis
##' of Leave-One-Out cross-validation.  The choice \code{which =
##' "diag"} is used to display diagnostics. Two types of diagnostic
##' plots are shown as sub-plots by default, but each can be discarded
##' if wanted.
##'
##' The choice \code{which = "evolution"} displays the evolution of the
##' quality of the configurations evaluated along the iterations, by the
##' model selection algorithm in the fgpm_factory function. For each
##' iteration, the performance statistic of all the evaluated models is
##' printed, along with the corresponding median of the group. The plot
##' also includes the global maximum, which corresponds to the best
##' performance statistic obtained up to the current iteration. In this
##' plot, it is typical to have some points falling relatively far from
##' the maximum, even after multiple iterations. This happens mainly
##' because we have multiple categorical features, whose alteration
##' might change the performance statistic in a nonsmooth way. On the
##' other hand, the points that fall bellow zero usually correspond to
##' models whose hyperparameters were hard to optimize. This occurs
##' sporadically during the log-likelihood optimization for Gaussian
##' processes, due to the non-linearity of the objective function.
##' As long as the maximum keeps improving and the median remains close
##' to it, none of the two aforementioned phenomena is matter for worries.
##' Both of them respond to the mechanism of exploration implemented in
##' the algorithm, which makes it able to progressively move towards better
##' model configurations.
##'
##'  The choice \code{which = "diag"} provides two plots for assessing the quality of the output delivered by the
##'  model selection algorithm in the \link[funGp]{fgpm_factory} function. The first one is a calibration
##'  plot similar to the one offered for \linkS4class{fgpm} objects by \link[funGp]{plot,fgpm-method}.
##'  This plot allows to validate the absolute quality of the selected model. The second one displays the
##'  performance statistic of all the models successfully evaluated by the model selection algorithm. This
##'  provides a notion of the relative quality of the selected model with respect to the other models that
##'  can be made using the same data.
##'
##' @title Plot method for the class \code{"Xfgpm"}
##'
##' @param x The \code{Xfgpm} object to plot.
##'
##' @param y Not used.
##'
##' @param which Character giving the type of plot wanted. Can take the value
##'     \code{"evol"} or \code{"diag"}. See \bold{Examples}.
##'
##' @param calib Logical. If \code{TRUE} the calibration plot of the
##'     selected model will be included in the display in its
##'     "diagnostic" part if \code{which} is set to \code{"diag"}.
##'
##' @param fitp Logical. If \code{TRUE} a scatter plot of the quality
##'     of all explored models will be included in the display in its
##'     "diagnostic" part if \code{which} is set to \code{"diag"}.
##'
##' @param horiz Logical. Used only when \code{which} is \code{"diag"}
##'     and when both \code{calib} and \code{fitp} are \code{TRUE}. If
##'     \code{horiz} is \code{TRUE} the two subplots are displayed
##'     horizontally (on a same row) rather than vertically which is
##'     the default.
##'
##' @param ... Other graphical parameters such as \code{main} of
##'     \code{xlab}. When \code{which} is \code{"diag"} and both
##'     \code{calib} and \code{fitp} are \code{TRUE}, the graphical parameters
##'     should be enclosed into a list and passed with the formal name
##'     \code{calib.gpars} or \code{fitp.gpars}.
##'
##' @return Nothing.
##'
##'
##' @export
##' @method plot Xfgpm
##'
##' @seealso \strong{*} \link[funGp]{fgpm_factory} for structural optimization of funGp models.
##'
##' @examples
##' # generating input and output data
##' set.seed(100)
##' n.tr <- 2^5
##' x1 <- x2 <- x3 <- x4 <- x5 <- seq(0, 1, length = n.tr^(1/5))
##' sIn <- expand.grid(x1 = x1, x2 = x2, x3 = x3, x4 = x4, x5 = x5)
##' fIn <- list(f1 = matrix(runif(n.tr * 10), ncol = 10),
##'             f2 = matrix(runif(n.tr * 22), ncol = 22))
##' sOut <- fgp_BB7(sIn, fIn, n.tr)
##' \dontrun{
##' # optimizing the model structure with 'fgpm_factory' (~5 seconds)
##' xm <- fgpm_factory(sIn = sIn, fIn = fIn, sOut = sOut)
##' # assessing the quality of the model - absolute and w.r.t. the other
##' # explored models
##' plot(xm, which = "evol")
##' # Diagnostics (two subplots)
##' plot(xm, which = "diag")
##' plot(xm, which = "diag", horiz = TRUE)
##' # Diagnostics (one plot)
##' plot(xm, which = "diag", fitp = FALSE)
##' plot(xm, which = "diag", calib = FALSE)
##' # customizing some graphical parameters
##' plotX(xm, calib.gpars = list(xlim = c(800,1000), ylim = c(600,1200)),
##'      fitp.gpars = list(main = "Relative quality", legends = FALSE))
##'
##' }
setMethod("plot", "Xfgpm",
          function(x, y = NULL, which = c("evol", "diag"),
              calib = TRUE, fitp = TRUE, horiz = FALSE, ...) {
                  which <- match.arg(which)
                  if (which == "evol") {
                      if (!missing(calib) || ! missing(fitp)) {
                          warning("the formal arguments 'calib' and 'fitp' are ",
                                  "used only when 'which' is \"diag\"")
                      }
                      plotEvol.Xfgpm(x.model = x,  ...)
                  } else if (which == "diag") {
                      plotX.Xfgpm(x.model = x, calib = calib, fitp = fitp,
                                 horiz = horiz, ...)
                  }
              })

## ==============================================================================
## plot method, S3 class "predict.fgpm"
## ==============================================================================
##' @title Plot method for the predictions of a \code{fgpm} model
##'
##' @description This method displays the predicted output values
##'     delivered by a funGp Gaussian process model.
##'
##' @param x An object with S3 class \code{"predict.fgpm"}. This is a
##'     list containing the predictions and confidence bands as
##'     created by \link[funGp]{predict,fgpm-method} for the S3 class
##'     \code{"fgpm"}.
##'
##' @param y An optional vector (or 1-column matrix) containing
##'     the true values of the scalar output at the prediction
##'     points. If provided, the method will display two figures: (i) a
##'     calibration plot with true vs predicted output values, and (ii)
##'     a plot including the true and predicted output along with the
##'     confidence bands, sorted according to the increasing order of
##'     the true output. If not provided, only the second plot will be
##'     made, and the predictions will be arranged according to the
##'     increasing order of the predicted output.
##'
##' @param sOut.pr Alias of \code{y}, used for compatibility reasons.
##'
##' @param calib An optional boolean indicating if the calibration
##'     plot should be displayed. Ignored if \code{sOut.pr} is not
##'     provided. Default is \code{TRUE}.
##'
##' @param sortp An optional boolean indicating if the plot of sorted
##'     output should be displayed. Default is TRUE.
##'
##' @param ... Additional arguments affecting the display. Since this
##'     method allows to generate two plots from a single function
##'     call, the extra arguments for each plot should be included in a
##'     list. For the calibration plot, the list should be called
##'     \emph{calib.gpars}. For the plot of the output in increasing
##'     order, the list should be called \emph{sortp.gpars}. The
##'     following typical graphics parameters are valid entries of both
##'     lists: \emph{xlim}, \emph{ylim}, \emph{xlab}, \emph{ylab},
##'     \emph{main}. The boolean argument \emph{legends} can also be
##'     included in any of the two lists in order to control the
##'     display of legends in the corresponding plot.
##'
##' @export
##' @method plot predict.fgpm
##'
##' @return None.
##'
##' @author José Betancourt, François Bachoc and Thierry Klein
##'
##' @seealso \strong{*} \link[funGp]{fgpm} for the construction of funGp models;
##' @seealso \strong{*} \link[funGp]{plot,fgpm-method} for model diagnostic plots;
##' @seealso \strong{*} \link[funGp]{simulate} for simulations based on a funGp model;
##' @seealso \strong{*} \link[funGp]{plot.simulate.fgpm} for simulation plots.
##'
##' @examples
##' # plotting predictions without the true output values_______________________
##' # building the model
##' set.seed(100)
##' n.tr <- 25
##' sIn <- expand.grid(x1 = seq(0, 1, length = sqrt(n.tr)),
##'                    x2 = seq(0, 1, length = sqrt(n.tr)))
##' fIn <- list(f1 = matrix(runif(n.tr * 10), ncol = 10),
##'             f2 = matrix(runif(n.tr * 22), ncol = 22))
##' sOut <- fgp_BB3(sIn, fIn, n.tr)
##' m1 <- fgpm(sIn = sIn, fIn = fIn, sOut = sOut)
##'
##' # making predictions
##' n.pr <- 100
##' sIn.pr <- as.matrix(expand.grid(x1 = seq(0,1,length = sqrt(n.pr)),
##'                                 x2 = seq(0,1,length = sqrt(n.pr))))
##' fIn.pr <- list(f1 = matrix(runif(n.pr * 10), ncol = 10),
##'                f2 = matrix(runif(n.pr * 22), ncol = 22))
##' m1.preds <- predict(m1, sIn.pr = sIn.pr, fIn.pr = fIn.pr)
##'
##' # plotting predictions
##' plot(m1.preds)
##'
##' # plotting predictions and true output values_______________________________
##' # building the model
##' set.seed(100)
##' n.tr <- 25
##' sIn <- expand.grid(x1 = seq(0, 1, length = sqrt(n.tr)),
##'                    x2 = seq(0, 1, length = sqrt(n.tr)))
##' fIn <- list(f1 = matrix(runif(n.tr * 10), ncol = 10),
##'             f2 = matrix(runif(n.tr * 22), ncol = 22))
##' sOut <- fgp_BB3(sIn, fIn, n.tr)
##' m1 <- fgpm(sIn = sIn, fIn = fIn, sOut = sOut)
##'
##' # making predictions
##' n.pr <- 100
##' sIn.pr <- as.matrix(expand.grid(x1 = seq(0,1,length = sqrt(n.pr)),
##'                                 x2 = seq(0,1,length = sqrt(n.pr))))
##' fIn.pr <- list(f1 = matrix(runif(n.pr*10), ncol = 10),
##'                f2 = matrix(runif(n.pr*22), ncol = 22))
##' m1.preds <- predict(m1, sIn.pr = sIn.pr, fIn.pr = fIn.pr)
##'
##' # generating output data for validation
##' sOut.pr <- fgp_BB3(sIn.pr, fIn.pr, n.pr)
##'
##' # plotting predictions. Note that the 2-nd argument is 'y'
##' plot(m1.preds, sOut.pr)
##'
##' # only calibration plot
##' plot(m1.preds, sOut.pr = sOut.pr, sortp = FALSE)
##'
##' # only sorted output plot
##' plot(m1.preds, sOut.pr = sOut.pr, calib = FALSE)
##'
plot.predict.fgpm <- function(x, y = NULL, sOut.pr = NULL,
                              calib = TRUE, sortp = TRUE, ...) {

    if (!is.null(y)) {
        if (!is.null(sOut.pr)) {
            warning("both 'y' and 'sOut.pr' are given. Only 'y' ",
                    "will be used")
        }
        sOut.pr <- y
    }

    plotPreds.fgpm(preds = x, sOut.pr = sOut.pr, calib = calib,
                   sortp = sortp, ...)

}

# ==============================================================================
## plot method, S3 class "simulate.fgpm"
## ==============================================================================
##' @title Plot method for the simulations of a \code{fgpm} model
##'
##' @description This method displays the simulated output values
##'     delivered by a funGp Gaussian process model.
##'
##' @param x An object with S3 class \code{simulate.fgpm} as created
##'     by \link[funGp]{simulate,fgpm-method}.
##'
##' @param y Not used.
##'
##' @param detail An optional character string specifying the data
##'     elements that should be included in the plot, to be chosen
##'     between \code{"light"} and \code{"full"}. A \emph{light} plot
##'     will include only the simulated values, while a a
##'     \emph{full} plot will also include the predicted mean and
##'     confidence bands at the simulation points. This argument will
##'     only be used if full simulations (including the mean and
##'     confidence bands) are provided, otherwise it will be
##'     ignored. See \link[funGp]{simulate} for more details on the
##'     generation of light and full simulations.
##'
##' @param ... Additional arguments affecting the display. The
##'     following typical graphics parameters are valid entries:
##'     \emph{xlim}, \emph{ylim}, \emph{xlab}, \emph{ylab},
##'     \emph{main}. The boolean argument \emph{legends} can also be
##'     included in any of the two lists in order to control the
##'     display of legends in the corresponding plot.
##'
##' @return None.
##'
##' @author José Betancourt, François Bachoc and Thierry Klein
##'
##' @export
##' @method plot simulate.fgpm
##'
##' @examples
##' # plotting light simulations________________________________________________
##' # building the model
##' set.seed(100)
##' n.tr <- 25
##' sIn <- expand.grid(x1 = seq(0, 1, length = sqrt(n.tr)),
##'                    x2 = seq(0, 1, length = sqrt(n.tr)))
##' fIn <- list(f1 = matrix(runif(n.tr * 10), ncol = 10),
##'             f2 = matrix(runif(n.tr * 22), ncol = 22))
##' sOut <- fgp_BB3(sIn, fIn, n.tr)
##' m1 <- fgpm(sIn = sIn, fIn = fIn, sOut = sOut)
##'
##' # making light simulations
##' n.sm <- 100
##' sIn.sm <- as.matrix(expand.grid(x1 = seq(0, 1, length = sqrt(n.sm)),
##'                                 x2 = seq(0, 1, length = sqrt(n.sm))))
##' fIn.sm <- list(f1 = matrix(runif(n.sm * 10), ncol = 10),
##'                f2 = matrix(runif(n.sm * 22), ncol = 22))
##' simsl <- simulate(m1, nsim = 10, sIn.sm = sIn.sm, fIn.sm = fIn.sm)
##'
##' # plotting light simulations
##' plot(simsl)
##'
##'
##' # plotting full simulations_________________________________________________
##' # building the model
##' set.seed(100)
##' n.tr <- 25
##' sIn <- expand.grid(x1 = seq(0, 1, length = sqrt(n.tr)),
##'                    x2 = seq(0, 1, length = sqrt(n.tr)))
##' fIn <- list(f1 = matrix(runif(n.tr * 10), ncol = 10),
##'             f2 = matrix(runif(n.tr * 22), ncol = 22))
##' sOut <- fgp_BB3(sIn, fIn, n.tr)
##' m1 <- fgpm(sIn = sIn, fIn = fIn, sOut = sOut)
##'
##' # making full simulations
##' n.sm <- 100
##' sIn.sm <- as.matrix(expand.grid(x1 = seq(0, 1, length = sqrt(n.sm)),
##'                                 x2 = seq(0, 1 ,length = sqrt(n.sm))))
##' fIn.sm <- list(f1 = matrix(runif(n.sm * 10), ncol = 10),
##'                f2 = matrix(runif(n.sm * 22), ncol = 22))
##' simsf <- simulate(m1, nsim = 10, sIn.sm = sIn.sm, fIn.sm = fIn.sm,
##'                   detail = "full")
##'
##' # plotting full simulations in "full" mode
##' plot(simsf)
##'
##' # plotting full simulations in "light" mode
##' plot(simsf, detail = "light")
##'
plot.simulate.fgpm <- function(x, y = NULL, detail = NA, ...) {

    if (is.na(detail)) {
        if (is.list(x)) detail <- "full"
        else detail <- "light"
    } else {
        detVal <- c("light", "full")
        if (is.na(i <- pmatch(detail, detVal))) {
            stop("invalid value for 'detail'")
        } else detail <- detVal[i]
        if ((detail == "full") && !is.list(x)) {
            warning("the value \"full\" for 'detail' can only ",
                    "be used when `x` is a list")
            detail <- "light"
        }
    }
    plotSims.fgpm(sims = x, detail = detail, ...)

}
