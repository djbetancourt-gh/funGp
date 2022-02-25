# ==========================================================================================================
# Diagnostic calibration plot for funGp models
# ==========================================================================================================
#' @title Leave-one-out calibration plot for regression models
#' @description This method provides a diagnostic plot for the validation of regression models. It displays
#'   a calibration plot based on the leave-one-out predictions of the output at the points used to train the
#'   model.
#'
#' @param model a model object for which the LOO calibration plot is to be made.
#' @param ... additional arguments affecting the plot.
#'
#' @return None.
#'
#' @seealso \strong{*} \link[funGp]{plotLOO} for the diagnostic plot of a funGp model.
#'
#' @examples
#' require(funGp) # a package with a plotLOO method implemented
#'
#' # generating input and output data for training
#' set.seed(100)
#' n.tr <- 25
#' sIn <- expand.grid(x1 = seq(0,1,length = sqrt(n.tr)), x2 = seq(0,1,length = sqrt(n.tr)))
#' fIn <- list(f1 = matrix(runif(n.tr*10), ncol = 10), f2 = matrix(runif(n.tr*22), ncol = 22))
#' sOut <- fgp_BB3(sIn, fIn, n.tr)
#'
#' # building the model
#' m1 <- fgpm(sIn = sIn, fIn = fIn, sOut = sOut)
#'
#' # plotting the model
#' plotLOO(m1)
#'
#' @name plotLOO
#' @export
#' @keywords internal
setGeneric("plotLOO", function(model, ...) standardGeneric("plotLOO"))

#' @title Leave-one-out calibration plot for a funGp model
#' @description This method provides a diagnostic plot for the validation of a funGp Gaussian process model.
#'   It displays a calibration plot based on the leave-one-out predictions of the output at the points
#'   used to train the model.
#'
#' @param model an object of class \linkS4class{fgpm} corresponding to the funGp model to validate.
#' @param ... additional arguments affecting the plot. The following typical graphics parameters are
#'   valid entries: \emph{xlim}, \emph{ylim}, \emph{xlab}, \emph{ylab}, \emph{main}.
#'
#' @return None.
#'
#' @author José Betancourt, François Bachoc and Thierry Klein
#'
#' @seealso \strong{*} \link[funGp]{fgpm} for the construction of funGp models;
#' @seealso \strong{*} \link[funGp]{plotPreds} for prediction plots;
#' @seealso \strong{*} \link[funGp]{plotSims} for simulation plots.
#'
#' @examples
#' # generating input and output data for training
#' set.seed(100)
#' n.tr <- 25
#' sIn <- expand.grid(x1 = seq(0,1,length = sqrt(n.tr)), x2 = seq(0,1,length = sqrt(n.tr)))
#' fIn <- list(f1 = matrix(runif(n.tr*10), ncol = 10), f2 = matrix(runif(n.tr*22), ncol = 22))
#' sOut <- fgp_BB3(sIn, fIn, n.tr)
#'
#' # building the model
#' m1 <- fgpm(sIn = sIn, fIn = fIn, sOut = sOut)
#'
#' # plotting the model
#' plotLOO(m1)
#'
#' @rdname plotLOO-method
setMethod("plotLOO", "fgpm", function(model, ...) {
  plotLOO.fgpm(model = model, ...)
})

plotLOO.fgpm <- function(model, ...) {
  # recover observed output
  y_obs <- model@sOut

  # compute loocv predictions
  R <- tcrossprod(model@preMats$L)/model@kern@varHyp + diag(model@nugget, nrow = model@n.tr, ncol = model@n.tr)
  Rinv <- solve(R)
  y_pre <- y_obs - diag(Rinv)^(-1) * Rinv %*% y_obs

  # compute LOO statistic
  q2 <- format(getFitness(model), digits = 3, nsmall = 3)

  # recover graphic parameters if provided
  gpars <- list(...)
  if (!is.null(gpars$xlim)) xlim <- gpars$xlim else xlim <- range(c(y_obs, y_pre))
  if (!is.null(gpars$ylim)) ylim <- gpars$ylim else ylim <- range(c(y_obs, y_pre))
  if (!is.null(gpars$pch)) pch <- gpars$pch else pch <- 21
  if (!is.null(gpars$pt.col)) pt.col <- gpars$pt.col else pt.col <- "red"
  if (!is.null(gpars$pt.bg)) pt.bg <- gpars$pt.bg else pt.bg <- "red"
  if (!is.null(gpars$pt.cex)) pt.cex <- gpars$pt.cex else pt.cex <- 1
  if (!is.null(gpars$line)) line <- gpars$line else line <- "blue"
  if (!is.null(gpars$xlab)) xlab <- gpars$xlab else xlab <- "Observed"
  if (!is.null(gpars$ylab)) ylab <- gpars$ylab else ylab <- "Predicted"
  if (!is.null(gpars$main)) main <- gpars$main else main <- "Model diagnostic by leave-one-out cross-validation"

  # save current par state
  opar <- par('mar', 'mfrow')
  on.exit(par(opar))

  # set up layout
  par(mar = c(5.1, 4.1, 4.1, 2.1), mfrow = c(1,1))

  # plot
  plot(y_obs, y_pre, xlim = xlim, ylim = ylim, pch = pch, col = pt.col, bg = pt.bg, cex = pt.cex,
       main = main, xlab = xlab, ylab = ylab)
  lnlims <- range(c(xlim, ylim))
  lines(lnlims, lnlims, col = line)
  legend("topleft", legend = paste("Q2loocv =", q2),
         xjust = 0.5, yjust = 0.5, x.intersp = -0.5, y.intersp = 0.3,
         adj = c(0, 0.5), inset = c(.02,.05))
}
# ==========================================================================================================



# ==========================================================================================================
# Plot of predictions in increasing order with confidence bands
# ==========================================================================================================
#' @title Plot for predictions of regression models
#' @description This method displays the predicted output values delivered by some regression model. The
#'   plot might be constituted differently, depending on the type of model at hand.
#'
#' @param model a model object for which the plot is to be made.
#' @param ... additional arguments affecting the plot.
#'
#' @return None.
#'
#' @seealso \strong{*} \link[funGp]{plotPreds} for the predictions plot of a funGp model.
#'
#' @examples
#' require(funGp) # a package with a plotPreds method implemented
#' # building the model
#' set.seed(100)
#' n.tr <- 25
#' sIn <- expand.grid(x1 = seq(0,1,length = sqrt(n.tr)), x2 = seq(0,1,length = sqrt(n.tr)))
#' fIn <- list(f1 = matrix(runif(n.tr*10), ncol = 10), f2 = matrix(runif(n.tr*22), ncol = 22))
#' sOut <- fgp_BB3(sIn, fIn, n.tr)
#' m1 <- fgpm(sIn = sIn, fIn = fIn, sOut = sOut)
#'
#' # making predictions
#' n.pr <- 100
#' sIn.pr <- as.matrix(expand.grid(x1 = seq(0,1,length = sqrt(n.pr)),
#'                                 x2 = seq(0,1,length = sqrt(n.pr))))
#' fIn.pr <- list(f1 = matrix(runif(n.pr*10), ncol = 10), matrix(runif(n.pr*22), ncol = 22))
#' m1.preds <- predict(m1, sIn.pr = sIn.pr, fIn.pr = fIn.pr)
#'
#' # plotting predictions
#' plotPreds(m1, preds = m1.preds)
#'
#' @name plotPreds
#' @export
#' @keywords internal
setGeneric("plotPreds", function(model, ...) standardGeneric("plotPreds"))

#' @title Plot for predictions of a funGp model
#' @description This method displays the predicted output values delivered by a funGp Gaussian process model.
#'
#' @param model a \linkS4class{fgpm} object for which the plot is to be made.
#' @param preds a list containing the predictions and confidence bands. In funGp, this argument is just the
#'   data structure delivered by the \link[funGp]{predict} method.
#' @param sOut.pr an optional vector (or 1-column matrix) containing the true values of the scalar output at
#'   the prediction points. If provided, the method will display two figures: (i) a calibration plot with
#'   true vs predicted output values, and (ii) a plot including the true and predicted output along with the
#'   confidence bands, sorted according to the increasing order of the true output. If not provided, only
#'   the second plot will be made, and the predictions will be arranged according to the increasing order of
#'   the predicted output.
#' @param calib an optional boolean indicating if the calibration plot should be displayed. Ignored if sOut.pr
#'   is not provided. Default is TRUE.
#' @param sortp an optional boolean indicating if the plot of sorted output should be displayed. Default is
#'   TRUE.
#' @param ... additional arguments affecting the display. Since this method allows to generate two plots from
#'   a single function call, the extra arguments for each plot should be included in a list. For the
#'   calibration plot, the list should be called \emph{calib.gpars}. For the plot of the output in increasing
#'   order, the list should be called \emph{sortp.gpars}. The following typical graphics parameters are valid
#'   entries of both lists: \emph{xlim}, \emph{ylim}, \emph{xlab}, \emph{ylab}, \emph{main}. The boolean
#'   argument \emph{legends} can also be included in any of the two lists in order to control the display of
#'   legends in the corresponding plot.
#'
#' @return None.
#'
#' @author José Betancourt, François Bachoc and Thierry Klein
#'
#' @seealso \strong{*} \link[funGp]{fgpm} for the construction of funGp models;
#' @seealso \strong{*} \link[funGp]{plotLOO} for model diagnostic plots;
#' @seealso \strong{*} \link[funGp]{simulate} for simulations based on a funGp model;
#' @seealso \strong{*} \link[funGp]{plotSims} for simulation plots.
#'
#' @examples
#' # plotting predictions without the true output values______________________________________
#' # building the model
#' set.seed(100)
#' n.tr <- 25
#' sIn <- expand.grid(x1 = seq(0,1,length = sqrt(n.tr)), x2 = seq(0,1,length = sqrt(n.tr)))
#' fIn <- list(f1 = matrix(runif(n.tr*10), ncol = 10), f2 = matrix(runif(n.tr*22), ncol = 22))
#' sOut <- fgp_BB3(sIn, fIn, n.tr)
#' m1 <- fgpm(sIn = sIn, fIn = fIn, sOut = sOut)
#'
#' # making predictions
#' n.pr <- 100
#' sIn.pr <- as.matrix(expand.grid(x1 = seq(0,1,length = sqrt(n.pr)),
#'                                 x2 = seq(0,1,length = sqrt(n.pr))))
#' fIn.pr <- list(f1 = matrix(runif(n.pr*10), ncol = 10), matrix(runif(n.pr*22), ncol = 22))
#' m1.preds <- predict(m1, sIn.pr = sIn.pr, fIn.pr = fIn.pr)
#'
#' # plotting predictions
#' plotPreds(m1, preds = m1.preds)
#'
#'
#' # plotting predictions and true output values______________________________________________
#' # building the model
#' set.seed(100)
#' n.tr <- 25
#' sIn <- expand.grid(x1 = seq(0,1,length = sqrt(n.tr)), x2 = seq(0,1,length = sqrt(n.tr)))
#' fIn <- list(f1 = matrix(runif(n.tr*10), ncol = 10), f2 = matrix(runif(n.tr*22), ncol = 22))
#' sOut <- fgp_BB3(sIn, fIn, n.tr)
#' m1 <- fgpm(sIn = sIn, fIn = fIn, sOut = sOut)
#'
#' # making predictions
#' n.pr <- 100
#' sIn.pr <- as.matrix(expand.grid(x1 = seq(0,1,length = sqrt(n.pr)),
#'                                 x2 = seq(0,1,length = sqrt(n.pr))))
#' fIn.pr <- list(f1 = matrix(runif(n.pr*10), ncol = 10), matrix(runif(n.pr*22), ncol = 22))
#' m1.preds <- predict(m1, sIn.pr = sIn.pr, fIn.pr = fIn.pr)
#'
#' # generating output data for validation
#' sOut.pr <- fgp_BB3(sIn.pr, fIn.pr, n.pr)
#'
#' # plotting predictions
#' plotPreds(m1, m1.preds, sOut.pr)
#'
#' # only calibration plot
#' plotPreds(m1, m1.preds, sOut.pr, sortp = FALSE)
#'
#' # only sorted output plot
#' plotPreds(m1, m1.preds, sOut.pr, calib = FALSE)
#'
#' @importFrom graphics lines plot polygon layout legend par mtext
#' @rdname plotPreds-method
setMethod("plotPreds", "fgpm",
          function(model, preds, sOut.pr = NULL, calib = TRUE, sortp = TRUE, ...) {
            plotPreds.fgpm(preds = preds, sOut.pr = sOut.pr,
                           calib = calib, sortp = sortp, ...)
          })

plotPreds.fgpm <- function(preds, sOut.pr, calib, sortp, ...) {
  # recover observed and predicted output
  if (!is.null(sOut.pr)) {
    ord <- order(sOut.pr)
    y_obs <- sOut.pr[ord]
  } else {
    ord <- order(preds$mean)
    y_obs <- NULL
  }
  y_pre <- preds$mean[ord]
  n.pr <- length(y_pre)

  # recover 95% confidence intervals
  ll <- preds$lower95[ord]
  ul <- preds$upper95[ord]

  # recover graphic parameters if provided
  gpars <- list(...)
  lin.gpars <- gpars$sortp.gpars
  cal.gpars <- gpars$calib.gpars

  # calibration plot _____________________________________________________
  plot.c <- function (tight = FALSE) {
    # <---> limits
    if (!is.null(cal.gpars$xlim)) xlim <- cal.gpars$xlim else xlim <- range(c(y_obs, y_pre))
    if (!is.null(cal.gpars$ylim)) ylim <- cal.gpars$ylim else ylim <- range(c(y_obs, y_pre))
    # <---> point colors
    if (!is.null(cal.gpars$pt.col)) pt.col <- cal.gpars$pt.col else pt.col <- "red"
    if (!is.null(cal.gpars$pt.bg)) pt.bg <- cal.gpars$pt.bg else pt.bg <- "red"
    # <---> point style
    if (!is.null(cal.gpars$pch)) pch <- cal.gpars$pch else pch <- 21
    if (!is.null(cal.gpars$pt.cex)) pt.cex <- cal.gpars$pt.cex else pt.cex <- 1
    # <---> line colors
    if (!is.null(cal.gpars$line)) line <- cal.gpars$line else line <- "blue"
    # <---> text
    if (!is.null(cal.gpars$xlab)) xlab <- cal.gpars$xlab else xlab <- "Observed"
    if (!is.null(cal.gpars$ylab)) ylab <- cal.gpars$ylab else ylab <- "Predicted"
    if (!is.null(cal.gpars$main)) main <- cal.gpars$main else main <- "Model predictions at new input points"
    # <---> legend
    if (!is.null(cal.gpars$legends)) legends <- cal.gpars$legends else legends <- T
    if (!is.null(cal.gpars$cex.leg)) cex.leg <- cal.gpars$cex.leg else cex.leg <- 1

    # add points
    if (tight) {
      plot(y_obs, y_pre, xlim = xlim, ylim = ylim, pch = pch, col = pt.col, bg = pt.bg, cex = pt.cex,
           main = main, xlab = "", ylab = ylab)
      mtext(side = 1, text = xlab, line = 2.0)
    } else {
      plot(y_obs, y_pre, xlim = xlim, ylim = ylim, pch = pch, col = pt.col, bg = pt.bg, cex = pt.cex,
           main = main, xlab = xlab, ylab = ylab)
    }

    # add reference line
    lines(xlim, ylim, col = line)

    if (legends) {
      # compute LOO statistic
      q2 <- format(1 - mean(sum((y_obs - y_pre)^2))/mean(sum((y_obs - mean(y_obs))^2)), digits = 3, nsmall = 3)

      # add hould-out Q2
      legend("topleft", legend = paste("Q2hout =", q2),
             xjust = 0.5, yjust = 0.5, x.intersp = -0.5, y.intersp = 0.3,
             adj = c(0, 0.5), inset = c(.02,.05), cex = cex.leg)
    }
  }
  # ______________________________________________________________________


  # line plot ____________________________________________________________
  plot.l <- function (tight = FALSE) {
    # <---> limits
    if (!is.null(lin.gpars$xlim)) xlim <- lin.gpars$xlim else xlim <- c(1, length(y_pre))
    if (!is.null(lin.gpars$ylim)) ylim <- lin.gpars$ylim else ylim <- range(c(ll, ul))
    # <---> line colors
    if (!is.null(lin.gpars$line.pr)) line.pr <- lin.gpars$line.pr else line.pr <- "red"
    if (!is.null(lin.gpars$line.ci)) line.ci <- lin.gpars$line.ci else line.ci <- "blue"
    if (!is.null(lin.gpars$line.true)) line.true <- lin.gpars$line.true else line.true <- "black"
    # <---> line styles
    if (!is.null(lin.gpars$lty.pr)) lty.pr <- lin.gpars$lty.pr else lty.pr <- 1
    if (!is.null(lin.gpars$lwd.pr)) lwd.pr <- lin.gpars$lwd.pr else lwd.pr <- 1
    if (!is.null(lin.gpars$lty.ci)) lty.ci <- lin.gpars$lty.ci else lty.ci <- 1
    if (!is.null(lin.gpars$lwd.ci)) lwd.ci <- lin.gpars$lwd.ci else lwd.ci <- 1
    if (!is.null(lin.gpars$lty.true)) lty.true <- lin.gpars$lty.true else lty.true <- 1
    if (!is.null(lin.gpars$lwd.true)) lwd.true <- lin.gpars$lwd.true else lwd.true <- 1
    # <---> polygon colors
    if (!is.null(lin.gpars$col.poly)) col.poly <- lin.gpars$col.poly else col.poly <- "grey85"
    # <---> text
    if (!is.null(lin.gpars$xlab)) xlab <- lin.gpars$xlab else xlab <- "Index"
    if (!is.null(lin.gpars$ylab)) ylab <- lin.gpars$ylab else ylab <- "Output"
    if (!is.null(lin.gpars$main)) main <- lin.gpars$main else main <- "Sorted predictions"
    # <---> legend
    if (!is.null(lin.gpars$legends)) legends <- lin.gpars$legends else legends <- T

    # set up frame
    if (tight) {
      plot(1, type = "n", xlim = xlim, ylim = ylim, main = main, xlab = "", xaxt = "n", ylab = ylab)
      mtext(side = 1, text = xlab, line = 2.0)
    } else {
      plot(1, type = "n", xlim = xlim, ylim = ylim, main = main, xlab = xlab, xaxt = "n", ylab = ylab)
    }
    axis(1, axtags(n.pr))

    # get device dimensions in inches to check arrows length
    units = par(c('usr', 'pin'))
    y_to_inches = with(units, pin[2L]/diff(usr[3:4]))

    # add confidence limits
    for (i in 1:n.pr) {
      # check arrow length in inches
      d = sqrt((y_to_inches * abs(ul[i]-ll[i]))^2)
      if (d >= .001) {
        arrows(x0 = i, x1 = i, y0 = ll[i], y1 = ul[i], length = .03, code = 3, angle = 90, lwd = 1.5, col = "blue")
      }
    }

    # add predicted mean
    points(y_pre, pch = 21, bg = "blue")

    # add true output if provided
    if (!is.null(y_obs)) {
      points(y_obs, pch = 21, bg = "green")
      if (legends) {
        legend("topleft", legend = c("True", "Pred. mean", "95% Conf."), pch = c(21, 21, NA),
               pt.bg = c("green", "blue", NA), inset = c(.02,.05), cex = .8)

        # arrow parameters in x
        limscr <- par('usr')[1:2]
        lcr <- limscr[1] + diff(limscr)*(0.02) # left coordinate of legend frame
        win <- par()$pin[1]*(1-0.02) # width of plot in inches (omitting left inset)
        pin <- 0.12 # location of point in inches
        lp <- pin/win # location as percentage of width
        wcr <- diff(limscr)*(1-0.02) # width of plot in coordinate system
        xcr <- lcr + wcr*lp # location of point in coordinate

        # arrow parameters in y
        limscr <- par('usr')[3:4]
        tcr <- limscr[2] - diff(limscr)*(0.05) # top coordinate of legend frame
        hin <- par()$pin[2]*(1-0.05) # height of plot in inches (omitting bottom inset)
        pin1 <- hin - 0.525 # end location of arrow in inches
        h <- .10 # arrow height in inches
        pin2 <- pin1+h # start location of arrow in inches
        lp1 <- 1-pin1/hin # location 1 as percentage of height
        lp2 <- 1-pin2/hin # location 2 as percentage of height
        hcr <- diff(limscr)*(1-0.05) # height of plot in coordinate system
        ycr1 <- tcr - hcr*lp1 # location 1 in coordinate
        ycr2 <- tcr - hcr*lp2 # location 1 in coordinate

        # arrow display
        arrows(x0 = xcr, x1 = xcr, y0 = ycr1, y1 = ycr2, length = 0.03, code = 3,
               angle = 90, lwd = 1.5, col = "blue")
      }

    } else {
      if (legends) {
        legend("topleft", legend = c("Pred. mean", "95% Conf."), pch = c(21, NA),
               pt.bg = c("blue", NA), inset = c(.02,.05))

        # arrow parameters in x
        limscr <- par('usr')[1:2]
        lcr <- limscr[1] + diff(limscr)*(0.02) # left coordinate of legend frame
        win <- par()$pin[1]*(1-0.02) # width of plot in inches (omitting left inset)
        pin <- 0.15 # location of point in inches
        lp <- pin/win # location as percentage of width
        wcr <- diff(limscr)*(1-0.02) # width of plot in coordinate system
        xcr <- lcr + wcr*lp # location of point in coordinate

        # arrow parameters in y
        limscr <- par('usr')[3:4]
        tcr <- limscr[2] - diff(limscr)*(0.05) # top coordinate of legend frame
        hin <- par()$pin[2]*(1-0.05) # height of plot in inches (omitting bottom inset)
        pin1 <- hin - 0.46 # end location of arrow in inches
        h <- .12 # arrow height in inches
        pin2 <- pin1+h # start location of arrow in inches
        lp1 <- 1-pin1/hin # location 1 as percentage of height
        lp2 <- 1-pin2/hin # location 2 as percentage of height
        hcr <- diff(limscr)*(1-0.05) # height of plot in coordinate system
        ycr1 <- tcr - hcr*lp1 # location 1 in coordinate
        ycr2 <- tcr - hcr*lp2 # location 1 in coordinate

        # arrow display
        arrows(x0 = xcr, x1 = xcr, y0 = ycr1, y1 = ycr2, length = 0.03, code = 3,
               angle = 90, lwd = 1.5, col = "blue")
      }
    }
  }
  # ______________________________________________________________________


  # save current par state
  opar <- par('mar', 'mfrow')
  on.exit(par(opar))

  # identify the case
  # case 1: plotPreds(model, preds) -> plot.l
  # case 2: plotPreds(model, preds, sOut.pr) -> [plot.c, plot.l]
  # case 3: plotPreds(model, preds, sOut.pr, calib = FALSE) -> [plot.l]
  if (is.null(sOut.pr)) {
    par(mar = c(5.1, 4.1, 4.1, 2.1), mfrow = c(1,1))
    plot.l()
  } else if (all(calib, sortp)) {
    par(mar = c(3.1, 4.1, 2.5, 2.1), mfrow = c(2,1))
    plot.c(T)
    plot.l(T)
  } else if (calib) {
    par(mar = c(5.1, 4.1, 4.1, 2.1), mfrow = c(1,1))
    plot.c()
  } else if (sortp) {
    par(mar = c(5.1, 4.1, 4.1, 2.1), mfrow = c(1,1))
    plot.l()
  }
}
# ==========================================================================================================



# ==========================================================================================================
# Plot of simulations by input index - option to add predicted mean and confidence bands
# ==========================================================================================================
#' @title Plot for simulations of random processes
#' @description This method displays the simulated output values delivered by some random process model.
#'   The plot might be constituted differently, depending on the type of model at hand.
#'
#' @param model a model object for which the plot is to be made.
#' @param sims data structure containing simulations Depending on the type of model and the data structure
#'   used, it might also contain, for instance, the mean and confidence bands at the simulation points.
#' @param ... additional arguments affecting the plot.
#'
#' @return None.
#'
#' @seealso \strong{*} \link[funGp]{plotSims} for the simulations plot of a funGp model.
#'
#' @examples
#' require(funGp) # a package with a plotSims method implemented
#' # building the model
#' set.seed(100)
#' n.tr <- 25
#' sIn <- expand.grid(x1 = seq(0,1,length = sqrt(n.tr)), x2 = seq(0,1,length = sqrt(n.tr)))
#' fIn <- list(f1 = matrix(runif(n.tr*10), ncol = 10), f2 = matrix(runif(n.tr*22), ncol = 22))
#' sOut <- fgp_BB3(sIn, fIn, n.tr)
#' m1 <- fgpm(sIn = sIn, fIn = fIn, sOut = sOut)
#'
#' # making simulations
#' n.sm <- 100
#' sIn.sm <- as.matrix(expand.grid(x1 = seq(0,1,length = sqrt(n.sm)),
#'                                 x2 = seq(0,1,length = sqrt(n.sm))))
#' fIn.sm <- list(f1 = matrix(runif(n.sm*10), ncol = 10), matrix(runif(n.sm*22), ncol = 22))
#' m1.sims <- simulate(m1, nsim = 10, sIn.sm = sIn.sm, fIn.sm = fIn.sm)
#'
#' # plotting simulations
#' plotSims(m1, m1.sims)
#'
#' @name plotSims
#' @export
#' @keywords internal
setGeneric("plotSims", function(model, sims, ...) standardGeneric("plotSims"))

#' @title Plot for simulations from a funGp model
#' @description This method displays the simulated output values delivered by a funGp Gaussian process model.
#'
#' @param model a \linkS4class{fgpm} object for which the plot is to be made.
#' @param sims a list containing the simulated output values. In funGp, this argument is just the data
#'   structure delivered by the \link[funGp]{simulate} method.
#' @param detail an optional character string specifying the data elements that should be included in the plot,
#'   to be chosen between "light" and "full". A \emph{light} plot will include only include the simulated
#'   values, while a a \emph{full} plot will also include the predicted mean and confidence bands at the
#'   simulation points. This argument will only be used if full simulations (including the mean and confidence
#'   bands) are provided, otherwise it will be dropped. See \link[funGp]{simulate} for more details on the
#'   generation of light and full simulations.
#' @param ... additional arguments affecting the display. The following typical graphics parameters are valid
#'   entries: \emph{xlim}, \emph{ylim}, \emph{xlab}, \emph{ylab}, \emph{main}. The boolean argument
#'   \emph{legends} can also be included in any of the two lists in order to control the display of legends
#'   in the corresponding plot.
#'
#' @return None.
#'
#' @author José Betancourt, François Bachoc and Thierry Klein
#'
#' @seealso \strong{*} \link[funGp]{fgpm} for the construction of funGp models;
#' @seealso \strong{*} \link[funGp]{plotLOO} for model diagnostic plots;
#' @seealso \strong{*} \link[funGp]{predict} for predictions based on a funGp model;
#' @seealso \strong{*} \link[funGp]{plotPreds} for prediction plots.
#'
#' @examples
#' # plotting light simulations_______________________________________________________________
#' # building the model
#' set.seed(100)
#' n.tr <- 25
#' sIn <- expand.grid(x1 = seq(0,1,length = sqrt(n.tr)), x2 = seq(0,1,length = sqrt(n.tr)))
#' fIn <- list(f1 = matrix(runif(n.tr*10), ncol = 10), f2 = matrix(runif(n.tr*22), ncol = 22))
#' sOut <- fgp_BB3(sIn, fIn, n.tr)
#' m1 <- fgpm(sIn = sIn, fIn = fIn, sOut = sOut)
#'
#' # making light simulations
#' n.sm <- 100
#' sIn.sm <- as.matrix(expand.grid(x1 = seq(0,1,length = sqrt(n.sm)),
#'                                 x2 = seq(0,1,length = sqrt(n.sm))))
#' fIn.sm <- list(f1 = matrix(runif(n.sm*10), ncol = 10), matrix(runif(n.sm*22), ncol = 22))
#' m1.sims <- simulate(m1, nsim = 10, sIn.sm = sIn.sm, fIn.sm = fIn.sm)
#'
#' # plotting light simulations
#' plotSims(m1, m1.sims)
#'
#'
#' # plotting full simulations________________________________________________________________
#' # building the model
#' set.seed(100)
#' n.tr <- 25
#' sIn <- expand.grid(x1 = seq(0,1,length = sqrt(n.tr)), x2 = seq(0,1,length = sqrt(n.tr)))
#' fIn <- list(f1 = matrix(runif(n.tr*10), ncol = 10), f2 = matrix(runif(n.tr*22), ncol = 22))
#' sOut <- fgp_BB3(sIn, fIn, n.tr)
#' m1 <- fgpm(sIn = sIn, fIn = fIn, sOut = sOut)
#'
#' # making full simulations
#' n.sm <- 100
#' sIn.sm <- as.matrix(expand.grid(x1 = seq(0,1,length = sqrt(n.sm)),
#'                                 x2 = seq(0,1,length = sqrt(n.sm))))
#' fIn.sm <- list(f1 = matrix(runif(n.sm*10), ncol = 10), matrix(runif(n.sm*22), ncol = 22))
#' m1.sims <- simulate(m1, nsim = 10, sIn.sm = sIn.sm, fIn.sm = fIn.sm, detail = "full")
#'
#' # plotting full simulations in full mode
#' plotSims(m1, m1.sims)
#'
#' # plotting full simulations in light mode
#' plotSims(m1, m1.sims, detail = "light")
#'
#' @importFrom graphics lines plot layout legend par matplot axis
#' @rdname plotSims-method
setMethod("plotSims", "fgpm",
          function(model, sims, detail = "full", ...) {
            plotSims.fgpm(sims = sims, detail = detail)
          })

plotSims.fgpm <- function(sims, detail, ...) {
  # recover realizations
  if (is.list(sims)) y_traj <- sims$sims else y_traj <- sims

  # recover graphic parameters if provided
  gpars <- list(...)
  # <---> limits
  if (!is.null(gpars$xlim)) xlim <- gpars$xlim else xlim <- c(1, ncol(y_traj))
  if (!is.null(gpars$ylim)) ylim <- gpars$ylim else ylim <- range(y_traj)
  # <---> line colors
  if (is.list(sims) & detail == "full") {
    if (!is.null(gpars$line.sm)) line.sm <- gpars$line.sm else line.sm <- "grey"
    if (!is.null(gpars$line.mean)) line.mean <- gpars$line.mean else line.mean <- "red"
    if (!is.null(gpars$line.ci)) line.ci <- gpars$line.ci else line.ci <- "blue"
  } else {
    if (!is.null(gpars$line.sm)) line.sm <- gpars$line.sm else line.sm <- "palegreen4"
  }
  # <---> line styles
  if (!is.null(gpars$lty.sm)) lty.sm <- gpars$lty.sm else lty.sm <- 1
  if (!is.null(gpars$lwd.sm)) lwd.sm <- gpars$lwd.sm else lwd.sm <- 1
  if (!is.null(gpars$lty.ci)) lty.ci <- gpars$lty.ci else lty.ci <- 1
  if (!is.null(gpars$lwd.ci)) lwd.ci <- gpars$lwd.ci else lwd.ci <- 1
  if (!is.null(gpars$lty.mean)) lty.mean <- gpars$lty.mean else lty.mean <- 1
  if (!is.null(gpars$lwd.mean)) lwd.mean <- gpars$lwd.mean else lwd.mean <- 1
  # <---> text
  if (!is.null(gpars$main)) main <- gpars$main else main <- "Simulations from a funGp model"
  if (!is.null(gpars$xlab)) xlab <- gpars$xlab else xlab <- "Sim. index"
  if (!is.null(gpars$ylab)) ylab <- gpars$ylab else ylab <- "Output"
  # <---> legend
  if (!is.null(gpars$legends)) legends <- gpars$legends else legends <- T

  # save current par state
  opar <- par('mar', 'mfrow')
  on.exit(par(opar))

  # set up layout
  par(mar = c(5.1, 4.1, 4.1, 2.1), mfrow = c(1,1))

  # plot
  matplot(t(y_traj), type = "l", col = line.sm, lty = lty.sm, lwd = lwd.sm, main = main, xlab = xlab, ylab = ylab, xaxt = "n")
  axis(1, axtags(ncol(y_traj)))
  if (all(is.list(sims), detail == "full")) {
    lines(sims$mean, col = line.mean)
    lines(sims$lower95, col = line.ci)
    lines(sims$upper95, col = line.ci)
    if (legends)
      legend("topleft", legend = c("Sims", "Mean", "95% CIs"), col = c(line.sm, line.mean, line.ci),
             lty = c(lty.sm, lty.mean, lty.ci), cex = 0.8)
  }
}
# ==========================================================================================================



# ==========================================================================================================
# Plot of decay functions used to set the pheromones for projection dimension in ACO
# ==========================================================================================================
#' @title Decay functions for ant colony optimization in funGp
#' @description This function is intended to aid the selection of the heuristic parameters \emph{tao0},
#'   \emph{delta} and \emph{dispr} in the call to the model selection function \link[funGp]{fgpm_factory}.
#'   The values computed by decay are the ones that would be used by the ant colony algorithm as initial
#'   pheromone load of the links pointing out to projection on each dimension. For more details, check the
#'   \href{https://hal.archives-ouvertes.fr/hal-02532713}{technical report}
#'   explaining the ant colony algorithm implemented in funGp, and the
#'   \href{https://hal.archives-ouvertes.fr/hal-02536624}{
#'   manual} of the package.
#'
#' @param k A number indicating the dimension of the functional input under analysis.
#' @param pmax An optional number specifying the hypothetical maximum projection dimension of this input. The
#'   user will be able to set this value later in the call to \link[funGp]{fgpm_factory} as a constraint. If
#'   not specified, it takes the value of k.
#' @param tao0 Explained in the description of \emph{dispr}.
#' @param delta Explained in the description of \emph{dispr}.
#' @param dispr The arguments \emph{tao0}, \emph{delta} and \emph{dispr}, are optional numbers specifying the
#'   loss function that determines the initial pheromone load on the links pointing out to projection
#'   dimensions. Such a function is defined as
#'
#'   \deqn{tao = tao0 * exp(-.5 * ((p - delta - 1)^2/(-dispr^2/(2*log(.5)),}
#'
#'   with p taking the values of the projection dimensions. The argument \emph{tao0} indicates the pheromone
#'   load in the links pointing out to the smallest dimensions; \emph{delta} specifies how many dimensions
#'   should preserve the maximum pheromone load; \emph{dispr} determines how fast the pheromone load drops
#'   in dimensions further than \eqn{delta + 1}. If \emph{pmax} = \emph{k}, then the dimension 0,
#'   representing no projection, receives a pheromone load identical to that of dimension \emph{k}. This, in
#'   order to represent the fact that both the representation of the function in its original dimension or
#'   a projection in a space of the same dimension, are equally heavy for the model. The default values of
#'   \emph{tao0}, \emph{delta} and \emph{dispr}, are 0.1, 2 and 1.4, respectively, which match the default
#'   values used by the \link[funGp]{fgpm_factory} function. Check
#'   \href{https://hal.archives-ouvertes.fr/hal-02532713}{this technical
#'   report} for more details.
#' @param doplot An optional boolean indicating if the pheromone loads should be plotted. Default = TRUE.
#' @param deliver An optional boolean indicating if the pheromone loads should be returned. Default = FALSE.
#'
#' @return If deliver is TRUE, an object of class \code{"numeric"} containing the initial pheromone values
#'   corresponding to the specified projection dimensions. Otherwise, the function plots the pheromones and
#'   nothing is returned.
#'
#' @author José Betancourt, François Bachoc and Thierry Klein
#'
#' @seealso \strong{*} \link[funGp]{fgpm_factory} for heuristic funGp model selection.
#'
#' @examples
#' # using default decay arguments____________________________________________________________
#' # input of dimension 15 projected maximum in dimension 15
#' decay(15)
#'
#' # input of dimension 15 projected maximum in dimension 8
#' decay(15, 8)
#'
#'
#' # playing with decay arguments_____________________________________________________________
#' # input of dimension 15 projected maximum in dimension 15
#' decay(15)
#'
#' # using a larger value of tao0
#' decay(15, tao0 = .3)
#'
#' # using a larger value of tao0, keeping it fixed up to higher dimensions
#' decay(15, tao0 = .3, delta = 5)
#'
#' # using a larger value of tao0, keeping it fixed up to higher dimensions, with slower decay
#' decay(15, tao0 = .3, delta = 5, dispr = 5.2)
#'
#'
#' # requesting pheromone values______________________________________________________________
#' # input of dimension 15 projected maximum in dimension 15
#' decay(15, deliver = TRUE)
#'
#' @importFrom graphics abline points axis
#' @export
decay <- function(k, pmax = NULL, tao0 = .1, delta = 2, dispr = 1.4, doplot = TRUE, deliver = FALSE) {
  if (is.null(pmax)) {
    pmin <- 0
    pmax <- k
  } else {
    if (pmax >= k) pmin <- 0 else pmin <- 1
    pmax <- min(k, pmax)
  }
  kvec <- 1:pmax
  tao <- sapply(kvec, function(i) {
    if (i <= delta) fact <- 1 else fact <- exp(-(.5) * ((abs(i-1) - delta)^2/(-dispr^2/(2*log(.5)))))
    tao0 * fact
    })
  if (pmax == k) tao <- c(tao[length(tao)], tao)
  if (doplot) {
    plot(1, type = "n", xlim = c(pmin, pmax), ylim = c(0,1), main = "Regularized pheromones",
         xlab = "Dimension", ylab = "Pheromone", xaxt = "n")
    if (pmax == k) axis(1, axtags(pmax+1)-1) else axis(1, axtags(pmax))
    if (length(kvec) <= 100) abline(v = pmin:pmax, lty = 2, col = "grey80")
    abline(h = c(0,1))
    points(pmin:pmax, tao, pch = 21, bg = "red")
  }

  if (deliver) {
    return(tao)
  }
}
# ==========================================================================================================



# ==========================================================================================================
# Plot of probability functions used for projection dimension in ACO
# ==========================================================================================================
#' @title Probability functions for ant colony optimization in funGp
#' @description This function is intended to aid the selection of the heuristic parameters \emph{tao0},
#'   \emph{delta} and \emph{dispr} in the call to the model selection function \link[funGp]{fgpm_factory}.
#'   The values computed by decay2probs are the ones that would be used by the ant colony algorithm as
#'   probability load of the links pointing out to projection on each dimension. These values result from
#'   the normalization of the initial pheromone loads delivered by the \link[funGp]{decay} function, which
#'   are made to sum 1. For more details, check the
#'   \href{https://hal.archives-ouvertes.fr/hal-02532713}{technical report}
#'   explaining the ant colony algorithm implemented in funGp, and the
#'   \href{https://hal.archives-ouvertes.fr/hal-02536624}{
#'   manual} of the package.
#'
#' @param k a number indicating the dimension of the functional input under analysis.
#' @param pmax an optional number specifying the hypothetical maximum projection dimension of this input. The
#'   user will be able to set this value later in the call to \link[funGp]{fgpm_factory} as a constraint. If
#'   not specified, it takes the value of k.
#' @param tao0 explained in the description of \emph{dispr}.
#' @param delta explained in the description of \emph{dispr}.
#' @param dispr the arguments \emph{tao0}, \emph{delta} and \emph{dispr}, are optional numbers specifying the
#'   loss function that determines the initial pheromone load on the links pointing out to projection
#'   dimensions. Such a function is defined as
#'
#'   \deqn{tao = tao0 * exp(-.5 * ((p - delta - 1)^2/(-dispr^2/(2*log(.5)),}
#'
#'   with p taking the values of the projection dimensions. The argument \emph{tao0} indicates the pheromone
#'   load in the links pointing out to the smallest dimensions; \emph{delta} specifies how many dimensions
#'   should preserve the maximum pheromone load; \emph{dispr} determines how fast the pheromone load drops
#'   in dimensions further than \eqn{delta + 1}. If \emph{pmax} = \emph{k}, then the dimension 0,
#'   representing no projection, receives a pheromone load identical to that of dimension \emph{k}. This, in
#'   order to represent the fact that both, the representation of the function in its original dimension or
#'   a projection in a space of the same dimension, are equally heavy for the model. In order to obtain the
#'   probability loads, the initial pheromone values are normalized to sum 1. Note that the normalization
#'   makes the value of tao0 become irrelevant in the initial probability load. This does not mean that the
#'   effect of tao0 is completely removed from the algorithm. Despite the fact that tao0 does not have
#'   influence on the selection of the projection dimension during the first iteration, it will be
#'   protagonist during the global pheromone update and will have an impact on every further iteration.
#'   The argument tao0 is left active in the input just for a better comprehension of the functioning of the
#'   mechanisms defining the initial pheromone and probability loads. The default values of
#'   \emph{tao0}, \emph{delta} and \emph{dispr}, are 0.1, 2 and 1.4, respectively, which match the default
#'   values used by the \link[funGp]{fgpm_factory} function. Check
#'   \href{https://hal.archives-ouvertes.fr/hal-02532713}{this technical
#'   report} for more details.
#' @param doplot an optional boolean indicating if the probability loads should be plotted. Default = TRUE.
#' @param deliver an optional boolean indicating if the probability loads should be returned. Default = FALSE.
#'
#' @return If deliver is TRUE, an object of class \code{"numeric"} containing the normalized initial pheromone values
#'   corresponding to the specified projection dimensions. Otherwise, the function plots the normalized
#'   pheromones and nothing is returned.
#'
#' @author José Betancourt, François Bachoc and Thierry Klein
#'
#' @seealso \strong{*} \link[funGp]{decay} for the function to generate the initial pheromone load;
#' @seealso \strong{*} \link[funGp]{fgpm_factory} for heuristic model selection in funGp.
#'
#' @examples
#' # using default decay arguments____________________________________________________________
#' # input of dimension 15 projected maximum in dimension 15
#' decay(15) # initial pheromone load
#' decay2probs(15) # initial probability load
#'
#' # input of dimension 15 projected maximum in dimension 8
#' decay(15, 8) # initial pheromone load
#' decay2probs(15, 8) # initial probability load
#'
#'
#' # playing with decay2probs arguments_______________________________________________________
#' # varying the initial pheromone load
#' decay(15) # input of dimension 15 projected maximum in dimension 15
#' decay(15, tao0 = .3) # larger value of tao0
#' decay(15, tao0 = .3, delta = 5) # larger tao0 kept to higher dimensions
#' decay(15, tao0 = .3, delta = 5, dispr = 5.2) # larger tao0 kept to higher dimensions
#'                                              # and slower decay
#'
#' # varying the initial probability load
#' decay2probs(15) # input of dimension 15 projected maximum in dimension 15
#' decay2probs(15, tao0 = .3) # larger value of tao0 (no effect whatsoever)
#' decay2probs(15, tao0 = .3, delta = 5) # larger tao0 kept to higher dimensions
#' decay2probs(15, tao0 = .3, delta = 5, dispr = 5.2) # larger tao0 kept to higher dimensions
#'                                                    # and slower decay
#'
#'
#' # requesting probability values____________________________________________________________
#' # input of dimension 15 projected maximum in dimension 15
#' decay2probs(15, deliver = TRUE)
#'
#' @importFrom graphics abline points axis
#' @export
decay2probs <- function(k, pmax = NULL, tao0 = .1, delta = 2, dispr = 1.4, doplot = TRUE, deliver = FALSE) {
  if (is.null(pmax)) {
    pmin <- 0
    pmax <- k
  } else {
    if (pmax >= k) pmin <- 0 else pmin <- 1
    pmax <- min(k, pmax)
  }
  v <- decay(k, pmax, tao0, delta, dispr, doplot = FALSE, deliver = TRUE)

  p <- v/sum(v)
  if (doplot) {
    plot(1, type = "n", xlim = c(pmin, pmax), ylim = c(0,1), main = "Regularized probabilities",
         xlab = "Dimension", ylab = "Probability", xaxt = "n")
    if (pmax == k) axis(1, axtags(pmax+1)-1) else axis(1, axtags(pmax))
    if (length(k) <= 100) abline(v = pmin:pmax, lty = 2, col = "grey80")
    abline(h = c(0,1))
    points(pmin:pmax, p, pch = 21, bg = "red")
  }

  if (deliver) {
    return(p)
  }
}
# ==========================================================================================================



# ==========================================================================================================
# Diagnostic calibration and fitness plot for Xfgpm objects
# ==========================================================================================================
#' @title Diagnostic plot for quality-enhanced models
#' @description This method provides plots for assessing the quality of regression models whose structure
#'   have been somehow optimized for predictability.
#'
#' @param x.model an object containing the model for which the quality plot is to be made.
#' @param ... additional arguments affecting the plot.
#'
#' @return None.
#'
#' @seealso \strong{*} \link[funGp]{plotX} for the diagnostic plot of a quality-enhanced funGp model.
#'
#' @examples
#' require(funGp) # a package with a plotX method implemented
#'
#' # generating input and output data
#' set.seed(100)
#' n.tr <- 2^5
#' sIn <- expand.grid(x1 = seq(0,1,length = n.tr^(1/5)), x2 = seq(0,1,length = n.tr^(1/5)),
#'                    x3 = seq(0,1,length = n.tr^(1/5)), x4 = seq(0,1,length = n.tr^(1/5)),
#'                    x5 = seq(0,1,length = n.tr^(1/5)))
#' fIn <- list(f1 = matrix(runif(n.tr*10), ncol = 10), f2 = matrix(runif(n.tr*22), ncol = 22))
#' sOut <- fgp_BB7(sIn, fIn, n.tr)
#' \dontrun{
#' # optimizing the model structure with fgpm_factory (~5 seconds)
#' xm <- fgpm_factory(sIn = sIn, fIn = fIn, sOut = sOut)
#'
#' # assessing the quality of the model - absolute and w.r.t. the other explored models
#' plotX(xm)
#' }
#'
#' @name plotX
#' @export
#' @keywords internal
setGeneric("plotX", function(x.model, ...) standardGeneric("plotX"))

#' @title Diagnostic plots for funGp factory output
#' @description This method provides two plots for assessing the quality of the output delivered by the
#'   model selection algorithm in the \link[funGp]{fgpm_factory} function. The first one is a calibration
#'   plot similar to the one offered for \linkS4class{fgpm} objects by the \link[funGp]{plotLOO} function.
#'   This plot allows to validate the absolute quality of the selected model. The second one displays the
#'   performance statistic of all the models successfully evaluated by the model selection algorithm. This
#'   provides a notion of the relative quality of the selected model with respect to the other models that
#'   can be made using the same data.
#'
#' @param x.model an object of class \linkS4class{Xfgpm} containing the output of the model selection
#'   algorithm in \link[funGp]{fgpm_factory}.
#' @param calib a boolean indicating whether the calibration plot of the selected model should be included
#'   in the display. Default is TRUE.
#' @param fitp a boolean indicating whether scatter plot of the quality of all explored models should be
#'   included in the display. Default is TRUE.
#'
#' @param horiz logical. If \code{TRUE} and if both \code{calib} and
#'     \code{fitp} the two subplots corresponding to the calibration
#'     and the fit quality are displayer horizontally (on a same row)
#'     rather than vertically as in the default behaviour.
#'
#' @param ... additional arguments affecting the display. Since this method allows to generate two plots
#'   from a single function call, the extra arguments for each plot should be included in a list. For the
#'   calibration plot, the list should be called \emph{calib.gpars}. For the plot of the fitness of
#'   explored models, the list should be called \emph{fitp.gpars}. The following typical graphics parameters
#'   are valid entries of both lists: \emph{xlim}, \emph{ylim}, \emph{xlab}, \emph{ylab}, \emph{main}. The
#'   boolean argument legends can also be included in any of the two lists in order to control the display
#'   of legends in the corresponding plot.
#'
#' @return None.
#'
#' @author José Betancourt, François Bachoc and Thierry Klein
#'
#' @references Betancourt, J., Bachoc, F., Klein, T., and Gamboa, F. (2020),
#' Technical Report: "Ant Colony Based Model Selection for Functional-Input Gaussian Process Regression. Ref. B3D-WP3.2".
#' \emph{RISCOPE project}.
#' \href{https://hal.archives-ouvertes.fr/hal-02532713}{[HAL]}
#'
#' @references Betancourt, J., Bachoc, F., and Klein, T. (2020),
#' R Package Manual: "Gaussian Process Regression for Scalar and Functional Inputs with funGp - The in-depth tour".
#' \emph{RISCOPE project}.
#' \href{https://hal.archives-ouvertes.fr/hal-02536624}{[HAL]}
#'
#' @seealso \strong{*} \link[funGp]{fgpm_factory} for structural optimization of funGp models;
#' @seealso \strong{*} \link[funGp]{plotEvol} for a plot on the evolution of the model selection algorithm
#'   in fgpm_factory.
#'
#' @examples
#' # generating input and output data
#' set.seed(100)
#' n.tr <- 2^5
#' sIn <- expand.grid(x1 = seq(0,1,length = n.tr^(1/5)), x2 = seq(0,1,length = n.tr^(1/5)),
#'                    x3 = seq(0,1,length = n.tr^(1/5)), x4 = seq(0,1,length = n.tr^(1/5)),
#'                    x5 = seq(0,1,length = n.tr^(1/5)))
#' fIn <- list(f1 = matrix(runif(n.tr*10), ncol = 10), f2 = matrix(runif(n.tr*22), ncol = 22))
#' sOut <- fgp_BB7(sIn, fIn, n.tr)
#' \dontrun{
#' # optimizing the model structure with fgpm_factory (~5 seconds)
#' xm <- fgpm_factory(sIn = sIn, fIn = fIn, sOut = sOut)
#'
#' # assessing the quality of the model - absolute and w.r.t. the other explored models
#' plotX(xm)
#'
#' # customizing some graphical parameters
#' plotX(xm, calib.gpars = list(xlim = c(800,1000), ylim = c(600,1200)),
#'           fitp.gpars = list(main = "Relative quality", legends = FALSE))
#' }
#'
#' @importFrom graphics lines points plot layout legend par arrows abline axis
#' @importFrom scales alpha
#' @rdname plotX-method
setMethod("plotX", "Xfgpm",
          function(x.model, calib = TRUE, fitp = TRUE, horiz = FALSE, ...) {
              plotX.Xfgpm(x.model = x.model, calib = calib, fitp = fitp,
                          horiz = horiz, ...)
          })

plotX.Xfgpm <- function(x.model, calib = TRUE, fitp = TRUE,
                        horiz = FALSE, ...) {

    ## recover graphic parameters if provided
    gpars <- list(...)
    cal.gpars <- gpars$calib.gpars
    fit.gpars <- gpars$fitp.gpars

  # loocv calibration plot _______________________________________________
  plot.c <- function(model) {
    # recover observed output
    y_obs <- model@sOut

    # compute loocv predictions
      R <- tcrossprod(model@preMats$L)/model@kern@varHyp +
          diag(model@nugget, nrow = model@n.tr, ncol = model@n.tr)
    Rinv <- solve(R)
    y_pre <- y_obs - diag(Rinv)^(-1) * Rinv %*% y_obs

    # compute LOO statistic
    q2 <- format(getFitness(x.model@model), digits = 3, nsmall = 3)

    # recover graphic parameters if provided
    # gpars <- list(...)
    if (!is.null(cal.gpars$xlim)) xlim <- cal.gpars$xlim else xlim <- range(c(y_obs, y_pre))
    if (!is.null(cal.gpars$ylim)) ylim <- cal.gpars$ylim else ylim <- range(c(y_obs, y_pre))
    if (!is.null(cal.gpars$pch)) pch <- cal.gpars$pch else pch <- 21
    if (!is.null(cal.gpars$pt.col)) pt.col <- cal.gpars$pt.col else pt.col <- "red"
    if (!is.null(cal.gpars$pt.bg)) pt.bg <- cal.gpars$pt.bg else pt.bg <- "red"
    if (!is.null(cal.gpars$pt.cex)) pt.cex <- cal.gpars$pt.cex else pt.cex <- 1
    if (!is.null(cal.gpars$line)) line <- cal.gpars$line else line <- "blue"
    if (!is.null(cal.gpars$xlab)) xlab <- cal.gpars$xlab else xlab <- "Observed"
    if (!is.null(cal.gpars$ylab)) ylab <- cal.gpars$ylab else ylab <- "Predicted"
    if (!is.null(cal.gpars$main)) main <- cal.gpars$main else main <- "Model diagnostic by leave-one-out cross-validation"
    if (!is.null(cal.gpars$legends)) legends <- cal.gpars$legends else legends <- T

    # add points
    plot(y_obs, y_pre, xlim = xlim, ylim = ylim, pch = pch, col = pt.col, bg = pt.bg, cex = pt.cex,
         main = main, xlab = "", ylab = ylab)
    mtext(side = 1, text = xlab, line = 2.0)

    # add reference line
    lines(xlim, ylim, col = line)

    # add Q2 based on x.model stat and fitness
    if (legends) {
      legend("topleft", legend = paste("Q2loocv =", q2),
             xjust = 0.5, yjust = 0.5, x.intersp = -0.5, y.intersp = 0.3,
             adj = c(0, 0.5), inset = c(.02, .05))
    }
  }
  # ______________________________________________________________________


  # fitness plot _________________________________________________________
  plot.f <- function(tight) {
    y <- x.model@log.success@fitness
    # gpars <- list(...)
    if (!is.null(fit.gpars$pt.cex)) pt.cex <- fit.gpars$pt.cex else pt.cex <- 1
    if (!is.null(fit.gpars$lwd)) lwd <- fit.gpars$lwd else lwd <- 2
    if (!is.null(fit.gpars$xlim)) xlim <- fit.gpars$xlim else xlim <- c(1,length(y))
    if (!is.null(fit.gpars$ylim)) ylim <- fit.gpars$ylim else ylim <- c(0.038,1)
    if (!is.null(fit.gpars$xlab)) xlab <- fit.gpars$xlab else xlab <- "Index @log.success"
    if (!is.null(fit.gpars$ylab)) ylab <- fit.gpars$ylab else ylab <- strsplit(x.model@stat, split = '\\(')[[1]][1]
    if (!is.null(fit.gpars$main)) main <- fit.gpars$main else main <- "Fitness of explored models"
    if (!is.null(fit.gpars$legends)) legends <- fit.gpars$legends else legends <- T

    plot(y, xlim = xlim, ylim = ylim, pch = 21, bg = "blue", main = main,
         xlab = "", ylab = ylab, xaxt = "n", yaxt = "n", cex = pt.cex)
    # axis(1, axtags(length(y)))
    axis(1, seq_len(length(y)))
    mtext(side = 1, text = xlab, line = 2.0)
    tag <- format(seq(0, 1, length.out = 4), digits = 2, nsmall = 2)
    axis(2, tag, tag, las = 1)
    abline(h = 1, lty = 2, col = "grey57")
    abline(h = x.model@fitness, lty = 2, col = "green", lwd = lwd)
    points(x.model@fitness, pch = 24, bg = "green", cex = pt.cex*1.1)
    points(x.model@fitness, pch = 25, bg = "green", cex = pt.cex*1.1)
    piv.arr <- which(y < 0)
    if (tight) y1 <- .078 else y1 <- .045
    for (i in seq_along(piv.arr)) {
      arrows(x0 = piv.arr[i], x1 = piv.arr[i], y0 = .01, y1 = y1, length = 0.06, code = 1,
             angle = 25, lwd = 1.5, col = alpha("blue", .7))
    }
    if (legends) {
      if (x.model@fitness > .7) {
        legend("bottomleft", legend = c("Sel. model", "All models", "fitness < 0"),
               pch = c(24, 21, NA), pt.bg = c("green", "blue", NA), inset = c(.02,.08))
        legend("bottomleft", legend = c("Sel. model", "All models", "fitness < 0"),
               pch = c(25, 21, NA), pt.bg = c("green", "blue", NA), inset = c(.02,.08), bg = NA)

        # arrow parameters in x
        limscr <- par('usr')[1:2]
        lcr <- limscr[1] + diff(limscr)*(0.02) # left coordinate of legend frame
        win <- par()$pin[1]*(1-0.02) # width of plot in inches (omitting left inset)
        pin <- 0.15 # location of point in inches
        lp <- pin/win # location as percentage of width
        wcr <- diff(limscr)*(1-0.02) # width of plot in coordinate system
        xcr <- lcr + wcr*lp # location of point in coordinate

        # arrow parameters in y
        limscr <- par('usr')[3:4]
        bcr <- limscr[1] + diff(limscr)*(0.08) # bottom coordinate of legend frame
        hin <- par()$pin[2]*(1-0.08) # height of plot in inches (omitting bottom inset)
        pin1 <- 0.14 # end location of arrow in inches
        h <- .12 # arrow height in inches
        pin2 <- pin1+h # start location of arrow in inches
        lp1 <- pin1/hin # location 1 as percentage of height
        lp2 <- pin2/hin # location 2 as percentage of height
        hcr <- diff(limscr)*(1-0.08) # height of plot in coordinate system
        ycr1 <- bcr + hcr*lp1 # location 1 in coordinate
        ycr2 <- bcr + hcr*lp2 # location 1 in coordinate

      } else {
        legend("topright", legend = c("Sel. model", "All models", "fitness < 0"),
               pch = c(24, 21, NA), pt.bg = c("green", "blue", NA), inset = c(.02,.08))
        legend("topright", legend = c("Sel. model", "All models", "fitness < 0"),
               pch = c(25, 21, NA), pt.bg = c("green", "blue", NA), inset = c(.02,.08), bg = NA)

        # arrow parameters in x
        limscr <- par('usr')[1:2]
        win <- par()$pin[1]*(1-0.02) # width of plot in inches (omitting right inset)
        pin <- win - 1.015 # location of point in inches
        lp <- pin/win # location as percentage of width
        wcr <- diff(limscr)*(1-0.02) # width of plot in coordinate system
        xcr <- limscr[1] + wcr*lp # location of point in coordinate

        # arrow parameters in x
        limscr <- par('usr')[3:4]
        hin <- par()$pin[2]*(1-0.08) # height of plot in inches (omitting bottom inset)
        pin1 <- hin - 0.659 # end location of arrow in inches (distance from top)
        h <- .12 # arrow height in inches
        pin2 <- pin1+h # start location of arrow in inches
        lp1 <- pin1/hin # location 1 as percentage of height
        lp2 <- pin2/hin # location 2 as percentage of height
        hcr <- diff(limscr)*(1-0.08) # height of plot in coordinate system
        ycr1 <- limscr[1] + hcr*lp1 # location 1 in coordinate
        ycr2 <- limscr[1] + hcr*lp2 # location 1 in coordinate
      }

      # arrow display
      arrows(x0 = xcr, x1 = xcr, y0 = ycr1, y1 = ycr2, length = 0.06, code = 1,
             angle = 25, lwd = 2.5, col = alpha("blue", .7))
    }
  }
  # ______________________________________________________________________

  # save current par state
  opar <- par('mar', 'mfrow')
  on.exit(par(opar))

  # plot
    if (all(calib, fitp)) {
        if (!horiz) par(mar = c(3.1, 4.1, 2.5, 2.1), mfrow = c(2, 1))
        else par(mar = c(3.1, 4.1, 2.5, 2.1), mfcol = c(1, 2))
        plot.c(x.model@model)
        plot.f(T)
    } else if (calib) {
        par(mar = c(5.1, 4.1, 4.1, 2.1), mfrow = c(1,1))
        plot.c(x.model@model)
    } else if (fitp) {
        par(mar = c(5.1, 4.1, 4.1, 2.1), mfrow = c(1,1))
        plot.f(F)
    }
}
# ==========================================================================================================



# ==========================================================================================================
# Diagnostic plot for ACO evolution along iterations
# ==========================================================================================================
#' @title Plot for the evolution of model selection algorithm
#' @description This method displays the evolution of an iterative algorithm for model selection.
#'
#' @param x.model an object containing the data structures returned by the model selection algorithm.
#' @param ... additional arguments affecting the plot.
#'
#' @return None.
#'
#' @seealso \strong{*} \link[funGp]{plotEvol} for a plot on the evolution of the model selection algorithm
#'   in fgpm_factory.
#'
#' @examples
#' require(funGp) # a package with a plotEvol method implemented
#'
#' # generating input and output data
#' set.seed(100)
#' n.tr <- 2^5
#' sIn <- expand.grid(x1 = seq(0,1,length = n.tr^(1/5)), x2 = seq(0,1,length = n.tr^(1/5)),
#'                    x3 = seq(0,1,length = n.tr^(1/5)), x4 = seq(0,1,length = n.tr^(1/5)),
#'                    x5 = seq(0,1,length = n.tr^(1/5)))
#' fIn <- list(f1 = matrix(runif(n.tr*10), ncol = 10), f2 = matrix(runif(n.tr*22), ncol = 22))
#' sOut <- fgp_BB7(sIn, fIn, n.tr)
#' \dontrun{
#' # optimizing the model structure with fgpm_factory (~5 seconds)
#' xm <- fgpm_factory(sIn = sIn, fIn = fIn, sOut = sOut)
#'
#' # assessing the evolution of the algorihtm
#' plotEvol(xm)
#' }
#'
#' @name plotEvol
#' @export
#' @keywords internal
setGeneric("plotEvol", function(x.model, ...) standardGeneric("plotEvol"))

#' @title Plot for the evolution of model selection algorithm in funGp
#' @description This method displays the evolution of the quality of the configurations evaluated along
#'   the iterations, by the model selection algorithm in the \link[funGp]{fgpm_factory} function. For
#'   each iteration, the performance statistic of all the evaluated models is printed, along with the
#'   corresponding median of the group. The plot also includes the global maximum, which corresponds
#'   to the best performance statistic obtained up to the current iteration. In this plot, it is
#'   typical to have some points falling relatively far from the maximum, even after multiple
#'   iterations. This happens mainly because we have multiple categorical features, whose alteration
#'   might change the performance statistic in a nonsmooth way. On the other hand, the points that fall
#'   bellow zero usually correspond to models whose hyperparameters were hard to optimize. This occurs
#'   sporadically during the log-likelihood optimization for Gaussian processes, due to the
#'   non-linearity of the objective function. As long as the maximum keeps improving and the median
#'   remains close to it, none of the two aforementioned phenomena is matter for worries. Both of them
#'   respond to the mechanism of exploration implemented in the algorithm, which makes it able to
#'   progressively move towards better model configurations.
#'
#' @param x.model an object of class \linkS4class{Xfgpm} containing the output of the model selection
#'   algorithm in \link[funGp]{fgpm_factory}.
#' @param ... additional arguments affecting the plot. The following typical graphics parameters are
#'   valid entries: \emph{xlim}, \emph{ylim}, \emph{xlab}, \emph{ylab}, \emph{main}.
#'
#' @return None.
#'
#' @author José Betancourt, François Bachoc and Thierry Klein
#'
#' @references Betancourt, J., Bachoc, F., and Klein, T., and Gamboa, F. (2020),
#' Technical Report: "Ant Colony Based Model Selection for Functional-Input Gaussian Process Regression. Ref. D3.b (WP3.2)".
#' \emph{RISCOPE project}.
#' \href{https://hal.archives-ouvertes.fr/hal-02532713}{[HAL]}
#'
#' @references Betancourt, J., Bachoc, F., and Klein, T. (2020),
#' R Package Manual: "Gaussian Process Regression for Scalar and Functional Inputs with funGp - The in-depth tour".
#' \emph{RISCOPE project}.
#' \href{https://hal.archives-ouvertes.fr/hal-02536624}{[HAL]}
#'
#' @seealso \strong{*} \link[funGp]{fgpm_factory} for structural optimization of funGp models;
#' @seealso \strong{*} \link[funGp]{plotX} for diagnostic plots for a fgpm_factory output and selected model.
#'
#' @examples
#' # generating input and output data
#' set.seed(100)
#' n.tr <- 2^5
#' sIn <- expand.grid(x1 = seq(0,1,length = n.tr^(1/5)), x2 = seq(0,1,length = n.tr^(1/5)),
#'                    x3 = seq(0,1,length = n.tr^(1/5)), x4 = seq(0,1,length = n.tr^(1/5)),
#'                    x5 = seq(0,1,length = n.tr^(1/5)))
#' fIn <- list(f1 = matrix(runif(n.tr*10), ncol = 10), f2 = matrix(runif(n.tr*22), ncol = 22))
#' sOut <- fgp_BB7(sIn, fIn, n.tr)
#' \dontrun{
#' # optimizing the model structure with fgpm_factory (~5 seconds)
#' xm <- fgpm_factory(sIn = sIn, fIn = fIn, sOut = sOut)
#'
#' # assessing the evolution of the algorithm
#' plotEvol(xm)
#' }
#'
#' @importFrom graphics lines points plot layout legend par arrows axis
#' @importFrom scales alpha
#' @importFrom stats median
#' @rdname plotEvol-method
setMethod("plotEvol", "Xfgpm",
          function(x.model, ...) {
            plotEvol.Xfgpm(x.model = x.model, ...)
          })

plotEvol.Xfgpm <- function(x.model, ...) {
  # recover fitness list
  all.fitness <- x.model@details$evolution
  n.gen <- length(all.fitness)

  # recover graphic parameters if provided
  gpars <- list(...)
  # <---> limits
  if (!is.null(gpars$xlim)) xlim <- gpars$xlim else xlim <- c(.7, (n.gen + .3))
  if (!is.null(gpars$ylim)) ylim <- gpars$ylim else ylim <- c(-.045, 1)
  # <---> point colors
  if (!is.null(gpars$pt.col)) pt.col <- gpars$pt.col else pt.col <- "blue"
  if (!is.null(gpars$pt.med)) pt.med <- gpars$pt.med else pt.med <- "red"
  # <---> line colors
  if (!is.null(gpars$col.arr)) col.arr <- gpars$col.arr else col.arr <- "red"
  # <---> text
  if (!is.null(gpars$main)) main <- gpars$main else main <- "ACO learning process"
  if (!is.null(gpars$xlab)) xlab <- gpars$xlab else xlab <- "Iteration"
  if (!is.null(gpars$ylab)) ylab <- gpars$ylab else ylab <- x.model@stat

  # save current par state
  opar <- par('mar', 'mfrow')
  on.exit(par(opar))

  # set up layout
  par(mar = c(5.1, 4.1, 4.1, 2.1), mfrow = c(1,1))

  # plot
  plot(1, type = "n", xlim = xlim, ylim = ylim, main = main, xlab = xlab, ylab = ylab, xaxt = "n", yaxs = "i")
  # axis(1, 1:n.gen)
  axis(1, axtags(n.gen))
  abline(h = 0, lty = 3, col = "grey75")
  if (n.gen <= 100) abline(v = seq_len(n.gen), lty = 3, col = "grey75")
  for (i in 1:n.gen) {
    ys <- all.fitness[[i]]
    points(rep(i, length(which(ys>=0))), ys[ys>=0], pch = 21, bg = "blue", col = "grey55", cex = 1)
    if (any(ys < 0)) {
      arrows(x0 = i, x1 = i, y0 = -.035, y1 = -0.005, length = 0.06, code = 1, angle = 25, lwd = 1.5,
             col = alpha(col.arr, .7))
    }
  }
  lines(sapply(all.fitness, median), col = "red")
  points(sapply(all.fitness, median), pch = 21, bg = "red")
  lines(cummax(sapply(all.fitness, max)), col = "green")
  points(cummax(sapply(all.fitness, max)), pch = 21, bg = "green")

  legend("bottomright", legend = c("All models", "Global max.", "Med. per gen.", "Points < 0"),
         pch = c(21, 21, 21, NA), col = c("blue", "green", "red", NA), pt.bg = c("blue", "green", "red", NA),
         lty = c(0, 1, 1, 0), pt.cex = c(1, 1, 1, NA), cex = .85, inset = c(.02,.07))

  # arrow parameters in x
  limscr <- par('usr')[1:2]
  win <- par()$pin[1]*(1-0.02) # width of plot in inches (omitting right inset)
  pin <- win - 1.24 # location of point in inches
  lp <- pin/win # location as percentage of width
  wcr <- diff(limscr)*(1-0.02) # width of plot in coordinate system
  xcr <- limscr[1] + wcr*lp # location of point in coordinate

  # arrow parameters in y
  limscr <- par('usr')[3:4]
  bcr <- limscr[1] + diff(limscr)*(0.07) # bottom coordinate of legend frame
  hin <- par()$pin[2]*(1-0.07) # height of plot in inches (omitting bottom inset)
  pin1 <- 0.11 # end location of arrow in inches
  h <- .12 # arrow height in inches
  pin2 <- pin1+h # start location of arrow in inches
  lp1 <- pin1/hin # location 1 as percentage of height
  lp2 <- pin2/hin # location 2 as percentage of height
  hcr <- diff(limscr)*(1-0.07) # height of plot in coordinate system
  ycr1 <- bcr + hcr*lp1 # location 1 in coordinate
  ycr2 <- bcr + hcr*lp2 # location 1 in coordinate

  # arrow display
  arrows(x0 = xcr, x1 = xcr, y0 = ycr1, y1 = ycr2, length = 0.06, code = 1,
         angle = 25, lwd = 1.5, col = alpha(col.arr, .7))
}

