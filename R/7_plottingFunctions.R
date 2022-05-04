## ==========================================================================================================
## Diagnostic calibration plot for funGp models
## ==========================================================================================================



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


#' @importFrom graphics lines plot polygon layout legend par mtext
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





#' @importFrom graphics lines plot layout legend par matplot axis
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
#' @author José Betancourt, François Bachoc, Thierry Klein and Jérémy Rohmer
#'
#' @seealso \strong{*} \link[funGp]{decay2probs} for the function to generate the initial probability load;
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
#' @param doplot An optional boolean indicating if the probability loads should be plotted. Default = TRUE.
#' @param deliver An optional boolean indicating if the probability loads should be returned. Default = FALSE.
#'
#' @return If deliver is TRUE, an object of class \code{"numeric"} containing the normalized initial pheromone values
#'   corresponding to the specified projection dimensions. Otherwise, the function plots the normalized
#'   pheromones and nothing is returned.
#'
#' @author José Betancourt, François Bachoc, Thierry Klein and Jérémy Rohmer
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



##' @importFrom graphics lines points plot layout legend par arrows abline axis
##'
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


#' @importFrom graphics lines points plot layout legend par arrows axis
#' @importFrom scales alpha
#' @importFrom stats median
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

