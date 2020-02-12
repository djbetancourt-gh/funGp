# funGp plotter
# ----------------------------------------------------------------------------------------------------------
#' @name plotLOO
#' @description This is my description
#' @rdname plotLOO-methods
#' @importFrom graphics lines plot
#' @param model An object to predict from.
#' @param ... fill
#'
#' @author José Betancourt, François Bachoc and Thierry Klein
#' @exportMethod plotLOO
if(!isGeneric("plotLOO")) {setGeneric("plotLOO", function(model, ...) standardGeneric("plotLOO"))}

#' @title Prediction Method for the apk Class
#' @name plotLOO
#' @rdname plotLOO-methods
#' @aliases plotLOO,funGp-method
setMethod("plotLOO", "funGp", function(model, ...) {
  plotLOO.funGp(model = model, ...)
})

plotLOO.funGp <- function(model, ...) {
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
  if (!is.null(gpars$main)) main <- gpars$main else main <- "Model diagnostic by leave-one-out cross-valitation"

  # save current par state
  opar <- par('mar', 'mfrow')
  on.exit(par(opar))

  # set up layout
  par(mar = c(5.1, 4.1, 4.1, 2.1), mfrow = c(1,1))

  # plot
  plot(y_obs, y_pre, xlim = xlim, ylim = ylim, pch = pch, col = pt.col, bg = pt.bg, cex = pt.cex,
       main = main, xlab = xlab, ylab = ylab)
  lines(xlim, ylim, col = line)
  legend("topleft", legend = paste("Q2loocv =", q2),
         xjust = 0.5, yjust = 0.5, x.intersp = -0.5, y.intersp = 0.3,
         adj = c(0, 0.5), inset = c(.02,.05))
}
# ----------------------------------------------------------------------------------------------------------


# plotter for sorted predictions
# ----------------------------------------------------------------------------------------------------------
#' @name plotPreds
#' @description This is my description
#' @rdname plotPreds-methods
#' @importFrom graphics lines plot polygon layout legend par mtext
#' @param model An object to predict from.
#' @param preds An object to predict from.
#' @param ... Further arguments for methods.
#'
#' @author José Betancourt, François Bachoc and Thierry Klein
#' @exportMethod plotPreds
if(!isGeneric("plotPreds")) {setGeneric("plotPreds", function(model, preds, ...) standardGeneric("plotPreds"))}

#' @title Prediction Method for the apk Class
#' @name plotPreds
#' @rdname plotPreds-methods
#' @aliases plotPreds,funGp-method
#'
#' @param model fill
#' @param preds fill
#' @param sOut.pr also
#' @param calib also
#' @param legends also
#' @param stat also
setMethod("plotPreds", "funGp",
          function(model, preds, sOut.pr = NULL, calib = T, legends = T, stat = T, ...) {
            plotPreds.funGp(preds = preds, sOut.pr = sOut.pr,
                            calib = calib, legends = legends, stat, ...)
          })

plotPreds.funGp <- function(preds, sOut.pr, calib, legends, stat, ...) {
  # recover observed and predicted output
  y_obs <- sOut.pr[order(preds$mean)]
  y_pre <- sort(preds$mean)

  # recover 95% confidence intervals
  ll <- preds$lower95[order(preds$mean)]
  ul <- preds$upper95[order(preds$mean)]

  # recover graphic parameters if provided
  gpars <- list(...)
  lin.gpars <- gpars$lin.gpars
  cal.gpars <- gpars$cal.gpars

  # calibration plot _____________________________________________________
  plot.c <- function (tight = F) {
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

    # compute LOO statistic
    q2 <- format(1 - mean(sum((y_obs - y_pre)^2))/mean(sum((y_obs - mean(y_obs))^2)), digits = 3, nsmall = 3)

    # add hould-out Q2
    legend("topleft", legend = paste("Q2hout =", q2),
           xjust = 0.5, yjust = 0.5, x.intersp = -0.5, y.intersp = 0.3,
           adj = c(0, 0.5), inset = c(.02,.05))
  }
  # ______________________________________________________________________


  # line plot ____________________________________________________________
  plot.l <- function (tight = F) {
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
    if (!is.null(lin.gpars$xlab.l)) xlab <- lin.gpars$xlab else xlab <- "Index"
    if (!is.null(lin.gpars$ylab.l)) ylab <- lin.gpars$ylab else ylab <- "Output"
    if (!is.null(lin.gpars$main.l)) main <- lin.gpars$main else main <- "Sorted predictions"

    # add confidence intervals
    if (tight) {
      plot(1, type = "n", xlim = xlim, ylim = ylim, main = main, xlab = "", ylab = ylab)
      mtext(side = 1, text = xlab, line = 2.0)
    } else {
      plot(1, type = "n", xlim = xlim, ylim = ylim, main = main, xlab = xlab, ylab = ylab)
    }
    x <- seq_along(y_pre)
    polygon(c(x, rev(x)), c(ul, rev(ll)), col = col.poly, border = NA)
    lines(ll, col = line.ci, lty = lty.ci, lwd = lwd.ci)
    lines(ul, col = line.ci, lty = lty.ci, lwd = lwd.ci)

    # add predicted mean
    lines(y_pre, col = line.pr, lty = lty.pr, lwd = lwd.pr)

    # add true output if provided
    if (!is.null(y_obs)) {
      lines(y_obs, col = line.true, lty = lty.true, lwd = lwd.true)
      if (legends) {
        if (tight) {
          legend("topleft", legend = c("True", "Pred. mean", "95% CIs"),
                 col = c(line.true, line.pr, line.ci), lty = c(lty.true, lty.pr, lty.ci),
                 cex = 0.6, inset = c(.02,.05))
        } else {
          legend("topleft", legend = c("True", "Pred. mean", "95% CIs"),
                 col = c(line.true, line.pr, line.ci), lty = c(lty.true, lty.pr, lty.ci),
                 cex = 0.8, inset = c(.02,.05))
        }
      }

    } else {
      if (legends)
        legend("topleft", legend = c("Pred. mean", "95% CIs"),
               col = c(line.pr, line.ci), lty = c(lty.pr, lty.ci), cex = 0.8, inset = c(.02,.05))
    }
  }
  # ______________________________________________________________________


  # save current par state
  opar <- par('mar', 'mfrow')
  on.exit(par(opar))

  # identify the case
  # case 1: plotPreds(model, preds) -> plot.l
  # case 2: plotPreds(model, preds, sOut.pr) -> [plot.c, plot.l]
  # case 3: plotPreds(model, preds, sOut.pr, calib = F) -> [plot.l]
  if (is.null(sOut.pr)) {
    par(mar = c(5.1, 4.1, 4.1, 2.1), mfrow = c(1,1))
    plot.l()
  } else if (calib) {
    par(mar = c(3.1, 4.1, 2.5, 2.1), mfrow = c(2,1))
    plot.c(T)
    plot.l(T)
  } else {
    par(mar = c(5.1, 4.1, 4.1, 2.1), mfrow = c(1,1))
    plot.l()
  }
}
# ----------------------------------------------------------------------------------------------------------


# plotter for sorted simulation
# ----------------------------------------------------------------------------------------------------------
#' @name plotSims
#' @description This is my description
#' @rdname plotSims-methods
#' @importFrom graphics lines plot layout legend par matplot
#' @param model An object to predict from.
#' @param sims An object to predict from.
#' @param ... Further arguments for methods.
#'
#' @author José Betancourt, François Bachoc and Thierry Klein
#' @exportMethod plotSims
if(!isGeneric("plotSims")) {setGeneric("plotSims", function(model, sims, ...) standardGeneric("plotSims"))}

#' @title Prediction Method for the apk Class
#' @name plotSims
#' @rdname plotSims-methods
#' @aliases plotSims,funGp-method
#'
#' @param detail fill!!!!
#' @param legends fill!!!!
setMethod("plotSims", "funGp",
          function(model, sims, detail = "full", legends = T, ...) {
            plotSims.funGp(sims = sims, detail = detail, legends = legends)
          })

plotSims.funGp <- function(sims, detail, legends, ...) {
  # recover realizations
  if (is.list(sims)) y_traj <- sims$obs else y_traj <- sims

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

  # save current par state
  opar <- par('mar', 'mfrow')
  on.exit(par(opar))

  # set up layout
  par(mar = c(5.1, 4.1, 4.1, 2.1), mfrow = c(1,1))

  # plot
  matplot(t(y_traj), type = "l", col = line.sm, lty = lty.sm, lwd = lwd.sm, main = main, xlab = xlab, ylab = ylab)
  if (all(is.list(sims), detail == "full")) {
    lines(sims$mean, col = line.mean)
    lines(sims$lower95, col = line.ci)
    lines(sims$upper95, col = line.ci)
    if (legends)
      legend("topleft", legend = c("Sims", "Mean", "95% CIs"), col = c(line.sm, line.mean, line.ci),
             lty = c(lty.sm, lty.mean, lty.ci), cex = 0.8)
  }
}
# ----------------------------------------------------------------------------------------------------------


# plotter for X-funGp
# ----------------------------------------------------------------------------------------------------------
#' @name plotX
#' @description This is my description
#' @rdname plotX-methods
#' @importFrom graphics lines plot layout legend par arrows
#' @param x.model An object to predict from.
#' @param ... Further arguments for methods.
#'
#' @author José Betancourt, François Bachoc and Thierry Klein
#' @exportMethod plotX
if(!isGeneric("plotX")) {setGeneric("plotX", function(x.model, ...) standardGeneric("plotX"))}

#' @title Prediction Method for the apk Class
#' @name plotX
#' @rdname plotX-methods
#' @aliases plotX,XfunGp-method
setMethod("plotX", "XfunGp",
          function(x.model, ...) {
            plotX.XfunGp(x.model = x.model, ...)
          })

plotX.XfunGp <- function(x.model, ...) {
  # loocv calibration plot _______________________________________________
  plot.c <- function(model) {
    # recover observed output
    y_obs <- model@sOut

    # compute loocv predictions
    R <- tcrossprod(model@preMats$L)/model@kern@varHyp + diag(model@nugget, nrow = model@n.tr, ncol = model@n.tr)
    Rinv <- solve(R)
    y_pre <- y_obs - diag(Rinv)^(-1) * Rinv %*% y_obs

    # compute LOO statistic
    q2 <- format(getFitness(x.model@model), digits = 3, nsmall = 3)

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
    if (!is.null(gpars$main)) main <- gpars$main else main <- "Model diagnostic by leave-one-out cross-valitation"

    # add points
    plot(y_obs, y_pre, xlim = xlim, ylim = ylim, pch = pch, col = pt.col, bg = pt.bg, cex = pt.cex,
         main = main, xlab = "", ylab = ylab)
    mtext(side = 1, text = xlab, line = 2.0)

    # add reference line
    lines(xlim, ylim, col = line)

    # add Q2 based on x.model stat and fitness
    legend("topleft", legend = paste("Q2loocv =", q2),
           xjust = 0.5, yjust = 0.5, x.intersp = -0.5, y.intersp = 0.3,
           adj = c(0, 0.5), inset = c(.02,.05))
  }
  # ______________________________________________________________________


  # fitness plot _________________________________________________________
  plot.f <- function() {
    y <- x.model@log.success@fitness
    ylab <- strsplit(x.model@stat, split = '\\(')[[1]][1]
    plot(y, ylim = c(0,1), pch = 21, bg = "blue", main = "Fitness of explored models",
         ylab = ylab, yaxt = "n", yaxs = "i")
    mtext(side = 1, text = "Index @log.success", line = 2.0)
    tag <- format(seq(0, 1, length.out = 4), digits = 2, nsmall = 2)
    axis(2, tag, tag, las = 1)
    abline(h = x.model@fitness, lty = 2, col = "green")
    points(x.model@fitness, pch = 24, bg = "green")
    points(x.model@fitness, pch = 25, bg = "green")
    piv.arr <- which(y < 0)
    for (i in seq_along(piv.arr)) {
      arrows(x0 = piv.arr[i], x1 = piv.arr[i], y0 = .01, y1 = .09, length = 0.06, code = 1,
             angle = 25, lwd = 1.5, col = alpha("blue", .7))
    }
  }
  # ______________________________________________________________________

  # save current par state
  opar <- par('mar', 'mfrow')
  on.exit(par(opar))

  # plot
  par(mar = c(3.1, 4.1, 2.5, 2.1), mfrow = c(2,1))
  plot.c(x.model@model)
  plot.f()
}
# ----------------------------------------------------------------------------------------------------------


# plotter for heuristic solution methods (only ant colony for now)
# ----------------------------------------------------------------------------------------------------------
#' @name plotEvol
#' @description This is my description
#' @rdname plotEvol-methods
#' @importFrom graphics lines plot layout legend par arrows boxplot
#' @param x.model An object to predict from.
#' @param ... Further arguments for methods.
#'
#' @author José Betancourt, François Bachoc and Thierry Klein
#' @exportMethod plotEvol
if(!isGeneric("plotEvol")) {setGeneric("plotEvol", function(x.model, ...) standardGeneric("plotEvol"))}

#' @title Prediction Method for the apk Class
#' @name plotEvol
#' @rdname plotEvol-methods
#' @aliases plotEvol,XfunGp-method
setMethod("plotEvol", "XfunGp",
          function(x.model, ...) {
            plotEvol.XfunGp(x.model = x.model, ...)
          })

plotEvol.XfunGp <- function(x.model, ...) {
  # recover fitness list
  all.fitness <- x.model@details$evolution
  n.gen <- length(all.fitness)

  # recover graphic parameters if provided
  gpars <- list(...)
  # <---> limits
  if (!is.null(gpars$xlim)) xlim <- gpars$xlim else xlim <- c(.7, (n.gen + .3))
  if (!is.null(gpars$ylim)) ylim <- gpars$ylim else ylim <- c(0, 1)
  # <---> point colors
  if (!is.null(gpars$pt.col)) pt.col <- gpars$pt.col else pt.col <- alpha("red", .4)
  if (!is.null(gpars$pt.med)) pt.med <- gpars$pt.med else pt.med <- "blue"
  # <---> line colors
  if (!is.null(gpars$col.arr)) col.arr <- gpars$col.arr else col.arr <- "red"
  # <---> text
  if (!is.null(gpars$main)) main <- gpars$main else main <- "ACO learning process"
  if (!is.null(gpars$xlab)) xlab <- gpars$xlab else xlab <- "Generation"
  if (!is.null(gpars$ylab)) ylab <- gpars$ylab else ylab <- x.model@stat

  # save current par state
  opar <- par('mar', 'mfrow')
  on.exit(par(opar))

  # set up layout
  par(mar = c(5.1, 4.1, 4.1, 2.1), mfrow = c(1,1))

  # plot
  plot(1, type = "n", xlim = xlim, ylim = ylim, main = main, xlab = xlab, ylab = ylab, xaxt = "n", yaxs = "i")
  axis(1, 1:n.gen)
  if (n.gen <= 6) {
    for (i in 1:n.gen) {
      ys <- all.fitness[[i]]
      boxplot(at = (i+.04), ys, medcol = pt.med, lty = 1, outline = F,
              boxwex = .2, boxlty = 0, boxfill = alpha(pt.med, .2), axes = F, add = T)
      points(rep((i-.06), length(ys)), ys, pch = "*", col = pt.col)
      arrows(x0 = i, x1 = i, y0 = .01, y1 = .04, length = 0.06, code = 1, angle = 25, lwd = 1.5,
             col = alpha(col.arr, (.3 + .3 * (sum(ys) < 0)/length(ys))))
    }
    abline(v = (1:n.gen)-.1, lty = 3, col = "grey75")
    abline(v = (1:n.gen)+.1, lty = 3, col = "grey75")
  } else {
    for (i in 1:n.gen) {
      ys <- all.fitness[[i]]
      abline(v = i, lty = 3, col = "grey75")
      points(rep(i, length(ys)), ys, pch = "*", col = pt.col, cex = 1.3)
      points(i, median(ys), pch = 21, bg = pt.med, col = pt.med, cex = .7)
      arrows(x0 = i, x1 = i, y0 = .01, y1 = .04, length = 0.06, code = 1, angle = 25, lwd = 1.5,
             col = alpha(col.arr, (.4 + .4 * (sum(ys) < 0)/length(ys))))
    }
  }
  if (T) {
    xpiv <- (xlim[1] + diff(range(xlim))*.77*par()$pin[1]/par()$pin[1])
    ypiv <- (0 + diff(c(0, ylim[2]))*.4*par()$pin[2]/par()$pin[2])
    xarr <- (xlim[1] + diff(range(xlim))*.794*par()$pin[1]/par()$pin[1])
    yarr <- (0 + diff(c(0, ylim[2]))*.213*par()$pin[2]/par()$pin[2])
    xscr0 <- (xlim[1] + diff(range(xlim))*.788*par()$pin[1]/par()$pin[1])
    xscr1 <- (xlim[1] + diff(range(xlim))*.8*par()$pin[1]/par()$pin[1])
    yscr <- (0 + diff(c(0, ylim[2]))*.28*par()$pin[2]/par()$pin[2])
    legend(xpiv, ypiv, legend = c("All models", "Median by gen.", "Points < 0"), pch = c("*", NA, NA),
           col = c("red", NA, NA), lty = c(0, 0, 0), pt.cex = c(1.5, NA, NA), cex = .85, inset = .2)
    if (n.gen <= 6) {
      arrows(x0 = xscr0, x1 = xscr1, y0 = yscr, y1 = yscr,
             length = 0.06, code = 0, angle = 0, lwd = 2.5, col = pt.med)
    } else {
      legend(xpiv, ypiv, legend = c("All models", "Median by gen.", "Points < 0"), pch = c(NA, 21, NA),
             col = c(NA, pt.med, NA), lty = c(0, 0, 0), pt.cex = c(NA, 1, NA), pt.bg = c(NA, pt.med, NA),
             cex = .85, inset = .2, bty = "n")
    }
    arrows(x0 = xarr, x1 = xarr, y0 = yarr, y1 = (yarr + .03),
           length = 0.06, code = 1, angle = 25, lwd = 1.5, col = alpha(col.arr, .7))
  }
}
# ----------------------------------------------------------------------------------------------------------
