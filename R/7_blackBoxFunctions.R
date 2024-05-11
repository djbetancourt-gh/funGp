# ==========================================================================================================
# Analytic black-box functions
# ==========================================================================================================

#' @title Analytic models for the exploration of the funGp package
#'
#' @description Set of analytic functions that take functional
#'     variables as inputs.  Since they run quickly, they can be used
#'     for testing of \pkg{funGp} functionalities as if they were black box
#'     computer models. They cover different situations (number of
#'     scalar inputs and complexity of the inputs-output mathematical
#'     relationship).
#'
#' @details
#'
#' For all the functions, the \eqn{d_s}{ds} scalar inputs
#' \eqn{x_i}{xi} are in the real interval \eqn{[0,\,1]}{[0, 1]} and
#' the \eqn{d_f}{df} functional inputs
#' \eqn{f_i(t_i)}{fi(ti)} are defined on the interval
#' \eqn{[0,\,1]}{[0, 1]}. Expressions for the values are as follows.
#'
#' \itemize{
#'
#' \item\bold{\code{fgp_BB1 }}With \eqn{d_s = 2}{ds = 2} \eqn{d_f = 2}{df = 2}
#'    \preformatted{
#'    x1 * sin(x2) + x1 * mean(f1) - x2^2 * diff(range(f2))}
#'
#' \item\bold{\code{fgp_BB2 }}With \eqn{d_s = 2}{ds = 2} and \eqn{d_f = 2}{df = 2}
#'    \preformatted{
#'    x1 * sin(x2) + mean(exp(x1 * t1) * f1) - x2^2 * mean(f2^2 * t2)}
#'
#' \item\bold{\code{fgp_BB3 }}With  \eqn{d_s = 2}{ds = 2} and \eqn{d_f = 2}{df = 2}
#'   is the first analytical example in Muehlenstaedt et al (2017)
#'    \preformatted{
#'    x1 + 2 * x2 + 4 * mean(t1 * f1) + mean(f2)}
#'
#' \item\bold{\code{fgp_BB4 }}With  \eqn{d_s = 2}{ds = 2} and \eqn{d_f = 2}{df = 2} is the
#'     second analytical example in \emph{preprint} of Muehlenstaedt et al (2017)
#'    \preformatted{
#'    (x2 - (5 / (4 * pi^2)) * x1^2 + (5 / pi) * x1 - 6)^2 +
#'        10 * (1 - (1 / (8 * pi))) * cos(x1) + 10 +
#'        (4 / 3) * pi * (42 * mean(f1 * (1 - t1)) +
#'                        pi * ((x1 + 5) / 5) + 15) * mean(t2 * f2))}
#'
#' \item\bold{\code{fgp_BB5 }}With  \eqn{d_s=2}{ds = 2} and \eqn{d_f=2}{df = 2} is
#'     inspired by the  second analytical example in \emph{final version} of Muehlenstaedt et al (2017)
#'    \preformatted{
#'    (x2 - (5 / (4 * pi^2)) * x1^2 + (5 / pi) * x1 - 6)^2 +
#'        10 * (1 - (1 / (8 * pi))) * cos(x1) + 10 +
#'        (4 / 3) * pi * (42 * mean(15 * f1 * (1 - t1) - 5) +
#'                        pi * ((x1 + 5) / 5) + 15) * mean(15 * t2 * f2))}
#'
#' \item\bold{\code{fgp_BB6 }}With  \eqn{d_s = 2}{ds = 2} and \eqn{d_f = 2}{df = 2}
#'     is inspired by the analytical example in Nanty et al (2016)
#'    \preformatted{
#'    2 * x1^2 + 2 * mean(f1 + t1) + 2 * mean(f2 + t2) + max(f2) + x2}
#'
#' \item\bold{\code{fgp_BB7 }}With \eqn{d_s = 5}{ds = 5} and \eqn{d_f = 2}{df = 2} is
#'    inspired by the second analytical example in \emph{final version} of Muehlenstaedt et al (2017)
#'    \preformatted{
#'    (x2 + 4 * x3 - (5 / (4 * pi^2)) * x1^2 + (5 / pi) * x1 - 6)^2 +
#'        10 * (1 - (1 / (8 * pi))) * cos(x1) * x2^2 * x5^3 + 10 +
#'        (4 / 3) * pi * (42 * sin(x4) * mean(15 * f1 * (1 - t1) - 5) +
#'                        pi * (((x1 * x5 + 5) / 5) + 15) * mean(15 * t2 * f2))}
#' }
#'
#' @param sIn Object with class \code{"matrix"}. The scalar input
#'         points. Variables are arranged by columns and coordinates by rows.
#'
#' @param fIn Object with class \code{"list"}. The functional inputs.
#' Each element of the list must be a matrix containing the set of curves
#' corresponding to one functional input.
#'
#' @param n.tr Object with class \code{"numeric"}. The number of
#'          input points provided and correspondingly, the number of observations
#'          to produce.
#'
#' @section Value:
#' An object of class \code{"matrix"} with the values of the output at the specified input coordinates.
#'
#' @section Note:
#' The functions listed above were used to validate the functionality and stability of this package.
#' Several tests involving \link[=funGp-package]{all main functions, plotters and getters} were run
#' for scalar-input, functional-input and hybrid-input models. In all cases, the output of the functions
#' were correct from the statistical and programmatic perspectives. For an example on the kind of tests
#' performed, the interested user is referred to the introductory funGp manual
#' (\doi{https://doi.org/10.18637/jss.v109.i05}).
#'
#' @references Muehlenstaedt, T., Fruth, J., and Roustant, O. (2017),
#' "Computer experiments with functional inputs and scalar outputs by a norm-based approach".
#' \emph{Statistics and Computing}, \strong{27}, 1083-1097.
#' \href{https://link.springer.com/article/10.1007/s11222-016-9672-z}{[SC]}
#'
#' @references Nanty, S., Helbert, C., Marrel, A., Pérot, N., and Prieur, C. (2016),
#' "Sampling, metamodeling, and sensitivity analysis of numerical simulators with functional stochastic inputs".
#' \emph{SIAM/ASA Journal on Uncertainty Quantification}, \strong{4}(1), 636-659.
#' \doi{10.1137/15M1033319}
#'
#' @name black-boxes
#' @rdname black-boxes
#'
NULL

#' BBK_1
#' @export
#' @rdname black-boxes
fgp_BB1 <- function(sIn, fIn, n.tr){
  sOut <- as.matrix(sapply(t(1:n.tr), function(i){
    x1 <- sIn[i,1]
    x2 <- sIn[i,2]
    f1 <- fIn[[1]][i,]
    f2 <- fIn[[2]][i,]
    as.numeric(x1 * sin(x2) + x1 * mean(f1) - x2^2 * diff(range(f2)))
  }))
  return(sOut)
}


#' BBK_2
#' @export
#' @rdname black-boxes
fgp_BB2 <- function(sIn, fIn, n.tr){
  sOut <- as.matrix(sapply(t(1:n.tr), function(i){
    x1 <- sIn[i,1]
    x2 <- sIn[i,2]
    f1 <- fIn[[1]][i,]
    f2 <- fIn[[2]][i,]
    t1 <- seq(0,1,length = length(f1))
    t2 <- seq(0,1,length = length(f2))
    as.numeric(x1 * sin(x2) + mean(exp(x1 * t1) * f1) - x2^2 * mean(f2^2 * t2))
  }))
  return(sOut)
}


#' MFR_1
#' @export
#' @rdname black-boxes
fgp_BB3 <- function(sIn, fIn, n.tr){
  sOut <- as.matrix(sapply(t(1:n.tr), function(i){
    x1 <- sIn[i,1]
    x2 <- sIn[i,2]
    f1 <- fIn[[1]][i,]
    f2 <- fIn[[2]][i,]
    t1 <- seq(0,1,length = length(f1))
    t2 <- seq(0,1,length = length(f2))
    as.numeric(x1 + 2 * x2 + 4 * mean(t1 * f1) + mean(f2))
  }))
  return(sOut)
}


#' MFR_2p
#' @export
#' @rdname black-boxes
fgp_BB4 <- function(sIn, fIn, n.tr){
  sOut <- as.matrix(sapply(t(1:n.tr), function(i){
    x1 <- sIn[i,1]
    x2 <- sIn[i,2]
    f1 <- fIn[[1]][i,]
    f2 <- fIn[[2]][i,]
    t1 <- seq(0,1,length = length(f1))
    t2 <- seq(0,1,length = length(f2))
    a <- (x2 - (5/(4*pi^2)) * x1^2 + (5/pi) * x1 - 6)^2
    b <- 10 * (1 - (1/(8*pi))) * cos(x1)
    c <- 10
    d1 <- 42 * mean(f1 * (1 - t1))
    d2 <- pi * (((x1 + 5)/5) + 15) * mean(t2 * f2)
    d <- (4/3) * pi * (d1+d2)
    as.numeric(sum(a, b, c, d))
  }))
  return(sOut)
}


#' MFR_2f
#' @export
#' @rdname black-boxes
fgp_BB5 <- function(sIn, fIn, n.tr){
  sOut <- as.matrix(sapply(t(1:n.tr), function(i){
    x1 <- sIn[i,1]
    x2 <- sIn[i,2]
    f1 <- fIn[[1]][i,]
    f2 <- fIn[[2]][i,]
    t1 <- seq(0,1,length = length(f1))
    t2 <- seq(0,1,length = length(f2))
    a <- (x2 - (5/(4*pi^2)) * x1^2 + (5/pi) * x1 - 6)^2
    b <- 10 * (1 - (1/(8*pi))) * cos(x1)
    c <- 10
    d1 <- 42 * mean(15 * f1 * (1 - t1) - 5)
    d2 <- pi * (((x1 + 5)/5) + 15) * mean(15 * t2 * f2)
    d <- (4/3) * pi * (d1+d2)
    as.numeric(sum(a, b, c, d))
  }))
  return(sOut)
}


#' NHMPP
#' @export
#' @rdname black-boxes
fgp_BB6 <- function(sIn, fIn, n.tr){
  sOut <- as.matrix(sapply(t(1:n.tr), function(i){
    x1 <- sIn[i,1]
    x2 <- sIn[i,2]
    f1 <- fIn[[1]][i,]
    f2 <- fIn[[2]][i,]
    t1 <- seq(0,1,length = length(f1))
    t2 <- seq(0,1,length = length(f2))
    as.numeric(2 * x1^2 + 2 * mean(f1 + t1) + 2 * mean(f2 + t2) + max(f2) + x2)
  }))
  return(sOut)
}


#' BBK_7
#' @export
#' @rdname black-boxes
fgp_BB7 <- function(sIn, fIn, n.tr){
  sOut <- as.matrix(sapply(t(1:n.tr), function(i){
    x1 <- sIn[i,1]
    x2 <- sIn[i,2]
    x3 <- sIn[i,3]
    x4 <- sIn[i,4]
    x5 <- sIn[i,5]
    f1 <- fIn[[1]][i,]
    f2 <- fIn[[2]][i,]
    t1 <- seq(0,1,length = length(f1))
    t2 <- seq(0,1,length = length(f2))
    a <- (x2 + 4*x3 - (5/(4*pi^2)) * x1^2 + (5/pi) * x1 - 6)^2
    b <- 10 * (1 - (1/(8*pi))) * cos(x1) * x2^2 * x5^3
    c <- 10
    d1 <- 42 * sin(x4) * mean(15 * f1 * (1 - t1) - 5)
    d2 <- pi * (((x1*x5 + 5)/5) + 15) * mean(15 * t2 * f2)
    d <- (4/3) * pi * (d1+d2)
    as.numeric(sum(a, b, c, d))
  }))
  return(sOut)
}
# ==========================================================================================================
