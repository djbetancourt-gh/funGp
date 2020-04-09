# ==========================================================================================================
# Analytic black-box functions
# ==========================================================================================================
#' @title Analytic black-boxes for the exploration of the funGp package
#' @description Set of black-box analytic functions for the discovering and testing of funGp functionalities.
#'
#' @section Usage:
#' \preformatted{
#' ## Own analytical function 1
#' ## -------------------------
#' ## x1 * sin(x2) + x1 * mean(f1) - x2^2 * diff(range(f2))
#' fgp_BB1(sIn, fIn, n.tr)
#'
#'
#' ## Own analytical function 2
#' ## -------------------------
#' ## x1 * sin(x2) + mean(exp(x1 * t1) * f1) - x2^2 * mean(f2^2 * t2)
#' fgp_BB2(sIn, fIn, n.tr)
#'
#'
#' ## First analytical example in Muehlenstaedt, Fruth & Roustant (2016)
#' ## ------------------------------------------------------------------
#' ## x1 + 2 * x2 + 4 * mean(t1 * f1) + mean(f2)
#' fgp_BB3(sIn, fIn, n.tr)
#'
#'
#' ## Second analytical example in preprint of Muehlenstaedt, Fruth & Roustant (2016)
#' ## -------------------------------------------------------------------------------
#' ## a = (x2 - (5/(4*pi^2)) * x1^2 + (5/pi) * x1 - 6)^2
#' ## b = 10 * (1 - (1/(8*pi))) * cos(x1)
#' ## c = 10
#' ## d = (4/3) * pi * (42 * mean(f1*(1-t1)) + pi * (((x1+5)/5) + 15) * mean(t2*f2))
#' ## a + b + c + d
#' fgp_BB4(sIn, fIn, n.tr)
#'
#'
#' ## Second analytical example in final version of Muehlenstaedt, Fruth & Roustant (2016)
#' ## ------------------------------------------------------------------------------------
#' ## a = (x2 - (5/(4*pi^2)) * x1^2 + (5/pi) * x1 - 6)^2
#' ## b = 10 * (1 - (1/(8*pi))) * cos(x1)
#' ## c = 10
#' ## d <- (4/3) * pi * (42 * mean(15*f1*(1-t1)-5) + pi * (((x1+5)/5) + 15) * mean(15*t2*f2))
#' ## a + b + c + d
#' fgp_BB5(sIn, fIn, n.tr)
#'
#'
#' ## Inspired by the analytical example in Nanty, Helbert, Marrel, Pérot, Prieur (2016)
#' ## ----------------------------------------------------------------------------------
#' ## 2 * x1^2 + 2 * mean(f1 + t1) + 2 * mean(f2 + t2) + max(f2) + x2
#' fgp_BB6(sIn, fIn, n.tr)
#'
#'
#' ## Inspired by the second analytical example in final version of Muehlenstaedt et al (2016)
#' ## ----------------------------------------------------------------------------------------
#' ## a = (x2 + 4*x3 - (5/(4*pi^2)) * x1^2 + (5/pi) * x1 - 6)^2
#' ## b = 10 * (1 - (1/(8*pi))) * cos(x1) * x2^2 * x5^3
#' ## c = 10
#' ## d <- (4/3) * pi * (42 * sin(x4) * mean(15*f1*(1-t1)-5) +
#'                                           pi * (((x1*x5+5)/5) + 15) * mean(15*t2*f2))
#' ## a + b + c + d
#' fgp_BB7(sIn, fIn, n.tr)
#' }
#'
#' @section Arguments:
#' \strong{*}\emph{sIn}: Object of class \code{"matrix"}. The scalar input points. Variables are arranged
#'   by columns and coordinates by rows. \cr\cr
#' \strong{*}\emph{fIn}: Object of class \code{"list"}. The functional input points. Each element of the
#'   list contains a functional input in the form of a matrix. In each matrix, curves representing
#'   functional coordinates are arranged by rows. \cr\cr
#' \strong{*}\emph{n.tr}: Object of class \code{"numeric"}. The number of input points provided and
#'   correspondingly, the number of observations to produce.
#'
#' @section Value:
#' An object of class \code{"matrix"} with the values of the output at the specified input coordinates.
#'
#' @section Note:
#' The functions listed above were used to validate the functionality and stability of this package.
#' Several tests involving \link[=funGp-package]{all main functions, plotters and getters} were run
#' for scalar-input, functional-input and hybrid-input models. In all cases, the output of the functions
#' were correct from the statistical and programmatic perspectives. For an example on the kind of tests
#' performed, the interested user is referred to
#' \href{https://hal.archives-ouvertes.fr/hal-02536624}{the introductory funGp manual}.
#'
#' @references Muehlenstaedt, T., Fruth, J., and Roustant, O. (2017),
#' "Computer experiments with functional inputs and scalar outputs by a norm-based approach".
#' \emph{Statistics and Computing}, \strong{27}, 1083-1097.
#' \href{https://link.springer.com/article/10.1007/s11222-016-9672-z}{[SC]}
#'
#' @references Nanty, S., Helbert, C., Marrel, A., Pérot, N., and Prieur, C. (2016),
#' "Sampling, metamodeling, and sensitivity analysis of numerical simulators with functional stochastic inputs".
#' \emph{SIAM/ASA Journal on Uncertainty Quantification}, \strong{4}(1), 636-659.
#' \href{https://epubs.siam.org/doi/10.1137/15M1033319}{[SA-JUQ]}
#'
#' @docType methods
#' @name black-boxes
#' @rdname black-boxes
NULL
# ==========================================================================================================

#' BBK_1
#' @export
#' @keywords internal
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
#' @keywords internal
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
#' @keywords internal
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
#' @keywords internal
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
#' @keywords internal
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
#' @keywords internal
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
#' @keywords internal
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
