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
#' ## Function inspired by the analytical example in Nanty, Helbert, Marrel, Pérot, Prieur (2016)
#' ## -------------------------------------------------------------------------------------------
#' ## 2 * x1^2 + 2 * mean(f1 + t1) + 2 * mean(f2 + t2) + max(f2) + x2
#' fgp_BB6(sIn, fIn, n.tr)
#' }
#'
#' @section Arguments:
#' \preformatted{
#' sIn  fill!!!
#' fIn  fill!!!
#' n.tr fill!!!
#' }
#'
#' @section Note:
#' The functions listed above were used to validate the functionality and stability of this package.
#' Several tests involving \link[=funGp-package]{all main functions, plotters and getters} were run
#' for scalar-input, functional-input and hybrid-input models. In all cases, the output of the functions
#' were correct from the statistical and programming perspectives. For an example on the kind of tests
#' performed, the interested user is referred to
#' \href{https://drive.google.com/open?id=0B5dR1D0AmTvsb01qcmJ0UWg1TlE}{the introductory funGp manual}.
#'
#' @docType methods
#' @name black-boxes
#' @rdname black-boxes
NULL

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


#' # MFR_1
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

# MFR_2p
# Second analytical example in preprint of Muehlenstaedt, Fruth & Roustant (2016)
# @export
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

# MFR_2f
# Second analytical example in final version of Muehlenstaedt, Fruth & Roustant (2016)
# @export
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

# NHMPP
# Function inspired by the analytical example in Nanty, Helbert, Marrel, Pérot, Prieur (2016)
# @export
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
