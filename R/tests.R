#' A tester
#'
#' Allows to make tests
#' @importFrom graphics abline
#' @export
tester <- function(){
  # generating input data for training
  n.tr <- 14^2 # should be a number with exact sqrt: {4, 9, 16, 25, 36, 49, 64, 81, 100, 121, 144, 169, 196, 225}
  sIn <- as.matrix(expand.grid(x1 = seq(0,1,length = sqrt(n.tr)), x2 = seq(0,1,length = sqrt(n.tr))))
  fIn <- list(f1 = matrix(runif(n.tr*10), ncol = 10), matrix(runif(n.tr*22), ncol = 22))

  # generating output data for training
  # sOut <- test1(sIn, fIn, n.tr) # R: both model optimizations get lost almost always.
  # sOut <- test2(sIn, fIn, n.tr) # R: ok for n.tr = 4, 9. from 16 to 49 it gets difficult to get good hypers, specially for the scalar meta.
  # sOut <- test3(sIn, fIn, n.tr) # R: ok for 4, 9, 225. This last with some difficulties for the scalar case.
  # sOut <- test4(sIn, fIn, n.tr) # R: ok for 4, 9, 16. Then 25, 36, 49 work some times.
  # sOut <- test5(sIn, fIn, n.tr) # R: ok for 4, 9. Hardly 16.
  sOut <- test6(sIn, fIn, n.tr) # R: ok for 4, 9, 16, 25, 36. Bit hard for 49.

  # creating the three types of model
  ms <- funGp(sIn = sIn, sOut = sOut)
  mf <- funGp(fIn = fIn, sOut = sOut)
  msf <- funGp(sIn = sIn, fIn = fIn, sOut = sOut)

  # plotting the models
  layout(matrix(1:3, nrow = 3))
  plotLOO(ms)
  abline(h = 0, col = "green")
  plotLOO(mf)
  abline(h = 0, col = "green")
  plotLOO(msf)
}

# O_1
# Our own analytical function 1
# @export
test1 <- function(sIn, fIn, n.tr){
  sOut <- as.matrix(sapply(t(1:n.tr), function(i){
    x1 <- sIn[i,1]
    x2 <- sIn[i,2]
    f1 <- fIn[[1]][i,]
    f2 <- fIn[[2]][i,]
    as.numeric(x1 * sin(x2) + x1 * mean(f1) - x2^2 * diff(range(f2)))
  }))
  return(sOut)
}

# O_2
# Our own analytical function 2
# @export
test2 <- function(sIn, fIn, n.tr){
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

# MFR_1
# First analytical example in Muehlenstaedt, Fruth & Roustant (2016)
# @export
test3 <- function(sIn, fIn, n.tr){
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
test4 <- function(sIn, fIn, n.tr){
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
test5 <- function(sIn, fIn, n.tr){
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
test6 <- function(sIn, fIn, n.tr){
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



# A tester
#
# Allows to make tests
# @export
testerKm <- function(){
  # generating input data for training
  # set.seed(100)
  # n.tr <- 100 # should be a number with exact sqrt: {4, 9, 16, 25, 36, 49, 64, 81, 100, 121, 144, 169, 196, 225}
  # sIn <- as.matrix(expand.grid(x1 = seq(0,1,length = sqrt(n.tr)), x2 = seq(0,1,length = sqrt(n.tr))))
  # fIn <- list(f1 = matrix(runif(n.tr*10), ncol = 10), matrix(runif(n.tr*22), ncol = 22))

  # generating output data for training
  # sOut <- test1(sIn, fIn, n.tr) # R: ok for 4, 9. Hardly 16. Both model optimizations get lost almost always for the others.
  # sOut <- test2(sIn, fIn, n.tr) # R: ok for n.tr = 4, 9. from 16 to 49 it gets difficult to get good hypers, specially for the scalar meta.
  # sOut <- test3(sIn, fIn, n.tr) # R: ok for 4, 9, 225. This last with some difficulties for the scalar case.
  # sOut <- test4(sIn, fIn, n.tr) # R: ok for 4, 9, 16. Then 25, 36, 49 work some times.
  # sOut <- test5(sIn, fIn, n.tr) # R: ok for 4, 9. Hardly 16.
  # sOut <- test6(sIn, fIn, n.tr) # R: ok for 4, 9, 16, 25, 36. Bit hard for 49.

  # creating the two types of model
  # ms <- funGp(sIn = sIn, sOut = sOut)
  # mk <- km(coef.trend = 0, design = data.frame(sIn), response = sOut, gr = FALSE)

  # ms <- funGp(sIn = sIn, sOut = sOut, n.presample = 1)
  # mk <- km(coef.trend = 0, design = data.frame(sIn), response = sOut, gr = FALSE, control = list(pop.size = 1))

  # plotting the models
  # layout(matrix(1:2, nrow = 2))
  # plotLOO(ms)
  # abline(h = 0, col = "green")
  # plot(mk)
  # plotLOO(msf)
}
