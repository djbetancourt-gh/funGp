#' A tester
#'
#' Allows to make tests
#' @importFrom graphics abline
#' @export
tester <- function(){
  # example funGp()
  # =========================================
  # generating input data for training
  set.seed(100)
  n.tr <- 25
  sIn <- expand.grid(x1 = seq(0,1,length = sqrt(n.tr)), x2 = seq(0,1,length = sqrt(n.tr)))
  fIn <- list(f1 = matrix(runif(n.tr*10), ncol = 10), f2 = matrix(runif(n.tr*22), ncol = 22))

  # generating output data for training
  sOut <- as.matrix(sapply(t(1:n.tr), function(i){
    x1 <- sIn[i,1]
    x2 <- sIn[i,2]
    f1 <- fIn[[1]][i,]
    f2 <- fIn[[2]][i,]
    t1 <- seq(0,1,length = length(f1))
    t2 <- seq(0,1,length = length(f2))
    as.numeric(x1 + 2 * x2 + 4 * mean(t1 * f1) + mean(f2))
  }))

  # creating a funGp model
  m1 <- funGp(sIn = sIn, fIn = fIn, sOut = sOut)

  # plotting the model
  plotLOO(m1)
  # =========================================


  # example show()
  # =========================================
  # generating input data for training
  set.seed(100)
  n.tr <- 25
  sIn <- expand.grid(x1 = seq(0,1,length = sqrt(n.tr)), x2 = seq(0,1,length = sqrt(n.tr)))
  fIn <- list(f1 = matrix(runif(n.tr*10), ncol = 10), f2 = matrix(runif(n.tr*22), ncol = 22))

  # generating output data for training
  sOut <- as.matrix(sapply(t(1:n.tr), function(i){
    x1 <- sIn[i,1]
    x2 <- sIn[i,2]
    f1 <- fIn[[1]][i,]
    f2 <- fIn[[2]][i,]
    t1 <- seq(0,1,length = length(f1))
    t2 <- seq(0,1,length = length(f2))
    as.numeric(x1 + 2 * x2 + 4 * mean(t1 * f1) + mean(f2))
  }))

  # creating three types of funGp models: scalar-input, functional-input and hybrid-input
  ms <- funGp(sIn = sIn, sOut = sOut)
  mf <- funGp(fIn = fIn, sOut = sOut)
  msf <- funGp(sIn = sIn, fIn = fIn, sOut = sOut)

  # printing each model
  show(ms)
  show(mf)
  show(msf)

  # printing the kernel structure of each model
  show(ms@kern)
  show(mf@kern)
  show(msf@kern)

  # printing the projection structure of each model
  show(ms@proj)
  show(mf@proj)
  show(msf@proj)

  # the show method is also called when the name of any object of class funGp, funGpKern or
  # funGpProj is typed followed by enter. For instance, 'show(ms)' produces the same output as 'ms'.
  # =========================================


  # example predict()
  # =========================================
  # defining blackbox to generate observations
  bbox <- as.matrix(sapply(t(1:n.tr), function(i){
    x1 <- sIn[i,1]
    x2 <- sIn[i,2]
    f1 <- fIn[[1]][i,]
    f2 <- fIn[[2]][i,]
    t1 <- seq(0,1,length = length(f1))
    t2 <- seq(0,1,length = length(f2))
    as.numeric(x1 + 2 * x2 + 4 * mean(t1 * f1) + mean(f2))
  }))

  # generating input data for training
  set.seed(100)
  n.tr <- 25
  sIn <- expand.grid(x1 = seq(0,1,length = sqrt(n.tr)), x2 = seq(0,1,length = sqrt(n.tr)))
  fIn <- list(f1 = matrix(runif(n.tr*10), ncol = 10), f2 = matrix(runif(n.tr*22), ncol = 22))

  # generating output data for training
  sOut <- ft(sIn, fIn, n.tr)

  # creating a funGp model
  m1 <- funGp(sIn = sIn, fIn = fIn, sOut = sOut)

  # generating input data for prediction
  n.pr <- 100
  sIn.pr <- expand.grid(x1 = seq(0,1,length = sqrt(n.pr)), x2 = seq(0,1,length = sqrt(n.pr)))
  fIn.pr <- list(f1 = matrix(runif(n.pr*10), ncol = 10), f2 = matrix(runif(n.pr*22), ncol = 22))

  # making predictions
  m1.preds <- predict(m1, sIn.pr = sIn.pr, fIn.pr = fIn.pr)

  # plotting predictions
  plotPreds(m1, preds = m1.preds)

  # It is also possible to compare against true output values
  sOut.pr <- ft(sIn.pr, fIn.pr, n.pr)
  plotPreds(m1, m1.preds, sOut.pr)
  # =========================================


  # example simulate()
  # =========================================
  # generating input data for training
  set.seed(100)
  n.tr <- 25
  sIn <- expand.grid(x1 = seq(0,1,length = sqrt(n.tr)), x2 = seq(0,1,length = sqrt(n.tr)))
  fIn <- list(f1 = matrix(runif(n.tr*10), ncol = 10), f2 = matrix(runif(n.tr*22), ncol = 22))

  # generating output data for training
  sOut <- as.matrix(sapply(t(1:n.tr), function(i){
    x1 <- sIn[i,1]
    x2 <- sIn[i,2]
    f1 <- fIn[[1]][i,]
    f2 <- fIn[[2]][i,]
    t1 <- seq(0,1,length = length(f1))
    t2 <- seq(0,1,length = length(f2))
    as.numeric(x1 + 2 * x2 + 4 * mean(t1 * f1) + mean(f2))
  }))

  # creating a funGp model
  m1 <- funGp(sIn = sIn, fIn = fIn, sOut = sOut)

  # generating input data for simulation
  set.seed(100)
  n.sm <- 100
  sIn.sm <- expand.grid(x1 = seq(0,1,length = sqrt(n.sm)), x2 = seq(0,1,length = sqrt(n.sm)))
  fIn.sm <- list(f1 = matrix(runif(n.sm*10), ncol = 10), f2 = matrix(runif(n.sm*22), ncol = 22))

  # making light simulations
  m1.sims_l <- simulate(m1, nsim = 10, sIn.sm = sIn.sm, fIn.sm = fIn.sm)

  # plotting light simulations
  plotSims(m1, m1.sims_l)

  # making full simulations
  m1.sims_f <- simulate(m1, nsim = 10, sIn.sm = sIn.sm, fIn.sm = fIn.sm, detail = "full")

  # plotting full simulations in full mode
  plotSims(m1, m1.sims_f)

  # plotting full simulations in light mode
  plotSims(m1, m1.sims_f, detail = "light")
  # =========================================


  # example update()
  # =========================================
  # defining blackbox to generate observations
  bbox <- as.matrix(sapply(t(1:n.tr), function(i){
    x1 <- sIn[i,1]
    x2 <- sIn[i,2]
    f1 <- fIn[[1]][i,]
    f2 <- fIn[[2]][i,]
    t1 <- seq(0,1,length = length(f1))
    t2 <- seq(0,1,length = length(f2))
    as.numeric(x1 + 2 * x2 + 4 * mean(t1 * f1) + mean(f2))
  }))

  # generating input data for training
  set.seed(100)
  n.tr <- 25
  sIn <- expand.grid(x1 = seq(0,1,length = sqrt(n.tr)), x2 = seq(0,1,length = sqrt(n.tr)))
  fIn <- list(f1 = matrix(runif(n.tr*10), ncol = 10), f2 = matrix(runif(n.tr*22), ncol = 22))

  # generating output data for training
  sOut <- ft(sIn, fIn, n.tr)

  # creating a funGp model
  m1 <- funGp(sIn = sIn, fIn = fIn, sOut = sOut)
  show(m1)
  plotLOO(m1)

  # generating complementary input data for updating
  n.nw <- 3
  sIn.nw <- matrix(runif(n.nw * m1@ds), nrow = n.nw)
  fIn.nw <- list(f1 = matrix(runif(n.nw*10), ncol = 10), f2 = matrix(runif(n.nw*22), ncol = 22))

  # generating complementary output data for updating
  sOut.nw <- ft(sIn.nw, fIn.nw, n.nw)

  # updating the model
  m1up1 <- update(m1, sIn.nw = sIn.nw, fIn.nw = fIn.nw, sOut.nw = sOut.nw)
  show(m1up1)
  plotLOO(m1up1)

  # generating substituting input data for updating
  n.sb <- 2
  sIn.sb <- matrix(runif(n.sb * m1@ds), nrow = n.sb)
  fIn.sb <- list(f1 = matrix(runif(n.sb*10), ncol = 10), f2 = matrix(runif(n.sb*22), ncol = 22))

  # generating substituting output data for updating
  sOut.sb <- ft(sIn.sb, fIn.sb, n.sb)

  # generating indices for substitution
  n.tot <- (m1up1@n.tr+n.nw) # should be replaced when the slot n.tot is added to the model !!!!!!!!!!!!!
  ind.sb <- sample(1:n.tot, n.sb)

  # updating the model
  m1up2 <- update(m1up1, sIn.sb = sIn.sb, fIn.sb = fIn.sb, sOut.sb = sOut.sb, ind.sb = ind.sb)
  show(m1up2)
  plotLOO(m1up2)
  # =========================================



  # choosing test function
  #xx ft <- test1 # R: both model optimizations get lost almost always.
  #x ft <- test2 # R: ok for n.tr = 4, 9. from 16 to 49 it gets difficult to get good hypers, specially for the scalar meta.
  ft <- test3 # R: ok for 4, 9, 225. This last with some difficulties for the scalar case.
  #x ft <- test4 # R: ok for 4, 9, 16. Then 25, 36, 49 work some times.
  #xx ft <- test5 # R: ok for 4, 9. Hardly 16.
  #x ft <- test6 # R: ok for 4, 9, 16, 25, 36. Bit hard for 49.

  # generating input data for training
  set.seed(100)
  n.tr <- 5^2 # should be a number with exact sqrt: {4, 9, 16, 25, 36, 49, 64, 81, 100, 121, 144, 169, 196, 225}
  sIn <- as.matrix(expand.grid(x1 = seq(0,1,length = sqrt(n.tr)), x2 = seq(0,1,length = sqrt(n.tr))))
  fIn <- list(f1 = matrix(runif(n.tr*10), ncol = 10), matrix(runif(n.tr*22), ncol = 22))

  # generating output data for training
  sOut <- ft(sIn, fIn, n.tr)

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

  # generating input data for prediction
  set.seed(100)
  n.pr <- 100
  sIn.pr <- as.matrix(expand.grid(x1 = seq(0,1,length = sqrt(n.pr)), x2 = seq(0,1,length = sqrt(n.pr))))
  fIn.pr <- list(f1 = matrix(runif(n.pr*10), ncol = 10), matrix(runif(n.pr*22), ncol = 22))

  # making predictions with each model
  ms.preds <- predict(ms, sIn.pr = sIn.pr)
  mf.preds <- predict(mf, fIn.pr = fIn.pr)
  msf.preds <- predict(msf, sIn.pr = sIn.pr, fIn.pr = fIn.pr)

  # plotting predictions
  plotPreds(ms, preds = ms.preds)
  plotPreds(mf, preds = mf.preds)
  plotPreds(msf, preds = msf.preds)

  # generating output data for prediction
  sOut.pr <- ft(sIn.pr, fIn.pr, n.pr)

  # comparing against true output values
  plotPreds(ms, ms.preds, sOut.pr)
  plotPreds(mf, mf.preds, sOut.pr)
  plotPreds(msf, msf.preds, sOut.pr)
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
