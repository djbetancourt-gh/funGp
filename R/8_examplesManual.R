# Creating a model
# ----------------------------------------------------------------------------------------------------------
section_1.1_create <- function(){
  # generating input data for training
  set.seed(100)
  n.tr <- 25
  sIn <- expand.grid(x1 = seq(0,1,length = sqrt(n.tr)), x2 = seq(0,1,length = sqrt(n.tr)))
  fIn <- list(f1 = matrix(runif(n.tr*10), ncol = 10), f2 = matrix(runif(n.tr*22), ncol = 22))

  # generating output data for training
  sOut <- fgp_BB3(sIn, fIn, n.tr)

  # creating a funGp model
  m1 <- funGp(sIn = sIn, fIn = fIn, sOut = sOut)

  # printing the model
  m1 # equivalent to show(m1)

  # plotting the model
  plotLOO(m1)
}
# ----------------------------------------------------------------------------------------------------------


# Making predictions
# ----------------------------------------------------------------------------------------------------------
section_1.2_predict <- function(){
  # building the model
  set.seed(100)
  n.tr <- 25
  sIn <- expand.grid(x1 = seq(0,1,length = sqrt(n.tr)), x2 = seq(0,1,length = sqrt(n.tr)))
  fIn <- list(f1 = matrix(runif(n.tr*10), ncol = 10), f2 = matrix(runif(n.tr*22), ncol = 22))
  sOut <- fgp_BB3(sIn, fIn, n.tr)
  m1 <- funGp(sIn = sIn, fIn = fIn, sOut = sOut)

  # generating input data for prediction
  set.seed(100)
  n.pr <- 100
  sIn.pr <- expand.grid(x1 = seq(0,1,length = sqrt(n.pr)), x2 = seq(0,1,length = sqrt(n.pr)))
  fIn.pr <- list(f1 = matrix(runif(n.pr*10), ncol = 10), f2 = matrix(runif(n.pr*22), ncol = 22))

  # making predictions
  m1.preds <- predict(m1, sIn.pr = sIn.pr, fIn.pr = fIn.pr)

  # checking content of the list
  summary(m1.preds)

  # plotting predictions
  plotPreds(m1, preds = m1.preds)

  # It is also possible to compare against true output values
  sOut.pr <- fgp_BB3(sIn.pr, fIn.pr, n.pr)
  plotPreds(m1, m1.preds, sOut.pr)

  # making full predictions
  m1.preds_f <- predict(m1, sIn.pr = sIn.pr, fIn.pr = fIn.pr, detail = "full")

  # checking content of the list
  summary(m1.preds_f)

  # plotting full predictions without true output
  plotPreds(m1, preds = m1.preds_f)

  # plotting full predictions along with true output
  plotPreds(m1, m1.preds_f, sOut.pr)
}
# ----------------------------------------------------------------------------------------------------------


# Making simulations
section_1.3_simulate <- function(){
  # building the model
  set.seed(100)
  n.tr <- 25
  sIn <- expand.grid(x1 = seq(0,1,length = sqrt(n.tr)), x2 = seq(0,1,length = sqrt(n.tr)))
  fIn <- list(f1 = matrix(runif(n.tr*10), ncol = 10), f2 = matrix(runif(n.tr*22), ncol = 22))
  sOut <- fgp_BB3(sIn, fIn, n.tr)
  m1 <- funGp(sIn = sIn, fIn = fIn, sOut = sOut)

  # generating input data for simulation
  set.seed(100)
  n.sm <- 100
  sIn.sm <- as.matrix(expand.grid(x1 = seq(0,1,length = sqrt(n.sm)), x2 = seq(0,1,length = sqrt(n.sm))))
  fIn.sm <- list(f1 = matrix(runif(n.sm*10), ncol = 10), matrix(runif(n.sm*22), ncol = 22))

  # making light simulations
  m1.sims_l <- simulate(m1, nsim = 10, sIn.sm = sIn.sm, fIn.sm = fIn.sm)

  # plotting light simulations
  plotSims(m1, m1.sims_l)

  # making full simulations
  m1.sims_f <- simulate(m1, nsim = 10, sIn.sm = sIn.sm, fIn.sm = fIn.sm, detail = "full")

  # checking content of the list
  summary(m1.sims_f)

  # plotting full simulations in full mode
  plotSims(m1, m1.sims_f)

  # plotting full simulations in light mode
  plotSims(m1, m1.sims_f, detail = "light")
}
# ----------------------------------------------------------------------------------------------------------


# Making updates: deletion
# ----------------------------------------------------------------------------------------------------------
section_1.4.a_update_delete <- function(){
  # building the initial model
  set.seed(100)
  n.tr <- 25
  sIn <- expand.grid(x1 = seq(0,1,length = sqrt(n.tr)), x2 = seq(0,1,length = sqrt(n.tr)))
  fIn <- list(f1 = matrix(runif(n.tr*10), ncol = 10), f2 = matrix(runif(n.tr*22), ncol = 22))
  sOut <- fgp_BB3(sIn, fIn, n.tr)
  m1 <- funGp(sIn = sIn, fIn = fIn, sOut = sOut)

  # generating deletion indices
  n.dl <- 2
  ind.dl <- sample(1:m1@n.tot, n.dl)

  # updating the model
  m1up <- update(m1, ind.dl = ind.dl)

  # checking model updates
  res <- matrix(c("m1", "m1up", m1@n.tr, m1up@n.tr, m1@n.tot, m1up@n.tot,
                  nrow(m1@sIn), nrow(m1up@sIn), nrow(m1@fIn[[1]]), nrow(m1up@fIn[[1]]),
                  nrow(m1@sOut), nrow(m1up@sOut)), nrow = 2)
  colnames(res) <- c("model", "n.tr", "n.tot", "nrow(sIn)", "nrow(fIn)", "nrow(sOut)")
  # knitr::kable(res)

  # printing m1's kern structure
  m1@kern

  # printing m1up's kern structure
  m1up@kern

  # generating data for prediction
  n.pr <- 100
  sIn.pr <- expand.grid(x1 = seq(0,1,length = sqrt(n.pr)), x2 = seq(0,1,length = sqrt(n.pr)))
  fIn.pr <- list(f1 = matrix(runif(n.pr*10), ncol = 10), f2 = matrix(runif(n.pr*22), ncol = 22))
  sOut.pr <- fgp_BB3(sIn.pr, fIn.pr, n.pr)

  # predicting with both models
  m1.preds <- predict(m1, sIn.pr = sIn.pr, fIn.pr = fIn.pr, detail = "full")
  m1up.preds <- predict(m1up, sIn.pr = sIn.pr, fIn.pr = fIn.pr, detail = "full")

  # comparing models' prediction plots
  plotPreds(m1, preds = m1.preds, sOut.pr, justCal = T)
  points(sOut.pr, m1up.preds$mean, pch = "*", cex = 1.5)
  legend("topleft", legend = c("m1", "m1up"), pch = c(21, NA), col = c("red", "red"), cex = 1,
         pt.bg = c("red", NA), inset = .02)
  legend("topleft", legend = c("m1", "m1up"), pch = c(NA, "*"), col = c(NA, "black"), cex = 1,
         pt.cex = 2, bg = "transparent", inset = .02)
}
# ----------------------------------------------------------------------------------------------------------


# Making updates: substitution
# ----------------------------------------------------------------------------------------------------------
section_1.4.b_update_substitute <- function(){
  # building the initial model
  set.seed(100)
  n.tr <- 25
  sIn <- expand.grid(x1 = seq(0,1,length = sqrt(n.tr)), x2 = seq(0,1,length = sqrt(n.tr)))
  fIn <- list(f1 = matrix(runif(n.tr*10), ncol = 10), f2 = matrix(runif(n.tr*22), ncol = 22))
  sOut <- fgp_BB3(sIn, fIn, n.tr)
  m1 <- funGp(sIn = sIn, fIn = fIn, sOut = sOut)

  # generating substituting input data for updating
  n.sb <- 2
  sIn.sb <- matrix(runif(n.sb * m1@ds), nrow = n.sb)
  fIn.sb <- list(f1 = matrix(runif(n.sb*10), ncol = 10), f2 = matrix(runif(n.sb*22), ncol = 22))

  # generating substituting output data for updating
  sOut.sb <- fgp_BB3(sIn.sb, fIn.sb, n.sb)

  # generating indices for substitution
  ind.sb <- sample(1:(m1@n.tot), n.sb)

  # updating only the scalar inputs
  m1up1 <- update(m1, sIn.sb = sIn.sb, ind.sb = ind.sb)

  # updating only the output
  m1up2 <- update(m1, sOut.sb = sOut.sb, ind.sb = ind.sb)

  # updating the scalar inputs and the output
  m1up3 <- update(m1, sIn.sb = sIn.sb, sOut.sb = sOut.sb, ind.sb = ind.sb)

  # updating all, the scalar inputs, functional inputs and the output
  m1up4 <- update(m1, sIn.sb = sIn.sb, fIn.sb = fIn.sb, sOut.sb = sOut.sb, ind.sb = ind.sb)

  # veryfing substitutions
  all(m1up1@sIn[ind.sb,] == sIn.sb)
  all(m1up2@sOut[ind.sb,] == sOut.sb)
  all(m1up3@sIn[ind.sb,] == sIn.sb, m1up3@sOut[ind.sb,] == sOut.sb)
  all(m1up4@sIn[ind.sb,] == sIn.sb, all(sapply(mapply(function(M, m) M[ind.sb,] == m, m1up4@fIn, fIn.sb), all)),
      m1up4@sOut[ind.sb,] == sOut.sb)
}
# ----------------------------------------------------------------------------------------------------------


# Making updates: addition
# ----------------------------------------------------------------------------------------------------------
section_1.4.c_update_addition <- function(){
  # building the initial model
  set.seed(100)
  n.tr <- 25
  sIn <- expand.grid(x1 = seq(0,1,length = sqrt(n.tr)), x2 = seq(0,1,length = sqrt(n.tr)))
  fIn <- list(f1 = matrix(runif(n.tr*10), ncol = 10), f2 = matrix(runif(n.tr*22), ncol = 22))
  sOut <- fgp_BB3(sIn, fIn, n.tr)
  m1 <- funGp(sIn = sIn, fIn = fIn, sOut = sOut)

  # generating complementary input data for updating
  n.nw <- 3
  sIn.nw <- matrix(runif(n.nw * m1@ds), nrow = n.nw)
  fIn.nw <- list(f1 = matrix(runif(n.nw*10), ncol = 10), f2 = matrix(runif(n.nw*22), ncol = 22))

  # generating complementary output data for updating
  sOut.nw <- fgp_BB3(sIn.nw, fIn.nw, n.nw)

  # updating the model
  m1up <- update(m1, sIn.nw = sIn.nw, fIn.nw = fIn.nw, sOut.nw = sOut.nw)

  # checking model updates
  res <- matrix(c("m1", "m1up", m1@n.tr, m1up@n.tr, m1@n.tot, m1up@n.tot,
                  nrow(m1@sIn), nrow(m1up@sIn), nrow(m1@fIn[[1]]), nrow(m1up@fIn[[1]]),
                  nrow(m1@sOut), nrow(m1up@sOut)), nrow = 2)
  colnames(res) <- c("model", "n.tr", "n.tot", "nrow(sIn)", "nrow(fIn)", "nrow(sOut)")
  # knitr::kable(res)
}
# ----------------------------------------------------------------------------------------------------------


# Making updates: substitution of hyperparameters
# ----------------------------------------------------------------------------------------------------------
section_1.4.d_update_subshypers <- function(){
  # building the initial model
  set.seed(100)
  n.tr <- 25
  sIn <- expand.grid(x1 = seq(0,1,length = sqrt(n.tr)), x2 = seq(0,1,length = sqrt(n.tr)))
  fIn <- list(f1 = matrix(runif(n.tr*10), ncol = 10), f2 = matrix(runif(n.tr*22), ncol = 22))
  sOut <- fgp_BB3(sIn, fIn, n.tr)
  m1 <- funGp(sIn = sIn, fIn = fIn, sOut = sOut)

  # defining var hyperparameter for substitution
  var.sb <- 3
  ls_s.sb <- c(2.44, 1.15)
  ls_f.sb <- c(5.83, 4.12)

  # updating the model
  m1up <- update(m1, var.sb = var.sb)
  m1up <- update(m1, ls_s.sb = ls_s.sb)
  m1up <- update(m1, ls_f.sb = ls_f.sb)
}


# Making updates: re-estimation of hyperparameters
# ----------------------------------------------------------------------------------------------------------
# building the initial model
section_1.4.e_update_reeshypers <- function(){
  set.seed(100)
  n.tr <- 25
  sIn <- expand.grid(x1 = seq(0,1,length = sqrt(n.tr)), x2 = seq(0,1,length = sqrt(n.tr)))
  fIn <- list(f1 = matrix(runif(n.tr*10), ncol = 10), f2 = matrix(runif(n.tr*22), ncol = 22))
  sOut <- fgp_BB3(sIn, fIn, n.tr)
  m1 <- funGp(sIn = sIn, fIn = fIn, sOut = sOut)

  # defining var hyperparameter for substitution
  var.sb <- 3
  ls_s.sb <- c(2.44, 1.15)
  ls_f.sb <- c(5.83, 4.12)

  # updating the model
  m1 <- update(m1, var.sb = var.sb)
  m1up <- update(m1, var.re = T)
  m1 <- update(m1, ls_s.sb = ls_s.sb)
  m1up <- update(m1, ls_s.re = T)
  m1 <- update(m1, ls_f.sb = ls_f.sb)
  m1up <- update(m1, ls_f.re = T)
  m1 <- update(m1, var.sb = var.sb, ls_s.sb = ls_s.sb, ls_f.sb = ls_f.sb)
  m1up <- update(m1, var.re = T, ls_s.re = T, ls_f.re = T)
}
# ----------------------------------------------------------------------------------------------------------

#' A tester
#'
#' Allows to make tests
#' @importFrom graphics abline points
examplesManual <- function(){

  # Making updates: deletion with different types of models
  # =========================================
  # building the initial model
  set.seed(100)
  n.tr <- 25
  sIn <- expand.grid(x1 = seq(0,1,length = sqrt(n.tr)), x2 = seq(0,1,length = sqrt(n.tr)))
  fIn <- list(f1 = matrix(runif(n.tr*10), ncol = 10), f2 = matrix(runif(n.tr*22), ncol = 22))
  sOut <- fgp_BB3(sIn, fIn, n.tr)
  ms <- funGp(sIn = sIn, sOut = sOut)
  mf <- funGp(fIn = fIn, sOut = sOut)
  msf <- funGp(sIn = sIn, fIn = fIn, sOut = sOut)

  # generating deletion indices
  n.dl <- 2
  ind.dl <- sample(1:n.tr, n.dl)

  # updating the model
  msup <- update(ms, ind.dl = ind.dl)
  mfup <- update(mf, ind.dl = ind.dl)
  msfup <- update(msf, ind.dl = ind.dl)
  # =========================================
}
