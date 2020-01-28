#' @author José Betancourt, François Bachoc and Thierry Klein
makePreds_S <- function(sMs.tp, sMs.pp, sig2, thetas_s, kerType, L, LInvY, detail){
  # create empty prediction list
  preds <- list()

  # compute and store conditional mean and standard deviation
  K.tp <- sig2 * setR(thetas_s, sMs.tp, kerType)
  LInvK <- backsolve(L, K.tp, upper.tri = F)
  preds$mean <- t(LInvK) %*% LInvY
  preds$sd <- sqrt(pmax(sig2 - apply(LInvK, 2, crossprod), 0))

  # if user requires details, provide K.tp and K.pp
  if (detail == "full") {
    preds$K.tp <- K.tp
    preds$K.pp <- sig2 * setR(thetas_s, sMs.pp, kerType)
  }

  return(preds)
}

#' @author José Betancourt, François Bachoc and Thierry Klein
preMats_S <- function(sMs, sOut, sig2, thetas_s, kerType){
  # precompute L and LInvY matrices
  K.tt <- sig2 * setR(thetas_s, sMs, kerType)
  L <- t(chol(K.tt))
  LInvY <- backsolve(L, sOut, upper.tri = F)

  return(list(L = L, LInvY = LInvY))
}
