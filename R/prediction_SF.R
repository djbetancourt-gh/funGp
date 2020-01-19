makePreds_SF <- function(sMs.tp, sMs.pp, fMs.tp, fMs.pp, sig2, thetas_s, thetas_f, kerType, L, LInvY, detail){
  # create empty prediction list
  preds <- list()

  # compute and store conditional mean and standard deviation
  K.tp <- sig2 * setR(thetas_s, sMs.tp, kerType) * setR(thetas_f, fMs.tp, kerType)
  LInvK <- backsolve(L, K.tp, upper.tri = F)
  preds$mean <- t(LInvK) %*% LInvY
  preds$sd <- sqrt(pmax(sig2 - apply(LInvK, 2, crossprod), 0))

  # if user requires details, provide K.tp and K.pp
  if (detail == "full") {
    preds$K.tp <- K.tp
    preds$K.pp <- sig2 * setR(thetas_s, sMs.pp, kerType) * setR(thetas_f, fMs.pp, kerType)
  }

  return(preds)
}

preMats_SF <- function(sMs, fMs, sOut, sig2, thetas_s, thetas_f, kerType){
  # precompute L and LInvY matrices
  K.tt <- sig2 * setR(thetas_s, sMs, kerType) * setR(thetas_f, fMs, kerType)
  L <- t(chol(K.tt))
  LInvY <- backsolve(L, sOut, upper.tri = F)

  return(list(L = L, LInvY = LInvY))
}
