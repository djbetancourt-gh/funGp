#' @author José Betancourt, François Bachoc and Thierry Klein
makePreds_F <- function(fMs.tp, fMs.pp, sig2, thetas_f, kerType, L, LInvY, detail, nugget){
  # create empty prediction list
  preds <- list()

  # compute and store conditional mean and standard deviation
  K.tp <- sig2 * setR(thetas_f, fMs.tp, kerType)
  LInvK <- backsolve(L, K.tp, upper.tri = F)
  preds$mean <- t(LInvK) %*% LInvY
  preds$sd <- sqrt(pmax(sig2 - apply(LInvK, 2, crossprod), 0))

  # if user requires details, provide K.tp and K.pp
  if (detail == "full") {
    preds$K.tp <- K.tp
    R <- setR(thetas_f, fMs.pp, kerType)
    preds$K.pp <- sig2 * (R + diag(nugget, nrow = nrow(R), ncol = ncol(R)))
  }

  return(preds)
}

#' @author José Betancourt, François Bachoc and Thierry Klein
preMats_F <- function(fMs, sOut, sig2, thetas_f, kerType, nugget){
  # precompute L and LInvY matrices
  R <- setR(thetas_f, fMs, kerType)
  K.tt <- sig2 * (R + diag(nugget, nrow = nrow(R), ncol = ncol(R)))
  L <- t(chol(K.tt))
  LInvY <- backsolve(L, sOut, upper.tri = F)

  return(list(L = L, LInvY = LInvY))
}
