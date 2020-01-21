#' @author José Betancourt, François Bachoc and Thierry Klein
makeSims_F <- function(fMs.ts, fMs.ss, sig2, thetas_f, kerType, L, LInvY, nsim, nug.sim, detail){
  # create empty prediction list
  sims <- list()

  # compute and store conditional realizations
  K.ts <- sig2 * setR(thetas_f, fMs.ts, kerType)
  K.ss <- sig2 * setR(thetas_f, fMs.ss, kerType)
  LInvK <- backsolve(L, K.ts, upper.tri = F)
  ys.mean <- t(LInvK) %*% LInvY
  n.sm <- nrow(K.ss)
  ys.noise <- t(chol(K.ss - t(LInvK) %*% LInvK + diag(nug.sim, n.sm))) %*% matrix(rnorm(n.sm * nsim), n.sm, nsim)
  sims$obs <- t(matrix(ys.mean, n.sm, nsim) + ys.noise)

  if (detail == "full") {
    sims$mean <- ys.mean
    sims$sd <- sqrt(pmax(sig2 - apply(LInvK, 2, crossprod), 0))
  }

  return(sims)
}
