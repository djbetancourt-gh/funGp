#' @importFrom stats rnorm
#' @author José Betancourt, François Bachoc and Thierry Klein
makeSims_SF <- function(sMs.ts, sMs.ss, fMs.ts, fMs.ss, sig2, thetas_s, thetas_f, kerType, L, LInvY, nsim, nugget.sm, detail, seed){
  # create empty prediction list
  sims <- list()

  # compute and store conditional realizations
  K.ts <- sig2 * setR(thetas_s, sMs.ts, kerType) * setR(thetas_f, fMs.ts, kerType)
  K.ss <- sig2 * setR(thetas_s, sMs.ss, kerType) * setR(thetas_f, fMs.ss, kerType)
  LInvK <- backsolve(L, K.ts, upper.tri = F)
  ys.mean <- t(LInvK) %*% LInvY
  n.sm <- nrow(K.ss)
  if (!is.null(seed)) set.seed(seed)
  ys.noise <- t(chol(K.ss - t(LInvK) %*% LInvK + diag(nugget.sm, n.sm))) %*% matrix(rnorm(n.sm * nsim), n.sm, nsim)
  sims$sims <- t(matrix(ys.mean, n.sm, nsim) + ys.noise)

  if (detail == "full") {
    sims$mean <- ys.mean
    sims$sd <- sqrt(pmax(sig2 - apply(LInvK, 2, crossprod), 0))
  }

  return(sims)
}
