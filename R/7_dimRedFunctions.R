# ==========================================================================================================
# Master function to request a projection
# ==========================================================================================================

#
# ----------------------------------------------------------------------------------------------------------
dimReduction <- function(fIn, df, fpDims, methvec) {
  basis <- coefs <- J <- list()
  for (i in 1:df) {
    if (fpDims[[i]] > 0) {
      switch(methvec[i],
             "B-splines" = {
               B <- proj_bsplines(fIn[[i]], fpDims[i])
             },
             "PCA" = {
               B <- proj_pca(fIn[[i]], fpDims[i])
             })
      Q <- crossprod(B)
      coefs[[i]] <- t(solve(Q, tcrossprod(t(B),fIn[[i]])))
      J[[i]] <- Q
      basis[[i]] <- B
    } else {
      basis[[i]] <- J[[i]] <- diag(ncol(fIn[[i]]))
      coefs[[i]] <- fIn[[i]]
    }
  }

  return(list(basis = basis, coefs = coefs, J = J))
}
# ----------------------------------------------------------------------------------------------------------


# ==========================================================================================================
# Basis family
# ==========================================================================================================

# B-Splines
# ----------------------------------------------------------------------------------------------------------
#' @importFrom splines splineDesign
proj_bsplines <- function(f, p){
  # ord <- 4 # order of the B-spline (degree of each polynomial - 1)
  if (p == 1) ord <- 3 else ord <- 4 # order of the B-spline (degree of each polynomial - 1)
  n.inner <- p - ord + 2 # number of inner knots
  n.outer <- ord - 1 # number of endpoint extra knots
  ll <- 1 # lower 'time' instant
  ul <- ncol(f) # upper 'time' instant
  if (n.inner < 0) {
    browser()
  }
  knots.inner <- seq(ll, ul, length.out = n.inner)
  knots.left <- rep(ll, n.outer)
  knots.right <- rep(ul, n.outer)
  knots <- c(knots.left, knots.inner, knots.right)
  return(splineDesign(knots = knots, x = ll:ul, outer.ok = T, ord = ord))
}
# ----------------------------------------------------------------------------------------------------------

# PCA
# ----------------------------------------------------------------------------------------------------------
#' @importFrom stats cov
proj_pca <- function(f, p){
  B <- (eigen(cov(f))$vectors)[,1:p]
  return(B)
}
# ----------------------------------------------------------------------------------------------------------
