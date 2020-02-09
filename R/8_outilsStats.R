getOut_loocv <- function(model){
  y_obs <- model@sOut
  R <- tcrossprod(model@preMats$L)/model@kern@varHyp + diag(model@nugget, nrow = model@n.tr, ncol = model@n.tr)
  Rinv <- solve(R)
  y_pre <- y_obs - diag(Rinv)^(-1) * Rinv %*% y_obs
  return(y_pre)
}
