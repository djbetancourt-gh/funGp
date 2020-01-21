#' @title Optimization of hyperparameters for funGp models
#' @description Sets good values for the hyperparameters of the Gaussian process model based on optimization.
#'
#' @param sIn a matrix of scalar input values to train the model. Each column must match an input variable and each row a training point.
#' @param fpIn a list with as many elements as functional inputs. The i-th element must be a matrix with the projection coefficients
#' for the i-th functional input.
#' @param J a list with as many elements as functional inputs. The i-th element must be the Gram matrix of the basis functions used for
#' the projection of the i-th functional input.
#' @param sMs a list with as many elements as scalar input variables. Each element of the list is a n times n matrix of differences
#' between the scalar observation coordinates.
#' @param fMs a list with as many elements as functional input variables. Each element of the list is a n times n matrix of differences
#' between the functional observation coordinates.
#' @param sOut a vector (or 1-column matrix) containing the values of the scalar output at the training points.
#' @param kerType a character string indicating the covariance structure to be used, to be choosen between "gauss", "matern5_2" or "matern3_2".
#' @param n.starts Fill!!!!!!!!!!
#' @param n.presample Fill!!!!!!!!!!
#' @return The hyperparameters.
#'
#' @keywords internal
#'
#' @author José Betancourt, François Bachoc and Thierry Klein
#' @export
setHypers_SF <- function(sIn, fpIn, J, sMs, fMs, sOut, kerType, n.starts, n.presample){
  # browser()
  # 1. set hypercube for solution space
  bnds <- setBounds_SF(sMs, fMs)

  # 2. set starting points
  cat("** Presampling...\n")
  spoints <- setSPoints_SF(bnds, sMs, fMs, sOut, kerType, n.starts, n.presample)

  # 3. Perform optimization
  cat("** Optimising...\n")
  return(optimHypers_SF(spoints, n.starts, bnds, sMs, fMs, sOut, kerType))
}
# -------------------------------------------------------------------------------------------------------------------------------------


#' @title Setting up of hyperparameters' boundaries for optimization
#' @description Fill this!!!!!!!!!!
#'
#' @param sMs a list with as many elements as scalar input variables. Each element of the list is a n times n matrix of differences
#' between the scalar observation coordinates.
#' @param fMs a list with as many elements as functional input variables. Each element of the list is a n times n matrix of differences
#' between the functional observation coordinates.
#' @return A matrix with two rows and as many columns as hyperparameters the metamodel requires. The first row indicates the lower limits
#' and the second row the upper limits for the optimization.
#'
#' @keywords internal
#'
#' @author José Betancourt, François Bachoc and Thierry Klein
#' @export
setBounds_SF <- function(sMs, fMs){
  # define lower and upper bounds for hypers

  # lower and upper bounds for length-scale hypers linked to scalar inputs
  mxs <- sapply(sMs, max)
  ll_s <- 10^-10
  ul_s <- 2 * mxs

  # lower and upper bounds for length-scale hypers linked to functional inputs
  mxf <- sapply(fMs, max)
  ll_f <- 10^-10
  ul_f <- 2 * mxf

  # grouping
  llims <- c(ll_s, ll_f)
  ulims <- c(ul_s, ul_f)
  wholeLims <- rbind(llims, ulims)

  return(wholeLims)
}
# -------------------------------------------------------------------------------------------------------------------------------------


#' @title Setting up of starting points for optimization of hyperparameters
#' @description Fill this!!!!!!!!!!
#' @param niter Fill this!!!!!!!!!!
#' @param bnds Fill this!!!!!!!!!!
#' @param sMs a list with as many elements as scalar input variables. Each element of the list is a n times n matrix of differences
#' between the scalar observation coordinates.
#' @param fMs a list with as many elements as functional input variables. Each element of the list is a n times n matrix of differences
#' between the functional observation coordinates.
#' @param sOut Fill this!!!!!!!!!!
#' @param kerType Fill this!!!!!!!!!!
#' @return Fill this!!!!!!!!!!
#'
#' @keywords internal
#'
#' @author José Betancourt, François Bachoc and Thierry Klein
#' @export
setSPoints_SF <- function(bnds, sMs, fMs, sOut, kerType, n.starts, n.presample){
  # recover lower and upper limits
  ll <- bnds[1,]
  ul <- bnds[2,]

  # generate random uniform points to test
  allspoints <- matrix(runif((ncol(bnds) * n.starts * n.presample)),
                       nrow = ncol(bnds), ncol = (n.starts * n.presample))
  allspoints <- ll + allspoints * (ul - ll)

  # compute fitness of each starting point
  fitvec <- apply(allspoints, 2, function(v) negLogLik_funGp_SF(v, sMs, fMs, sOut, kerType))

  # get the best n.starts points
  spoints <- as.matrix(allspoints[,order(fitvec, decreasing = T)[1:n.starts]])

  return(spoints)
}
# -------------------------------------------------------------------------------------------------------------------------------------


#' @title Optimizing the hyperparameters
#' @description Fill this!!!!!!!!!!
#'
#' @param spoints Fill this!!!!!!!!!!
#' @param bnds Fill this!!!!!!!!!!
#' @param sMs a list with as many elements as scalar input variables. Each element of the list is a n times n matrix of differences
#' between the scalar observation coordinates.
#' @param fMs a list with as many elements as functional input variables. Each element of the list is a n times n matrix of differences
#' between the functional observation coordinates.
#' @param sOut Fill this!!!!!!!!!!
#' @param kerType Fill this!!!!!!!!!!
#' @return Fill this!!!!!!!!!!
#'
#' @keywords internal
#'
#' @importFrom foreach foreach
#' @importFrom foreach "%do%"
#' @importFrom foreach "%dopar%"
#' @importFrom foreach getDoParRegistered
#'
#' @author José Betancourt, François Bachoc and Thierry Klein
#' @export
optimHypers_SF <- function(spoints, n.starts, bnds, sMs, fMs, sOut, kerType){
  # if multistart is required then parallelize, else run single optimization
  if (ncol(spoints) == 1){
    optOut <- optim(par = as.numeric(spoints), fn = negLogLik_funGp_SF, method = "L-BFGS-B",
                   lower = bnds[1,], upper = bnds[2,], control = list(trace = T),
                   sMs = sMs, fMs = fMs, sOut = sOut, kerType = kerType)
  } else {
    # if (!requireNamespace("foreach", quietly = TRUE)){
    if (!getDoParRegistered()){
      cat("Parallel backend register not found. Multistart optimizations done in sequence.\n\n")
      optOutList <- "%do%"(foreach(i = 1:n.starts, .errorhandling = 'remove'), {
        optim(par = as.numeric(spoints[,i]), fn = negLogLik_funGp_SF, method = "L-BFGS-B",
              lower = bnds[1,], upper = bnds[2,], control = list(trace = T),
              sMs = sMs, fMs = fMs, sOut = sOut, kerType = kerType)})
    } else {
      cat("Parallel backend register found. Multistart optimizations done in parallel.\n")
      optOutList <- "%dopar%"(foreach(i = 1:n.starts, .errorhandling = 'remove'), {
        optim(par = as.numeric(spoints[,i]), fn = negLogLik_funGp_SF, method = "L-BFGS-B",
              lower = bnds[1,], upper = bnds[2,], control = list(trace = T),
              sMs = sMs, fMs = fMs, sOut = sOut, kerType = kerType)})
    }

    # check if there are usable results
    if (length(optOutList) == 0) stop("All model optimizations crashed!")

    # extract fitness from usable results
    fitvec <- lapply(optOutList, function(sol) sol$value)

    # recover best solution
    optOut <- optOutList[[which.max(fitvec)]]
  }

  # recovering length-scale hypers linked to scalar and functional inputs
  ds <- length(sMs)
  df <- length(fMs)
  thetas_s <- optOut$par[1:ds]
  thetas_f <- optOut$par[(ds+1):(ds+df)]

  # recovering relevant information for the estimation of the process a priori variance
  n.tr <- length(sOut)
  R <- setR(thetas_s, sMs, kerType) * setR(thetas_f, fMs, kerType) # + diag(10^-8, nrow = n.tr, ncol = n.tr)
  U <- chol(R)

  # estimation of the variance
  UInvY <- backsolve(t(U), sOut, upper.tri = F)
  sig2 <- crossprod(UInvY)/n.tr

  return(c(sig2, c(thetas_s, thetas_f)))
}
# -------------------------------------------------------------------------------------------------------------------------------------


# Potential objective functions for optimization of hyperparameters
# =========================================================================================================
#' @title Setting up of starting points for optimization of hyperparameters
#' @description Fill this!!!!!!!!!!
#'
#' @param thetas Fill this!!!!!!!!!!
#' @param sMs a list with as many elements as scalar input variables. Each element of the list is a n times n matrix of differences
#' between the scalar observation coordinates.
#' @param fMs a list with as many elements as functional input variables. Each element of the list is a n times n matrix of differences
#' between the functional observation coordinates.
#' @param sOut Fill this!!!!!!!!!!
#' @param kerType Fill this!!!!!!!!!!
#' @return Fill this!!!!!!!!!!
#'
#' @keywords internal
#'
#' @author José Betancourt, François Bachoc and Thierry Klein
#' @export
negLogLik_funGp_SF <- function(thetas, sMs, fMs, sOut, kerType){
  # recovering length-scale hypers linked to scalar and functional inputs
  ds <- length(sMs)
  df <- length(fMs)
  thetas_s <- thetas[1:ds]
  thetas_f <- thetas[(ds+1):(ds+df)]

  # Estimation of the correlation matrix
  n.tr <- length(sOut)
  R <- setR(thetas_s, sMs, kerType) * setR(thetas_f, fMs, kerType) # + diag(10^-8, nrow = n.tr, ncol = n.tr)
  U <- chol(R)

  # Estimation of the a priori process variance
  UInvY <- backsolve(t(U), sOut, upper.tri = F)
  sig2 <- crossprod(UInvY)/n.tr

  # compute loglikelihood
  # DiceKgiging 2108
  llik <- -0.5 * (n.tr * log(2*pi*sig2) + 2*sum(log(diag(U))) + n.tr)

  return(-llik)
}
# =========================================================================================================
