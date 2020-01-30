#' @title Optimization of hyperparameters for funGp models
#' @description Sets good values for the hyperparameters of the Gaussian process model based on optimization.
#'
#' @param sIn a matrix of scalar input values to train the model. Each column must match an input variable and each row a training point.
#' @param sMs a list with as many elements as scalar input variables. Each element of the list is a n times n matrix of differences
#' between the scalar observation coordinates.
#' @param sOut a vector (or 1-column matrix) containing the values of the scalar output at the training points.
#' @param kerType a character string indicating the covariance structure to be used, to be choosen between "gauss", "matern5_2" or "matern3_2".
#' @param var.known Fill!!!!!!!!!!
#' @param ls_s.known Fill!!!!!!!!!!
#' @param n.starts Fill!!!!!!!!!!
#' @param n.presample Fill!!!!!!!!!!
#' @param nugget Fill!!!!!!!!!!
#' @return The hyperparameters.
#'
#' @keywords internal
#'
#' @author José Betancourt, François Bachoc and Thierry Klein
#' @export
setHypers_S <- function(sIn, sMs, sOut, kerType, var.known, ls_s.known, n.starts, n.presample, nugget){
  # if the length-scale coefficients are known, skip optim and compute var analytically. Else optimize
  if (!is.null(ls_s.known)) {
    # 1. estimation of the correlation matrix
    n.tr <- length(sOut)
    R <- setR(ls_s.known, sMs, kerType) + diag(nugget, nrow = n.tr, ncol = n.tr)
    U <- chol(R)

    # 2. estimate the a priori process variance
    cat("** Computing optimal variance...\n")
    sig2 <- analyticVar_llik(U, sOut, n.tr)

    # 3. merge hyperparameters and return
    return(c(sig2, ls_s.known))

  } else {
    # 1. set hypercube for solution space
    bnds <- setBounds_S(sMs)

    # 2. set up variance function
    if (is.null(var.known)) { # the variance is computed based on the analytic formula for optimal var given ls with loglikelihood
      varfun <- analyticVar_llik
    } else { # the variance is set fixed at its known value using a closure
      g <- function(var.known) function(...) var.known
      varfun <- g(var.known)
    }

    # 3. set starting points
    cat("** Presampling...\n")
    spoints <- setSPoints_S(bnds, sMs, sOut, kerType, varfun, n.starts, n.presample, nugget)

    # 4. Perform optimization
    cat("** Optimising...\n")
    return(optimHypers_S(spoints, n.starts, bnds, sMs, sOut, kerType, varfun, nugget))
  }
}
# -------------------------------------------------------------------------------------------------------------------------------------


#' @title Setting up of hyperparameters' boundaries for optimization
#' @description Fill this!!!!!!!!!!
#'
#' @param sMs a list with as many elements as scalar input variables. Each element of the list is a n times n matrix of differences
#'            between the scalar observation coordinates.
#' @return A matrix with two rows and as many columns as hyperparameters the metamodel requires. The first row indicates the lower limits
#'         and the second row the upper limits for the optimization.
#'
#' @keywords internal
#'
#' @author José Betancourt, François Bachoc and Thierry Klein
#' @export
setBounds_S <- function(sMs){
  # define lower and upper bounds for length-scale hypers linked to scalar inputs
  mxs <- sapply(sMs, max)
  ll_s <- 10^-10
  ul_s <- 2 * mxs

  # grouping
  wholeLims <- rbind(ll_s, ul_s)

  return(wholeLims)
}
# -------------------------------------------------------------------------------------------------------------------------------------


#' @title Setting up of starting points for optimization of hyperparameters
#' @description Fill this!!!!!!!!!!
#' @param niter Fill this!!!!!!!!!!
#' @param bnds Fill this!!!!!!!!!!
#' @param sMs a list with as many elements as scalar input variables. Each element of the list is a n times n matrix of differences
#' between the scalar observation coordinates.
#' @param sOut Fill this!!!!!!!!!!
#' @param kerType Fill this!!!!!!!!!!
#' @return Fill this!!!!!!!!!!
#'
#' @keywords internal
#'
#' @author José Betancourt, François Bachoc and Thierry Klein
#' @export
setSPoints_S <- function(bnds, sMs, sOut, kerType, varfun, n.starts, n.presample, nugget){
  # recover lower and upper limits
  ll <- bnds[1,]
  ul <- bnds[2,]
  n.ls <- ncol(bnds)

  # generate random uniform points to test
  allspoints <- matrix(runif(n.ls * n.presample), nrow = n.ls, ncol = n.presample)
  allspoints <- ll + allspoints * (ul - ll)

  # compute fitness of each starting point
  fitvec <- apply(allspoints, 2, negLogLik_funGp_S, sMs, sOut, kerType, varfun)

  # get the best n.starts points
  spoints <- allspoints[,order(fitvec)[1:n.starts], drop = F]

  return(spoints)
}
# -------------------------------------------------------------------------------------------------------------------------------------


globalVariables('i')
#' @title Optimizing the hyperparameters
#' @description Fill this!!!!!!!!!!
#'
#' @param spoints Fill this!!!!!!!!!!
#' @param bnds Fill this!!!!!!!!!!
#' @param sMs a list with as many elements as scalar input variables. Each element of the list is a n times n matrix of differences
#' between the scalar observation coordinates.
#' @param sOut Fill this!!!!!!!!!!
#' @param kerType Fill this!!!!!!!!!!
#' @return Fill this!!!!!!!!!!
#'
#' @keywords internal
#'
#' @author José Betancourt, François Bachoc and Thierry Klein
#' @export
optimHypers_S <- function(spoints, n.starts, bnds, sMs, sOut, kerType, varfun, nugget){
  # if multistart is required then parallelize, else run single optimization
  if (n.starts == 1){
    optOut <- optim(par = as.numeric(spoints), fn = negLogLik_funGp_S, method = "L-BFGS-B",
                    lower = bnds[1,], upper = bnds[2,], control = list(trace = T),
                    sMs = sMs, sOut = sOut, kerType = kerType, varfun = varfun, nugget = nugget)
  } else {
    # if (!requireNamespace("foreach", quietly = TRUE)){
    if (!getDoParRegistered()){
      cat("Parallel backend register not found. Multistart optimizations done in sequence.\n\n")
      optOutList <- "%do%"(foreach(i = 1:n.starts, .errorhandling = 'remove'), {
        optim(par = as.numeric(spoints[,i]), fn = negLogLik_funGp_S, method = "L-BFGS-B",
              lower = bnds[1,], upper = bnds[2,], control = list(trace = T),
              sMs = sMs, sOut = sOut, kerType = kerType, varfun = varfun, nugget = nugget)})
    } else {
      cat("Parallel backend register found. Multistart optimizations done in parallel.\n")
      optOutList <- "%dopar%"(foreach(i = 1:n.starts, .errorhandling = 'remove'), {
        optim(par = as.numeric(spoints[,i]), fn = negLogLik_funGp_S, method = "L-BFGS-B",
              lower = bnds[1,], upper = bnds[2,], control = list(trace = T),
              sMs = sMs, sOut = sOut, kerType = kerType, varfun = varfun, nugget = nugget)})
    }

    # check if there are usable results
    if (length(optOutList) == 0) stop("All model optimizations crashed!")

    # extract fitness from usable results
    fitvec <- lapply(optOutList, function(sol) sol$value)

    # recover best solution
    optOut <- optOutList[[which.min(fitvec)]]
  }

  # recovering relevant information for the estimation of the process a priori variance
  thetas_s <- optOut$par
  n.tr <- length(sOut)
  R <- setR(thetas_s, sMs, kerType) + diag(nugget, nrow = n.tr, ncol = n.tr)
  U <- chol(R)

  # estimation of the variance
  sig2 <- varfun(U, sOut, n.tr)

  return(c(sig2, thetas_s))
}
# -------------------------------------------------------------------------------------------------------------------------------------


# Potential objective functions for optimization of hyperparameters
# =========================================================================================================
#' @title Setting up of starting points for optimization of hyperparameters
#' @description Fill this!!!!!!!!!!
#'
#' @param theta_s Fill this!!!!!!!!!!
#' @param sMs a list with as many elements as scalar input variables. Each element of the list is a n times n matrix of differences
#' between the scalar observation coordinates.
#' @param sOut Fill this!!!!!!!!!!
#' @param kerType Fill this!!!!!!!!!!
#' @return Fill this!!!!!!!!!!
#'
#' @keywords internal
#'
#' @author José Betancourt, François Bachoc and Thierry Klein
#' @export
negLogLik_funGp_S <- function(thetas_s, sMs, sOut, kerType, varfun, nugget){
  # Estimation of the correlation matrix
  n.tr <- length(sOut)
  R <- setR(thetas_s, sMs, kerType) + diag(nugget, nrow = n.tr, ncol = n.tr)
  U <- chol(R)

  # Estimation of the a priori process variance
  sig2 <- varfun(U, sOut, n.tr)

  # compute loglikelihood
  llik <- -0.5 * (n.tr * log(2*pi*sig2) + 2*sum(log(diag(U))) + n.tr)

  return(-llik)
}
# =========================================================================================================
