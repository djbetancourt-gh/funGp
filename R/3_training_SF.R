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
#' @param var.known Fill!!!!!!!!!!
#' @param ls_s.known Fill!!!!!!!!!!
#' @param ls_f.known Fill!!!!!!!!!!
#' @param n.starts Fill!!!!!!!!!!
#' @param n.presample Fill!!!!!!!!!!
#' @return The hyperparameters.
#'
#' @keywords internal
#'
#' @author José Betancourt, François Bachoc and Thierry Klein
#' @export
setHypers_SF <- function(sIn, fpIn, J, sMs, fMs, sOut, kerType, var.known, ls_s.known, ls_f.known, n.starts, n.presample){
  # if all the length-scale coefficients are known, skip optim and compute var analytically. Else optimize
  if (all(!is.null(ls_s.known), !is.null(ls_f.known))) {
    # 1. estimation of the correlation matrix
    n.tr <- length(sOut)
    R <- setR(ls_s.known, sMs, kerType) * setR(ls_f.known, fMs, kerType) # + diag(10^-8, nrow = n.tr, ncol = n.tr)!!!!!!!!!!!
    U <- chol(R)

    # 2. estimate the a priori process variance
    cat("** Computing optimal variance...\n")
    sig2 <- analyticVar_llik(U, sOut, n.tr)

    # 3. merge hyperparameters and return
    return(c(sig2, ls_s.known, ls_f.known))

  } else {
    # 1. set hypercube for solution space
    if (all(is.null(ls_s.known), is.null(ls_f.known))) { # case 1: all ls coefficients are unknown
      bnds <- setBounds_SF(sMs, fMs)
    } else if (is.null(ls_f.known)) { # case 2: only the functional ls coefficients are unknown
      bnds <- setBounds_F(fMs)
    } else { # case 3: only the scalar ls coefficients are unknown
      bnds <- setBounds_S(sMs)
    }

    # 2. set up variance function
    if (is.null(var.known)) { # the variance is computed based on the analytic formula for optimal var given ls with loglikelihood
      varfun <- analyticVar_llik
    } else { # the variance is set fixed at its known value using a closure
      g <- function(var.known) function(...) var.known
      varfun <- g(var.known)
    }

    # 3. set starting points
    cat("** Presampling...\n")
    spoints <- setSPoints_SF(bnds, sMs, fMs, sOut, kerType, varfun, ls_s.known, ls_f.known, n.starts, n.presample)

    # 4. Perform optimization
    cat("** Optimising...\n")
    return(optimHypers_SF(spoints, n.starts, bnds, sMs, fMs, sOut, kerType, varfun, ls_s.known, ls_f.known))
  }
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
setSPoints_SF <- function(bnds, sMs, fMs, sOut, kerType, varfun, ls_s.known, ls_f.known, n.starts, n.presample){
  # recover lower and upper limits
  ll <- bnds[1,]
  ul <- bnds[2,]
  n.ls <- ncol(bnds)

  # generate random uniform points to test
  allspoints <- matrix(runif(n.ls * n.presample), nrow = n.ls, ncol = n.presample)
  allspoints <- ll + allspoints * (ul - ll)

  # complement the candidate points with pre-fixed ls coefficients if available
  # if (!is.null(ls_f.known)) {
  #   fullspoints <- rbind(allspoints, matrix(ls_f.known, nrow = length(ls_f.known), ncol = n.presample))
  # } else if (!is.null(ls_s.known)) {
  #   fullspoints <- rbind(matrix(ls_s.known, nrow = length(ls_s.known), ncol = n.presample), allspoints)
  # }

  # compute fitness of each starting point
  # fitvec <- apply(fullspoints, 2, negLogLik_funGp_SF, sMs, fMs, sOut, kerType, varfun, ls_s.known, ls_f.known)
  fitvec <- apply(allspoints, 2, negLogLik_funGp_SF, sMs, fMs, sOut, kerType, varfun, ls_s.known, ls_f.known)

  # get the best n.starts points
  spoints <- allspoints[,order(fitvec)[1:n.starts], drop = F]

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
optimHypers_SF <- function(spoints, n.starts, bnds, sMs, fMs, sOut, kerType, varfun, ls_s.known, ls_f.known){
  # if multistart is required then parallelize, else run single optimization
  if (n.starts == 1){
    optOut <- optim(par = as.numeric(spoints), fn = negLogLik_funGp_SF, method = "L-BFGS-B",
                   lower = bnds[1,], upper = bnds[2,], control = list(trace = T),
                   sMs = sMs, fMs = fMs, sOut = sOut, kerType = kerType,
                   varfun = varfun, ls_s.known = ls_s.known, ls_f.known = ls_f.known)
  } else {
    # if (!requireNamespace("foreach", quietly = TRUE)){
    if (!getDoParRegistered()){
      cat("Parallel backend register not found. Multistart optimizations done in sequence.\n\n")
      optOutList <- "%do%"(foreach(i = 1:n.starts, .errorhandling = 'remove'), {
        optim(par = as.numeric(spoints[,i]), fn = negLogLik_funGp_SF, method = "L-BFGS-B",
              lower = bnds[1,], upper = bnds[2,], control = list(trace = T),
              sMs = sMs, fMs = fMs, sOut = sOut, kerType = kerType,
              varfun = varfun, ls_s.known = ls_s.known, ls_f.known = ls_f.known)})
    } else {
      cat("Parallel backend register found. Multistart optimizations done in parallel.\n")
      optOutList <- "%dopar%"(foreach(i = 1:n.starts, .errorhandling = 'remove'), {
        optim(par = as.numeric(spoints[,i]), fn = negLogLik_funGp_SF, method = "L-BFGS-B",
              lower = bnds[1,], upper = bnds[2,], control = list(trace = T),
              sMs = sMs, fMs = fMs, sOut = sOut, kerType = kerType,
              varfun = varfun, ls_s.known = ls_s.known, ls_f.known = ls_f.known)})
    }

    # check if there are usable results
    if (length(optOutList) == 0) stop("All model optimizations crashed!")

    # extract fitness from usable results
    fitvec <- lapply(optOutList, function(sol) sol$value)

    # recover best solution
    optOut <- optOutList[[which.min(fitvec)]]
  }

  # recovering length-scale hypers linked to scalar and functional inputs
  ds <- length(sMs)
  df <- length(fMs)
  if (all(is.null(ls_s.known), is.null(ls_f.known))) { # case 1: all ls coefficients are unknown
    thetas_s <- optOut$par[1:ds]
    thetas_f <- optOut$par[(ds+1):(ds+df)]
  } else if (is.null(ls_f.known)) { # case 2: only the functional ls coefficients are unknown
    thetas_s <- ls_s.known
    thetas_f <- optOut$par
  } else { # case 3: only the scalar ls coefficients are unknown
    thetas_s <- optOut$par
    thetas_f <- ls_f.known
  }

  # recovering relevant information for the estimation of the process a priori variance
  n.tr <- length(sOut)
  R <- setR(thetas_s, sMs, kerType) * setR(thetas_f, fMs, kerType) # + diag(10^-8, nrow = n.tr, ncol = n.tr)
  U <- chol(R)

  # estimation of the variance
  sig2 <- varfun(U, sOut, n.tr)

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
negLogLik_funGp_SF <- function(thetas, sMs, fMs, sOut, kerType, varfun, ls_s.known, ls_f.known){
  # recovering length-scale hypers linked to scalar and functional inputs
  ds <- length(sMs)
  df <- length(fMs)
  if (all(is.null(ls_s.known), is.null(ls_f.known))) { # case 1: all ls coefficients are unknown
    thetas_s <- thetas[1:ds]
    thetas_f <- thetas[(ds+1):(ds+df)]
  } else if (is.null(ls_f.known)) { # case 2: only the functional ls coefficients are unknown
    thetas_s <- ls_s.known
    thetas_f <- thetas
  } else { # case 3: only the scalar ls coefficients are unknown
    thetas_s <- thetas
    thetas_f <- ls_f.known
  }

  # Estimation of the correlation matrix
  n.tr <- length(sOut)
  R <- setR(thetas_s, sMs, kerType) * setR(thetas_f, fMs, kerType) # + diag(10^-8, nrow = n.tr, ncol = n.tr)
  U <- chol(R)

  # Estimation of the a priori process variance
  sig2 <- varfun(U, sOut, n.tr)

  # compute loglikelihood
  # DiceKgiging 2108
  llik <- -0.5 * (n.tr * log(2*pi*sig2) + 2*sum(log(diag(U))) + n.tr)

  return(-llik)
}
# =========================================================================================================


# =========================================================================================================
analyticVar_llik <- function(U, sOut, n.tr) {
  UInvY <- backsolve(t(U), sOut, upper.tri = F)
  return(crossprod(UInvY)/n.tr)
}
# =========================================================================================================
