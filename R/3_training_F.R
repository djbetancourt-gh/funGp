#' @title Optimization of hyperparameters for funGp models
#' @description Sets good values for the hyperparameters of the Gaussian process model based on optimization.
#'
#' @param fMs a list with as many elements as functional input variables. Each element of the list is a n times n matrix of differences
#' between the functional observation coordinates.
#' @param sOut a vector (or 1-column matrix) containing the values of the scalar output at the training points.
#' @param kerType a character string indicating the covariance structure to be used, to be choosen between "gauss", "matern5_2" or "matern3_2".
#' @param var.known Fill!!!!!!!!!!
#' @param ls_f.known Fill!!!!!!!!!!
#' @param n.starts Fill!!!!!!!!!!
#' @param n.presample Fill!!!!!!!!!!
#' @param nugget Fill!!!!!!!!!!
#' @return The hyperparameters.
#'
#' @keywords internal
#'
#' @author José Betancourt, François Bachoc and Thierry Klein
#' @export
setHypers_F <- function(fMs, sOut, kerType, var.known, ls_f.known, n.starts, n.presample, nugget, par.clust){
  # if the length-scale coefficients are known, skip optim and compute var analytically. Else optimize
  if (!is.null(ls_f.known)) {
    # 1. estimation of the correlation matrix
    n.tr <- length(sOut)
    R <- setR(ls_f.known, fMs, kerType) + diag(nugget, nrow = n.tr, ncol = n.tr)
    U <- chol(R)

    # 2. estimate the a priori process variance
    cat("** Computing optimal variance...\n")
    sig2 <- analyticVar_llik(U, sOut, n.tr)

    # 3. merge hyperparameters and return
    return(c(sig2, ls_f.known))

  } else {
    # 1. set hypercube for solution space
    bnds <- setBounds_F(fMs)

    # 2. set up variance function
    if (is.null(var.known)) { # the variance is computed based on the analytic formula for optimal var given ls with loglikelihood
      varfun <- analyticVar_llik
    } else { # the variance is set fixed at its known value using a closure
      g <- function(var.known) function(...) var.known
      varfun <- g(var.known)
    }

    # 3. set starting points
    cat("** Presampling...\n")
    spoints <- setSPoints_F(bnds, fMs, sOut, kerType, varfun, n.starts, n.presample, nugget)

    # 4. Perform optimization
    cat("** Optimising...\n")
    return(optimHypers_F(spoints, n.starts, bnds, fMs, sOut, kerType, varfun, nugget, par.clust))
  }
}
# -------------------------------------------------------------------------------------------------------------------------------------


#' @title Setting up of hyperparameters' boundaries for optimization
#' @description Fill this!!!!!!!!!!
#'
#' @param fMs a list with as many elements as functional input variables. Each element of the list is a n times n matrix of differences
#' between the functional observation coordinates.
#' @return A matrix with two rows and as many columns as hyperparameters the metamodel requires. The first row indicates the lower limits
#' and the second row the upper limits for the optimization.
#'
#' @keywords internal
#'
#' @author José Betancourt, François Bachoc and Thierry Klein
#' @export
setBounds_F <- function(fMs){
  # define lower and upper bounds for length-scale hypers linked to functional inputs
  mxf <- sapply(fMs, max)
  ll_f <- rep(10^-10, length(fMs))
  ul_f <- 2 * mxf

  # grouping
  wholeLims <- rbind(ll_f, ul_f)

  return(wholeLims)
}
# -------------------------------------------------------------------------------------------------------------------------------------


#' @title Setting up of starting points for optimization of hyperparameters
#' @description Fill this!!!!!!!!!!
#' @param niter Fill this!!!!!!!!!!
#' @param bnds Fill this!!!!!!!!!!
#' @param fMs a list with as many elements as functional input variables. Each element of the list is a n times n matrix of differences
#' between the functional observation coordinates.
#' @param sOut Fill this!!!!!!!!!!
#' @param kerType Fill this!!!!!!!!!!
#' @return Fill this!!!!!!!!!!
#'
#' @keywords internal
#'
#' @importFrom stats runif
#' @author José Betancourt, François Bachoc and Thierry Klein
#' @export
setSPoints_F <- function(bnds, fMs, sOut, kerType, varfun, n.starts, n.presample, nugget){
  # recover lower and upper limits
  ll <- bnds[1,]
  ul <- bnds[2,]
  n.ls <- ncol(bnds)

  # generate random uniform points to test
  allspoints <- matrix(runif(n.ls * n.presample), nrow = n.ls, ncol = n.presample)
  allspoints <- ll + allspoints * (ul - ll)

  # compute fitness of each starting point
  fitvec <- apply(allspoints, 2, negLogLik_funGp_F, fMs, sOut, kerType, varfun, nugget)

  # get the best n.starts points
  spoints <- allspoints[,order(fitvec)[1:n.starts], drop = F]

  return(spoints)
}
# -------------------------------------------------------------------------------------------------------------------------------------


#' @title Optimizing the hyperparameters
#' @description Fill this!!!!!!!!!!
#'
#' @param spoints Fill this!!!!!!!!!!
#' @param n.starts Fill this!!!!!!!!!!
#' @param bnds Fill this!!!!!!!!!!
#' @param fMs a list with as many elements as functional input variables. Each element of the list is a n times n matrix of differences
#' between the functional observation coordinates.
#' @param sOut Fill this!!!!!!!!!!
#' @param kerType Fill this!!!!!!!!!!
#' @return Fill this!!!!!!!!!!
#'
#' @keywords internal
#'
#' @importFrom foreach foreach `%dopar%`
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom doSNOW registerDoSNOW
#' @importFrom stats optim
#'
#' @author José Betancourt, François Bachoc and Thierry Klein
#' @export
optimHypers_F <- function(spoints, n.starts, bnds, fMs, sOut, kerType, varfun, nugget, par.clust){
  # if multistart is required then parallelize, else run single optimization
  if (n.starts == 1){
    optOut <- optim(par = as.numeric(spoints), fn = negLogLik_funGp_F, method = "L-BFGS-B",
                    lower = bnds[1,], upper = bnds[2,], control = list(trace = T),
                    fMs = fMs, sOut = sOut, kerType = kerType, varfun = varfun, nugget = nugget)
  } else {
    # if (!requireNamespace("foreach", quietly = TRUE)){
    # if (!getDoParRegistered()){
    if (is.null(par.clust)) {
      cat("Parallel backend register not found. Multistart optimizations done in sequence.\n\n")

      # set up progress bar
      pb <- txtProgressBar(min = 0, max = n.starts, style = 3)
      cat("\n")

      # set up progress controller
      progress <- function(n) setTxtProgressBar(pb, n)
      opts <- list(progress = progress)

      optOutList <- list()
      for (i in 1:n.starts) {
        modeval <- tryCatch(
          {
            optim(par = as.numeric(spoints[,i]), fn = negLogLik_funGp_F, method = "L-BFGS-B",
                          lower = bnds[1,], upper = bnds[2,], control = list(trace = T),
                          fMs = fMs, sOut = sOut, kerType = kerType, varfun = varfun, nugget = nugget)
          },
          error = function(e) e
        )

        if (!inherits(modeval, "error")) {
          optOutList[[i]] <- modeval
        }
        setTxtProgressBar(pb, i)
        cat("\n")
      }
      close(pb)
      # optOutList <- "%do%"(foreach(i = 1:n.starts, .errorhandling = 'remove'), {
      #   optim(par = as.numeric(spoints[,i]), fn = negLogLik_funGp_F, method = "L-BFGS-B",
      #         lower = bnds[1,], upper = bnds[2,], control = list(trace = T),
      #         fMs = fMs, sOut = sOut, kerType = kerType, varfun = varfun, nugget = nugget)})

    } else {
      cat("Parallel backend register found. Multistart optimizations done in parallel.\n")

      # set up progress bar
      pb <- txtProgressBar(min = 0, max = n.starts, style = 3)

      # set up progress controller
      progress <- function(n) setTxtProgressBar(pb, n)
      opts <- list(progress = progress)

      registerDoSNOW(par.clust)
      optOutList <- foreach(i = 1:n.starts, .errorhandling = "remove", .options.snow = opts) %dopar%
      {
          optim(par = as.numeric(spoints[,i]), fn = negLogLik_funGp_F, method = "L-BFGS-B",
                lower = bnds[1,], upper = bnds[2,], control = list(trace = T),
                fMs = fMs, sOut = sOut, kerType = kerType, varfun = varfun, nugget = nugget)
      }
      close(pb)
      # optOutList <- "%dopar%"(foreach(i = 1:n.starts, .errorhandling = 'remove'), {
      #   optim(par = as.numeric(spoints[,i]), fn = negLogLik_funGp_F, method = "L-BFGS-B",
      #         lower = bnds[1,], upper = bnds[2,], control = list(trace = T),
      #         fMs = fMs, sOut = sOut, kerType = kerType, varfun = varfun, nugget = nugget)})
    }

    # check if there are usable results
    if (length(optOutList) == 0) stop("All model optimizations crashed!")

    # extract fitness from usable results
    fitvec <- lapply(optOutList, function(sol) sol$value)

    # recover best solution
    optOut <- optOutList[[which.min(fitvec)]]
  }

  # recovering relevant information for the estimation of the process a priori variance
  thetas_f <- optOut$par
  n.tr <- length(sOut)
  R <- setR(thetas_f, fMs, kerType) + diag(nugget, nrow = n.tr, ncol = n.tr)
  U <- chol(R)

  # estimation of the variance
  sig2 <- varfun(U, sOut, n.tr)

  return(c(sig2, thetas_f))
}
# -------------------------------------------------------------------------------------------------------------------------------------


# Potential objective functions for optimization of hyperparameters
# =========================================================================================================
#' @title Setting up of starting points for optimization of hyperparameters
#' @description Fill this!!!!!!!!!!
#'
#' @param thetas Fill this!!!!!!!!!!
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
negLogLik_funGp_F <- function(thetas_f, fMs, sOut, kerType, varfun, nugget){
  # Estimation of the correlation matrix
  n.tr <- length(sOut)
  R <- setR(thetas_f, fMs, kerType) + diag(nugget, nrow = n.tr, ncol = n.tr)
  U <- chol(R)

  # Estimation of the a priori process variance
  sig2 <- varfun(U, sOut, n.tr)

  # compute loglikelihood
  llik <- -0.5 * (n.tr * log(2*pi*sig2) + 2*sum(log(diag(U))) + n.tr)

  return(-llik)
}
# =========================================================================================================
