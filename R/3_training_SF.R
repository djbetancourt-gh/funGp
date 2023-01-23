# ==========================================================================================================
# Master function to manage the optimization of hybrid-input models
# ==========================================================================================================
setHypers_SF <- function(sMs, fMs, sOut, kerType, var.known, ls_s.known, ls_f.known, n.starts, n.presample, nugget, par.clust, trace, pbars, control.optim){
  # if all the length-scale coefficients are known, skip optim and compute var analytically. Else optimize
  if (all(!is.null(ls_s.known), !is.null(ls_f.known))) {
    # 1. estimation of the correlation matrix
    n.tr <- length(sOut)
    R <- setR(ls_s.known, sMs, kerType) * setR(ls_f.known, fMs, kerType) + diag(nugget, nrow = n.tr, ncol = n.tr)
    U <- chol(R)

    # 2. estimate the a priori process variance
    if (trace) message("** Computing optimal variance...")
    sig2 <- analyticVar_llik(U, sOut, n.tr)

    # 3. merge hyperparameters and return
    return(list(hypers = c(sig2, ls_s.known, ls_f.known), convg = as.numeric(NA), nllik = as.numeric(NA)))

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
    if (trace) message("** Presampling...")
    spoints <- setSPoints_SF(bnds, sMs, fMs, sOut, kerType, varfun, ls_s.known, ls_f.known, n.starts, n.presample, nugget)

    # 4. Perform optimization
    if (trace) message("** Optimising hyperparameters...")
    optResult <- optimHypers_SF(spoints, n.starts, bnds, sMs, fMs, sOut, kerType, varfun, ls_s.known, ls_f.known, nugget, par.clust, trace, pbars, control.optim)
    if (trace) message("** Hyperparameters done!")
    return(optResult)
  }
}
# ==========================================================================================================



# ==========================================================================================================
# Function to set the boundaries for hyperparameters optimization - hybrid inputs
# ==========================================================================================================
setBounds_SF <- function(sMs, fMs){
  # lower and upper bounds for length-scale hypers linked to scalar inputs
  mxs <- sapply(sMs, max)
  ll_s <- rep(10^-10, length(sMs))
  ul_s <- 2 * mxs

  # lower and upper bounds for length-scale hypers linked to functional inputs
  mxf <- sapply(fMs, max)
  ll_f <- rep(10^-10, length(fMs))
  ul_f <- 2 * mxf

  # grouping
  llims <- c(ll_s, ll_f)
  ulims <- c(ul_s, ul_f)
  wholeLims <- rbind(llims, ulims)

  return(wholeLims)
}
# ==========================================================================================================



# ==========================================================================================================
# Function to set the starting points for hyperparameters optimization - hybrid inputs
# ==========================================================================================================
#' @importFrom stats runif
setSPoints_SF <- function(bnds, sMs, fMs, sOut, kerType, varfun, ls_s.known, ls_f.known, n.starts, n.presample, nugget){
  # recover lower and upper limits
  ll <- bnds[1,]
  ul <- bnds[2,]
  n.ls <- ncol(bnds)

  # generate random uniform points to test
  allspoints <- matrix(runif(n.ls * n.presample), nrow = n.ls, ncol = n.presample)
  allspoints <- ll + allspoints * (ul - ll)

  # compute fitness of each starting point
  fitvec <- apply(allspoints, 2, negLogLik_funGp_SF, sMs, fMs, sOut, kerType, varfun, ls_s.known, ls_f.known, nugget)

  # get the best n.starts points
  spoints <- allspoints[,order(fitvec)[1:n.starts], drop = FALSE]

  return(spoints)
}
# ==========================================================================================================



# ==========================================================================================================
# Function optimize the hyperparameters of hybrid-input models
# ==========================================================================================================
#' @importFrom stats optim
#' @importFrom foreach setDoPar %dopar%
#' @importFrom doFuture registerDoFuture
#' @importFrom doRNG registerDoRNG %dorng%
#' @importFrom future plan cluster
#' @importFrom progressr with_progress progressor
optimHypers_SF <- function(spoints, n.starts, bnds, sMs, fMs, sOut, kerType, varfun, ls_s.known, ls_f.known, nugget, par.clust, trace, pbars, control.optim){
  # if multistart is required then parallelize, else run single optimization
  if (n.starts == 1){
    optOut <- optim(par = as.numeric(spoints), fn = negLogLik_funGp_SF, method = "L-BFGS-B",
                    lower = bnds[1,], upper = bnds[2,], control = control.optim,
                    sMs = sMs, fMs = fMs, sOut = sOut, kerType = kerType,
                    varfun = varfun, ls_s.known = ls_s.known, ls_f.known = ls_f.known, nugget = nugget)
  } else {
    if (is.null(par.clust)) {
      if (trace) message("** Parallel backend register not found. Multistart optimizations done in sequence.")

      # set up progress bar
      if (pbars) {
        pb <- txtProgressBar(min = 0, max = n.starts, style = 3)
        ## cat("\n")
      }

      optOutList <- list()
      for (i in 1:n.starts) {
        modeval <- tryCatch(
          {
            optim(par = as.numeric(spoints[,i]), fn = negLogLik_funGp_SF, method = "L-BFGS-B",
                  lower = bnds[1,], upper = bnds[2,], control = control.optim,
                  sMs = sMs, fMs = fMs, sOut = sOut, kerType = kerType,
                  varfun = varfun, ls_s.known = ls_s.known, ls_f.known = ls_f.known, nugget = nugget)
          },
          error = function(e) e
        )

        if (!inherits(modeval, "error")) {
          optOutList[[i]] <- modeval
        }
        if (pbars) {
          setTxtProgressBar(pb, i)
          ## cat("\n")
        }
      }
      if (pbars) close(pb)

    } else {
      if (trace) message("** Parallel backend register found. Multistart optimizations done in parallel.")

      # register parallel backend
      oldDoPar <- registerDoFuture()
      on.exit(with(oldDoPar, setDoPar(fun = fun, data = data, info = info)), add = TRUE)
      # registerDoFuture()
      # registerDoRNG()

      # plan(cluster, workers = par.clust)
      oplan <- plan(cluster, workers = par.clust)
      on.exit(plan(oplan), add = TRUE)

      with_progress({
        if (pbars) p <- progressor(along = 1:n.starts, auto_finish = FALSE)
        optOutList <- foreach(i = 1:n.starts, .errorhandling = "remove") %dorng% {
          o <- optim(par = as.numeric(spoints[,i]), fn = negLogLik_funGp_SF, method = "L-BFGS-B",
                     lower = bnds[1,], upper = bnds[2,], control = control.optim,
                     sMs = sMs, fMs = fMs, sOut = sOut, kerType = kerType,
                     varfun = varfun, ls_s.known = ls_s.known, ls_f.known = ls_f.known, nugget = nugget)
          ## cat("\n")
          if (pbars) p()
          return(o)
        }
      })
    }

    # check if there are usable results
    if (length(optOutList) == 0) stop("All model optimizations crashed!")

    # extract fitness from usable results
    fitvec <- lapply(optOutList, function(sol) sol$value)

    # recover best solution
    optOut <- optOutList[[which.min(fitvec)]]
  }

  if (isTRUE(control.optim$trace)) message("The function value is the negated log-likelihood")

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
  R <- setR(thetas_s, sMs, kerType) * setR(thetas_f, fMs, kerType) + diag(nugget, nrow = n.tr, ncol = n.tr)
  U <- chol(R)

  # estimation of the variance
  sig2 <- varfun(U, sOut, n.tr)

  return(list(hypers = c(sig2, c(thetas_s, thetas_f)), convg = optOut$convergence, nllik = optOut$value))
}
# ==========================================================================================================



# ==========================================================================================================
# Function to compute the negative log likelihood - hybrid inputs
# ==========================================================================================================
negLogLik_funGp_SF <- function(thetas, sMs, fMs, sOut, kerType, varfun, ls_s.known, ls_f.known, nugget){
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
  R <- setR(thetas_s, sMs, kerType) * setR(thetas_f, fMs, kerType) + diag(nugget, nrow = n.tr, ncol = n.tr)
  U <- chol(R)

  # Estimation of the a priori process variance
  sig2 <- varfun(U, sOut, n.tr)

  # compute loglikelihood
  llik <- -0.5 * (n.tr * log(2*pi*sig2) + 2*sum(log(diag(U))) + n.tr)

  return(-llik)
}
# ==========================================================================================================



# ==========================================================================================================
# Analytic function for the optimal variance parameter in log likelihood - hybrid inputs
# ==========================================================================================================
analyticVar_llik <- function(U, sOut, n.tr) {
  UInvY <- backsolve(t(U), sOut, upper.tri = FALSE)
  return(crossprod(UInvY)/n.tr)
}
# ==========================================================================================================
