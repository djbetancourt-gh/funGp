# ==========================================================================================================
# Master function to manage the optimization of functional-input models
# ==========================================================================================================
setHypers_F <- function(fMs, sOut, kerType, var.known, ls_f.known, n.starts, n.presample, spoints.usr,
                        nugget, par.clust, trace, pbars, control.optim){
  # if the length-scale coefficients are known, skip optim and compute var analytically. Else optimize
  if (!is.null(ls_f.known)) {
    # 1. estimation of the correlation matrix
    n.tr <- length(sOut)
    R <- setR(ls_f.known, fMs, kerType) + diag(nugget, nrow = n.tr, ncol = n.tr)
    U <- chol(R)

    # 2. estimate the a priori process variance
    if (trace) message("** Computing optimal variance...")
    sig2 <- analyticVar_llik(U, sOut, n.tr)

    # 3. merge hyperparameters and return
    return(list(hypers = c(sig2, ls_f.known), convg = as.numeric(NA), nllik = as.numeric(NA)))

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
    if (trace) message("** Presampling...")
    spoints <- setSPoints_F(bnds, fMs, sOut, kerType, varfun, n.starts, n.presample, spoints.usr, nugget)

    # 4. Perform optimization
    if (trace) message("** Optimising hyperparameters...")
    optResult <- optimHypers_F(spoints, n.starts, bnds, fMs, sOut, kerType, varfun, nugget, par.clust, trace, pbars, control.optim)
    if (trace) message("** Hyperparameters done!")
    return(optResult)
  }
}
# ==========================================================================================================



# ==========================================================================================================
# Function to set the boundaries for hyperparameters optimization - functional inputs
# ==========================================================================================================
setBounds_F <- function(fMs){
  # define lower and upper bounds for length-scale hypers linked to functional inputs
  mxf <- sapply(fMs, max)
  ll_f <- rep(10^-10, length(fMs))
  ul_f <- 2 * mxf

  # grouping
  wholeLims <- rbind(ll_f, ul_f)

  return(wholeLims)
}
# ==========================================================================================================



# ==========================================================================================================
# Function to set the starting points for hyperparameters optimization - functional inputs
# ==========================================================================================================
#' @importFrom stats runif
setSPoints_F <- function(bnds, fMs, sOut, kerType, varfun, n.starts, n.presample, spoints.usr, nugget){
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
  spoints <- allspoints[,order(fitvec)[1:n.starts], drop = FALSE]

  # if starting points were provided by the user, use them in place of the worst randomly-generated ones
  if (!is.null(spoints.usr)) {
    n.sp.usr <- ncol(spoints.usr)
    n.toRep <- min(n.starts, n.sp.usr)
    spoints[,(n.starts - n.toRep + 1):n.starts] <- spoints.usr
  }

  return(spoints)
}
# ==========================================================================================================



# ==========================================================================================================
# Function optimize the hyperparameters of functional-input models
# ==========================================================================================================
#' @importFrom stats optim
#' @importFrom foreach setDoPar %dopar%
#' @importFrom doRNG %dorng%
optimHypers_F <- function(spoints, n.starts, bnds, fMs, sOut, kerType, varfun, nugget, par.clust, trace, pbars, control.optim){
  # if multistart is required then parallelize, else run single optimization
  if (n.starts == 1){
    optOut <- optim(par = as.numeric(spoints), fn = negLogLik_funGp_F, method = "L-BFGS-B",
                    lower = bnds[1,], upper = bnds[2,], control = control.optim,
                    fMs = fMs, sOut = sOut, kerType = kerType, varfun = varfun, nugget = nugget)
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
            optim(par = as.numeric(spoints[,i]), fn = negLogLik_funGp_F, method = "L-BFGS-B",
                  lower = bnds[1,], upper = bnds[2,], control = control.optim,
                  fMs = fMs, sOut = sOut, kerType = kerType, varfun = varfun, nugget = nugget)
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
          o <- optim(par = as.numeric(spoints[,i]), fn = negLogLik_funGp_F, method = "L-BFGS-B",
                     lower = bnds[1,], upper = bnds[2,], control = control.optim,
                     fMs = fMs, sOut = sOut, kerType = kerType, varfun = varfun, nugget = nugget)
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

  # recovering relevant information for the estimation of the process a priori variance
  thetas_f <- optOut$par
  n.tr <- length(sOut)
  R <- setR(thetas_f, fMs, kerType) + diag(nugget, nrow = n.tr, ncol = n.tr)
  U <- chol(R)

  # estimation of the variance
  sig2 <- varfun(U, sOut, n.tr)

  return(list(hypers = c(sig2, thetas_f), convg = optOut$convergence, nllik = optOut$value))
}
# ==========================================================================================================



# ==========================================================================================================
# Function to compute the negative log likelihood - functional inputs
# ==========================================================================================================
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
# ==========================================================================================================
