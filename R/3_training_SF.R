# ==========================================================================================================
# Master function to manage the optimization of hybrid-input models
# ==========================================================================================================
setHypers_SF <- function(sMs, fMs, sOut, kerType, var.known, ls_s.known, ls_f.known, n.starts, n.presample, nugget, par.clust, quietly){
  # if all the length-scale coefficients are known, skip optim and compute var analytically. Else optimize
  if (all(!is.null(ls_s.known), !is.null(ls_f.known))) {
    # 1. estimation of the correlation matrix
    n.tr <- length(sOut)
    R <- setR(ls_s.known, sMs, kerType) * setR(ls_f.known, fMs, kerType) + diag(nugget, nrow = n.tr, ncol = n.tr)
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
    spoints <- setSPoints_SF(bnds, sMs, fMs, sOut, kerType, varfun, ls_s.known, ls_f.known, n.starts, n.presample, nugget)

    # 4. Perform optimization
    cat("** Optimising...\n")
    # if (quietly) {
    #   hypers <- quiet(optimHypers_SF(spoints, n.starts, bnds, sMs, fMs, sOut, kerType, varfun, ls_s.known, ls_f.known, nugget, par.clust))
    # } else {
    #   hypers <- optimHypers_SF(spoints, n.starts, bnds, sMs, fMs, sOut, kerType, varfun, ls_s.known, ls_f.known, nugget, par.clust)
    # }
    hypers <- optimHypers_SF(spoints, n.starts, bnds, sMs, fMs, sOut, kerType, varfun, ls_s.known, ls_f.known, nugget, par.clust, quietly)
    cat("** Hyperparameters done!\n")
    return(hypers)
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
#' @importFrom doFuture registerDoFuture
#' @importFrom future plan cluster
#' @importFrom progressr with_progress progressor
optimHypers_SF <- function(spoints, n.starts, bnds, sMs, fMs, sOut, kerType, varfun, ls_s.known, ls_f.known, nugget, par.clust, quietly){
  # if multistart is required then parallelize, else run single optimization
  if (n.starts == 1){
    if (quietly) {
      optOut <- quiet(optim(par = as.numeric(spoints), fn = negLogLik_funGp_SF, method = "L-BFGS-B",
                            lower = bnds[1,], upper = bnds[2,], control = list(trace = TRUE),
                            sMs = sMs, fMs = fMs, sOut = sOut, kerType = kerType,
                            varfun = varfun, ls_s.known = ls_s.known, ls_f.known = ls_f.known, nugget = nugget))
    } else {
      optOut <- optim(par = as.numeric(spoints), fn = negLogLik_funGp_SF, method = "L-BFGS-B",
                      lower = bnds[1,], upper = bnds[2,], control = list(trace = TRUE),
                      sMs = sMs, fMs = fMs, sOut = sOut, kerType = kerType,
                      varfun = varfun, ls_s.known = ls_s.known, ls_f.known = ls_f.known, nugget = nugget)
    }

  } else {
    if (is.null(par.clust)) {
      cat("** Parallel backend register not found. Multistart optimizations done in sequence.\n\n")

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
            if (quietly) {
              quiet(optim(par = as.numeric(spoints[,i]), fn = negLogLik_funGp_SF, method = "L-BFGS-B",
                          lower = bnds[1,], upper = bnds[2,], control = list(trace = TRUE),
                          sMs = sMs, fMs = fMs, sOut = sOut, kerType = kerType,
                          varfun = varfun, ls_s.known = ls_s.known, ls_f.known = ls_f.known, nugget = nugget))
            } else {
              optim(par = as.numeric(spoints[,i]), fn = negLogLik_funGp_SF, method = "L-BFGS-B",
                    lower = bnds[1,], upper = bnds[2,], control = list(trace = TRUE),
                    sMs = sMs, fMs = fMs, sOut = sOut, kerType = kerType,
                    varfun = varfun, ls_s.known = ls_s.known, ls_f.known = ls_f.known, nugget = nugget)
            }
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

    } else {
      cat("** Parallel backend register found. Multistart optimizations done in parallel.\n")

      # # set up progress bar
      # pb <- txtProgressBar(min = 0, max = n.starts, style = 3)

      # register parallel backend
      registerDoFuture()
      plan(cluster, workers = par.clust)

      # # set up progress controller
      # progress <- function(n) setTxtProgressBar(pb, n)
      # opts <- list(progress = progress)

      # registerDoSNOW(par.clust)
      # optOutList <- foreach(i = 1:n.starts, .errorhandling = "remove", .options.snow = opts) %dopar%
      # {
      with_progress({
        p <- progressor(along = 1:n.starts, auto_finish = FALSE)
        optOutList <- foreach(i = 1:n.starts, .errorhandling = "remove") %dopar% {
          if (quietly) {
            o <- quiet(optim(par = as.numeric(spoints[,i]), fn = negLogLik_funGp_SF, method = "L-BFGS-B",
                             lower = bnds[1,], upper = bnds[2,], control = list(trace = TRUE),
                             sMs = sMs, fMs = fMs, sOut = sOut, kerType = kerType,
                             varfun = varfun, ls_s.known = ls_s.known, ls_f.known = ls_f.known, nugget = nugget))
          } else {
            o <- optim(par = as.numeric(spoints[,i]), fn = negLogLik_funGp_SF, method = "L-BFGS-B",
                       lower = bnds[1,], upper = bnds[2,], control = list(trace = TRUE),
                       sMs = sMs, fMs = fMs, sOut = sOut, kerType = kerType,
                       varfun = varfun, ls_s.known = ls_s.known, ls_f.known = ls_f.known, nugget = nugget)
          }
          p()
          return(o)
        }
        # close(pb)
      })
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
  R <- setR(thetas_s, sMs, kerType) * setR(thetas_f, fMs, kerType) + diag(nugget, nrow = n.tr, ncol = n.tr)
  U <- chol(R)

  # estimation of the variance
  sig2 <- varfun(U, sOut, n.tr)

  return(c(sig2, c(thetas_s, thetas_f)))
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
