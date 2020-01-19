# -------------------------------------------------------------------------------------------------------------------------------------
# funGp master function: used for construction and training of a funGP model
# -------------------------------------------------------------------------------------------------------------------------------------
#' @title Fitting of functional-input Gaussian process models
#' @description Creates a Gaussian process model based on the nature of the inputs which could be scalar, functional or hybrid.
#'              For functional inputs, the user might specify a projection method, projection dimension and distance type seeking
#'              for optimal processing time or metamodel predictability.
#'
#' @param sIn a matrix of scalar input values to train the model. Each column must match an input variable and each row a training point.
#' @param fIn a list of functional inputs values to fit the model. Each element of the list must contain a matrix.
#' @param sOut a vector (or 1-column matrix) containing the values of the scalar output at the training points.
#' @param doProj a boolean indicating whether a projection of the inputs should be done for dimension reduction.
#' @param fpDims an optional array with the projection dimension for each functional input.
#' @param kerType an optional character specifying the covariance structure to be used. To be chosen between "gauss", "matern5_2" and
#' "matern3_2". Default is "matern5_2".
#' @param disType an optional character specifying the distance function to use for the functional inputs within the covariance
#' function. To be chosen between "scalar" and "functional". Default is "functional".
#' @param n.starts Fill!!!!!!!!!!
#' @param n.presample Fill!!!!!!!!!!
#'
#' @importFrom methods new
#' @importFrom stats optim
#' @importFrom stats runif
#'
#' @examples
#'
#' # generating input data
#' n.tr <- 16
#' sIn <- as.matrix(expand.grid(x1 = seq(0,1,length = sqrt(n.tr)), x2 = seq(0,1,length = sqrt(n.tr))))
#' fIn <- list(f1 = matrix(runif(n.tr*10), ncol = 10), matrix(runif(n.tr*22), ncol = 22))
#'
#' # generating output data
#' sOut <- as.matrix(sapply(t(1:n.tr), function(i){
#'   x1 <- sIn[i,1]
#'   x2 <- sIn[i,2]
#'   f1 <- fIn[[1]][i,]
#'   f2 <- fIn[[2]][i,]
#'   as.numeric(x1 * sin(x2) + x1 * mean(f1) - x2^2 * diff(range(f2)))
#' }))
#'
#' # creating a funGp model
#' m1 <- funGp(sIn = sIn, fIn = fIn, sOut = sOut)
#'
#' # plotting the model
#'
#' @author José Betancourt
#' @export
funGp <- function(sIn = NULL, fIn = NULL, sOut, doProj = T, fpDims = NULL, kerType = "matern5_2", disType = "functional",
                  n.starts = 1, n.presample = 1) {
  # =====================================================================================================
  # Attributes checklist
  # =====================================================================================================
  # 0.  * funCall .......... call ................ functional call
  # 1.  * ds ............... scalar .............. number of scalar inputs
  # 2.  * df ............... scalar .............. number of functional inputs
  # 3.  * fDims ............ array (df) .......... dimension of each functional input
  # 4.  * sIn .............. matrix (n x ds) ..... scalar inputs
  # 5.  * fIn .............. list (df) ........... each element (n x fDims_i) contains a functional input
  # 6.  * sOut ............. matrix (n x 1) ...... scalar output
  # 7.  * n.tr ............. scalar .............. number of training points
  # 8.  * proj ............. proj ................ structures related to the projection
  # 9.    - doProj ......... boolean ............. should projection of functional inputs be done?
  # 10.   - fpDims ......... array (df) .......... projection dimension of each functional input
  # 11.   - basis .......... list (df) ........... each element (fDims_i x fpDims_i) contains the basis
  #                                                functions used for the projection of one fun. input
  # 12.   - coefs .......... list (df) ........... each element (n x fpDims_i) contains the coefficients
  #                                                used for the projection of one fun. input
  # 13. * kern ............. kernel .............. structures related to the kernel
  # 14.   - kerType ........ char ................ kernel type from {"gauss", "matern5_2", "matern3_2"}
  # 15.   - disType ........ char ................ distance type from {"scalar", "functional"}
  # 16.   - varHyp ......... scalar .............. estimated variance parameter
  # 17.   - lsHyps ......... array (ds + df) ..... estimated length-scale parameters
  # 18. * preMats .......... list (2) ............ KttInv and KttInv.sOut matrices for prediction
  # =====================================================================================================
  # browser()
  checkVal(as.list(environment()))

  # create objects of class funGpProj, funGpKern and funGp
  proj <- new("funGpProj")
  kern <- new("funGpKern")
  model <- new("funGp")

  # extract generic information from user inputs
  sOut <- as.matrix(sOut)
  n.tr <- length(sOut)

  # 3 possible cases
  # Case 1: scalar and functional
  # Case 2: functional only
  # Case 3: scalar only
  if (all(!is.null(sIn), !is.null(fIn))) { # Hybrid-input case *******************************************
    # extract information from user inputs specific to the hybrid-input case
    ds <- ncol(sIn)
    df <- length(fIn)
    fDims <- sapply(fIn, ncol)

    # Extend to other possible cases!!!!!!!!!!!!!!!!!!
    if (all(doProj, is.null(fpDims))) {
      # fpDims <- rep(3, df)
      fpDims <- c(3,2)
    }

    # project functional inputs
    basis <- fpIn <- J <- list()
    for (i in 1:df) {
      B <- (eigen(cov(fIn[[i]]))$vectors)[,1:fpDims[i]]
      fpIn[[i]] <- t(solve(t(B) %*% B) %*% t(B) %*% t(fIn[[i]]))
      J[[i]] <- t(B) %*% B
      basis[[i]] <- B
    }

    # compute scalar distance matrices
    sMs <- setScalDistance(sIn, sIn)

    # compute functional distance matrices
    fMs <- setFunDistance(fpIn, fpIn, J)

    # optimize hyperparameters
    hypers <- setHypers_SF(sIn, fpIn, J, sMs, fMs, sOut, kerType, n.starts, n.presample)
    varHyp <- hypers[1]
    lsHyps <- hypers[-1]

    # pre-commpute KttInv and KttInv.sOut matrices for prediction and add them to the model
    model@preMats <- preMats_SF(sMs, fMs, sOut, varHyp, lsHyps[1:ds], lsHyps[(ds+1):(ds+df)], kerType)

    # fill funGpProj slots specific to the hybrid-input case
    proj@doProj <- doProj
    proj@fpDims <- fpDims
    proj@basis <- basis
    proj@coefs <- fpIn

    # fill funGpModel slots specific to the hybrid-input case
    model@ds <- ds
    model@df <- df
    model@fDims <- fDims
    model@sIn <- sIn
    model@fIn <- fIn

  } else if(!is.null(fIn)) { # functional-input case ***************************************
    # extract information from user inputs specific to the functional-input case
    df <- length(fIn)
    fDims <- sapply(fIn, ncol)

    # Extend to other possible cases!!!!!!!!!!!!!!!!!!
    if (all(doProj, is.null(fpDims))) {
      # fpDims <- rep(3, df)
      fpDims <- c(3,2)
    }

    # project functional inputs
    basis <- fpIn <- J <- list()
    for (i in 1:df) { # functional-input case
      B <- (eigen(cov(fIn[[i]]))$vectors)[,1:fpDims[i]]
      fpIn[[i]] <- t(solve(t(B) %*% B) %*% t(B) %*% t(fIn[[i]]))
      J[[i]] <- t(B) %*% B
      basis[[i]] <- B
    }

    # compute functional distance matrices
    fMs <- setFunDistance(fpIn, fpIn, J)

    # optimize hyperparameters
    hypers <- setHypers_F(fpIn, J, fMs, sOut, kerType, n.starts, n.presample)
    varHyp <- hypers[1]
    lsHyps <- hypers[-1]

    # pre-commpute KttInv and KttInv.sOut matrices for prediction and add them to the model
    model@preMats <- preMats_F(fMs, sOut, varHyp, lsHyps, kerType)

    # fill funGpProj slots specific to the hybrid-input case
    proj@doProj <- doProj
    proj@fpDims <- fpDims
    proj@basis <- basis
    proj@coefs <- fpIn

    # fill funGpModel slots specific to the hybrid-input case
    model@ds <- 0
    model@df <- df
    model@fDims <- fDims
    model@fIn <- fIn

  } else if(!is.null(sIn)) { # scalar-input case *******************************************
    # extract information from user inputs specific to the scalar-input case
    ds <- ncol(sIn)

    # compute scalar distance matrices
    sMs <- setScalDistance(sIn, sIn)

    # optimize hyperparameters
    hypers <- setHypers_S(sIn, sMs, sOut, kerType, n.starts, n.presample)
    varHyp <- hypers[1]
    lsHyps <- hypers[-1]

    # pre-commpute KttInv and KttInv.sOut matrices for prediction and add them to the model
    model@preMats <- preMats_S(sMs, sOut, varHyp, lsHyps, kerType)

    # fill funGpModel slots specific to the hybrid-input case
    model@ds <- ds
    model@df <- 0
    model@sIn <- sIn

  } else { # error: no inputs were provided
    stop("User must provide either a scalar-input matrix, a functional-input list or both of them. None has been detected.")
  }

  # fill general funGpKern slots
  kern@kerType <- kerType
  kern@disType <- disType
  kern@varHyp <- varHyp
  kern@lsHyps <- lsHyps

  # fill general funGpModel slots
  model@call <- match.call()
  model@sOut <- sOut
  model@n.tr <- n.tr
  model@proj <- proj
  model@kern <- kern

  return(model)
}
# -------------------------------------------------------------------------------------------------------------------------------------



#' @title Fill!!!!!!!!!!!
#' @description Fill!!!!!!!!!!!
#'
#' @param env Fill!!!!!!!!!!!
#'
#' @author José Betancourt
checkVal <- function(env){
  # browser()
  if (all(!is.null(env$sIn), !is.null(env$fIn))) { # Hybrid-input case *******************************************

    # consistency in number of points
    if (length(unique(c(nrow(env$sIn), as.numeric(sapply(env$fIn, nrow)), length(env$sOut)))) > 1) {
      stop("Inconsistent number of points. Please check that sIn, sOut and each matrix in fIn have all the same number of rows.")
    }

    # consistency in number of points
    if (!is.null(env$fpDims)) {
      if (length(env$fpDims) != length(env$fIn)) {
        stop(paste("Inconsistent number of projection dimensions. The functional input list has", length(env$fIn), "elements, but",
                   length(env$fpDims), "projection dimensions were specified."))
      }
    }

  } else if(!is.null(env$fIn)) { # functional-input case ***************************************
    # check validity and consistency of user inputs
    if (length(unique(c(as.numeric(sapply(env$fIn, nrow)), length(env$sOut)))) > 1) {
      stop("Inconsistent number of points. Please check that sOut and each matrix in fIn have all the same number of rows.")
    }
    if (!is.null(env$fpDims)) {
      if (length(env$fpDims) != length(env$fIn)) {
        stop(paste("Inconsistent number of projection dimensions. The functional input list has", length(env$fIn), "elements, but",
                   length(env$fpDims), "projection dimensions were specified."))
      }
    }

  } else if(!is.null(env$sIn)) { # scalar-input case *******************************************
    # check validity and consistency of user inputs
    if (nrow(env$sIn) != length(env$sOut)) {
      stop("Inconsistent number of points. Please check that sIn and sOut have the same number of rows.")
    }
  }
}



# @title Just to run tests
# @description Fill!!!!
#
# @importFrom graphics layout lines par plot points
# @importFrom methods show
#
# @author José Betancourt
# @export
# torun <- function(){
#   set.seed(100)
#   prind <- sort(sample(1:length(sOut), round(length(sOut) * .2)))
#   trind <- (1:length(sOut))[-prind]
#   sOut <- as.matrix(sOut)
#
#   sIn.tr <- sIn[trind,]
#   fIn.tr <- lapply(fIn, function(M) M[trind,])
#   sOut.tr <- sOut[trind,]
#
#   sIn.pr <- sIn[prind,]
#   fIn.pr <- lapply(fIn, function(M) M[prind,])
#   sOut.pr <- sOut[prind,]
#
#   # Test 1 ms vs mf vs msf
#   # =====================================================================================================
#   ms <- funGp(sIn = sIn.tr, sOut = sOut.tr)
#   c(ms@kern@varHyp, ms@kern@lsHyps)
#   mf <- funGp(fIn = fIn.tr, sOut = sOut.tr)
#   c(mf@kern@varHyp, mf@kern@lsHyps)
#   msf <- funGp(sIn = sIn.tr, fIn = fIn.tr, sOut = sOut.tr)
#   c(msf@kern@varHyp, msf@kern@lsHyps)
#
#   show(ms)
#   show(mf)
#   show(msf)
#
#   ps <- predict(ms, sIn.pr = sIn.pr)
#   pf <- predict(mf, fIn.pr = fIn.pr)
#   psf <- predict(msf, sIn.pr = sIn.pr, fIn.pr = fIn.pr)
#
#   yps <- ps$mean
#   lls <- ps$lower95
#   uls <- ps$upper95
#
#   ypf <- pf$mean
#   llf <- pf$lower95
#   ulf <- pf$upper95
#
#   ypsf <- psf$mean
#   llsf <- psf$lower95
#   ulsf <- psf$upper95
#
#   par(mar = rep(2.5, 4))
#   layout(mat = matrix(1:3, nrow = 3))
#
#   yr <- range(c(sOut.pr, yps))
#   plot(1:length(sOut.pr), sort(sOut.pr), type = "l", col = "blue", ylim = yr)
#   lines(lls[order(sOut.pr)], col = "grey55")
#   lines(uls[order(sOut.pr)], col = "grey55")
#   lines(yps[order(sOut.pr)], col = "red")
#
#   yr <- range(c(sOut.pr, ypf))
#   plot(1:length(sOut.pr), sort(sOut.pr), type = "l", col = "blue", ylim = yr)
#   lines(llf[order(sOut.pr)], col = "grey55")
#   lines(ulf[order(sOut.pr)], col = "grey55")
#   lines(ypf[order(sOut.pr)], col = "red")
#
#   yr <- range(c(sOut.pr, ypsf))
#   plot(1:length(sOut.pr), sort(sOut.pr), type = "l", col = "blue", ylim = yr)
#   lines(llsf[order(sOut.pr)], col = "grey55")
#   lines(ulsf[order(sOut.pr)], col = "grey55")
#   lines(ypsf[order(sOut.pr)], col = "red")
#   # =====================================================================================================
#
#
#   # Test 1 against DiceKriging
#   # =====================================================================================================
#   # require(DiceKriging)
#
#   ###Test 1
#   mk <- km(design = data.frame(sIn.tr), response = sOut.tr, nugget = 10^-8, coef.trend = 0)
#   ms <- funGp(sIn = sIn.tr, sOut = sOut.tr, n.starts = 1, n.presample = 1)
#   "km:"
#   c(mk@covariance@sd2, mk@covariance@range.val)
#   "sm"
#   c(ms@kern@varHyp, ms@kern@lsHyps)
#
#   pk <- predict(mk, sIn.pr, "SK")
#   ps <- predict(ms, sIn.pr = sIn.pr)
#
#   ypk <- pk$mean
#   llk <- pk$lower95
#   ulk <- pk$upper95
#
#   yps <- ps$mean
#   lls <- ps$lower95
#   uls <- ps$upper95
#
#   yr <- range(c(sOut.pr, yps))
#   plot(sOut.pr, sOut.pr, type = "l", col = "blue", ylim = yr)
#   lines(sort(sOut.pr), lls[order(sOut.pr)], col = "grey55")
#   lines(sort(sOut.pr), uls[order(sOut.pr)], col = "grey55")
#   points(sOut.pr, yps, pch = 21, bg = "red")
#
#   lines(sort(sOut.pr), llk[order(sOut.pr)], col = "green")
#   lines(sort(sOut.pr), ulk[order(sOut.pr)], col = "green")
#   points(sOut.pr, ypk, pch = 21, bg = "green")
#   # =====================================================================================================
#
#
#   ###Test 2
#   design.fact <- expand.grid(x1=seq(0,1,length=4), x2=seq(0,1,length=4))
#   y <- apply(design.fact, 1, branin)
#   mk <- km(design = data.frame(design.fact), response = y, nugget = 10^-8)
#   ms <- funGp(sIn = as.matrix(design.fact), sOut = y, n.starts = 1, n.presample = 1)
#
#   "km:"
#   c(mk@covariance@sd2, mk@covariance@range.val)
#   "sm"
#   c(ms@kern@varHyp, ms@kern@lsHyps)
#
#   n.grid <- 50
#   x.grid <- y.grid <- seq(0,1,length=n.grid)
#   design.grid <- expand.grid(x1=x.grid, x2=y.grid)
#   response.grid <- apply(design.grid, 1, branin)
#
#   pk <- predict(mk, design.grid, "SK")
#   ps <- predict(ms, sIn.pr = as.matrix(design.grid))
#
#   ypk <- pk$mean
#   llk <- pk$lower95
#   ulk <- pk$upper95
#
#   yps <- ps$mean
#   lls <- ps$lower95
#   uls <- ps$upper95
#
#   yr <- range(c(sOut.pr, yps))
#   plot(sOut.pr, sOut.pr, type = "l", col = "blue", ylim = yr)
#   lines(sort(sOut.pr), lls[order(sOut.pr)], col = "grey55")
#   lines(sort(sOut.pr), uls[order(sOut.pr)], col = "grey55")
#   points(sOut.pr, yps, pch = 21, bg = "red")
#
#   lines(sort(sOut.pr), llk[order(sOut.pr)], col = "green")
#   lines(sort(sOut.pr), ulk[order(sOut.pr)], col = "green")
#   points(sOut.pr, ypk, pch = 21, bg = "green")
#   ###Out
#
#
#
#
#   ###In
#   pk <- predict(mk, sIn.pr, "SK")
#   ps <- predict(ms, sIn.pr = sIn.pr)
#   ###Out
#
#   ###In
#   ypk <- pk$mean
#   llk <- pk$lower95
#   ulk <- pk$upper95
#
#   yps <- ps$mean
#   lls <- ps$lower95
#   uls <- ps$upper95
#   ###Out
#
#
#
#   ###In
#   yr <- range(c(sOut.pr, yps))
#   plot(sOut.pr, sOut.pr, type = "l", col = "blue", ylim = yr)
#   lines(sort(sOut.pr), lls[order(sOut.pr)], col = "grey55")
#   lines(sort(sOut.pr), uls[order(sOut.pr)], col = "grey55")
#   points(sOut.pr, yps, pch = 21, bg = "red")
#
#   lines(sort(sOut.pr), llk[order(sOut.pr)], col = "green")
#   lines(sort(sOut.pr), ulk[order(sOut.pr)], col = "green")
#   points(sOut.pr, ypk, pch = 21, bg = "green")
#   ###Out
#
#   yr <- range(c(sOut.pr, ypf))
#   plot(sOut.pr, sOut.pr, type = "l", col = "blue", ylim = yr)
#   lines(sort(sOut.pr), llf[order(sOut.pr)], col = "grey55")
#   lines(sort(sOut.pr), ulf[order(sOut.pr)], col = "grey55")
#   points(sOut.pr, ypf, pch = 21, bg = "red")
#
#   yr <- range(c(sOut.pr, ypsf))
#   plot(sOut.pr, sOut.pr, type = "l", col = "blue", ylim = yr)
#   lines(sort(sOut.pr), llsf[order(sOut.pr)], col = "grey55")
#   lines(sort(sOut.pr), ulsf[order(sOut.pr)], col = "grey55")
#   points(sOut.pr, ypsf, pch = 21, bg = "red")
#
#   # ----
#
#   yr <- range(c(sOut.pr, yps))
#   plot(1:length(sOut.pr), sort(sOut.pr), type = "l", col = "blue", ylim = yr)
#   lines(lls[order(sOut.pr)], col = "grey55")
#   lines(uls[order(sOut.pr)], col = "grey55")
#   lines(yps[order(sOut.pr)], col = "red")
#
#   yr <- range(c(sOut.pr, ypf))
#   plot(1:length(sOut.pr), sort(sOut.pr), type = "l", col = "blue", ylim = yr)
#   lines(llf[order(sOut.pr)], col = "grey55")
#   lines(ulf[order(sOut.pr)], col = "grey55")
#   lines(ypf[order(sOut.pr)], col = "red")
#
#   yr <- range(c(sOut.pr, ypsf))
#   plot(1:length(sOut.pr), sort(sOut.pr), type = "l", col = "blue", ylim = yr)
#   lines(llsf[order(sOut.pr)], col = "grey55")
#   lines(ulsf[order(sOut.pr)], col = "grey55")
#   lines(ypsf[order(sOut.pr)], col = "red")
#
#   # require(doParallel)
#   # registerDoParallel(4)
#   #
#   # env <- foreach:::.foreachGlobals
#   # rm(list=ls(name=env), pos=env)
#   #
#   # Filter(isGeneric,ls(all.names=TRUE, env = baseenv()))
#   # 19. * opt .............. optimization ........ structures related to hyperparameters optimization
# }
