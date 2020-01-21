# ==========================================================================================================
# Class for funGp models
# ==========================================================================================================

# Constructor of the class
# ----------------------------------------------------------------------------------------------------------
#' @title Class: functional input Gaussian process model
#' @description To create a funGp object, use \link[funGp]{funGp} . See also this function for mode details.
#'
#' @slot call Object of class \code{"language"}. User call reminder.
#' @slot ds Object of class \code{"numeric"}. Number of scalar inputs.
#' @slot df Object of class \code{"numeric"}. Number of functional inputs.
#' @slot fDims Object of class \code{"numeric"}. Dimension of each functional input.
#' @slot sIn Object of class \code{"matrix"}. Scalar inputs.
#' @slot fIn Object of class \code{"list"}. Functional inputs. Each element of the list contains a functional input in the form of a matrix.
#' @slot sOut Object of class \code{"matrix"}. Scalar output.
#' @slot n.tr Object of class \code{"integer"}. Number of training points.
#' @slot proj Object of class \code{"funGpProj"}. Data structures related to the projection.
#' @slot kern Object of class \code{"funGpKern"}. Data structures related to the kernel.
#' @slot preMats Object of class \code{"list"}. L and LInvY matrices pre-computed for prediction. L is a lower diagonal matrix such that
#' \eqn{L'L} equals the training cross covariance matrix \eqn{K.tt}. On the other hand, \eqn{LInvY = L^(-1) * sOut}.
#'
#' @rdname funGp-class
#' @include funGpProj_Class.R
#' @include funGpKern_Class.R
#'
#' @author José Betancourt, François Bachoc and Thierry Klein
#' @export
setClass("funGp",
         representation(
           call = "language",          # user call reminder
           ds = "numeric",             # number of scalar inputs
           df = "numeric",             # number of functional inputs
           fDims = "numeric",          # dimension of each functional input
           sIn = "matrix",             # scalar inputs
           fIn = "list",               # each element (n x fDims_i) contains a functional input
           sOut = "matrix",            # scalar output
           n.tr = "integer",           # number of training points
           proj = "funGpProj",         # structures related to the projection
           kern = "funGpKern",         # structures related to the kernel
           preMats = "list"            # Pre-computed KttInv and KttInv.sOut matrices
         ),
         validity = function(object) {T})
# ----------------------------------------------------------------------------------------------------------



# ==========================================================================================================
# User oriented methods.
# ==========================================================================================================

# funGp master function: used for construction and training of a funGP model
# ----------------------------------------------------------------------------------------------------------
#' @title Fitting of functional-input Gaussian process models
#' @description Creates a Gaussian process model based on the nature of the inputs which could be scalar,
#' functional or hybrid. For functional inputs, the user might specify a projection method, projection
#' dimension and distance type seeking for optimal processing time or metamodel predictability.
#' @param sIn a matrix of scalar input values to train the model. Each column must match an input variable
#' and each row a training point.
#' @param fIn a list of functional inputs values to fit the model. Each element of the list must contain a
#' matrix.
#' @param sOut a vector (or 1-column matrix) containing the values of the scalar output at the training
#' points.
#' @param doProj a boolean indicating whether a projection of the inputs should be done for dimension
#' reduction.
#' @param fpDims an optional array with the projection dimension for each functional input.
#' @param kerType an optional character specifying the covariance structure to be used. To be chosen between
#' "gauss", "matern5_2" and "matern3_2". Default is "matern5_2".
#' @param disType an optional character specifying the distance function to use for the functional inputs
#' within the covariance function. To be chosen between "scalar" and "functional". Default is "functional".
#' @param n.starts Fill!!!!!!!!!!
#' @param n.presample Fill!!!!!!!!!!
#'
#' @importFrom methods new
#' @importFrom stats optim
#' @importFrom stats runif
#'
#' @examples
#' # generating input data
#' n.tr <- 25
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
#' plotLOO(m1)
#'
#' @author José Betancourt, François Bachoc and Thierry Klein
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
  checkVal_funGp(as.list(environment()))

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
    if (doProj) {
      if (is.null(fpDims)) {
        # fpDims <- rep(3, df)
        fpDims <- c(3,2)
      }

      # project functional inputs
      basis <- fpIn <- J <- list()
      for (i in 1:df) {
        if (fpDims[i] > 0) {
          B <- (eigen(cov(fIn[[i]]))$vectors)[,1:fpDims[i]]
          fpIn[[i]] <- t(solve(t(B) %*% B) %*% t(B) %*% t(fIn[[i]]))
          J[[i]] <- t(B) %*% B
        } else {
          J[[i]] <- B <- diag(ncol(fIn[[i]]))
          fpIn[[i]] <- fIn[[i]]
        }
        basis[[i]] <- B
      }
    } else {
      fpDims <- rep(0, df)
      basis <- J <- lapply(fIn, function(m) diag(ncol(m)))
      fpIn <- fIn
    }

    # compute scalar distance matrices
    sMs <- setScalDistance(sIn, sIn)

    # compute functional distance matrices
    fMs <- setFunDistance(fpIn, fpIn, J)

    # optimize hyperparameters
    hypers <- setHypers_SF(sIn, fpIn, J, sMs, fMs, sOut, kerType, n.starts, n.presample)
    varHyp <- hypers[1]
    lsHyps <- hypers[-1]

    # fill funGpKern slots specific to the functional-input case
    kern@s_lsHyps <- lsHyps[1:ds]
    kern@f_lsHyps <- lsHyps[-c(1:ds)]

    # pre-commpute KttInv and KttInv.sOut matrices for prediction and add them to the model
    model@preMats <- preMats_SF(sMs, fMs, sOut, varHyp, lsHyps[1:ds], lsHyps[(ds+1):(ds+df)], kerType)

    # fill funGpProj slots specific to the hybrid-input case
    proj@doProj <- doProj
    proj@fpDims <- fpDims
    proj@basis <- basis
    proj@coefs <- fpIn

    # fill funGp slots specific to the hybrid-input case
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

    # fill funGpKern slots specific to the functional-input case
    kern@f_lsHyps <- lsHyps

    # pre-commpute KttInv and KttInv.sOut matrices for prediction and add them to the model
    model@preMats <- preMats_F(fMs, sOut, varHyp, lsHyps, kerType)

    # fill funGpProj slots specific to the hybrid-input case
    proj@doProj <- doProj
    proj@fpDims <- fpDims
    proj@basis <- basis
    proj@coefs <- fpIn

    # fill funGp slots specific to the functional-input case
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

    # fill funGpKern slots specific to the scalar-input case
    kern@s_lsHyps <- lsHyps

    # pre-commpute KttInv and KttInv.sOut matrices for prediction and add them to the model
    model@preMats <- preMats_S(sMs, sOut, varHyp, lsHyps, kerType)

    # fill funGp slots specific to the scalar-input case
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
#' @author José Betancourt, François Bachoc and Thierry Klein
checkVal_funGp <- function(env){
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
# ----------------------------------------------------------------------------------------------------------


# Method to print a funGp model
# ----------------------------------------------------------------------------------------------------------
#' @name show
#' @description This is my description
#' @rdname show-methods
#' @importFrom methods show
#' @param object An object to show.
#'
#' @author José Betancourt, François Bachoc and Thierry Klein
if(!isGeneric("show")) {setGeneric(name = "show", def = function(object) standardGeneric("show"))}

#' @title Fill!!!!!!!!!!!
#' @name show
#' @rdname show-methods
#' @aliases show,funGp-method
setMethod("show", "funGp", function(object) show.funGp(object))

show.funGp <- function(object) {
  mainTxt <- "Gaussian Process Model"
  callTxt <- paste("* Call: ", as.expression(object@call), sep = "")
  cat(paste("\n", mainTxt, paste(rep("_", min(30, (nchar(callTxt) - nchar(mainTxt) - 1))), collapse = "")))

  cat(paste("\n\n", callTxt, "\n\n", sep = ""))

  cat(paste("* Scalar inputs: ", object@ds, "\n", sep = ""))
  cat(paste("* Functional inputs: ", object@df, "\n", sep = ""))
  if (object@df > 0) {
    cat("  -> Dimension:\n")
    for (i in 1:object@df) {
      cat(paste("\t F", i, ": ", object@fDims[i], "\n", sep = ""))
    }
  }
  cat(paste("* Training points: ", object@n.tr, "\n\n", sep = ""))

  cat(paste("* Kernel type: ", object@kern@kerType, "\n", sep = ""))
  cat(paste("* Distance type: ", object@kern@disType, "\n\n", sep = ""))

  if (object@df > 0) {
    cat(paste("* Do projection: ", object@proj@doProj, "\n", sep = ""))
    if (object@proj@doProj) {
      cat("  -> Proj. dimension:\n")
      for (i in 1:object@df) {
        if (object@proj@fpDims[i] > 0) {
          cat(paste("\t F", i, ": ", object@proj@fpDims[i], "\n", sep = ""))
        } else {
          cat(paste("\t F", i, ": not required\n", sep = ""))
        }
      }
    }
  }

  cat("\n* Hyperparameters:\n")
  cat(paste("  -> variance: ", format(object@kern@varHyp, digits = 3, nsmall = 4), "\n", sep = ""))
  cat("  -> length-scale:\n")
  if (object@ds > 0) {
    for (i in 1:object@ds) {
      cat(paste("\t ls(X", i, "): ", format(object@kern@s_lsHyps[i], digits = 3, nsmall = 4), "\n", sep = ""))
    }
  }
  if (object@df > 0) {
    for (i in 1:object@df) {
      cat(paste("\t ls(F", i, "): ", format(object@kern@f_lsHyps[i], digits = 3, nsmall = 4), "\n", sep = ""))
    }
  }
  cat(paste(rep("_", max(30, (nchar(callTxt)))), collapse = ""))
}
# ----------------------------------------------------------------------------------------------------------

# Method to make predictions with afunGp model
# ----------------------------------------------------------------------------------------------------------
#' @name predict
#' @rdname predict-methods
#' @importFrom stats predict qnorm
#' @param object An object to predict from.
#' @param sIn.pr fill!!
#' @param fIn.pr fill!!
#' @param detail fill!!
#' @param ... Further arguments for methods.
#'
#' @examples
#' # generating input data for training
#' n.tr <- 25
#' sIn.tr <- expand.grid(x1 = seq(0,1,length = sqrt(n.tr)), x2 = seq(0,1,length = sqrt(n.tr)))
#' sIn.tr <- as.matrix(sIn.tr)
#' fIn.tr <- list(f1 = matrix(runif(n.tr*10), ncol = 10), matrix(runif(n.tr*22), ncol = 22))
#'
#' # generating output data for training
#' sOut.tr <- as.matrix(sapply(t(1:n.tr), function(i){
#'   x1 <- sIn.tr[i,1]
#'   x2 <- sIn.tr[i,2]
#'   f1 <- fIn.tr[[1]][i,]
#'   f2 <- fIn.tr[[2]][i,]
#'   as.numeric(x1 * sin(x2) + x1 * mean(f1) - x2^2 * diff(range(f2)))
#' }))
#'
#' # creating a funGp model
#' m1 <- funGp(sIn = sIn.tr, fIn = fIn.tr, sOut = sOut.tr)
#'
#' # generating input data for prediction
#' n.pr <- 100
#' sIn.pr <- expand.grid(x1 = seq(0,1,length = sqrt(n.pr)), x2 = seq(0,1,length = sqrt(n.pr)))
#' sIn.pr <- as.matrix(sIn.pr)
#' fIn.pr <- list(f1 = matrix(runif(n.pr*10), ncol = 10), matrix(runif(n.pr*22), ncol = 22))
#'
#' # making predictions
#' m1.preds <- predict(m1, sIn.pr = sIn.pr, fIn.pr = fIn.pr)
#'
#' # plotting predictions
#' plotPreds(m1, preds = m1.preds)
#'
#' # also possible to compare against true output values
#' sOut.pr <- as.matrix(sapply(t(1:n.pr), function(i){
#'   x1 <- sIn.pr[i,1]
#'   x2 <- sIn.pr[i,2]
#'   f1 <- fIn.pr[[1]][i,]
#'   f2 <- fIn.pr[[2]][i,]
#'   as.numeric(x1 * sin(x2) + x1 * mean(f1) - x2^2 * diff(range(f2)))
#' }))
#' plotPreds(m1, m1.preds, sOut.pr)
#'
#' @author José Betancourt, François Bachoc and Thierry Klein
#' @export predict
setGeneric(name = "predict", def = function(object, ...) standardGeneric("predict"))

#' @title Prediction method for the funGp Class
#' @name predict
#' @rdname predict-methods
#' @aliases predict,funGp-method
setMethod("predict", "funGp",
          function(object, sIn.pr = NULL, fIn.pr = NULL, detail = "light", ...){
            predict.funGp(object = object, sIn.pr = sIn.pr, fIn.pr = fIn.pr, detail = detail, ...)
          })

predict.funGp <- function(object, sIn.pr, fIn.pr, detail = "light", ...) {
  # =====================================================================================================
  # Prediction output checklist
  # =====================================================================================================
  # 1.  * mean ............... array (n.pr) .............. predicted mean
  # 2.  * sd ................. array (n.pr) .............. predicted standard deviation
  # 3.  * lower95 ............ array (n.pr) .............. lower bounds of 95% confidence intervals
  # 4.  * upper95 ............ array (n.pr) .............. upper bounds of 95% confidence intervals
  # 5.  * K.pp ............... matrix(n.pr x n.pr) ....... conditional covariance matrix
  # 6.  * K.tp ............... matrix(n.tr x n.pr) ....... training vs prediction cross covariance matrix
  # =====================================================================================================

  if (all(object@ds > 0, object@df > 0)) { # Hybrid-input case *******************************************
    print("I'm hybrid!")

    # project functional inputs
    fpIn.pr <- J <- list()
    for (i in 1:object@df) {
      B <- object@proj@basis[[i]]
      fpIn.pr[[i]] <- t(solve(t(B) %*% B) %*% t(B) %*% t(fIn.pr[[i]]))
      J[[i]] <- t(B) %*% B
    }

    # compute scalar distance matrices
    sMs.tp <- setScalDistance(object@sIn, sIn.pr)
    sMs.pp <- setScalDistance(sIn.pr, sIn.pr)

    # compute functional distance matrices
    fMs.tp <- setFunDistance(object@proj@coefs, fpIn.pr, J)
    fMs.pp <- setFunDistance(fpIn.pr, fpIn.pr, J)

    # make predictions based on the Gaussian Conditioning Theorem
    preds <- makePreds_SF(sMs.tp, sMs.pp, fMs.tp, fMs.pp,
                          object@kern@varHyp, object@kern@s_lsHyps, object@kern@f_lsHyps,
                          object@kern@kerType, object@preMats$L, object@preMats$LInvY, detail)

  } else if (object@df > 0) { # functional-input case *******************************************
    print("I'm functional!")

    # project functional inputs
    fpIn.pr <- J <- list()
    for (i in 1:object@df) {
      B <- object@proj@basis[[i]]
      fpIn.pr[[i]] <- t(solve(t(B) %*% B) %*% t(B) %*% t(fIn.pr[[i]]))
      J[[i]] <- t(B) %*% B
    }

    # compute functional distance matrices
    fMs.tp <- setFunDistance(object@proj@coefs, fpIn.pr, J)
    fMs.pp <- setFunDistance(fpIn.pr, fpIn.pr, J)

    # make predictions based on the Gaussian Conditioning Theorem
    preds <- makePreds_F(fMs.tp, fMs.pp, object@kern@varHyp, object@kern@f_lsHyps, object@kern@kerType,
                         object@preMats$L, object@preMats$LInvY, detail)

  } else { # scalar-input case *******************************************
    print("I'm scalar!")

    # compute scalar distance matrices
    sMs.tp <- setScalDistance(object@sIn, sIn.pr)
    sMs.pp <- setScalDistance(sIn.pr, sIn.pr)

    # make predictions based on the Gaussian Conditioning Theorem
    preds <- makePreds_S(sMs.tp, sMs.pp, object@kern@varHyp, object@kern@s_lsHyps, object@kern@kerType,
                         object@preMats$L, object@preMats$LInvY, detail)
  }

  # compute confidence intervals
  preds$lower95 <- preds$mean - qnorm(0.975) * preds$sd
  preds$upper95 <- preds$mean + qnorm(0.975) * preds$sd

  return(preds)
}
# ----------------------------------------------------------------------------------------------------------


# Method to make simulations from afunGp model
# ----------------------------------------------------------------------------------------------------------
#' @name simulate
#' @rdname simulate-methods
#' @importFrom stats simulate rnorm
#' @param object An object to simulate from.
#' @param nsim fill!!!
#' @param seed fill!!!
#' @param sIn.sm fill!!
#' @param fIn.sm fill!!
#' @param nug.sim fill!!!
#' @param detail fill!!!
#' @param ... Further arguments for methods.
#'
#' @examples
#' # generating input data for training
#' n.tr <- 25
#' sIn.tr <- expand.grid(x1 = seq(0,1,length = sqrt(n.tr)), x2 = seq(0,1,length = sqrt(n.tr)))
#' sIn.tr <- as.matrix(sIn.tr)
#' fIn.tr <- list(f1 = matrix(runif(n.tr*10), ncol = 10), matrix(runif(n.tr*22), ncol = 22))
#'
#' # generating output data for training
#' sOut.tr <- as.matrix(sapply(t(1:n.tr), function(i){
#'   x1 <- sIn.tr[i,1]
#'   x2 <- sIn.tr[i,2]
#'   f1 <- fIn.tr[[1]][i,]
#'   f2 <- fIn.tr[[2]][i,]
#'   as.numeric(x1 * sin(x2) + x1 * mean(f1) - x2^2 * diff(range(f2)))
#' }))
#'
#' # creating a funGp model
#' m1 <- funGp(sIn = sIn.tr, fIn = fIn.tr, sOut = sOut.tr)
#'
#' # generating input data for simulation
#' n.sm <- 100
#' sIn.sm <- expand.grid(x1 = seq(0,1,length = sqrt(n.sm)), x2 = seq(0,1,length = sqrt(n.sm)))
#' sIn.sm <- as.matrix(sIn.sm)
#' fIn.sm <- list(f1 = matrix(runif(n.sm*10), ncol = 10), matrix(runif(n.sm*22), ncol = 22))
#'
#' # making simulations
#' m1.sims <- simulate(m1, nsim = 10, sIn.sm = sIn.sm, fIn.sm = fIn.sm, detail = "full")
#' matplot(t(m1.sims$obs), type = 'l', col = "grey")
#' lines(m1.sims$mean, col = "red")
#' lines(m1.sims$lower95, col = "blue")
#' lines(m1.sims$upper95, col = "blue")
#'
#' @author José Betancourt, François Bachoc and Thierry Klein
#' @export simulate
setGeneric(name = "simulate", def = function(object, nsim = 1, seed = NULL, ...) standardGeneric("simulate"))

#' @title Simulation method for the funGp Class
#' @name simulate
#' @rdname simulate-methods
#' @aliases simulate,funGp-method
setMethod("simulate", "funGp",
          function(object, nsim = 1, seed = NULL, sIn.sm = NULL, fIn.sm = NULL, nug.sim = 0, detail = "light", ...) {
            simulate.funGp(object = object, nsim = nsim, seed = seed, sIn.sm = sIn.sm, fIn.sm = fIn.sm,
                           nug.sim = nug.sim, detail = detail, ...)
          })

simulate.funGp <- function(object, nsim, seed, sIn.sm, fIn.sm, nug.sim, detail, ...) {
  checkVal_simulate(as.list(environment()))

  # check which type of model it is
  if (all(object@ds > 0, object@df > 0)) { # Hybrid-input case *******************************************
    print("I'm hybrid!")

    # project functional inputs
    fpIn.sm <- J <- list()
    for (i in 1:object@df) {
      B <- object@proj@basis[[i]]
      fpIn.sm[[i]] <- t(solve(t(B) %*% B) %*% t(B) %*% t(fIn.sm[[i]]))
      J[[i]] <- t(B) %*% B
    }

    # compute scalar distance matrices
    sMs.ts <- setScalDistance(object@sIn, sIn.sm)
    sMs.ss <- setScalDistance(sIn.sm, sIn.sm)

    # compute functional distance matrices
    fMs.ts <- setFunDistance(object@proj@coefs, fpIn.sm, J)
    fMs.ss <- setFunDistance(fpIn.sm, fpIn.sm, J)

    # make simulations based on the Gaussian Conditioning Theorem
    sims <- makeSims_SF(sMs.ts, sMs.ss, fMs.ts, fMs.ss,
                        object@kern@varHyp, object@kern@s_lsHyps, object@kern@f_lsHyps,
                        object@kern@kerType, object@preMats$L, object@preMats$LInvY, nsim, nug.sim, detail)

  } else if (object@df > 0) { # functional-input case *******************************************
    print("I'm functional!")

    # compute functional distance matrices
    fMs.ts <- setFunDistance(object@proj@coefs, fpIn.sm, J)
    fMs.ss <- setFunDistance(fpIn.sm, fpIn.sm, J)

    # make simulations based on the Gaussian Conditioning Theorem
    sims <- makeSims_F(fMs.ts, fMs.ss, object@kern@varHyp, object@kern@f_lsHyps, object@kern@kerType,
                       object@preMats$L, object@preMats$LInvY, nsim, nug.sim, detail)

  } else { # scalar-input case *******************************************
    print("I'm scalar!")

    # compute scalar distance matrices
    sMs.ts <- setScalDistance(object@sIn, sIn.sm)
    sMs.ss <- setScalDistance(sIn.sm, sIn.sm)

    # make simulations based on the Gaussian Conditioning Theorem
    sims <- makeSims_S(sMs.ts, sMs.ss, object@kern@varHyp, object@kern@s_lsHyps, object@kern@kerType,
                       object@preMats$L, object@preMats$LInvY, nsim, nug.sim, detail)
  }

  # if detail == 'full', confidence intervals at simulation points are provided,
  # else the sims list is dropped to a matrix with the observations only
  if (detail == "full") {
    # compute confidence intervals
    sims$lower95 <- sims$mean - qnorm(0.975) * sims$sd
    sims$upper95 <- sims$mean + qnorm(0.975) * sims$sd
  } else {
    sims <- sims$obs
  }

  return(sims)
}
# ----------------------------------------------------------------------------------------------------------



# ==========================================================================================================
# Plotters
# ==========================================================================================================

# Method to plot a funGp model
# ----------------------------------------------------------------------------------------------------------
#' @name plotLOO
#' @description This is my description
#' @rdname plotLOO-methods
#' @importFrom graphics lines plot
#' @param object An object to predict from.
#'
#' @author José Betancourt, François Bachoc and Thierry Klein
#' @exportMethod plotLOO
if(!isGeneric("plotLOO")) {setGeneric("plotLOO", function(object) standardGeneric("plotLOO"))}

#' @title Prediction Method for the apk Class
#' @name plotLOO
#' @rdname plotLOO-methods
#' @aliases plotLOO,funGp-method
setMethod("plotLOO", "funGp", function(object) plotLOO.funGp(object))

plotLOO.funGp <- function(object) {
  y_obs <- object@sOut
  R <- tcrossprod(object@preMats$L)/object@kern@varHyp
  Rinv <- solve(R)
  y_pre <- y_obs - diag(Rinv)^(-1) * Rinv %*% y_obs
  yr <- range(c(y_obs, y_pre))
  plot(y_obs, y_pre, xlim = yr, ylim = yr, pch = 21, col = "red", bg = "red", xlab = "Observed", ylab = "Predicted")
  lines(y_obs, y_obs, col = "blue")
}
# ----------------------------------------------------------------------------------------------------------


# Method to plot predictions of a funGp model
# ----------------------------------------------------------------------------------------------------------
#' @name plotPreds
#' @description This is my description
#' @rdname plotPreds-methods
#' @importFrom graphics lines plot polygon layout legend par
#' @param object An object to predict from.
#' @param ... Further arguments for methods.
#'
#' @author José Betancourt, François Bachoc and Thierry Klein
#' @exportMethod plotPreds
if(!isGeneric("plotPreds")) {setGeneric("plotPreds", function(object, ...) standardGeneric("plotPreds"))}

#' @title Prediction Method for the apk Class
#' @name plotPreds
#' @rdname plotPreds-methods
#' @aliases plotPreds,funGp-method
#' @param preds something
#' @param sOut.pr also
setMethod("plotPreds", "funGp",
          function(object, preds, sOut.pr = NULL, ...) {
            plotPreds.funGp(preds = preds, sOut.pr = sOut.pr)
          })

plotPreds.funGp <- function(preds, sOut.pr) {
  if (!is.null(sOut.pr)) {
    layout(matrix(1:2, nrow = 2))
    par(mar = c(4.2,4.1,1,2.1))
    # calibration plot
    y_obs <- sOut.pr
    y_pre <- preds$mean
    yr <- range(c(y_obs, y_pre))
    plot(y_obs, y_pre, xlim = yr, ylim = yr, pch = 21, col = "red", bg = "red", xlab = "Observed", ylab = "Predicted")
    lines(y_obs, y_obs, col = "blue")

    # prediction, limits and true curve
    y <- sort(preds$mean)
    n.pr <- length(y)
    ll <- (preds$lower95)[order(preds$mean)]
    ul <- (preds$upper95)[order(preds$mean)]

    plot(1, type = "n", xlab = "Index", ylab = "Predicted", xlim = c(1, n.pr), ylim = range(y))
    x <- 1:n.pr
    polygon(c(x, rev(x)), c(ul, rev(ll)), col = "grey85", border = NA)
    lines(sOut.pr[order(preds$mean)], col = "black")
    lines(y, col = "red")
    lines(ll, col = "blue")
    lines(ul, col = "blue")
  } else {
    # prediction, limits and true curve
    y <- sort(preds$mean)
    n.pr <- length(y)
    ll <- (preds$lower95)[order(preds$mean)]
    ul <- (preds$upper95)[order(preds$mean)]

    plot(1, type = "n", xlab = "Index", ylab = "Predicted", xlim = c(1, n.pr), ylim = range(y))
    x <- 1:n.pr
    polygon(c(x, rev(x)), c(ul, rev(ll)), col = "grey85", border = NA)
    lines(y, col = "red")
    lines(ll, col = "blue")
    lines(ul, col = "blue")
  }
  legend("topleft", legend = c("True", "Predicted", "95% CIs"), col = c("black", "red", "blue"), lty = 1, cex = 0.8)

  # reset plotting default setup
  par(mar = c(5.1, 4.1, 4.1, 2.1), mfrow = c(1,1))
}
# ----------------------------------------------------------------------------------------------------------


# Method to plot simulations of a funGp model
# ----------------------------------------------------------------------------------------------------------
#' @name plotSims
#' @description This is my description
#' @rdname plotSims-methods
#' @param object An object to predict from.
#' @param ... Further arguments for methods.
#'
#' @importFrom graphics matplot
#' @author José Betancourt, François Bachoc and Thierry Klein
#' @exportMethod plotSims
if(!isGeneric("plotSims")) {setGeneric("plotSims", function(object, ...) standardGeneric("plotSims"))}

#' @title Fill Method for the apk Class
#' @name plotSims
#' @rdname plotSims-methods
#' @aliases plotSims,funGp-method
#' @param sims something
#' @param detail fill!!!!
setMethod("plotSims", "funGp",
          function(object, sims = sims, detail = "full", ...) {
            plotSims.funGp(sims = sims, detail = detail)
          })

plotSims.funGp <- function(sims, detail) {
  if (!is.list(sims)) {
    detail <- "light"
  }

  if (detail == "light") {
    matplot(t(sims$obs), type = "l", lty = 1, col = "palegreen4", xlab = "Sim. index", ylab = "Output")
  } else {
    matplot(t(sims$obs), type = "l", lty = 1, col = "grey", xlab = "Sim. index", ylab = "Output")
    lines(sims$mean, col = "red")
    lines(sims$lower95, col = "blue")
    lines(sims$upper95, col = "blue")
    legend("topleft", legend = c("Sims", "Mean", "95% CIs"), col = c("grey50", "red", "blue"), lty = 1, cex = 0.8)
  }
}
# ----------------------------------------------------------------------------------------------------------



# ==========================================================================================================
# Getters
# ==========================================================================================================

# Method to get the list of hyperparameters of a funGp model
# ----------------------------------------------------------------------------------------------------------
#' @name getCoef
#' @description This is my description
#' @rdname getCoef-methods
#' @param object An object to predict from.
#'
#' @author José Betancourt, François Bachoc and Thierry Klein
#' @exportMethod getCoef
if(!isGeneric("getCoef")) {setGeneric(name = "getCoef", def = function(object) standardGeneric("getCoef"))}

#' @title Prediction Method for the apk Class
#' @name getCoef
#' @rdname getCoef-methods
#' @aliases getCoef,funGp-method
setMethod("getCoef", "funGp", function(object) getCoef.funGp(object))

getCoef.funGp <- function(object) {
  coefs <- object@kern@varHyp
  names_ls_s <- c()
  if (object@ds > 0) {
    names_ls_s <- paste("ls(X", 1:object@ds, ")", sep = "")
    coefs <- c(coefs, object@kern@s_lsHyps)
  }
  names_ls_f <- c()
  if (object@df > 0) {
    names_ls_f <- paste("ls(F", 1:object@df, ")", sep = "")
    coefs <- c(coefs, object@kern@f_lsHyps)
  }
  names(coefs) <- c("var", names_ls_s, names_ls_f)
  return(coefs)
}


# Method to get the list of basis functions used for the projection of the functional inputs
# ----------------------------------------------------------------------------------------------------------
#' @name getBasis
#' @description This is my description
#' @rdname getBasis-methods
#' @param object An object to predict from.
#'
#' @author José Betancourt, François Bachoc and Thierry Klein
#' @exportMethod getBasis
if(!isGeneric("getBasis")) {setGeneric(name = "getBasis", def = function(object) standardGeneric("getBasis"))}

#' @title Prediction Method for the apk Class
#' @name getBasis
#' @rdname getBasis-methods
#' @aliases getBasis,funGp-method
setMethod("getBasis", "funGp", function(object) getBasis.funGp(object))

getBasis.funGp <- function(object) {
  if (object@df > 0) {
    return(object@proj@basis)
  } else {
    cat("The provided funGp model does not have functional inputs. Basis functions are not defined for it.")
  }
}


# Method to get the list of reduced functional inputs used in the model (coefficents of the projection)
# ----------------------------------------------------------------------------------------------------------
#' @name getRedfIn
#' @description This is my description
#' @rdname getRedfIn-methods
#' @param object An object to predict from.
#'
#' @author José Betancourt, François Bachoc and Thierry Klein
#' @exportMethod getRedfIn
if(!isGeneric("getRedfIn")) {setGeneric(name = "getRedfIn", def = function(object) standardGeneric("getRedfIn"))}

#' @title Prediction Method for the apk Class
#' @name getRedfIn
#' @rdname getRedfIn-methods
#' @aliases getRedfIn,funGp-method
setMethod("getRedfIn", "funGp", function(object) getRedfIn.funGp(object))

getRedfIn.funGp <- function(object) {
  if (object@df > 0) {
    rfIn <- object@proj@coefs
    names(rfIn) <- paste("alpha(F", 1:object@df, ")", sep = "")
    return(rfIn)
  } else {
    cat("The provided funGp model does not have functional inputs. Reduced functional inputs are not defined for it.")
  }
}


# Method to get the list of projected functional inputs (not used directly in the model)
# ----------------------------------------------------------------------------------------------------------
#' @name getProjfIn
#' @description This is my description
#' @rdname getProjfIn-methods
#' @param object An object to predict from.
#'
#' @author José Betancourt, François Bachoc and Thierry Klein
#' @exportMethod getProjfIn
if(!isGeneric("getProjfIn")) {setGeneric(name = "getProjfIn", def = function(object) standardGeneric("getProjfIn"))}

#' @title Prediction Method for the apk Class
#' @name getProjfIn
#' @rdname getProjfIn-methods
#' @aliases getProjfIn,funGp-method
setMethod("getProjfIn", "funGp", function(object) getProjfIn.funGp(object))

getProjfIn.funGp <- function(object) {
  if (object@df > 0) {
    fpIn <- mapply(function(m1, m2) m1 %*% t(m2), m1 = object@proj@coefs, m2 = object@proj@basis)
    names(fpIn) <- paste("Fp", 1:object@df, sep = "")
    return(fpIn)
  } else {
    cat("The provided funGp model does not have functional inputs. Projected functional inputs are not defined for it.")
  }
}


# Method to get the list of Gram matrices computed from the basis functions
# ----------------------------------------------------------------------------------------------------------
#' @name getProjgram
#' @description This is my description
#' @rdname getProjgram-methods
#' @param object An object to predict from.
#'
#' @author José Betancourt, François Bachoc and Thierry Klein
#' @exportMethod getProjgram
if(!isGeneric("getProjgram")) {setGeneric("getProjgram", function(object) standardGeneric("getProjgram"))}

#' @title Prediction Method for the apk Class
#' @name getProjgram
#' @rdname getProjgram-methods
#' @aliases getProjgram,funGp-method
setMethod("getProjgram", "funGp", function(object) getProjgram.funGp(object))

getProjgram.funGp <- function(object) {
  if (object@df > 0) {
    gram <- list()
    for (i in 1:object@df) {
      gram[[i]] <- t(object@proj@basis[[i]]) %*% object@proj@basis[[i]]
    }
    names(gram) <- paste("G(F", 1:object@df, ")", sep = "")
    return(gram)
  } else {
    cat("The provided funGp model does not have functional inputs. The projection Gram matrix is not defined for it.")
  }
}
# ----------------------------------------------------------------------------------------------------------


# Method to get the training covariance matrix
# ----------------------------------------------------------------------------------------------------------
#' @name getTrainCov
#' @description This is my description
#' @rdname getTrainCov-methods
#' @param object An object to predict from.
#'
#' @author José Betancourt, François Bachoc and Thierry Klein
#' @exportMethod getTrainCov
if(!isGeneric("getTrainCov")) {setGeneric("getTrainCov", function(object) standardGeneric("getTrainCov"))}

#' @title Prediction Method for the apk Class
#' @name getTrainCov
#' @rdname getTrainCov-methods
#' @aliases getTrainCov,funGp-method
setMethod("getTrainCov", "funGp", function(object) getTrainCov.funGp(object))

getTrainCov.funGp <- function(object) {
  return(tcrossprod(object@preMats$L))
}
# ----------------------------------------------------------------------------------------------------------



# ==========================================================================================================
# Validators
# ==========================================================================================================

# Function to check that user inputs for 'predict' method are ok
# ----------------------------------------------------------------------------------------------------------
#' @title Fill!!!!!!!!!!!
#' @description Fill!!!!!!!!!!!
#'
#' @param env Fill!!!!!!!!!!!
#'
#' @author José Betancourt, François Bachoc and Thierry Klein
checkVal_simulate <- function(env){
  # recover the model
  model <- env$object

  if (all(!is.null(model@sIn), !is.null(model@fIn))) { # Hybrid-input case *******************************************
    # consistency in data structures
    if (all(is.null(env$sIn.sm), is.null(env$fIn.sm))) {
      stop(paste("Invalid input. The model has both, scalar and functional inputs. Please provide valid new scalar and\n",
                 "functional points to proceed with the simulation."))
    } else if (is.null(env$fIn.sm)) {
      stop(paste("Inconsistent data structures. The model has both, scalar and functional inputs, but only new scalar points\n",
                 "for simulation were specified. Please provide also valid new functional points to proceed with the simulation."))
    } else if (is.null(env$sIn.sm)) {
      stop(paste("Inconsistent data structures. The model has both, scalar and functional inputs, but only new functional points\n",
                 "for simulation were specified. Please provide also valid new scalar points to proceed with the simulation."))
    }

    # consistency in number of points
    if (!all(nrow(env$sIn.sm) == c(sapply(env$fIn.sm, nrow)))) {
      stop("Inconsistent number of points. Please check that sIn.sm and each matrix in fIn have all the same number of rows.")
    }

  } else if(!is.null(model@fIn)) { # functional-input case ***************************************
    # consistency in data structures
    if (!is.null(env$sIn.sm)) {
      stop(paste("Inconsistent data structures. The model has only functional inputs, but new scalar points\n",
                 "for simulation were specified. Please provide new functional points instead to proceed with the simulation."))
    }
    if (is.null(env$fIn.sm)) {
      stop(paste("Invalid input. The model has scalar inputs. Please provide valid new scalar points to proceed with\n",
                 "the simulation."))
    }

  } else if(!is.null(model@sIn)) { # scalar-input case *******************************************
    # consistency in data structures
    if (!is.null(env$fIn.sm)) {
      stop(paste("Inconsistent data structures. The model has only scalar inputs, but new functional points\n",
                 "for simulation were specified. Please provide new scalar points instead to proceed with the simulation."))
    }
    if (is.null(env$sIn.sm)) {
      stop(paste("Invalid input. The model has functional inputs. Please provide valid new functional points to proceed with\n",
                 "the simulation."))
    }
  }
}
# ----------------------------------------------------------------------------------------------------------
