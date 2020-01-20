# ==========================================================================================================
# Class for funGp models
# ==========================================================================================================



# ==========================================================================================================
# Developer oriented methods
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
#' @author José Betancourt
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
# User oriented methods. For documentation of generic methods check the extraDoc.R file
# ==========================================================================================================

# Method to make predictions with afunGp model
# ----------------------------------------------------------------------------------------------------------
#' @name predict
#' @rdname predict-methods
#' @importFrom stats predict
#' @param object An object to predict from.
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
#' plotCalib(m1, m1.preds, sOut.pr)
#'
#' @export predict
setGeneric(name = "predict", def = function(object, ...) standardGeneric("predict"))

#' @importFrom stats qnorm
predict.funGp <- function(object, sIn.pr = NULL, fIn.pr = NULL, detail = "light", ...) {
  # =====================================================================================================
  # Prediction checklist
  # =====================================================================================================
  # 1.  * mean ............... array (n.pr) .............. predicted mean
  # 2.  * sd ................. array (n.pr) .............. predicted standard deviation
  # 3.  * lower95 ............ array (n.pr) .............. lower bounds of 95% confidence intervals
  # 4.  * upper95 ............ array (n.pr) .............. upper bounds of 95% confidence intervals
  # 5.  * K.pp ............... matrix(n.pr x n.pr) ....... conditional covariance matrix
  # 6.  * K.tp ............... matrix(n.tr x n.pr) ....... training vs prediction cross covariance matrix
  # =====================================================================================================

  if (all(object@ds > 0, object@df > 0)) {
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

  } else if (object@df > 0) {
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

  } else {
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

#' @title Prediction method for the funGp Class
#' @name predict
#' @rdname predict-methods
#' @aliases predict,funGp-method
setMethod("predict", "funGp", function(object, ...) predict.funGp(object, ...))
# ----------------------------------------------------------------------------------------------------------


# Method to print a funGp model
# ----------------------------------------------------------------------------------------------------------
#' @name show
#' @description This is my description
#' @rdname show-methods
#' @importFrom methods show
#' @param object An object to show.
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


# Method to get the list of hyperparameters of a funGp model
# ----------------------------------------------------------------------------------------------------------
#' @name getCoef
#' @description This is my description
#' @rdname getCoef-methods
#' @exportMethod getCoef
#' @param object An object to predict from.
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
#' @exportMethod getBasis
#' @param object An object to predict from.
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
#' @exportMethod getRedfIn
#' @param object An object to predict from.
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
#' @exportMethod getProjfIn
#' @param object An object to predict from.
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
#' @exportMethod getProjgram
#' @param object An object to predict from.
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
#' @exportMethod getTrainCov
#' @param object An object to predict from.
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


# Method to plot a funGp model
# ----------------------------------------------------------------------------------------------------------
#' @name plotLOO
#' @description This is my description
#' @rdname plotLOO-methods
#' @exportMethod plotLOO
#' @importFrom graphics lines plot
#' @param object An object to predict from.
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


# Method to plot predictions of a funGp model
# ----------------------------------------------------------------------------------------------------------
#' @name plotPreds
#' @description This is my description
#' @rdname plotPreds-methods
#' @exportMethod plotPreds
#' @importFrom graphics lines plot polygon
#' @param object An object to predict from.
#' @param ... Further arguments for methods.
if(!isGeneric("plotPreds")) {setGeneric("plotPreds", function(object, ...) standardGeneric("plotPreds"))}

#' @title Prediction Method for the apk Class
#' @name plotPreds
#' @rdname plotPreds-methods
#' @aliases plotPreds,funGp-method
#' @param preds something
#' @param sOut.pr also
setMethod("plotPreds", "funGp", function(object, preds, sOut.pr = NULL, ...) plotPreds.funGp(object, preds, sOut.pr))

plotPreds.funGp <- function(object, preds, sOut.pr) {
  y <- sort(preds$mean)
  n.pr <- length(y)
  ll <- (preds$lower95)[order(preds$mean)]
  ul <- (preds$upper95)[order(preds$mean)]

  plot(1, type = "n", xlab = "Index", ylab = "Predicted", xlim = c(1, n.pr), ylim = range(y))
  x <- 1:n.pr
  polygon(c(x, rev(x)), c(ul, rev(ll)), col = "grey85", border = NA)
  if (!is.null(sOut.pr)) {
    lines(sOut.pr[order(preds$mean)], col = "black")
  }
  lines(y, col = "red")
  lines(ll, col = "blue")
  lines(ul, col = "blue")
}


# Method to plot calibration plot from predictions of a funGp model
# ----------------------------------------------------------------------------------------------------------
#' @name plotCalib
#' @description This is my description
#' @rdname plotCalib-methods
#' @exportMethod plotCalib
#' @importFrom graphics lines plot
#' @param object An object to predict from.
#' @param ... Further arguments for methods.
if(!isGeneric("plotCalib")) {setGeneric("plotCalib", function(object, ...) standardGeneric("plotCalib"))}

#' @title Prediction Method for the apk Class
#' @name plotCalib
#' @rdname plotCalib-methods
#' @aliases plotCalib,funGp-method
#' @param preds something
#' @param sOut.pr also
setMethod("plotCalib", "funGp", function(object, preds, sOut.pr, ...) plotCalib.funGp(object, preds, sOut.pr))

plotCalib.funGp <- function(object, preds, sOut.pr) {
  y_obs <- sOut.pr
  y_pre <- preds$mean
  yr <- range(c(y_obs, y_pre))
  plot(y_obs, y_pre, xlim = yr, ylim = yr, pch = 21, col = "red", bg = "red", xlab = "Observed", ylab = "Predicted")
  lines(y_obs, y_obs, col = "blue")
}
