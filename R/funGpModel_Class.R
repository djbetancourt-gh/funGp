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
    preds <- makePreds_SF(sMs.tp, sMs.pp, fMs.tp, fMs.pp, object@kern@varHyp, object@kern@lsHyps[1:object@ds],
                          object@kern@lsHyps[(object@ds+1):(object@ds+object@df)], object@kern@kerType,
                          object@preMats$L, object@preMats$LInvY, detail)

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
    preds <- makePreds_F(fMs.tp, fMs.pp, object@kern@varHyp, object@kern@lsHyps, object@kern@kerType,
                         object@preMats$L, object@preMats$LInvY, detail)

  } else {
    print("I'm scalar!")

    # compute scalar distance matrices
    sMs.tp <- setScalDistance(object@sIn, sIn.pr)
    sMs.pp <- setScalDistance(sIn.pr, sIn.pr)

    # make predictions based on the Gaussian Conditioning Theorem
    preds <- makePreds_S(sMs.tp, sMs.pp, object@kern@varHyp, object@kern@lsHyps, object@kern@kerType,
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
  cat(paste("\nGaussian Process Model", paste(rep("_", min(30, (nchar(callTxt) - nchar(mainTxt) - 1))), collapse = "")))
  # cat("----------------------\n")
  # cat("======================\n\n")

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

  cat(paste("* Do projection: ", object@proj@doProj, "\n", sep = ""))
  if (object@proj@doProj) {
    cat("  -> Proj. dimension:\n")
    for (i in 1:object@df) {
      cat(paste("\t F", i, ": ", object@proj@fpDims[i], "\n", sep = ""))
    }
  }

  cat("\n* Hyperparameters:\n")
  cat(paste("  -> variance: ", format(object@kern@varHyp, digits = 4, nsmall = 4), "\n", sep = ""))
  cat("  -> length-scale:\n")
  ls_s <- object@kern@lsHyps[1:object@ds]
  for (i in 1:object@ds) {
    cat(paste("\t ls(X", i, "): ", format(ls_s[i], digits = 4, nsmall = 4), "\n", sep = ""))
  }
  ls_f <- object@kern@lsHyps[-c(1:object@ds)]
  for (i in 1:object@df) {
    cat(paste("\t ls(F", i, "): ", format(ls_f[i], digits = 4, nsmall = 4), "\n", sep = ""))
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
  coefs <- c(object@kern@varHyp, object@kern@lsHyps)
  names_ls_s <- paste("ls(X", 1:object@ds, ")", sep = "")
  names_ls_f <- paste("ls(F", 1:object@df, ")", sep = "")
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
  return(object@proj@basis)
}


# Method to get the list of gram matrices based on the basis functions used for the projection of the
# functional inputs
# ----------------------------------------------------------------------------------------------------------
#' @name getGram
#' @description This is my description
#' @rdname getGram-methods
#' @exportMethod getGram
#' @param object An object to predict from.
if(!isGeneric("getGram")) {setGeneric(name = "getGram", def = function(object) standardGeneric("getGram"))}

#' @title Prediction Method for the apk Class
#' @name getGram
#' @rdname getGram-methods
#' @aliases getGram,funGp-method
setMethod("getGram", "funGp", function(object) getGram.funGp(object))

getGram.funGp <- function(object) {
  return(lapply(object@proj@basis, crossprod))
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
  rfIn <- object@proj@coefs
  names(rfIn) <- paste("alpha(F", 1:object@df, ")", sep = "")
  return(rfIn)
}


# Method to get the list of projected functional (not used directly in the model)
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
  fpIn <- mapply(function(m1, m2) m1 %*% t(m2), m1 = object@proj@coefs, m2 = object@proj@basis)
  names(fpIn) <- paste("Fp", 1:object@df, sep = "")
  return(fpIn)
}


# Method to get the list of projected functional inputs
# ----------------------------------------------------------------------------------------------------------
if(!isGeneric("getProjfuns")) {setGeneric(name = "getProjfuns", def = function(object) standardGeneric("getProjfuns"))}
#' @rdname getProjfuns
#' @export
setMethod("getProjfuns", "funGp", function(object) getProjfuns.funGp(object))

getProjfuns.funGp <- function(object) {
  if (object@df > 0) {
    fpIn <- list()
    for (i in 1:object@df) {
      fpIn[[i]] <- object@proj@coefs[[i]] %*% t(object@proj@basis[[i]])
    }
    return(fpIn)
  } else {
    cat("The provided funGp model does not have functional inputs. Projected functional inputs are not defined for it.")
  }
}
# ----------------------------------------------------------------------------------------------------------

# Method to get the list of Gram matrices computed from the basis functions
# ----------------------------------------------------------------------------------------------------------
if(!isGeneric("getProjgram")) {setGeneric("getProjgram", function(object) standardGeneric("getProjgram"))}

#' @rdname getProjgram
#' @export
setMethod("getProjgram", "funGp", function(object) getProjgram.funGp(object))

getProjgram.funGp <- function(object) {
  if (object@df > 0) {
    gram <- list()
    for (i in 1:object@df) {
      gram[[i]] <- t(object@proj@basis[[i]]) %*% object@proj@basis[[i]]
    }
    return(gram)
  } else {
    cat("The provided funGp model does not have functional inputs. Thus, the projection Gram matrix is not defined for it.")
  }
}
# ----------------------------------------------------------------------------------------------------------
